#include "lp_backend.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include "resource_scope.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <limits>

namespace {

std::atomic<long long> highs_solves(0);
std::atomic<long long> pivotcopy_solves(0);
std::atomic<long long> pivotcopy_pivots(0);
std::atomic<long long> pivotcopy_warm_solves(0);
std::atomic<long long> pivotcopy_warm_retries(0);
std::atomic<long long> pivotcopy_failures(0);
std::atomic<long long> pivotcopy_fallbacks(0);
std::atomic<long long> cached_solves(0);
std::atomic<long long> highs_nanoseconds(0);
std::atomic<long long> pivotcopy_nanoseconds(0);
std::atomic<int> max_columns(0);
std::atomic<int> max_rows(0);
std::atomic<int> max_highs_columns(0);
std::atomic<int> max_highs_rows(0);
std::atomic<int> max_pivotcopy_columns(0);
std::atomic<int> max_pivotcopy_rows(0);
std::recursive_mutex highs_mutex;

void update_max(std::atomic<int> &target, int value)
{
    int current = target.load(std::memory_order_relaxed);
    while (current < value &&
           !target.compare_exchange_weak(current, value,
                                         std::memory_order_relaxed)) {}
}

long long elapsed_nanoseconds(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start).count();
}

void invalidate(LpModel *model)
{
    model->last_status = NUMFAILURE;
    model->pivotcopy_solution = false;
    model->needs_solve = true;
}

void invalidate_structure(LpModel *model)
{
    invalidate(model);
    ++model->structure_revision;
    // Zero is reserved for callers that cannot provide a trusted revision.
    if (model->structure_revision == 0) model->structure_revision = 1;
}

double row_lower(const LpModel *model, int row)
{
    return model->types[row] == kLpLe ? -kHighsInf : model->rhs[row];
}

double row_upper(const LpModel *model, int row)
{
    return model->types[row] == kLpGe ? kHighsInf : model->rhs[row];
}

} // namespace

std::recursive_mutex& highs_api_mutex() { return highs_mutex; }

void reset_lp_backend_stats()
{
    highs_solves = 0;
    pivotcopy_solves = 0;
    pivotcopy_pivots = 0;
    pivotcopy_warm_solves = 0;
    pivotcopy_warm_retries = 0;
    pivotcopy_failures = 0;
    pivotcopy_fallbacks = 0;
    cached_solves = 0;
    highs_nanoseconds = 0;
    pivotcopy_nanoseconds = 0;
    max_columns = 0;
    max_rows = 0;
    max_highs_columns = 0;
    max_highs_rows = 0;
    max_pivotcopy_columns = 0;
    max_pivotcopy_rows = 0;
}

LpBackendStats last_lp_backend_stats()
{
    LpBackendStats result;
    result.highs_solves = highs_solves.load();
    result.pivotcopy_solves = pivotcopy_solves.load();
    result.pivotcopy_pivots = pivotcopy_pivots.load();
    result.pivotcopy_warm_solves = pivotcopy_warm_solves.load();
    result.pivotcopy_warm_retries = pivotcopy_warm_retries.load();
    result.pivotcopy_failures = pivotcopy_failures.load();
    result.pivotcopy_fallbacks = pivotcopy_fallbacks.load();
    result.cached_solves = cached_solves.load();
    result.max_columns = max_columns.load();
    result.max_rows = max_rows.load();
    result.max_highs_columns = max_highs_columns.load();
    result.max_highs_rows = max_highs_rows.load();
    result.max_pivotcopy_columns = max_pivotcopy_columns.load();
    result.max_pivotcopy_rows = max_pivotcopy_rows.load();
    result.highs_s = highs_nanoseconds.load() * 1e-9;
    result.pivotcopy_s = pivotcopy_nanoseconds.load() * 1e-9;
    return result;
}

std::string selected_lp_backend()
{
    const char *value = std::getenv("DECIDEIT_LP_BACKEND");
    return (!value || std::string(value) == "auto") ? "highs" : value;
}

void highs_require(HighsStatus status, const char *operation)
{
    if (status == HighsStatus::kError)
        throw std::runtime_error(std::string("HiGHS: ") + operation);
}

void configure_highs(Highs &highs, bool output)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    highs_require(highs.setOptionValue("output_flag", output), "output_flag");
    highs_require(highs.setOptionValue("solver", "simplex"), "solver");
    highs_require(highs.setOptionValue("threads", 1), "threads");
    highs_require(highs.setOptionValue("primal_feasibility_tolerance", 1e-8),
                  "primal tolerance");
    highs_require(highs.setOptionValue("dual_feasibility_tolerance", 1e-8),
                  "dual tolerance");
    highs_require(highs.setOptionValue("simplex_iteration_limit", 1000000),
                  "iteration limit");
}

bool highs_solution_valid(const Highs &highs)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    const HighsInfo &info = highs.getInfo();
    const HighsSolution &solution = highs.getSolution();
    if (highs.getModelStatus() != HighsModelStatus::kOptimal ||
        !solution.value_valid || !info.valid ||
        info.primal_solution_status != kSolutionStatusFeasible ||
        info.dual_solution_status != kSolutionStatusFeasible ||
        !std::isfinite(highs.getObjectiveValue()) ||
        !std::isfinite(info.max_primal_infeasibility) ||
        info.max_primal_infeasibility > 1e-8 ||
        !std::isfinite(info.max_dual_infeasibility) ||
        info.max_dual_infeasibility > 1e-8)
        return false;
    return std::all_of(solution.col_value.begin(), solution.col_value.end(),
                       [](double value) { return std::isfinite(value); });
}

LpModel::LpModel(int row_count, int column_count)
    : columns(column_count),
      rows(row_count, std::vector<double>(column_count, 0.0)),
      cost(column_count, 0.0), lower(column_count, 0.0),
      upper(column_count, kHighsInf), rhs(row_count, 0.0),
      types(row_count, kLpLe), backend(selected_lp_backend())
{
    if (backend != "highs")
        throw std::invalid_argument("This build requires the HiGHS local LP backend");
    ensure_highs();
}

LpModel::~LpModel()
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    highs.reset();
}

Highs &LpModel::ensure_highs()
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    if (!highs)
    {
        highs.reset(new Highs());
        configure_highs(*highs);
        highs_require(highs->setOptionValue("presolve", "off"), "presolve");
    }
    return *highs;
}

void LpModel::require_solution() const
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    if (last_status != OPTIMAL)
        throw std::runtime_error("LP solution unavailable");
    if (!pivotcopy_solution && (!highs || !highs_solution_valid(*highs)))
        throw std::runtime_error("HiGHS LP solution unavailable");
}

LpModel *make_lp(int rows, int columns)
{
    if (rows < 0 || columns <= 0)
        throw std::invalid_argument("Invalid LP dimensions");
    return managed::own<LpModel, managed::delete_object<LpModel>>(
        new LpModel(rows, columns));
}

void delete_lp(LpModel *model) { managed::destroy(model); }
void free_lp(LpModel **model) { delete_lp(*model); *model = nullptr; }
int get_Ncolumns(const LpModel *model) { return model->columns; }
int get_Nrows(const LpModel *model) { return static_cast<int>(model->rows.size()); }
double get_infinite(const LpModel *) { return kHighsInf; }

void set_verbose(LpModel *model, int level)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    if (model->highs)
        highs_require(model->highs->setOptionValue("output_flag", level != 0),
                      "output_flag");
}

static void set_sense(LpModel *model, ObjSense sense)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    if (model->sense == sense) return;
    model->sense = sense;
    invalidate(model);
    if (model->loaded)
        highs_require(model->highs->changeObjectiveSense(sense),
                      "objective sense");
}

void set_minim(LpModel *model) { set_sense(model, ObjSense::kMinimize); }
void set_maxim(LpModel *model) { set_sense(model, ObjSense::kMaximize); }

static void update_column_bounds(LpModel *model, int column)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    invalidate_structure(model);
    if (model->loaded)
        highs_require(model->highs->changeColBounds(
                          column, model->lower[column], model->upper[column]),
                      "column bounds");
}

void set_lowbo(LpModel *model, int column, double value)
{
    if (model->lower.at(column - 1) == value) return;
    model->lower.at(column - 1) = value;
    update_column_bounds(model, column - 1);
}

void set_upbo(LpModel *model, int column, double value)
{
    if (model->upper.at(column - 1) == value) return;
    model->upper.at(column - 1) = value;
    update_column_bounds(model, column - 1);
}

static void update_row_bounds(LpModel *model, int row)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    invalidate_structure(model);
    if (model->loaded)
        highs_require(model->highs->changeRowBounds(
                          row, row_lower(model, row), row_upper(model, row)),
                      "row bounds");
}

void set_constr_type(LpModel *model, int row, int type)
{
    if (type != kLpLe && type != kLpGe && type != kLpEq)
        throw std::invalid_argument("row type");
    if (model->types.at(row - 1) == type) return;
    model->types.at(row - 1) = type;
    update_row_bounds(model, row - 1);
}

void set_rh(LpModel *model, int row, double value)
{
    if (!std::isfinite(value)) throw std::invalid_argument("nonfinite rhs");
    if (model->rhs.at(row - 1) == value) return;
    model->rhs.at(row - 1) = value;
    update_row_bounds(model, row - 1);
}

void set_row(LpModel *model, int row, const double *values)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    std::vector<double> &destination = model->rows.at(row - 1);
    bool changed = false;
    for (int column = 0; column < model->columns; ++column)
    {
        if (!std::isfinite(values[column + 1]))
            throw std::invalid_argument("nonfinite coefficient");
        changed = changed || destination[column] != values[column + 1];
    }
    if (!changed) return;
    invalidate_structure(model);
    for (int column = 0; column < model->columns; ++column)
    {
        if (model->loaded && destination[column] != values[column + 1])
            highs_require(model->highs->changeCoeff(
                              row - 1, column, values[column + 1]),
                          "coefficient");
        destination[column] = values[column + 1];
    }
}

void set_rowex(LpModel *model, int row, int count, const double *values,
               const int *indices)
{
    std::vector<double> dense(model->columns + 1, 0.0);
    for (int i = 0; i < count; ++i) dense.at(indices[i]) += values[i];
    set_row(model, row, dense.data());
}

void set_obj_fn(LpModel *model, const double *values)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    bool changed = false;
    for (int column = 0; column < model->columns; ++column)
    {
        if (!std::isfinite(values[column + 1]))
            throw std::invalid_argument("nonfinite objective");
        changed = changed || model->cost[column] != values[column + 1];
    }
    if (!changed) return;
    invalidate(model);
    for (int column = 0; column < model->columns; ++column)
        model->cost[column] = values[column + 1];
    if (model->loaded)
        highs_require(model->highs->changeColsCost(
                          0, model->columns - 1, model->cost.data()),
                      "objective");
}

void set_obj_fnex(LpModel *model, int count, const double *values,
                  const int *indices)
{
    std::vector<double> dense(model->columns + 1, 0.0);
    for (int i = 0; i < count; ++i) dense.at(indices[i]) += values[i];
    set_obj_fn(model, dense.data());
}

void add_constraint(LpModel *model, const double *values, int type, double rhs)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    invalidate_structure(model);
    const bool loaded = model->loaded;
    model->loaded = false;
    model->rows.emplace_back(model->columns, 0.0);
    model->rhs.push_back(rhs);
    model->types.push_back(type);
    set_row(model, get_Nrows(model), values);
    model->loaded = loaded;
    if (loaded)
    {
        std::vector<HighsInt> indices;
        std::vector<double> nonzeros;
        for (int column = 0; column < model->columns; ++column)
            if (values[column + 1] != 0.0)
            {
                indices.push_back(column);
                nonzeros.push_back(values[column + 1]);
            }
        const int row = get_Nrows(model) - 1;
        highs_require(model->highs->addRow(
                          row_lower(model, row), row_upper(model, row),
                          static_cast<HighsInt>(indices.size()), indices.data(),
                          nonzeros.data()),
                      "addRow");
    }
}

static int solve_highs(LpModel *model)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    auto load_model = [model](Highs &highs)
    {
        HighsLp lp;
        lp.num_col_ = model->columns;
        lp.num_row_ = get_Nrows(model);
        lp.col_cost_ = model->cost;
        lp.col_lower_ = model->lower;
        lp.col_upper_ = model->upper;
        lp.sense_ = model->sense;
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.start_.assign(1, 0);
        for (int row = 0; row < lp.num_row_; ++row)
        {
            lp.row_lower_.push_back(row_lower(model, row));
            lp.row_upper_.push_back(row_upper(model, row));
            for (int column = 0; column < lp.num_col_; ++column)
                if (model->rows[row][column] != 0.0)
                {
                    lp.a_matrix_.index_.push_back(column);
                    lp.a_matrix_.value_.push_back(model->rows[row][column]);
                }
            lp.a_matrix_.start_.push_back(lp.a_matrix_.value_.size());
        }
        highs_require(highs.passModel(lp), "passModel");
        model->loaded = true;
    };
    Highs *highs = &model->ensure_highs();
    if (!model->loaded)
    {
        load_model(*highs);
    }
    HighsStatus status = highs->run();
    const HighsModelStatus first_model_status = highs->getModelStatus();
    if ((status == HighsStatus::kError || !highs_solution_valid(*highs)) &&
        first_model_status != HighsModelStatus::kInfeasible &&
        first_model_status != HighsModelStatus::kUnbounded)
    {
        // clearSolver() retains the incrementally modified model and, on a
        // rare numerically difficult reoptimization, can retain the same
        // unusable state.  Reconstructing from our authoritative dense copy
        // is a true cold retry and makes the exported/reloaded LP behaviour
        // available automatically without penalizing successful solves.
        model->highs.reset();
        model->loaded = false;
        highs = &model->ensure_highs();
        highs_require(highs->setOptionValue("presolve", "on"),
                      "numeric retry presolve");
        load_model(*highs);
        status = highs->run();
    }
    model->last_status = NUMFAILURE;
    model->pivotcopy_solution = false;
    if (status != HighsStatus::kError && highs_solution_valid(*highs))
        model->last_status = OPTIMAL;
    else if (highs->getModelStatus() == HighsModelStatus::kInfeasible)
        model->last_status = INFEASIBLE;
    else if (highs->getModelStatus() == HighsModelStatus::kUnbounded)
        model->last_status = UNBOUNDED;
    if (model->last_status == NUMFAILURE && std::getenv("DECIDEIT_TRACE"))
        std::cerr << "TRACE {\"event\":\"highs_numfailure\",\"run_status\":"
                  << static_cast<int>(status) << ",\"model_status\":"
                  << static_cast<int>(highs->getModelStatus())
                  << ",\"rows\":" << get_Nrows(model) << "}\n";
    if (model->last_status != OPTIMAL && std::getenv("DECIDEIT_FAILED_LP"))
        highs_require(highs->writeModel(std::getenv("DECIDEIT_FAILED_LP")),
                      "failed LP diagnostic export");
    return model->last_status;
}

int solve(LpModel *model)
{
    if (!model->needs_solve)
    {
        cached_solves++;
        return model->last_status;
    }
    const int row_count = get_Nrows(model);
    update_max(max_columns, model->columns);
    update_max(max_rows, row_count);
    update_max(max_highs_columns, model->columns);
    update_max(max_highs_rows, row_count);
    const auto start = std::chrono::steady_clock::now();
    const int status = solve_highs(model);
    highs_nanoseconds += elapsed_nanoseconds(start);
    highs_solves++;
    model->needs_solve = false;
    return status;
}

double get_objective(const LpModel *model)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    model->require_solution();
    return model->pivotcopy_solution
        ? model->objective : model->highs->getObjectiveValue();
}

void get_variables(const LpModel *model, double *values)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    model->require_solution();
    if (model->pivotcopy_solution)
        std::copy(model->solution.begin(), model->solution.end(), values);
    else
        std::copy(model->highs->getSolution().col_value.begin(),
                  model->highs->getSolution().col_value.end(), values);
}

void get_constraints(const LpModel *model, double *values)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    model->require_solution();
    if (model->pivotcopy_solution)
        std::copy(model->row_activity.begin(), model->row_activity.end(), values);
    else
        std::copy(model->highs->getSolution().row_value.begin(),
                  model->highs->getSolution().row_value.end(), values);
}

void get_basis(LpModel *model, int *basis, bool include_nonbasic)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    model->require_solution();
    if (include_nonbasic)
        throw std::invalid_argument("nonbasic output not supported");
    basis[0] = 0;
    if (model->pivotcopy_solution)
    {
        if (static_cast<int>(model->row_basis.size()) != get_Nrows(model))
            throw std::runtime_error("Pivot-Copy basis dimension mismatch");
        std::copy(model->row_basis.begin(), model->row_basis.end(), basis + 1);
        return;
    }
    std::vector<HighsInt> variables(get_Nrows(model));
    highs_require(model->highs->getBasicVariables(variables.data()), "basis");
    // HiGHS: j>=0 means structural column j; -(r+1) means logical row r.
    // Geometry convention: row slacks 1..m, then columns m+1..m+n.
    for (size_t i = 0; i < variables.size(); ++i)
        basis[i + 1] = variables[i] >= 0
            ? get_Nrows(model) + variables[i] + 1 : -variables[i];
}

void get_column(const LpModel *model, int column, double *values)
{
    values[0] = model->cost.at(column - 1);
    for (int row = 0; row < get_Nrows(model); ++row)
        values[row + 1] = model->rows[row].at(column - 1);
}

double get_rh(const LpModel *model, int row) { return model->rhs.at(row - 1); }

void clear_lp_solver(LpModel *model)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    if (model->highs) highs_require(model->highs->clearSolver(), "clearSolver");
    invalidate(model);
}

long long get_simplex_iterations(const LpModel *model)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    if (model->pivotcopy_solution) return 0;
    return model->highs
        ? static_cast<long long>(model->highs->getInfo().simplex_iteration_count)
        : 0;
}

void print_lp(const LpModel *model)
{
    std::cout << (model->backend == "highs" ? "HiGHS" : "Pivot-Copy")
              << " LP: " << get_Nrows(model) << " rows, "
              << model->columns << " columns\n";
}
