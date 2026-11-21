#include "dblp.h"
#include "rlt_audit.h"
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <queue>
#include <set>
#include <sstream>

void recover_cut(gsl_vector *onecut, gsl_vector *normal_vector,
                 gsl_vector *xl, gsl_vector_int *fixed_variables);
void linearization_onecut2go(gsl_vector *onecut,
                             LinearizationParams *LinPara,
                             d_row P_lobo, d_row P_upbo,
                             int Ai_begin, int Aj_begin, int n_alt1,
                             const std::vector<LinCons> *opposite_rows = nullptr,
                             const double *separation_values = nullptr,
                             size_t separation_size = 0,
                             const std::vector<int> *forced_cross_rows = nullptr);
void build_linearization(
    int n_alt1, int n_alt2,
    std::vector<LinCons> &linear_cons_P,
    std::vector<LinCons> &linear_cons_V,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    int mX, int mY, LinearizationParams *LinPara);
void putto_global_x_solution(gsl_vector_int *fixed_variables,
                             gsl_vector *whole_solution,
                             gsl_vector *x_solution);
void getfrom_x_solution(gsl_vector_int *fixed_variables,
                        gsl_vector *whole_solution, gsl_vector *xl,
                        gsl_vector *x_solution);

namespace {
thread_local DecideItRunStats decideit_stats;

struct ConePilotResult
{
    bool attempted = false;
    bool certified = false;
    int children = 0;
    int infeasible = 0;
    long long simplex_iterations = 0;
    double lower_bound = std::numeric_limits<double>::infinity();
    double wall_s = 0.0;
};

uint64_t cone_geometry_signature(const ConeRayGeometry &geometry)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    const auto add_bytes = [&](const void *data, size_t length) {
        const unsigned char *bytes = static_cast<const unsigned char *>(data);
        for (size_t i = 0; i < length; ++i)
        {
            hash ^= bytes[i];
            hash *= UINT64_C(1099511628211);
        }
    };
    add_bytes(&geometry.dimension, sizeof(geometry.dimension));
    add_bytes(&geometry.ray_count, sizeof(geometry.ray_count));
    add_bytes(&geometry.alpha, sizeof(geometry.alpha));
    for (double value : geometry.apex) add_bytes(&value, sizeof(value));
    for (double value : geometry.rays) add_bytes(&value, sizeof(value));
    for (double value : geometry.gauges) add_bytes(&value, sizeof(value));
    return hash;
}

int cone_pilot_limit(const char *name, int fallback, int maximum)
{
    const char *text = std::getenv(name);
    if (!text) return fallback;
    char *end = nullptr;
    errno = 0;
    const long value = std::strtol(text, &end, 10);
    if (errno || end == text || *end || value < 1 || value > maximum)
        throw std::invalid_argument(std::string(name) +
                                    " must be in [1," +
                                    std::to_string(maximum) + "]");
    return static_cast<int>(value);
}

double bounded_environment_double(const char *name, double fallback,
                                  double minimum, double maximum)
{
    const char *text = std::getenv(name);
    if (!text) return fallback;
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(text, &end);
    if (errno || end == text || *end || !std::isfinite(value) ||
        value < minimum || value > maximum)
        throw std::invalid_argument(std::string(name) +
                                    " is outside its valid range");
    return value;
}

bool binary_environment_option(const char *name, bool fallback = false)
{
    const char *value = std::getenv(name);
    if (!value) return fallback;
    if (std::strcmp(value, "0") == 0) return false;
    if (std::strcmp(value, "1") == 0) return true;
    throw std::invalid_argument(std::string(name) + " must be 0 or 1");
}

enum class CrossRowRltMode { disabled, violated, all };

CrossRowRltMode selected_cross_row_rlt_mode()
{
    const char *value = std::getenv("DECIDEIT_RLT_CROSS_ROWS");
    if (!value || std::strcmp(value, "0") == 0)
        return CrossRowRltMode::disabled;
    if (std::strcmp(value, "1") == 0 ||
        std::strcmp(value, "violated") == 0)
        return CrossRowRltMode::violated;
    if (std::strcmp(value, "all") == 0)
        return CrossRowRltMode::all;
    throw std::invalid_argument(
        "DECIDEIT_RLT_CROSS_ROWS must be 0, 1, violated, or all");
}

enum class PlateauTightenMode {
    disabled,
    propagation,
    obbt,
    cutoff_obbt,
    reduced_cost,
    reduced_cost_cutoff
};

PlateauTightenMode selected_plateau_tighten_mode()
{
    const char *value = std::getenv("DECIDEIT_RLT_PLATEAU_TIGHTEN");
    if (!value || std::strcmp(value, "0") == 0)
        return PlateauTightenMode::disabled;
    if (std::strcmp(value, "propagate") == 0)
        return PlateauTightenMode::propagation;
    if (std::strcmp(value, "1") == 0 || std::strcmp(value, "obbt") == 0)
        return PlateauTightenMode::obbt;
    if (std::strcmp(value, "cutoff") == 0 ||
        std::strcmp(value, "cutoff-obbt") == 0)
        return PlateauTightenMode::cutoff_obbt;
    if (std::strcmp(value, "rc") == 0 ||
        std::strcmp(value, "reduced-cost") == 0)
        return PlateauTightenMode::reduced_cost;
    if (std::strcmp(value, "rc-cutoff") == 0 ||
        std::strcmp(value, "reduced-cost-cutoff") == 0)
        return PlateauTightenMode::reduced_cost_cutoff;
    throw std::invalid_argument(
        "DECIDEIT_RLT_PLATEAU_TIGHTEN must be 0, propagate, obbt, "
        "cutoff, reduced-cost, or reduced-cost-cutoff");
}

bool reduced_cost_snapshots_requested()
{
    const char *plateau = std::getenv("DECIDEIT_RLT_PLATEAU_TIGHTEN");
    const bool plateau_reduced_cost = plateau &&
        (std::strcmp(plateau, "rc") == 0 ||
         std::strcmp(plateau, "reduced-cost") == 0 ||
         std::strcmp(plateau, "rc-cutoff") == 0 ||
         std::strcmp(plateau, "reduced-cost-cutoff") == 0);
    const char *branch = std::getenv(
        "DECIDEIT_BRANCH_RLT_REDUCED_COST_TIGHTEN");
    return plateau_reduced_cost || (branch && std::strcmp(branch, "1") == 0);
}

void configure_rlt_highs(Highs &solver)
{
    configure_highs(solver);
    if (const char *value = std::getenv("DECIDEIT_RLT_SOLVER"))
        highs_require(solver.setOptionValue("solver", value), "RLT solver");
    if (const char *value = std::getenv("DECIDEIT_RLT_PRESOLVE"))
        highs_require(solver.setOptionValue("presolve", value),
                      "RLT presolve");
    if (const char *value = std::getenv("DECIDEIT_RLT_SIMPLEX_STRATEGY"))
        highs_require(solver.setOptionValue(
                          "simplex_strategy", std::stoi(value)),
                      "RLT simplex strategy");
    if (const char *value = std::getenv(
            "DECIDEIT_RLT_SIMPLEX_SCALE_STRATEGY"))
        if (std::strcmp(value, "auto") != 0 &&
            std::strcmp(value, "adaptive") != 0)
            highs_require(solver.setOptionValue(
                              "simplex_scale_strategy", std::stoi(value)),
                          "RLT simplex scale strategy");
    if (const char *value = std::getenv(
            "DECIDEIT_RLT_DUAL_EDGE_WEIGHT_STRATEGY"))
        highs_require(solver.setOptionValue(
                          "simplex_dual_edge_weight_strategy",
                          std::stoi(value)),
                      "RLT dual edge weight strategy");
    if (const char *value = std::getenv(
            "DECIDEIT_RLT_PRIMAL_EDGE_WEIGHT_STRATEGY"))
        highs_require(solver.setOptionValue(
                          "simplex_primal_edge_weight_strategy",
                          std::stoi(value)),
                      "RLT primal edge weight strategy");
}

void configure_rlt_scale_for_columns(Highs &solver, int columns)
{
    const char *value = std::getenv(
        "DECIDEIT_RLT_SIMPLEX_SCALE_STRATEGY");
    const bool default_auto = !value || std::strcmp(value, "auto") == 0;
    const bool large_only = value && std::strcmp(value, "adaptive") == 0;
    if (!default_auto && !large_only) return;
    const int minimum = large_only ? cone_pilot_limit(
        "DECIDEIT_RLT_SCALE4_MIN_COLUMNS", 700, 4000000) : 1;
    const int maximum = cone_pilot_limit(
        "DECIDEIT_RLT_SCALE4_MAX_COLUMNS", 2500, 4000000);
    highs_require(solver.setOptionValue(
                      "simplex_scale_strategy",
                      columns >= minimum && columns <= maximum ? 4 : 2),
                  "size-aware RLT simplex scale strategy");
}

struct BoundRow
{
    std::vector<double> coefficient;
    double rhs = 0.0;
};

struct PlateauTighteningResult
{
    bool feasible = true;
    int propagation_updates = 0;
    int obbt_solves = 0;
    int obbt_updates = 0;
    int reduced_cost_updates = 0;
    int selected_x = 0;
    int selected_y = 0;
    bool reduced_cost_available = false;
    bool cutoff_certified = false;
    bool cutoff_budget_exhausted = false;
    long long simplex_iterations = 0;
    double width_reduction = 0.0;
    double wall_s = 0.0;
};

BoundRow dense_bound_row(const LinCons &source, int variables,
                         double multiplier = 1.0)
{
    BoundRow row;
    row.coefficient.assign(variables, 0.0);
    row.rhs = multiplier * source.b;
    for (size_t k = 0; k < source.index.size(); ++k)
    {
        const int index = source.index[k] - 1;
        if (index < 0 || index >= variables)
            throw std::runtime_error("Bound-tightening row index out of range");
        row.coefficient[index] += multiplier * source.coefficient[k];
    }
    return row;
}

std::vector<BoundRow> bound_rows_from_constraints(
    const std::vector<LinCons> &constraints, int variables)
{
    std::vector<BoundRow> rows;
    rows.reserve(constraints.size());
    for (const LinCons &source : constraints)
        rows.push_back(dense_bound_row(source, variables));
    return rows;
}

void append_probability_normalizations(std::vector<BoundRow> *rows,
                                       int first_block, int variables)
{
    for (int block = 0; block < 2; ++block)
    {
        BoundRow positive, negative;
        positive.coefficient.assign(variables, 0.0);
        negative.coefficient.assign(variables, 0.0);
        const int begin = block == 0 ? 0 : first_block;
        const int end = block == 0 ? first_block : variables;
        for (int column = begin; column < end; ++column)
        {
            positive.coefficient[column] = 1.0;
            negative.coefficient[column] = -1.0;
        }
        positive.rhs = 1.0;
        negative.rhs = -1.0;
        rows->push_back(positive);
        rows->push_back(negative);
    }
}

bool propagate_linear_bounds(const std::vector<BoundRow> &rows,
                             std::vector<double> *lower,
                             std::vector<double> *upper,
                             int maximum_passes, int *updates,
                             double *width_reduction)
{
    const double tolerance = 1e-10;
    for (int pass = 0; pass < maximum_passes; ++pass)
    {
        bool changed = false;
        for (const BoundRow &row : rows)
        {
            double minimum_activity = 0.0;
            for (size_t column = 0; column < row.coefficient.size(); ++column)
            {
                const double value = row.coefficient[column];
                minimum_activity += value >= 0.0
                    ? value * (*lower)[column]
                    : value * (*upper)[column];
            }
            if (!std::isfinite(minimum_activity)) continue;
            if (minimum_activity > row.rhs + 1e-8)
                return false;
            for (size_t column = 0; column < row.coefficient.size(); ++column)
            {
                const double value = row.coefficient[column];
                if (std::fabs(value) <= 1e-14) continue;
                const double selected = value >= 0.0
                    ? (*lower)[column] : (*upper)[column];
                const double other_minimum =
                    minimum_activity - value * selected;
                double candidate = (row.rhs - other_minimum) / value;
                const double guard = tolerance * (1.0 + std::fabs(candidate));
                if (value > 0.0)
                {
                    candidate += guard;
                    if (candidate < (*upper)[column] - tolerance)
                    {
                        if (candidate < (*lower)[column] - 1e-8) return false;
                        const double old = (*upper)[column];
                        (*upper)[column] = (std::max)(candidate,
                                                     (*lower)[column]);
                        *width_reduction += old - (*upper)[column];
                        ++*updates;
                        changed = true;
                    }
                }
                else
                {
                    candidate -= guard;
                    if (candidate > (*lower)[column] + tolerance)
                    {
                        if (candidate > (*upper)[column] + 1e-8) return false;
                        const double old = (*lower)[column];
                        (*lower)[column] = (std::min)(candidate,
                                                     (*upper)[column]);
                        *width_reduction += (*lower)[column] - old;
                        ++*updates;
                        changed = true;
                    }
                }
            }
        }
        if (!changed) break;
    }
    return true;
}

struct LpModelDeleter
{
    void operator()(LpModel *model) const { delete_lp(model); }
};

bool obbt_selected_variables(const std::vector<BoundRow> &rows,
                             const std::vector<int> &candidates,
                             std::vector<double> *lower,
                             std::vector<double> *upper,
                             PlateauTighteningResult *result)
{
    if (candidates.empty()) return true;
    const int variables = static_cast<int>(lower->size());
    std::unique_ptr<LpModel, LpModelDeleter> model(
        make_lp(static_cast<int>(rows.size()), variables));
    set_verbose(model.get(), kLpNeutral);
    std::vector<double> dense(variables + 1, 0.0);
    for (size_t row = 0; row < rows.size(); ++row)
    {
        std::fill(dense.begin(), dense.end(), 0.0);
        for (int column = 0; column < variables; ++column)
            dense[column + 1] = rows[row].coefficient[column];
        set_row(model.get(), static_cast<int>(row + 1), dense.data());
        set_constr_type(model.get(), static_cast<int>(row + 1), kLpLe);
        set_rh(model.get(), static_cast<int>(row + 1), rows[row].rhs);
    }
    for (int column = 0; column < variables; ++column)
    {
        set_lowbo(model.get(), column + 1, (*lower)[column]);
        set_upbo(model.get(), column + 1, (*upper)[column]);
    }

    for (int column : candidates)
    {
        if (column < 0 || column >= variables)
            throw std::runtime_error("OBBT candidate index out of range");
        if ((*upper)[column] - (*lower)[column] <= 1e-10) continue;
        std::fill(dense.begin(), dense.end(), 0.0);
        dense[column + 1] = 1.0;
        set_obj_fn(model.get(), dense.data());

        set_minim(model.get());
        int status = solve(model.get());
        result->obbt_solves++;
        result->simplex_iterations += get_simplex_iterations(model.get());
        if (status == INFEASIBLE) return false;
        if (status == OPTIMAL)
        {
            const double optimum = get_objective(model.get());
            const double candidate = optimum -
                1e-9 * (1.0 + std::fabs(optimum));
            if (candidate > (*lower)[column] + 1e-10)
            {
                const double old = (*lower)[column];
                (*lower)[column] = (std::min)(candidate, (*upper)[column]);
                result->width_reduction += (*lower)[column] - old;
                result->obbt_updates++;
                set_lowbo(model.get(), column + 1, (*lower)[column]);
            }
        }

        set_maxim(model.get());
        status = solve(model.get());
        result->obbt_solves++;
        result->simplex_iterations += get_simplex_iterations(model.get());
        if (status == INFEASIBLE) return false;
        if (status == OPTIMAL)
        {
            const double optimum = get_objective(model.get());
            const double candidate = optimum +
                1e-9 * (1.0 + std::fabs(optimum));
            if (candidate < (*upper)[column] - 1e-10)
            {
                const double old = (*upper)[column];
                (*upper)[column] = (std::max)(candidate, (*lower)[column]);
                result->width_reduction += old - (*upper)[column];
                result->obbt_updates++;
                set_upbo(model.get(), column + 1, (*upper)[column]);
            }
        }
    }
    return true;
}

std::pair<std::vector<int>, std::vector<int> > residual_candidates(
    const LinearizationParams &parameters, const gsl_vector *solution,
    int variables, int limit)
{
    if (limit >= 2 * variables)
    {
        std::vector<int> all(variables);
        for (int index = 0; index < variables; ++index) all[index] = index;
        return {all, all};
    }
    struct Candidate { double residual; int x; int y; };
    std::vector<Candidate> pairs;
    pairs.reserve(static_cast<size_t>(variables) * variables);
    for (int i = 0; i < variables; ++i)
        for (int j = 0; j < variables; ++j)
        {
            const int product = 2 * variables + i * variables + j;
            if (product >= parameters.num_col) continue;
            const double residual = std::fabs(
                gsl_vector_get(solution, product) -
                gsl_vector_get(solution, i) *
                gsl_vector_get(solution, variables + j));
            if (residual > 10.0 * DBLP_EPSILON)
                pairs.push_back(Candidate{residual, i, j});
        }
    std::sort(pairs.begin(), pairs.end(),
              [](const Candidate &left, const Candidate &right) {
                  return left.residual > right.residual;
              });
    std::set<int> selected_x, selected_y;
    for (const Candidate &candidate : pairs)
    {
        if (static_cast<int>(selected_x.size() + selected_y.size()) >= limit)
            break;
        selected_x.insert(candidate.x);
        if (static_cast<int>(selected_x.size() + selected_y.size()) < limit)
            selected_y.insert(candidate.y);
    }
    return {std::vector<int>(selected_x.begin(), selected_x.end()),
            std::vector<int>(selected_y.begin(), selected_y.end())};
}

bool objective_cutoff_rlt_obbt(
    const LinearizationParams &parameters,
    const std::vector<int> &x_candidates,
    const std::vector<int> &y_candidates,
    int variables, double incumbent,
    std::vector<double> *x_lower, std::vector<double> *x_upper,
    std::vector<double> *y_lower, std::vector<double> *y_upper,
    PlateauTighteningResult *result)
{
    std::vector<int> candidates = x_candidates;
    for (int column : y_candidates)
        candidates.push_back(variables + column);
    if (candidates.empty()) return true;
    const int total_limit_ms = cone_pilot_limit(
        "DECIDEIT_RLT_CUTOFF_OBBT_TIME_LIMIT_MS", 3000, 60000);
    const auto started = std::chrono::steady_clock::now();
    const auto remaining_seconds = [&]() {
        const long long elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - started).count();
        return (total_limit_ms - elapsed) / 1000.0;
    };

    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    Highs solver;
    configure_rlt_highs(solver);
    highs_require(solver.setOptionValue("output_flag", false),
                  "cutoff OBBT output flag");
    HighsLp lp;
    lp.num_col_ = parameters.num_col;
    lp.num_row_ = parameters.num_row;
    lp.col_cost_.assign(parameters.num_col, 0.0);
    lp.col_lower_.assign(parameters.col_lower->data,
                         parameters.col_lower->data + parameters.num_col);
    lp.col_upper_.assign(parameters.col_upper->data,
                         parameters.col_upper->data + parameters.num_col);
    // Interval propagation precedes cutoff OBBT.  Apply those bounds to the
    // copied relaxation immediately rather than waiting until a selected
    // variable happens to be optimized.
    for (int column = 0; column < variables; ++column)
    {
        lp.col_lower_[column] = (*x_lower)[column];
        lp.col_upper_[column] = (*x_upper)[column];
        lp.col_lower_[variables + column] = (*y_lower)[column];
        lp.col_upper_[variables + column] = (*y_upper)[column];
    }
    lp.row_lower_.assign(parameters.row_lower->data,
                         parameters.row_lower->data + parameters.num_row);
    lp.row_upper_.assign(parameters.row_upper->data,
                         parameters.row_upper->data + parameters.num_row);
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.start_.assign(parameters.a_start->data,
        parameters.a_start->data + parameters.num_row + 1);
    lp.a_matrix_.index_.assign(parameters.a_index->data,
        parameters.a_index->data + parameters.num_nz);
    lp.a_matrix_.value_.assign(parameters.a_value->data,
        parameters.a_value->data + parameters.num_nz);
    lp.sense_ = ObjSense::kMinimize;
    if (solver.passModel(lp) == HighsStatus::kError) return false;
    configure_rlt_scale_for_columns(solver, parameters.num_col);

    if (parameters.persistent)
    {
        const HighsBasis &basis = parameters.persistent->highs.getBasis();
        if (basis.valid && basis.col_status.size() ==
                static_cast<size_t>(parameters.num_col) &&
            basis.row_status.size() ==
                static_cast<size_t>(parameters.num_row))
            solver.setBasis(basis);
    }
    std::vector<HighsInt> cutoff_index;
    std::vector<double> cutoff_value;
    for (int column = 0; column < parameters.num_col; ++column)
    {
        const double value = gsl_vector_get(parameters.col_cost, column);
        if (std::fabs(value) > 1e-15)
        {
            cutoff_index.push_back(column);
            cutoff_value.push_back(value);
        }
    }
    const double cutoff = incumbent - DBLP_EPSILON - parameters.offset;
    if (solver.addRow(-1.0e30, cutoff,
                      static_cast<HighsInt>(cutoff_index.size()),
                      cutoff_index.data(), cutoff_value.data()) ==
        HighsStatus::kError)
        return false;

    int previous_column = -1;
    for (int column : candidates)
    {
        if (column < 0 || column >= 2 * variables)
            throw std::runtime_error(
                "Cutoff RLT OBBT candidate index out of range");
        if (remaining_seconds() <= 0.0)
        {
            result->cutoff_budget_exhausted = true;
            break;
        }
        if (previous_column >= 0)
            highs_require(solver.changeColCost(previous_column, 0.0),
                          "cutoff OBBT clear objective");
        highs_require(solver.changeColCost(column, 1.0),
                      "cutoff OBBT set objective");
        previous_column = column;

        for (int sense = 0; sense < 2; ++sense)
        {
            const double remaining = remaining_seconds();
            if (remaining <= 0.0)
            {
                result->cutoff_budget_exhausted = true;
                break;
            }
            highs_require(solver.setOptionValue("time_limit", remaining),
                          "cutoff OBBT time limit");
            highs_require(solver.changeObjectiveSense(
                sense == 0 ? ObjSense::kMinimize : ObjSense::kMaximize),
                "cutoff OBBT objective sense");
            const HighsStatus run_status = solver.run();
            result->obbt_solves++;
            const HighsInt iterations =
                solver.getInfo().simplex_iteration_count;
            if (iterations > 0)
                result->simplex_iterations +=
                    static_cast<long long>(iterations);
            const HighsModelStatus model_status = solver.getModelStatus();
            if (model_status == HighsModelStatus::kInfeasible)
            {
                result->cutoff_certified = true;
                return true;
            }
            if (run_status == HighsStatus::kError ||
                model_status != HighsModelStatus::kOptimal ||
                !highs_solution_valid(solver))
            {
                result->cutoff_budget_exhausted = true;
                return true;
            }
            const double optimum = solver.getObjectiveValue();
            const double guard = 1e-9 * (1.0 + std::fabs(optimum));
            std::vector<double> *lower = column < variables
                ? x_lower : y_lower;
            std::vector<double> *upper = column < variables
                ? x_upper : y_upper;
            const int local = column < variables
                ? column : column - variables;
            if (sense == 0)
            {
                const double candidate = optimum - guard;
                if (candidate > (*lower)[local] + 1e-10)
                {
                    const double old = (*lower)[local];
                    (*lower)[local] = (std::min)(candidate, (*upper)[local]);
                    result->width_reduction += (*lower)[local] - old;
                    result->obbt_updates++;
                    highs_require(solver.changeColBounds(
                        column, (*lower)[local], (*upper)[local]),
                        "cutoff OBBT lower bound");
                }
            }
            else
            {
                const double candidate = optimum + guard;
                if (candidate < (*upper)[local] - 1e-10)
                {
                    const double old = (*upper)[local];
                    (*upper)[local] = (std::max)(candidate, (*lower)[local]);
                    result->width_reduction += old - (*upper)[local];
                    result->obbt_updates++;
                    highs_require(solver.changeColBounds(
                        column, (*lower)[local], (*upper)[local]),
                        "cutoff OBBT upper bound");
                }
            }
        }
        if (result->cutoff_budget_exhausted) break;
    }
    return true;
}

bool reduced_cost_rlt_tightening(
    const LinearizationParams &parameters,
    const gsl_vector *solution, int variables, int candidate_limit,
    double incumbent,
    std::vector<double> *x_lower, std::vector<double> *x_upper,
    std::vector<double> *y_lower, std::vector<double> *y_upper,
    PlateauTighteningResult *result)
{
    if (!parameters.persistent) return true;
    const HighsSolver &solver = *parameters.persistent;
    if (!solver.reduced_cost_snapshot_valid ||
        !std::isfinite(solver.reduced_cost_objective) ||
        solver.reduced_cost_primal.size() !=
            static_cast<size_t>(parameters.num_col) ||
        solver.reduced_cost_dual.size() !=
            static_cast<size_t>(parameters.num_col) ||
        !solver.reduced_cost_basis.valid ||
        solver.reduced_cost_basis.col_status.size() !=
            static_cast<size_t>(parameters.num_col))
        return true;

    result->reduced_cost_available = true;
    const double cutoff = incumbent - DBLP_EPSILON - parameters.offset;
    const double objective_guard = 1e-7 *
        (1.0 + std::fabs(cutoff) +
         std::fabs(solver.reduced_cost_objective));
    // A slightly enlarged allowance and outward-rounded variable bounds make
    // the inference conservative under the LP solver's dual tolerance.
    const double allowance = (std::max)(
        0.0, cutoff - solver.reduced_cost_objective) + objective_guard;
    const double dual_tolerance = 1e-8;

    // Residual-only selection often chooses basic columns, for which reduced
    // costs cannot imply a bound.  Rank every eligible nonbasic marginal by
    // the shrink it can actually deliver.  The default hybrid score also
    // favours variables participating in the largest rank-one violations,
    // since tightening an inactive part of the McCormick envelope can make a
    // rebuilt RLT LP slower without reducing the spatial-search tree.
    std::vector<double> x_residual(variables, 0.0);
    std::vector<double> y_residual(variables, 0.0);
    double maximum_residual = 0.0;
    for (int i = 0; i < variables; ++i)
        for (int j = 0; j < variables; ++j)
        {
            const int product = 2 * variables + i * variables + j;
            if (product >= parameters.num_col) continue;
            const double residual = std::fabs(
                gsl_vector_get(solution, product) -
                gsl_vector_get(solution, i) *
                gsl_vector_get(solution, variables + j));
            x_residual[i] = (std::max)(x_residual[i], residual);
            y_residual[j] = (std::max)(y_residual[j], residual);
            maximum_residual = (std::max)(maximum_residual, residual);
        }
    double x_width_sum = 0.0, y_width_sum = 0.0;
    for (int i = 0; i < variables; ++i)
    {
        x_width_sum += (std::max)(0.0, (*x_upper)[i] - (*x_lower)[i]);
        y_width_sum += (std::max)(0.0, (*y_upper)[i] - (*y_lower)[i]);
    }
    const char *rank_value = std::getenv(
        "DECIDEIT_RLT_REDUCED_COST_RANK");
    const std::string rank_mode = rank_value ? rank_value : "hybrid";
    if (rank_mode != "width" && rank_mode != "residual" &&
        rank_mode != "hybrid")
        throw std::invalid_argument(
            "DECIDEIT_RLT_REDUCED_COST_RANK must be width, residual, or "
            "hybrid");

    struct ReducedCostCandidate
    {
        int column;
        double bound;
        double reduction;
        double score;
        bool upper;
    };
    std::vector<ReducedCostCandidate> candidates;
    candidates.reserve(2 * variables);
    for (int column = 0; column < 2 * variables; ++column)
    {
        std::vector<double> *lower = column < variables
            ? x_lower : y_lower;
        std::vector<double> *upper = column < variables
            ? x_upper : y_upper;
        const int local = column < variables
            ? column : column - variables;
        if ((*upper)[local] - (*lower)[local] <= 1e-10) continue;

        const double primal = solver.reduced_cost_primal[column];
        const double dual = solver.reduced_cost_dual[column];
        const double model_lower = gsl_vector_get(
            parameters.col_lower, column);
        const double model_upper = gsl_vector_get(
            parameters.col_upper, column);
        const HighsBasisStatus status =
            solver.reduced_cost_basis.col_status[column];
        double candidate = 0.0;
        double reduction = 0.0;
        bool tighten_upper = false;
        if (status == HighsBasisStatus::kLower &&
            std::isfinite(model_lower) && dual > dual_tolerance &&
            std::fabs(primal - model_lower) <=
                1e-7 * (1.0 + std::fabs(model_lower)))
        {
            candidate = model_lower + allowance / dual;
            candidate += 1e-9 * (1.0 + std::fabs(candidate));
            candidate = (std::max)(candidate, (*lower)[local]);
            reduction = (*upper)[local] - candidate;
            tighten_upper = true;
        }
        else if (status == HighsBasisStatus::kUpper &&
                 std::isfinite(model_upper) && dual < -dual_tolerance &&
                 std::fabs(primal - model_upper) <=
                    1e-7 * (1.0 + std::fabs(model_upper)))
        {
            candidate = model_upper - allowance / (-dual);
            candidate -= 1e-9 * (1.0 + std::fabs(candidate));
            candidate = (std::min)(candidate, (*upper)[local]);
            reduction = candidate - (*lower)[local];
        }
        if (reduction <= 1e-10) continue;

        const double residual = column < variables
            ? x_residual[local] : y_residual[local];
        const double envelope_scale = column < variables
            ? y_width_sum : x_width_sum;
        double score = reduction * envelope_scale;
        if (rank_mode == "residual")
            score = reduction * residual;
        else if (rank_mode == "hybrid")
        {
            const double residual_weight = maximum_residual > 0.0
                ? residual / maximum_residual : 0.0;
            score *= 0.25 + 0.75 * residual_weight;
        }
        candidates.push_back(ReducedCostCandidate{
            column, candidate, reduction, score, tighten_upper});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const ReducedCostCandidate &left,
                 const ReducedCostCandidate &right) {
                  return left.score != right.score
                      ? left.score > right.score
                      : left.reduction > right.reduction;
              });
    if (candidate_limit < static_cast<int>(candidates.size()))
        candidates.resize(candidate_limit);
    for (const ReducedCostCandidate &candidate : candidates)
    {
        const int column = candidate.column;
        std::vector<double> *lower = column < variables
            ? x_lower : y_lower;
        std::vector<double> *upper = column < variables
            ? x_upper : y_upper;
        const int local = column < variables
            ? column : column - variables;
        if (candidate.upper)
            (*upper)[local] = candidate.bound;
        else
            (*lower)[local] = candidate.bound;
        result->width_reduction += candidate.reduction;
        result->reduced_cost_updates++;
        if (column < variables) result->selected_x++;
        else result->selected_y++;
    }
    return true;
}

PlateauTighteningResult tighten_plateau_bounds(
    PlateauTightenMode mode, const LinearizationParams &parameters,
    const gsl_vector *solution, int variables, int first_block,
    const std::vector<LinCons> &x_constraints,
    const std::vector<LinCons> &y_constraints,
    const dataset *dataset_in, int noc,
    gsl_vector_int *fixed_variables,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Aj_begin, int candidate_limit, double incumbent)
{
    const auto started = std::chrono::steady_clock::now();
    PlateauTighteningResult result;
    const int propagation_passes = cone_pilot_limit(
        "DECIDEIT_RLT_PROPAGATION_PASSES", 8, 1000);
    const int followup_propagation_passes = cone_pilot_limit(
        "DECIDEIT_RLT_FOLLOWUP_PROPAGATION_PASSES", 4, 1000);
    std::vector<double> x_lower(variables), x_upper(variables);
    std::vector<double> y_lower(variables), y_upper(variables);
    for (int i = 0; i < variables; ++i)
    {
        const int original = location(i, Ai_begin, Aj_begin, first_block);
        x_lower[i] = V_lobo[original]; x_upper[i] = V_upbo[original];
        y_lower[i] = P_lobo[original]; y_upper[i] = P_upbo[original];
    }

    std::vector<BoundRow> x_rows =
        bound_rows_from_constraints(x_constraints, variables);
    gsl_vector *cut = managed::vector_calloc(variables + 1);
    for (int row = dataset_in->ori_num_cons; row < noc; ++row)
    {
        gsl_vector_view source = gsl_matrix_row(dataset_in->valid_A, row);
        recover_cut(cut, &source.vector, nullptr, fixed_variables);
        BoundRow recovered;
        recovered.coefficient.resize(variables);
        recovered.rhs = gsl_vector_get(cut, 0);
        for (int column = 0; column < variables; ++column)
            recovered.coefficient[column] = gsl_vector_get(cut, column + 1);
        x_rows.push_back(recovered);
    }
    managed::vector_free(cut);
    std::vector<BoundRow> y_rows =
        bound_rows_from_constraints(y_constraints, variables);
    append_probability_normalizations(&y_rows, first_block, variables);

    result.feasible = propagate_linear_bounds(
        x_rows, &x_lower, &x_upper, propagation_passes,
        &result.propagation_updates,
        &result.width_reduction) &&
        propagate_linear_bounds(
            y_rows, &y_lower, &y_upper, propagation_passes,
            &result.propagation_updates, &result.width_reduction);

    if (result.feasible &&
        (mode == PlateauTightenMode::obbt ||
         mode == PlateauTightenMode::cutoff_obbt ||
         mode == PlateauTightenMode::reduced_cost ||
         mode == PlateauTightenMode::reduced_cost_cutoff))
    {
        if (mode == PlateauTightenMode::obbt)
        {
            const auto candidates = residual_candidates(
                parameters, solution, variables, candidate_limit);
            result.selected_x = static_cast<int>(candidates.first.size());
            result.selected_y = static_cast<int>(candidates.second.size());
            result.feasible = obbt_selected_variables(
                x_rows, candidates.first, &x_lower, &x_upper, &result) &&
                obbt_selected_variables(
                    y_rows, candidates.second, &y_lower, &y_upper, &result);
        }
        else if (mode == PlateauTightenMode::cutoff_obbt)
        {
            const auto candidates = residual_candidates(
                parameters, solution, variables, candidate_limit);
            result.selected_x = static_cast<int>(candidates.first.size());
            result.selected_y = static_cast<int>(candidates.second.size());
            result.feasible = objective_cutoff_rlt_obbt(
                parameters, candidates.first, candidates.second,
                variables, incumbent,
                &x_lower, &x_upper, &y_lower, &y_upper, &result);
        }
        else
        {
            result.feasible = reduced_cost_rlt_tightening(
                parameters, solution, variables, candidate_limit, incumbent,
                &x_lower, &x_upper, &y_lower, &y_upper, &result);
            if (result.feasible &&
                mode == PlateauTightenMode::reduced_cost_cutoff)
            {
                const auto candidates = residual_candidates(
                    parameters, solution, variables, candidate_limit);
                result.feasible = objective_cutoff_rlt_obbt(
                    parameters, candidates.first, candidates.second,
                    variables, incumbent,
                    &x_lower, &x_upper, &y_lower, &y_upper, &result);
            }
        }
        if (result.feasible && !result.cutoff_certified)
            result.feasible = propagate_linear_bounds(
                x_rows, &x_lower, &x_upper, followup_propagation_passes,
                &result.propagation_updates, &result.width_reduction) &&
                propagate_linear_bounds(
                    y_rows, &y_lower, &y_upper,
                    followup_propagation_passes,
                    &result.propagation_updates, &result.width_reduction);
    }

    if (result.feasible)
        for (int i = 0; i < variables; ++i)
        {
            const int original = location(i, Ai_begin, Aj_begin, first_block);
            V_lobo[original] = x_lower[i]; V_upbo[original] = x_upper[i];
            P_lobo[original] = y_lower[i]; P_upbo[original] = y_upper[i];
        }
    result.wall_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return result;
}

bool propagate_branch_rlt_bounds(
    std::vector<double> *combined_lower,
    std::vector<double> *combined_upper,
    int variables, int first_block,
    const std::vector<LinCons> &x_constraints,
    const std::vector<LinCons> &y_constraints,
    const dataset *dataset_in, int noc,
    gsl_vector_int *fixed_variables)
{
    if (combined_lower->size() != static_cast<size_t>(2 * variables) ||
        combined_upper->size() != static_cast<size_t>(2 * variables))
        throw std::runtime_error("Branch RLT bound dimension mismatch");
    std::vector<double> x_lower(combined_lower->begin(),
                                combined_lower->begin() + variables);
    std::vector<double> x_upper(combined_upper->begin(),
                                combined_upper->begin() + variables);
    std::vector<double> y_lower(combined_lower->begin() + variables,
                                combined_lower->end());
    std::vector<double> y_upper(combined_upper->begin() + variables,
                                combined_upper->end());
    std::vector<BoundRow> x_rows =
        bound_rows_from_constraints(x_constraints, variables);
    gsl_vector *cut = managed::vector_calloc(variables + 1);
    for (int row = dataset_in->ori_num_cons; row < noc; ++row)
    {
        gsl_vector_view source = gsl_matrix_row(dataset_in->valid_A, row);
        recover_cut(cut, &source.vector, nullptr, fixed_variables);
        BoundRow recovered;
        recovered.coefficient.resize(variables);
        recovered.rhs = gsl_vector_get(cut, 0);
        for (int column = 0; column < variables; ++column)
            recovered.coefficient[column] =
                gsl_vector_get(cut, column + 1);
        x_rows.push_back(recovered);
    }
    managed::vector_free(cut);
    std::vector<BoundRow> y_rows =
        bound_rows_from_constraints(y_constraints, variables);
    append_probability_normalizations(&y_rows, first_block, variables);
    int updates = 0;
    double width_reduction = 0.0;
    if (!propagate_linear_bounds(x_rows, &x_lower, &x_upper, 8,
                                 &updates, &width_reduction) ||
        !propagate_linear_bounds(y_rows, &y_lower, &y_upper, 8,
                                 &updates, &width_reduction))
        return false;
    for (int column = 0; column < variables; ++column)
    {
        (*combined_lower)[column] = x_lower[column];
        (*combined_upper)[column] = x_upper[column];
        (*combined_lower)[variables + column] = y_lower[column];
        (*combined_upper)[variables + column] = y_upper[column];
    }
    return true;
}

bool branch_rlt_propagation_enabled()
{
    const char *value = std::getenv("DECIDEIT_BRANCH_RLT_PROPAGATE");
    if (!value || std::strcmp(value, "0") == 0) return false;
    if (std::strcmp(value, "1") == 0) return true;
    throw std::invalid_argument(
        "DECIDEIT_BRANCH_RLT_PROPAGATE must be 0 or 1");
}

std::string cone_pilot_file(const char *directory, bool min_flag,
                            const char *extension)
{
    std::string path(directory ? directory : "");
    if (!path.empty() && path.back() != '/') path.push_back('/');
    path += min_flag ? "min" : "max";
    path += extension;
    return path;
}

bool export_cone_geometry(const ConeRayGeometry &geometry, bool min_flag,
                          const char *directory)
{
    if (!directory || !geometry.valid) return false;
    std::ofstream stream(cone_pilot_file(directory, min_flag, ".cone"));
    if (!stream) return false;
    stream << std::setprecision(17);
    stream << "DECIDEIT_CONE 1\n";
    stream << "dimension " << geometry.dimension << '\n';
    stream << "ray_count " << geometry.ray_count << '\n';
    stream << "alpha " << geometry.alpha << '\n';
    stream << "signature " << cone_geometry_signature(geometry) << '\n';
    stream << "apex";
    for (double value : geometry.apex) stream << ' ' << value;
    stream << '\n';
    for (int ray = 0; ray < geometry.ray_count; ++ray)
    {
        stream << "ray " << geometry.gauges[ray];
        for (int coordinate = 0; coordinate < geometry.dimension; ++coordinate)
            stream << ' ' << geometry.rays[
                static_cast<size_t>(ray) * geometry.dimension + coordinate];
        stream << '\n';
    }
    return static_cast<bool>(stream);
}

bool load_cone_triangulation(const ConeRayGeometry &geometry, bool min_flag,
                             const char *directory,
                             std::vector<std::vector<int>> &children,
                             double *coverage_error,
                             double *maximum_condition)
{
    std::ifstream stream(cone_pilot_file(directory, min_flag, ".tri"));
    if (!stream) return false;
    std::string token;
    int version = 0, dimension = 0, ray_count = 0, child_count = 0;
    uint64_t signature = 0;
    int validated = 0;
    if (!(stream >> token >> version) || token != "DECIDEIT_TRIANGULATION" ||
        version != 1 || !(stream >> token >> dimension) ||
        token != "dimension" || !(stream >> token >> ray_count) ||
        token != "ray_count" || !(stream >> token >> signature) ||
        token != "signature" || !(stream >> token >> validated) ||
        token != "validated" || !(stream >> token >> *coverage_error) ||
        token != "coverage_error" ||
        !(stream >> token >> *maximum_condition) ||
        token != "maximum_condition" || !(stream >> token >> child_count) ||
        token != "children")
        throw std::runtime_error("invalid cone triangulation header");
    const int child_limit = cone_pilot_limit(
        "DECIDEIT_CONE_CHILD_MAX", 64, 10000);
    if (!validated || dimension != geometry.dimension ||
        ray_count != geometry.ray_count ||
        signature != cone_geometry_signature(geometry) || child_count < 1 ||
        child_count > child_limit || !std::isfinite(*coverage_error) ||
        *coverage_error < 0.0 || *coverage_error > 1e-7 ||
        !std::isfinite(*maximum_condition) ||
        *maximum_condition > 1e8)
        throw std::runtime_error(
            "cone triangulation failed validation, signature, or child limit");
    children.assign(child_count, std::vector<int>(dimension));
    std::vector<bool> covered_ray(ray_count, false);
    std::set<std::vector<int>> unique_children;
    for (int child = 0; child < child_count; ++child)
    {
        if (!(stream >> token) || token != "child")
            throw std::runtime_error("missing cone child row");
        std::vector<bool> used(ray_count, false);
        for (int j = 0; j < dimension; ++j)
        {
            int index = -1;
            if (!(stream >> index) || index < 0 || index >= ray_count ||
                used[index])
                throw std::runtime_error("invalid or duplicate child ray");
            used[index] = true;
            covered_ray[index] = true;
            children[child][j] = index;
        }
        std::vector<int> canonical = children[child];
        std::sort(canonical.begin(), canonical.end());
        if (!unique_children.insert(canonical).second)
            throw std::runtime_error("duplicate cone child");
    }
    if (std::find(covered_ray.begin(), covered_ray.end(), false) !=
        covered_ray.end())
        throw std::runtime_error("cone triangulation omitted a true ray");
    if (stream >> token)
        throw std::runtime_error("unexpected trailing cone triangulation data");
    return true;
}

void copy_vector_prefix(gsl_vector *destination, const gsl_vector *source,
                        size_t count)
{
    gsl_vector_const_view from = gsl_vector_const_subvector(source, 0, count);
    gsl_vector_view to = gsl_vector_subvector(destination, 0, count);
    gsl_vector_memcpy(&to.vector, &from.vector);
}

void copy_vector_prefix(gsl_vector_int *destination,
                        const gsl_vector_int *source, size_t count)
{
    gsl_vector_int_const_view from =
        gsl_vector_int_const_subvector(source, 0, count);
    gsl_vector_int_view to = gsl_vector_int_subvector(destination, 0, count);
    gsl_vector_int_memcpy(&to.vector, &from.vector);
}

void clone_rlt_for_cone_child(const LinearizationParams &source,
                              int added_cuts,
                              LinearizationParams *destination)
{
    const int variables = source.n_XY / 2;
    const size_t row_capacity = source.num_row +
        static_cast<size_t>(added_cuts) * (2 * variables + 2) + 2;
    const size_t nz_capacity = source.num_nz +
        static_cast<size_t>(added_cuts) *
            (4ULL * variables * variables + 3ULL * variables + 4ULL) + 2;
    destination->m_total = source.m_total;
    destination->n_total = source.n_total;
    destination->n_XY = source.n_XY;
    destination->n_XX = source.n_XX;
    destination->n_YY = source.n_YY;
    destination->num_col = source.num_col;
    destination->num_row = source.num_row;
    destination->num_nz = source.num_nz;
    destination->offset = source.offset;
    destination->a_start_column = source.a_start_column;
    destination->a_index_column = source.a_index_column;
    destination->col_cost = managed::vector_calloc(source.num_col);
    destination->col_lower = managed::vector_calloc(source.num_col);
    destination->col_upper = managed::vector_calloc(source.num_col);
    destination->row_lower = managed::vector_calloc(row_capacity);
    destination->row_upper = managed::vector_calloc(row_capacity);
    destination->a_start = managed::vector_int_calloc(row_capacity + 1);
    destination->a_value = managed::vector_calloc(nz_capacity);
    destination->a_index = managed::vector_int_calloc(nz_capacity);
    copy_vector_prefix(destination->col_cost, source.col_cost, source.num_col);
    copy_vector_prefix(destination->col_lower, source.col_lower, source.num_col);
    copy_vector_prefix(destination->col_upper, source.col_upper, source.num_col);
    copy_vector_prefix(destination->row_lower, source.row_lower, source.num_row);
    copy_vector_prefix(destination->row_upper, source.row_upper, source.num_row);
    copy_vector_prefix(destination->a_start, source.a_start,
                       source.num_row + 1);
    copy_vector_prefix(destination->a_value, source.a_value, source.num_nz);
    copy_vector_prefix(destination->a_index, source.a_index, source.num_nz);
}

enum class IsolatedRltStatus { optimal, infeasible, interrupted, failed };

IsolatedRltStatus solve_isolated_rlt(
    LinearizationParams *parameters, double time_limit, double *objective,
    long long *iterations, std::vector<double> *primal = nullptr,
    std::vector<double> *dual = nullptr,
    const HighsBasis *warm_basis = nullptr,
    HighsBasis *resulting_basis = nullptr, bool *warm_basis_used = nullptr,
    std::shared_ptr<Highs> *continuation = nullptr)
{
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    const bool resumed = continuation && static_cast<bool>(*continuation);
    std::shared_ptr<Highs> holder = resumed
        ? *continuation : std::make_shared<Highs>();
    Highs &solver = *holder;
    if (warm_basis_used) *warm_basis_used = false;
    if (!resumed) configure_rlt_highs(solver);
    highs_require(solver.setOptionValue("output_flag", false),
                  "isolated RLT output flag");
    // HiGHS interprets time_limit against the cumulative runtime of a model.
    // A continuation therefore needs an absolute deadline, not another copy
    // of the per-call allowance, or every resumed run stops immediately.
    const double absolute_time_limit = resumed
        ? solver.getRunTime() + time_limit : time_limit;
    highs_require(solver.setOptionValue("time_limit", absolute_time_limit),
                  "isolated RLT time limit");
    if (!resumed)
    {
        HighsLp lp;
        lp.num_col_ = parameters->num_col;
        lp.num_row_ = parameters->num_row;
        lp.col_cost_.assign(parameters->col_cost->data,
                            parameters->col_cost->data + parameters->num_col);
        lp.col_lower_.assign(parameters->col_lower->data,
                             parameters->col_lower->data + parameters->num_col);
        lp.col_upper_.assign(parameters->col_upper->data,
                             parameters->col_upper->data + parameters->num_col);
        lp.row_lower_.assign(parameters->row_lower->data,
                             parameters->row_lower->data + parameters->num_row);
        lp.row_upper_.assign(parameters->row_upper->data,
                             parameters->row_upper->data + parameters->num_row);
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.start_.assign(parameters->a_start->data,
            parameters->a_start->data + parameters->num_row + 1);
        lp.a_matrix_.index_.assign(parameters->a_index->data,
            parameters->a_index->data + parameters->num_nz);
        lp.a_matrix_.value_.assign(parameters->a_value->data,
            parameters->a_value->data + parameters->num_nz);
        lp.sense_ = ObjSense::kMinimize;
        if (solver.passModel(lp) == HighsStatus::kError)
            return IsolatedRltStatus::failed;
        configure_rlt_scale_for_columns(solver, parameters->num_col);
        if (warm_basis && warm_basis->valid &&
            warm_basis->col_status.size() ==
                static_cast<size_t>(parameters->num_col) &&
            warm_basis->row_status.size() ==
                static_cast<size_t>(parameters->num_row) &&
            solver.setBasis(*warm_basis) != HighsStatus::kError)
        {
            if (warm_basis_used) *warm_basis_used = true;
        }
    }
    else if (warm_basis_used) *warm_basis_used = true;
    if (solver.run() == HighsStatus::kError)
    {
        if (continuation) continuation->reset();
        return IsolatedRltStatus::failed;
    }
    *iterations = solver.getInfo().simplex_iteration_count;
    const HighsModelStatus status = solver.getModelStatus();
    if (status == HighsModelStatus::kInfeasible)
    {
        if (continuation) continuation->reset();
        return IsolatedRltStatus::infeasible;
    }
    if (status == HighsModelStatus::kTimeLimit ||
        status == HighsModelStatus::kIterationLimit)
    {
        if (resulting_basis) *resulting_basis = solver.getBasis();
        if (primal)
        {
            const HighsSolution &partial = solver.getSolution();
            if (partial.col_value.size() ==
                    static_cast<size_t>(parameters->num_col))
                *primal = partial.col_value;
        }
        if (continuation) *continuation = holder;
        return IsolatedRltStatus::interrupted;
    }
    if (status != HighsModelStatus::kOptimal || !highs_solution_valid(solver))
    {
        if (continuation) continuation->reset();
        return IsolatedRltStatus::failed;
    }
    *objective = solver.getObjectiveValue();
    const HighsSolution &solution = solver.getSolution();
    if (primal) *primal = solution.col_value;
    if (dual) *dual = solution.col_dual;
    if (resulting_basis) *resulting_basis = solver.getBasis();
    if (continuation) continuation->reset();
    return IsolatedRltStatus::optimal;
}

bool append_cone_row(const std::vector<double> &coefficients, double rhs,
                     gsl_vector *local_row, gsl_vector *full_row,
                     LinearizationParams *parameters,
                     const dblp *instance, gsl_vector_int *fixed_variables,
                     d_row P_lobo, d_row P_upbo,
                     int Ai_begin, int Aj_begin, int n_alt_i)
{
    double norm_squared = 0.0;
    for (double value : coefficients) norm_squared += value * value;
    const double norm = std::sqrt(norm_squared);
    if (!std::isfinite(norm) || norm <= 1e-12) return false;
    gsl_vector_set_zero(local_row);
    gsl_vector_set(local_row, 0, rhs / norm);
    for (size_t i = 0; i < coefficients.size(); ++i)
        gsl_vector_set(local_row, i + 1, coefficients[i] / norm);
    recover_cut(full_row, local_row, instance->xl, fixed_variables);
    linearization_onecut2go(full_row, parameters, P_lobo, P_upbo,
                            Ai_begin, Aj_begin, n_alt_i);
    return true;
}

ConePilotResult run_cone_rlt_pilot(
    const ConeRayGeometry &geometry, const LinearizationParams *base,
    const dblp *instance, gsl_vector_int *fixed_variables,
    d_row P_lobo, d_row P_upbo, int Ai_begin, int Aj_begin, int n_alt_i,
    double incumbent, bool min_flag)
{
    ConePilotResult result;
    const char *directory = std::getenv("DECIDEIT_CONE_PILOT_DIR");
    if (!directory || !geometry.valid) return result;
    const auto started = std::chrono::steady_clock::now();
    const bool exported = export_cone_geometry(geometry, min_flag, directory);
    const int delta = geometry.ray_count - geometry.dimension;
    const int delta_limit = cone_pilot_limit(
        "DECIDEIT_CONE_DELTA_MAX", 10, 1000);
    if (!exported || delta < 1 || delta > delta_limit || !base ||
        !base->col_cost)
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"cone_rlt_pilot\",\"status\":\"export_only\",\"sense\":\"%s\",\"dimension\":%d,\"rays\":%d,\"delta\":%d,\"exported\":%s,\"rlt_ready\":%s}\n",
                    min_flag ? "min" : "negated_max", geometry.dimension,
                    geometry.ray_count, delta, exported ? "true" : "false",
                    base && base->col_cost ? "true" : "false");
        return result;
    }
    std::vector<std::vector<int>> children;
    double coverage_error = 0.0, triangulation_condition = 0.0;
    if (!load_cone_triangulation(geometry, min_flag, directory, children,
                                 &coverage_error,
                                 &triangulation_condition))
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"cone_rlt_pilot\",\"status\":\"awaiting_triangulation\",\"sense\":\"%s\",\"dimension\":%d,\"rays\":%d,\"delta\":%d}\n",
                    min_flag ? "min" : "negated_max", geometry.dimension,
                    geometry.ray_count, delta);
        return result;
    }
    result.attempted = true;
    result.children = static_cast<int>(children.size());
    const double child_time_limit = cone_pilot_limit(
        "DECIDEIT_CONE_CHILD_TIME_LIMIT", 10, 3600);
    const double condition_limit = 1e8;
    managed::Scope scope;
    gsl_matrix *ray_matrix = managed::matrix_alloc(
        geometry.dimension, geometry.dimension);
    gsl_matrix *lu = managed::matrix_alloc(
        geometry.dimension, geometry.dimension);
    gsl_matrix *inverse = managed::matrix_alloc(
        geometry.dimension, geometry.dimension);
    gsl_permutation *permutation = managed::permutation_alloc(
        geometry.dimension);
    gsl_vector *local_row = managed::vector_calloc(geometry.dimension + 1);
    gsl_vector *full_row = managed::vector_calloc(base->n_XY / 2 + 1);
    for (const auto &child_indices : children)
    {
        for (int column = 0; column < geometry.dimension; ++column)
            for (int row = 0; row < geometry.dimension; ++row)
                gsl_matrix_set(ray_matrix, row, column,
                    geometry.rays[static_cast<size_t>(child_indices[column]) *
                                      geometry.dimension + row]);
        gsl_matrix_memcpy(lu, ray_matrix);
        int sign = 0;
        if (gsl_linalg_LU_decomp(lu, permutation, &sign) != GSL_SUCCESS ||
            std::fabs(gsl_linalg_LU_det(lu, sign)) <= 1e-14 ||
            gsl_linalg_LU_invert(lu, permutation, inverse) != GSL_SUCCESS)
            throw std::runtime_error("singular cone child basis");
        double ray_norm = 0.0, inverse_norm = 0.0;
        for (int row = 0; row < geometry.dimension; ++row)
        {
            double ray_sum = 0.0, inverse_sum = 0.0;
            for (int column = 0; column < geometry.dimension; ++column)
            {
                ray_sum += std::fabs(gsl_matrix_get(ray_matrix, row, column));
                inverse_sum += std::fabs(gsl_matrix_get(inverse, row, column));
            }
            ray_norm = (std::max)(ray_norm, ray_sum);
            inverse_norm = (std::max)(inverse_norm, inverse_sum);
        }
        if (!std::isfinite(ray_norm * inverse_norm) ||
            ray_norm * inverse_norm > condition_limit)
            throw std::runtime_error("ill-conditioned cone child basis");

        managed::Scope child_scope;
        LinearizationParams child{};
        clone_rlt_for_cone_child(*base, geometry.dimension + 1, &child);
        std::vector<double> coefficients(geometry.dimension, 0.0);
        for (int coordinate = 0; coordinate < geometry.dimension;
             ++coordinate)
        {
            double rhs = 0.0;
            for (int j = 0; j < geometry.dimension; ++j)
            {
                coefficients[j] = -gsl_matrix_get(inverse, coordinate, j);
                rhs += gsl_matrix_get(inverse, coordinate, j) *
                       geometry.apex[j];
            }
            if (!append_cone_row(coefficients, -rhs, local_row, full_row,
                                 &child, instance, fixed_variables,
                                 P_lobo, P_upbo, Ai_begin, Aj_begin, n_alt_i))
                throw std::runtime_error("invalid cone membership row");
        }
        std::fill(coefficients.begin(), coefficients.end(), 0.0);
        for (int coordinate = 0; coordinate < geometry.dimension;
             ++coordinate)
            for (int ray = 0; ray < geometry.dimension; ++ray)
                coefficients[coordinate] -=
                    geometry.gauges[child_indices[ray]] *
                    gsl_matrix_get(inverse, ray, coordinate);
        double polar_norm_squared = 0.0, polar_at_apex = 0.0;
        for (int coordinate = 0; coordinate < geometry.dimension;
             ++coordinate)
        {
            polar_norm_squared += coefficients[coordinate] *
                                  coefficients[coordinate];
            polar_at_apex -= coefficients[coordinate] *
                             geometry.apex[coordinate];
        }
        const double polar_norm = std::sqrt(polar_norm_squared);
        if (!std::isfinite(polar_norm) || polar_norm <= 1e-12 ||
            1.0 / polar_norm <= 8 * DBLP_EPSILON)
            throw std::runtime_error("numerically shallow child polar cut");
        const double polar_rhs = -1.0 - polar_at_apex;
        if (!append_cone_row(coefficients, polar_rhs, local_row, full_row,
                             &child, instance, fixed_variables,
                             P_lobo, P_upbo, Ai_begin, Aj_begin, n_alt_i))
            throw std::runtime_error("invalid child polar row");

        double bound = 0.0;
        long long iterations = 0;
        const IsolatedRltStatus status = solve_isolated_rlt(
            &child, child_time_limit, &bound, &iterations);
        result.simplex_iterations += iterations;
        if (status == IsolatedRltStatus::failed ||
            status == IsolatedRltStatus::interrupted)
            throw std::runtime_error("cone child RLT did not finish");
        if (status == IsolatedRltStatus::infeasible)
        {
            result.infeasible++;
            continue;
        }
        result.lower_bound = (std::min)(result.lower_bound, bound);
    }
    result.certified = result.infeasible == result.children ||
        result.lower_bound >= incumbent - DBLP_EPSILON;
    result.wall_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    if (std::getenv("DECIDEIT_TRACE"))
        fprintf(stderr, "TRACE {\"event\":\"cone_rlt_pilot\",\"status\":\"%s\",\"sense\":\"%s\",\"dimension\":%d,\"rays\":%d,\"delta\":%d,\"children\":%d,\"infeasible\":%d,\"lower_bound\":%.17g,\"incumbent\":%.17g,\"coverage_error\":%.12g,\"triangulation_condition\":%.12g,\"simplex_iterations\":%lld,\"wall_s\":%.9g}\n",
                result.certified ? "certified" : "open",
                min_flag ? "min" : "negated_max", geometry.dimension,
                geometry.ray_count, delta, result.children,
                result.infeasible, result.lower_bound, incumbent,
                coverage_error, triangulation_condition,
                result.simplex_iterations, result.wall_s);
    managed::matrix_free(ray_matrix);
    managed::matrix_free(lu);
    managed::matrix_free(inverse);
    managed::permutation_free(permutation);
    managed::vector_free(local_row);
    managed::vector_free(full_row);
    return result;
}

bool branch_rlt_tail_enabled()
{
    return binary_environment_option("DECIDEIT_BRANCH_RLT_TAIL");
}

struct BranchRltLimits
{
    int nodes = 128;
    int total_ms = 20000;
    int lp_ms = 8000;
    bool borrow_remaining = false;
};

BranchRltLimits selected_branch_rlt_limits(bool early)
{
    BranchRltLimits limits;
    if (early)
    {
        limits.nodes = cone_pilot_limit(
            "DECIDEIT_EARLY_BRANCH_RLT_NODE_LIMIT", 16, 4096);
        limits.total_ms = cone_pilot_limit(
            "DECIDEIT_EARLY_BRANCH_RLT_TIME_LIMIT_MS", 3000, 3600000);
        limits.lp_ms = cone_pilot_limit(
            "DECIDEIT_EARLY_BRANCH_RLT_LP_TIME_LIMIT_MS", 1500, 3600000);
        limits.borrow_remaining = binary_environment_option(
            "DECIDEIT_EARLY_BRANCH_RLT_BORROW_REMAINING");
    }
    else
    {
        limits.nodes = cone_pilot_limit(
            "DECIDEIT_BRANCH_RLT_NODE_LIMIT", 128, 4096);
        limits.total_ms = cone_pilot_limit(
            "DECIDEIT_BRANCH_RLT_TIME_LIMIT_MS", 20000, 3600000);
        limits.lp_ms = cone_pilot_limit(
            "DECIDEIT_BRANCH_RLT_LP_TIME_LIMIT_MS", 8000, 3600000);
        limits.borrow_remaining = binary_environment_option(
            "DECIDEIT_BRANCH_RLT_BORROW_REMAINING");
    }
    return limits;
}

struct BranchRltAnalysis
{
    bool feasible_marginals = false;
    bool rank_one_optimum = false;
    double feasible_objective = std::numeric_limits<double>::infinity();
    double maximum_residual = std::numeric_limits<double>::infinity();
    int residual_i = -1;
    int residual_j = -1;
};

BranchRltAnalysis analyze_branch_rlt_solution(
    const LinearizationParams &parameters,
    const std::vector<double> &primal, int n, double lower_bound)
{
    BranchRltAnalysis analysis;
    if (primal.size() != static_cast<size_t>(parameters.num_col) ||
        parameters.num_col != 2 * n + n * n)
        return analysis;

    std::vector<double> rank_one(parameters.num_col, 0.0);
    for (int column = 0; column < 2 * n; ++column)
    {
        if (!std::isfinite(primal[column])) return analysis;
        rank_one[column] = primal[column];
    }
    const int product_offset = 2 * n;
    analysis.maximum_residual = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            const int column = product_offset + i * n + j;
            rank_one[column] = primal[i] * primal[n + j];
            const double residual =
                std::fabs(primal[column] - rank_one[column]);
            if (residual > analysis.maximum_residual)
            {
                analysis.maximum_residual = residual;
                analysis.residual_i = i;
                analysis.residual_j = j;
            }
        }

    const double feasibility_tolerance = DBLP_EPSILON;
    for (int column = 0; column < parameters.num_col; ++column)
    {
        const double value = rank_one[column];
        if (!std::isfinite(value) ||
            value < gsl_vector_get(parameters.col_lower, column) -
                        feasibility_tolerance ||
            value > gsl_vector_get(parameters.col_upper, column) +
                        feasibility_tolerance)
            return analysis;
    }
    for (int row = 0; row < parameters.num_row; ++row)
    {
        double activity = 0.0;
        const int begin = gsl_vector_int_get(parameters.a_start, row);
        const int end = gsl_vector_int_get(parameters.a_start, row + 1);
        for (int entry = begin; entry < end; ++entry)
            activity += gsl_vector_get(parameters.a_value, entry) *
                rank_one[gsl_vector_int_get(parameters.a_index, entry)];
        if (!std::isfinite(activity) ||
            activity < gsl_vector_get(parameters.row_lower, row) -
                           feasibility_tolerance ||
            activity > gsl_vector_get(parameters.row_upper, row) +
                           feasibility_tolerance)
            return analysis;
    }
    analysis.feasible_marginals = true;
    analysis.feasible_objective = parameters.offset;
    for (int column = 0; column < parameters.num_col; ++column)
        analysis.feasible_objective +=
            gsl_vector_get(parameters.col_cost, column) * rank_one[column];
    analysis.rank_one_optimum =
        analysis.maximum_residual <= 10.0 * DBLP_EPSILON &&
        std::fabs(analysis.feasible_objective - lower_bound) <=
            DBLP_EPSILON;
    return analysis;
}

struct BranchRltNode
{
    std::vector<double> lower;
    std::vector<double> upper;
    std::vector<double> primal;
    std::vector<double> dual;
    std::vector<std::vector<int> > cut_cross_rows;
    HighsBasis basis;
    std::shared_ptr<Highs> continuation;
    BranchRltAnalysis analysis;
    double bound = -std::numeric_limits<double>::infinity();
    int depth = 0;
    uint64_t serial = 0;
};

struct BranchRltNodeOrder
{
    bool operator()(const BranchRltNode &left,
                    const BranchRltNode &right) const
    {
        if (left.bound != right.bound) return left.bound > right.bound;
        if (left.depth != right.depth) return left.depth < right.depth;
        return left.serial > right.serial;
    }
};

struct BranchRltTailResult
{
    bool attempted = false;
    bool certified = false;
    bool budget_exhausted = false;
    bool lp_failed = false;
    bool initial_incumbent_valid = false;
    bool final_incumbent_valid = false;
    int nodes = 0;
    int lp_solves = 0;
    int infeasible = 0;
    int bound_pruned = 0;
    int exact = 0;
    int incumbent_updates = 0;
    int polish_attempts = 0;
    int polish_lp_solves = 0;
    int polish_updates = 0;
    int resumed_lp_solves = 0;
    int warm_basis_uses = 0;
    int reduced_cost_updates = 0;
    int maximum_depth = 0;
    int open_nodes = 0;
    long long simplex_iterations = 0;
    double best_open_bound = std::numeric_limits<double>::infinity();
    double maximum_residual = 0.0;
    double reduced_cost_width_reduction = 0.0;
    double node_model_build_s = 0.0;
    double node_lp_s = 0.0;
    double polish_s = 0.0;
    double wall_s = 0.0;
};

void tighten_branch_rlt_node_from_reduced_costs(
    BranchRltNode *node, int n, double incumbent, int candidate_limit,
    const std::string &rank_mode, BranchRltTailResult *result)
{
    const size_t columns = static_cast<size_t>(2 * n + n * n);
    if (!node || !result || !std::isfinite(incumbent) ||
        !std::isfinite(node->bound) || node->primal.size() != columns ||
        node->dual.size() != columns || !node->basis.valid ||
        node->basis.col_status.size() != columns)
        return;

    std::vector<double> x_residual(n, 0.0), y_residual(n, 0.0);
    double maximum_residual = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            const int product = 2 * n + i * n + j;
            const double residual = std::fabs(
                node->primal[product] -
                node->primal[i] * node->primal[n + j]);
            x_residual[i] = (std::max)(x_residual[i], residual);
            y_residual[j] = (std::max)(y_residual[j], residual);
            maximum_residual = (std::max)(maximum_residual, residual);
        }
    double x_width_sum = 0.0, y_width_sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        x_width_sum += (std::max)(0.0, node->upper[i] - node->lower[i]);
        y_width_sum += (std::max)(
            0.0, node->upper[n + i] - node->lower[n + i]);
    }
    const double cutoff = incumbent - DBLP_EPSILON;
    const double objective_guard = 1e-7 *
        (1.0 + std::fabs(cutoff) + std::fabs(node->bound));
    const double allowance = (std::max)(0.0, cutoff - node->bound) +
        objective_guard;
    const double dual_tolerance = 1e-8;
    struct Candidate
    {
        int column;
        double bound;
        double reduction;
        double score;
        bool upper;
    };
    std::vector<Candidate> candidates;
    candidates.reserve(2 * n);
    for (int column = 0; column < 2 * n; ++column)
    {
        const double lower = node->lower[column];
        const double upper = node->upper[column];
        if (upper - lower <= 1e-10) continue;
        const double primal = node->primal[column];
        const double dual = node->dual[column];
        const HighsBasisStatus status = node->basis.col_status[column];
        double bound = 0.0, reduction = 0.0;
        bool tighten_upper = false;
        if (status == HighsBasisStatus::kLower && dual > dual_tolerance &&
            std::fabs(primal - lower) <=
                1e-7 * (1.0 + std::fabs(lower)))
        {
            bound = lower + allowance / dual;
            bound += 1e-9 * (1.0 + std::fabs(bound));
            bound = (std::max)(bound, lower);
            reduction = upper - bound;
            tighten_upper = true;
        }
        else if (status == HighsBasisStatus::kUpper &&
                 dual < -dual_tolerance &&
                 std::fabs(primal - upper) <=
                    1e-7 * (1.0 + std::fabs(upper)))
        {
            bound = upper - allowance / (-dual);
            bound -= 1e-9 * (1.0 + std::fabs(bound));
            bound = (std::min)(bound, upper);
            reduction = bound - lower;
        }
        if (reduction <= 1e-10) continue;
        const int local = column < n ? column : column - n;
        const double residual = column < n
            ? x_residual[local] : y_residual[local];
        const double envelope_scale = column < n
            ? y_width_sum : x_width_sum;
        double score = reduction * envelope_scale;
        if (rank_mode == "residual") score = reduction * residual;
        else if (rank_mode == "hybrid")
        {
            const double residual_weight = maximum_residual > 0.0
                ? residual / maximum_residual : 0.0;
            score *= 0.25 + 0.75 * residual_weight;
        }
        candidates.push_back(Candidate{
            column, bound, reduction, score, tighten_upper});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate &left, const Candidate &right) {
                  return left.score != right.score
                      ? left.score > right.score
                      : left.reduction > right.reduction;
              });
    if (candidate_limit < static_cast<int>(candidates.size()))
        candidates.resize(candidate_limit);
    for (const Candidate &candidate : candidates)
    {
        if (candidate.upper)
            node->upper[candidate.column] = candidate.bound;
        else
            node->lower[candidate.column] = candidate.bound;
        result->reduced_cost_updates++;
        result->reduced_cost_width_reduction += candidate.reduction;
    }
}

bool original_incumbent_feasible(
    const global_solution &solution, const dblp &instance,
    gsl_vector_int *fixed_variables,
    const std::vector<LinCons> &original_y_constraints,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Aj_begin, int n_alt_i, int n_alt_j,
    bool min_flag, double *recomputed_objective)
{
    const int variables = n_alt_i + n_alt_j;
    gsl_vector *whole_x = managed::vector_calloc(variables);
    getfrom_x_solution(fixed_variables, whole_x, instance.xl,
                       solution.x_solution);
    std::vector<double> x(whole_x->data, whole_x->data + variables);
    std::vector<double> y(solution.y_solution->data,
                          solution.y_solution->data + variables);
    managed::vector_free(whole_x);
    // This gate controls a global-pruning threshold, so it is deliberately
    // tighter than the ordinary LP feasibility tolerance.  A merely
    // near-feasible AMC point may remain a useful warm start, but must not be
    // allowed to certify a branch tree.
    const double tolerance = 1e-9;
    auto within = [&](double value, double lower, double upper) {
        const double scale = 1.0 + std::fabs(value) +
            (std::isfinite(lower) ? std::fabs(lower) : 0.0) +
            (std::isfinite(upper) ? std::fabs(upper) : 0.0);
        return value >= lower - tolerance * scale &&
               value <= upper + tolerance * scale;
    };
    for (int column = 0; column < variables; ++column)
    {
        const int original = location(
            column, Ai_begin, Aj_begin, n_alt_i);
        if (!std::isfinite(x[column]) || !std::isfinite(y[column]) ||
            !within(x[column], V_lobo[original], V_upbo[original]) ||
            !within(y[column], P_lobo[original], P_upbo[original]))
            return false;
    }
    auto rows_feasible = [&](const std::vector<LinCons> &rows,
                             const std::vector<double> &point) {
        for (const LinCons &row : rows)
        {
            double activity = 0.0;
            double scale = 1.0 + std::fabs(row.b);
            for (size_t k = 0; k < row.index.size(); ++k)
            {
                const int column = row.index[k] - 1;
                if (column < 0 || column >= variables) return false;
                activity += row.coefficient[k] * point[column];
                scale += std::fabs(row.coefficient[k] * point[column]);
            }
            if (!std::isfinite(activity) ||
                activity > row.b + tolerance * scale)
                return false;
        }
        return true;
    };
    // Validate X against every original row, not only the reduced set used by
    // the fast solver.  A row removed as numerically redundant may still be
    // the one that distinguishes a near-feasible incumbent from a safe global
    // pruning threshold.
    for (size_t row = 0; row < instance.Ax->size1; ++row)
    {
        double activity = 0.0;
        double scale = 1.0 + std::fabs(gsl_vector_get(instance.bx, row));
        for (int column = 0; column < variables; ++column)
        {
            const double term = gsl_matrix_get(instance.Ax, row, column) *
                x[column];
            activity += term;
            scale += std::fabs(term);
        }
        if (!std::isfinite(activity) ||
            activity > gsl_vector_get(instance.bx, row) + tolerance * scale)
            return false;
    }
    if (!rows_feasible(original_y_constraints, y))
        return false;
    double first_sum = 0.0, second_sum = 0.0;
    for (int column = 0; column < variables; ++column)
        if (column < n_alt_i) first_sum += y[column];
        else second_sum += y[column];
    if (std::fabs(first_sum - 1.0) > tolerance * (1.0 + n_alt_i) ||
        std::fabs(second_sum - 1.0) > tolerance * (1.0 + n_alt_j))
        return false;
    double objective = 0.0;
    for (int column = 0; column < variables; ++column)
        objective += (column < n_alt_i ? 1.0 : -1.0) *
            x[column] * y[column];
    if (!min_flag) objective = -objective;
    if (recomputed_objective) *recomputed_objective = objective;
    return std::isfinite(objective) &&
        std::fabs(objective - solution.alpha) <=
            tolerance * (1.0 + std::fabs(objective));
}

IsolatedRltStatus evaluate_branch_rlt_node(
    BranchRltNode *node, const dblp *instance, const dataset *dataset_in,
    gsl_vector_int *fixed_variables,
    std::vector<LinCons> &linear_cons_P,
    std::vector<LinCons> &linear_cons_V,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    int n_alt_i, int n_alt_j, int mX, int mY, int noc,
    bool min_flag, double time_limit,
    const std::vector<double> *separation_primal,
    const std::vector<std::vector<int> > *forced_cross_plan,
    const HighsBasis *parent_basis, bool capture_reduced_costs,
    BranchRltTailResult *result)
{
    const auto node_started = std::chrono::steady_clock::now();
    const int n = n_alt_i + n_alt_j;
    if (branch_rlt_propagation_enabled() &&
        !propagate_branch_rlt_bounds(
            &node->lower, &node->upper, n, n_alt_i,
            linear_cons_V, linear_cons_P,
            dataset_in, noc, fixed_variables))
        return IsolatedRltStatus::infeasible;
    d_row node_P_lobo, node_P_upbo, node_V_lobo, node_V_upbo;
    for (int i = 0; i <= VX_MAXDIM; ++i)
    {
        node_P_lobo[i] = P_lobo[i]; node_P_upbo[i] = P_upbo[i];
        node_V_lobo[i] = V_lobo[i]; node_V_upbo[i] = V_upbo[i];
    }
    for (int i = 0; i < n; ++i)
    {
        const size_t original = location(i, Ai_begin, Aj_begin, n_alt_i);
        node_V_lobo[original] = node->lower[i];
        node_V_upbo[original] = node->upper[i];
        node_P_lobo[original] = node->lower[n + i];
        node_P_upbo[original] = node->upper[n + i];
    }

    LinearizationParams parameters;
    construct_LinearizationParameters(mX, mY, n, n, &parameters);
    construct_linearization_objective(n_alt_i, n_alt_j, &parameters);
    if (!min_flag) gsl_vector_scale(parameters.col_cost, -1.0);
    build_linearization(n_alt_i, n_alt_j, linear_cons_P, linear_cons_V,
                        node_P_lobo, node_P_upbo,
                        node_V_lobo, node_V_upbo,
                        Ai_begin, Ai_end, Aj_begin, Aj_end,
                        mX, mY, &parameters);
    gsl_vector *onecut = managed::vector_calloc(n + 1);
    int cut_index = 0;
    for (int row = dataset_in->ori_num_cons; row < noc;
         ++row, ++cut_index)
    {
        gsl_vector_view cut = gsl_matrix_row(dataset_in->valid_A, row);
        recover_cut(onecut, &cut.vector, instance->xl, fixed_variables);
        linearization_onecut2go(onecut, &parameters,
                                node_P_lobo, node_P_upbo,
                                Ai_begin, Aj_begin, n_alt_i,
                                &linear_cons_P,
                                separation_primal
                                    ? separation_primal->data() : nullptr,
                                separation_primal
                                    ? separation_primal->size() : 0,
                                forced_cross_plan && cut_index <
                                        static_cast<int>(forced_cross_plan->size())
                                    ? &(*forced_cross_plan)[cut_index]
                                    : nullptr);
    }
    if (forced_cross_plan && cut_index !=
            static_cast<int>(forced_cross_plan->size()))
        throw std::runtime_error("Branch RLT cross-row plan size mismatch");
    bool warm_used = false;
    long long iterations = 0;
    const double setup_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - node_started).count();
    result->node_model_build_s += setup_s;
    const double remaining_s = time_limit - setup_s;
    if (remaining_s <= 0.0)
    {
        managed::vector_free(onecut);
        return IsolatedRltStatus::failed;
    }
    const auto solve_started = std::chrono::steady_clock::now();
    const IsolatedRltStatus status = solve_isolated_rlt(
        &parameters, remaining_s, &node->bound, &iterations,
        &node->primal, capture_reduced_costs ? &node->dual : nullptr,
        parent_basis, &node->basis, &warm_used,
        &node->continuation);
    result->node_lp_s += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solve_started).count();
    result->simplex_iterations += iterations;
    if (warm_used) result->warm_basis_uses++;
    if (status == IsolatedRltStatus::optimal)
    {
        node->analysis = analyze_branch_rlt_solution(
            parameters, node->primal, n, node->bound);
    }
    if (status == IsolatedRltStatus::optimal ||
        status == IsolatedRltStatus::interrupted)
        node->cut_cross_rows = parameters.cut_cross_rows;
    managed::vector_free(onecut);
    return status;
}

void install_branch_rlt_incumbent(const BranchRltNode &node, int n,
                                  gsl_vector_int *fixed_variables,
                                  global_solution *solution)
{
    gsl_vector_view first = gsl_vector_view_array(
        const_cast<double *>(node.primal.data()), n);
    putto_global_x_solution(fixed_variables, &first.vector,
                            solution->x_solution);
    gsl_vector_view second = gsl_vector_view_array(
        const_cast<double *>(node.primal.data()) + n, n);
    gsl_vector_memcpy(solution->y_solution, &second.vector);
    solution->alpha = node.analysis.feasible_objective;
}

bool try_install_branch_rlt_incumbent(
    const BranchRltNode &node, int n, const dblp &instance,
    gsl_vector_int *fixed_variables,
    const std::vector<LinCons> &original_y_constraints,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Aj_begin, int n_alt_i, int n_alt_j,
    bool min_flag, global_solution *solution)
{
    gsl_vector *saved_x = managed::vector_alloc(solution->x_solution->size);
    gsl_vector *saved_y = managed::vector_alloc(solution->y_solution->size);
    gsl_vector_memcpy(saved_x, solution->x_solution);
    gsl_vector_memcpy(saved_y, solution->y_solution);
    const double saved_alpha = solution->alpha;

    install_branch_rlt_incumbent(node, n, fixed_variables, solution);
    double recomputed_objective = solution->alpha;
    const bool valid = original_incumbent_feasible(
        *solution, instance, fixed_variables, original_y_constraints,
        P_lobo, P_upbo, V_lobo, V_upbo,
        Ai_begin, Aj_begin, n_alt_i, n_alt_j,
        min_flag, &recomputed_objective);
    if (valid)
        solution->alpha = recomputed_objective;
    else
    {
        gsl_vector_memcpy(solution->x_solution, saved_x);
        gsl_vector_memcpy(solution->y_solution, saved_y);
        solution->alpha = saved_alpha;
    }
    managed::vector_free(saved_x);
    managed::vector_free(saved_y);
    return valid;
}

LpModel *clone_local_lp(const LpModel &source)
{
    LpModel *copy = make_lp(get_Nrows(&source), source.columns);
    copy->sense = source.sense;
    copy->rows = source.rows;
    copy->cost = source.cost;
    copy->lower = source.lower;
    copy->upper = source.upper;
    copy->rhs = source.rhs;
    copy->types = source.types;
    copy->backend = source.backend;
    copy->loaded = false;
    copy->needs_solve = true;
    return copy;
}

bool polish_branch_rlt_marginals(
    const BranchRltNode &node, int n, int maximum_rounds,
    LpModel *lp_X_std, LpModel *lp_Y, const dataset *dataset_in,
    const dblp &instance, gsl_vector_int *fixed_variables,
    const std::vector<LinCons> &original_y_constraints,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Aj_begin, int n_alt_i, int n_alt_j,
    bool min_flag, bool *incumbent_valid, global_solution *solution,
    BranchRltTailResult *result)
{
    if (!lp_X_std || !lp_Y || !dataset_in || maximum_rounds <= 0 ||
        node.primal.size() < static_cast<size_t>(2 * n))
        return false;

    const auto started = std::chrono::steady_clock::now();
    result->polish_attempts++;
    const int nx = get_Ncolumns(lp_X_std);
    const int ny = get_Ncolumns(lp_Y);
    if (ny != n || dataset_in->Q->size1 != static_cast<size_t>(nx) ||
        dataset_in->Q->size2 != static_cast<size_t>(ny))
        return false;

    gsl_vector *objective_x = managed::vector_calloc(nx + 1);
    gsl_vector *objective_y = managed::vector_calloc(ny + 1);
    gsl_vector *x = managed::vector_calloc(nx);
    gsl_vector *y = managed::vector_alloc(ny);
    gsl_vector *previous_y = managed::vector_alloc(ny);
    gsl_vector_view objective_x_body = gsl_vector_subvector(
        objective_x, 1, nx);
    gsl_vector_view objective_y_body = gsl_vector_subvector(
        objective_y, 1, ny);
    for (int column = 0; column < ny; ++column)
        gsl_vector_set(y, column, node.primal[n + column]);

    // Polish private copies inside the current branch box.  X is represented
    // by z=x-dataset.xl in the small LP, whereas Y uses original coordinates.
    // Private copies make a failed tail observationally harmless when ODBA
    // resumes with its existing bases and incumbent path.
    LpModel *polish_x = clone_local_lp(*lp_X_std);
    LpModel *polish_y = clone_local_lp(*lp_Y);
    int reduced_x = 0;
    for (int full = 0; full < n; ++full)
    {
        if (gsl_vector_int_get(fixed_variables, full + 1) != 0) continue;
        set_lowbo(polish_x, reduced_x + 1,
                  node.lower[full] - gsl_vector_get(dataset_in->xl,
                                                     reduced_x));
        set_upbo(polish_x, reduced_x + 1,
                 node.upper[full] - gsl_vector_get(dataset_in->xl,
                                                    reduced_x));
        reduced_x++;
    }
    for (int column = 0; column < ny; ++column)
    {
        set_lowbo(polish_y, column + 1, node.lower[n + column]);
        set_upbo(polish_y, column + 1, node.upper[n + column]);
    }

    std::vector<double> best_x;
    std::vector<double> best_y;
    double best_objective = std::numeric_limits<double>::infinity();
    for (int round = 0; round < maximum_rounds; ++round)
    {
        gsl_vector_memcpy(previous_y, y);
        gsl_blas_dgemv(CblasNoTrans, 1.0, dataset_in->Q, y, 0.0,
                       &objective_x_body.vector);
        set_obj_fn(polish_x, objective_x->data);
        result->polish_lp_solves++;
        if (solve(polish_x) != OPTIMAL) break;
        get_variables(polish_x, x->data);
        gsl_vector_add(x, dataset_in->xl);

        gsl_blas_dgemv(CblasTrans, 1.0, dataset_in->Q, x, 0.0,
                       &objective_y_body.vector);
        gsl_vector_add(&objective_y_body.vector, dataset_in->dy);
        set_obj_fn(polish_y, objective_y->data);
        result->polish_lp_solves++;
        if (solve(polish_y) != OPTIMAL) break;
        get_variables(polish_y, y->data);
        const double objective = get_objective(polish_y);
        if (std::isfinite(objective) && objective < best_objective)
        {
            best_objective = objective;
            best_x.assign(x->data, x->data + nx);
            best_y.assign(y->data, y->data + ny);
        }
        gsl_vector_sub(previous_y, y);
        if (gsl_blas_dnrm2(previous_y) <= DBLP_EPSILON) break;
    }

    bool improved = false;
    if (!best_x.empty() && (!*incumbent_valid ||
            best_objective < solution->alpha - DBLP_EPSILON))
    {
        gsl_vector *saved_x = managed::vector_alloc(
            solution->x_solution->size);
        gsl_vector *saved_y = managed::vector_alloc(
            solution->y_solution->size);
        gsl_vector_memcpy(saved_x, solution->x_solution);
        gsl_vector_memcpy(saved_y, solution->y_solution);
        const double saved_alpha = solution->alpha;
        for (int column = 0; column < nx; ++column)
            gsl_vector_set(solution->x_solution, column, best_x[column]);
        for (int column = 0; column < ny; ++column)
            gsl_vector_set(solution->y_solution, column, best_y[column]);
        solution->alpha = best_objective;
        double recomputed = best_objective;
        const bool valid = original_incumbent_feasible(
            *solution, instance, fixed_variables, original_y_constraints,
            P_lobo, P_upbo, V_lobo, V_upbo,
            Ai_begin, Aj_begin, n_alt_i, n_alt_j,
            min_flag, &recomputed);
        if (valid && (!*incumbent_valid ||
                recomputed < saved_alpha - DBLP_EPSILON))
        {
            solution->alpha = recomputed;
            *incumbent_valid = true;
            result->polish_updates++;
            improved = true;
        }
        else
        {
            gsl_vector_memcpy(solution->x_solution, saved_x);
            gsl_vector_memcpy(solution->y_solution, saved_y);
            solution->alpha = saved_alpha;
        }
        managed::vector_free(saved_x);
        managed::vector_free(saved_y);
    }

    managed::vector_free(objective_x);
    managed::vector_free(objective_y);
    managed::vector_free(x);
    managed::vector_free(y);
    managed::vector_free(previous_y);
    free_lp(&polish_x);
    free_lp(&polish_y);
    result->polish_s += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return improved;
}

BranchRltTailResult run_branch_rlt_tail(
    const LinearizationParams *root_parameters,
    const gsl_vector *root_solution, double root_bound,
    const dblp *instance, const dataset *dataset_in,
    gsl_vector_int *fixed_variables,
    std::vector<LinCons> &linear_cons_P,
    std::vector<LinCons> &linear_cons_V,
    const std::vector<LinCons> &original_P_rows,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    d_row original_P_lobo, d_row original_P_upbo,
    d_row original_V_lobo, d_row original_V_upbo,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    int n_alt_i, int n_alt_j, int mX, int mY, int noc,
    bool min_flag, global_solution *solution,
    const BranchRltLimits &limits,
    LpModel *lp_X_std, LpModel *lp_Y)
{
    BranchRltTailResult result;
    result.attempted = true;
    const int node_limit = limits.nodes;
    const int total_limit_ms = limits.total_ms;
    const int lp_limit_ms = limits.lp_ms;
    const double minimum_width = 1e-8;
    const int n = n_alt_i + n_alt_j;
    const bool polish_marginals_requested = binary_environment_option(
        "DECIDEIT_BRANCH_RLT_POLISH");
    const int polish_max_columns = polish_marginals_requested
        ? cone_pilot_limit(
            "DECIDEIT_BRANCH_RLT_POLISH_MAX_COLUMNS", 2500, 4000000)
        : 0;
    // Alternating marginal polishing is intentionally restricted to modest
    // RLTs.  On larger relaxations a better incumbent can change the later
    // cutting-plane path without being strong enough to certify the node,
    // making the complete solve slower.  The limit is tunable for ablations.
    const bool polish_marginals = polish_marginals_requested &&
        root_parameters && root_parameters->num_col <= polish_max_columns;
    const bool resume_interrupted = binary_environment_option(
        "DECIDEIT_BRANCH_RLT_RESUME_INTERRUPTED");
    const bool reduced_cost_tightening = binary_environment_option(
        "DECIDEIT_BRANCH_RLT_REDUCED_COST_TIGHTEN");
    const int reduced_cost_candidate_limit = reduced_cost_tightening
        ? cone_pilot_limit(
            "DECIDEIT_BRANCH_RLT_REDUCED_COST_VARIABLES", 4, 2 * VX_MAXDIM)
        : 0;
    const char *reduced_cost_rank_value = std::getenv(
        "DECIDEIT_BRANCH_RLT_REDUCED_COST_RANK");
    const std::string reduced_cost_rank = reduced_cost_rank_value
        ? reduced_cost_rank_value : "residual";
    if (reduced_cost_rank != "width" &&
        reduced_cost_rank != "residual" &&
        reduced_cost_rank != "hybrid")
        throw std::invalid_argument(
            "DECIDEIT_BRANCH_RLT_REDUCED_COST_RANK must be width, "
            "residual, or hybrid");
    const int polish_rounds = polish_marginals
        ? cone_pilot_limit("DECIDEIT_BRANCH_RLT_POLISH_ROUNDS", 2, 20)
        : 0;
    const int polish_node_limit = polish_marginals
        ? cone_pilot_limit("DECIDEIT_BRANCH_RLT_POLISH_NODE_LIMIT", 4, 4096)
        : 0;
    const char *branching_value = std::getenv(
        "DECIDEIT_BRANCH_RLT_BRANCHING");
    const std::string branching = branching_value
        ? branching_value : "residual";
    if (branching != "residual" && branching != "objective")
        throw std::invalid_argument(
            "DECIDEIT_BRANCH_RLT_BRANCHING must be residual or objective");
    double recomputed_incumbent = solution->alpha;
    bool incumbent_valid = original_incumbent_feasible(
        *solution, *instance, fixed_variables,
        original_P_rows,
        original_P_lobo, original_P_upbo,
        original_V_lobo, original_V_upbo,
        Ai_begin, Aj_begin, n_alt_i, n_alt_j,
        min_flag, &recomputed_incumbent);
    result.initial_incumbent_valid = incumbent_valid;
    if (incumbent_valid) solution->alpha = recomputed_incumbent;
    const auto started = std::chrono::steady_clock::now();
    const auto elapsed_ms = [&]() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - started).count();
    };

    BranchRltNode root;
    root.lower.resize(2 * n);
    root.upper.resize(2 * n);
    if (root_parameters && root_solution)
    {
        root.primal.assign(root_solution->data,
                           root_solution->data + root_parameters->num_col);
        root.bound = root_bound;
        for (int column = 0; column < 2 * n; ++column)
        {
            root.lower[column] = gsl_vector_get(
                root_parameters->col_lower, column);
            root.upper[column] = gsl_vector_get(
                root_parameters->col_upper, column);
        }
        root.analysis = analyze_branch_rlt_solution(
            *root_parameters, root.primal, n, root.bound);
        if (root_parameters->persistent)
        {
            root.basis = root_parameters->persistent->highs.getBasis();
            const HighsSolver &persistent = *root_parameters->persistent;
            if (reduced_cost_tightening &&
                persistent.reduced_cost_snapshot_valid &&
                persistent.reduced_cost_dual.size() ==
                    static_cast<size_t>(root_parameters->num_col))
                root.dual = persistent.reduced_cost_dual;
        }
        root.cut_cross_rows = root_parameters->cut_cross_rows;
        result.nodes = 1;
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            const size_t original = location(
                i, Ai_begin, Aj_begin, n_alt_i);
            root.lower[i] = V_lobo[original];
            root.upper[i] = V_upbo[original];
            root.lower[n + i] = P_lobo[original];
            root.upper[n + i] = P_upbo[original];
        }
        const double root_seconds = (limits.borrow_remaining
            ? total_limit_ms
            : (std::min)(lp_limit_ms, total_limit_ms)) / 1000.0;
        const IsolatedRltStatus root_status = evaluate_branch_rlt_node(
            &root, instance, dataset_in, fixed_variables,
            linear_cons_P, linear_cons_V,
            P_lobo, P_upbo, V_lobo, V_upbo,
            Ai_begin, Ai_end, Aj_begin, Aj_end,
            n_alt_i, n_alt_j, mX, mY, noc, min_flag,
            root_seconds, nullptr, nullptr, nullptr,
            reduced_cost_tightening, &result);
        result.nodes = 1;
        result.lp_solves = 1;
        if (root_status == IsolatedRltStatus::infeasible)
        {
            result.infeasible = 1;
            result.certified = incumbent_valid;
            result.final_incumbent_valid = incumbent_valid;
            result.wall_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            return result;
        }
        if (root_status == IsolatedRltStatus::failed)
        {
            result.lp_failed = true;
            result.budget_exhausted = true;
            result.open_nodes = 1;
            result.wall_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            return result;
        }
        if (root_status == IsolatedRltStatus::interrupted)
        {
            result.budget_exhausted = true;
            result.open_nodes = 1;
            result.wall_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            return result;
        }
    }
    root.serial = 1;

    std::priority_queue<BranchRltNode, std::vector<BranchRltNode>,
                        BranchRltNodeOrder> open;
    open.push(std::move(root));
    uint64_t serial = 1;
    int detached_open_nodes = 0;
    double detached_open_bound = std::numeric_limits<double>::infinity();
    while (!open.empty())
    {
        if (elapsed_ms() >= total_limit_ms)
        {
            // Solved nodes already carry rigorous lower bounds.  A time
            // budget must prevent new LP work, but it must not hide a free
            // certificate that is already sitting in the queue.  Since this
            // priority queue exposes the smallest bound first, one comparison
            // can prune the entire remaining queue.
            if (incumbent_valid &&
                open.top().bound >= solution->alpha - DBLP_EPSILON)
            {
                result.bound_pruned += static_cast<int>(open.size());
                while (!open.empty()) open.pop();
            }
            else
                result.budget_exhausted = true;
            break;
        }
        BranchRltNode node = open.top();
        open.pop();
        result.maximum_depth = (std::max)(result.maximum_depth, node.depth);
        result.maximum_residual = (std::max)(
            result.maximum_residual, node.analysis.maximum_residual);

        if (polish_marginals &&
            result.polish_attempts < polish_node_limit &&
            (!incumbent_valid ||
             node.bound < solution->alpha - DBLP_EPSILON))
            polish_branch_rlt_marginals(
                node, n, polish_rounds, lp_X_std, lp_Y, dataset_in,
                *instance, fixed_variables, original_P_rows,
                original_P_lobo, original_P_upbo,
                original_V_lobo, original_V_upbo,
                Ai_begin, Aj_begin, n_alt_i, n_alt_j,
                min_flag, &incumbent_valid, solution, &result);

        if (node.analysis.feasible_marginals &&
            (!incumbent_valid || node.analysis.feasible_objective <
                solution->alpha - DBLP_EPSILON))
        {
            if (try_install_branch_rlt_incumbent(
                    node, n, *instance, fixed_variables, original_P_rows,
                    original_P_lobo, original_P_upbo,
                    original_V_lobo, original_V_upbo,
                    Ai_begin, Aj_begin, n_alt_i, n_alt_j,
                    min_flag, solution))
            {
                incumbent_valid = true;
                result.incumbent_updates++;
            }
        }
        if (incumbent_valid &&
            node.bound >= solution->alpha - DBLP_EPSILON)
        {
            result.bound_pruned++;
            continue;
        }
        if (node.analysis.rank_one_optimum)
        {
            result.exact++;
            continue;
        }
        if (reduced_cost_tightening && incumbent_valid)
            tighten_branch_rlt_node_from_reduced_costs(
                &node, n, solution->alpha, reduced_cost_candidate_limit,
                reduced_cost_rank, &result);

        int i = node.analysis.residual_i;
        int j = node.analysis.residual_j;
        if (branching == "objective")
        {
            double best_score = -1.0;
            for (int candidate_i = 0; candidate_i < n; ++candidate_i)
                for (int candidate_j = 0; candidate_j < n; ++candidate_j)
                {
                    const int product = 2 * n + candidate_i * n +
                        candidate_j;
                    const double residual = std::fabs(
                        node.primal[product] -
                        node.primal[candidate_i] *
                            node.primal[n + candidate_j]);
                    const double first_width = (std::max)(0.0,
                        node.upper[candidate_i] -
                        node.lower[candidate_i]);
                    const double second_width = (std::max)(0.0,
                        node.upper[n + candidate_j] -
                        node.lower[n + candidate_j]);
                    // Only diagonal W columns appear directly in the DBLP
                    // objective.  Retain a positive score for off-diagonal
                    // residuals because constraints can still couple them.
                    const double objective_weight =
                        candidate_i == candidate_j ? 5.0 : 1.0;
                    const double score = residual * objective_weight *
                        std::sqrt(first_width * second_width);
                    if (score > best_score)
                    {
                        best_score = score;
                        i = candidate_i;
                        j = candidate_j;
                    }
                }
        }
        if (i < 0 || j < 0)
        {
            result.lp_failed = true;
            detached_open_nodes++;
            detached_open_bound = (std::min)(detached_open_bound, node.bound);
            break;
        }
        const int first_column = i;
        const int second_column = n + j;
        const double first_width =
            node.upper[first_column] - node.lower[first_column];
        const double second_width =
            node.upper[second_column] - node.lower[second_column];
        int branch_column = first_width >= second_width
            ? first_column : second_column;
        if (node.upper[branch_column] - node.lower[branch_column] <=
            minimum_width)
            branch_column = branch_column == first_column
                ? second_column : first_column;
        const double width =
            node.upper[branch_column] - node.lower[branch_column];
        if (!std::isfinite(width) || width <= minimum_width)
        {
            result.budget_exhausted = true;
            detached_open_nodes++;
            detached_open_bound = (std::min)(detached_open_bound, node.bound);
            break;
        }
        const double lower = node.lower[branch_column];
        const double upper = node.upper[branch_column];
        double split = node.primal[branch_column];
        if (!std::isfinite(split) ||
            split <= lower + 0.2 * width || split >= upper - 0.2 * width)
            split = lower + 0.5 * width;
        if (split <= lower + minimum_width ||
            split >= upper - minimum_width)
        {
            result.budget_exhausted = true;
            detached_open_nodes++;
            detached_open_bound = (std::min)(detached_open_bound, node.bound);
            break;
        }

        for (int side = 0; side < 2; ++side)
        {
            if (result.nodes >= node_limit ||
                elapsed_ms() >= total_limit_ms)
            {
                result.budget_exhausted = true;
                detached_open_nodes++;
                detached_open_bound =
                    (std::min)(detached_open_bound, node.bound);
                break;
            }
            BranchRltNode child;
            child.lower = node.lower;
            child.upper = node.upper;
            if (side == 0) child.upper[branch_column] = split;
            else child.lower[branch_column] = split;
            child.depth = node.depth + 1;
            child.serial = ++serial;
            const HighsBasis *warm_basis = node.basis.valid
                ? &node.basis : nullptr;
            const std::vector<std::vector<int> > *cross_plan =
                selected_cross_row_rlt_mode() !=
                    CrossRowRltMode::disabled
                ? &node.cut_cross_rows : nullptr;
            const std::vector<double> *separation_primal = &node.primal;
            IsolatedRltStatus status = IsolatedRltStatus::failed;
            HighsBasis resumed_basis;
            while (true)
            {
                const int remaining_ms = total_limit_ms -
                    static_cast<int>(elapsed_ms());
                if (remaining_ms <= 0)
                {
                    status = IsolatedRltStatus::interrupted;
                    break;
                }
                const double solve_seconds = (limits.borrow_remaining
                    ? remaining_ms
                    : (std::min)(lp_limit_ms, remaining_ms)) / 1000.0;
                if (!separation_primal) result.resumed_lp_solves++;
                status = evaluate_branch_rlt_node(
                    &child, instance, dataset_in, fixed_variables,
                    linear_cons_P, linear_cons_V,
                    P_lobo, P_upbo, V_lobo, V_upbo,
                    Ai_begin, Ai_end, Aj_begin, Aj_end,
                    n_alt_i, n_alt_j, mX, mY, noc, min_flag,
                    solve_seconds, separation_primal, cross_plan,
                    warm_basis, reduced_cost_tightening, &result);
                result.lp_solves++;
                if (status != IsolatedRltStatus::interrupted ||
                    !resume_interrupted || limits.borrow_remaining ||
                    (!child.continuation && !child.basis.valid))
                    break;
                resumed_basis = child.basis;
                warm_basis = &resumed_basis;
                separation_primal = nullptr;
                cross_plan = selected_cross_row_rlt_mode() !=
                    CrossRowRltMode::disabled
                    ? &child.cut_cross_rows : nullptr;
            }
            result.nodes++;
            if (status == IsolatedRltStatus::infeasible)
            {
                result.infeasible++;
                continue;
            }
            if (status == IsolatedRltStatus::interrupted)
            {
                result.budget_exhausted = true;
                detached_open_nodes++;
                detached_open_bound =
                    (std::min)(detached_open_bound, node.bound);
                break;
            }
            if (status == IsolatedRltStatus::failed)
            {
                result.lp_failed = true;
                result.budget_exhausted = true;
                detached_open_nodes++;
                detached_open_bound =
                    (std::min)(detached_open_bound, node.bound);
                break;
            }
            open.push(std::move(child));
        }
        if (result.budget_exhausted || result.lp_failed) break;
    }
    result.final_incumbent_valid = incumbent_valid;
    result.certified = incumbent_valid && open.empty() &&
        !result.budget_exhausted && !result.lp_failed;
    result.open_nodes = static_cast<int>(open.size()) + detached_open_nodes;
    if (!open.empty()) result.best_open_bound = open.top().bound;
    result.best_open_bound = (std::min)(result.best_open_bound,
                                        detached_open_bound);
    result.wall_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return result;
}
}

void reset_decideit_stats() { decideit_stats = DecideItRunStats{}; }
void set_decideit_stats(const DecideItRunStats &stats) { decideit_stats = stats; }
const DecideItRunStats &last_decideit_stats() { return decideit_stats; }

std::string selected_rlt_backend()
{
    const char *value = std::getenv("DECIDEIT_RLT_BACKEND");
    return value ? value : "highs";
}

extern "C"
{
    HighsSolver *createHighsSolver()
    {
        HighsSolver* solver = new HighsSolver();
        configure_rlt_highs(solver->highs);
        return solver;
    }

    void destroyHighsSolver(HighsSolver *solver)
    {
        delete solver;
    }

    void setLP(HighsSolver *solver,
               int num_col, int num_row,
               const int *start, const int *index, const double *value,
               const double *col_cost, const double *col_lower, const double *col_upper,
               const double *row_lower, const double *row_upper)
    {
        HighsLp lp;

        lp.num_col_ = num_col;
        lp.num_row_ = num_row;

        lp.col_cost_.assign(col_cost, col_cost + num_col);

        lp.col_lower_.assign(col_lower, col_lower + num_col);
        lp.col_upper_.assign(col_upper, col_upper + num_col);

        lp.row_lower_.assign(row_lower, row_lower + num_row);
        lp.row_upper_.assign(row_upper, row_upper + num_row);

        // lp.a_matrix_.format_ = MatrixFormat::kColwise;
        // lp.a_matrix_.start_.assign(start, start + num_col + 1);
        // lp.a_matrix_.index_.assign(index, index + start[num_col]);
        // lp.a_matrix_.value_.assign(value, value + start[num_col]);

        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.start_.assign(start, start + num_row + 1);
        lp.a_matrix_.index_.assign(index, index + start[num_row]);
        lp.a_matrix_.value_.assign(value, value + start[num_row]);

        // for (int i = 0; i < num_col; ++i)
        // {
        //     printf("col[%d]: [%g,%g]\n", i, lp.col_lower_[i], lp.col_upper_[i]);
        // }

        // for (int i = 0; i < num_row; ++i)
        // {
        //     printf("row[%d]: [%g,%g]\n", i, lp.row_lower_[i], lp.row_upper_[i]);
        // }

        lp.sense_ = ObjSense::kMinimize;
        highs_require(solver->highs.passModel(lp), "RLT passModel");
        configure_rlt_scale_for_columns(solver->highs, num_col);
        solver->loaded = true;
        solver->loaded_rows = num_row;
        solver->solved = false;
        solver->reduced_cost_snapshot_valid = false;
        solver->reduced_cost_objective =
            std::numeric_limits<double>::quiet_NaN();
        solver->reduced_cost_primal.clear();
        solver->reduced_cost_dual.clear();
        solver->reduced_cost_basis = HighsBasis{};
    }

    void setWarmStart(HighsSolver *solver, const double *x0, int num_col)
    {
        HighsSolution sol;
        sol.col_value.assign(x0, x0 + num_col);
        highs_require(solver->highs.setSolution(sol), "RLT initial solution");
    }

    void solve_with_highs(HighsSolver *solver, bool outputFlag)
    {
        if (outputFlag)
        {
            solver->highs.setOptionValue("output_flag", true);
        }
        else
        {
            solver->highs.setOptionValue("output_flag", false);
        }

        rlt_audit_snapshot(solver->highs);
        const bool audit=std::getenv("DECIDEIT_RLT_DIAGNOSTICS")!=nullptr;
        std::chrono::steady_clock::time_point audit_begin;
        if(audit){
            audit_begin=std::chrono::steady_clock::now();
            fprintf(stderr,"RLT_ATTEMPT {\"basis_valid_before\":%s,\"primal_valid_before\":%s,\"rows\":%d,\"columns\":%d}\n",
                solver->highs.getBasis().valid?"true":"false",solver->highs.getSolution().value_valid?"true":"false",int(solver->highs.getNumRow()),int(solver->highs.getNumCol()));
            fflush(stderr);
        }
        bool audit_fallback=false;
        auto status = solver->highs.run();
        if (status == HighsStatus::kError ||
            (solver->highs.getModelStatus() == HighsModelStatus::kOptimal && !highs_solution_valid(solver->highs)))
        {
            audit_fallback=true;
            if(audit){fprintf(stderr,"RLT_FALLBACK {\"event\":\"cold_retry\"}\n");fflush(stderr);}
            highs_require(solver->highs.clearSolver(), "RLT cold fallback");
            status = solver->highs.run();
        }
        highs_require(status, "RLT run");
        if (!highs_solution_valid(solver->highs))
            throw std::runtime_error("RLT has no certified LP solution: " + solver->highs.modelStatusToString(solver->highs.getModelStatus()));
        if(audit)rlt_audit_result(solver->highs,std::chrono::duration<double>(std::chrono::steady_clock::now()-audit_begin).count(),audit_fallback);
        solver->solved = true;
        if (!reduced_cost_snapshots_requested())
        {
            solver->reduced_cost_snapshot_valid = false;
            solver->reduced_cost_objective =
                std::numeric_limits<double>::quiet_NaN();
            solver->reduced_cost_primal.clear();
            solver->reduced_cost_dual.clear();
            solver->reduced_cost_basis = HighsBasis{};
            return;
        }
        const HighsSolution &solution = solver->highs.getSolution();
        const HighsBasis basis = solver->highs.getBasis();
        solver->reduced_cost_snapshot_valid =
            solver->highs.getModelStatus() == HighsModelStatus::kOptimal &&
            solution.value_valid && solution.dual_valid && basis.valid &&
            solution.col_value.size() == solution.col_dual.size() &&
            basis.col_status.size() == solution.col_value.size();
        if (solver->reduced_cost_snapshot_valid)
        {
            solver->reduced_cost_objective =
                solver->highs.getObjectiveValue();
            solver->reduced_cost_primal = solution.col_value;
            solver->reduced_cost_dual = solution.col_dual;
            solver->reduced_cost_basis = basis;
        }
        else
        {
            solver->reduced_cost_objective =
                std::numeric_limits<double>::quiet_NaN();
            solver->reduced_cost_primal.clear();
            solver->reduced_cost_dual.clear();
            solver->reduced_cost_basis = HighsBasis{};
        }
    }

    void getSolution(HighsSolver *solver, double *x_out, int num_col)
    {
        if (!solver->solved || !highs_solution_valid(solver->highs))
            throw std::runtime_error("Cannot read unsuccessful RLT solve");
        HighsSolution sol = solver->highs.getSolution();

        if (sol.col_value.size() == static_cast<size_t>(num_col))
        {
            std::memcpy(x_out, sol.col_value.data(), sizeof(double) * num_col);
        }
        else throw std::runtime_error("RLT solution dimension mismatch");
    }

    double getObjective(HighsSolver *solver)
    {
        if (!solver->solved || !highs_solution_valid(solver->highs))
            throw std::runtime_error("Cannot read unsuccessful RLT objective");
        return solver->highs.getObjectiveValue();
    }

} // extern "C"

/*
    tableau_construction() constructs the standard simplex tableau for pivoting,
    which includes the original linear constraints, finite upper bounds imposed on variables, and slacks.
    as it is "standard", each variable >= 0.

    inputs
        lp_X_std: the lp model in its standard form without slack variables。
        tableau: the tableau to be constructed.

    output
        tableau: the constructed tableau.
*/

void tableau_construction(LpModel *lp_X_std,
                          gsl_matrix *tableau)
{
    const int m = get_Nrows(lp_X_std) + 1;
    int i, j;
    std::vector<REAL> column(m);

    gsl_vector_view src_v;
    gsl_matrix_view dest_m;

    for (i = 1; i <= get_Ncolumns(lp_X_std); i++)
    {
        get_column(lp_X_std, i, column.data());
        src_v = gsl_vector_view_array(column.data() + 1, get_Nrows(lp_X_std));
        gsl_matrix_set_col(tableau, i, &src_v.vector);
    }

    j = 0;

    for (i = 1; i <= get_Nrows(lp_X_std); i++)
    {
        gsl_matrix_set(tableau, j++, 0, get_rh(lp_X_std, i));
    }

    dest_m = gsl_matrix_submatrix(tableau, 0, 1 + get_Ncolumns(lp_X_std), get_Nrows(lp_X_std), get_Nrows(lp_X_std));
    gsl_matrix_set_identity(&dest_m.matrix);
}

/*
    amc() implements ordinary alternating mountain climbing.  It does not
    inspect adjacent vertices; that optional, deliberately expensive step is
    kept separate so the two variants can be measured independently.

    inputs
        local_amc: the structure to save the identified heuristic solution;
        lp_X_std: the constructed standard LP for V
        lp_Y: the constructed LP for P;
        dataset_in: the structure containing the actual linear constraints, bounds and variables in X;
        noc: the # of linear constraints in X/V.

    outputs
        local_amc: saving the information concerning the PGM, its adjacent vertices, etc.
*/

bool amc(struct local_solution *local_amc,
         LpModel *lp_X_std, LpModel *lp_Y,
         struct dataset *dataset_in,
         const gsl_vector *initial_y,
         DecideItDirectionStats *stats,
         int maximum_rounds = 0)
{
    const int n_variable_x = get_Ncolumns(lp_X_std), n_variable_y = get_Ncolumns(lp_Y);

    REAL cbov = 1e7;

    gsl_vector *coef_x = managed::vector_calloc(n_variable_x + 1);
    gsl_vector *coef_y = managed::vector_calloc(n_variable_y + 1);
    gsl_vector *vertex_x = managed::vector_calloc(n_variable_x);
    gsl_vector *vertex_y = managed::vector_calloc(n_variable_y);
    gsl_vector *temp_y = managed::vector_calloc(n_variable_y);

    gsl_vector_view sub_coef_x = gsl_vector_subvector(coef_x, 1, n_variable_x);
    gsl_vector_view sub_coef_y = gsl_vector_subvector(coef_y, 1, n_variable_y);

    gsl_vector_set(vertex_y, 0, 1);
    if (initial_y && initial_y->size == temp_y->size)
        gsl_vector_memcpy(temp_y, initial_y);
    else
        gsl_vector_set_all(temp_y, 1);

    int amc_steps=0;
    while (gsl_blas_dnrm2(vertex_y) > DBLP_EPSILON)
    {
        if (++amc_steps>10000) throw std::runtime_error("AMC iteration limit exceeded");
        gsl_vector_memcpy(vertex_y, temp_y); // save the solution for comparison in the next iteration

        gsl_blas_dgemv(CblasNoTrans, 1, dataset_in->Q, temp_y, 0, vertex_x); // the coefficient of x, (Q * y)'
        gsl_vector_memcpy(&sub_coef_x.vector, vertex_x);
        set_obj_fn(lp_X_std, coef_x->data);
        solve(lp_X_std);
        if (stats) stats->block_lp_solves++;

        get_variables(lp_X_std, vertex_x->data);
        gsl_vector_add(vertex_x, dataset_in->xl);

        gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, vertex_x, 0, temp_y); // the coefficient of y, Q' * x
        gsl_vector_add(temp_y, dataset_in->dy);                            // the final coefficient of y, Q' * x + dy
        gsl_vector_memcpy(&sub_coef_y.vector, temp_y);
        set_obj_fn(lp_Y, coef_y->data);
        solve(lp_Y);
        if (stats)
        {
            stats->block_lp_solves++;
            stats->amc_rounds++;
        }
        get_variables(lp_Y, temp_y->data);

        gsl_vector_sub(vertex_y, temp_y);
        if (maximum_rounds > 0 && amc_steps >= maximum_rounds) break;
    }

    const bool converged = gsl_blas_dnrm2(vertex_y) <= DBLP_EPSILON;

    cbov = get_objective(lp_Y);
    resolve_zero_vector(vertex_x);
    resolve_zero_vector(temp_y);

    /* save results and free memory */
    local_amc->cbov = cbov;
    gsl_vector_memcpy(local_amc->x_solution, vertex_x);
    gsl_vector_memcpy(local_amc->y_solution, temp_y);

    managed::vector_free(temp_y);
    managed::vector_free(vertex_x);
    managed::vector_free(vertex_y);
    managed::vector_free(coef_x);
    managed::vector_free(coef_y);
    return converged;
}

/*
   Evaluate adjacent X vertices for the historical augmented-climb ablation
   and for guarded certification immediately before adding a polar cut.  The
   ordinary plain-climb path never calls this routine.  Returning on the first
   improvement reproduces the old better-first policy and makes the LP cost
   visible in DecideItDirectionStats.
*/
static bool improving_neighbor(
    struct dataset *dataset_in,
    struct local_solution *local_amc,
    LpModel *lp_Y,
    gsl_vector *better_y,
    DecideItDirectionStats *stats)
{
    const int n_variable_y = get_Ncolumns(lp_Y);
    gsl_vector *coef = managed::vector_calloc(1 + n_variable_y);
    gsl_vector_view objective = gsl_vector_subvector(coef, 1, n_variable_y);

    for (int j = 0; j < local_amc->n_adj_vertices; ++j)
    {
        gsl_vector_view neighbor = gsl_matrix_column(local_amc->adj_vertices, j);
        gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, &neighbor.vector,
                       0, &objective.vector);
        gsl_vector_add(&objective.vector, dataset_in->dy);
        set_obj_fn(lp_Y, coef->data);
        solve(lp_Y);
        if (stats) stats->neighbor_lp_solves++;
        if (local_amc->cbov - get_objective(lp_Y) > DBLP_EPSILON)
        {
            get_variables(lp_Y, better_y->data);
            managed::vector_free(coef);
            return true;
        }
    }

    managed::vector_free(coef);
    return false;
}

void packed2normal(
    const gsl_vector_int *a_start,
    size_t pointer,
    const gsl_vector *a_value,
    const gsl_vector_int *a_index,
    size_t rows,
    size_t column,
    size_t n_variable)
{
    size_t i, j;
    gsl_vector *v = managed::vector_alloc(n_variable);

    pointer--;

    for (i = 0; i < rows; i++)
    {
        column = column - (gsl_vector_int_get(a_start, pointer) - gsl_vector_int_get(a_start, pointer - 1));
        pointer--;
    }

    size_t ceiling;

    for (i = 0; i < rows; i++)
    {
        printf("\n\n");

        gsl_vector_set_zero(v);
        ceiling = gsl_vector_int_get(a_start, pointer + 1) - gsl_vector_int_get(a_start, pointer);

        for (j = 0; j < ceiling; j++)
        {
            gsl_vector_set(v, gsl_vector_int_get(a_index, column), gsl_vector_get(a_value, column));
            column++;
        }

        pointer++;

        for (j = 0; j < n_variable; j++)
        {
            printf("   %8.4f   ", gsl_vector_get(v, j));
        }
    }

    printf("\n");

    managed::vector_free(v);
}

int previousExistance(const size_t *records_pos, size_t current_pos, size_t num)
{
    size_t i;

    for (i = 0; i < num; i++)
    {
        if (records_pos[i] == current_pos)
        {
            break;
        }
    }

    return i;
}

void sparse_outer_product_simple(
    const size_t *u_idx, const double *u_val, size_t u_nnz,
    const size_t *v_idx, const double *v_val, size_t v_nnz,
    size_t *result_i, size_t *result_j, double *result_val,
    size_t *result_nnz)
{

    size_t count = 0, i, j;

    for (i = 0; i < u_nnz; i++)
    {
        for (j = 0; j < v_nnz; j++)
        {
            result_i[count] = u_idx[i];
            result_j[count] = v_idx[j];
            result_val[count] = u_val[i] * v_val[j];

            count++;
        }
    }

    *result_nnz = count;
}

void print_OuterProductResult(size_t actual_nnz, size_t *result_i, size_t *result_j, double *result_val)
{
    size_t i;

    printf("Simplified Version - Outer Product Result (nnz=%zu):\n", actual_nnz);

    for (i = 0; i < actual_nnz; i++)
    {
        printf("  (%zu, %zu) = %.4f\n", result_i[i], result_j[i], result_val[i]);
    }
}

#if 0 // Archived XX/YY full RLT; not part of the active paper algorithm.
void build_full_linearization(
    int n_alt1, int n_alt2,
    std::vector<LinCons> &linear_cons_P, std::vector<LinCons> &linear_cons_V,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    const int mX, const int mY,
    LinearizationParams *LinPara)
{
    LinPara->persistent.reset(); // A new min/max query has its own RLT and basis.
    const int nX = n_alt1 + n_alt2, nY = nX;
    int i, j, t, k, row, column, r, pointer, pos, loc, loc_i, loc_j, n_XY, n_XX, n_YY, n_total;

    n_total = LinPara->n_total;

    n_XY = LinPara->n_XY;
    n_XX = LinPara->n_XX;
    n_YY = LinPara->n_YY;

    gsl_vector *row_lower = LinPara->row_lower;
    gsl_vector *row_upper = LinPara->row_upper;
    gsl_vector *a_value = LinPara->a_value;
    gsl_vector_int *a_index = LinPara->a_index;
    gsl_vector_int *a_start = LinPara->a_start;

    gsl_vector *col_lower = LinPara->col_lower;
    gsl_vector *col_upper = LinPara->col_upper;

    /* set bounds */
    j = 0;

    for (i = Ai_begin; i <= Ai_end; i++, j++)
    {
        gsl_vector_set(col_lower, j, V_lobo[i]);
        gsl_vector_set(col_upper, j, V_upbo[i]);
        gsl_vector_set(col_lower, j + nX, P_lobo[i]);
        gsl_vector_set(col_upper, j + nX, P_upbo[i]);
    }

    for (i = Aj_begin; i <= Aj_end; i++, j++)
    {
        gsl_vector_set(col_lower, j, V_lobo[i]);
        gsl_vector_set(col_upper, j, V_upbo[i]);
        gsl_vector_set(col_lower, j + nX, P_lobo[i]);
        gsl_vector_set(col_upper, j + nX, P_upbo[i]);
    }

    for (j = nX + nY; j < n_total; j++)
    {
        gsl_vector_set(col_lower, j, 0);
        gsl_vector_set(col_upper, j, 1);
    }

    pointer = 0; // the position in a_start
    gsl_vector_int_set(a_start, pointer++, 0);
    row = 0;    // the # of rows
    column = 0; // the position in a_index and a_value

    /* (2.2a) X: V */
    for (r = 0; r < mX; ++r)
    {
        gsl_vector_int_set(a_index, column, linear_cons_V[r].j1 - 1);
        gsl_vector_set(a_value, column++, -1.0);
        gsl_vector_int_set(a_index, column, linear_cons_V[r].j2 - 1);
        gsl_vector_set(a_value, column++, 1.0);
        gsl_vector_set(row_lower, row, -linear_cons_V[r].b);
        gsl_vector_set(row_upper, row++, 1.0e30);
        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);

        pointer++;
    }

    /* (2.2a) Y: P */
    for (r = 0; r < mY; ++r)
    {
        gsl_vector_int_set(a_index, column, nX + linear_cons_P[r].j1 - 1);
        gsl_vector_set(a_value, column++, -1.0);
        gsl_vector_int_set(a_index, column, nX + linear_cons_P[r].j2 - 1);
        gsl_vector_set(a_value, column++, 1.0);
        gsl_vector_set(row_lower, row, -linear_cons_P[r].b);
        gsl_vector_set(row_upper, row++, 1.0e30);
        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);

        pointer++;
    }

    int entry_1, entry_2;

    /* (2.2b) */
    for (i = 0; i < nX; ++i)
    {
        loc = location(i, Ai_begin, Aj_begin, n_alt1);

        for (r = 0; r < mX; ++r)
        {
            j = linear_cons_V[r].j1;
            pos = sym_index1(j, i + 1, nX) - 1;
            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, -1.0);

            j = linear_cons_V[r].j2;
            pos = sym_index1(j, i + 1, nX) - 1;
            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, 1.0);

            entry_1 = linear_cons_V[r].j1 - 1;
            gsl_vector_int_set(a_index, column, entry_1);
            gsl_vector_set(a_value, column++, V_lobo[loc]);

            entry_2 = linear_cons_V[r].j2 - 1;
            gsl_vector_int_set(a_index, column, entry_2);
            gsl_vector_set(a_value, column++, -V_lobo[loc]);

            if (i != entry_1 && i != entry_2)
            {
                gsl_vector_int_set(a_index, column, i);
                gsl_vector_set(a_value, column++, linear_cons_V[r].b);

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            }
            else
            {
                if (i == entry_1)
                {
                    gsl_vector_set(a_value, column - 2, gsl_vector_get(a_value, column - 2) + linear_cons_V[r].b);
                }
                else
                {
                    gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) + linear_cons_V[r].b);
                }

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 4);
            }

            pointer++;

            gsl_vector_set(row_lower, row, V_lobo[loc] * linear_cons_V[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }

        for (r = 0; r < mY; ++r)
        {
            j = linear_cons_P[r].j1 - 1;
            pos = xy_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            j = linear_cons_P[r].j2 - 1;
            pos = xy_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            gsl_vector_int_set(a_index, column, nX + linear_cons_P[r].j1 - 1);
            gsl_vector_set(a_value, column++, V_lobo[loc]);

            gsl_vector_int_set(a_index, column, nX + linear_cons_P[r].j2 - 1);
            gsl_vector_set(a_value, column++, -V_lobo[loc]);

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, linear_cons_P[r].b);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            pointer++;

            gsl_vector_set(row_lower, row, V_lobo[loc] * linear_cons_P[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2c) */

    for (i = 0; i < nX; ++i)
    {
        loc = location(i, Ai_begin, Aj_begin, n_alt1);

        for (r = 0; r < mX; ++r)
        {
            j = linear_cons_V[r].j1;
            pos = sym_index1(j, i + 1, nX) - 1;
            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, 1.0);

            j = linear_cons_V[r].j2;
            pos = sym_index1(j, i + 1, nX) - 1;
            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, -1.0);

            entry_1 = linear_cons_V[r].j1 - 1;
            gsl_vector_int_set(a_index, column, entry_1);
            gsl_vector_set(a_value, column++, -V_upbo[loc]);

            entry_2 = linear_cons_V[r].j2 - 1;
            gsl_vector_int_set(a_index, column, entry_2);
            gsl_vector_set(a_value, column++, V_upbo[loc]);

            if (i != entry_1 && i != entry_2)
            {
                gsl_vector_int_set(a_index, column, i);
                gsl_vector_set(a_value, column++, -linear_cons_V[r].b);

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            }
            else
            {
                if (i == entry_1)
                {
                    gsl_vector_set(a_value, column - 2, gsl_vector_get(a_value, column - 2) - linear_cons_V[r].b);
                }
                else
                {
                    gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) - linear_cons_V[r].b);
                }

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 4);
            }

            pointer++;

            gsl_vector_set(row_lower, row, -V_upbo[loc] * linear_cons_V[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }

        for (r = 0; r < mY; ++r)
        {
            j = linear_cons_P[r].j1 - 1;
            pos = xy_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            j = linear_cons_P[r].j2 - 1;
            pos = xy_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            gsl_vector_int_set(a_index, column, nX + linear_cons_P[r].j1 - 1);
            gsl_vector_set(a_value, column++, -V_upbo[loc]);

            gsl_vector_int_set(a_index, column, nX + linear_cons_P[r].j2 - 1);
            gsl_vector_set(a_value, column++, V_upbo[loc]);

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, -linear_cons_P[r].b);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            pointer++;

            gsl_vector_set(row_lower, row, -V_upbo[loc] * linear_cons_P[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2d) */
    for (j = 0; j < nY; ++j)
    {
        loc = location(j, Ai_begin, Aj_begin, n_alt1);

        for (r = 0; r < mX; ++r)
        {
            t = linear_cons_V[r].j1 - 1;
            pos = xy_index0(t, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            t = linear_cons_V[r].j2 - 1;
            pos = xy_index0(t, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            gsl_vector_int_set(a_index, column, linear_cons_V[r].j1 - 1);
            gsl_vector_set(a_value, column++, P_lobo[loc]);

            gsl_vector_int_set(a_index, column, linear_cons_V[r].j2 - 1);
            gsl_vector_set(a_value, column++, -P_lobo[loc]);

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, linear_cons_V[r].b);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            pointer++;

            gsl_vector_set(row_lower, row, P_lobo[loc] * linear_cons_V[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }

        for (r = 0; r < mY; ++r)
        {
            k = linear_cons_P[r].j1 - 1;
            pos = sym_index0(k, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            k = linear_cons_P[r].j2 - 1;
            pos = sym_index0(k, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            entry_1 = linear_cons_P[r].j1 - 1;
            gsl_vector_int_set(a_index, column, nX + entry_1);
            gsl_vector_set(a_value, column++, P_lobo[loc]);

            entry_2 = linear_cons_P[r].j2 - 1;
            gsl_vector_int_set(a_index, column, nX + entry_2);
            gsl_vector_set(a_value, column++, -P_lobo[loc]);

            if (j != entry_1 && j != entry_2)
            {
                gsl_vector_int_set(a_index, column, nX + j);
                gsl_vector_set(a_value, column++, linear_cons_P[r].b);

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            }
            else
            {
                if (j == entry_1)
                {
                    gsl_vector_set(a_value, column - 2, gsl_vector_get(a_value, column - 2) + linear_cons_P[r].b);
                }
                else
                {
                    gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) + linear_cons_P[r].b);
                }

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 4);
            }

            pointer++;

            gsl_vector_set(row_lower, row, P_lobo[loc] * linear_cons_P[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2e) */
    for (j = 0; j < nY; ++j)
    {
        loc = location(j, Ai_begin, Aj_begin, n_alt1);

        for (r = 0; r < mX; ++r)
        {
            t = linear_cons_V[r].j1 - 1;
            pos = xy_index0(t, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            t = linear_cons_V[r].j2 - 1;
            pos = xy_index0(t, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            gsl_vector_int_set(a_index, column, linear_cons_V[r].j1 - 1);
            gsl_vector_set(a_value, column++, -P_upbo[loc]);

            gsl_vector_int_set(a_index, column, linear_cons_V[r].j2 - 1);
            gsl_vector_set(a_value, column++, P_upbo[loc]);

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, -linear_cons_V[r].b);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            pointer++;

            gsl_vector_set(row_lower, row, -P_upbo[loc] * linear_cons_V[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }

        for (r = 0; r < mY; ++r)
        {
            k = linear_cons_P[r].j1 - 1;
            pos = sym_index0(k, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            k = linear_cons_P[r].j2 - 1;
            pos = sym_index0(k, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            entry_1 = linear_cons_P[r].j1 - 1;
            gsl_vector_int_set(a_index, column, nX + entry_1);
            gsl_vector_set(a_value, column++, -P_upbo[loc]);

            entry_2 = linear_cons_P[r].j2 - 1;
            gsl_vector_int_set(a_index, column, nX + entry_2);
            gsl_vector_set(a_value, column++, P_upbo[loc]);

            if (j != entry_1 && j != entry_2)
            {
                gsl_vector_int_set(a_index, column, nX + j);
                gsl_vector_set(a_value, column++, -linear_cons_P[r].b);

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 5);
            }
            else
            {
                if (j == entry_1)
                {
                    gsl_vector_set(a_value, column - 2, gsl_vector_get(a_value, column - 2) - linear_cons_P[r].b);
                }
                else
                {
                    gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) - linear_cons_P[r].b);
                }

                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 4);
            }

            pointer++;

            gsl_vector_set(row_lower, row, -P_upbo[loc] * linear_cons_P[r].b);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2f) XX, lower xl */
    pos = 0;

    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = i; j < nX; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, 1.0);

            pos++;

            gsl_vector_int_set(a_index, column, j);
            gsl_vector_set(a_value, column++, -V_lobo[loc_i]);

            if (i != j)
            {
                gsl_vector_int_set(a_index, column, i);
                gsl_vector_set(a_value, column++, -V_lobo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            }
            else
            {
                gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) - V_lobo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);
            }

            pointer++;

            gsl_vector_set(row_lower, row, -V_lobo[loc_i] * V_lobo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2g) YY, lower yl */
    pos = 0;

    for (i = 0; i < nY; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = i; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            pos++;

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, -P_lobo[loc_i]);

            if (i != j)
            {
                gsl_vector_int_set(a_index, column, nX + i);
                gsl_vector_set(a_value, column++, -P_lobo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            }
            else
            {
                gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) - P_lobo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);
            }

            pointer++;

            gsl_vector_set(row_lower, row, -P_lobo[loc_i] * P_lobo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2h) XX, upper xu */
    pos = 0;

    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = i; j < nX; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, 1.0);

            pos++;

            gsl_vector_int_set(a_index, column, j);
            gsl_vector_set(a_value, column++, -V_upbo[loc_i]);

            if (i != j)
            {
                gsl_vector_int_set(a_index, column, i);
                gsl_vector_set(a_value, column++, -V_upbo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            }
            else
            {
                gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) - V_upbo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);
            }

            pointer++;

            gsl_vector_set(row_lower, row, -V_upbo[loc_i] * V_upbo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2i) YY, upper yu */
    pos = 0;

    for (i = 0; i < nY; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = i; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            pos++;

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, -P_upbo[loc_i]);

            if (i != j)
            {
                gsl_vector_int_set(a_index, column, nX + i);
                gsl_vector_set(a_value, column++, -P_upbo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            }
            else
            {
                gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) - P_upbo[loc_j]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);
            }

            pointer++;

            gsl_vector_set(row_lower, row, -P_upbo[loc_i] * P_upbo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2j) */
    pos = 0;

    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = 0; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            pos++;

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, -V_lobo[loc_i]);

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, -P_lobo[loc_j]);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            pointer++;

            gsl_vector_set(row_lower, row, -V_lobo[loc_i] * P_lobo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2k) */
    pos = 0;

    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = 0; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);

            pos++;

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, -V_upbo[loc_i]);

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, -P_upbo[loc_j]);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            pointer++;

            gsl_vector_set(row_lower, row, -V_upbo[loc_i] * P_upbo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2l) */
    pos = 0;

    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = 0; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            pos++;

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, V_lobo[loc_i]);

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, P_upbo[loc_j]);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            pointer++;

            gsl_vector_set(row_lower, row, V_lobo[loc_i] * P_upbo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2m) */
    pos = 0;

    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = 0; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            pos++;

            gsl_vector_int_set(a_index, column, nX + j);
            gsl_vector_set(a_value, column++, V_upbo[loc_i]);

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, P_lobo[loc_j]);

            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            pointer++;

            gsl_vector_set(row_lower, row, V_upbo[loc_i] * P_lobo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2n) */
    for (i = 0; i < nX; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = 0; j < nX; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            pos = sym_index0(i, j, nX);
            gsl_vector_int_set(a_index, column, n_XX + pos);
            gsl_vector_set(a_value, column++, -1.0);

            pos++;

            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, V_upbo[loc_j]);

            if (i != j)
            {
                gsl_vector_int_set(a_index, column, j);
                gsl_vector_set(a_value, column++, V_lobo[loc_i]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            }
            else
            {
                gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) + V_lobo[loc_i]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);
            }

            pointer++;

            gsl_vector_set(row_lower, row, V_lobo[loc_i] * V_upbo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2o) */
    for (i = 0; i < nY; ++i)
    {
        loc_i = location(i, Ai_begin, Aj_begin, n_alt1);

        for (j = 0; j < nY; ++j)
        {
            loc_j = location(j, Ai_begin, Aj_begin, n_alt1);

            pos = sym_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, -1.0);

            gsl_vector_int_set(a_index, column, nX + i);
            gsl_vector_set(a_value, column++, P_upbo[loc_j]);

            if (i != j)
            {
                gsl_vector_int_set(a_index, column, nX + j);
                gsl_vector_set(a_value, column++, P_lobo[loc_i]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 3);
            }
            else
            {
                gsl_vector_set(a_value, column - 1, gsl_vector_get(a_value, column - 1) + P_lobo[loc_i]);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 2);
            }

            pointer++;

            gsl_vector_set(row_lower, row, P_lobo[loc_i] * P_upbo[loc_j]);
            gsl_vector_set(row_upper, row++, 1.0e30);
        }
    }

    /* (2.2p) */

    size_t ii, jj, actual_nnz, indicator, temp, counter = 0;
    size_t u_idx[2], v_idx[2], result_i[4], result_j[4], records_pos[4];
    const double u_val[2] = {1.0, -1.0}, v_val[2] = {1.0, -1.0};
    double result_val[4];

    for (i = 0; i < (mX + mY); ++i)
    {
        for (j = i; j < (mX + mY); ++j)
        {
            if (i < mX)
            {
                if (j < mX)
                {
                    /* x'[A'(i)A(j)]x */
                    u_idx[0] = linear_cons_V[i].j1 - 1;
                    u_idx[1] = linear_cons_V[i].j2 - 1;
                    v_idx[0] = linear_cons_V[j].j1 - 1;
                    v_idx[1] = linear_cons_V[j].j2 - 1;

                    sparse_outer_product_simple(u_idx, u_val, 2, v_idx, v_val, 2, result_i, result_j, result_val, &actual_nnz);

                    // print_OuterProductResult(actual_nnz, result_i, result_j, result_val);
                    // printf("\n");

                    counter = 0;

                    for (size_t k = 0; k < actual_nnz; k++)
                    {
                        pos = sym_index0(result_i[k], result_j[k], nX);

                        if (k > 0)
                        {
                            indicator = previousExistance(records_pos, pos, counter);

                            if (indicator == counter)
                            {
                                records_pos[counter++] = pos;

                                gsl_vector_int_set(a_index, column, n_XX + pos);
                                gsl_vector_set(a_value, column++, result_val[k]);
                                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 1);
                            }
                            else
                            {
                                temp = column - counter + indicator;
                                gsl_vector_set(a_value, temp, gsl_vector_get(a_value, temp) + result_val[k]);
                            }
                        }
                        else
                        {
                            records_pos[counter++] = pos;

                            gsl_vector_int_set(a_index, column, n_XX + pos);
                            gsl_vector_set(a_value, column++, result_val[k]);
                            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 1);
                        }
                    }

                    // -[b(1i)A(1j) + b(1j)A(1i)]x
                    gsl_vector_int_set(a_index, column, linear_cons_V[j].j1 - 1);
                    gsl_vector_set(a_value, column++, -linear_cons_V[i].b);
                    gsl_vector_int_set(a_index, column, linear_cons_V[j].j2 - 1);
                    gsl_vector_set(a_value, column++, linear_cons_V[i].b);
                    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 2);

                    counter = 2;
                    records_pos[0] = linear_cons_V[j].j1;
                    records_pos[1] = linear_cons_V[j].j2;

                    indicator = previousExistance(records_pos, linear_cons_V[i].j1, counter);

                    if (indicator == counter)
                    {
                        records_pos[counter++] = linear_cons_V[i].j1;

                        gsl_vector_int_set(a_index, column, linear_cons_V[i].j1 - 1);
                        gsl_vector_set(a_value, column++, -linear_cons_V[j].b);
                        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 1);
                    }
                    else
                    {
                        temp = column - counter + indicator;
                        gsl_vector_set(a_value, temp, gsl_vector_get(a_value, temp) - linear_cons_V[j].b);
                    }

                    indicator = previousExistance(records_pos, linear_cons_V[i].j2, counter);

                    if (indicator == counter)
                    {
                        gsl_vector_int_set(a_index, column, linear_cons_V[i].j2 - 1);
                        gsl_vector_set(a_value, column++, linear_cons_V[j].b);
                        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 1);
                    }
                    else
                    {
                        temp = column - counter + indicator;
                        gsl_vector_set(a_value, temp, gsl_vector_get(a_value, temp) + linear_cons_V[j].b);
                    }

                    pointer++;

                    // - b(1i)b(1j)
                    gsl_vector_set(row_lower, row, -linear_cons_V[i].b * linear_cons_V[j].b);
                    gsl_vector_set(row_upper, row++, 1.0e30);
                }
                else
                {
                    jj = j - mX;

                    // x'[A'(1i)D(1j) + A'(1j)D(1i)]y
                    u_idx[0] = linear_cons_V[i].j1 - 1;
                    u_idx[1] = linear_cons_V[i].j2 - 1;
                    v_idx[0] = linear_cons_P[jj].j1 - 1;
                    v_idx[1] = linear_cons_P[jj].j2 - 1;

                    sparse_outer_product_simple(u_idx, u_val, 2, v_idx, v_val, 2, result_i, result_j, result_val, &actual_nnz);

                    // print_OuterProductResult(actual_nnz, result_i, result_j, result_val);
                    // printf("\n");

                    counter = 0;

                    for (size_t k = 0; k < actual_nnz; k++)
                    {
                        pos = xy_index0(result_i[k], result_j[k], nY);

                        gsl_vector_int_set(a_index, column, n_XY + pos);
                        gsl_vector_set(a_value, column++, result_val[k]);
                    }

                    // -[b(1j)A(1i)]x
                    gsl_vector_int_set(a_index, column, linear_cons_V[i].j1 - 1);
                    gsl_vector_set(a_value, column++, -linear_cons_P[jj].b);
                    gsl_vector_int_set(a_index, column, linear_cons_V[i].j2 - 1);
                    gsl_vector_set(a_value, column++, linear_cons_P[jj].b);

                    // -[b(1i)D(1j)]y
                    gsl_vector_int_set(a_index, column, nX + linear_cons_P[jj].j1 - 1);
                    gsl_vector_set(a_value, column++, -linear_cons_V[i].b);
                    gsl_vector_int_set(a_index, column, nX + linear_cons_P[jj].j2 - 1);
                    gsl_vector_set(a_value, column++, linear_cons_V[i].b);

                    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 8);
                    pointer++;

                    // -b(1i)b(1j)
                    gsl_vector_set(row_lower, row, -linear_cons_V[i].b * linear_cons_P[jj].b);
                    gsl_vector_set(row_upper, row++, 1.0e30);
                }
            }
            else
            {
                ii = i - mX;
                jj = j - mX;

                // y'[D'(1i)D(1j)]y
                u_idx[0] = linear_cons_P[ii].j1 - 1;
                u_idx[1] = linear_cons_P[ii].j2 - 1;
                v_idx[0] = linear_cons_P[jj].j1 - 1;
                v_idx[1] = linear_cons_P[jj].j2 - 1;

                sparse_outer_product_simple(u_idx, u_val, 2, v_idx, v_val, 2, result_i, result_j, result_val, &actual_nnz);

                // print_OuterProductResult(actual_nnz, result_i, result_j, result_val);
                // printf("\n");

                counter = 0;

                for (size_t k = 0; k < actual_nnz; k++)
                {
                    pos = sym_index0(result_i[k], result_j[k], nY);

                    if (k > 0)
                    {
                        indicator = previousExistance(records_pos, pos, counter);

                        if (indicator == counter)
                        {
                            records_pos[counter++] = pos;

                            gsl_vector_int_set(a_index, column, n_YY + pos);
                            gsl_vector_set(a_value, column++, result_val[k]);
                            gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 1);
                        }
                        else
                        {
                            temp = column - counter + indicator;
                            gsl_vector_set(a_value, temp, gsl_vector_get(a_value, temp) + result_val[k]);
                        }
                    }
                    else
                    {
                        records_pos[counter++] = pos;

                        gsl_vector_int_set(a_index, column, n_YY + pos);
                        gsl_vector_set(a_value, column++, result_val[k]);
                        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + 1);
                    }
                }

                // -[b(1i)D(1j) + b(1j)D(1i)]y
                gsl_vector_int_set(a_index, column, nX + linear_cons_P[jj].j1 - 1);
                gsl_vector_set(a_value, column++, -linear_cons_P[ii].b);
                gsl_vector_int_set(a_index, column, nX + linear_cons_P[jj].j2 - 1);
                gsl_vector_set(a_value, column++, linear_cons_P[ii].b);
                gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 2);

                counter = 2;
                records_pos[0] = linear_cons_P[jj].j1;
                records_pos[1] = linear_cons_P[jj].j2;

                indicator = previousExistance(records_pos, linear_cons_P[ii].j1, counter);

                if (indicator == counter)
                {
                    records_pos[counter++] = linear_cons_P[ii].j1;

                    gsl_vector_int_set(a_index, column, nX + linear_cons_P[ii].j1 - 1);
                    gsl_vector_set(a_value, column++, -linear_cons_P[jj].b);
                    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 1);
                }
                else
                {
                    temp = column - counter + indicator;
                    gsl_vector_set(a_value, temp, gsl_vector_get(a_value, temp) - linear_cons_P[jj].b);
                }

                indicator = previousExistance(records_pos, linear_cons_P[ii].j2, counter);

                if (indicator == counter)
                {
                    gsl_vector_int_set(a_index, column, nX + linear_cons_P[ii].j2 - 1);
                    gsl_vector_set(a_value, column++, linear_cons_P[jj].b);
                    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer) + 1);
                }
                else
                {
                    temp = column - counter + indicator;
                    gsl_vector_set(a_value, temp, gsl_vector_get(a_value, temp) + linear_cons_P[jj].b);
                }

                pointer++;

                // -b(1i)b(1j)
                gsl_vector_set(row_lower, row, -linear_cons_P[ii].b * linear_cons_P[jj].b);
                gsl_vector_set(row_upper, row++, 1.0e30);
            }

            // packed2normal(a_start, pointer, a_value, a_index, 1, column, n_total);
        }
    }

    /* (2.2q) */
    for (i = 0; i < nX; ++i)
    {
        pos = i * nY;

        for (r = 0; r < n_alt1; ++r, ++pos)
        {
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);
        }

        gsl_vector_int_set(a_index, column, i);
        gsl_vector_set(a_value, column++, -1.0);

        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + n_alt1 + 1);
        pointer++;

        gsl_vector_set(row_lower, row, 0);
        gsl_vector_set(row_upper, row++, 0);
    }

    for (i = 0; i < nX; ++i)
    {
        pos = i * nY + n_alt1;

        for (r = 0; r < n_alt2; ++r, ++pos)
        {
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, 1.0);
        }

        gsl_vector_int_set(a_index, column, i);
        gsl_vector_set(a_value, column++, -1.0);

        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + n_alt2 + 1);
        pointer++;

        gsl_vector_set(row_lower, row, 0);
        gsl_vector_set(row_upper, row++, 0);
    }

    /* (2.2r) */
    for (i = 0; i < nY; ++i)
    {
        for (j = 0; j < n_alt1; ++j)
        {
            pos = sym_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, 1.0);
        }

        gsl_vector_int_set(a_index, column, nX + i);
        gsl_vector_set(a_value, column++, -1.0);

        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + n_alt1 + 1);
        pointer++;

        gsl_vector_set(row_lower, row, 0);
        gsl_vector_set(row_upper, row++, 0);
    }

    for (i = 0; i < nY; ++i)
    {
        for (j = n_alt1; j < nY; ++j)
        {
            pos = sym_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_YY + pos);
            gsl_vector_set(a_value, column++, 1.0);
        }

        gsl_vector_int_set(a_index, column, nX + i);
        gsl_vector_set(a_value, column++, -1.0);

        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + n_alt2 + 1);
        pointer++;

        gsl_vector_set(row_lower, row, 0);
        gsl_vector_set(row_upper, row++, 0);
    }

    /* set two normalizations as inequalities */
    for (r = 0; r < n_alt1; ++r)
    {
        gsl_vector_int_set(a_index, column, nX + r);
        gsl_vector_set(a_value, column++, 1.0);
    }

    gsl_vector_set(row_lower, row, 1.0);
    gsl_vector_set(row_upper, row++, 1.0);
    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + n_alt1);

    pointer++;

    for (r = 0; r < n_alt2; ++r)
    {
        gsl_vector_int_set(a_index, column, nX + n_alt1 + r);
        gsl_vector_set(a_value, column++, 1.0);
    }

    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + n_alt2);

    pointer++;

    gsl_vector_set(row_lower, row, 1.0);
    gsl_vector_set(row_upper, row++, 1.0);

    LinPara->num_nz = gsl_vector_int_get(a_start, pointer - 1);
    LinPara->num_col = n_total;
    LinPara->num_row = row;
    LinPara->offset = 0.0;

    LinPara->a_start_column = pointer;
    LinPara->a_index_column = column;

    row_lower = NULL;
    row_upper = NULL;
    a_value = NULL;
    a_index = NULL;
    a_start = NULL;
    col_lower = NULL;
    col_upper = NULL;
}
#endif

namespace
{
struct RltBuilder
{
    LinearizationParams *parameters;
    size_t row = 0;
    size_t nonzero = 0;

    explicit RltBuilder(LinearizationParams *value) : parameters(value) {}

    void append(std::vector<std::pair<int, double> > terms, double lower, double upper)
    {
        std::sort(terms.begin(), terms.end());
        std::vector<std::pair<int, double> > merged;
        for (const auto &term : terms)
        {
            if (std::fabs(term.second) <= 1e-15) continue;
            if (!merged.empty() && merged.back().first == term.first)
                merged.back().second += term.second;
            else
                merged.push_back(term);
        }
        merged.erase(std::remove_if(merged.begin(), merged.end(),
            [](const std::pair<int, double> &term) {
                return std::fabs(term.second) <= 1e-15;
            }), merged.end());
        ensure_rlt_capacity(parameters, row + 2, nonzero + merged.size() + 1);
        gsl_vector_int_set(parameters->a_start, row, nonzero);
        for (const auto &term : merged)
        {
            if (term.first < 0 || term.first >= parameters->n_total ||
                !std::isfinite(term.second))
                throw std::runtime_error("Invalid cross-term RLT coefficient");
            gsl_vector_int_set(parameters->a_index, nonzero, term.first);
            gsl_vector_set(parameters->a_value, nonzero++, term.second);
        }
        gsl_vector_set(parameters->row_lower, row, lower);
        gsl_vector_set(parameters->row_upper, row, upper);
        ++row;
        gsl_vector_int_set(parameters->a_start, row, nonzero);
    }
};

bool append_cross_row_product(
    RltBuilder &builder,
    const std::vector<std::pair<int, double> > &x_terms, double x_rhs,
    const LinCons &y_row, int n, const double *separation_values = nullptr,
    size_t separation_size = 0, bool violated_only = false)
{
    std::vector<std::pair<int, double> > terms;
    terms.reserve(x_terms.size() * y_row.index.size() +
                  x_terms.size() + y_row.index.size());
    const int product_offset = 2 * n;
    for (const auto &x_term : x_terms)
    {
        terms.push_back({x_term.first, -y_row.b * x_term.second});
        for (size_t k = 0; k < y_row.index.size(); ++k)
        {
            const int y = y_row.index[k] - 1;
            if (x_term.first < 0 || x_term.first >= n || y < 0 || y >= n)
                throw std::runtime_error("Cross-row RLT index out of range");
            terms.push_back({product_offset + x_term.first * n + y,
                             x_term.second * y_row.coefficient[k]});
        }
    }
    for (size_t k = 0; k < y_row.index.size(); ++k)
        terms.push_back({n + y_row.index[k] - 1,
                         -x_rhs * y_row.coefficient[k]});
    const double lower = -x_rhs * y_row.b;
    if (violated_only)
    {
        if (!separation_values || separation_size <
                static_cast<size_t>(builder.parameters->num_col))
            return false;
        double activity = 0.0;
        for (const auto &term : terms)
            activity += term.second * separation_values[term.first];
        if (!std::isfinite(activity) || activity >= lower - 1e-8)
            return false;
    }
    builder.append(terms, lower, 1.0e30);
    return true;
}
} // namespace

void build_linearization(
    int n_alt1, int n_alt2,
    std::vector<LinCons> &linear_cons_P, std::vector<LinCons> &linear_cons_V,
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    const int mX, const int mY,
    LinearizationParams *LinPara)
{
    LinPara->persistent.reset();
    const int n = n_alt1 + n_alt2;
    if (mX != static_cast<int>(linear_cons_V.size()) ||
        mY != static_cast<int>(linear_cons_P.size()))
        throw std::runtime_error("RLT source-row count mismatch");
    const double infinity = 1.0e30;
    const int w_offset = 2 * n;
    auto w = [=](int i, int j) { return w_offset + i * n + j; };

    std::vector<double> xl(n), xu(n), yl(n), yu(n);
    for (int i = 0; i < n; ++i)
    {
        const int original = location(i, Ai_begin, Aj_begin, n_alt1);
        xl[i] = V_lobo[original]; xu[i] = V_upbo[original];
        yl[i] = P_lobo[original]; yu[i] = P_upbo[original];
        gsl_vector_set(LinPara->col_lower, i, xl[i]);
        gsl_vector_set(LinPara->col_upper, i, xu[i]);
        gsl_vector_set(LinPara->col_lower, n + i, yl[i]);
        gsl_vector_set(LinPara->col_upper, n + i, yu[i]);
    }
    for (int column = w_offset; column < LinPara->n_total; ++column)
    {
        gsl_vector_set(LinPara->col_lower, column, 0.0);
        gsl_vector_set(LinPara->col_upper, column, 1.0);
    }

    const CrossRowRltMode cross_mode = selected_cross_row_rlt_mode();
    const bool add_cross_rows = cross_mode != CrossRowRltMode::disabled;
    if (add_cross_rows)
    {
        size_t cross_nonzeros = 0;
        for (const LinCons &x_row : linear_cons_V)
            for (const LinCons &y_row : linear_cons_P)
                cross_nonzeros += x_row.index.size() * y_row.index.size() +
                                  x_row.index.size() + y_row.index.size();
        ensure_rlt_capacity(
            LinPara,
            LinPara->row_lower->size +
                linear_cons_V.size() * linear_cons_P.size(),
            LinPara->a_value->size + cross_nonzeros);
    }

    RltBuilder builder{LinPara};
    for (const LinCons &source : linear_cons_V)
    {
        std::vector<std::pair<int, double> > terms;
        for (size_t k = 0; k < source.index.size(); ++k)
            terms.push_back({source.index[k] - 1, source.coefficient[k]});
        builder.append(terms, -infinity, source.b);
    }
    for (const LinCons &source : linear_cons_P)
    {
        std::vector<std::pair<int, double> > terms;
        for (size_t k = 0; k < source.index.size(); ++k)
            terms.push_back({n + source.index[k] - 1, source.coefficient[k]});
        builder.append(terms, -infinity, source.b);
    }

    if (add_cross_rows)
    {
        for (const LinCons &x_row : linear_cons_V)
        {
            std::vector<std::pair<int, double> > x_terms;
            x_terms.reserve(x_row.index.size());
            for (size_t k = 0; k < x_row.index.size(); ++k)
                x_terms.push_back(
                    {x_row.index[k] - 1, x_row.coefficient[k]});
            for (const LinCons &y_row : linear_cons_P)
                append_cross_row_product(builder, x_terms, x_row.b,
                                         y_row, n);
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            builder.append({{w(i,j), 1.0}, {i, -yl[j]}, {n+j, -xl[i]}},
                           -xl[i] * yl[j], infinity);
            builder.append({{w(i,j), 1.0}, {i, -yu[j]}, {n+j, -xu[i]}},
                           -xu[i] * yu[j], infinity);
            builder.append({{w(i,j), 1.0}, {i, -yl[j]}, {n+j, -xu[i]}},
                           -infinity, -xu[i] * yl[j]);
            builder.append({{w(i,j), 1.0}, {i, -yu[j]}, {n+j, -xl[i]}},
                           -infinity, -xl[i] * yu[j]);
        }
    }

    for (const LinCons &source : linear_cons_V)
    {
        for (int j = 0; j < n; ++j)
        {
            std::vector<std::pair<int, double> > low, high;
            for (size_t k = 0; k < source.index.size(); ++k)
            {
                const int i = source.index[k] - 1;
                const double a = source.coefficient[k];
                low.push_back({w(i,j), -a}); low.push_back({i, yl[j] * a});
                high.push_back({w(i,j), a}); high.push_back({i, -yu[j] * a});
            }
            low.push_back({n+j, source.b});
            high.push_back({n+j, -source.b});
            builder.append(low, source.b * yl[j], infinity);
            builder.append(high, -source.b * yu[j], infinity);
        }
    }
    for (const LinCons &source : linear_cons_P)
    {
        for (int i = 0; i < n; ++i)
        {
            std::vector<std::pair<int, double> > low, high;
            for (size_t k = 0; k < source.index.size(); ++k)
            {
                const int j = source.index[k] - 1;
                const double a = source.coefficient[k];
                low.push_back({w(i,j), -a}); low.push_back({n+j, xl[i] * a});
                high.push_back({w(i,j), a}); high.push_back({n+j, -xu[i] * a});
            }
            low.push_back({i, source.b});
            high.push_back({i, -source.b});
            builder.append(low, source.b * xl[i], infinity);
            builder.append(high, -source.b * xu[i], infinity);
        }
    }

    for (int block = 0; block < 2; ++block)
    {
        const int begin = block == 0 ? 0 : n_alt1;
        const int end = block == 0 ? n_alt1 : n;
        std::vector<std::pair<int, double> > normalization;
        for (int j = begin; j < end; ++j) normalization.push_back({n+j, 1.0});
        builder.append(normalization, 1.0, 1.0);
        for (int i = 0; i < n; ++i)
        {
            std::vector<std::pair<int, double> > product{{i, -1.0}};
            for (int j = begin; j < end; ++j) product.push_back({w(i,j), 1.0});
            builder.append(product, 0.0, 0.0);
        }
    }

    LinPara->num_row = static_cast<int>(builder.row);
    LinPara->num_col = LinPara->n_total;
    LinPara->num_nz = static_cast<int>(builder.nonzero);
    LinPara->offset = 0.0;
    LinPara->a_start_column = static_cast<int>(builder.row + 1);
    LinPara->a_index_column = static_cast<int>(builder.nonzero);
}

void putto_global_x_solution(gsl_vector_int *fixed_variables, gsl_vector *whole_solution, gsl_vector *x_solution)
{
    size_t i, j;

    if (gsl_vector_int_get(fixed_variables, 0) > 0)
    {
        j = 0;

        for (i = 1; i < fixed_variables->size; i++)
        {
            if (gsl_vector_int_get(fixed_variables, i) == 0)
            {
                gsl_vector_set(x_solution, j, gsl_vector_get(whole_solution, i - 1));
                j++;
            }
        }
    }
    else
    {
        gsl_vector_memcpy(x_solution, whole_solution);
    }
}

void getfrom_x_solution(gsl_vector_int *fixed_variables, gsl_vector *whole_solution, gsl_vector *xl, gsl_vector *x_solution)
{
    size_t i, j;

    if (gsl_vector_int_get(fixed_variables, 0) > 0)
    {
        j = 0;

        for (i = 1; i <= fixed_variables->size - 1; i++)
        {
            if (gsl_vector_int_get(fixed_variables, i) == 1)
            {
                gsl_vector_set(whole_solution, i - 1, gsl_vector_get(xl, i - 1));
            }
            else
            {
                gsl_vector_set(whole_solution, i - 1, gsl_vector_get(x_solution, j));
                j++;
            }
        }
    }
    else
    {
        gsl_vector_memcpy(whole_solution, x_solution);
    }
}

void construct_initial_solution(LinearizationParams *LinPara,
                                gsl_vector *x_solution,
                                gsl_vector *y_solution,
                                const int n_variable,
                                gsl_vector *initial_solution)
{
    int i, j, k;

    gsl_vector_view temp_v;

    temp_v = gsl_vector_subvector(initial_solution, 0, n_variable);
    gsl_vector_memcpy(&temp_v.vector, x_solution);

    temp_v = gsl_vector_subvector(initial_solution, n_variable, n_variable);
    gsl_vector_memcpy(&temp_v.vector, y_solution);

    k = LinPara->n_XY;

    for (i = 0; i < n_variable; i++)
    {
        for (j = 0; j < n_variable; j++)
        {
            gsl_vector_set(initial_solution, k, gsl_vector_get(x_solution, i) * gsl_vector_get(y_solution, j));
            k++;
        }
    }
    if (LinPara->n_YY < LinPara->n_total)
    {
        k = LinPara->n_XX;
        for (i = 0; i < n_variable; ++i)
            for (j = i; j < n_variable; ++j)
                gsl_vector_set(initial_solution, k++,
                    gsl_vector_get(x_solution, i) * gsl_vector_get(x_solution, j));
        k = LinPara->n_YY;
        for (i = 0; i < n_variable; ++i)
            for (j = i; j < n_variable; ++j)
                gsl_vector_set(initial_solution, k++,
                    gsl_vector_get(y_solution, i) * gsl_vector_get(y_solution, j));
    }
}

double solve_rlt_lp(LinearizationParams *p, gsl_vector *initial_solution,
                    bool use_initial, bool outputFlag,
                    RltCrossLevelState *cross_level_state,
                    DecideItDirectionStats *stats)
{
    if (selected_rlt_backend() != "highs")
        throw std::invalid_argument("This build requires the HiGHS RLT backend");

    // The server's packaged libhighs can corrupt its process-global internal
    // state when two independent instances enter passModel/run concurrently.
    // Serialize only the HiGHS-backed portion; Pivot-Copy and all geometry
    // above remain concurrent between the two objective directions.
    std::lock_guard<std::recursive_mutex> highs_lock(highs_api_mutex());

    if (!p->persistent) p->persistent.reset(createHighsSolver());
    HighsSolver* solver=p->persistent.get();
    if (!solver->loaded) {
        setLP(solver,p->num_col,p->num_row,p->a_start->data,p->a_index->data,p->a_value->data,
              p->col_cost->data,p->col_lower->data,p->col_upper->data,p->row_lower->data,p->row_upper->data);
        // A previous contraction level has the same variables and often the
        // same row layout.  Its optimal basis is a much richer warm start
        // than a primal point.  HiGHS refactorizes the supplied statuses for
        // the new coefficients and bounds; any incompatible dimensions fall
        // back to the feasible AMC rank-one point.
        bool reused_cross_level_basis = false;
        const char *cross_level_mode =
            std::getenv("DECIDEIT_CROSS_LEVEL_RLT_START");
        const bool basis_requested =
            (cross_level_mode && std::string(cross_level_mode) == "basis") ||
            std::getenv("DECIDEIT_CROSS_DIRECTION_RLT_BASIS");
        if (cross_level_state && cross_level_state->basis_valid &&
            cross_level_state->rows == p->num_row &&
            cross_level_state->columns == p->num_col &&
            basis_requested &&
            !std::getenv("DECIDEIT_COLD_RLT"))
        {
            const HighsStatus basis_status = solver->highs.setBasis(
                cross_level_state->basis, "preceding contraction");
            reused_cross_level_basis = basis_status != HighsStatus::kError;
            if (reused_cross_level_basis && stats)
                stats->cross_level_basis_uses++;
        }
        bool reused_cross_level_primal = false;
        if (!reused_cross_level_basis && cross_level_state &&
            cross_level_state->valid &&
            cross_level_state->columns == p->num_col &&
            static_cast<int>(cross_level_state->solution.size()) == p->num_col &&
            static_cast<int>(cross_level_state->col_lower.size()) == p->num_col &&
            static_cast<int>(cross_level_state->col_upper.size()) == p->num_col &&
            cross_level_mode && std::string(cross_level_mode) == "transport" &&
            !std::getenv("DECIDEIT_COLD_RLT"))
        {
            // Transport the preceding lifted solution through the same
            // homothety used by contraction.  For x'=r*x+a and y'=r*y+b,
            // xy is mapped exactly as r^2*xy+r*b*x+r*a*y+a*b.  This retains
            // the useful non-rank-one RLT information instead of discarding
            // it when only P is transported between levels.
            std::vector<double> candidate = cross_level_state->solution;
            const int n = p->n_XY / 2;
            double contraction_ratio = -1.0;
            for (int column = 0; column < 2 * n; ++column)
            {
                const double old_width =
                    cross_level_state->col_upper[column] -
                    cross_level_state->col_lower[column];
                const double new_width = gsl_vector_get(p->col_upper, column) -
                    gsl_vector_get(p->col_lower, column);
                if (old_width > 1e-12)
                {
                    contraction_ratio = new_width / old_width;
                    break;
                }
            }
            contraction_ratio = std::fmax(0.0,
                std::fmin(1.0, contraction_ratio < 0.0 ? 1.0 :
                          contraction_ratio));
            std::vector<double> offset(2 * n, 0.0);
            for (int column = 0; column < 2 * n; ++column)
            {
                offset[column] = gsl_vector_get(p->col_lower, column) -
                    contraction_ratio * cross_level_state->col_lower[column];
                candidate[column] = contraction_ratio *
                    cross_level_state->solution[column] + offset[column];
            }
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    const int w = 2 * n + i * n + j;
                    candidate[w] = contraction_ratio * contraction_ratio *
                        cross_level_state->solution[w] +
                        contraction_ratio * offset[n + j] *
                        cross_level_state->solution[i] +
                        contraction_ratio * offset[i] *
                        cross_level_state->solution[n + j] +
                        offset[i] * offset[n + j];
                }

            // The AMC rank-one point is feasible for the current relaxation.
            // Guard against changed active rows and rounding by clipping the
            // segment to the transported point at the first current row or
            // bound intersection.  This remains an O(nnz) calculation.
            double scale = 1.0;
            auto restrict_scale = [&](double direction, double slack)
            {
                if (direction > 1e-14)
                    scale = std::fmin(scale, std::fmax(0.0, slack / direction));
            };
            for (int column = 0; column < p->num_col; ++column)
            {
                const double anchor = gsl_vector_get(initial_solution, column);
                const double direction = candidate[column] - anchor;
                restrict_scale(direction,
                    gsl_vector_get(p->col_upper, column) - anchor);
                restrict_scale(-direction,
                    anchor - gsl_vector_get(p->col_lower, column));
            }
            for (int row = 0; row < p->num_row; ++row)
            {
                double anchor_lhs = 0.0;
                double direction_lhs = 0.0;
                const int begin = gsl_vector_int_get(p->a_start, row);
                const int end = gsl_vector_int_get(p->a_start, row + 1);
                for (int entry = begin; entry < end; ++entry)
                {
                    const int column = gsl_vector_int_get(p->a_index, entry);
                    const double coefficient = gsl_vector_get(p->a_value, entry);
                    const double anchor = gsl_vector_get(initial_solution, column);
                    anchor_lhs += coefficient * anchor;
                    direction_lhs += coefficient *
                        (candidate[column] - anchor);
                }
                restrict_scale(direction_lhs,
                    gsl_vector_get(p->row_upper, row) - anchor_lhs);
                restrict_scale(-direction_lhs,
                    anchor_lhs - gsl_vector_get(p->row_lower, row));
            }
            scale = std::fmax(0.0, std::fmin(1.0, scale));
            if (scale < 1.0) scale *= 1.0 - 1e-10;
            for (int column = 0; column < p->num_col; ++column)
            {
                const double anchor = gsl_vector_get(initial_solution, column);
                gsl_vector_set(initial_solution, column,
                    anchor + scale *
                    (candidate[column] - anchor));
            }
            setWarmStart(solver, initial_solution->data, p->num_col);
            reused_cross_level_primal = true;
            if (stats)
            {
                stats->cross_level_primal_uses++;
                stats->cross_level_primal_scale_sum += scale;
            }
        }
        // Use the current feasible AMC point at activation when no compatible
        // cross-level basis is available, then retain the resulting basis for
        // all row additions within this DBLP.
        if (!reused_cross_level_basis && !reused_cross_level_primal &&
            use_initial &&
            !std::getenv("DECIDEIT_COLD_RLT"))
            setWarmStart(solver,initial_solution->data,p->num_col);
        solver->cold_solves++;
    } else {
        const int begin=solver->loaded_rows, count=p->num_row-begin;
        if(count<0) throw std::runtime_error("RLT rows were removed without resetting model");
        if(count) {
            const int first=gsl_vector_int_get(p->a_start,begin);
            const int last=gsl_vector_int_get(p->a_start,p->num_row);
            std::vector<HighsInt> starts(count+1);
            for(int i=0;i<=count;i++) starts[i]=gsl_vector_int_get(p->a_start,begin+i)-first;
            highs_require(solver->highs.addRows(count,p->row_lower->data+begin,p->row_upper->data+begin,
                last-first,starts.data(),p->a_index->data+first,p->a_value->data+first),"RLT addRows");
        }
        solver->loaded_rows=p->num_row;
        if(std::getenv("DECIDEIT_COLD_RLT")) {
            highs_require(solver->highs.clearSolver(),"RLT requested cold solve");
            solver->cold_solves++;
        } else solver->warm_solves++;
        solver->solved=false;
    }
    solve_with_highs(solver,outputFlag);
    getSolution(solver,initial_solution->data,p->num_col);
    if (cross_level_state)
    {
        const HighsBasis &basis = solver->highs.getBasis();
        cross_level_state->valid = true;
        cross_level_state->basis_valid = basis.valid;
        cross_level_state->rows = p->num_row;
        cross_level_state->columns = p->num_col;
        cross_level_state->basis = basis;
        cross_level_state->solution.assign(
            initial_solution->data, initial_solution->data + p->num_col);
        cross_level_state->col_lower.assign(
            p->col_lower->data, p->col_lower->data + p->num_col);
        cross_level_state->col_upper.assign(
            p->col_upper->data, p->col_upper->data + p->num_col);
    }
    return getObjective(solver);
}

// Add every row accumulated in LinearizationParams to the persistent HiGHS
// model without optimizing it.  The controller calls this only after checking
// that the preceding certified optimum satisfies every new row.  That point
// is then still feasible and still attains the old lower bound, so it remains
// optimal for the restricted LP.  Keeping the model synchronized here is
// essential for a later warm solve or branch-and-RLT hand-off.
void append_rlt_rows_without_solving(LinearizationParams *p)
{
    if (selected_rlt_backend() != "highs" || !p->persistent ||
        !p->persistent->loaded)
        throw std::runtime_error(
            "Cannot synchronize unsolved or non-HiGHS RLT model");
    std::lock_guard<std::recursive_mutex> highs_lock(highs_api_mutex());
    HighsSolver *solver = p->persistent.get();
    const int begin = solver->loaded_rows;
    const int count = p->num_row - begin;
    if (count < 0)
        throw std::runtime_error(
            "RLT rows were removed without resetting model");
    if (count)
    {
        const int first = gsl_vector_int_get(p->a_start, begin);
        const int last = gsl_vector_int_get(p->a_start, p->num_row);
        std::vector<HighsInt> starts(count + 1);
        for (int i = 0; i <= count; ++i)
            starts[i] = gsl_vector_int_get(p->a_start, begin + i) - first;
        highs_require(solver->highs.addRows(
            count, p->row_lower->data + begin,
            p->row_upper->data + begin, last - first, starts.data(),
            p->a_index->data + first, p->a_value->data + first),
            "RLT add satisfied rows");
        solver->loaded_rows = p->num_row;
        // HiGHS may invalidate its solution status after addRows.  The
        // controller retains the independently checked primal and objective;
        // a future solve will reuse the expanded basis but must call run().
        solver->solved = false;
    }
}

void recover_cut(gsl_vector *onecut, gsl_vector *normal_vector, gsl_vector *, gsl_vector_int *fixed_variables)
{
    gsl_vector_set_zero(onecut);
    gsl_vector_set(onecut,0,gsl_vector_get(normal_vector,0));
    size_t source=1;
    for(size_t i=1;i<fixed_variables->size;i++)
        if(!gsl_vector_int_get(fixed_variables,i))
            gsl_vector_set(onecut,i,gsl_vector_get(normal_vector,source++));
}

/*
    record_results() saves the final results

    inputs
        solution: the structure containing the global solution information;
        P_point, V_point: final P and V results starting with appropriate indices;
        Ai_begin, Ai_end: positiion for alternavtive i;
        Aj_begin, Aj_end: positiion for alternavtive j;
        n_cuts: the # of cuts;
        solverFlag: to indicate the reason to terminate the solver;
        localFlag: to indicate whether the solver is terminated during the local phase;
        globalFlag: to indicate whether the solver is terminated during the global phase;
        fixed_varibles: the # and indices of fixed variables;
        instance: the structure containing the original data.

    outputs
        P_point, V_point: final P and V results starting with appropriate indices;
        P_point[0]: the objective value of minimizing problem;
        V_point[0]: the # of cuts, a negative number cannot guarantee the global optimality.
*/

void record_results(struct global_solution *solution,
                    d_row P_point,
                    d_row V_point,
                    int Ai_begin,
                    int Ai_end,
                    int Aj_begin,
                    int Aj_end,
                    int n_cuts,
                    int solverFlag,
                    bool localFlag,
                    bool globalFlag,
                    gsl_vector_int *fixed_variables,
                    struct dblp *instance,
                    gsl_vector *whole_solution)
{
    size_t i, j;

    /* construct the whole_solution if some variables are fixed */
    if (gsl_vector_int_get(fixed_variables, 0) > 0)
    {
        j = 0;

        for (i = 1; i <= fixed_variables->size - 1; i++)
        {
            if (gsl_vector_int_get(fixed_variables, i) == 1)
            {
                gsl_vector_set(whole_solution, i - 1, gsl_vector_get(instance->xl, i - 1));
            }
            else
            {
                gsl_vector_set(whole_solution, i - 1, gsl_vector_get(solution->x_solution, j));
                j++;
            }
        }
    }
    else
    {
        gsl_vector_memcpy(whole_solution, solution->x_solution);
    }

    P_point[0] = solution->alpha;

    if (localFlag)
    {
        switch (solverFlag)
        {
        case OPTIMAL:
            V_point[0] = n_cuts;
            break;

        case GSL_LU:
            V_point[0] = -n_cuts;
            break;

        case DEGENERACY: // neighbor or visited-basis capacity reached
        case NO_T_COL: // the starting basis has no transition column
            V_point[0] = -n_cuts;
            break;

        default:
            break;
        }
    }
    else if (globalFlag)
    {
        switch (solverFlag)
        {
        case GSL_LU:
            V_point[0] = -n_cuts;
            break;

        case SOLUTION_LP:
            V_point[0] = -n_cuts;
            break;

        case CUT_SEPARATION_FAILED:
            V_point[0] = -n_cuts;
            break;

        case ALL_NEGATIVE: // due to all negative cutting points
            V_point[0] = n_cuts;
            break;

        default:
            break;
        }
    }
    else
    {
        switch (solverFlag)
        {
        case INFEASIBLE: // due to an empty feasible region
            V_point[0] = n_cuts;
            break;

        case CUT_LIMIT: // due to the pre-sepecified maximal # of cuts
            V_point[0] = -n_cuts;
            break;

        default:
            break;
        }
    }

    j = 0;

    if ((!localFlag) || (localFlag && solverFlag == OPTIMAL) || (localFlag && (solverFlag == DEGENERACY || solverFlag == NO_T_COL)) || (localFlag && solverFlag == GSL_LU && n_cuts > 0))
    {
        for (i = Ai_begin; i <= (size_t)Ai_end; i++)
        {
            V_point[i] = gsl_vector_get(whole_solution, j);
            P_point[i] = gsl_vector_get(solution->y_solution, j++);
        }

        for (i = Aj_begin; i <= (size_t)Aj_end; i++)
        {
            V_point[i] = gsl_vector_get(whole_solution, j);
            P_point[i] = gsl_vector_get(solution->y_solution, j++);
        }
    }
}

bool check_optimality(gsl_vector *solution, int n_variable)
{
    int i, j, n_XY = 2 * n_variable, position;
    bool Flag = TRUE;

    for (i = 0; i < n_variable; i++)
    {
        position = n_XY + i * n_variable;

        for (j = 0; j < n_variable; j++)
        {
            if (!std::isfinite(gsl_vector_get(solution,i)) || !std::isfinite(gsl_vector_get(solution,n_variable+j)) || !std::isfinite(gsl_vector_get(solution,position)) || fabs(gsl_vector_get(solution, i) * gsl_vector_get(solution, n_variable + j) - gsl_vector_get(solution, position++)) > 10 * DBLP_EPSILON)
            {
                Flag = FALSE;
                break;
            }
        }
    }

    return Flag;
}

void linearization_onecut2go(gsl_vector *onecut,
                             LinearizationParams *LinPara,
                             d_row P_lobo, d_row P_upbo,
                             int Ai_begin, int Aj_begin, int n_alt1,
                             const std::vector<LinCons> *opposite_rows,
                             const double *separation_values,
                             size_t separation_size,
                             const std::vector<int> *forced_cross_rows)
{
    size_t i, j, pos, loc;
    const size_t n_variable = onecut->size - 1;
    if (binary_environment_option("DECIDEIT_RLT_NORMALIZE_CUTS"))
    {
        double scale = 0.0;
        for (i = 0; i <= n_variable; ++i)
            scale = (std::max)(scale,
                               std::fabs(gsl_vector_get(onecut, i)));
        if (std::isfinite(scale) && scale > 0.0 &&
            std::fabs(scale - 1.0) > 1e-15)
            gsl_vector_scale(onecut, 1.0 / scale);
    }
    size_t cut_nonzeros = 0, cross_nonzeros = 0;
    for (i = 1; i <= n_variable; ++i)
        if (std::fabs(gsl_vector_get(onecut, i)) > 1e-15)
            ++cut_nonzeros;
    const CrossRowRltMode cross_mode = selected_cross_row_rlt_mode();
    const bool add_cross_rows =
        opposite_rows && (cross_mode != CrossRowRltMode::disabled ||
                          forced_cross_rows);
    if (add_cross_rows)
        for (const LinCons &row : *opposite_rows)
            cross_nonzeros += cut_nonzeros * row.index.size() +
                              cut_nonzeros + row.index.size();
    ensure_rlt_capacity(
        LinPara,
        LinPara->num_row + 2*n_variable + 2 +
            (add_cross_rows ? opposite_rows->size() : 0),
        LinPara->num_nz + 4*n_variable*n_variable + 3*n_variable +
            cross_nonzeros);
    const size_t nX = n_variable;
    const size_t nY = n_variable;

    size_t pointer = LinPara->a_start_column;
    size_t column = LinPara->a_index_column;
    size_t row = LinPara->num_row;

    size_t n_XY = (size_t)LinPara->n_XY;

    gsl_vector_int *a_index = LinPara->a_index;
    gsl_vector *a_value = LinPara->a_value;
    gsl_vector_int *a_start = LinPara->a_start;
    gsl_vector *row_lower = LinPara->row_lower;
    gsl_vector *row_upper = LinPara->row_upper;

    /* (2.2a) */
    for (i = 1; i <= nX; i++)
    {
        gsl_vector_int_set(a_index, column, i - 1);
        gsl_vector_set(a_value, column++, gsl_vector_get(onecut, i));
    }

    gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + nX);
    pointer++;

    gsl_vector_set(row_lower, row, -1.0e30);
    gsl_vector_set(row_upper, row, gsl_vector_get(onecut, 0));
    row++;

    /* (2.2d) */
    for (j = 0; j < nX; j++)
    {
        loc = location(j, Ai_begin, Aj_begin, n_alt1);

        // -[A(1r) * x] * y(j)
        for (i = 0; i < nX; i++)
        {
            pos = xy_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, -gsl_vector_get(onecut, i + 1));
        }

        // L(j) * [A(1r) * x]
        for (i = 0; i < nX; i++)
        {
            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, P_lobo[loc] * gsl_vector_get(onecut, i + 1));
        }

        // b(1r) * y(j)
        gsl_vector_int_set(a_index, column, nX + j);
        gsl_vector_set(a_value, column++, gsl_vector_get(onecut, 0));
        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + nX + nX + 1);

        // std::cout << gsl_vector_int_get(a_start, pointer) << std::endl;

        pointer++;

        gsl_vector_set(row_lower, row, P_lobo[loc] * gsl_vector_get(onecut, 0));
        gsl_vector_set(row_upper, row++, 1.0e30);
    }

    /* (2.2e) */
    for (j = 0; j < nY; ++j)
    {
        loc = location(j, Ai_begin, Aj_begin, n_alt1);

        for (i = 0; i < nX; i++)
        {
            pos = xy_index0(i, j, nY);
            gsl_vector_int_set(a_index, column, n_XY + pos);
            gsl_vector_set(a_value, column++, gsl_vector_get(onecut, i + 1));
        }

        for (i = 0; i < nX; i++)
        {
            gsl_vector_int_set(a_index, column, i);
            gsl_vector_set(a_value, column++, -P_upbo[loc] * gsl_vector_get(onecut, i + 1));
        }

        gsl_vector_int_set(a_index, column, nX + j);
        gsl_vector_set(a_value, column++, -gsl_vector_get(onecut, 0));

        gsl_vector_int_set(a_start, pointer, gsl_vector_int_get(a_start, pointer - 1) + nX + nX + 1);
        pointer++;

        gsl_vector_set(row_lower, row, -P_upbo[loc] * gsl_vector_get(onecut, 0));
        gsl_vector_set(row_upper, row++, 1.0e30);
    }

    if (add_cross_rows)
    {
        std::vector<std::pair<int, double> > cut_terms;
        cut_terms.reserve(cut_nonzeros);
        for (i = 1; i <= n_variable; ++i)
            if (std::fabs(gsl_vector_get(onecut, i)) > 1e-15)
                cut_terms.push_back(
                    {static_cast<int>(i - 1), gsl_vector_get(onecut, i)});
        RltBuilder builder{LinPara};
        builder.row = row;
        builder.nonzero = column;
        std::vector<int> selected;
        std::vector<bool> forced(opposite_rows->size(), false);
        if (forced_cross_rows)
            for (int index : *forced_cross_rows)
            {
                if (index < 0 || index >=
                        static_cast<int>(opposite_rows->size()) ||
                    forced[index])
                    throw std::runtime_error(
                        "Invalid branch RLT cross-row replay index");
                forced[index] = true;
            }
        for (size_t index = 0; index < opposite_rows->size(); ++index)
        {
            const bool replay = forced_cross_rows && forced[index];
            if (forced_cross_rows && !replay) continue;
            const bool appended = append_cross_row_product(
                builder, cut_terms, gsl_vector_get(onecut, 0),
                (*opposite_rows)[index], static_cast<int>(n_variable),
                separation_values, separation_size,
                !forced_cross_rows &&
                    cross_mode == CrossRowRltMode::violated);
            if (appended) selected.push_back(static_cast<int>(index));
        }
        LinPara->cut_cross_rows.push_back(selected);
        row = builder.row;
        column = builder.nonzero;
        pointer = row + 1;
    }

    LinPara->num_row = row;

    LinPara->num_nz = gsl_vector_int_get(LinPara->a_start, pointer - 1);
    LinPara->a_start_column = pointer;
    LinPara->a_index_column = column;

    a_index = NULL;
    a_value = NULL;
    a_start = NULL;
    row_lower = NULL;
    row_upper = NULL;
}

bool identify_basics(
    LpModel *lp_X,
    gsl_vector_int *basic_index,
    gsl_vector_int *nonbasic_index,
    gsl_vector *whole_solution,
    std::vector<int> &basic,
    std::vector<int> &nonbasic)
{
    const int m = get_Nrows(lp_X), n = get_Ncolumns(lp_X);
    std::vector<REAL> constr(m);
    std::vector<REAL> var(n);
    REAL value_r;
    std::vector<int> bascolumn(1 + m);
    int i, value_i;
    std::vector<int> total_index(1 + m + n);

    gsl_vector_view src_v, dest_v;

    get_basis(lp_X, bascolumn.data(), FALSE);

    basic[0] = 0;

    for (i = 1; i <= m; i++)
    {
        value_i = abs(bascolumn[i]);

        if (value_i <= m)
        {
            basic[i] = value_i + n;
        }
        else
        {
            basic[i] = value_i - m;
        }
    }

    // qsort(basic.data(), 1 + m, sizeof(int), compare_ints);
    std::sort(basic.data(), basic.data() + 1 + m);

    basic[0] = m;
    // show_array_int(basic);

    // src_v = gsl_vector_int_view_array(swap_array, 1 + m);
    // gsl_vector_int_memcpy(basic_index, &src_v.vector);
    // gsl_vector_int_set(basic_index, 0, m);

    get_variables(lp_X, var.data());

    src_v = gsl_vector_view_array(var.data(), n);
    dest_v = gsl_vector_subvector(whole_solution, 1, n);
    gsl_vector_memcpy(&dest_v.vector, &src_v.vector);

    get_constraints(lp_X, constr.data());

    for (i = 0; i < m; i++)
    {
        value_r = get_rh(lp_X, i + 1) - constr[i];

        if (fabs(value_r) <= DBLP_EPSILON)
        {
            gsl_vector_set(whole_solution, n + i + 1, 0);
        }
        else
        {
            gsl_vector_set(whole_solution, n + i + 1, value_r);
        }
    }

    gsl_vector_set(whole_solution, 0, m + n);

    gsl_vector_int_set(basic_index, 0, 0);
    gsl_vector_int_set(nonbasic_index, 0, 0);

    for (i = 1; i <= gsl_vector_get(whole_solution, 0); i++)
    {
        if (gsl_vector_get(whole_solution, i) > DBLP_EPSILON)
        {
            gsl_vector_int_set(basic_index, 0, gsl_vector_int_get(basic_index, 0) + 1);
            gsl_vector_int_set(basic_index, gsl_vector_int_get(basic_index, 0), i);
        }
        else
        {
            gsl_vector_int_set(nonbasic_index, 0, gsl_vector_int_get(nonbasic_index, 0) + 1);
            gsl_vector_int_set(nonbasic_index, gsl_vector_int_get(nonbasic_index, 0), i);
        }
    }

    if (gsl_vector_int_get(basic_index, 0) == m)
    {
        return FALSE;
    }
    else
    {
        for (i = 1; i <= m + n; i++)
        {
            total_index[i] = i;
        }

        total_index[0] = m + n;

        setDifference(total_index, basic, nonbasic);

        return TRUE;
    }

    // show_gsl_vector_int(basic_index);
    // show_gsl_vector_int(nonbasic_index);
    // printf("\n");
}

/*
    solution_from_tableau_new() constructs the solution to regular variables.

    inputs
        solution: to save the constructed solution;
        n_variable: the # of regular variables;
        basic_index: the indices of basic variables;
        rhs: the right-hand side.

    output
        solution: the constructed solution.
*/

void solution_from_tableau_new(
    gsl_vector *solution,
    int n_variable,
    gsl_vector_int *basic_index,
    gsl_vector *rhs)
{
    int j;

    gsl_vector_set_zero(solution);

    for (j = 1; j <= gsl_vector_int_get(basic_index, 0); j++)
    {
        if (gsl_vector_int_get(basic_index, j) <= n_variable) // it is a regular variable
        {
            gsl_vector_set(solution, gsl_vector_int_get(basic_index, j) - 1, gsl_vector_get(rhs, j - 1));
        }
    }
}

/*
    adj_vertices_normal() locates all adjacent vertices of a non-degenerate local solution in X.

    inputs
        tableau: the tableau containing the transformed effective linear constraints, intervals (0 <= x <= ub - lb) and slack variables;
        basic_index and nonbasic_index: the indices of identified basic and non-basic variables;
        noc: the # of effective linear constraints in X;
        adj_vertices: the matrix containing all vertices adjacent to x;
        lb: lower bounds for the effective variables in X.

    outputs
        true for a successful identification, false otherwise;
        adj_vertices: the constructed matrix containing the identified vertices (in columns) adjacent to x.
*/

bool adj_vertices_normal(gsl_matrix *tableau,
                         gsl_vector_int *basic_index,
                         gsl_vector_int *nonbasic_index,
                         struct local_solution *local_amc,
                         struct dataset *dataset_in)
{
    const int n_variable = tableau->size2 - tableau->size1 - 1; // the # of regular variables in X
    int j;
    size_t i, row, column;
    REAL min_value, temp_value, theta;

    gsl_matrix *inverse = managed::matrix_calloc(tableau->size1, tableau->size1);
    gsl_matrix *temp_inverse = managed::matrix_calloc(tableau->size1, tableau->size1);
    gsl_permutation *p = managed::permutation_calloc(tableau->size1);

    gsl_vector_view src_v;

    /* find the inverse of the current basis */
    for (j = 1; j <= gsl_vector_int_get(basic_index, 0); j++)
    {
        src_v = gsl_matrix_column(tableau, gsl_vector_int_get(basic_index, j));
        gsl_matrix_set_col(inverse, j - 1, &src_v.vector);
    }

    // if (!condition_number(inverse, 1e12)) // if satisfied, the matrix is close to singular, thus not invertible
    // {
    //     managed::matrix_free(inverse);
    //     managed::permutation_free(p);

    //     return (false);
    // }

    gsl_matrix_memcpy(temp_inverse, inverse);

    if (safe_inverse(temp_inverse, inverse, 1e-12) != 0)
    {
        managed::matrix_free(inverse);
        managed::matrix_free(temp_inverse);
        managed::permutation_free(p);

        return (false);
    }

    // if (gsl_linalg_LU_decomp(inverse, p, &s) != GSL_SUCCESS)
    // {
    //     managed::matrix_free(inverse);
    //     managed::permutation_free(p);

    //     return (false);
    // }

    // if (gsl_linalg_LU_invert(inverse, p, inverse) != GSL_SUCCESS)
    // {
    //     managed::matrix_free(inverse);
    //     managed::permutation_free(p);

    //     return (false);
    // }

    gsl_matrix *tableau_new = managed::matrix_calloc(tableau->size1, tableau->size2);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, inverse, tableau, 0, tableau_new); // inverse * tableau

    // gsl_vector_view temp = gsl_matrix_column(tableau_new, 0);

    /* figure out all vertices adjacent to x by pivoting each non-basic variable*/
    gsl_vector *new_solution = managed::vector_calloc(n_variable);
    gsl_vector_int *basic_copy = managed::vector_int_calloc(basic_index->size);
    gsl_vector *a = managed::vector_calloc(tableau_new->size1);
    gsl_vector *c = managed::vector_calloc(tableau_new->size1);
    bool valid_neighbors = true;

    // gsl_matrix_get_col(a, tableau_new, 0);
    // solution_from_tableau_new(new_solution, n_variable, basic_index, a);
    // gsl_vector_add(new_solution, lb); // transform the derived vertex back to the original space
    // show_gsl_vector(new_solution);

    for (j = 1; j <= gsl_vector_int_get(nonbasic_index, 0); j++)
    {
        column = gsl_vector_int_get(nonbasic_index, j);
        min_value = std::numeric_limits<double>::infinity();
        row = tableau_new->size1;

        for (i = 0; i < tableau_new->size1; i++)
        {
            if (gsl_matrix_get(tableau_new, i, column) > DBLP_EPSILON)
            {
                temp_value = gsl_matrix_get(tableau_new, i, 0) / gsl_matrix_get(tableau_new, i, column);

                if (std::isfinite(temp_value) && temp_value < min_value)
                {
                    min_value = temp_value;
                    row = i;
                }
            }
        }

        // A bounded X should offer a leaving row for every entering column.
        // If numerical filtering removes all pivots, do not manufacture a
        // neighbor by dividing the arbitrary first row by zero.
        if (row == tableau_new->size1)
        {
            valid_neighbors = false;
            break;
        }

        gsl_vector_int_memcpy(basic_copy, basic_index);
        gsl_vector_int_set(basic_copy, row + 1, column); // change the basic indices

        theta = gsl_matrix_get(tableau_new, row, 0) / gsl_matrix_get(tableau_new, row, column);
        src_v = gsl_matrix_column(tableau_new, 0);
        gsl_vector_memcpy(a, &src_v.vector);
        src_v = gsl_matrix_column(tableau_new, column);
        gsl_vector_memcpy(c, &src_v.vector);
        gsl_vector_scale(c, theta);
        gsl_vector_sub(a, c); // current values of basic variables
        gsl_vector_set(a, row, theta);

        solution_from_tableau_new(new_solution, n_variable, basic_copy, a);
        gsl_vector_add(new_solution, dataset_in->xl);
        gsl_matrix_set_col(local_amc->adj_vertices, j - 1, new_solution);
    }
    
    managed::matrix_free(inverse);
    managed::matrix_free(temp_inverse);
    managed::permutation_free(p);
    managed::matrix_free(tableau_new);
    managed::vector_free(new_solution);
    managed::vector_int_free(basic_copy);
    managed::vector_free(a);
    managed::vector_free(c);

    return valid_neighbors;
}

int neighborhood(int noc, struct dataset *dataset_in, LpModel *lp_X_std, struct local_solution *local_amc)
{
    const int n_variable_x = get_Ncolumns(lp_X_std), n_constraint_x = noc;
    const int n_constraint_tableau = n_constraint_x + n_variable_x;
    const int n_variable_tableau = n_constraint_tableau + n_variable_x;
    int i, num_neighbors;
    bool degeneracyFlag;
    local_amc->degenerate_vertex = false;
    local_amc->neighbor_states.clear();
    local_amc->cone_geometry = ConeRayGeometry{};
    std::vector<int> basic(1 + get_Nrows(lp_X_std));
    std::vector<int> nonbasic(1 + get_Ncolumns(lp_X_std));

    gsl_matrix *adj_vertices = managed::matrix_calloc(
        n_variable_x, adjacency_capacity(n_variable_x));
    gsl_matrix *tableau = managed::matrix_calloc(n_constraint_tableau, 1 + n_variable_tableau);

    gsl_vector *vertex_x = managed::vector_calloc(n_variable_x);
    gsl_vector *coef_x_std = managed::vector_calloc(1 + get_Ncolumns(lp_X_std));
    gsl_vector *x_solution_std = managed::vector_calloc(get_Ncolumns(lp_X_std));
    gsl_vector *x_solution = managed::vector_calloc(n_variable_x);
    gsl_vector *whole_solution = managed::vector_calloc(1 + n_constraint_x + 2 * n_variable_x);

    gsl_vector_int *basic_index = managed::vector_int_calloc(1 + n_constraint_x + n_variable_x);
    gsl_vector_int *nonbasic_index = managed::vector_int_calloc(1 + n_constraint_x + 2 * n_variable_x);

    gsl_vector_view src_v, dest_v;    
    gsl_vector_int_view src_v_int, dest_v_int;

    if (n_constraint_x == 0)
    {
        /*
            this is a very special case, which may occur in the first iteration with only box constraints.
            in particular, it is IMPOSSIBLE for x to be degenerate.
            for such a case, it is unncessary to utilize a simplex tableau to derive the neighbors of x.
            after the first iteration, a cut will definitely be incorporated into Valid_A so that n_constraint_x > 0.
        */

        for (i = 0; i < n_variable_x; i++)
        {
            gsl_vector_memcpy(vertex_x, local_amc->x_solution);

            if (fabs(gsl_vector_get(vertex_x, i) - gsl_vector_get(dataset_in->xl, i)) < DBLP_EPSILON)
            {
                gsl_vector_set(vertex_x, i, gsl_vector_get(dataset_in->xu, i));
            }
            else
            {
                gsl_vector_set(vertex_x, i, gsl_vector_get(dataset_in->xl, i));
            }

            gsl_matrix_set_col(local_amc->adj_vertices, i, vertex_x);
        }

        num_neighbors = n_variable_x;
    }
    else
    {
        degeneracyFlag = identify_basics(lp_X_std, basic_index, nonbasic_index, whole_solution, basic, nonbasic);
        local_amc->degenerate_vertex = degeneracyFlag;
        // show_gsl_vector_int(basic_index);
        // show_gsl_vector(whole_solution);

        tableau_construction(lp_X_std, tableau);
        
        // gsl_vector_view temp = gsl_matrix_column(tableau_copy, 0);

        if (!degeneracyFlag) // the solution is non-degenerate
        {
            /* get the adjacent vertices by pivoting the simplex tableau */
            if (!adj_vertices_normal(tableau, basic_index, nonbasic_index, local_amc, dataset_in))
            {
                num_neighbors = GSL_LU;
            }
            else
            {
                num_neighbors = n_variable_x;
            }
        }
        else // the solution is degenerate
        {
            src_v_int = gsl_vector_int_view_array(basic.data(), 1 + get_Nrows(lp_X_std));
            gsl_vector_int_memcpy(basic_index, &src_v_int.vector);
            src_v_int = gsl_vector_int_view_array(nonbasic.data(), 1 + get_Ncolumns(lp_X_std));
            dest_v_int = gsl_vector_int_subvector(nonbasic_index, 0, 1 + get_Ncolumns(lp_X_std));
            gsl_vector_int_memcpy(&dest_v_int.vector, &src_v_int.vector);

            const char *store_setting =
                std::getenv("DECIDEIT_STORE_NEIGHBOR_STATE");
            const bool store_states = store_setting &&
                std::strcmp(store_setting, "1") == 0;
            const char *boundary_setting =
                std::getenv("DECIDEIT_COMPLETE_BOUNDARY_DIAGNOSTIC");
            const char *boundary_cut_setting =
                std::getenv("DECIDEIT_COMPLETE_BOUNDARY_CUT");
            const bool boundary_diagnostic =
                (boundary_setting && std::strcmp(boundary_setting, "1") == 0) ||
                (boundary_cut_setting &&
                 std::strcmp(boundary_cut_setting, "1") == 0);
            num_neighbors = reverse_search(
                tableau, basic_index, nonbasic_index, adj_vertices,
                (store_states || boundary_diagnostic) ?
                    &local_amc->neighbor_states : nullptr);

            if (num_neighbors < 0 ||
                ((store_states || boundary_diagnostic) && static_cast<int>(
                    local_amc->neighbor_states.size()) != num_neighbors))
                local_amc->neighbor_states.clear();

            if (num_neighbors > 0)
                grow_neighbor_matrix(local_amc->adj_vertices,
                                     static_cast<size_t>(num_neighbors));
            for (i = 0; i < num_neighbors; i++)
            {
                dest_v = gsl_matrix_column(local_amc->adj_vertices, i);
                src_v = gsl_matrix_column(adj_vertices, i);
                gsl_vector_memcpy(&dest_v.vector, &src_v.vector);
                gsl_vector_add(&dest_v.vector, dataset_in->xl);
            }
        }
    }

    local_amc->n_adj_vertices = num_neighbors;
    
    managed::vector_free(vertex_x);
    managed::vector_free(x_solution_std);
    managed::vector_free(coef_x_std);
    managed::matrix_free(adj_vertices);
    managed::matrix_free(tableau);
    managed::vector_free(x_solution);
    managed::vector_int_free(basic_index);
    managed::vector_int_free(nonbasic_index);
    managed::vector_free(whole_solution);

    return num_neighbors;
}

int cpm_controller(struct dblp *instance,
                   struct local_solution *local_amc,
                   LpModel *lp_X_std,
                   LpModel *lp_Y,
                   struct dataset *dataset_in,
                   int *noc,
                   d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
                   struct global_solution *solution,
                   int CUTS_MAX,
                   d_row P_point, d_row V_point,
                   int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
                   gsl_vector_int *fixed_variables,
                   std::vector<LinCons> &LC_V, std::vector<LinCons> &LC_P,
                   const std::vector<LinCons> &original_P_rows,
                   int mX, int mY, LinearizationParams *LinPara,
                   bool minFlag,
                   const gsl_vector *initial_y,
                   RltCrossLevelState *cross_level_rlt_state)
{
    const int n_variable = get_Ncolumns(lp_X_std);
    const int n_alt_i = Ai_end - Ai_begin + 1, n_alt_j = Aj_end - Aj_begin + 1;
    const int actual_n_variables = n_alt_i + n_alt_j;

    bool loop_flag = true, localFlag, globalFlag;
    int exitFlag = OK, solverFlag, n_cuts = 0;
    double temp_obj = -std::numeric_limits<double>::infinity();
    // Tightened bounds are private to this objective direction and to the
    // residual region left by its cuts.  Never write them back to the shared
    // problem arrays: min/max directions may run concurrently and their cuts
    // are not interchangeable.
    d_row rlt_P_lobo, rlt_P_upbo, rlt_V_lobo, rlt_V_upbo;
    for (int index = 0; index <= VX_MAXDIM; ++index)
    {
        rlt_P_lobo[index] = P_lobo[index];
        rlt_P_upbo[index] = P_upbo[index];
        rlt_V_lobo[index] = V_lobo[index];
        rlt_V_upbo[index] = V_upbo[index];
    }

    gsl_vector *coef_x = managed::vector_calloc(n_variable + 1);
    gsl_vector *initial_solution;
    gsl_vector_view temp_v;
    gsl_vector *whole_solution = managed::vector_calloc(actual_n_variables);
    gsl_vector *onecut = managed::vector_calloc(actual_n_variables + 1);

    const char *configured_method = std::getenv("DECIDEIT_METHOD");
    const std::string method = configured_method ? configured_method : "bcm";
    const char *configured_trigger = std::getenv("DECIDEIT_TRIGGER");
    const int trigger = configured_trigger ? std::stoi(configured_trigger) : 3;
    const char *configured_climb = std::getenv("DECIDEIT_CLIMB");
    const std::string climb = configured_climb ? configured_climb : "guarded";
    const char *configured_rlt_start = std::getenv("DECIDEIT_RLT_START");
    const std::string rlt_start = configured_rlt_start
        ? configured_rlt_start : "staged";
    const char *configured_failure_rescue =
        std::getenv("DECIDEIT_RLT_FAILURE_RESCUE");
    if (configured_failure_rescue &&
        std::strcmp(configured_failure_rescue, "0") != 0 &&
        std::strcmp(configured_failure_rescue, "1") != 0)
        throw std::invalid_argument(
            "DECIDEIT_RLT_FAILURE_RESCUE must be 0 or 1");
    // ODBA is meant to activate RLT when its cut path cannot progress.
    // Preserve pure-CP ablations unless rescue is explicitly requested.
    const bool failure_rescue = configured_failure_rescue
        ? std::strcmp(configured_failure_rescue, "1") == 0
        : method == "odba";
    const bool branch_rlt_tail = branch_rlt_tail_enabled();
    const bool early_branch_rlt = binary_environment_option(
        "DECIDEIT_EARLY_BRANCH_RLT");
    const bool adaptive_branch_rlt_escalation = binary_environment_option(
        "DECIDEIT_BRANCH_RLT_ADAPTIVE_ESCALATION");
    const double branch_rlt_escalation_relative_gap =
        adaptive_branch_rlt_escalation
        ? bounded_environment_double(
            "DECIDEIT_BRANCH_RLT_ESCALATION_REL_GAP", 0.005, 0.0, 1.0)
        : 0.0;
    const int branch_rlt_escalation_max_columns =
        adaptive_branch_rlt_escalation
        ? cone_pilot_limit(
            "DECIDEIT_BRANCH_RLT_ESCALATION_MAX_COLUMNS", 2500, 4000000)
        : 0;
    const bool extend_small_rlt_early_lp = std::getenv(
        "DECIDEIT_EARLY_BRANCH_RLT_SMALL_LP_TIME_LIMIT_MS") != nullptr;
    const int small_rlt_early_lp_time_limit_ms = extend_small_rlt_early_lp
        ? cone_pilot_limit(
            "DECIDEIT_EARLY_BRANCH_RLT_SMALL_LP_TIME_LIMIT_MS", 2500,
            3600000)
        : 0;
    const int small_rlt_early_lp_max_columns = extend_small_rlt_early_lp
        ? cone_pilot_limit(
            "DECIDEIT_EARLY_BRANCH_RLT_SMALL_LP_MAX_COLUMNS", 2500,
            4000000)
        : 0;
    const int early_branch_rlt_window = early_branch_rlt
        ? cone_pilot_limit("DECIDEIT_EARLY_BRANCH_RLT_WINDOW", 4, 100)
        : 0;
    const bool skip_satisfied_rlt = binary_environment_option(
        "DECIDEIT_RLT_SKIP_SATISFIED");
    const PlateauTightenMode plateau_tighten_mode =
        selected_plateau_tighten_mode();
    const int plateau_window = plateau_tighten_mode ==
        PlateauTightenMode::disabled ? 0 :
        cone_pilot_limit("DECIDEIT_RLT_PLATEAU_WINDOW", 3, 100);
    const int plateau_candidate_limit = plateau_tighten_mode ==
        PlateauTightenMode::disabled ? 0 :
        ((plateau_tighten_mode == PlateauTightenMode::cutoff_obbt ||
          plateau_tighten_mode == PlateauTightenMode::reduced_cost_cutoff)
            ? cone_pilot_limit("DECIDEIT_RLT_CUTOFF_OBBT_VARIABLES", 4,
                               2 * VX_MAXDIM)
            : (plateau_tighten_mode == PlateauTightenMode::reduced_cost
                ? cone_pilot_limit(
                    "DECIDEIT_RLT_REDUCED_COST_VARIABLES", 8,
                    2 * VX_MAXDIM)
                : cone_pilot_limit("DECIDEIT_RLT_OBBT_VARIABLES", 8,
                                   2 * VX_MAXDIM)));
    const int plateau_round_limit = plateau_tighten_mode ==
        PlateauTightenMode::disabled ? 0 :
        cone_pilot_limit("DECIDEIT_RLT_TIGHTEN_ROUNDS", 1, 10);
    const bool early_bcm_rlt = method == "bcm" && climb == "plain" &&
        std::getenv("DECIDEIT_BCM_EARLY_RLT");
    bool rlt_active = method == "bcm";
    bool rlt_built = false;
    int pending_failure = OK;
    bool first_climb = true;
    int rlt_noc = dataset_in->ori_num_cons;
    initial_solution = nullptr;
    DecideItDirectionStats &direction_stats = minFlag
        ? decideit_stats.minimum : decideit_stats.maximum;
    gsl_vector *neighbor_seed = climb != "plain"
        ? managed::vector_calloc(get_Ncolumns(lp_Y)) : nullptr;
    bool guarded_restart = false;
    int guarded_restarts = 0;
    bool branch_rlt_tail_attempted = false;
    bool early_branch_rlt_attempted = false;
    bool rlt_solution_certified = false;
    int plateau_count = 0;
    int plateau_tightening_rounds = 0;
    int rebuilt_rlt_cold_solves = 0;
    int rebuilt_rlt_warm_solves = 0;
    double preceding_rlt_bound = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> preceding_rlt_solution;
    auto account_branch_rlt = [&](const BranchRltTailResult &tail,
                                  const char *phase, int failure) {
        direction_stats.lower_bound_calls += tail.lp_solves;
        direction_stats.lower_bound_s += tail.wall_s;
        direction_stats.simplex_iterations += tail.simplex_iterations;
        direction_stats.block_lp_solves += tail.polish_lp_solves;
        if (!std::getenv("DECIDEIT_TRACE")) return;
        char best_open_bound[64];
        if (std::isfinite(tail.best_open_bound))
            std::snprintf(best_open_bound, sizeof(best_open_bound),
                          "%.17g", tail.best_open_bound);
        else
            std::strcpy(best_open_bound, "null");
        fprintf(stderr,
            "TRACE {\"event\":\"branch_rlt_tail\",\"phase\":\"%s\",\"status\":\"%s\",\"sense\":\"%s\",\"failure\":%d,\"initial_incumbent_valid\":%s,\"final_incumbent_valid\":%s,\"nodes\":%d,\"lp_solves\":%d,\"resumed_lp_solves\":%d,\"infeasible\":%d,\"bound_pruned\":%d,\"exact\":%d,\"incumbent_updates\":%d,\"polish_attempts\":%d,\"polish_lp_solves\":%d,\"polish_updates\":%d,\"polish_s\":%.9g,\"warm_basis_uses\":%d,\"reduced_cost_updates\":%d,\"reduced_cost_width_reduction\":%.17g,\"node_model_build_s\":%.9g,\"node_lp_s\":%.9g,\"max_depth\":%d,\"open_nodes\":%d,\"best_open_bound\":%s,\"incumbent\":%.17g,\"max_residual\":%.17g,\"simplex_iterations\":%lld,\"wall_s\":%.9g}\n",
            phase,
            tail.certified ? "certified" :
                (tail.budget_exhausted ? "budget" :
                 (tail.lp_failed ? "lp_failed" : "open")),
            minFlag ? "min" : "negated_max", failure,
            tail.initial_incumbent_valid ? "true" : "false",
            tail.final_incumbent_valid ? "true" : "false",
            tail.nodes, tail.lp_solves, tail.resumed_lp_solves,
            tail.infeasible,
            tail.bound_pruned, tail.exact,
            tail.incumbent_updates, tail.polish_attempts,
            tail.polish_lp_solves, tail.polish_updates,
            tail.polish_s, tail.warm_basis_uses,
            tail.reduced_cost_updates,
            tail.reduced_cost_width_reduction,
            tail.node_model_build_s, tail.node_lp_s,
            tail.maximum_depth, tail.open_nodes,
            best_open_bound, solution->alpha,
            tail.maximum_residual, tail.simplex_iterations,
            tail.wall_s);
    };
    auto run_controlled_branch_rlt = [&](bool early, int failure) {
        bool *attempted = early ? &early_branch_rlt_attempted
                                : &branch_rlt_tail_attempted;
        if (*attempted || (early ? !early_branch_rlt : !branch_rlt_tail))
            return false;
        *attempted = true;
        try
        {
            const double initial_relative_gap = std::isfinite(temp_obj) &&
                std::isfinite(solution->alpha)
                ? (std::max)(0.0, solution->alpha - temp_obj) /
                    (1.0 + std::fabs(solution->alpha))
                : std::numeric_limits<double>::infinity();
            const bool pre_escalate = early && branch_rlt_tail &&
                adaptive_branch_rlt_escalation &&
                !branch_rlt_tail_attempted && rlt_built &&
                LinPara->num_col <= branch_rlt_escalation_max_columns &&
                initial_relative_gap <=
                    branch_rlt_escalation_relative_gap;
            if (pre_escalate)
            {
                branch_rlt_tail_attempted = true;
                if (std::getenv("DECIDEIT_TRACE"))
                    fprintf(stderr,
                        "TRACE {\"event\":\"branch_rlt_escalation\",\"sense\":\"%s\",\"relative_gap\":%.17g,\"threshold\":%.17g,\"rlt_columns\":%d}\n",
                        minFlag ? "min" : "negated_max",
                        initial_relative_gap,
                        branch_rlt_escalation_relative_gap,
                        LinPara->num_col);
            }
            BranchRltLimits limits = selected_branch_rlt_limits(
                pre_escalate ? false : early);
            if (early && !pre_escalate && rlt_built &&
                LinPara->num_col <= small_rlt_early_lp_max_columns &&
                small_rlt_early_lp_time_limit_ms > limits.lp_ms)
            {
                limits.lp_ms = (std::min)(
                    limits.total_ms, small_rlt_early_lp_time_limit_ms);
                if (std::getenv("DECIDEIT_TRACE"))
                    fprintf(stderr,
                        "TRACE {\"event\":\"branch_rlt_small_lp_budget\",\"sense\":\"%s\",\"rlt_columns\":%d,\"lp_time_limit_ms\":%d,\"total_time_limit_ms\":%d}\n",
                        minFlag ? "min" : "negated_max",
                        LinPara->num_col, limits.lp_ms, limits.total_ms);
            }
            const BranchRltTailResult tail = run_branch_rlt_tail(
                rlt_built ? LinPara : nullptr,
                rlt_built ? initial_solution : nullptr, temp_obj,
                instance, dataset_in, fixed_variables,
                LC_P, LC_V, original_P_rows,
                rlt_P_lobo, rlt_P_upbo, rlt_V_lobo, rlt_V_upbo,
                P_lobo, P_upbo, V_lobo, V_upbo,
                Ai_begin, Ai_end, Aj_begin, Aj_end,
                n_alt_i, n_alt_j, mX, mY, *noc,
                minFlag, solution, limits,
                lp_X_std, lp_Y);
            account_branch_rlt(
                tail, pre_escalate ? "escalated" :
                    (early ? "early" : "final"), failure);
            return tail.certified;
        }
        catch (const std::exception &error)
        {
            if (std::getenv("DECIDEIT_TRACE"))
                fprintf(stderr,
                    "TRACE {\"event\":\"branch_rlt_tail\",\"phase\":\"%s\",\"status\":\"discarded\",\"sense\":\"%s\",\"failure\":%d,\"reason\":\"controlled_exception\"}\n",
                    early ? "early" : "final",
                    minFlag ? "min" : "negated_max", failure);
            (void)error;
            return false;
        }
    };
    auto try_branch_rlt_tail = [&](int failure) {
        if (failure != DEGENERACY && failure != CUT_SEPARATION_FAILED)
            return false;
        return run_controlled_branch_rlt(false, failure);
    };
    auto activate_failure_rescue = [&](int failure) {
        if (branch_rlt_tail_attempted ||
            (!failure_rescue && !branch_rlt_tail) || rlt_active ||
            (failure != DEGENERACY && failure != CUT_SEPARATION_FAILED))
            return false;
        // No cut was appended by either failed operation.  The next outer
        // iteration builds/synchronizes the current RLT and performs its
        // ordinary bound-versus-incumbent certificate check without repeating
        // the failed neighborhood search. This is a relaxation of the
        // unexplored region: earlier polar cuts may exclude the incumbent,
        // so its lower bound can legitimately be higher than that incumbent.
        rlt_active = true;
        pending_failure = failure;
        if (std::getenv("DECIDEIT_TRACE"))
            printf("TRACE {\"event\":\"rlt_failure_rescue\",\"method\":\"%s\",\"sense\":\"%s\",\"failure\":%d,\"cuts\":%d}\n",
                   method.c_str(), minFlag ? "min" : "negated_max",
                   failure, n_cuts);
        return true;
    };

    while (loop_flag)
    {
        localFlag = false;
        globalFlag = false;
        exitFlag = OK;
        solverFlag = -100;

        bool neighborhood_ready = false;
        const bool shortened_climb = early_bcm_rlt && first_climb;
        bool climb_converged = true;
        const gsl_vector *climb_seed = guarded_restart ? neighbor_seed
            : (first_climb ? initial_y : nullptr);
        guarded_restart = false;
        first_climb = false;
        int augmented_restarts = 0;
        while (true)
        {
            climb_converged = amc(
                local_amc, lp_X_std, lp_Y, dataset_in,
                climb_seed, &direction_stats, shortened_climb ? 1 : 0);
            if (climb != "augmented") break;

            exitFlag = neighborhood(*noc, dataset_in, lp_X_std, local_amc);
            direction_stats.neighborhood_calls++;
            if (exitFlag < 0) break;
            neighborhood_ready = true;
            if (!improving_neighbor(dataset_in, local_amc, lp_Y,
                                    neighbor_seed, &direction_stats))
                break;
            if (++augmented_restarts > 10000)
                throw std::runtime_error("Augmented climbing restart limit exceeded");
            direction_stats.climb_restarts++;
            neighborhood_ready = false;
            climb_seed = neighbor_seed;
        }

        // Only augmented climbing performs a neighborhood call in the loop
        // above.  Plain/guarded climbing must not inspect a status left by an
        // earlier outer iteration.
        if (climb == "augmented" && exitFlag < 0)
        {
            if (try_branch_rlt_tail(exitFlag))
            {
                if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                localFlag = true;
                solverFlag = OPTIMAL;
                break;
            }
            if (activate_failure_rescue(exitFlag)) continue;
            if (minFlag) P_lobo[0] = exitFlag;
            else P_upbo[0] = exitFlag;
            localFlag = true;
            solverFlag = exitFlag;
            if (exitFlag == GSL_LU && n_cuts == 0) V_point[0] = 0;
            break;
        }

        if (local_amc->cbov < solution->alpha)
        {
            solution->alpha = local_amc->cbov;
            gsl_vector_memcpy(solution->x_solution, local_amc->x_solution);
            gsl_vector_memcpy(solution->y_solution, local_amc->y_solution);
        }

        if (method == "odba" && n_cuts >= trigger)
            rlt_active = true;
        if (rlt_active)
        {
            if (!rlt_built)
            {
                if (LinPara->col_cost)
                    destroy_LinearizationParameters(LinPara);
                construct_LinearizationParameters(
                    mX, mY, actual_n_variables, actual_n_variables, LinPara);
                construct_linearization_objective(n_alt_i, n_alt_j, LinPara);
                if (!minFlag) gsl_vector_scale(LinPara->col_cost, -1.0);
                build_linearization(n_alt_i, n_alt_j, LC_P, LC_V,
                    rlt_P_lobo, rlt_P_upbo, rlt_V_lobo, rlt_V_upbo,
                    Ai_begin, Ai_end, Aj_begin, Aj_end, mX, mY, LinPara);
                initial_solution = managed::vector_calloc(LinPara->n_total);
                rlt_built = true;

                // Optional staged activation: first optimize the smaller base
                // RLT, then append all cuts accumulated before activation and
                // reoptimize from that exact LP basis.  This is a genuine
                // basis warm start, unlike merely supplying a rank-one point.
                if (rlt_start == "staged" && rlt_noc < *noc &&
                    !std::getenv("DECIDEIT_COLD_RLT"))
                {
                    getfrom_x_solution(fixed_variables, whole_solution,
                                       instance->xl, local_amc->x_solution);
                    construct_initial_solution(LinPara, whole_solution,
                                               local_amc->y_solution,
                                               actual_n_variables,
                                               initial_solution);
                    const auto auxiliary_start = std::chrono::steady_clock::now();
                    temp_obj = solve_rlt_lp(
                        LinPara, initial_solution, true,
                        std::getenv("DECIDEIT_LP_LOG") != nullptr,
                        cross_level_rlt_state, &direction_stats);
                    rlt_solution_certified = true;
                    direction_stats.rlt_auxiliary_solves++;
                    direction_stats.rlt_auxiliary_s +=
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - auxiliary_start).count();
                    if (LinPara->persistent)
                        direction_stats.simplex_iterations += static_cast<long long>(
                            LinPara->persistent->highs.getInfo().simplex_iteration_count);
                }
            }
            int newly_appended_cuts = 0;
            double maximum_new_cut_violation = 0.0;
            double maximum_new_rlt_row_violation = 0.0;
            double maximum_new_rlt_relative_violation = 0.0;
            while (rlt_noc < *noc)
            {
                temp_v = gsl_matrix_row(dataset_in->valid_A, rlt_noc++);
                recover_cut(onecut, &temp_v.vector, instance->xl, fixed_variables);
                const bool have_relaxation = rlt_solution_certified;
                if (have_relaxation)
                {
                    double activity = 0.0;
                    for (int column = 0; column < actual_n_variables; ++column)
                        activity += gsl_vector_get(onecut, column + 1) *
                            gsl_vector_get(initial_solution, column);
                    maximum_new_cut_violation = (std::max)(
                        maximum_new_cut_violation,
                        activity - gsl_vector_get(onecut, 0));
                }
                const int first_new_rlt_row = LinPara->num_row;
                linearization_onecut2go(onecut, LinPara,
                                        rlt_P_lobo, rlt_P_upbo,
                                        Ai_begin, Aj_begin, n_alt_i,
                                        &LC_P,
                                        have_relaxation
                                            ? initial_solution->data : nullptr,
                                        have_relaxation
                                            ? initial_solution->size : 0);
                if (have_relaxation)
                    for (int row = first_new_rlt_row;
                         row < LinPara->num_row; ++row)
                    {
                        double activity = 0.0;
                        const int begin = gsl_vector_int_get(
                            LinPara->a_start, row);
                        const int end = gsl_vector_int_get(
                            LinPara->a_start, row + 1);
                        for (int entry = begin; entry < end; ++entry)
                            activity += gsl_vector_get(
                                LinPara->a_value, entry) *
                                gsl_vector_get(
                                    initial_solution,
                                    gsl_vector_int_get(
                                        LinPara->a_index, entry));
                        maximum_new_rlt_row_violation = (std::max)(
                            maximum_new_rlt_row_violation,
                            (std::max)(
                                gsl_vector_get(LinPara->row_lower, row) -
                                    activity,
                                activity - gsl_vector_get(
                                    LinPara->row_upper, row)));
                        const double lower = gsl_vector_get(
                            LinPara->row_lower, row);
                        const double upper = gsl_vector_get(
                            LinPara->row_upper, row);
                        double scale = 1.0 + std::fabs(activity);
                        if (std::isfinite(lower) && lower > -1.0e29)
                            scale = (std::max)(scale,
                                               1.0 + std::fabs(lower));
                        if (std::isfinite(upper) && upper < 1.0e29)
                            scale = (std::max)(scale,
                                               1.0 + std::fabs(upper));
                        maximum_new_rlt_relative_violation = (std::max)(
                            maximum_new_rlt_relative_violation,
                            (std::max)(0.0, (std::max)(
                                lower - activity, activity - upper)) / scale);
                    }
                newly_appended_cuts++;
            }

            const bool skipped_rlt_solve = skip_satisfied_rlt &&
                selected_rlt_backend() == "highs" &&
                rlt_solution_certified &&
                maximum_new_rlt_relative_violation <= 1e-10;
            if (skipped_rlt_solve)
            {
                const int synchronized_rows = LinPara->persistent
                    ? LinPara->num_row - LinPara->persistent->loaded_rows : 0;
                append_rlt_rows_without_solving(LinPara);
                direction_stats.rlt_skipped_solves++;
                if (std::getenv("DECIDEIT_TRACE"))
                    printf("TRACE {\"event\":\"rlt_solve_skip\",\"sense\":\"%s\",\"new_cuts\":%d,\"new_rows\":%d,\"relative_violation\":%.17g,\"lower_bound\":%.17g}\n",
                           minFlag ? "min" : "negated_max",
                           newly_appended_cuts, synchronized_rows,
                           maximum_new_rlt_relative_violation, temp_obj);
            }
            else
            {
                getfrom_x_solution(fixed_variables, whole_solution,
                                   instance->xl, local_amc->x_solution);
                construct_initial_solution(
                    LinPara, whole_solution, local_amc->y_solution,
                    actual_n_variables, initial_solution);
                const auto lower_bound_start =
                    std::chrono::steady_clock::now();
                temp_obj = solve_rlt_lp(
                    LinPara, initial_solution, rlt_start != "none",
                    std::getenv("DECIDEIT_LP_LOG") != nullptr,
                    cross_level_rlt_state, &direction_stats);
                rlt_solution_certified = true;
                direction_stats.lower_bound_calls++;
                direction_stats.lower_bound_s +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() -
                        lower_bound_start).count();
            }
            if (LinPara->persistent)
            {
                direction_stats.cold_solves = rebuilt_rlt_cold_solves +
                    LinPara->persistent->cold_solves;
                direction_stats.warm_solves = rebuilt_rlt_warm_solves +
                    LinPara->persistent->warm_solves;
                if (!skipped_rlt_solve)
                {
                    const HighsInt iterations = LinPara->persistent->highs
                        .getInfo().simplex_iteration_count;
                    if (iterations > 0)
                        direction_stats.simplex_iterations +=
                            static_cast<long long>(iterations);
                }
            }
            if (LinPara->num_row > direction_stats.max_rows)
                direction_stats.max_rows = LinPara->num_row;
            direction_stats.columns = LinPara->num_col;

            auto maximum_rank_one_residual = [&]() {
                double maximum = 0.0;
                for (int i = 0; i < actual_n_variables; ++i)
                    for (int j = 0; j < actual_n_variables; ++j)
                    {
                        const int product = 2 * actual_n_variables +
                            i * actual_n_variables + j;
                        maximum = (std::max)(maximum, std::fabs(
                            gsl_vector_get(initial_solution, product) -
                            gsl_vector_get(initial_solution, i) *
                            gsl_vector_get(initial_solution,
                                           actual_n_variables + j)));
                    }
                return maximum;
            };
            double rank_one_residual = maximum_rank_one_residual();
            const bool have_preceding_bound =
                std::isfinite(preceding_rlt_bound);
            const double bound_improvement = have_preceding_bound
                ? temp_obj - preceding_rlt_bound
                : std::numeric_limits<double>::quiet_NaN();
            double solution_movement = 0.0;
            if (preceding_rlt_solution.size() == initial_solution->size)
                for (size_t column = 0; column < initial_solution->size;
                     ++column)
                    solution_movement = (std::max)(solution_movement,
                        std::fabs(gsl_vector_get(initial_solution, column) -
                                  preceding_rlt_solution[column]));
            const double plateau_tolerance = 1e-8 *
                (1.0 + (have_preceding_bound
                    ? std::fabs(preceding_rlt_bound) : 0.0));
            if (!have_preceding_bound)
                plateau_count = 0;
            else if (newly_appended_cuts > 0)
            {
                if (bound_improvement <= plateau_tolerance)
                    plateau_count++;
                else
                    plateau_count = 0;
            }
            const int plateau_streak = plateau_count;

            if (std::getenv("DECIDEIT_TRACE"))
            {
                char improvement_text[64];
                if (have_preceding_bound)
                    std::snprintf(improvement_text,
                                  sizeof(improvement_text), "%.17g",
                                  bound_improvement);
                else
                    std::strcpy(improvement_text, "null");
                printf("TRACE {\"event\":\"rlt_plateau_profile\",\"sense\":\"%s\",\"cuts\":%d,\"new_cuts\":%d,\"solve_skipped\":%s,\"new_cut_max_violation\":%.17g,\"new_rlt_row_max_violation\":%.17g,\"new_rlt_row_relative_violation\":%.17g,\"lower_bound\":%.17g,\"bound_improvement\":%s,\"solution_linf\":%.17g,\"rank_one_residual\":%.17g,\"plateau_count\":%d}\n",
                       minFlag ? "min" : "negated_max", n_cuts,
                       newly_appended_cuts,
                       skipped_rlt_solve ? "true" : "false",
                       (std::max)(0.0, maximum_new_cut_violation),
                       (std::max)(0.0, maximum_new_rlt_row_violation),
                       maximum_new_rlt_relative_violation,
                       temp_obj, improvement_text, solution_movement,
                       rank_one_residual, plateau_count);
            }

            if (plateau_tighten_mode != PlateauTightenMode::disabled &&
                newly_appended_cuts > 0 &&
                plateau_count >= plateau_window &&
                plateau_tightening_rounds < plateau_round_limit)
            {
                plateau_tightening_rounds++;
                const double bound_before_tightening = temp_obj;
                const PlateauTighteningResult tightening =
                    tighten_plateau_bounds(
                        plateau_tighten_mode, *LinPara, initial_solution,
                        actual_n_variables, n_alt_i, LC_V, LC_P,
                        dataset_in, *noc, fixed_variables,
                        rlt_P_lobo, rlt_P_upbo,
                        rlt_V_lobo, rlt_V_upbo,
                        Ai_begin, Aj_begin, plateau_candidate_limit,
                        solution->alpha);
                direction_stats.simplex_iterations +=
                    tightening.simplex_iterations;
                const bool tightened = tightening.feasible &&
                    !tightening.cutoff_certified &&
                    tightening.width_reduction > 1e-10;
                if (tightened)
                {
                    if (LinPara->persistent)
                    {
                        rebuilt_rlt_cold_solves +=
                            LinPara->persistent->cold_solves;
                        rebuilt_rlt_warm_solves +=
                            LinPara->persistent->warm_solves;
                    }
                    const std::vector<double> separation_primal(
                        initial_solution->data,
                        initial_solution->data + initial_solution->size);
                    destroy_LinearizationParameters(LinPara);
                    rlt_solution_certified = false;
                    construct_LinearizationParameters(
                        mX, mY, actual_n_variables,
                        actual_n_variables, LinPara);
                    construct_linearization_objective(
                        n_alt_i, n_alt_j, LinPara);
                    if (!minFlag) gsl_vector_scale(LinPara->col_cost, -1.0);
                    build_linearization(
                        n_alt_i, n_alt_j, LC_P, LC_V,
                        rlt_P_lobo, rlt_P_upbo,
                        rlt_V_lobo, rlt_V_upbo,
                        Ai_begin, Ai_end, Aj_begin, Aj_end,
                        mX, mY, LinPara);
                    rlt_noc = dataset_in->ori_num_cons;
                    while (rlt_noc < *noc)
                    {
                        temp_v = gsl_matrix_row(
                            dataset_in->valid_A, rlt_noc++);
                        recover_cut(onecut, &temp_v.vector,
                                    instance->xl, fixed_variables);
                        linearization_onecut2go(
                            onecut, LinPara,
                            rlt_P_lobo, rlt_P_upbo,
                            Ai_begin, Aj_begin, n_alt_i, &LC_P,
                            separation_primal.data(),
                            separation_primal.size());
                    }
                    const auto tightened_solve_start =
                        std::chrono::steady_clock::now();
                    temp_obj = solve_rlt_lp(
                        LinPara, initial_solution, false,
                        std::getenv("DECIDEIT_LP_LOG") != nullptr,
                        selected_rlt_backend() == "pivotcopy"
                            ? cross_level_rlt_state : nullptr,
                        &direction_stats);
                    rlt_solution_certified = true;
                    direction_stats.lower_bound_calls++;
                    direction_stats.lower_bound_s +=
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() -
                            tightened_solve_start).count();
                    if (LinPara->persistent)
                    {
                        direction_stats.cold_solves =
                            rebuilt_rlt_cold_solves +
                            LinPara->persistent->cold_solves;
                        direction_stats.warm_solves =
                            rebuilt_rlt_warm_solves +
                            LinPara->persistent->warm_solves;
                        direction_stats.simplex_iterations +=
                            static_cast<long long>(
                                LinPara->persistent->highs.getInfo()
                                    .simplex_iteration_count);
                        if (cross_level_rlt_state)
                        {
                            const HighsBasis &basis =
                                LinPara->persistent->highs.getBasis();
                            cross_level_rlt_state->valid = true;
                            cross_level_rlt_state->basis_valid = basis.valid;
                            cross_level_rlt_state->rows = LinPara->num_row;
                            cross_level_rlt_state->columns = LinPara->num_col;
                            cross_level_rlt_state->basis = basis;
                            cross_level_rlt_state->solution.assign(
                                initial_solution->data,
                                initial_solution->data + LinPara->num_col);
                            cross_level_rlt_state->col_lower.assign(
                                LinPara->col_lower->data,
                                LinPara->col_lower->data + LinPara->num_col);
                            cross_level_rlt_state->col_upper.assign(
                                LinPara->col_upper->data,
                                LinPara->col_upper->data + LinPara->num_col);
                        }
                    }
                    direction_stats.max_rows = (std::max)(
                        direction_stats.max_rows, LinPara->num_row);
                    direction_stats.columns = LinPara->num_col;
                    rank_one_residual = maximum_rank_one_residual();
                    if (temp_obj < bound_before_tightening - 1e-7 *
                            (1.0 + std::fabs(bound_before_tightening)))
                        throw std::runtime_error(
                            "Tightened RLT produced a weaker lower bound");
                }
                if (std::getenv("DECIDEIT_TRACE"))
                    printf("TRACE {\"event\":\"rlt_plateau_tightening\",\"sense\":\"%s\",\"mode\":\"%s\",\"round\":%d,\"feasible\":%s,\"cutoff_certified\":%s,\"cutoff_budget\":%s,\"rebuilt\":%s,\"propagation_updates\":%d,\"selected_x\":%d,\"selected_y\":%d,\"obbt_solves\":%d,\"obbt_updates\":%d,\"reduced_cost_available\":%s,\"reduced_cost_updates\":%d,\"width_reduction\":%.17g,\"bound_before\":%.17g,\"bound_after\":%.17g,\"rank_one_residual\":%.17g,\"simplex_iterations\":%lld,\"wall_s\":%.9g}\n",
                           minFlag ? "min" : "negated_max",
                           plateau_tighten_mode == PlateauTightenMode::obbt
                               ? "obbt" :
                           (plateau_tighten_mode ==
                                PlateauTightenMode::cutoff_obbt
                                ? "cutoff" :
                            (plateau_tighten_mode ==
                                PlateauTightenMode::reduced_cost
                                ? "reduced-cost" :
                            (plateau_tighten_mode ==
                                PlateauTightenMode::reduced_cost_cutoff
                                ? "reduced-cost-cutoff" : "propagate"))),
                           plateau_tightening_rounds,
                           tightening.feasible ? "true" : "false",
                           tightening.cutoff_certified ? "true" : "false",
                           tightening.cutoff_budget_exhausted
                               ? "true" : "false",
                           tightened ? "true" : "false",
                           tightening.propagation_updates,
                           tightening.selected_x, tightening.selected_y,
                           tightening.obbt_solves,
                           tightening.obbt_updates,
                           tightening.reduced_cost_available
                               ? "true" : "false",
                           tightening.reduced_cost_updates,
                           tightening.width_reduction,
                           bound_before_tightening, temp_obj,
                           rank_one_residual,
                           tightening.simplex_iterations,
                           tightening.wall_s);
                if (tightening.cutoff_certified)
                {
                    localFlag = true;
                    if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                    solverFlag = OPTIMAL;
                    break;
                }
                plateau_count = 0;
            }

            preceding_rlt_bound = temp_obj;
            preceding_rlt_solution.assign(
                initial_solution->data,
                initial_solution->data + initial_solution->size);
            double lifted_obj = 0.0;
            for (int z = 0; z < actual_n_variables; ++z)
                lifted_obj += (z < n_alt_i ? 1.0 : -1.0) *
                    gsl_vector_get(initial_solution, z) *
                    gsl_vector_get(initial_solution, actual_n_variables + z);
            if (!minFlag) lifted_obj = -lifted_obj;
            const bool objective_consistent =
                check_optimality(initial_solution, actual_n_variables) &&
                std::fabs(lifted_obj - temp_obj) <= DBLP_EPSILON;

            if (std::getenv("DECIDEIT_TRACE"))
                printf("TRACE {\"method\":\"%s\",\"sense\":\"%s\",\"cuts\":%d,\"rlt_active\":true,\"current_search_value\":%.17g,\"incumbent\":%.17g,\"lower_bound\":%.17g,\"gap\":%.17g,\"rlt_rows\":%d,\"rlt_columns\":%d,\"cold_solves\":%d,\"warm_solves\":%d,\"simplex_iterations\":%d}\n",
                    method.c_str(), minFlag ? "min" : "negated_max", n_cuts,
                    local_amc->cbov, solution->alpha, temp_obj, std::fmax(0.0, solution->alpha-temp_obj),
                    LinPara->num_row, LinPara->num_col,
                    direction_stats.cold_solves, direction_stats.warm_solves,
                    LinPara->persistent
                        ? (std::max)(0, static_cast<int>(
                              LinPara->persistent->highs.getInfo()
                                  .simplex_iteration_count))
                        : static_cast<int>(direction_stats.rlt_pivotcopy_pivots));

            // The RLT bounds the region left after all valid polar cuts.
            // Those cuts may remove the incumbent itself while preserving
            // every point that could improve it.  Thus a residual-region
            // lower bound above the incumbent is the intended certificate.
            if (objective_consistent || temp_obj > solution->alpha - DBLP_EPSILON)
            {
                direction_stats.certificate = objective_consistent ? "rlt_rank_one" : "rlt_residual_gap";
                if (objective_consistent && lifted_obj < solution->alpha)
                {
                    solution->alpha = lifted_obj;
                    temp_v = gsl_vector_subvector(initial_solution, 0, actual_n_variables);
                    putto_global_x_solution(fixed_variables, &temp_v.vector,
                                            solution->x_solution);
                    temp_v = gsl_vector_subvector(initial_solution,
                                                  actual_n_variables,
                                                  actual_n_variables);
                    gsl_vector_memcpy(solution->y_solution, &temp_v.vector);
                }
                localFlag = true;
                if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                solverFlag = OPTIMAL;
                break;
            }

            // Once several consecutive cuts leave the certified RLT optimum
            // feasible and do not move its bound, spend a small, strictly
            // bounded branch-and-RLT budget before generating more geometry.
            // Failure to close the tree is harmless: ODBA resumes and retains
            // any improved feasible incumbent found by the pilot.
            if (early_branch_rlt && !early_branch_rlt_attempted &&
                newly_appended_cuts > 0 &&
                plateau_streak >= early_branch_rlt_window &&
                maximum_new_rlt_relative_violation <= 1e-10 &&
                rank_one_residual > 10.0 * DBLP_EPSILON &&
                run_controlled_branch_rlt(true, OK))
            {
                if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                localFlag = true;
                solverFlag = OPTIMAL;
                break;
            }

            // A failed rescue cannot make the incomplete neighborhood or
            // failed separation valid.  Report the original failure without
            // paying for a second identical geometry search.
            if (pending_failure != OK)
            {
                if (try_branch_rlt_tail(pending_failure))
                {
                    if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                    localFlag = true;
                    solverFlag = OPTIMAL;
                    break;
                }
                if (minFlag) P_lobo[0] = pending_failure;
                else P_upbo[0] = pending_failure;
                localFlag = pending_failure == DEGENERACY;
                globalFlag = pending_failure == CUT_SEPARATION_FAILED;
                solverFlag = pending_failure;
                break;
            }

            /* The one-round BCM seed is enough for a rigorous RLT
               certificate.  If the base RLT is not exact, finish ordinary
               mountain climbing before constructing a polar cut. */
            if (shortened_climb && !climb_converged)
            {
                amc(local_amc, lp_X_std, lp_Y, dataset_in,
                    local_amc->y_solution, &direction_stats);
                if (local_amc->cbov < solution->alpha)
                {
                    solution->alpha = local_amc->cbov;
                    gsl_vector_memcpy(solution->x_solution,
                                      local_amc->x_solution);
                    gsl_vector_memcpy(solution->y_solution,
                                      local_amc->y_solution);
                }
            }
        }
        else if (std::getenv("DECIDEIT_TRACE"))
        {
            printf("TRACE {\"method\":\"%s\",\"sense\":\"%s\",\"cuts\":%d,\"rlt_active\":false,\"current_search_value\":%.17g,\"incumbent\":%.17g}\n",
                   method.c_str(), minFlag ? "min" : "negated_max",
                   n_cuts, local_amc->cbov, solution->alpha);
        }

        // exitFlag = local_phase(local_amc, lp_X_std, lp_Y, dataset_in, *noc, n_cuts); // exitFlag now records the # of vertices adjacent to a local solution

        if (!neighborhood_ready)
        {
            exitFlag = neighborhood(*noc, dataset_in, lp_X_std, local_amc);
            direction_stats.neighborhood_calls++;
        }

        // A negative neighborhood result is a failed enumeration, never a
        // count of adjacent vertices.  In particular NO_T_COL must not pass
        // the positive-count test below and be reported as optimal.
        if (exitFlag < 0)
        {
            if (try_branch_rlt_tail(exitFlag))
            {
                if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                localFlag = true;
                solverFlag = OPTIMAL;
                break;
            }
            if (activate_failure_rescue(exitFlag)) continue;
            if (minFlag) P_lobo[0] = exitFlag; else P_upbo[0] = exitFlag;
            localFlag = true;
            solverFlag = exitFlag;
            if (exitFlag == GSL_LU && n_cuts == 0) V_point[0] = 0;
            break;
        }

        /* exitFlag may be abnormal, but a local solution may exist. update if necessary */
        if (solution->alpha - local_amc->cbov > DBLP_EPSILON)
        {
            solution->alpha = local_amc->cbov;
            gsl_vector_memcpy(solution->x_solution, local_amc->x_solution);
            gsl_vector_memcpy(solution->y_solution, local_amc->y_solution);
        }

        /* A polar cut is safe only after the AMC point has no improving
           adjacent X vertex.  The historical augmented mode paid for this
           scan before every RLT solve.  Guarded mode defers it until the RLT
           relaxation is non-exact and a cut is actually about to be added.
           Thus exact base-RLT cases pay zero neighbor LPs, while cut-producing
           cases retain the local-optimality condition required by the cut. */
        if (climb == "guarded" &&
            improving_neighbor(dataset_in, local_amc, lp_Y,
                               neighbor_seed, &direction_stats))
        {
            if (++guarded_restarts > 10000)
                throw std::runtime_error("Guarded climbing restart limit exceeded");
            direction_stats.climb_restarts++;
            guarded_restart = true;
            continue;
        }

        // printf("\nthe objective value is %f\n", solution->alpha);

        if (exitFlag < n_variable)
        {
            // The surviving polyhedron has lower local dimension.  The AMC
            // point is final on that face; this is the termination rule used
            // by the archived CP experiment (the reverse-search truncation
            // signal is handled separately above).
            localFlag = true;
            solverFlag = OPTIMAL;
            if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
            break;
        }

        /* start the global optimization phase */
        exitFlag = global_phase(local_amc, solution->alpha, lp_X_std, lp_Y, dataset_in, *noc);

        // Failed cut construction did not append a constraint.  Do not count
        // it as a cut or advance the active-row count in the result record.
        if (exitFlag >= 0 || exitFlag == ALL_NEGATIVE)
        {
            n_cuts++;
            (*noc)++;
        }

        // if(!consistency_check(dataset_in, local_amc->x_solution, *noc))
        // {
        //     printf("the local solution is infeasible with the new introduced cut\n");
        // }

        if (exitFlag == DEGENERACY)
        {
            if (try_branch_rlt_tail(exitFlag))
            {
                if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                localFlag = true;
                solverFlag = OPTIMAL;
                break;
            }
            if (activate_failure_rescue(exitFlag)) continue;
            if (minFlag)
            {
                P_lobo[0] = -6;
            }
            else
            {
                P_upbo[0] = -6;
            }

            localFlag = true;
            solverFlag = DEGENERACY;
            break;
        }

        if (exitFlag == CUT_SEPARATION_FAILED)
        {
            // Experimental cone-decomposition certificate.  The first pass
            // only exports the exact ray geometry.  A separately validated
            // triangulation enables a second pass that solves one child RLT
            // per simplicial cone.  Invalid or incomplete pilot input is
            // deliberately non-fatal and falls through to the ordinary ODBA
            // failure/rescue path.
            try
            {
                const ConePilotResult cone_result = run_cone_rlt_pilot(
                    local_amc->cone_geometry,
                    rlt_built ? LinPara : nullptr,
                    instance, fixed_variables, P_lobo, P_upbo,
                    Ai_begin, Aj_begin, n_alt_i, solution->alpha, minFlag);
                if (cone_result.attempted)
                {
                    direction_stats.lower_bound_calls += cone_result.children;
                    direction_stats.lower_bound_s += cone_result.wall_s;
                    direction_stats.simplex_iterations +=
                        cone_result.simplex_iterations;
                }
                if (cone_result.certified)
                {
                    localFlag = true;
                    solverFlag = OPTIMAL;
                    if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                    break;
                }
            }
            catch (const std::exception &error)
            {
                if (std::getenv("DECIDEIT_TRACE"))
                    fprintf(stderr,
                            "TRACE {\"event\":\"cone_rlt_pilot\",\"status\":\"discarded\",\"sense\":\"%s\",\"reason\":\"%s\"}\n",
                            minFlag ? "min" : "negated_max", error.what());
            }
            if (try_branch_rlt_tail(exitFlag))
            {
                if (minFlag) P_lobo[0] = 0; else P_upbo[0] = 0;
                localFlag = true;
                solverFlag = OPTIMAL;
                break;
            }
            if (activate_failure_rescue(exitFlag)) continue;
            if (minFlag) P_lobo[0] = CUT_SEPARATION_FAILED;
            else P_upbo[0] = CUT_SEPARATION_FAILED;
            globalFlag = true;
            solverFlag = CUT_SEPARATION_FAILED;
            break;
        }

        if (exitFlag == GSL_LU)
        {
            globalFlag = true;
            solverFlag = GSL_LU;

            if (minFlag)
            {
                P_lobo[0] = -4;
            }
            else
            {
                P_upbo[0] = -4;
            }

            break;
        }

        if (exitFlag == SOLUTION_LP)
        {
            globalFlag = true;
            solverFlag = SOLUTION_LP;

            if (minFlag)
            {
                P_lobo[0] = -3;
            }
            else
            {
                P_upbo[0] = -3;
            }

            break;
        }

        if (exitFlag == ALL_NEGATIVE) // all negative cutting points
        {
            globalFlag = true;

            if (minFlag)
            {
                P_lobo[0] = 0;
            }
            else
            {
                P_upbo[0] = 0;
            }

            solverFlag = ALL_NEGATIVE;
            break;
        }

        // print_lp(lp_X_std);

        gsl_vector_set_all(coef_x, 1);
        set_obj_fn(lp_X_std, coef_x->data);

        solverFlag = solve(lp_X_std); // the new cut has already been introduced into lp_X_std
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"post_cut_lp\",\"status\":%d,\"cuts\":%d,\"rows\":%d}\n",
                    solverFlag, n_cuts, get_Nrows(lp_X_std));

        switch (solverFlag) // solver_flag already exists
        {
        case INFEASIBLE: // an empty feasible region
            if (minFlag)
            {
                P_lobo[0] = 0;
            }
            else
            {
                P_upbo[0] = 0;
            }

            loop_flag = FALSE;
            solverFlag = INFEASIBLE;
            break;

        case OPTIMAL: // the feasible region is non-empty
            if (n_cuts >= CUTS_MAX)
            {
                if (minFlag)
                {
                    P_lobo[0] = -5;
                }
                else
                {
                    P_upbo[0] = -5;
                }

                loop_flag = FALSE;
                solverFlag = CUT_LIMIT;
            }

            break;
        default:
            loop_flag = FALSE;
            solverFlag = SOLUTION_LP;
            break;
        }
    }

    gsl_vector_set_all(whole_solution, -1);

    record_results(solution, P_point, V_point, Ai_begin, Ai_end, Aj_begin, Aj_end, n_cuts, solverFlag, localFlag, globalFlag, fixed_variables, instance, whole_solution);

    DecideItDirectionStats &final_stats = minFlag
        ? decideit_stats.minimum : decideit_stats.maximum;
    final_stats.cuts = n_cuts;
    final_stats.termination_flag = solverFlag;
    final_stats.residual_bound = temp_obj;
    if (final_stats.certificate == "not_recorded")
        final_stats.certificate = (solverFlag == OPTIMAL || solverFlag == INFEASIBLE || solverFlag == ALL_NEGATIVE)
            ? "controller_geometric_or_tail" : "uncertified";
    if (LinPara->persistent && LinPara->persistent->solved) {
        std::lock_guard<std::recursive_mutex> info_lock(highs_api_mutex());
        const HighsInfo& info = LinPara->persistent->highs.getInfo();
        final_stats.rlt_primal_infeasibility = info.max_primal_infeasibility;
        final_stats.rlt_dual_infeasibility = info.max_dual_infeasibility;
    }

    if (std::getenv("DECIDEIT_TRACE"))
        printf("TRACE {\"event\":\"termination\",\"method\":\"%s\",\"sense\":\"%s\",\"cuts\":%d,\"solver_flag\":%d,\"local\":%s,\"global\":%s}\n",
               method.c_str(), minFlag ? "min" : "negated_max", n_cuts,
               solverFlag, localFlag ? "true" : "false",
               globalFlag ? "true" : "false");

    // consistency_check(instance, whole_solution);

    managed::vector_free(whole_solution);
    managed::vector_free(coef_x);
    managed::vector_free(onecut);
    managed::vector_free(initial_solution);
    managed::vector_free(neighbor_seed);

    if (solverFlag != OPTIMAL && solverFlag != INFEASIBLE && solverFlag != ALL_NEGATIVE)
    {
        exitFlag = SLOW_CONVERGENCE;
    }
    else
    {
        exitFlag = OK;
    }

    return exitFlag;
}
