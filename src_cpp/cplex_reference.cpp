#include "cplex_reference.h"
#include <ilcplex/cplex.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <numeric>
#include <sstream>
#include <vector>

namespace
{
thread_local CplexRunStats stats;

std::string cplex_error(CPXENVptr env, int code)
{
    char buffer[CPXMESSAGEBUFSIZE] = {};
    if (CPXgeterrorstring(env, code, buffer))
        return "CPLEX error " + std::to_string(code);
    return std::string(buffer);
}

void check(CPXENVptr env, int code, const char *operation)
{
    if (code)
        throw std::runtime_error(std::string(operation) + ": " + cplex_error(env, code));
}

struct Environment
{
    CPXENVptr value = nullptr;
    Environment()
    {
        int status = 0;
        value = CPXopenCPLEX(&status);
        if (!value)
            throw std::runtime_error("Cannot open CPLEX environment: " + cplex_error(value, status));
    }
    ~Environment() { if (value) CPXcloseCPLEX(&value); }
};

struct Problem
{
    CPXENVptr env;
    CPXLPptr value = nullptr;
    explicit Problem(CPXENVptr environment) : env(environment)
    {
        int status = 0;
        value = CPXcreateprob(env, &status, "DecideIT_DBLP");
        check(env, status, "CPXcreateprob");
        if (!value) throw std::runtime_error("CPXcreateprob returned a null model");
    }
    ~Problem() { if (value) CPXfreeprob(env, &value); }
};

bool included(int i, int ai, int ae, int bi, int be)
{
    return (i >= ai && i <= ae) || (i >= bi && i <= be);
}

int local_index(int original, int ai, int ae, int bi)
{
    return original <= ae ? original - ai : (ae - ai + 1) + original - bi;
}

void append_linear_rows(
    const std::vector<LinCons> &rows, int column_offset,
    int ai, int ae, int bi, int be,
    std::vector<int> &beg, std::vector<int> &ind,
    std::vector<double> &val, std::vector<double> &rhs,
    std::vector<char> &sense)
{
    for (const LinCons &row : rows)
    {
        bool use = true;
        for (int index : row.index)
            use = use && included(index, ai, ae, bi, be);
        if (!use) continue;
        beg.push_back(static_cast<int>(ind.size()));
        for (size_t k = 0; k < row.index.size(); ++k)
        {
            ind.push_back(column_offset + local_index(row.index[k], ai, ae, bi));
            val.push_back(row.coefficient[k]);
        }
        rhs.push_back(row.b);
        sense.push_back('L');
    }
}

bool optimal_status(int status)
{
    return status == CPX_STAT_OPTIMAL ||
           status == CPXMIP_OPTIMAL ||
           status == CPXMIP_OPTIMAL_TOL;
}

std::string status_text(CPXENVptr env, int status)
{
    char buffer[CPXMESSAGEBUFSIZE] = {};
    CPXgetstatstring(env, status, buffer);
    return std::string(buffer);
}

struct OneSolve
{
    int status = 0;
    double seconds = 0.0;
    double objective = std::numeric_limits<double>::quiet_NaN();
    double bound = std::numeric_limits<double>::quiet_NaN();
    long long nodes = 0;
    std::string text;
    std::vector<double> x;
};

OneSolve optimize(CPXENVptr env, CPXLPptr lp, int sense, int columns)
{
    check(env, CPXchgobjsen(env, lp, sense), "CPXchgobjsen");
    const auto begin = std::chrono::steady_clock::now();
    check(env, CPXqpopt(env, lp), "CPXqpopt");
    OneSolve result;
    result.seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    result.status = CPXgetstat(env, lp);
    result.text = status_text(env, result.status);
    result.nodes = static_cast<long long>(CPXgetnodecnt(env, lp));
    (void)CPXgetbestobjval(env, lp, &result.bound);

    if (optimal_status(result.status))
    {
        check(env, CPXgetobjval(env, lp, &result.objective), "CPXgetobjval");
        result.x.resize(columns);
        check(env, CPXgetx(env, lp, result.x.data(), 0, columns - 1), "CPXgetx");
    }
    return result;
}

void verify_solution(
    const std::vector<double> &x, double objective,
    d_row pl, d_row pu, const std::vector<LinCons> &pc,
    d_row vl, d_row vu, const std::vector<LinCons> &vc,
    int ai, int ae, int bi, int be)
{
    const int n = ae - ai + 1 + be - bi + 1;
    double recomputed = 0.0, sum_a = 0.0, sum_b = 0.0;
    for (int original = ai; original <= be; ++original)
    {
        if (!included(original, ai, ae, bi, be)) continue;
        const int j = local_index(original, ai, ae, bi);
        const double p = x[j], v = x[n + j];
        if (!std::isfinite(p) || !std::isfinite(v) ||
            p < pl[original] - 1e-7 || p > pu[original] + 1e-7 ||
            v < vl[original] - 1e-7 || v > vu[original] + 1e-7)
            throw std::runtime_error("CPLEX returned a point outside an original bound");
        recomputed += (original <= ae ? 1.0 : -1.0) * p * v;
        (original <= ae ? sum_a : sum_b) += p;
    }
    if (std::fabs(sum_a - 1.0) > 1e-7 || std::fabs(sum_b - 1.0) > 1e-7)
        throw std::runtime_error("CPLEX probabilities do not sum to one");
    for (int base = 0; base < 2; ++base)
    {
        const std::vector<LinCons> &rows = base ? vc : pc;
        const int offset = base ? n : 0;
        for (const LinCons &row : rows)
        {
            bool use = true;
            for (int index : row.index)
                use = use && included(index, ai, ae, bi, be);
            if (!use) continue;
            double lhs = 0.0;
            for (size_t k = 0; k < row.index.size(); ++k)
            {
                lhs += row.coefficient[k] *
                    x[offset + local_index(row.index[k], ai, ae, bi)];
            }
            if (lhs > row.b + 1e-7)
                throw std::runtime_error("CPLEX returned a point violating a linear statement");
        }
    }
    if (!std::isfinite(objective) || std::fabs(recomputed - objective) > 1e-7)
        throw std::runtime_error("CPLEX objective disagrees with its returned point");
}
} // namespace

const CplexRunStats &last_cplex_stats() { return stats; }

int cplex_solve_pair_rows(
    d_row pl, d_row pu, const std::vector<LinCons> &pc,
    d_row vl, d_row vu, const std::vector<LinCons> &vc,
    d_row pmin, d_row vmin, d_row pmax, d_row vmax,
    int ai, int ae, int bi, int be,
    double time_limit_seconds, int threads, bool display_log)
{
    stats = CplexRunStats{};
    if (threads < 1 || time_limit_seconds < 0.0)
        throw std::invalid_argument("Invalid CPLEX thread count or time limit");
    const auto build_begin = std::chrono::steady_clock::now();
    Environment holder;
    CPXENVptr env = holder.value;
    stats.version = CPXversion(env);

    check(env, CPXsetintparam(env, CPXPARAM_ScreenOutput, display_log ? CPX_ON : CPX_OFF), "CPLEX ScreenOutput");
    check(env, CPXsetintparam(env, CPXPARAM_OptimalityTarget, CPX_OPTIMALITYTARGET_OPTIMALGLOBAL), "CPLEX OptimalityTarget");
    check(env, CPXsetintparam(env, CPXPARAM_Threads, threads), "CPLEX Threads");
    check(env, CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC), "CPLEX Parallel");
    check(env, CPXsetintparam(env, CPXPARAM_RandomSeed, 1), "CPLEX RandomSeed");
    if (std::getenv("DECIDEIT_CPLEX_NUMERICAL"))
    {
        check(env, CPXsetintparam(env, CPXPARAM_Emphasis_Numerical, CPX_ON),
              "CPLEX numerical emphasis");
        check(env, CPXsetdblparam(env,
                                  CPXPARAM_Simplex_Tolerances_Optimality,
                                  1e-9),
              "CPLEX optimality tolerance");
    }
    check(env, CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9), "CPLEX MIPGap");
    check(env, CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_AbsMIPGap, 1e-9), "CPLEX AbsMIPGap");
    check(env, CPXsetdblparam(env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9), "CPLEX feasibility tolerance");
    if (std::getenv("DECIDEIT_CPLEX_NO_PRESOLVE"))
    {
        check(env, CPXsetintparam(env, CPXPARAM_Preprocessing_Presolve,
                                  CPX_OFF),
              "CPLEX presolve");
        check(env, CPXsetintparam(env, CPXPARAM_Preprocessing_Aggregator, 0),
              "CPLEX aggregator");
        check(env, CPXsetintparam(env, CPXPARAM_Preprocessing_Reduce, 0),
              "CPLEX reductions");
    }
    if (time_limit_seconds > 0.0)
        check(env, CPXsetdblparam(env, CPXPARAM_TimeLimit, time_limit_seconds), "CPLEX TimeLimit");

    Problem model(env);
    const int n_a = ae - ai + 1, n_b = be - bi + 1, n = n_a + n_b;
    const int columns = 2 * n;
    std::vector<double> obj(columns, 0.0), lower(columns), upper(columns);
    for (int original = ai; original <= be; ++original)
    {
        if (!included(original, ai, ae, bi, be)) continue;
        const int j = local_index(original, ai, ae, bi);
        lower[j] = pl[original]; upper[j] = pu[original];
        lower[n + j] = vl[original]; upper[n + j] = vu[original];
    }
    check(env, CPXnewcols(env, model.value, columns, obj.data(), lower.data(), upper.data(), nullptr, nullptr), "CPXnewcols");

    std::vector<int> beg, ind;
    std::vector<double> val, rhs;
    std::vector<char> row_sense;
    beg.push_back(0);
    for (int j = 0; j < n_a; ++j) { ind.push_back(j); val.push_back(1.0); }
    rhs.push_back(1.0); row_sense.push_back('E');
    beg.push_back(static_cast<int>(ind.size()));
    for (int j = n_a; j < n; ++j) { ind.push_back(j); val.push_back(1.0); }
    rhs.push_back(1.0); row_sense.push_back('E');
    append_linear_rows(pc, 0, ai, ae, bi, be, beg, ind, val, rhs, row_sense);
    append_linear_rows(vc, n, ai, ae, bi, be, beg, ind, val, rhs, row_sense);
    check(env, CPXaddrows(env, model.value, 0, static_cast<int>(rhs.size()),
                          static_cast<int>(val.size()), rhs.data(), row_sense.data(),
                          beg.data(), ind.data(), val.data(), nullptr, nullptr), "CPXaddrows");

    for (int j = 0; j < n; ++j)
    {
        const double sign = j < n_a ? 1.0 : -1.0;
        check(env, CPXchgqpcoef(env, model.value, j, n + j, sign), "CPXchgqpcoef");
    }
    stats.build_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - build_begin).count();

    OneSolve lo = optimize(env, model.value, CPX_MIN, columns);
    stats.min_s = lo.seconds; stats.min_bound = lo.bound;
    stats.min_nodes = lo.nodes; stats.min_status = lo.status;
    stats.min_status_text = lo.text;
    if (!optimal_status(lo.status))
        return (lo.status == CPX_STAT_INFEASIBLE || lo.status == CPXMIP_INFEASIBLE)
                   ? INCONSISTENT : SOLVER_ABORTED;
    verify_solution(lo.x, lo.objective, pl, pu, pc, vl, vu, vc, ai, ae, bi, be);

    // Optional diagnostic start: whitespace-separated P values followed by V
    // values in the local alternative-pair order.  CPLEX checks feasibility
    // before accepting it, which is useful when auditing another solver's
    // incumbent against the independent global model.
    if (const char *start_path = std::getenv("DECIDEIT_CPLEX_START_POINT"))
    {
        std::ifstream source(start_path);
        std::vector<double> start(columns);
        for (double &value : start)
            if (!(source >> value))
                throw std::runtime_error("Invalid CPLEX diagnostic start point");
        double surplus = 0.0;
        if (source >> surplus)
            throw std::runtime_error("Too many values in CPLEX diagnostic start point");
        std::vector<int> indices(columns);
        std::iota(indices.begin(), indices.end(), 0);
        const int begin = 0;
        const int effort = CPX_MIPSTART_CHECKFEAS;
        check(env, CPXaddmipstarts(env, model.value, 1, columns, &begin,
                                  indices.data(), start.data(), &effort,
                                  nullptr),
              "CPXaddmipstarts diagnostic incumbent");
    }

    OneSolve hi = optimize(env, model.value, CPX_MAX, columns);
    stats.max_s = hi.seconds; stats.max_bound = hi.bound;
    stats.max_nodes = hi.nodes; stats.max_status = hi.status;
    stats.max_status_text = hi.text;
    if (!optimal_status(hi.status))
        return (hi.status == CPX_STAT_INFEASIBLE || hi.status == CPXMIP_INFEASIBLE)
                   ? INCONSISTENT : SOLVER_ABORTED;
    verify_solution(hi.x, hi.objective, pl, pu, pc, vl, vu, vc, ai, ae, bi, be);

    for (int original = ai; original <= be; ++original)
    {
        if (!included(original, ai, ae, bi, be)) continue;
        const int j = local_index(original, ai, ae, bi);
        pmin[original] = lo.x[j]; vmin[original] = lo.x[n + j];
        pmax[original] = hi.x[j]; vmax[original] = hi.x[n + j];
    }
    pmin[0] = lo.objective;
    pmax[0] = -hi.objective;
    vmin[0] = 0.0;
    vmax[0] = 0.0;
    return OK;
}

int cplex_solve_pair(
    d_row pl, d_row pu, link_set pc, d_row vl, d_row vu, link_set vc,
    d_row pmin, d_row vmin, d_row pmax, d_row vmax,
    int ai, int ae, int bi, int be,
    double time_limit_seconds, int threads, bool display_log)
{
    std::vector<LinCons> probability_rows, value_rows;
    for (int k = 1; k <= pc[0].head; ++k)
        probability_rows.push_back(pc[k].cmp == '<'
            ? difference_row(pc[k].head, pc[k].tail, pc[k].value)
            : difference_row(pc[k].tail, pc[k].head, -pc[k].value));
    for (int k = 1; k <= vc[0].head; ++k)
        value_rows.push_back(vc[k].cmp == '<'
            ? difference_row(vc[k].head, vc[k].tail, vc[k].value)
            : difference_row(vc[k].tail, vc[k].head, -vc[k].value));
    return cplex_solve_pair_rows(pl, pu, probability_rows, vl, vu, value_rows,
        pmin, vmin, pmax, vmax, ai, ae, bi, be,
        time_limit_seconds, threads, display_log);
}
