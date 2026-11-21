#ifndef DECIDEIT_CPLEX_REFERENCE_H
#define DECIDEIT_CPLEX_REFERENCE_H

#include "dblp.h"
#include <string>

struct CplexRunStats
{
    double build_s = 0.0;
    double min_s = 0.0;
    double max_s = 0.0;
    double min_bound = std::numeric_limits<double>::quiet_NaN();
    double max_bound = std::numeric_limits<double>::quiet_NaN();
    long long min_nodes = 0;
    long long max_nodes = 0;
    int min_status = 0;
    int max_status = 0;
    std::string min_status_text;
    std::string max_status_text;
    std::string version;
};

const CplexRunStats &last_cplex_stats();

int cplex_solve_pair(
    d_row P_lobo, d_row P_upbo, link_set P_links,
    d_row V_lobo, d_row V_upbo, link_set V_links,
    d_row P_min_point, d_row V_min_point,
    d_row P_max_point, d_row V_max_point,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    double time_limit_seconds, int threads, bool display_log);

int cplex_solve_pair_rows(
    d_row P_lobo, d_row P_upbo, const std::vector<LinCons> &P_rows,
    d_row V_lobo, d_row V_upbo, const std::vector<LinCons> &V_rows,
    d_row P_min_point, d_row V_min_point,
    d_row P_max_point, d_row V_max_point,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    double time_limit_seconds, int threads, bool display_log);

#endif
