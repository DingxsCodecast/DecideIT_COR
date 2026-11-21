#include "dblp.h"
#include <cmath>
#include <iostream>

int main()
{
    d_row pl = {}, pu = {}, vl = {}, vu = {};
    d_row pmin = {}, vmin = {}, pmax = {}, vmax = {};
    link_set pc = {}, vc = {};
    for (int i = 1; i <= 4; ++i) { pl[i] = 0.0; pu[i] = 1.0; }
    vl[1] = vu[1] = 0.5;
    vl[2] = 0.0; vu[2] = 1.0;
    vl[3] = vu[3] = 0.4;
    vl[4] = 0.0; vu[4] = 1.0;

    /* Row 0 is slack and removed by redundancy().  Row 1 becomes a
       redundant zero row after fixed-variable substitution.  Row 2 remains. */
    vc[0].head = 3;
    vc[1] = {1, 3, '<', 0.2};
    vc[2] = {1, 3, '<', 0.1};
    vc[3] = {2, 4, '<', 0.0};

    std::vector<LinCons> lp(1), lv(4);
    int mx = -1, my = -1;
    setenv("DECIDEIT_ROW_FILTER", "legacy", 1);
    const int status = cpm_solver(
        pl, pu, pc, vl, vu, vc, pmin, vmin, pmax, vmax,
        1, 2, 3, 4, 30, lv, lp, &mx, &my);
    const double legacy_min = pmin[0], legacy_max = -pmax[0];
    const RowReductionStats legacy = last_row_reduction_stats();

    mx = -1;
    my = -1;
    setenv("DECIDEIT_ROW_FILTER", "hybrid", 1);
    const int hybrid_status = cpm_solver(
        pl, pu, pc, vl, vu, vc, pmin, vmin, pmax, vmax,
        1, 2, 3, 4, 30, lv, lp, &mx, &my);
    const double hybrid_min = pmin[0], hybrid_max = -pmax[0];
    const RowReductionStats hybrid = last_row_reduction_stats();

    mx = -1;
    my = -1;
    unsetenv("DECIDEIT_ROW_FILTER");
    const int cheap_status = cpm_solver(
        pl, pu, pc, vl, vu, vc, pmin, vmin, pmax, vmax,
        1, 2, 3, 4, 30, lv, lp, &mx, &my);
    const RowReductionStats cheap = last_row_reduction_stats();

    if (status != OK || hybrid_status != OK || cheap_status != OK ||
        std::fabs(hybrid_min - legacy_min) > 1e-10 ||
        std::fabs(hybrid_max - legacy_max) > 1e-10 ||
        std::fabs(pmin[0] - legacy_min) > 1e-10 ||
        std::fabs((-pmax[0]) - legacy_max) > 1e-10 ||
        std::fabs(pmin[0] + 1.0) > 1e-7 ||
        std::fabs((-pmax[0]) - 0.6) > 1e-7)
    {
        std::cerr << "fixed reduction failed: legacy status=" << status
                  << " hybrid status=" << hybrid_status
                  << " cheap status=" << cheap_status
                  << " min=" << pmin[0] << " max=" << -pmax[0] << "\n";
        return 1;
    }
    if (hybrid.lp_calls >= legacy.lp_calls || hybrid.cheap_removed == 0)
    {
        std::cerr << "hybrid reduction did not reduce LP work: legacy="
                  << legacy.lp_calls << " hybrid=" << hybrid.lp_calls
                  << " cheap_removed=" << hybrid.cheap_removed << "\n";
        return 1;
    }
    if (cheap.lp_calls != 0 || cheap.cheap_removed == 0 ||
        cheap.retained_rows < hybrid.retained_rows)
    {
        std::cerr << "cheap reduction invariant failed: LP calls="
                  << cheap.lp_calls << " cheap_removed=" << cheap.cheap_removed
                  << " retained=" << cheap.retained_rows
                  << " hybrid_retained=" << hybrid.retained_rows << "\n";
        return 1;
    }
    std::cout << "Fixed-variable reduction passed: legacy LP calls="
              << legacy.lp_calls << ", hybrid LP calls=" << hybrid.lp_calls
              << ", cheap LP calls=" << cheap.lp_calls
              << "\n";
    return 0;
}
