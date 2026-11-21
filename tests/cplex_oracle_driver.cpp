#include "cplex_reference.h"
#include <fstream>
#include <iomanip>

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    std::ifstream in(argv[1]);
    int a, b; in >> a >> b; const int n = a + b;
    d_row pl = {}, pu = {}, vl = {}, vu = {};
    d_row pmin = {}, vmin = {}, pmax = {}, vmax = {};
    link_set pc = {}, vc = {};
    for (int i = 1; i <= n; ++i) in >> pl[i] >> pu[i] >> vl[i] >> vu[i];
    int np, nv; in >> np; pc[0].head = np;
    for (int i = 1; i <= np; ++i) {
        in >> pc[i].head >> pc[i].tail >> pc[i].value; pc[i].cmp = '<';
    }
    in >> nv; vc[0].head = nv;
    for (int i = 1; i <= nv; ++i) {
        in >> vc[i].head >> vc[i].tail >> vc[i].value; vc[i].cmp = '<';
    }
    if (!in) return 3;
    const int status = cplex_solve_pair(
        pl, pu, pc, vl, vu, vc, pmin, vmin, pmax, vmax,
        1, a, a + 1, n, 60.0, 1, false);
    std::cout << std::setprecision(17) << "AUDIT " << status << " "
              << pmin[0] << " " << -pmax[0] << " " << pl[0] << " " << pu[0]
              << " " << vmin[0] << " " << vmax[0];
    for (int i = 1; i <= n; ++i)
        std::cout << " " << pmin[i] << " " << vmin[i]
                  << " " << pmax[i] << " " << vmax[i];
    std::cout << "\n";
    return status == OK ? 0 : 1;
}
