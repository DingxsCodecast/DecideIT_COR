#ifndef DECIDEIT_LP_BACKEND_H
#define DECIDEIT_LP_BACKEND_H
#include "Highs.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

using REAL = double;
constexpr bool TRUE=true, FALSE=false;
constexpr int OPTIMAL=0, INFEASIBLE=2, UNBOUNDED=3, NUMFAILURE=5;
// Keep the compatibility API values out of the global names used by solver
// headers.  HiGHS 1.12+ exposes ConstraintType::EQ in the global namespace.
constexpr int kLpLe = 1;
constexpr int kLpGe = 2;
constexpr int kLpEq = 3;
constexpr int kLpNeutral = 0;

// Project-owned LP interface. Row/column arguments are one-based, and objective
// arrays have an unused element zero, matching the existing geometry routines.
// The public build uses HiGHS for all LPs.
struct LpModel {
    std::unique_ptr<Highs> highs;
    int columns, last_status=NUMFAILURE;
    bool loaded=false;
    ObjSense sense=ObjSense::kMinimize;
    std::vector<std::vector<double>> rows;
    std::vector<double> cost, lower, upper, rhs;
    std::vector<int> types;
    std::vector<double> solution, row_activity;
    std::vector<int> row_basis;
    double objective=0.0;
    bool pivotcopy_solution=false;
    bool needs_solve=true;
    std::uint64_t structure_revision=1;
    std::string backend;
    LpModel(int m,int n);
    ~LpModel();
    void require_solution() const;
    Highs& ensure_highs();
};

struct LpBackendStats {
    long long highs_solves=0;
    long long pivotcopy_solves=0;
    long long pivotcopy_pivots=0;
    long long pivotcopy_warm_solves=0;
    long long pivotcopy_warm_retries=0;
    long long pivotcopy_failures=0;
    long long pivotcopy_fallbacks=0;
    long long cached_solves=0;
    int max_columns=0;
    int max_rows=0;
    int max_highs_columns=0;
    int max_highs_rows=0;
    int max_pivotcopy_columns=0;
    int max_pivotcopy_rows=0;
    double highs_s=0.0;
    double pivotcopy_s=0.0;
};

void reset_lp_backend_stats();
LpBackendStats last_lp_backend_stats();
std::string selected_lp_backend();
void clear_lp_solver(LpModel* lp);
long long get_simplex_iterations(const LpModel* lp);
void highs_require(HighsStatus status, const char* operation);
// Some packaged HiGHS builds are not safe when independent Highs instances
// are entered concurrently. All project HiGHS API calls share this mutex;
// Pivot-Copy and the surrounding geometry remain fully parallel.
std::recursive_mutex& highs_api_mutex();
void configure_highs(Highs& highs, bool output=false);
bool highs_solution_valid(const Highs& highs);
LpModel* make_lp(int rows,int columns);
void delete_lp(LpModel* lp);
void free_lp(LpModel** lp);
int get_Ncolumns(const LpModel* lp);
int get_Nrows(const LpModel* lp);
double get_infinite(const LpModel* lp);
void set_verbose(LpModel* lp,int level);
void set_minim(LpModel* lp);
void set_maxim(LpModel* lp);
void set_lowbo(LpModel* lp,int col,double value);
void set_upbo(LpModel* lp,int col,double value);
void set_row(LpModel* lp,int row,const double* values);
void set_rowex(LpModel* lp,int row,int count,const double* values,const int* indices);
void set_constr_type(LpModel* lp,int row,int type);
void set_rh(LpModel* lp,int row,double value);
void set_obj_fn(LpModel* lp,const double* values);
void set_obj_fnex(LpModel* lp,int count,const double* values,const int* indices);
void add_constraint(LpModel* lp,const double* values,int type,double rhs);
int solve(LpModel* lp);
double get_objective(const LpModel* lp);
void get_variables(const LpModel* lp,double* values);
void get_constraints(const LpModel* lp,double* values);
void get_basis(LpModel* lp,int* basis,bool include_nonbasic);
void get_column(const LpModel* lp,int column,double* values);
double get_rh(const LpModel* lp,int row);
void print_lp(const LpModel* lp);
#endif
