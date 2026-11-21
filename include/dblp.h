#ifndef DECIDEIT_DBLP_H
#define DECIDEIT_DBLP_H

// #define VX_MAXDIM  100    // max nbr of prob vars
#define VX_MAXVER 100000 // max nbr of vertices inside hyper-pyramide (the rest outside)
#define WARP_EPS 1.0E-6  // minimum mass difference
#define WARP_EPS10 10.0 * WARP_EPS
#define NORM_LIMIT 1.0E-180

// #define NEIGHBORs_MAX 500 // the maximal # of neighbors adjacent to x
#define DBLP_EPSILON 1e-6
#define DBLP_sum_EPSILON 1e-12
#define DBLP_OBJ_EPSILON 1e-4 / 2
#define VX_MAXDIM 300 // max nbr of prob vars
#define DBLP_INFINITY 1e7
#define mem_alloc(size, type, source, tag, perm) malloc(size)
#define mem_free(mem_ptr, source, tag) free(mem_ptr)

#define OK 0
#define INCONSISTENT 1
#define INPUT_ERROR 2
#define DETACHED 3
#define ATTACHED 4
#define OUT_OF_MEMORY 5
#define TOO_MANY_CONS 6
#define TOO_MANY_ALTS 7
#define UNLIMITED 8
#define NOT_IMPLEMENTED 9
#define CORRUPTED 10
#define TOO_FEW_ALTS 11
#define ILLEGAL_NODE 12
#define TOO_MANY_STMTS 13
#define NO_FILE 14
#define TOO_NARROW_STMT 15
#define RELEASE_EXPIRED 16
#define CONFIG_MISMATCH 17
#define NOT_INITIATED 18
#define TREE_ERROR 19
#define SLOW_CONVERGENCE 20
#define SYS_CORRUPT 21
#define SOLVER_ABORTED 22

// the following pre-defined signals are only used as internal signals of cpm_solver()
#define SOLUTION_LP -3
#define GSL_LU -4
#define CUT_LIMIT -5
#define DEGENERACY -6
#define ALL_NEGATIVE -7
#define INSUFFICIENT -8
#define NO_T_COL -9
#define CUT_SEPARATION_FAILED -10
#define NEIGHBOR_TIME_LIMIT -11

typedef double h_vector[VX_MAXDIM + 1];
typedef double d_row[VX_MAXDIM + 1];

#include <cstddef>
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <ctime>
#include <cmath>
#include <cassert>
#include <limits>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <future>

#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <fenv.h>

#include <gsl/gsl_matrix.h>      
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_errno.h>
#include "resource_scope.h"

#include "lp_backend.h"
#include "Highs.h"
#include "highs_solver.h"

#define max(a, b) ((a) > (b) ? (a) : (b))

#define MAT_SET(M, rows, cols, r, c, val)                \
    do                                                   \
    {                                                    \
        (M)[(size_t)(r) * (cols) + (size_t)(c)] = (val); \
    } while (0)

#define MAT_ADD(M, rows, cols, r, c, val)                 \
    do                                                    \
    {                                                     \
        (M)[(size_t)(r) * (cols) + (size_t)(c)] += (val); \
    } while (0)

static inline size_t nCk2(size_t n) { return (n > 1) ? (n * (n - 1)) / 2 : 0; }
static inline size_t nSym(size_t n) { return n + nCk2(n); }

/* 1-based packing index for the upper triangle including the diagonal (as in the MATLAB code) */
static inline size_t sym_index1(size_t i1, size_t j1, size_t n)
{
    if (j1 < i1)
    {
        return (size_t)((j1 - 1) * (2 * n - j1 + 2) / 2 + i1 - j1 + 1);
    }
    else
    {
        return (size_t)((i1 - 1) * (2 * n - i1 + 2) / 2 + j1 - i1 + 1);
    }
}

static inline size_t location(size_t j, size_t Ai_begin, size_t Aj_begin, size_t n_alt1)
{
    if (j < n_alt1)
    {
        return (size_t)Ai_begin + j;
    }
    else
    {
        return (size_t)Aj_begin + j - n_alt1;
    }
}

static inline size_t sym_index0(size_t i0, size_t j0, size_t n) { return sym_index1(i0 + 1, j0 + 1, n) - 1; }
static inline size_t xy_index0(size_t i0, size_t j0, size_t nY) { return i0 * nY + j0; }
#ifndef DECIDEIT_ADJACENCY_MULTIPLIER
#define DECIDEIT_ADJACENCY_MULTIPLIER 5
#endif
static inline size_t reverse_limit_override(const char *name, size_t fallback)
{
    const char *value = std::getenv(name);
    if (!value) return fallback;
    if (!*value || *value == '-')
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    char *end = nullptr;
    errno = 0;
    const unsigned long parsed = std::strtoul(value, &end, 10);
    if (errno || *end || parsed == 0)
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return parsed < VX_MAXVER ? static_cast<size_t>(parsed) : VX_MAXVER;
}
static inline size_t adjacency_capacity(size_t dimension)
{
    // Exhaust small degenerate faces exactly.  For larger problems retain a
    // bounded enumeration.  The multiplier can be overridden in a controlled
    // benchmark; the production build uses five neighbors per dimension.
    const size_t multiplier = reverse_limit_override(
        "DECIDEIT_NEIGHBOR_CAP_MULTIPLIER", DECIDEIT_ADJACENCY_MULTIPLIER);
    const size_t scaled = dimension <= 8 &&
        !std::getenv("DECIDEIT_NEIGHBOR_CAP_MULTIPLIER") ? 256 :
        (dimension && multiplier > VX_MAXVER / dimension ? VX_MAXVER :
         multiplier * dimension);
    return scaled < VX_MAXVER ? scaled : VX_MAXVER;
}
static inline size_t reverse_node_capacity(size_t dimension)
{
    const size_t scaled = dimension <= 16 ? 2048 : 7 * dimension;
    return reverse_limit_override("DECIDEIT_REVERSE_NODE_CAP",
                                  scaled < VX_MAXVER ? scaled : VX_MAXVER);
}

struct HighsSolver
{
    Highs highs;
    bool loaded = false;
    bool solved = false;
    // Retain the last certified optimal primal/dual basis.  When a newly
    // appended RLT row is satisfied by that primal, extending the dual by a
    // zero row multiplier preserves a valid reduced-cost certificate even
    // though HiGHS invalidates its public solution status after addRows().
    bool reduced_cost_snapshot_valid = false;
    double reduced_cost_objective =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<double> reduced_cost_primal;
    std::vector<double> reduced_cost_dual;
    HighsBasis reduced_cost_basis;
    int loaded_rows = 0;
    int cold_solves = 0;
    int warm_solves = 0;
};

struct DecideItDirectionStats
{
    int cuts = 0;
    int termination_flag = 0;
    std::string certificate = "not_recorded";
    double residual_bound = std::numeric_limits<double>::quiet_NaN();
    double rlt_primal_infeasibility = std::numeric_limits<double>::quiet_NaN();
    double rlt_dual_infeasibility = std::numeric_limits<double>::quiet_NaN();
    int amc_rounds = 0;
    int block_lp_solves = 0;
    int neighbor_lp_solves = 0;
    int neighborhood_calls = 0;
    int climb_restarts = 0;
    int rlt_auxiliary_solves = 0;
    int lower_bound_calls = 0;
    int rlt_skipped_solves = 0;
    int cold_solves = 0;
    int warm_solves = 0;
    int rlt_pivotcopy_solves = 0;
    long long rlt_pivotcopy_pivots = 0;
    int cross_level_basis_uses = 0;
    int cross_level_primal_uses = 0;
    double cross_level_primal_scale_sum = 0.0;
    long long simplex_iterations = 0;
    int max_rows = 0;
    int columns = 0;
    double rlt_auxiliary_s = 0.0;
    double lower_bound_s = 0.0;
};

struct RltCrossLevelState
{
    bool valid = false;
    bool basis_valid = false;
    int rows = 0;
    int columns = 0;
    HighsBasis basis;
    std::vector<double> solution;
    std::vector<double> col_lower;
    std::vector<double> col_upper;
};

struct DecideItRunStats
{
    DecideItDirectionStats minimum;
    DecideItDirectionStats maximum;
};

struct RowReductionStats
{
    int input_rows = 0;
    int retained_rows = 0;
    int cheap_removed = 0;
    int cheap_kept = 0;
    int conservative_kept = 0;
    int witness_kept = 0;
    int lp_calls = 0;
    long long simplex_iterations = 0;
    double wall_s = 0.0;
};

void reset_decideit_stats();
void set_decideit_stats(const DecideItRunStats &stats);
const DecideItRunStats &last_decideit_stats();
std::string selected_rlt_backend();
void reset_row_reduction_stats();
const RowReductionStats &last_row_reduction_stats();

struct LinearizationParams
{
    LinearizationParams();
    ~LinearizationParams();
    int m_total;
    int n_total;

    int n_XY;
    int n_XX;
    int n_YY;

    int num_col;
    int num_row;
    int num_nz;
    double offset;

    gsl_vector *col_cost;

    gsl_vector *row_lower;
    gsl_vector *row_upper;
    gsl_vector *a_value;
    gsl_vector_int *a_index;
    gsl_vector_int *a_start;

    gsl_vector *col_lower;
    gsl_vector *col_upper;

    int a_start_column;
    int a_index_column;
    // For each incremental X cut, record the Y source-row products actually
    // present in the model.  Spatial branch children replay this plan in the
    // same row order so the parent's HiGHS basis remains compatible.
    std::vector<std::vector<int> > cut_cross_rows;
    std::unique_ptr<HighsSolver> persistent;
};

void ensure_rlt_capacity(LinearizationParams*, size_t rows, size_t nonzeros);
void construct_LinearizationParameters(int, int, int, int, LinearizationParams*);
void destroy_LinearizationParameters(LinearizationParams*);
void construct_linearization_objective(int, int, LinearizationParams*);

struct LinCons
{
    /* Sparse row sum(value[k] * z[index[k]]) <= b.  Indices are 1-based.
       The old implementation stored only z[j1]-z[j2] <= b; using a sparse
       row here also covers the compound Delta statements in data/B. */
    std::vector<int> index;
    std::vector<double> coefficient;
    double b = 0.0;
};

typedef struct link_rec
{
    int head;
    int tail;
    char cmp;
    double value;
} link_rec;

typedef struct link_rec link_set[5 * VX_MAXDIM + 1];
void read_delta(const char*, std::vector<int>&, d_row, d_row, link_set,
                std::vector<LinCons>&, d_row, d_row, link_set,
                std::vector<LinCons>&);
LinCons difference_row(int head, int tail, double rhs);
int probability_focal(int n,h_vector lo,h_vector mid,h_vector hi);

struct Solution_Pool
{
    REAL vertex_obj;
    REAL adj_vertex_obj;
    int num_adj;
    int out_index;
    gsl_vector *vertex_x;
    gsl_vector *vertex_y;
    gsl_matrix *adj_vertices;
};

struct dataset
{
    gsl_matrix *valid_A;
    gsl_vector *xl;
    gsl_vector *xu;
    gsl_matrix *Q;
    gsl_vector *dy;
    int ori_num_cons;
};

struct raw_data
{
    gsl_matrix *Ax;
    gsl_vector *bx;
    gsl_vector *xl;
    gsl_vector *xu;
    gsl_matrix *Ay;
    gsl_vector *by;
    gsl_vector *yl;
    gsl_vector *yu;
    gsl_matrix *Q;
};

struct dblp
{
    gsl_matrix *Ax;
    gsl_vector *bx;
    gsl_vector *xl;
    gsl_vector *xu;
    gsl_matrix *Q;
};

// Optional exact search state for one geometric neighbor.  Coordinates alone
// do not identify a basis on an ill-conditioned degenerate face; keeping the
// pivoted dictionary lets a later boundary search resume from the state that
// actually produced the neighbor.  These records are populated only under
// DECIDEIT_STORE_NEIGHBOR_STATE, so the production chain pays no tableau-
// allocation or copy cost.
struct NeighborDictionaryState
{
    int rows = 0;
    int columns = 0;
    std::vector<double> tableau;
    std::vector<int> basic_index;
    std::vector<int> nonbasic_index;
    std::vector<double> point;
};

// Certified first-layer geometry retained only by the opt-in cone pilot.
// rays are unit directions stored row-major (ray, coordinate); gauges are
// reciprocal safe forward radii and are zero on recession rays.
struct ConeRayGeometry
{
    bool valid = false;
    int dimension = 0;
    int ray_count = 0;
    double alpha = 0.0;
    std::vector<double> apex;
    std::vector<double> rays;
    std::vector<double> gauges;
};

struct local_solution
{
    int n_adj_vertices;
    bool degenerate_vertex;
    gsl_matrix *adj_vertices;
    gsl_vector *x_solution;
    gsl_vector *y_solution;
    double cbov;
    std::vector<NeighborDictionaryState> neighbor_states;
    ConeRayGeometry cone_geometry;
};

struct global_solution
{
    gsl_vector *x_solution;
    gsl_vector *y_solution;
    double alpha;
};

void show_gsl_sparse_matrix(gsl_matrix *p);

bool DELTA_data(
    struct dblp *original,
    char *file_name,
    int dim_Alt1,
    int dim_Alt2,
    int num_P,
    int num_V);

bool data(struct dblp *a, char *b);
void show_matrix(REAL **p, int rows, int columns, int indicator);

int local_phase(
    struct local_solution *local_amc,
    LpModel *lp_X_std,
    LpModel *lp_Y,
    struct dataset *dataset_in,
    int noc,
    int n_cuts);

void show_gsl_matrix_columns(gsl_matrix *p, int columns);
void show_gsl_matrix(gsl_matrix *p);
void show_gsl_matrix_int(gsl_matrix_int *p);
void show_gsl_vector(gsl_vector *p);
void show_array_int(int *p);
void show_array_int_all(int *p, int bound);
void show_array_float(REAL *p);

// void dblp_constructor(struct dblp *p, int *dimension);
// void dblp_destructor(struct dblp *p);
// void local_solution_constructor(struct local_solution *p, int *dimension, int CUTS_MAX);
// void local_solution_destructor(struct local_solution *p);
// void global_solution_constructor(struct global_solution *p, int *dimension);
// void global_solution_destructor(struct global_solution *p);
// bool adj_vertices_normal(gsl_matrix *adj_vertices, gsl_matrix *tableau, int *basic, int *non_basic, int *var_status_new, gsl_vector *xl, gsl_vector *xu);

int global_phase(
    struct local_solution *local_amc,
    double alpha,
    LpModel *lp_X_std,
    LpModel *lp_Y,
    struct dataset *dataset_in,
    int noc);

// points contains all adjacent vertices in rows and the current vertex last.
// A successful cut retains every neighbor and strictly removes the last row.
bool conservative_neighbor_cut(gsl_matrix *points, gsl_matrix *valid_A,
                               int valid_rows, const gsl_vector *xl,
                               const gsl_vector *xu, gsl_vector *normal,
                               double *rhs);
// Experimental, opt-in true-cone positive-ray parent cut.  Unlike the
// conservative cut, it may remove locally dominated neighboring vertices.
bool positive_ray_parent_cut(gsl_matrix *points, LpModel *lp_y,
                             const gsl_matrix *Q, const gsl_vector *dy,
                             gsl_vector *normal, double *rhs,
                             ConeRayGeometry *geometry = nullptr);

void show_cutting_values(gsl_vector *b, struct dataset *dataset_in, LpModel *lp_Y, int cuts_indicator, int No);
void show_gsl_vector_int(gsl_vector_int *p);
void show_gsl_matrix_cuts(gsl_matrix *p, int n_cuts);
void check_objective(gsl_vector *x, gsl_vector *y, struct dblp *instance);
void resolve_zero_matrix(gsl_matrix *m);
void resolve_zero_vector(gsl_vector *v);

int reverse_search(gsl_matrix *tableau,
                   gsl_vector_int *basic_index,
                   gsl_vector_int *nonbasic_index,
                   gsl_matrix *&adj_vertices,
                   std::vector<NeighborDictionaryState> *neighbor_states = nullptr);
void grow_neighbor_matrix(gsl_matrix *&neighbors, size_t columns);

void show_gsl_matrix_neighbors(gsl_matrix *neighbors, size_t num_neighbors);

int generate_mp(
    int n_cons,
    double target,
    h_vector lobo,
    h_vector mid,
    h_vector upbo);

bool consistency_check(
    struct dataset *dataset_in,
    gsl_vector *vertex_x,
    int noc);

int cpm_solver(
    d_row P_lobo,
    d_row P_upbo,
    link_set P_links,
    d_row V_lobo,
    d_row V_upbo,
    link_set V_links,
    d_row P_min_point,
    d_row V_min_point,
    d_row P_max_point,
    d_row V_max_point,
    int Ai_begin,
    int Ai_end,
    int Aj_begin,
    int Aj_end,
    int CUTS_MAX,
    std::vector<LinCons> &LC_V,
    std::vector<LinCons> &LV_P,
    int *mX,
    int *mY);

int cpm_solver_rows(
    d_row P_lobo,
    d_row P_upbo,
    const std::vector<LinCons> &P_rows,
    d_row V_lobo,
    d_row V_upbo,
    const std::vector<LinCons> &V_rows,
    d_row P_min_point,
    d_row V_min_point,
    d_row P_max_point,
    d_row V_max_point,
    int Ai_begin,
    int Ai_end,
    int Aj_begin,
    int Aj_end,
    int CUTS_MAX,
    std::vector<LinCons> &LC_V,
    std::vector<LinCons> &LC_P,
    int *mX,
    int *mY,
    const double *P_min_start = nullptr,
    const double *P_max_start = nullptr,
    RltCrossLevelState *min_rlt_state = nullptr,
    RltCrossLevelState *max_rlt_state = nullptr);

void show_gsl_matrix_original(gsl_matrix *p);
void fpe_handler(int signum);
int matrix_rank(const gsl_matrix *A, double tol);
void show_gsl_matrix_columns_int(gsl_matrix_int *p, int columns);
void show_gsl_vector_int_all(gsl_vector_int *p);
void show_gsl_matrix_int_rows(gsl_matrix_int *p, int rows);
void show_LC(LinCons *LC, int num);

// gsl_matrix *construct_constraints(struct dblp *instance, gsl_matrix *cuts, gsl_vector_int *effective_constraints, gsl_vector_int *effective_cuts, gsl_matrix *valid_A);
unsigned long long nchoosek(int n, int k);
void nchoosek_enum(int n, int k, int i, std::vector<int> &result);
void setDifference(std::vector<int> &set1, std::vector<int> &set2, std::vector<int> &set3);
bool condition_number(gsl_matrix *A, double threshold);

int cpm_controller(
    struct dblp *instance,
    struct local_solution *local_amc,
    LpModel *lp_X_std, LpModel *lp_Y,
    struct dataset *dataset_in_use,
    int *noc,
    d_row P_bo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    struct global_solution *solution,
    int CUTS_MAX,
    d_row P_point, d_row V_point,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    gsl_vector_int *fixed_variables,
    std::vector<LinCons> &LC_V, std::vector<LinCons> &LC_P,
    const std::vector<LinCons> &original_P_rows,
    int mX, int mY, LinearizationParams *LinPara,
    bool minFlag,
    const gsl_vector *initial_y = nullptr,
    RltCrossLevelState *cross_level_rlt_state = nullptr);

void construct_complete_solution(gsl_vector *vector_new,
                                 gsl_vector_int *basic_index,
                                 gsl_vector *column_solution);

void show_array_int(int my_array[]);

bool check_consistency_reverse_search(gsl_matrix *A,
                                      gsl_vector *b,
                                      gsl_vector *solution);

void identify_lu_indicator(gsl_matrix *A,
                           gsl_vector *b,
                           gsl_vector *solution,
                           gsl_vector_int *basic_index,
                           gsl_vector_int *nonbasic_index);

int results_BARON(d_row p_point, d_row v_point);

int svd_solve_safe(gsl_matrix *A, gsl_vector *b, gsl_vector *x, double tol);
int safe_inverse(const gsl_matrix *A, gsl_matrix *A_inv, double tol);

int heuristic_solver(
    d_row P_lobo,
    d_row P_upbo,
    link_set P_links,
    d_row V_lobo,
    d_row V_upbo,
    link_set V_links,
    d_row P_min_point,
    d_row V_min_point,
    d_row P_max_point,
    d_row V_max_point,
    int Ai_begin,
    int Ai_end,
    int Aj_begin,
    int Aj_end,
    int max_cuts);

void DataPreparation(
    const std::vector<LinCons> &V_rows,
    gsl_vector_int *effective_constraints,
    const std::vector<LinCons> &P_rows,
    d_row P_lobo,
    d_row P_upbo,
    int Ai_begin,
    int Ai_end,
    int Aj_begin,
    int Aj_end,
    std::vector<LinCons> &LC_V,
    std::vector<LinCons> &LC_P,
    int *mX,
    int *mY);

int redundancy(struct dblp *instance,
               gsl_vector_int *eff_cons,
               gsl_vector_int *fixed_variables);

void show_gsl_vector_counter(gsl_vector *p, int start, int end);
void show_gsl_vector_int_counter(gsl_vector_int *p, int start, int end);

void linearization(
    d_row P_lobo, d_row P_upbo, d_row V_lobo, d_row V_upbo,
    d_row P_point, d_row V_point,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    bool *successFlag,
    LinCons *LC_V, LinCons *LC_P, int mX, int mY, LinearizationParams *LinPara);

// void global_solution_destructor(struct global_solution *p);
// void global_solution_constructor(
//     struct global_solution *p,
//     int dim_x,
//     int dim_y);

// void local_solution_destructor(struct local_solution *p);
// void local_solution_constructor(
//     struct local_solution *p,
//     int dim_x,
//     int dim_y,
//     int CUTS_MAX);

// void dblp_destructor(struct dblp *p);
// void dblp_constructor(
//     struct dblp *p,
//     int *dimension);

// void dataset_destructor(struct dataset *p);
// void dataset_constructor(
//     struct dataset *p,
//     int noc,
//     int nov_x,
//     int nov_y,
//     int CUTS_MAX);

#endif
