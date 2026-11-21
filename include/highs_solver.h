#ifndef HIGHS_SOLVER_H
#define HIGHS_SOLVER_H

#include "Highs.h"

#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct HighsSolver HighsSolver;

    HighsSolver *createHighsSolver();

    void destroyHighsSolver(HighsSolver *solver);

    void setLP(HighsSolver *solver,
               int num_col, int num_row,
               const int *start, const int *index, const double *value,
               const double *col_cost, const double *col_lower, const double *col_upper,
               const double *row_lower, const double *row_upper);

    void setWarmStart(HighsSolver *solver, const double *x0, int num_col);

    void solve_with_highs(HighsSolver *solver, bool outputFlag);

    void getSolution(HighsSolver *solver, double *x_out, int num_col);
    
    double getObjective(HighsSolver *solver);

#ifdef __cplusplus
}
#endif

#endif
