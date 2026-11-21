#include "dblp.h"
#include "decision_path.h"
#ifdef DECIDEIT_WITH_CPLEX
#include "cplex_reference.h"
#endif

// void full_api()
// {
//     printf("\nHiGHS version %s\n", Highs_version());
//     printf("      Major version %" HIGHSINT_FORMAT "\n", Highs_versionMajor());
//     printf("      Minor version %" HIGHSINT_FORMAT "\n", Highs_versionMinor());
//     printf("      Patch version %" HIGHSINT_FORMAT "\n", Highs_versionPatch());
//     printf("      Githash %s\n", Highs_githash());

//     // Compilation date is deprecated.
//     // printf("      compilation date %s\n", Highs_compilationDate());

//     // This illustrates the use of Highs_lpCall, the simple C interface to
//     // HiGHS. It's designed to solve the general LP problem
//     //
//     // Min c^Tx + d subject to L <= Ax <= U; l <= x <= u
//     //
//     // where A is a matrix with m rows and n columns
//     //
//     // The scalar n is num_col
//     // The scalar m is num_row
//     //
//     // The vector c is col_cost
//     // The scalar d is offset
//     // The vector l is col_lower
//     // The vector u is col_upper
//     // The vector L is row_lower
//     // The vector U is row_upper
//     //
//     // The matrix A is represented in packed vector form, either
//     // row-wise or column-wise: only its nonzeros are stored
//     //
//     // * The number of nonzeros in A is num_nz
//     //
//     // * The indices of the nonnzeros in the vectors of A are stored in a_index
//     //
//     // * The values of the nonnzeros in the vectors of A are stored in a_value
//     //
//     // * The position in a_index/a_value of the index/value of the first
//     // nonzero in each vector is stored in a_start
//     //
//     // Note that a_start[0] must be zero
//     //
//     // After a successful call to Highs_lpCall, the primal and dual
//     // solution, and the simplex basis are returned as follows
//     //
//     // The vector x is col_value
//     // The vector Ax is row_value
//     // The vector of dual values for the variables x is col_dual
//     // The vector of dual values for the variables Ax is row_dual
//     // The basic/nonbasic status of the variables x is col_basis_status
//     // The basic/nonbasic status of the variables Ax is row_basis_status
//     //
//     // The status of the solution obtained is model_status

//     // This example does exactly the same as the minimal example above,
//     // but illustrates the full C API.  It first forms and solves the LP
//     //
//     // Min    f  =  x_0 +  x_1 + 3
//     // s.t.                x_1 <= 7
//     //        5 <=  x_0 + 2x_1 <= 15
//     //        6 <= 3x_0 + 2x_1
//     // 0 <= x_0 <= 4; 1 <= x_1
//     //
//     // It then solves it as a maximization, then as a MIP.
//     //

//     const HighsInt num_col = 2;
//     const HighsInt num_row = 3;
//     const HighsInt num_nz = 5;

//     // Define the optimization sense and objective offset
//     HighsInt sense = kHighsObjSenseMinimize;
//     const double offset = 3;

//     // Define the column costs, lower bounds and upper bounds
//     const double col_cost[2] = {1.0, 1.0};
//     const double col_lower[2] = {0.0, 1.0};
//     const double col_upper[2] = {4.0, 1.0e30};

//     // Define the row lower bounds and upper bounds
//     const double row_lower[3] = {-1.0e30, 5.0, 6.0};
//     const double row_upper[3] = {7.0, 15.0, 1.0e30};

//     // Define the constraint matrix column-wise
//     const HighsInt a_format = kHighsMatrixFormatColwise;
//     const HighsInt a_start[2] = {0, 2};
//     const HighsInt a_index[5] = {1, 2, 0, 1, 2};
//     const double a_value[5] = {1.0, 3.0, 1.0, 2.0, 2.0};

//     HighsInt run_status;
//     HighsInt model_status;
//     double objective_function_value;
//     HighsInt simplex_iteration_count;
//     // int64_t mip_node_count;
//     HighsInt primal_solution_status;
//     HighsInt dual_solution_status;
//     HighsInt basis_validity;

//     // Create a Highs instance
//     void *highs = Highs_create();

//     // Pass the LP to HiGHS
//     run_status = Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense,
//                             offset, col_cost, col_lower, col_upper, row_lower,
//                             row_upper, a_start, a_index, a_value);
//     assert(run_status == kHighsStatusOk);

//     // Solve the incumbent model
//     run_status = Highs_run(highs);

//     // The run must be successful
//     assert(run_status == kHighsStatusOk);

//     // Get the model status - which must be optimal
//     model_status = Highs_getModelStatus(highs);
//     assert(model_status == kHighsModelStatusOptimal);

//     printf("Run status = %" HIGHSINT_FORMAT "; Model status = %" HIGHSINT_FORMAT "\n", run_status, model_status);

//     // Get scalar information about the solution
//     Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
//     Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
//     Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
//     Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
//     Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);

//     // The primal and dual solution status values should indicate feasibility
//     assert(primal_solution_status == kHighsSolutionStatusFeasible);
//     assert(dual_solution_status == kHighsSolutionStatusFeasible);

//     // The basis status should indicate validity
//     assert(basis_validity == kHighsBasisValidityValid);

//     double *col_value = (double *)malloc(sizeof(double) * num_col);
//     double *col_dual = (double *)malloc(sizeof(double) * num_col);
//     double *row_value = (double *)malloc(sizeof(double) * num_row);
//     double *row_dual = (double *)malloc(sizeof(double) * num_row);

//     HighsInt *col_basis_status = (HighsInt *)malloc(sizeof(HighsInt) * num_col);
//     HighsInt *row_basis_status = (HighsInt *)malloc(sizeof(HighsInt) * num_row);

//     // Get the primal and dual solution
//     Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);

//     // Get the basis
//     Highs_getBasis(highs, col_basis_status, row_basis_status);

//     // Report the column primal and dual values, and basis status
//     for (HighsInt i = 0; i < num_col; i++)
//     {
//         printf("Col%" HIGHSINT_FORMAT " = %lf; dual = %lf; status = %" HIGHSINT_FORMAT "; \n",
//             i, col_value[i], col_dual[i], col_basis_status[i]);
//     }

//     // Report the row primal and dual values, and basis status
//     for (HighsInt i = 0; i < num_row; i++)
//     {
//         printf("Row%" HIGHSINT_FORMAT " = %lf; dual = %lf; status = %" HIGHSINT_FORMAT "; \n",
//             i, row_value[i], row_dual[i], row_basis_status[i]);
//     }

//     printf("Objective value = %g; Iteration count = %" HIGHSINT_FORMAT "\n", objective_function_value, simplex_iteration_count);

//     // Illustrate extraction of model data
//     HighsInt check_sense;
//     run_status = Highs_getObjectiveSense(highs, &check_sense);
//     assert(run_status == 0);
//     printf("LP problem has objective sense = %" HIGHSINT_FORMAT "\n", check_sense);
//     assert(check_sense == sense);

//     // Illustrate change of model data
//     Highs_changeObjectiveSense(highs, kHighsObjSenseMaximize);

//     // Solve the incumbent model
//     run_status = Highs_run(highs);
//     // The run must be successful
//     assert(run_status == kHighsStatusOk);
//     // Get the model status - which must be optimal
//     model_status = Highs_getModelStatus(highs);
//     assert(model_status == kHighsModelStatusOptimal);

//     printf("Run status = %" HIGHSINT_FORMAT "; Model status = %" HIGHSINT_FORMAT "\n", run_status, model_status);

//     // Get scalar information about the solution
//     Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
//     Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
//     Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
//     Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
//     Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);

//     // The primal and dual solution status values should indicate feasibility
//     assert(primal_solution_status == kHighsSolutionStatusFeasible);
//     assert(dual_solution_status == kHighsSolutionStatusFeasible);
//     // The basis status should indicate validity
//     assert(basis_validity == kHighsBasisValidityValid);

//     // Get the primal and dual solution
//     Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);
//     // Get the basis
//     Highs_getBasis(highs, col_basis_status, row_basis_status);

//     // Report the column primal and dual values, and basis status
//     for (HighsInt i = 0; i < num_col; i++)
//     {
//         printf("Col%" HIGHSINT_FORMAT " = %lf; dual = %lf; status = %" HIGHSINT_FORMAT "; \n",
//             i, col_value[i], col_dual[i], col_basis_status[i]);
//     }

//     // Report the row primal and dual values, and basis status
//     for (HighsInt i = 0; i < num_row; i++)
//     {
//         printf("Row%" HIGHSINT_FORMAT " = %lf; dual = %lf; status = %" HIGHSINT_FORMAT "; \n",
//             i, row_value[i], row_dual[i], row_basis_status[i]);
//     }

//     printf("Objective value = %g; Iteration count = %" HIGHSINT_FORMAT "\n", objective_function_value, simplex_iteration_count);

//     // Now illustrate how LPs can be built within HiGHS by constructing
//     // the same maximization problem. First clear the incumbent model
//     Highs_clearModel(highs);

//     // Demonstrate that the incumbent model is empty
//     const HighsInt check_num_col = Highs_getNumCol(highs);
//     const HighsInt check_num_row = Highs_getNumRow(highs);
//     const HighsInt check_num_nz = Highs_getNumNz(highs);

//     assert(check_num_col == 0);
//     assert(check_num_row == 0);
//     assert(check_num_nz == 0);

//     printf("\nCleared model has %" HIGHSINT_FORMAT " columns, %" HIGHSINT_FORMAT " rows and %" HIGHSINT_FORMAT " nonzeros\n",
//         check_num_col, check_num_row, check_num_nz);

//     // This example does exactly the same as the minimal example above. It forms and solves the LP

//     // Min    f  =  x_0 +  x_1 + 3
//     // s.t.                x_1 <= 7
//     //        5 <=  x_0 + 2x_1 <= 15
//     //        6 <= 3x_0 + 2x_1
//     // 0 <= x_0 <= 4; 1 <= x_1

//     // Define the constraint matrix row-wise, as it is added to the LP with the rows
//     const HighsInt ar_start[3] = {0, 1, 3};
//     const HighsInt ar_index[5] = {1, 0, 1, 0, 1};
//     const double ar_value[5] = {1.0, 1.0, 2.0, 3.0, 2.0};

//     // Add two columns to the empty LP
//     run_status = Highs_addCols(highs, num_col, col_cost, col_lower, col_upper, 0, NULL, NULL, NULL);
//     assert(run_status == 0);

//     // Add three rows to the 2-column LP
//     run_status = Highs_addRows(highs, num_row, row_lower, row_upper, num_nz, ar_start, ar_index, ar_value);
//     assert(run_status == 0);

//     // By default, the optimization sense is minimization, and the objective offset is zero, so these need to be changed
//     Highs_changeObjectiveSense(highs, kHighsObjSenseMaximize);
//     Highs_changeObjectiveOffset(highs, offset);

//     // Illustrate extraction of model data
//     run_status = Highs_getObjectiveSense(highs, &check_sense);
//     assert(run_status == 0);
//     printf("LP problem has objective sense = %" HIGHSINT_FORMAT "\n", check_sense);
//     assert(check_sense == kHighsObjSenseMaximize);

//     // Illustrate a change of option

//     // Check what type of option value you should provide
//     HighsInt option_type;
//     const char *option_string = "primal_feasibility_tolerance";
//     run_status = Highs_getOptionType(highs, option_string, &option_type);
//     printf("Option %s is of type %" HIGHSINT_FORMAT "\n", option_string, option_type);
//     assert(run_status == kHighsStatusOk);
//     assert(option_type == 2);

//     double primal_feasibility_tolerance;
//     Highs_getDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
//     printf("primal_feasibility_tolerance = %g: setting it to 1e-6\n", primal_feasibility_tolerance);
//     primal_feasibility_tolerance = 1e-6;
//     Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance", primal_feasibility_tolerance);

//     // Illustrate how HiGHS can run quietly
//     Highs_setBoolOptionValue(highs, "output_flag", 0);
//     printf("Running quietly...\n");
//     // Solve the incumbent model
//     run_status = Highs_run(highs);
//     printf("Running loudly...\n");
//     Highs_setBoolOptionValue(highs, "output_flag", 1);
//     assert(run_status == 0);

//     // Get the model status - which must be optimal
//     model_status = Highs_getModelStatus(highs);
//     assert(model_status == kHighsModelStatusOptimal);

//     printf("\nRun status = %" HIGHSINT_FORMAT "; Model status = %" HIGHSINT_FORMAT "\n", run_status, model_status);

//     // Check what type of info values are
//     HighsInt info_type;
//     const char *info_string = "objective_function_value";
//     run_status = Highs_getInfoType(highs, info_string, &info_type);
//     printf("Info %s is of type %" HIGHSINT_FORMAT "\n", info_string, info_type);
//     assert(run_status == kHighsStatusOk);
//     assert(info_type == kHighsInfoTypeDouble);

//     Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
//     Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
//     Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
//     Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
//     Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);

//     // The primal and dual solution status values should indicate feasibility
//     assert(primal_solution_status == kHighsSolutionStatusFeasible);
//     assert(dual_solution_status == kHighsSolutionStatusFeasible);
//     // The basis status should indicate validity
//     assert(basis_validity == kHighsBasisValidityValid);

//     // Get the primal
//     Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);
//     Highs_getBasis(highs, col_basis_status, row_basis_status);

//     // Report the column primal and dual values, and basis status
//     for (HighsInt i = 0; i < num_col; i++)
//     {
//         printf("Col%" HIGHSINT_FORMAT " = %lf; dual = %lf; status = %" HIGHSINT_FORMAT "; \n",
//             i, col_value[i], col_dual[i], col_basis_status[i]);
//     }

//     // Report the row primal and dual values, and basis status
//     for (HighsInt i = 0; i < num_row; i++)
//     {
//         printf("Row%" HIGHSINT_FORMAT " = %lf; dual = %lf; status = %" HIGHSINT_FORMAT "; \n",
//             i, row_value[i], row_dual[i], row_basis_status[i]);
//     }

//     printf("Objective value = %g; Iteration count = %" HIGHSINT_FORMAT "\n", objective_function_value, simplex_iteration_count);

//     free(col_value);
//     free(col_dual);
//     free(row_value);
//     free(row_dual);
//     free(col_basis_status);
//     free(row_basis_status);

//     Highs_destroy(highs);
// }

void final_consistency_check(struct raw_data *instance, gsl_vector *vertex_x)
{
    int i, j, counter;
    REAL value, sum;
    const int m = instance->Ax->size1, n = instance->Ax->size2;

    counter = 0;

    for (i = 0; i < m; i++)
    {
        sum = 0;

        for (j = 0; j < n; j++)
        {
            value = gsl_matrix_get(instance->Ax, i, j);

            if (value != 0)
                sum += value * gsl_vector_get(vertex_x, j);
        }

        if (sum > gsl_vector_get(instance->bx, i) + 1e-4)
        {
            if (counter == 0)
            {
                printf("\n");
            }

            printf("the %dth linear constraint in Ax <= b is violated, lhs = %f, rhs = %f\n",
                   i, sum, gsl_vector_get(instance->bx, i));

            counter++;
        }
    }

    counter = 0;

    for (j = 0; j < n; j++)
    {
        value = gsl_vector_get(vertex_x, j);

        if ((value < gsl_vector_get(instance->xl, j) - 1e-6) || (value > gsl_vector_get(instance->xu, j) + 1e-6))
        {
            if (counter == 0)
            {
                printf("\n");
            }

            printf("the box constraints for variable %d is violated, v[%d]_l = %f, v[%d]_u = %f, v[%d] = %f\n", j + 1, j + 1, gsl_vector_get(instance->xl, j), j + 1, gsl_vector_get(instance->xu, j), j + 1, gsl_vector_get(vertex_x, j));

            counter++;
        }
    }
}

void rawData_constructor(struct raw_data *p, int *dimension)
{
    p->Ax = managed::matrix_calloc(dimension[2], dimension[3]);
    p->bx = managed::vector_calloc(dimension[2]);
    p->xl = managed::vector_calloc(dimension[3]);
    p->xu = managed::vector_calloc(dimension[3]);
    gsl_vector_set_all(p->xl, 0);
    gsl_vector_set_all(p->xu, 1);

    p->Ay = managed::matrix_calloc(dimension[0], dimension[1]);
    p->by = managed::vector_calloc(dimension[0]);
    p->yl = managed::vector_calloc(dimension[1]);
    p->yu = managed::vector_calloc(dimension[1]);
    gsl_vector_set_all(p->yl, 0);
    gsl_vector_set_all(p->yu, 1);

    p->Q = managed::matrix_calloc(dimension[1], dimension[3]);
}

void rawData_destructor(struct raw_data *p)
{
    managed::matrix_free(p->Ax);
    managed::vector_free(p->bx);
    managed::vector_free(p->xl);
    managed::vector_free(p->xu);
    managed::matrix_free(p->Ay);
    managed::vector_free(p->by);
    managed::vector_free(p->yl);
    managed::vector_free(p->yu);
    managed::matrix_free(p->Q);
    mem_free(p, "exit_cpm", 1);
}

void duplicate_rawData(
    struct raw_data *instance_copy,
    struct raw_data *instance)
{
    gsl_matrix_memcpy(instance_copy->Q, instance->Q);

    gsl_matrix_memcpy(instance_copy->Ax, instance->Ax);
    gsl_vector_memcpy(instance_copy->bx, instance->bx);
    gsl_vector_memcpy(instance_copy->xl, instance->xl);
    gsl_vector_memcpy(instance_copy->xu, instance->xu);

    gsl_matrix_memcpy(instance_copy->Ay, instance->Ay);
    gsl_vector_memcpy(instance_copy->by, instance->by);
    gsl_vector_memcpy(instance_copy->yl, instance->yl);
    gsl_vector_memcpy(instance_copy->yu, instance->yu);
}

/*
    get_FocalPoint() computes the focal points of P-/V-base, respectively.

    inputs
        instance: the structure saving the original data;
        dim_Alt1, dim_Alt2: the dimensions of alternative 1 and 2, respectively;
        focal_P, focal_V: the arrays to save the focal points of P-/V-base, respectively.

    output
        focal_P, focal_V: the arrays saving the focal points of P-/V-base, respectively.
*/

void get_FocalPoint(struct raw_data *instance, int dim_Alt1, int dim_Alt2, h_vector focal_P, h_vector focal_V)
{
    int i, j;
    REAL sum;

    h_vector mid_alt1, mid_alt2;
    h_vector lobo_alt1, lobo_alt2;
    h_vector upbo_alt1, upbo_alt2;

    lobo_alt1[0] = dim_Alt1;
    mid_alt1[0] = dim_Alt1;
    upbo_alt1[0] = dim_Alt1;

    for (i = 1; i <= dim_Alt1; i++)
    {
        lobo_alt1[i] = gsl_vector_get(instance->yl, i - 1);
        mid_alt1[i] = 0;
        upbo_alt1[i] = gsl_vector_get(instance->yu, i - 1);
    }

    // show_array_float(lobo_alt1);
    // show_array_float(upbo_alt1);

    if (probability_focal(dim_Alt1, lobo_alt1, mid_alt1, upbo_alt1)!=dim_Alt1)
        throw std::runtime_error("Cannot construct probability focal point for alternative 1");

    sum = 0;

    for (i = 1; i <= dim_Alt1; i++)
    {
        sum = sum + mid_alt1[i];
    }

    // show_array_float(mid_alt1);
    // printf("the sum of %d figures is %f\n", (int)mid_alt1[0], sum);

    lobo_alt2[0] = dim_Alt2;
    mid_alt2[0] = dim_Alt2;
    upbo_alt2[0] = dim_Alt2;

    for (i = dim_Alt1 + 1; i <= dim_Alt1 + dim_Alt2; i++)
    {
        lobo_alt2[i - dim_Alt1] = gsl_vector_get(instance->yl, i - 1);
        mid_alt2[i - dim_Alt1] = 0;
        upbo_alt2[i - dim_Alt1] = gsl_vector_get(instance->yu, i - 1);
    }

    // show_array_float(lobo_alt2);
    // show_array_float(upbo_alt2);

    if (probability_focal(dim_Alt2, lobo_alt2, mid_alt2, upbo_alt2)!=dim_Alt2)
        throw std::runtime_error("Cannot construct probability focal point for alternative 2");

    sum = 0;

    for (i = 1; i <= dim_Alt2; i++)
    {
        sum = sum + mid_alt2[i];
    }

    // show_array_float(mid_alt2);
    // printf("the sum of %d figures is %f\n", (int)mid_alt2[0], sum);

    focal_P[0] = mid_alt1[0] + mid_alt2[0];

    for (i = 1; i <= mid_alt1[0]; i++)
    {
        focal_P[i] = mid_alt1[i];
    }

    for (j = 1; j <= mid_alt2[0]; j++)
    {
        focal_P[i] = mid_alt2[j];
        i++;
    }

    focal_V[0] = focal_P[0];

    for (i = 1; i <= focal_V[0]; i++)
    {
        focal_V[i] = (gsl_vector_get(instance->xl, i - 1) + gsl_vector_get(instance->xu, i - 1)) / 2;
    }

    // show_array_float(focal_P);
    // show_array_float(focal_V);
}

/*
    adjust_Constraints() adjusts the rhs of linear constraints in V-base and P-base with their focal points.

    inputs
        instance: the structure saving the original data and to save the updated data;
        dim_Alt1, dim_Alt2: the dimensions of alternative 1 and 2, respectively;
        focal_Alt1, focal_Alt2: the arrays to save the focal points of alternative 1 and 2, respectively.

    output
        instance: the structure saving the updated data.
*/

void adjust_constraints(struct raw_data *instance, int, int,
                        double *mid_VBase, double *mid_PBase)
{
    int i, j;
    double sum;

    // show_array_float(mid_VBase);

    for (i = 0; i < instance->Ax->size1; i++) // adjust the right-hand sides of constraints in V-Base
    {
        sum = 0;
        // printf("\n");

        for (j = 0; j < instance->Ax->size2; j++)
        {
            if (gsl_matrix_get(instance->Ax, i, j) != 0)
            {
                sum = sum + gsl_matrix_get(instance->Ax, i, j) * mid_VBase[j + 1];
                // printf("Ax[%d][%d] = %f, mass[%d] = %f, xl[%d] = %f, xu[%d] = %f\n", i, j, gsl_matrix_get(instance->Ax, i, j), j + 1, mid_VBase[j + 1], j, gsl_vector_get(instance->xl, j), j, gsl_vector_get(instance->xu, j));
            }
        }

        // printf("sum = %f, current right-hand side = %f\n", sum, gsl_vector_get(instance->bx, i));

        if (sum > gsl_vector_get(instance->bx, i))
        {
            // printf("\nbx[%d] is changed from %f to %f\n\n", i, gsl_vector_get(instance->bx, i), sum);
            gsl_vector_set(instance->bx, i, sum);
        }
    }

    // show_array_float(mid_PBase);

    // The final four rows encode the two probability normalizations as pairs
    // of inequalities.  Every earlier row is a genuine (possibly sparse)
    // user constraint and must be adjusted independently of its coefficient
    // sum; a general row may otherwise be mistaken for a normalization row.
    const int probability_statement_rows =
        static_cast<int>(instance->Ay->size1) - 4;
    for (i = 0; i < probability_statement_rows; i++)
    {
        sum = 0;
        for (j = 0; j < instance->Ay->size2; j++)
        {
            if (gsl_matrix_get(instance->Ay, i, j) != 0)
            {
                sum = sum + gsl_matrix_get(instance->Ay, i, j) * mid_PBase[j + 1];
            }
        }
        if (sum > gsl_vector_get(instance->by, i))
        {
            gsl_vector_set(instance->by, i, sum);
        }
    }
}

void clean_results(
    d_row P_min_point,
    d_row V_min_point,
    d_row P_max_point,
    d_row V_max_point,
    int dim)
{
    int i;

    for (i = 1; i <= dim; i++)
    {
        P_min_point[i] = 0;
        P_max_point[i] = 0;
        V_min_point[i] = 0;
        V_max_point[i] = 0;
    }
}

void display_solution(d_row P_min_point, d_row V_min_point, d_row P_max_point, d_row V_max_point, int Ai_begin, int Ai_end, int Aj_begin, int Aj_end)
{
    int i;
    REAL obj;

    obj = 0;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        obj = obj + P_min_point[i] * V_min_point[i];
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        obj = obj - P_min_point[i] * V_min_point[i];
    }

    printf("\nthe computed min: %f, P_min_point[0] = %f\n", obj, P_min_point[0]);

    if (fabs(obj - P_min_point[0]) > DBLP_EPSILON)
    {
        printf("\n\ntwo values differ up to %f\n\n\n", DBLP_EPSILON);
    }

    printf("\nvariables in P-base (min):\n");

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        printf("% .3lf", P_min_point[i]);
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        printf("% .3lf", P_min_point[i]);
    }

    printf("\n");
    printf("variables in V-base (min):\n");

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        printf("% .3lf", V_min_point[i]);
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        printf("% .3lf", V_min_point[i]);
    }

    printf("\n");

    obj = 0;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        obj = obj + P_max_point[i] * V_max_point[i];
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        obj = obj - P_max_point[i] * V_max_point[i];
    }

    printf("\nthe computed max: %f, P_max_point[0] = %f\n", obj, -P_max_point[0]);

    if (fabs(obj + P_max_point[0]) > DBLP_EPSILON)
    {
        printf("\n\ntwo values differ up to %f\n\n\n", DBLP_EPSILON);
    }

    printf("variables in P-base (max):\n");

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        printf("% .3lf", P_max_point[i]);
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        printf("% .3lf", P_max_point[i]);
    }

    printf("\n");
    printf("variables in V-base (max):\n");

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        printf("% .3lf", V_max_point[i]);
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        printf("% .3lf", V_max_point[i]);
    }

    printf("\n");
}

void output_datafile(struct raw_data *instance_copy, double percentage, char *file_name, int num_alt_i, int num_alt_j)
{
    std::string filename(file_name);
    auto dot=filename.find_last_of('.');
    if(dot!=std::string::npos)filename.resize(dot);
    char suffix[32];snprintf(suffix,sizeof(suffix),"_%03dp.m",static_cast<int>(percentage));
    filename+=suffix;
    int i,j,ceiling,pointer_1=0,pointer_2=0;
    FILE* fid=fopen(filename.c_str(),"w");
    if(!fid)throw std::runtime_error("Cannot create export file: "+filename);
    fprintf(fid," %3d %3d\\n",num_alt_i,num_alt_j);
    fprintf(fid," %3d %3d\\n",static_cast<int>(instance_copy->Ay->size1)-4,static_cast<int>(instance_copy->Ax->size1));
    ceiling = num_alt_i + num_alt_j;

    for (i = 0; i < ceiling; i++)
    {
        fprintf(fid, " %3d: [%.12f %.12f] [%.12f %.12f]\n", i + 1, gsl_vector_get(instance_copy->yl, i), gsl_vector_get(instance_copy->yu, i), gsl_vector_get(instance_copy->xl, i), gsl_vector_get(instance_copy->xu, i));
    }

    ceiling = instance_copy->Ay->size1 - 4;

    for (i = 0; i < ceiling; i++)
    {
        for (j = 0; j < instance_copy->Ay->size2; j++)
        {
            if (gsl_matrix_get(instance_copy->Ay, i, j) == 1)
            {
                pointer_1 = j + 1;
            }

            if (gsl_matrix_get(instance_copy->Ay, i, j) == -1)
            {
                pointer_2 = j + 1;
            }
        }

        fprintf(fid, " %3d: P%d-P%d < %.12f\n", i + 1, pointer_1, pointer_2, gsl_vector_get(instance_copy->by, i));
    }

    ceiling = instance_copy->Ax->size1;

    for (i = 0; i < ceiling; i++)
    {
        for (j = 0; j < instance_copy->Ax->size2; j++)
        {
            if (gsl_matrix_get(instance_copy->Ax, i, j) == 1)
            {
                pointer_1 = j + 1;
            }

            if (gsl_matrix_get(instance_copy->Ax, i, j) == -1)
            {
                pointer_2 = j + 1;
            }
        }

        fprintf(fid, " %3d: V%d-V%d < %.12f\n", i + 1, pointer_1, pointer_2, gsl_vector_get(instance_copy->bx, i));
    }

    fclose(fid);

}

/*
    this program tries to solve dblp specially designed for imprecise decision analysis.

    input
        argv[1]: file path for the input dataset
        argv[2]: the # of contractions
        argv[3]: 0 for debug, and 1 otherwise

    output
        results are in P_point and V_point
*/

static bool selected_cplex()
{
    const char *solver = std::getenv("DECIDEIT_SOLVER");
    return solver && std::string(solver) == "cplex";
}

static std::string selected_method()
{
    const char *value = std::getenv("DECIDEIT_METHOD");
    return value ? value : "bcm";
}

static int selected_trigger()
{
    const char *value = std::getenv("DECIDEIT_TRIGGER");
    return value ? std::stoi(value) : 3;
}

static int selected_cut_limit()
{
    const char *value = std::getenv("DECIDEIT_CUT_LIMIT");
    return value ? std::stoi(value) : 5000;
}

static int selected_bcm_through()
{
    const char *value = std::getenv("DECIDEIT_BCM_THROUGH");
    return value ? std::stoi(value) : -1;
}

static std::string selected_climb()
{
    const char *value = std::getenv("DECIDEIT_CLIMB");
    return value ? value : "guarded";
}

static std::string selected_contraction_start()
{
    const char *value = std::getenv("DECIDEIT_CONTRACTION_START");
    // Map the preceding same-sense solution into the new homothetic box.  The
    // map is a convex combination with the focal point, hence it also
    // preserves normalization and every general linear P constraint.
    return value ? value : "transport";
}

static std::string selected_contraction_path()
{
    const char *value = std::getenv("DECIDEIT_CONTRACTION_PATH");
    // The paper's complete comparison chain includes the uncontracted model:
    // 0,1/K,...,(K-1)/K plus direct evaluation of the focal point.
    return value ? value : "inclusive";
}

static std::string selected_contraction_reduction()
{
    const char *value = std::getenv("DECIDEIT_CONTRACTION_REDUCTION");
    return value ? value : "refresh";
}

static std::string selected_row_filter()
{
    const char *value = std::getenv("DECIDEIT_ROW_FILTER");
    return value ? value : "cheap";
}

static void print_row_reduction_stats(int level)
{
    const RowReductionStats &stats = last_row_reduction_stats();
    printf("REDUCTION {\"level\":%d,\"row_filter\":\"%s\","
           "\"input_rows\":%d,\"retained_rows\":%d,"
           "\"cheap_removed\":%d,\"cheap_kept\":%d,"
           "\"conservative_kept\":%d,"
           "\"witness_kept\":%d,\"lp_calls\":%d,"
           "\"simplex_iterations\":%lld,\"wall_s\":%.9g}\n",
           level, selected_row_filter().c_str(), stats.input_rows,
           stats.retained_rows, stats.cheap_removed, stats.cheap_kept,
           stats.conservative_kept, stats.witness_kept, stats.lp_calls,
           stats.simplex_iterations, stats.wall_s);
}

static void print_lp_backend_stats(int level)
{
    const LpBackendStats stats = last_lp_backend_stats();
    printf("LP_BACKEND {\"level\":%d,\"mode\":\"%s\","
           "\"highs_solves\":%lld,\"highs_s\":%.9g,"
           "\"pivotcopy_solves\":%lld,\"pivotcopy_pivots\":%lld,"
           "\"pivotcopy_warm_solves\":%lld,"
           "\"pivotcopy_warm_retries\":%lld,"
           "\"pivotcopy_failures\":%lld,\"pivotcopy_fallbacks\":%lld,"
           "\"cached_solves\":%lld,"
           "\"pivotcopy_s\":%.9g,"
           "\"max_columns\":%d,\"max_rows\":%d,"
           "\"max_highs_columns\":%d,\"max_highs_rows\":%d,"
           "\"max_pivotcopy_columns\":%d,\"max_pivotcopy_rows\":%d}\n",
           level, selected_lp_backend().c_str(), stats.highs_solves,
           stats.highs_s, stats.pivotcopy_solves, stats.pivotcopy_pivots,
           stats.pivotcopy_warm_solves,
           stats.pivotcopy_warm_retries,
           stats.pivotcopy_failures, stats.pivotcopy_fallbacks,
           stats.cached_solves, stats.pivotcopy_s,
           stats.max_columns, stats.max_rows,
           stats.max_highs_columns, stats.max_highs_rows,
           stats.max_pivotcopy_columns, stats.max_pivotcopy_rows);
}

static void print_rlt_backend_stats(int level)
{
    const DecideItRunStats &stats = last_decideit_stats();
    const char *scale = std::getenv("DECIDEIT_RLT_SIMPLEX_SCALE_STRATEGY");
    printf("RLT_BACKEND {\"level\":%d,\"mode\":\"%s\","
           "\"scale_strategy\":\"%s\","
           "\"pivotcopy_solves\":%d,\"pivotcopy_pivots\":%lld}\n",
           level, selected_rlt_backend().c_str(),
           scale ? scale : "auto",
           stats.minimum.rlt_pivotcopy_solves +
               stats.maximum.rlt_pivotcopy_solves,
           stats.minimum.rlt_pivotcopy_pivots +
               stats.maximum.rlt_pivotcopy_pivots);
}

static std::string selected_contraction_order()
{
    const char *value = std::getenv("DECIDEIT_CONTRACTION_ORDER");
    return value ? value : "inward";
}

static std::string selected_contraction_precheck()
{
    const char *value = std::getenv("DECIDEIT_CONTRACTION_PRECHECK");
    // Contracted levels contain the focal point by construction, so their
    // zero-objective feasibility LP is redundant.  Inclusive level zero is
    // the unmodified user model and keeps the explicit feasibility check.
    return value ? value : "solve";
}

static std::string selected_rlt_start()
{
    const char *value = std::getenv("DECIDEIT_RLT_START");
    // Solve the smaller base RLT once, then add the accumulated cuts and
    // reoptimize from its exact simplex basis.
    return value ? value : "staged";
}

#ifdef DECIDEIT_WITH_CPLEX
static double cplex_time_limit()
{
    const char *value = std::getenv("DECIDEIT_CPLEX_TIMELIMIT");
    if (!value) return 0.0;
    size_t used = 0;
    const double result = std::stod(value, &used);
    if (used != std::strlen(value) || !std::isfinite(result) || result < 0.0)
        throw std::invalid_argument("CPLEX time limit must be a finite nonnegative number");
    return result;
}

static int cplex_threads()
{
    const char *value = std::getenv("DECIDEIT_CPLEX_THREADS");
    if (!value) return 1;
    size_t used = 0;
    const int result = std::stoi(value, &used);
    if (used != std::strlen(value) || result < 1 || result > 1024)
        throw std::invalid_argument("CPLEX threads must be 1..1024");
    return result;
}
#endif

static int solve_selected(
    d_row pl, d_row pu, link_set pc, const std::vector<LinCons> &p_rows,
    d_row vl, d_row vu, link_set vc, const std::vector<LinCons> &v_rows,
    d_row pmin, d_row vmin, d_row pmax, d_row vmax,
    int ai, int ae, int bi, int be, int cuts,
    std::vector<LinCons> &lc_v, std::vector<LinCons> &lc_p, int *mx, int *my,
    const double *pmin_start, const double *pmax_start,
    RltCrossLevelState *min_rlt_state = nullptr,
    RltCrossLevelState *max_rlt_state = nullptr)
{
    if (selected_cplex())
    {
#ifdef DECIDEIT_WITH_CPLEX
        return cplex_solve_pair_rows(
            pl, pu, p_rows, vl, vu, v_rows, pmin, vmin, pmax, vmax,
            ai, ae, bi, be, cplex_time_limit(), cplex_threads(),
            std::getenv("DECIDEIT_CPLEX_LOG") != nullptr);
#else
        throw std::runtime_error(
            "This binary was built without CPLEX; run 'make cplex' and use bin/DELTA_solver_cplex");
#endif
    }
    return cpm_solver_rows(
        pl, pu, p_rows, vl, vu, v_rows, pmin, vmin, pmax, vmax,
        ai, ae, bi, be, cuts, lc_v, lc_p, mx, my,
        pmin_start, pmax_start, min_rlt_state, max_rlt_state);
}

int run_cli(int, char *argv[])
{
    namespace dp = decision_path;
    const auto chain_begin = std::chrono::steady_clock::now();
    auto option = [](const char* name) {const char* x=std::getenv(name);return x?std::string(x):std::string();};
    dp::init_log(option("DECIDEIT_DIAGNOSTICS"));
    dp::stage("input");
    managed::Scope invocation_resources;
    const bool debug=std::string(argv[3])=="1";
    char* file_name=argv[1];
    link_set P_links={},V_links={};
    std::vector<LinCons> P_rows, V_rows;
    d_row P_lobo={},P_upbo={},V_lobo={},V_upbo={},P_min_point={},V_min_point={},P_max_point={},V_max_point={};
    int i,j;
    std::vector<int> dimension;
    read_delta(file_name, dimension, P_lobo, P_upbo, P_links, P_rows,
               V_lobo, V_upbo, V_links, V_rows);
    if (std::string(argv[4]) == "validate")
    {
        printf("VALID {\"alternatives\":2,\"first_dimension\":%d,\"second_dimension\":%d,\"probability_rows\":%zu,\"value_rows\":%zu,\"contraction_compatible\":null,\"validation_scope\":\"syntax_only\"}\n",
               dimension[1], dimension[2], P_rows.size(), V_rows.size());
        return 0;
    }
    int Ai_begin, Ai_end, Aj_begin, Aj_end, exit_flag;

    Ai_begin = 1;
    Ai_end = dimension[1];

    Aj_begin = dimension[1] + 1;
    Aj_end = dimension[1] + dimension[2];

    int new_dimension[4];
    new_dimension[3] = Ai_end - Ai_begin + 1 + Aj_end - Aj_begin + 1;
    new_dimension[2] = static_cast<int>(V_rows.size());
    new_dimension[1] = new_dimension[3];
    new_dimension[0] = static_cast<int>(P_rows.size()) + 4;

    struct raw_data *instance = (struct raw_data *)mem_alloc(sizeof(struct raw_data), "raw_data", "DELTA_solver", 1, FALSE);
    struct raw_data *instance_copy = (struct raw_data *)mem_alloc(sizeof(struct raw_data), "raw_data", "DELTA_solver", 1, FALSE);

    rawData_constructor(instance, new_dimension);
    rawData_constructor(instance_copy, new_dimension);
    std::unique_ptr<raw_data,void(*)(raw_data*)> original_owner(instance,rawData_destructor);
    std::unique_ptr<raw_data,void(*)(raw_data*)> copy_owner(instance_copy,rawData_destructor);

    const int num_alt_i = Ai_end - Ai_begin + 1, num_alt_j = Aj_end - Aj_begin + 1;
    // set Q
    gsl_matrix_set_identity(instance->Q);

    j = num_alt_i;

    for (i = 0; i < num_alt_j; i++)
    {
        gsl_matrix_set(instance->Q, j, j, -1);
        j++;
    }

    // set lower and upper bounds in V-base
    j = 0;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        gsl_vector_set(instance->xl, j, V_lobo[i]);
        gsl_vector_set(instance->xu, j, V_upbo[i]);
        j++;
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        gsl_vector_set(instance->xl, j, V_lobo[i]);
        gsl_vector_set(instance->xu, j, V_upbo[i]);
        j++;
    }

    // set lower and upper bounds in P-base
    j = 0;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        gsl_vector_set(instance->yl, j, P_lobo[i]);
        gsl_vector_set(instance->yu, j, P_upbo[i]);
        j++;
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        gsl_vector_set(instance->yl, j, P_lobo[i]);
        gsl_vector_set(instance->yu, j, P_upbo[i]);
        j++;
    }

    auto set_sparse_rows = [](gsl_matrix *matrix, gsl_vector *rhs,
                              const std::vector<LinCons> &rows)
    {
        for (size_t row = 0; row < rows.size(); ++row)
        {
            for (size_t term = 0; term < rows[row].index.size(); ++term)
                gsl_matrix_set(matrix, row, rows[row].index[term] - 1,
                               rows[row].coefficient[term]);
            gsl_vector_set(rhs, row, rows[row].b);
        }
    };
    set_sparse_rows(instance->Ax, instance->bx, V_rows);
    set_sparse_rows(instance->Ay, instance->by, P_rows);

    // show_gsl_matrix(instance->Ay);
    // set four normalization constraints in P-base
    i = new_dimension[0] - 4;

    for (j = 0; j < num_alt_i; j++)
    {
        gsl_matrix_set(instance->Ay, i, j, 1);
        gsl_matrix_set(instance->Ay, i + 1, j, -1);
    }

    gsl_vector_set(instance->by, i, 1);
    gsl_vector_set(instance->by, i + 1, -1);

    i = i + 2;
    for (j = num_alt_i; j < num_alt_i + num_alt_j; j++)
    {
        gsl_matrix_set(instance->Ay, i, j, 1);
        gsl_matrix_set(instance->Ay, i + 1, j, -1);
    }

    gsl_vector_set(instance->by, i, 1);
    gsl_vector_set(instance->by, i + 1, -1);

    // show_gsl_matrix(instance->Ay);
    // show_gsl_vector(instance->by);

    int mX=-1, mY=-1;
    const int CUTS_MAX = selected_cut_limit();
    clock_t start, end;
    double cpu_time_used;
    std::vector<LinCons> LC_V, LC_P;

    const bool fixed_statements = option("DECIDEIT_CONTRACTION_MODEL") != "legacy";
    const std::string input_hash = dp::file_hash(file_name);
    const std::string original_rows_hash = dp::rows_hash(instance);
    dp::log("\"event\":\"input\",\"input_sha256\":"+dp::quote(input_hash)+
            ",\"original_rows_sha256\":"+dp::quote(original_rows_hash)+
            ",\"first_dimension\":"+std::to_string(num_alt_i)+
            ",\"second_dimension\":"+std::to_string(num_alt_j)+
            ",\"dimension_per_block\":"+std::to_string(num_alt_i+num_alt_j)+
            ",\"contraction_model\":"+dp::quote(fixed_statements?"fixed-statements":"legacy"));
    auto emit_directions = [&](int level,const raw_data* model,int status) {
        double min_bound=last_decideit_stats().minimum.residual_bound;
        double max_bound=-last_decideit_stats().maximum.residual_bound;
        int min_status=0,max_status=0;
#ifdef DECIDEIT_WITH_CPLEX
        if(selected_cplex()) {
            const CplexRunStats& cs=last_cplex_stats();
            min_bound=cs.min_bound;max_bound=cs.max_bound;
            min_status=cs.min_status;max_status=cs.max_status;
        }
#endif
        dp::direction(level,"min",model,num_alt_i,P_min_point,V_min_point,
                      last_decideit_stats().minimum,status,selected_cplex(),min_bound,min_status);
        dp::direction(level,"max",model,num_alt_i,P_max_point,V_max_point,
                      last_decideit_stats().maximum,status,selected_cplex(),max_bound,max_status);
    };
    h_vector focal_V={}, focal_P={};
    dp::FocalInfo focal_info;
    if(std::string(argv[4])=="normal" || !option("DECIDEIT_PREPARE_FOCAL").empty()) {
        dp::stage("focal_preparation");
        get_FocalPoint(instance,num_alt_i,num_alt_j,focal_P,focal_V);
        if(fixed_statements) {
            focal_info=dp::feasible_focal(instance,num_alt_i,num_alt_j,focal_P,focal_V,
                input_hash,option("DECIDEIT_FOCAL_FILE"),option("DECIDEIT_PREPARE_FOCAL"));
        } else {
            focal_info.candidate_violation=dp::point_violation(instance,focal_P,focal_V,num_alt_i);
            focal_info.final_violation=focal_info.candidate_violation;
            focal_info.hash=dp::sha256(dp::array(focal_P,num_alt_i+num_alt_j)+dp::array(focal_V,num_alt_i+num_alt_j));
        }
        dp::log("\"event\":\"focal\",\"focal_sha256\":"+dp::quote(focal_info.hash)+
            ",\"loaded\":"+(focal_info.loaded?"true":"false")+
            ",\"projection_lps\":"+std::to_string(focal_info.projection_lps)+
            ",\"preparation_s\":"+dp::number(focal_info.preparation_s)+
            ",\"candidate_violation\":"+dp::number(focal_info.candidate_violation)+
            ",\"original_model_violation\":"+dp::number(focal_info.final_violation)+
            ",\"p\":"+dp::array(focal_P,num_alt_i+num_alt_j)+
            ",\"v\":"+dp::array(focal_V,num_alt_i+num_alt_j));
        if(!option("DECIDEIT_PREPARE_FOCAL").empty()) return 0;
    }

    if(std::string(argv[4])=="export") {
        output_datafile(instance,0,file_name,num_alt_i,num_alt_j);
        return 0;
    }
    if(std::string(argv[4])=="single") {
        dp::stage("dblp_solve",0);
        auto start_wall=std::chrono::steady_clock::now();
        exit_flag=solve_selected(P_lobo,P_upbo,P_links,P_rows,V_lobo,V_upbo,V_links,V_rows,
            P_min_point,V_min_point,P_max_point,V_max_point,Ai_begin,Ai_end,Aj_begin,Aj_end,
            CUTS_MAX,LC_V,LC_P,&mX,&mY,nullptr,nullptr);
        emit_directions(0,instance,exit_flag);
        if(exit_flag==OK && (dp::point_violation(instance,P_min_point,V_min_point,num_alt_i)>1e-8 ||
                            dp::point_violation(instance,P_max_point,V_max_point,num_alt_i)>1e-8))
            throw std::runtime_error("Returned single-level witness is infeasible");
        if(exit_flag!=OK) {
#ifdef DECIDEIT_WITH_CPLEX
            if (selected_cplex()) {
                const CplexRunStats &cs = last_cplex_stats();
                throw std::runtime_error("CPLEX query not certified: min=" + cs.min_status_text +
                                         ", max=" + cs.max_status_text);
            }
#endif
            const DecideItRunStats &ds = last_decideit_stats();
            throw std::runtime_error(
                "Single query not certified, status=" + std::to_string(exit_flag) +
                ", min_termination=" + std::to_string(ds.minimum.termination_flag) +
                ", max_termination=" + std::to_string(ds.maximum.termination_flag) +
                ", cuts_min=" + std::to_string(ds.minimum.cuts) +
                ", cuts_max=" + std::to_string(ds.maximum.cuts));
        }
        double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start_wall).count();
        if (selected_cplex()) {
#ifdef DECIDEIT_WITH_CPLEX
            const CplexRunStats &cs = last_cplex_stats();
            printf("RESULT {\"level\":0,\"solver\":\"cplex-global\",\"cplex_version\":\"%s\",\"threads\":%d,\"time_limit_s\":%.9g,\"min\":%.17g,\"max\":%.17g,\"wall_s\":%.9g,\"build_s\":%.9g,\"solve_s\":%.9g,\"min_bound\":%.17g,\"max_bound\":%.17g,\"min_nodes\":%lld,\"max_nodes\":%lld,\"min_status\":%d,\"max_status\":%d,\"status\":\"optimal\"}\n",
                   cs.version.c_str(),cplex_threads(),cplex_time_limit(),
                   P_min_point[0],-P_max_point[0],wall,cs.build_s,cs.min_s+cs.max_s,
                   cs.min_bound,cs.max_bound,cs.min_nodes,cs.max_nodes,cs.min_status,cs.max_status);
#endif
        } else {
            const DecideItRunStats &ds = last_decideit_stats();
            printf("RESULT {\"level\":0,\"solver\":\"decideit-highs\",\"method\":\"%s\",\"climb\":\"%s\",\"parallel_directions\":%s,\"rlt_start\":\"%s\",\"trigger\":%d,\"min\":%.17g,\"max\":%.17g,\"wall_s\":%.9g,\"cuts_min\":%d,\"cuts_max\":%d,\"amc_rounds\":%d,\"block_lp_solves\":%d,\"neighbor_lp_solves\":%d,\"neighborhood_calls\":%d,\"climb_restarts\":%d,\"rlt_auxiliary_solves\":%d,\"rlt_skipped_solves\":%d,\"rlt_auxiliary_s\":%.9g,\"lb_calls_min\":%d,\"lb_calls_max\":%d,\"lb_s_min\":%.9g,\"lb_s_max\":%.9g,\"cold_solves\":%d,\"warm_solves\":%d,\"cross_level_basis_uses\":%d,\"cross_level_primal_uses\":%d,\"cross_level_primal_scale_sum\":%.9g,\"simplex_iterations\":%lld,\"rlt_max_rows\":%d,\"rlt_columns\":%d,\"status\":\"optimal\"}\n",
                   selected_method().c_str(), selected_climb().c_str(),
                   std::getenv("DECIDEIT_PARALLEL_DIRECTIONS") ? "true" : "false",
                   selected_rlt_start().c_str(), selected_trigger(),
                   P_min_point[0],-P_max_point[0],wall,(int)V_min_point[0],(int)V_max_point[0],
                   ds.minimum.amc_rounds + ds.maximum.amc_rounds,
                   ds.minimum.block_lp_solves + ds.maximum.block_lp_solves,
                   ds.minimum.neighbor_lp_solves + ds.maximum.neighbor_lp_solves,
                   ds.minimum.neighborhood_calls + ds.maximum.neighborhood_calls,
                   ds.minimum.climb_restarts + ds.maximum.climb_restarts,
                   ds.minimum.rlt_auxiliary_solves + ds.maximum.rlt_auxiliary_solves,
                   ds.minimum.rlt_skipped_solves + ds.maximum.rlt_skipped_solves,
                   ds.minimum.rlt_auxiliary_s + ds.maximum.rlt_auxiliary_s,
                   ds.minimum.lower_bound_calls, ds.maximum.lower_bound_calls,
                   ds.minimum.lower_bound_s, ds.maximum.lower_bound_s,
                   ds.minimum.cold_solves + ds.maximum.cold_solves,
                   ds.minimum.warm_solves + ds.maximum.warm_solves,
                   ds.minimum.cross_level_basis_uses +
                       ds.maximum.cross_level_basis_uses,
                   ds.minimum.cross_level_primal_uses +
                       ds.maximum.cross_level_primal_uses,
                   ds.minimum.cross_level_primal_scale_sum +
                       ds.maximum.cross_level_primal_scale_sum,
                   ds.minimum.simplex_iterations + ds.maximum.simplex_iterations,
                   ds.minimum.max_rows > ds.maximum.max_rows
                       ? ds.minimum.max_rows : ds.maximum.max_rows,
                   ds.minimum.columns > ds.maximum.columns
                       ? ds.minimum.columns : ds.maximum.columns);
            print_lp_backend_stats(0);
            print_rlt_backend_stats(0);
            print_row_reduction_stats(0);
        }
        return 0;
    }

    // start = clock(); // the cutting plane method starts

    // exit_flag = cpm_solver(
    //     P_lobo,
    //     P_upbo,
    //     P_links,
    //     V_lobo,
    //     V_upbo,
    //     V_links,
    //     P_min_point,
    //     V_min_point,
    //     P_max_point,
    //     V_max_point,
    //     Ai_begin,
    //     Ai_end,
    //     Aj_begin,
    //     Aj_end,
    //     CUTS_MAX,
    //     LC_V,
    //     LC_P,
    //     &mX,
    //     &mY);


    // end = clock(); // the cutting plane method ends

    // cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;

    // printf("\nexitFlag: %d,   min: % 7.6f,   max: % 7.6f,   time: %9.3f ms,   cuts_min: % 3d,   cuts_max: % 3d,   P_lobo[0]: % 2d,   P_upbo[0]: % 2d\n", exit_flag, P_min_point[0], -P_max_point[0], cpu_time_used, (int)V_min_point[0], (int)V_max_point[0], (int)P_lobo[0], (int)P_upbo[0]);

    // fflush(stdout);
    
    const std::vector<LinCons> original_V_rows = V_rows;
    const std::vector<LinCons> original_P_rows = P_rows;
    duplicate_rawData(instance_copy, instance);
    if(!fixed_statements)
        adjust_constraints(instance_copy, num_alt_i, num_alt_j, focal_V, focal_P);

    // Keep the solver rows aligned with any right-hand-side relaxation made
    // to include the focal point.  This works for both legacy differences and
    // arbitrary sparse linear constraints.
    for (size_t row = 0; row < V_rows.size(); ++row)
        V_rows[row].b = gsl_vector_get(instance_copy->bx, row);
    for (size_t row = 0; row < P_rows.size(); ++row)
        P_rows[row].b = gsl_vector_get(instance_copy->by, row);

    gsl_vector *mass_V = managed::vector_calloc(focal_V[0]);
    gsl_vector *mass_P = managed::vector_calloc(focal_P[0]);

    for (i = 0; i < focal_V[0]; i++)
    {
        gsl_vector_set(mass_V, i, focal_V[i + 1]);
    }

    for (i = 0; i < focal_P[0]; i++)
    {
        gsl_vector_set(mass_P, i, focal_P[i + 1]);
    }

    gsl_vector *increment_V = managed::vector_calloc(focal_V[0]);
    gsl_vector *increment_V_right = managed::vector_calloc(focal_V[0]);
    gsl_vector *increment_P_left = managed::vector_calloc(focal_P[0]);
    gsl_vector *increment_P_right = managed::vector_calloc(focal_P[0]);
    gsl_vector *base_xl = managed::vector_alloc(instance_copy->xl->size);
    gsl_vector *base_xu = managed::vector_alloc(instance_copy->xu->size);
    gsl_vector *base_yl = managed::vector_alloc(instance_copy->yl->size);
    gsl_vector *base_yu = managed::vector_alloc(instance_copy->yu->size);
    gsl_vector_memcpy(base_xl, instance_copy->xl);
    gsl_vector_memcpy(base_xu, instance_copy->xu);
    gsl_vector_memcpy(base_yl, instance_copy->yl);
    gsl_vector_memcpy(base_yu, instance_copy->yu);

    const int counter_contraction = atoi(argv[2]);
    const double delta = 1.0 / counter_contraction;

    /* derive the increments for V-/P-base, respectively */
    gsl_vector_memcpy(increment_V, mass_V);
    gsl_vector_sub(increment_V, instance_copy->xl);
    gsl_vector_scale(increment_V, delta);
    gsl_vector_memcpy(increment_V_right, mass_V);
    gsl_vector_sub(increment_V_right, instance_copy->xu);
    gsl_vector_scale(increment_V_right, delta);

    gsl_vector_memcpy(increment_P_left, mass_P);
    gsl_vector_sub(increment_P_left, instance_copy->yl);
    gsl_vector_scale(increment_P_left, delta);

    gsl_vector_memcpy(increment_P_right, mass_P);
    gsl_vector_sub(increment_P_right, instance_copy->yu);
    gsl_vector_scale(increment_P_right, -delta);

    gsl_matrix *final_results = managed::matrix_calloc(4, counter_contraction + 1);
    gsl_matrix *final_results_heuristic = managed::matrix_calloc(4, counter_contraction + 1);
    gsl_vector *whole_solution = managed::vector_calloc(focal_V[0]);
    gsl_vector_view src_v1, src_v2;
    d_row previous_P_min = {}, previous_V_min = {};
    d_row previous_P_max = {}, previous_V_max = {};
    d_row penultimate_P_min = {}, penultimate_P_max = {};
    d_row default_start = {}, midpoint_start = {};
    d_row transport_min_start = {}, transport_max_start = {};
    d_row transport_midpoint_start = {};
    d_row ray_min_start = {}, ray_max_start = {};
    d_row secant_min_candidate = {}, secant_max_candidate = {};
    d_row secant_min_start = {}, secant_max_start = {};
    for (j = 1; j <= focal_P[0]; ++j) default_start[j] = 1.0;

    // LinCons LC_V[V_links[0].head + 1], LC_P[P_links[0].head + 1];

    

    const std::string contraction_path = selected_contraction_path();
    const int first_level = contraction_path == "inclusive" ? 0 : 1;
    const std::string contraction_order = selected_contraction_order();
    const std::string contraction_reduction =
        selected_contraction_reduction() == "reuse" &&
        contraction_order == "inward" ? "reuse" : "refresh";
    const std::string contraction_precheck = selected_contraction_precheck();
    const std::string base_method = selected_method();
    const int bcm_through = selected_bcm_through();
    std::vector<int> contraction_levels;
    if (contraction_order == "inward")
        for (int level = first_level; level < counter_contraction; ++level)
            contraction_levels.push_back(level);
    else
        for (int level = counter_contraction - 1; level >= first_level; --level)
            contraction_levels.push_back(level);
    bool previous_solution_available = false;
    bool penultimate_solution_available = false;
    int previous_level = -1;
    int penultimate_level = -1;
    RltCrossLevelState min_rlt_state, max_rlt_state;
    const bool cross_level_rlt_basis =
        std::getenv("DECIDEIT_CROSS_LEVEL_RLT_START") != nullptr;
    if (contraction_reduction == "reuse")
        setenv("DECIDEIT_REUSE_REDUCTION", "1", 1);
    else
        unsetenv("DECIDEIT_REUSE_REDUCTION");
    /* Solve only non-point levels.  The focal endpoint is evaluated directly
       below, so the default inclusive path with K contractions solves K
       levels, or exactly 2*K DBLP directions.  The optional interior path
       omits level zero when it was solved separately. */
    for (size_t level_index = 0; level_index < contraction_levels.size(); ++level_index)
    {
        i = contraction_levels[level_index];
        const bool original_problem = i == 0;
        if (!selected_cplex())
            setenv("DECIDEIT_METHOD",
                   i <= bcm_through ? "bcm" : base_method.c_str(), 1);
        const std::vector<LinCons> &level_P_rows =
            original_problem ? original_P_rows : P_rows;
        const std::vector<LinCons> &level_V_rows =
            original_problem ? original_V_rows : V_rows;
        struct raw_data *level_instance = original_problem ? instance : instance_copy;
        const bool skip_precheck = contraction_precheck == "skip" ||
            (contraction_precheck == "auto" && !original_problem);
        if (skip_precheck)
            setenv("DECIDEIT_SKIP_INITIAL_X_SOLVE", "1", 1);
        else
            unsetenv("DECIDEIT_SKIP_INITIAL_X_SOLVE");

        /* Level zero uses the unrelaxed input rows.  Reduction state from
           that model cannot be carried into level one, whose rows may have
           been relaxed to include the focal point. */
        const bool crossed_original_boundary =
            level_index > 0 && (i == 0 || previous_level == 0);
        if (level_index == 0 || contraction_order == "outward" ||
            crossed_original_boundary)
        {
            mX = -1;
            mY = -1;
        }

        /* Set this level directly so the path may be traversed in either
           direction without accumulating floating-point or state errors. */
        gsl_vector_memcpy(instance_copy->xl, base_xl);
        gsl_vector_memcpy(instance_copy->xu, base_xu);
        gsl_vector_memcpy(instance_copy->yl, base_yl);
        gsl_vector_memcpy(instance_copy->yu, base_yu);
        gsl_blas_daxpy(static_cast<double>(i), increment_V, instance_copy->xl);
        // Legacy is retained byte-for-byte in its historical midpoint path.
        // A projected (non-midpoint) focus needs its own upper increment.
        if(fixed_statements)
            gsl_blas_daxpy(static_cast<double>(i),increment_V_right,instance_copy->xu);
        else
            gsl_blas_daxpy(-static_cast<double>(i),increment_V,instance_copy->xu);
        gsl_blas_daxpy(static_cast<double>(i), increment_P_left, instance_copy->yl);
        gsl_blas_daxpy(-static_cast<double>(i), increment_P_right, instance_copy->yu);

        for (j = 1; j <= instance_copy->xl->size; j++)
        {
            V_lobo[j] = gsl_vector_get(instance_copy->xl, j - 1);
            V_upbo[j] = gsl_vector_get(instance_copy->xu, j - 1);
        }

        for (j = 1; j <= instance_copy->yl->size; j++)
        {
            P_lobo[j] = gsl_vector_get(instance_copy->yl, j - 1);
            P_upbo[j] = gsl_vector_get(instance_copy->yu, j - 1);
        }

        clean_results(P_min_point, V_min_point, P_max_point, V_max_point, focal_V[0]);
        dp::stage("dblp_solve",i);
        dp::log("\"event\":\"level_model\",\"level\":"+std::to_string(i)+
            ",\"fraction\":"+dp::number(static_cast<double>(i)/counter_contraction)+
            ",\"rows_sha256\":"+dp::quote(dp::rows_hash(level_instance))+
            ",\"bounds_sha256\":"+dp::quote(dp::bounds_hash(level_instance))+
            ",\"focal_sha256\":"+dp::quote(focal_info.hash)+
            ",\"p_lower\":"+dp::array(P_lobo,Aj_end)+",\"p_upper\":"+dp::array(P_upbo,Aj_end)+
            ",\"v_lower\":"+dp::array(V_lobo,Aj_end)+",\"v_upper\":"+dp::array(V_upbo,Aj_end));

        // if (i == 7)
        // {
        //     printf("\n");
        // }

        // output_datafile(instance_copy, delta * i * 100, file_name, num_alt_i, num_alt_j);

        auto start_wall=std::chrono::steady_clock::now();
        start = clock();

        const std::string contraction_start = selected_contraction_start();
        const double *min_start = default_start, *max_start = default_start;
        double transport_min_scale = 0.0, transport_max_scale = 0.0;
        double transport_midpoint_scale = 0.0;
        double ray_min_scale = 0.0, ray_max_scale = 0.0;
        double secant_min_scale = 0.0, secant_max_scale = 0.0;
        if (previous_solution_available)
        {
            const double transport_ratio = (counter_contraction - i) /
                static_cast<double>(counter_contraction - previous_level);
            for (j = 1; j <= focal_P[0]; ++j)
                midpoint_start[j] =
                    0.5 * (previous_P_min[j] + previous_P_max[j]);

            auto feasible_ray_scale = [&](const double *source)
            {
                double scale = transport_ratio;
                if (scale <= 1.0) return scale;
                for (int k = Ai_begin; k <= Aj_end; ++k)
                {
                    if (k > Ai_end && k < Aj_begin) continue;
                    const double direction = source[k] - focal_P[k];
                    if (direction > 1e-14)
                        scale = std::fmin(scale,
                            (P_upbo[k] - focal_P[k]) / direction);
                    else if (direction < -1e-14)
                        scale = std::fmin(scale,
                            (P_lobo[k] - focal_P[k]) / direction);
                }
                for (const LinCons &row : level_P_rows)
                {
                    double focal_lhs = 0.0, direction_lhs = 0.0;
                    for (size_t z = 0; z < row.index.size(); ++z)
                    {
                        const int k = row.index[z];
                        focal_lhs += row.coefficient[z] * focal_P[k];
                        direction_lhs += row.coefficient[z] *
                            (source[k] - focal_P[k]);
                    }
                    if (direction_lhs > 1e-14)
                        scale = std::fmin(scale,
                            (row.b - focal_lhs) / direction_lhs);
                }
                return std::fmax(1.0, std::fmin(transport_ratio, scale));
            };
            transport_min_scale = feasible_ray_scale(previous_P_min);
            transport_max_scale = feasible_ray_scale(previous_P_max);
            transport_midpoint_scale = feasible_ray_scale(midpoint_start);

            /* A transported point is guaranteed feasible but can be much
               closer to the focal point than necessary.  Clip the segment
               from the focal point to a normalized candidate at the last
               bound or general-row intersection instead. */
            auto maximal_feasible_scale = [&](const double *candidate)
            {
                double scale = 1.0;
                for (int k = Ai_begin; k <= Aj_end; ++k)
                {
                    if (k > Ai_end && k < Aj_begin) continue;
                    const double direction = candidate[k] - focal_P[k];
                    if (direction > 1e-14)
                        scale = std::fmin(scale,
                            (P_upbo[k] - focal_P[k]) / direction);
                    else if (direction < -1e-14)
                        scale = std::fmin(scale,
                            (P_lobo[k] - focal_P[k]) / direction);
                }
                for (const LinCons &row : level_P_rows)
                {
                    double focal_lhs = 0.0, direction_lhs = 0.0;
                    for (size_t z = 0; z < row.index.size(); ++z)
                    {
                        const int k = row.index[z];
                        focal_lhs += row.coefficient[z] * focal_P[k];
                        direction_lhs += row.coefficient[z] *
                            (candidate[k] - focal_P[k]);
                    }
                    if (direction_lhs > 1e-14)
                        scale = std::fmin(scale,
                            (row.b - focal_lhs) / direction_lhs);
                }
                return std::fmax(0.0, std::fmin(1.0, scale));
            };

            ray_min_scale = maximal_feasible_scale(previous_P_min);
            ray_max_scale = maximal_feasible_scale(previous_P_max);
            for (j = 1; j <= focal_P[0]; ++j)
            {
                ray_min_start[j] = focal_P[j] + ray_min_scale *
                    (previous_P_min[j] - focal_P[j]);
                ray_max_start[j] = focal_P[j] + ray_max_scale *
                    (previous_P_max[j] - focal_P[j]);
            }

            if (penultimate_solution_available)
            {
                const double step = (i - previous_level) /
                    static_cast<double>(previous_level - penultimate_level);
                for (j = 1; j <= focal_P[0]; ++j)
                {
                    secant_min_candidate[j] = previous_P_min[j] + step *
                        (previous_P_min[j] - penultimate_P_min[j]);
                    secant_max_candidate[j] = previous_P_max[j] + step *
                        (previous_P_max[j] - penultimate_P_max[j]);
                }
            }
            else
            {
                std::memcpy(secant_min_candidate, previous_P_min, sizeof(d_row));
                std::memcpy(secant_max_candidate, previous_P_max, sizeof(d_row));
            }
            secant_min_scale = maximal_feasible_scale(secant_min_candidate);
            secant_max_scale = maximal_feasible_scale(secant_max_candidate);
            for (j = 1; j <= focal_P[0]; ++j)
            {
                transport_min_start[j] = focal_P[j] + transport_min_scale *
                    (previous_P_min[j] - focal_P[j]);
                transport_max_start[j] = focal_P[j] + transport_max_scale *
                    (previous_P_max[j] - focal_P[j]);
                transport_midpoint_start[j] = focal_P[j] +
                    transport_midpoint_scale *
                    (midpoint_start[j] - focal_P[j]);
                secant_min_start[j] = focal_P[j] + secant_min_scale *
                    (secant_min_candidate[j] - focal_P[j]);
                secant_max_start[j] = focal_P[j] + secant_max_scale *
                    (secant_max_candidate[j] - focal_P[j]);
            }
        }
        if (contraction_start == "previous" && previous_solution_available)
        {
            min_start = previous_P_min;
            max_start = previous_P_max;
        }
        else if (contraction_start == "focal")
        {
            min_start = focal_P;
            max_start = focal_P;
        }
        else if (contraction_start == "midpoint" && previous_solution_available)
        {
            min_start = midpoint_start;
            max_start = midpoint_start;
        }
        else if (contraction_start == "transport" && previous_solution_available)
        {
            min_start = transport_min_start;
            max_start = transport_max_start;
        }
        else if (contraction_start == "transport-midpoint" && previous_solution_available)
        {
            min_start = transport_midpoint_start;
            max_start = transport_midpoint_start;
        }
        else if (contraction_start == "ray" && previous_solution_available)
        {
            min_start = ray_min_start;
            max_start = ray_max_start;
        }
        else if (contraction_start == "secant" && previous_solution_available)
        {
            min_start = secant_min_start;
            max_start = secant_max_start;
        }
        else if (contraction_start == "transport" ||
                 contraction_start == "transport-midpoint" ||
                 contraction_start == "ray" || contraction_start == "secant")
        {
            // There is no preceding contracted solution at the first level.
            // The focal point is always feasible for arbitrary linear bases.
            min_start = focal_P;
            max_start = focal_P;
        }
        /* The focal point need not satisfy the user's original rows.  A null
           seed reproduces standalone --single behavior at 0% contraction. */
        if (original_problem)
        {
            min_start = nullptr;
            max_start = nullptr;
        }

        dp::log("\"event\":\"start_reuse\",\"level\":"+std::to_string(i)+
            ",\"strategy\":"+dp::quote(contraction_start)+
            ",\"previous_available\":"+(previous_solution_available?"true":"false")+
            ",\"previous_min_full_violation\":"+dp::number(previous_solution_available?
                dp::point_violation(level_instance,previous_P_min,previous_V_min,num_alt_i):std::numeric_limits<double>::quiet_NaN())+
            ",\"previous_max_full_violation\":"+dp::number(previous_solution_available?
                dp::point_violation(level_instance,previous_P_max,previous_V_max,num_alt_i):std::numeric_limits<double>::quiet_NaN())+
            ",\"min_seed_p\":"+(min_start?dp::array(min_start,Aj_end):"null")+
            ",\"max_seed_p\":"+(max_start?dp::array(max_start,Aj_end):"null")+
            ",\"reuse_scope\":"+dp::quote(selected_cplex()
                ? "candidate only; CPLEX baseline does not consume the AMC seed"
                : "P-block AMC seed only; no inherited optimality certificate"));
        exit_flag = solve_selected(
            P_lobo,
            P_upbo,
            P_links,
            level_P_rows,
            V_lobo,
            V_upbo,
            V_links,
            level_V_rows,
            P_min_point,
            V_min_point,
            P_max_point,
            V_max_point,
            Ai_begin,
            Ai_end,
            Aj_begin,
            Aj_end,
            CUTS_MAX,
            LC_V,
            LC_P,
            &mX,
            &mY,
            min_start,
            max_start,
            cross_level_rlt_basis ? &min_rlt_state : nullptr,
            cross_level_rlt_basis ? &max_rlt_state : nullptr);

        end = clock();

        cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;
        emit_directions(i,level_instance,exit_flag);
        if(exit_flag==OK && (dp::point_violation(level_instance,P_min_point,V_min_point,num_alt_i)>1e-8 ||
                            dp::point_violation(level_instance,P_max_point,V_max_point,num_alt_i)>1e-8))
            throw std::runtime_error("Returned contraction witness is infeasible");

        // printf("exitFlag: %d,   min: % 7.6f,   max: % 7.6f,   time: %9.3f ms,   cuts_min: % 4d,   cuts_max: % 4d,   P_lobo[0]: % 2d,   P_upbo[0]: % 2d\n", exit_flag, P_min_point[0], -P_max_point[0], cpu_time_used, (int)V_min_point[0], (int)V_max_point[0], (int)P_lobo[0], (int)P_upbo[0]);

        if ((exit_flag != OK))
        {
#ifdef DECIDEIT_WITH_CPLEX
            if (selected_cplex()) {
                const CplexRunStats &cs = last_cplex_stats();
                throw std::runtime_error("CPLEX contraction " + std::to_string(i) +
                                         " not certified: min=" + cs.min_status_text +
                                         ", max=" + cs.max_status_text);
            }
#endif
            const DecideItRunStats& ds=last_decideit_stats();
            throw std::runtime_error("Contraction "+std::to_string(i)+" not certified, status="+
                std::to_string(exit_flag)+", min_termination="+std::to_string(ds.minimum.termination_flag)+
                ", max_termination="+std::to_string(ds.maximum.termination_flag));
        }

        // printf("\n< cpm > ");
        // display_solution(P_min_point, V_min_point, P_max_point, V_max_point, Ai_begin, Ai_end, Aj_begin, Aj_end);

        src_v1 = gsl_vector_view_array(V_min_point, focal_V[0] + 1);
        src_v2 = gsl_vector_subvector(&src_v1.vector, 1, focal_V[0]);
        gsl_vector_memcpy(whole_solution, &src_v2.vector);
        final_consistency_check(level_instance, whole_solution);

        auto point_objective = [&](const double *p, const double *v)
        {
            double value = 0.0;
            for (int k = Ai_begin; k <= Ai_end; ++k) value += p[k] * v[k];
            for (int k = Aj_begin; k <= Aj_end; ++k) value -= p[k] * v[k];
            return value;
        };
        auto point_in_contracted_box = [&](const double *p, const double *v)
        {
            double sum_i = 0.0, sum_j = 0.0;
            for (int k = Ai_begin; k <= Aj_end; ++k)
            {
                if (k > Ai_end && k < Aj_begin) continue;
                if (!std::isfinite(p[k]) || !std::isfinite(v[k]) ||
                    p[k] < P_lobo[k] - 1e-8 || p[k] > P_upbo[k] + 1e-8 ||
                    v[k] < V_lobo[k] - 1e-8 || v[k] > V_upbo[k] + 1e-8)
                    return false;
                (k <= Ai_end ? sum_i : sum_j) += p[k];
            }
            if (std::fabs(sum_i - 1.0) > 1e-8 ||
                std::fabs(sum_j - 1.0) > 1e-8) return false;
            auto rows_hold = [&](const std::vector<LinCons> &rows,
                                 const double *point)
            {
                for (const LinCons &row : rows)
                {
                    double lhs = 0.0;
                    for (size_t z = 0; z < row.index.size(); ++z)
                        lhs += row.coefficient[z] * point[row.index[z]];
                    if (lhs > row.b + 1e-8) return false;
                }
                return true;
            };
            if (!rows_hold(level_P_rows, p) || !rows_hold(level_V_rows, v))
                return false;
            return true;
        };
        auto point_distance = [&](const double *p1, const double *v1,
                                  const double *p2, const double *v2)
        {
            double result = 0.0;
            for (int k = Ai_begin; k <= Aj_end; ++k)
            {
                if (k > Ai_end && k < Aj_begin) continue;
                result = std::fmax(result, std::fabs(p1[k] - p2[k]));
                result = std::fmax(result, std::fabs(v1[k] - v2[k]));
            }
            return result;
        };
        const bool previous_available = previous_solution_available;
        const bool previous_min_feasible = previous_available &&
            point_in_contracted_box(previous_P_min, previous_V_min);
        const bool previous_max_feasible = previous_available &&
            point_in_contracted_box(previous_P_max, previous_V_max);
        const bool previous_min_still_optimal = previous_min_feasible &&
            std::fabs(point_objective(previous_P_min, previous_V_min) -
                      P_min_point[0]) <= DBLP_EPSILON;
        const bool previous_max_still_optimal = previous_max_feasible &&
            std::fabs(point_objective(previous_P_max, previous_V_max) +
                      P_max_point[0]) <= DBLP_EPSILON;
        const double previous_min_linf = previous_available
            ? point_distance(previous_P_min, previous_V_min,
                             P_min_point, V_min_point) : 0.0;
        const double previous_max_linf = previous_available
            ? point_distance(previous_P_max, previous_V_max,
                             P_max_point, V_max_point) : 0.0;
        const double focal_min_linf = point_distance(
            focal_P, focal_V, P_min_point, V_min_point);
        const double focal_max_linf = point_distance(
            focal_P, focal_V, P_max_point, V_max_point);
        auto probability_seed_feasible = [&](const double *p)
        {
            if (!p) return false;
            double sum_i = 0.0, sum_j = 0.0;
            for (int k = Ai_begin; k <= Aj_end; ++k)
            {
                if (k > Ai_end && k < Aj_begin) continue;
                if (!std::isfinite(p[k]) || p[k] < P_lobo[k] - 1e-8 || p[k] > P_upbo[k] + 1e-8)
                    return false;
                (k <= Ai_end ? sum_i : sum_j) += p[k];
            }
            if (std::fabs(sum_i - 1.0) > 1e-8 ||
                std::fabs(sum_j - 1.0) > 1e-8) return false;
            for (const LinCons &row : level_P_rows)
            {
                double lhs = 0.0;
                for (size_t z = 0; z < row.index.size(); ++z)
                    lhs += row.coefficient[z] * p[row.index[z]];
                if (lhs > row.b + 1e-8) return false;
            }
            return true;
        };
        auto probability_distance = [&](const double *p1, const double *p2)
        {
            if (!p1 || !p2) return 0.0;
            double result = 0.0;
            for (int k = Ai_begin; k <= Aj_end; ++k)
            {
                if (k > Ai_end && k < Aj_begin) continue;
                result = std::fmax(result, std::fabs(p1[k] - p2[k]));
            }
            return result;
        };
        const bool start_min_feasible = probability_seed_feasible(min_start);
        const bool start_max_feasible = probability_seed_feasible(max_start);
        const double start_min_linf = probability_distance(min_start, P_min_point);
        const double start_max_linf = probability_distance(max_start, P_max_point);

        src_v1 = gsl_vector_view_array(V_max_point, focal_V[0] + 1);
        src_v2 = gsl_vector_subvector(&src_v1.vector, 1, focal_V[0]);
        gsl_vector_memcpy(whole_solution, &src_v2.vector);
        final_consistency_check(level_instance, whole_solution);

        double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start_wall).count();
        if (selected_cplex()) {
#ifdef DECIDEIT_WITH_CPLEX
            const CplexRunStats &cs = last_cplex_stats();
            printf("RESULT {\"level\":%d,\"original_problem\":%s,\"solver\":\"cplex-global\",\"contraction_path\":\"%s\",\"contraction_order\":\"%s\",\"cplex_version\":\"%s\",\"threads\":%d,\"time_limit_s\":%.9g,\"min\":%.17g,\"max\":%.17g,\"cpu_ms\":%.9g,\"wall_s\":%.9g,\"build_s\":%.9g,\"solve_s\":%.9g,\"min_bound\":%.17g,\"max_bound\":%.17g,\"min_nodes\":%lld,\"max_nodes\":%lld,\"min_status\":%d,\"max_status\":%d,\"status\":\"optimal\"}\n",
                   i,original_problem ? "true" : "false",
                   contraction_path.c_str(),contraction_order.c_str(),
                   cs.version.c_str(),cplex_threads(),cplex_time_limit(),
                   P_min_point[0],-P_max_point[0],cpu_time_used,wall,cs.build_s,
                   cs.min_s+cs.max_s,cs.min_bound,cs.max_bound,cs.min_nodes,cs.max_nodes,
                   cs.min_status,cs.max_status);
#endif
        } else {
            const DecideItRunStats &ds = last_decideit_stats();
            printf("RESULT {\"level\":%d,\"original_problem\":%s,\"solver\":\"decideit-highs\",\"method\":\"%s\",\"climb\":\"%s\",\"parallel_directions\":%s,\"contraction_path\":\"%s\",\"contraction_order\":\"%s\",\"contraction_reduction\":\"%s\",\"contraction_precheck\":\"%s\",\"contraction_start\":\"%s\",\"rlt_start\":\"%s\",\"trigger\":%d,\"min\":%.17g,\"max\":%.17g,\"cpu_ms\":%.9g,\"wall_s\":%.9g,\"cuts_min\":%d,\"cuts_max\":%d,\"amc_rounds\":%d,\"block_lp_solves\":%d,\"neighbor_lp_solves\":%d,\"neighborhood_calls\":%d,\"climb_restarts\":%d,\"rlt_auxiliary_solves\":%d,\"rlt_skipped_solves\":%d,\"rlt_auxiliary_s\":%.9g,\"lb_calls_min\":%d,\"lb_calls_max\":%d,\"lb_s_min\":%.9g,\"lb_s_max\":%.9g,\"cold_solves\":%d,\"warm_solves\":%d,\"cross_level_basis_uses\":%d,\"cross_level_primal_uses\":%d,\"cross_level_primal_scale_sum\":%.9g,\"simplex_iterations\":%lld,\"rlt_max_rows\":%d,\"rlt_columns\":%d,\"previous_available\":%s,\"previous_min_feasible\":%s,\"previous_max_feasible\":%s,\"previous_min_still_optimal\":%s,\"previous_max_still_optimal\":%s,\"previous_min_linf\":%.9g,\"previous_max_linf\":%.9g,\"focal_min_linf\":%.9g,\"focal_max_linf\":%.9g,\"start_min_feasible\":%s,\"start_max_feasible\":%s,\"start_min_linf\":%.9g,\"start_max_linf\":%.9g,\"transport_min_scale\":%.9g,\"transport_max_scale\":%.9g,\"transport_midpoint_scale\":%.9g,\"status\":\"optimal\"}\n",
                   i,original_problem ? "true" : "false",
                   selected_method().c_str(),selected_climb().c_str(),
                   std::getenv("DECIDEIT_PARALLEL_DIRECTIONS") ? "true" : "false",
                   contraction_path.c_str(),contraction_order.c_str(),
                   contraction_reduction.c_str(),
                   contraction_precheck.c_str(),
                   contraction_start.c_str(),selected_rlt_start().c_str(),
                   selected_trigger(),P_min_point[0],-P_max_point[0],cpu_time_used,wall,
                   (int)V_min_point[0],(int)V_max_point[0],
                   ds.minimum.amc_rounds + ds.maximum.amc_rounds,
                   ds.minimum.block_lp_solves + ds.maximum.block_lp_solves,
                   ds.minimum.neighbor_lp_solves + ds.maximum.neighbor_lp_solves,
                   ds.minimum.neighborhood_calls + ds.maximum.neighborhood_calls,
                   ds.minimum.climb_restarts + ds.maximum.climb_restarts,
                   ds.minimum.rlt_auxiliary_solves + ds.maximum.rlt_auxiliary_solves,
                   ds.minimum.rlt_skipped_solves + ds.maximum.rlt_skipped_solves,
                   ds.minimum.rlt_auxiliary_s + ds.maximum.rlt_auxiliary_s,
                   ds.minimum.lower_bound_calls, ds.maximum.lower_bound_calls,
                   ds.minimum.lower_bound_s, ds.maximum.lower_bound_s,
                   ds.minimum.cold_solves + ds.maximum.cold_solves,
                   ds.minimum.warm_solves + ds.maximum.warm_solves,
                   ds.minimum.cross_level_basis_uses +
                       ds.maximum.cross_level_basis_uses,
                   ds.minimum.cross_level_primal_uses +
                       ds.maximum.cross_level_primal_uses,
                   ds.minimum.cross_level_primal_scale_sum +
                       ds.maximum.cross_level_primal_scale_sum,
                   ds.minimum.simplex_iterations + ds.maximum.simplex_iterations,
                   ds.minimum.max_rows > ds.maximum.max_rows
                       ? ds.minimum.max_rows : ds.maximum.max_rows,
                   ds.minimum.columns > ds.maximum.columns
                       ? ds.minimum.columns : ds.maximum.columns,
                   previous_available ? "true" : "false",
                   previous_min_feasible ? "true" : "false",
                   previous_max_feasible ? "true" : "false",
                   previous_min_still_optimal ? "true" : "false",
                   previous_max_still_optimal ? "true" : "false",
                   previous_min_linf, previous_max_linf,
                   focal_min_linf, focal_max_linf,
                   start_min_feasible ? "true" : "false",
                   start_max_feasible ? "true" : "false",
                   start_min_linf, start_max_linf,
                   transport_min_scale, transport_max_scale,
                   transport_midpoint_scale);
            print_lp_backend_stats(i);
            print_rlt_backend_stats(i);
            print_row_reduction_stats(i);
        }
        gsl_matrix_set(final_results, 0, i, P_min_point[0]);
        gsl_matrix_set(final_results, 1, i, (P_min_point[0] - P_max_point[0]) / 2);
        gsl_matrix_set(final_results, 2, i, -P_max_point[0]);
        gsl_matrix_set(final_results, 3, i, cpu_time_used);

        if (previous_solution_available)
        {
            std::memcpy(penultimate_P_min, previous_P_min, sizeof(d_row));
            std::memcpy(penultimate_P_max, previous_P_max, sizeof(d_row));
            penultimate_solution_available = true;
            penultimate_level = previous_level;
        }
        std::memcpy(previous_P_min, P_min_point, sizeof(d_row));
        std::memcpy(previous_V_min, V_min_point, sizeof(d_row));
        std::memcpy(previous_P_max, P_max_point, sizeof(d_row));
        std::memcpy(previous_V_max, V_max_point, sizeof(d_row));
        previous_solution_available = true;
        previous_level = i;

        // clean_results(P_min_point, V_min_point, P_max_point, V_max_point, focal_V[0]);

        // start = clock();

        // exit_flag = heuristic_solver(P_lobo, P_upbo, P_links, V_lobo, V_upbo, V_links, P_min_point, V_min_point, P_max_point, V_max_point, Ai_begin, Ai_end, Aj_begin, Aj_end, CUTS_MAX);

        // end = clock();

        // cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;

        // printf("exitFlag: %d,   min: % 7.6f,   max: % 7.6f,   time: %9.3f ms,   cuts_min: % 4d,   cuts_max: % 4d,   P_lobo[0]: % 2d,   P_upbo[0]: % 2d\n", exit_flag, P_min_point[0], -P_max_point[0], cpu_time_used, (int)V_min_point[0], (int)V_max_point[0], (int)P_lobo[0], (int)P_upbo[0]);

        // if ((exit_flag != OK))
        // {
        //     printf("\n\n\nthe %dth contraction failed\n\n\n", i);
        //     return (FALSE);
        // }

        // printf("\n< heuristic > ");
        // display_solution(P_min_point, V_min_point, P_max_point, V_max_point, Ai_begin, Ai_end, Aj_begin, Aj_end);
        // printf("*******\n");

        // src_v1 = gsl_vector_view_array(V_min_point, focal_V[0] + 1);
        // src_v2 = gsl_vector_subvector(&src_v1.vector, 1, focal_V[0]);
        // gsl_vector_memcpy(whole_solution, &src_v2.vector);
        // final_consistency_check(instance_copy, whole_solution);

        // src_v1 = gsl_vector_view_array(V_max_point, focal_V[0] + 1);
        // src_v2 = gsl_vector_subvector(&src_v1.vector, 1, focal_V[0]);
        // gsl_vector_memcpy(whole_solution, &src_v2.vector);
        // final_consistency_check(instance_copy, whole_solution);

        // gsl_matrix_set(final_results_heuristic, 0, i, P_min_point[0]);
        // gsl_matrix_set(final_results_heuristic, 1, i, (P_min_point[0] - P_max_point[0]) / 2);
        // gsl_matrix_set(final_results_heuristic, 2, i, -P_max_point[0]);
        // gsl_matrix_set(final_results_heuristic, 3, i, cpu_time_used);

        // results_BARON(focal_P[0], P_max_point, V_max_point);
        // src_v1 = gsl_vector_view_array(V_max_point, focal_V[0] + 1);
        // src_v2 = gsl_vector_subvector(&src_v1.vector, 1, focal_V[0]);
        // gsl_vector_memcpy(whole_solution, &src_v2.vector);
        // final_consistency_check(instance_copy, whole_solution);

        // REAL obj = 0;

        // for (j = Ai_begin; j <= Ai_end; j++)
        // {
        //     obj = obj + P_max_point[j] * V_max_point[j];
        // }

        // for (j = Aj_begin; j <= Aj_end; j++)
        // {
        //     obj = obj - P_max_point[j] * V_max_point[j];
        // }

        // printf("\nthe computed max: %f\n", obj);
    }
    if (!selected_cplex())
        setenv("DECIDEIT_METHOD", base_method.c_str(), 1);

    // At 100% contraction both bases are the focal point.  Direct evaluation
    // is exact and deliberately bypasses the DBLP solver.
    double focal_objective = 0.0;
    dp::stage("direct_point",counter_contraction);
    if(fixed_statements && dp::point_violation(instance,focal_P,focal_V,num_alt_i)>1e-8)
        throw std::runtime_error("Direct focal endpoint is infeasible in original model");
    for (i = Ai_begin; i <= Ai_end; ++i)
        focal_objective += focal_P[i] * focal_V[i];
    for (i = Aj_begin; i <= Aj_end; ++i)
        focal_objective -= focal_P[i] * focal_V[i];
    gsl_matrix_set(final_results, 0, counter_contraction, focal_objective);
    gsl_matrix_set(final_results, 1, counter_contraction, focal_objective);
    gsl_matrix_set(final_results, 2, counter_contraction, focal_objective);
    gsl_matrix_set(final_results, 3, counter_contraction, 0.0);
    printf("RESULT {\"level\":%d,\"solver\":\"direct-point\",\"contraction_path\":\"%s\",\"contraction_order\":\"%s\",\"contraction_reduction\":\"%s\",\"contraction_precheck\":\"%s\",\"contraction_start\":\"%s\",\"min\":%.17g,\"max\":%.17g,\"wall_s\":0,\"direct_point\":true,\"status\":\"optimal\"}\n",
           counter_contraction, contraction_path.c_str(),
           contraction_order.c_str(), contraction_reduction.c_str(),
           contraction_precheck.c_str(),
           selected_contraction_start().c_str(), focal_objective, focal_objective);
    dp::log("\"event\":\"chain_complete\",\"nonpoint_levels\":"+std::to_string(contraction_levels.size())+
        ",\"dblp_queries\":"+std::to_string(2*contraction_levels.size())+
        ",\"direct_point\":true,\"focal_objective\":"+dp::number(focal_objective)+
        ",\"internal_chain_wall_s\":"+dp::number(std::chrono::duration<double>(std::chrono::steady_clock::now()-chain_begin).count())+
        ",\"focal_sha256\":"+dp::quote(focal_info.hash));

    if (debug)
    {
        gsl_vector_view src_v;
        gsl_vector_view shown;

        // printf("filename: %s\n  %%  ", file_name);

        for (i = first_level; i <= counter_contraction; i++)
        {
            printf("% 8.1f%%", delta * i * 100);
        }

        printf("\n< cpm >\n");
        printf(" min: ");
        src_v = gsl_matrix_row(final_results, 0);
        shown = gsl_vector_subvector(&src_v.vector, first_level,
                                     counter_contraction - first_level + 1);
        show_gsl_vector(&shown.vector);
        printf(" mid: ");
        src_v = gsl_matrix_row(final_results, 1);
        shown = gsl_vector_subvector(&src_v.vector, first_level,
                                     counter_contraction - first_level + 1);
        show_gsl_vector(&shown.vector);
        printf(" max: ");
        src_v = gsl_matrix_row(final_results, 2);
        shown = gsl_vector_subvector(&src_v.vector, first_level,
                                     counter_contraction - first_level + 1);
        show_gsl_vector(&shown.vector);
        printf("\ntime: ");
        src_v = gsl_matrix_row(final_results, 3);
        shown = gsl_vector_subvector(&src_v.vector, first_level,
                                     counter_contraction - first_level + 1);
        show_gsl_vector(&shown.vector);
        printf("\n");

        // printf("\n< heuristic >\n");
        // printf(" min: ");
        // src_v = gsl_matrix_row(final_results_heuristic, 0);
        // show_gsl_vector(&src_v.vector);
        // printf(" mid: ");
        // src_v = gsl_matrix_row(final_results_heuristic, 1);
        // show_gsl_vector(&src_v.vector);
        // printf(" max: ");
        // src_v = gsl_matrix_row(final_results_heuristic, 2);
        // show_gsl_vector(&src_v.vector);
        // printf("\ntime: ");
        // src_v = gsl_matrix_row(final_results_heuristic, 3);
        // show_gsl_vector(&src_v.vector);
    }

    // int counter_cpm = 0, counter_heuristic = 0;

    // counter = 0;

    // FILE *fp = NULL;

    // for (i = 0; i <= counter_contraction; i++)
    // {
    //     if (fabs(gsl_matrix_get(final_results, 1, i) - gsl_matrix_get(final_results_heuristic, 1, i)) > 10 * DBLP_EPSILON)
    //     {
    //         // gsl_vector_view src_v;

    //         // printf("\n  %%  ");

    //         // for (j = 0; j <= counter_contraction; j++)
    //         // {
    //         //     printf("% 8.1f%%", delta * j * 100);
    //         // }

    //         // printf("\n< cpm >\n");
    //         // printf(" min: ");
    //         // src_v = gsl_matrix_row(final_results, 0);
    //         // show_gsl_vector(&src_v.vector);
    //         // printf(" mid: ");
    //         // src_v = gsl_matrix_row(final_results, 1);
    //         // show_gsl_vector(&src_v.vector);
    //         // printf(" max: ");
    //         // src_v = gsl_matrix_row(final_results, 2);
    //         // show_gsl_vector(&src_v.vector);
    //         // printf("\ntime: ");
    //         // src_v = gsl_matrix_row(final_results, 3);
    //         // show_gsl_vector(&src_v.vector);

    //         // printf("\n< heuristic >\n");
    //         // printf(" min: ");
    //         // src_v = gsl_matrix_row(final_results_heuristic, 0);
    //         // show_gsl_vector(&src_v.vector);
    //         // printf(" mid: ");
    //         // src_v = gsl_matrix_row(final_results_heuristic, 1);
    //         // show_gsl_vector(&src_v.vector);
    //         // printf(" max: ");
    //         // src_v = gsl_matrix_row(final_results_heuristic, 2);
    //         // show_gsl_vector(&src_v.vector);
    //         // printf("\ntime: ");
    //         // src_v = gsl_matrix_row(final_results_heuristic, 3);
    //         // show_gsl_vector(&src_v.vector);

    //         counter++;

    //         if (fabs(gsl_matrix_get(final_results, 0, i) - gsl_matrix_get(final_results_heuristic, 0, i)) > 10 * DBLP_EPSILON || fabs(gsl_matrix_get(final_results, 2, i) - gsl_matrix_get(final_results_heuristic, 2, i)) > 10 * DBLP_EPSILON)
    //         {
    //             if (gsl_matrix_get(final_results, 0, i) > gsl_matrix_get(final_results_heuristic, 0, i) + DBLP_EPSILON || gsl_matrix_get(final_results, 2, i) < gsl_matrix_get(final_results_heuristic, 2, i) - DBLP_EPSILON)
    //             {
    //                 counter_heuristic++;
    //             }
    //             else
    //             {
    //                 counter_cpm++;
    //             }
    //         }
    //     }
    // }

    // if (counter == 0)
    // {
    //     printf("\nooooooo two sets of mids almost coincide up to tol = %f ooooooo\n\n", 10 * DBLP_EPSILON);
    // }
    // else
    // {
    //     if (counter_heuristic > 0)
    //     {

    //         fp = fopen("./bin/heuristics.log", "a");

    //         fprintf(fp, "\n%s", file_name);

    //         fclose(fp);

    //         printf("\nooooooo statistics (tol = %f): %d different results in mid, cpm outperforms %d times in both min and max, heuristic outperforms %d times in either min or max ooooooo\n\n", 10 * DBLP_EPSILON, counter, counter_cpm, counter_heuristic);
    //     }
    //     else
    //     {
    //         printf("\n******* statistics (tol = %f): %d different results in mid, cpm outperforms %d times in both min and max, heuristic outperforms %d times in either min or max *******\n\n", 10 * DBLP_EPSILON, counter, counter_cpm, counter_heuristic);
    //     }
    // }

    managed::vector_free(increment_P_left);
    managed::vector_free(increment_P_right);
    managed::vector_free(increment_V);
    managed::vector_free(base_xl);
    managed::vector_free(base_xu);
    managed::vector_free(base_yl);
    managed::vector_free(base_yu);
    managed::vector_free(mass_P);
    managed::vector_free(mass_V);
    managed::matrix_free(final_results);
    managed::matrix_free(final_results_heuristic);
    managed::vector_free(whole_solution);

    return 0;
}
int main(int argc,char** argv) {
    try {
        if(argc<2 || std::string(argv[1])=="--help") {
            std::cout<<"Usage: DELTA_solver file.ddt [contractions=5] [verbose=0] [--single|--export|--validate-only] [--fast-chain] [--method=cp|bcm|odba] [--lp-backend=highs] [--rlt-backend=highs] [--bcm-through=LEVEL] [--bcm-early-rlt] [--parallel-directions|--sequential-directions] [--cross-direction-rlt-basis] [--cross-level-rlt-start=basis|transport] [--climb=plain|guarded|augmented] [--row-filter=cheap|hybrid|legacy] [--contraction-path=interior|inclusive] [--contraction-order=inward|outward] [--contraction-reduction=reuse|refresh] [--contraction-precheck=auto|skip|solve] [--contraction-start=default|previous|focal|midpoint|transport|transport-midpoint|ray|secant] [--rlt-start=amc|none|staged] [--rlt-solver=simplex|ipm|choose] [--rlt-presolve=choose|on|off] [--rlt-simplex-strategy=0..4] [--trigger=N] [--cut-limit=N] [--cold-rlt|--trace] [--cplex] [--cplex-time-limit=SECONDS] [--cplex-threads=N] [--cplex-log]\n";
            std::cout << "Path audit: [--contraction-model=fixed-statements|legacy] [--prepare-focal=NEW_FILE] [--focal-file=FILE] [--diagnostics=NEW_JSONL] [--rlt-diagnostics] [--rlt-snapshot-dir=EXISTING_DIR]\n";
            return argc<2?2:0;
        }
        if (!std::getenv("DECIDEIT_METHOD")) setenv("DECIDEIT_METHOD","bcm",1);
        if (!std::getenv("DECIDEIT_PARALLEL_DIRECTIONS") &&
            !std::getenv("DECIDEIT_SEQUENTIAL_DIRECTIONS"))
            setenv("DECIDEIT_PARALLEL_DIRECTIONS","1",1);
        // The production comparison uses five non-point contraction levels:
        // 0%, 20%, 40%, 60%, and 80%; the 100% focal point is evaluated
        // directly.  Keeping the historical default of ten here silently
        // doubled the number of DBLP solves when only a file name was given.
        std::vector<std::string> values={argv[0],argv[1],"5","0","normal"};
        int positional=2;
        bool fast_chain_requested=false;
        for(int k=2;k<argc;k++) {
            std::string arg=argv[k];
            if(arg=="--single")values[4]="single";
            else if(arg.compare(0,20,"--contraction-model=")==0)
                setenv("DECIDEIT_CONTRACTION_MODEL",arg.substr(20).c_str(),1);
            else if(arg.compare(0,16,"--prepare-focal=")==0)
                setenv("DECIDEIT_PREPARE_FOCAL",arg.substr(16).c_str(),1);
            else if(arg.compare(0,13,"--focal-file=")==0)
                setenv("DECIDEIT_FOCAL_FILE",arg.substr(13).c_str(),1);
            else if(arg.compare(0,14,"--diagnostics=")==0)
                setenv("DECIDEIT_DIAGNOSTICS",arg.substr(14).c_str(),1);
            else if(arg=="--export")values[4]="export";
            else if(arg=="--validate-only")values[4]="validate";
            else if(arg=="--trace")setenv("DECIDEIT_TRACE","1",1);
            else if(arg=="--rlt-diagnostics")setenv("DECIDEIT_RLT_DIAGNOSTICS","1",1);
            else if(arg.find("--rlt-snapshot-dir=")==0)setenv("DECIDEIT_RLT_SNAPSHOT_DIR",arg.substr(19).c_str(),1);
            else if(arg=="--cold-rlt")setenv("DECIDEIT_COLD_RLT","1",1);
            else if(arg.compare(0,9,"--method=")==0)
                setenv("DECIDEIT_METHOD",arg.substr(9).c_str(),1);
            else if(arg.compare(0,13,"--lp-backend=")==0)
                setenv("DECIDEIT_LP_BACKEND",arg.substr(13).c_str(),1);
            else if(arg.compare(0,14,"--rlt-backend=")==0)
                setenv("DECIDEIT_RLT_BACKEND",arg.substr(14).c_str(),1);
            else if(arg.compare(0,14,"--bcm-through=")==0)
                setenv("DECIDEIT_BCM_THROUGH",arg.substr(14).c_str(),1);
            else if(arg=="--bcm-early-rlt")
                setenv("DECIDEIT_BCM_EARLY_RLT","1",1);
            else if(arg=="--parallel-directions")
            {
                setenv("DECIDEIT_PARALLEL_DIRECTIONS","1",1);
                unsetenv("DECIDEIT_SEQUENTIAL_DIRECTIONS");
            }
            else if(arg=="--sequential-directions")
            {
                unsetenv("DECIDEIT_PARALLEL_DIRECTIONS");
                setenv("DECIDEIT_SEQUENTIAL_DIRECTIONS","1",1);
            }
            else if(arg=="--fast-chain")fast_chain_requested=true;
            else if(arg=="--cross-direction-rlt-basis")
                setenv("DECIDEIT_CROSS_DIRECTION_RLT_BASIS","1",1);
            else if(arg.compare(0,24,"--cross-level-rlt-start=")==0)
                setenv("DECIDEIT_CROSS_LEVEL_RLT_START",arg.substr(24).c_str(),1);
            else if(arg.compare(0,8,"--climb=")==0)
                setenv("DECIDEIT_CLIMB",arg.substr(8).c_str(),1);
            else if(arg.compare(0,13,"--row-filter=")==0)
                setenv("DECIDEIT_ROW_FILTER",arg.substr(13).c_str(),1);
            else if(arg.compare(0,20,"--contraction-start=")==0)
                setenv("DECIDEIT_CONTRACTION_START",arg.substr(20).c_str(),1);
            else if(arg.compare(0,19,"--contraction-path=")==0)
                setenv("DECIDEIT_CONTRACTION_PATH",arg.substr(19).c_str(),1);
            else if(arg.compare(0,24,"--contraction-reduction=")==0)
                setenv("DECIDEIT_CONTRACTION_REDUCTION",arg.substr(24).c_str(),1);
            else if(arg.compare(0,20,"--contraction-order=")==0)
                setenv("DECIDEIT_CONTRACTION_ORDER",arg.substr(20).c_str(),1);
            else if(arg.compare(0,23,"--contraction-precheck=")==0)
                setenv("DECIDEIT_CONTRACTION_PRECHECK",arg.substr(23).c_str(),1);
            else if(arg.compare(0,12,"--rlt-start=")==0)
                setenv("DECIDEIT_RLT_START",arg.substr(12).c_str(),1);
            else if(arg.compare(0,13,"--rlt-solver=")==0)
                setenv("DECIDEIT_RLT_SOLVER",arg.substr(13).c_str(),1);
            else if(arg.compare(0,15,"--rlt-presolve=")==0)
                setenv("DECIDEIT_RLT_PRESOLVE",arg.substr(15).c_str(),1);
            else if(arg.compare(0,23,"--rlt-simplex-strategy=")==0)
                setenv("DECIDEIT_RLT_SIMPLEX_STRATEGY",arg.substr(23).c_str(),1);
            else if(arg.compare(0,29,"--rlt-simplex-scale-strategy=")==0)
                setenv("DECIDEIT_RLT_SIMPLEX_SCALE_STRATEGY",
                       arg.substr(29).c_str(),1);
            else if(arg.compare(0,25,"--rlt-scale4-min-columns=")==0)
                setenv("DECIDEIT_RLT_SCALE4_MIN_COLUMNS",
                       arg.substr(25).c_str(),1);
            else if(arg.compare(0,25,"--rlt-scale4-max-columns=")==0)
                setenv("DECIDEIT_RLT_SCALE4_MAX_COLUMNS",
                       arg.substr(25).c_str(),1);
            else if(arg.compare(0,10,"--trigger=")==0)
                setenv("DECIDEIT_TRIGGER",arg.substr(10).c_str(),1);
            else if(arg.compare(0,12,"--cut-limit=")==0)
                setenv("DECIDEIT_CUT_LIMIT",arg.substr(12).c_str(),1);
            else if(arg=="--cplex")setenv("DECIDEIT_SOLVER","cplex",1);
            else if(arg=="--cplex-log")setenv("DECIDEIT_CPLEX_LOG","1",1);
            else if(arg.compare(0,19,"--cplex-time-limit=")==0)
                setenv("DECIDEIT_CPLEX_TIMELIMIT",arg.substr(19).c_str(),1);
            else if(arg.compare(0,16,"--cplex-threads=")==0)
                setenv("DECIDEIT_CPLEX_THREADS",arg.substr(16).c_str(),1);
            else if(arg.compare(0,2,"--")==0)throw std::invalid_argument("Unexpected option: "+arg);
            else if(positional<4)values[positional++]=arg;
            else throw std::invalid_argument("Unexpected argument: "+arg);
        }
        // Apply the production bundle after parsing so --fast-chain has the
        // same meaning regardless of its position among experimental flags.
        if(fast_chain_requested)
        {
            setenv("DECIDEIT_METHOD","bcm",1);
            setenv("DECIDEIT_LP_BACKEND","auto",1);
            setenv("DECIDEIT_RLT_BACKEND","highs",1);
            setenv("DECIDEIT_RLT_SIMPLEX_SCALE_STRATEGY","auto",1);
            setenv("DECIDEIT_RLT_SCALE4_MAX_COLUMNS","2500",1);
            unsetenv("DECIDEIT_RLT_SCALE4_MIN_COLUMNS");
            setenv("DECIDEIT_PARALLEL_DIRECTIONS","1",1);
            unsetenv("DECIDEIT_SEQUENTIAL_DIRECTIONS");
            setenv("DECIDEIT_CLIMB","guarded",1);
            setenv("DECIDEIT_ROW_FILTER","cheap",1);
            setenv("DECIDEIT_CONTRACTION_START","transport",1);
            setenv("DECIDEIT_CONTRACTION_PATH","inclusive",1);
            setenv("DECIDEIT_CONTRACTION_ORDER","inward",1);
            setenv("DECIDEIT_CONTRACTION_REDUCTION","refresh",1);
            setenv("DECIDEIT_CONTRACTION_PRECHECK","solve",1);
            setenv("DECIDEIT_RLT_START","staged",1);
            unsetenv("DECIDEIT_COLD_RLT");
            unsetenv("DECIDEIT_BCM_THROUGH");
            unsetenv("DECIDEIT_BCM_EARLY_RLT");
            unsetenv("DECIDEIT_CROSS_DIRECTION_RLT_BASIS");
            unsetenv("DECIDEIT_CROSS_LEVEL_RLT_START");
            unsetenv("DECIDEIT_RLT_SOLVER");
            unsetenv("DECIDEIT_RLT_PRESOLVE");
            unsetenv("DECIDEIT_RLT_SIMPLEX_STRATEGY");
        }
        size_t used=0;int steps=std::stoi(values[2],&used);
        const char* cm=std::getenv("DECIDEIT_CONTRACTION_MODEL");
        if(cm && std::string(cm)!="legacy" && std::string(cm)!="fixed-statements")
            throw std::invalid_argument("contraction model must be fixed-statements or legacy");
        if(cm && std::string(cm)=="legacy" && (std::getenv("DECIDEIT_FOCAL_FILE") || std::getenv("DECIDEIT_PREPARE_FOCAL")))
            throw std::invalid_argument("shared focal artifacts require fixed-statements model");
        if(used!=values[2].size()||steps<1||steps>10000)throw std::invalid_argument("contractions must be 1..10000");
        if(values[3]!="0"&&values[3]!="1")throw std::invalid_argument("verbose must be 0 or 1");
        const std::string method=selected_method();
        if(method!="cp"&&method!="bcm"&&method!="odba")
            throw std::invalid_argument("method must be cp, bcm, or odba");
        const std::string lp_backend=selected_lp_backend();
        if(lp_backend!="highs")
            throw std::invalid_argument("This build requires --lp-backend=highs");
        const std::string rlt_backend=selected_rlt_backend();
        if(rlt_backend!="highs")
            throw std::invalid_argument("This build requires --rlt-backend=highs");
        if(selected_bcm_through() < -1 || selected_bcm_through() >= steps)
            throw std::invalid_argument("bcm through level must be -1..contractions-1");
        const std::string climb=selected_climb();
        if(climb!="plain"&&climb!="guarded"&&climb!="augmented")
            throw std::invalid_argument("climb must be plain, guarded, or augmented");
        const std::string contraction_start=selected_contraction_start();
        if(contraction_start!="default"&&contraction_start!="previous"&&
           contraction_start!="focal"&&contraction_start!="midpoint"&&
           contraction_start!="transport"&&contraction_start!="transport-midpoint"&&
           contraction_start!="ray"&&contraction_start!="secant")
            throw std::invalid_argument("contraction start must be default, previous, focal, midpoint, transport, transport-midpoint, ray, or secant");
        const std::string contraction_path=selected_contraction_path();
        if(contraction_path!="interior"&&contraction_path!="inclusive")
            throw std::invalid_argument("contraction path must be interior or inclusive");
        const std::string contraction_reduction=selected_contraction_reduction();
        if(contraction_reduction!="reuse"&&contraction_reduction!="refresh")
            throw std::invalid_argument("contraction reduction must be reuse or refresh");
        const std::string row_filter=selected_row_filter();
        if(row_filter!="cheap"&&row_filter!="hybrid"&&row_filter!="legacy")
            throw std::invalid_argument("row filter must be cheap, hybrid, or legacy");
        const std::string contraction_order=selected_contraction_order();
        if(contraction_order!="inward"&&contraction_order!="outward")
            throw std::invalid_argument("contraction order must be inward or outward");
        const std::string contraction_precheck=selected_contraction_precheck();
        if(contraction_precheck!="auto"&&contraction_precheck!="skip"&&
           contraction_precheck!="solve")
            throw std::invalid_argument("contraction precheck must be auto, skip, or solve");
        const std::string rlt_start=selected_rlt_start();
        if(rlt_start!="amc"&&rlt_start!="none"&&rlt_start!="staged")
            throw std::invalid_argument("RLT start must be amc, none, or staged");
        if(const char *value=std::getenv("DECIDEIT_RLT_SOLVER")) {
            const std::string solver=value;
            if(solver!="simplex"&&solver!="ipm"&&solver!="choose")
                throw std::invalid_argument("RLT solver must be simplex, ipm, or choose");
        }
        if(const char *value=std::getenv("DECIDEIT_RLT_PRESOLVE")) {
            const std::string presolve=value;
            if(presolve!="choose"&&presolve!="on"&&presolve!="off")
                throw std::invalid_argument("RLT presolve must be choose, on, or off");
        }
        if(const char *value=std::getenv("DECIDEIT_RLT_SIMPLEX_STRATEGY")) {
            const int strategy=std::stoi(value);
            if(strategy<0||strategy>4)
                throw std::invalid_argument("RLT simplex strategy must be 0..4");
        }
        if(const char *value=std::getenv(
               "DECIDEIT_RLT_SIMPLEX_SCALE_STRATEGY")) {
            const std::string scale=value;
            if(scale!="auto"&&scale!="adaptive") {
                size_t parsed=0;
                const int strategy=std::stoi(scale,&parsed);
                if(parsed!=scale.size()||strategy<0||strategy>4)
                    throw std::invalid_argument(
                        "RLT simplex scale strategy must be auto, adaptive, or 0..4");
            }
        }
        for(const char *name : {"DECIDEIT_RLT_SCALE4_MIN_COLUMNS",
                                "DECIDEIT_RLT_SCALE4_MAX_COLUMNS"}) {
            if(const char *value=std::getenv(name)) {
                size_t parsed=0;
                const long long columns=std::stoll(value,&parsed);
                if(parsed!=std::strlen(value)||columns<1||columns>4000000)
                    throw std::invalid_argument(
                        std::string(name)+" must be 1..4000000");
            }
        }
        if(const char *value=std::getenv("DECIDEIT_CROSS_LEVEL_RLT_START")) {
            const std::string start=value;
            if(start!="basis"&&start!="transport")
                throw std::invalid_argument("cross-level RLT start must be basis or transport");
        }
        if(std::getenv("DECIDEIT_PARALLEL_DIRECTIONS") &&
           std::getenv("DECIDEIT_CROSS_DIRECTION_RLT_BASIS"))
            throw std::invalid_argument(
                "cross-direction RLT basis requires sequential directions");
        if(selected_trigger()<1||selected_trigger()>10000)
            throw std::invalid_argument("trigger must be 1..10000");
        if(selected_cut_limit()<1||selected_cut_limit()>10000)
            throw std::invalid_argument("cut limit must be 1..10000");
        std::vector<char*> args;for(auto& v:values)args.push_back(&v[0]);
        gsl_set_error_handler_off();
        return run_cli(args.size(),args.data());
    } catch(const std::exception& e) {
        decision_path::log_error(e.what());
        std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;
    }
}
