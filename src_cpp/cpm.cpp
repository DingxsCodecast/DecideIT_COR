#include "dblp.h"

static void grow(gsl_vector*& v,size_t size) {
    if(size<=v->size)return;
    gsl_vector* next=managed::vector_calloc(size);
    if(!next)throw std::bad_alloc();
    gsl_vector_view part=gsl_vector_subvector(next,0,v->size);
    gsl_vector_memcpy(&part.vector,v);managed::vector_free(v);v=next;
}
static void grow(gsl_vector_int*& v,size_t size) {
    if(size<=v->size)return;
    gsl_vector_int* next=managed::vector_int_calloc(size);
    if(!next)throw std::bad_alloc();
    gsl_vector_int_view part=gsl_vector_int_subvector(next,0,v->size);
    gsl_vector_int_memcpy(&part.vector,v);managed::vector_int_free(v);v=next;
}
void ensure_rlt_capacity(LinearizationParams* p,size_t rows,size_t nz) {
    if(rows>2000000 || nz>32000000)throw std::runtime_error("RLT capacity limit reached");
    grow(p->row_lower,rows);grow(p->row_upper,rows);grow(p->a_start,rows+1);
    grow(p->a_value,nz);grow(p->a_index,nz);
}

/* set the objective function */
void construct_linearization_objective(int n_alt1, int n_alt2, LinearizationParams *LinPara)
{
    int n_variable = n_alt1 + n_alt2;
    int i, t = 2 * n_variable;

    for (i = 0; i < n_alt1; i++)
    {
        gsl_vector_set(LinPara->col_cost, t, 1);
        t = t + n_variable + 1;
    }

    for (i = 0; i < n_alt2; i++)
    {
        gsl_vector_set(LinPara->col_cost, t, -1);
        t = t + n_variable + 1;
    }
}

void initialize_LinearizationParameters(LinearizationParams *LinPara)
{
    LinPara->m_total = 0;
    LinPara->n_total = 0;

    LinPara->n_XY = 0;
    LinPara->n_XX = 0;
    LinPara->n_YY = 0;

    LinPara->num_col = 0;
    LinPara->num_row = 0;
    LinPara->num_nz = 0;
    LinPara->offset = 0;

    LinPara->col_cost = NULL;

    LinPara->row_lower = NULL;
    LinPara->row_upper = NULL;
    LinPara->a_value = NULL;
    LinPara->a_index = NULL;
    LinPara->a_start = NULL;

    LinPara->col_lower = NULL;
    LinPara->col_upper = NULL;

    LinPara->a_start_column = 0;
    LinPara->a_index_column = 0;
    LinPara->cut_cross_rows.clear();
}

void destroy_LinearizationParameters(LinearizationParams *p)
{
    if(p->col_cost) managed::vector_free(p->col_cost);
    if(p->row_lower) managed::vector_free(p->row_lower);
    if(p->row_upper) managed::vector_free(p->row_upper);
    if(p->a_value) managed::vector_free(p->a_value);
    if(p->a_index) managed::vector_int_free(p->a_index);
    if(p->a_start) managed::vector_int_free(p->a_start);
    if(p->col_lower) managed::vector_free(p->col_lower);
    if(p->col_upper) managed::vector_free(p->col_upper);
    p->persistent.reset();
    initialize_LinearizationParameters(p);
}
LinearizationParams::LinearizationParams(){initialize_LinearizationParameters(this);}
LinearizationParams::~LinearizationParams(){destroy_LinearizationParameters(this);}

#if 0 // Archived XX/YY full RLT; the paper implementation is cross-term-only.
void construct_full_LinearizationParameters(
    const int mX,
    const int mY,
    const int nX,
    const int nY,
    LinearizationParams *linPara)
{
    int t;

    linPara->m_total += mX;
    linPara->m_total += mY;
    linPara->m_total += nX * (mX + mY);
    linPara->m_total += nX * (mX + mY);
    linPara->m_total += nY * (mX + mY);
    linPara->m_total += nY * (mX + mY);
    linPara->m_total += nSym(nX);
    linPara->m_total += nSym(nY);
    linPara->m_total += nSym(nX);
    linPara->m_total += nSym(nY);
    linPara->m_total += nX * nY;
    linPara->m_total += nX * nY;
    linPara->m_total += nX * nY;
    linPara->m_total += nX * nY;
    linPara->m_total += nX * nX;
    linPara->m_total += nY * nY;

    if (mX + mY >= 2)
    {
        t = mX + mY;
        linPara->m_total += mX + mY + (t * (t - 1)) / 2;
    }
    else
    {
        linPara->m_total += mX + mY;
    }

    int m_eq_P = 2;

    linPara->m_total += m_eq_P + nX * m_eq_P + nY * m_eq_P;

    linPara->n_total = 0;

    linPara->n_total += nX;
    linPara->n_total += nY;
    linPara->n_total += nX * nY;
    linPara->n_total += nSym(nX);
    linPara->n_total += nSym(nY);

    linPara->n_XY = nX + nY;
    linPara->n_XX = linPara->n_XY + nX * nY;
    linPara->n_YY = linPara->n_XX + nSym(nX);

    linPara->col_cost = managed::vector_calloc(linPara->n_total);

    // Initial RLT rows have at most 16 nonzeros except normalization
    // products (covered by 4*nX*nY). New cuts add 2*nY+1 rows.
    const size_t rows = static_cast<size_t>(linPara->m_total) + 1;
    const size_t nz = 16 * rows + 4 * static_cast<size_t>(nX) * nY;
    if (rows > 2000000 || nz > 32000000)
        throw std::runtime_error("RLT resource limit exceeded before allocation");
    linPara->row_lower = managed::vector_calloc(rows);
    linPara->row_upper = managed::vector_calloc(rows);
    linPara->a_value = managed::vector_calloc(nz);
    linPara->a_index = managed::vector_int_calloc(nz);
    linPara->a_start = managed::vector_int_calloc(rows + 1);

    linPara->col_lower = managed::vector_calloc(linPara->n_total);
    linPara->col_upper = managed::vector_calloc(linPara->n_total);
}
#endif
void construct_LinearizationParameters(
    const int mX, const int mY, const int nX, const int nY,
    LinearizationParams *linPara)
{
    if (nX != nY || nX < 1)
        throw std::invalid_argument("Cross-term RLT requires equal positive block dimensions");
    linPara->n_XY = nX + nY;
    linPara->n_XX = linPara->n_XY + nX * nY;
    linPara->n_YY = linPara->n_XX;
    linPara->n_total = linPara->n_XX;

    const size_t rows = static_cast<size_t>(mX + mY) +
        4ULL * nX * nY + 2ULL * nY * mX + 2ULL * nX * mY +
        2ULL + 2ULL * nX + 1ULL;
    const size_t nonzeros = 16ULL * rows + 4ULL * nX * nY;
    if (rows > 2000000 || nonzeros > 32000000)
        throw std::runtime_error("Cross-term RLT resource limit exceeded before allocation");

    linPara->m_total = static_cast<int>(rows - 1);
    linPara->col_cost = managed::vector_calloc(linPara->n_total);
    linPara->row_lower = managed::vector_calloc(rows);
    linPara->row_upper = managed::vector_calloc(rows);
    linPara->a_value = managed::vector_calloc(nonzeros);
    linPara->a_index = managed::vector_int_calloc(nonzeros);
    linPara->a_start = managed::vector_int_calloc(rows + 1);
    linPara->col_lower = managed::vector_calloc(linPara->n_total);
    linPara->col_upper = managed::vector_calloc(linPara->n_total);
}

int check_PBase(const d_row P_lobo,
                const d_row P_upbo,
                const int n_alt1,
                const int n_alt2,
                const int Ai_begin,
                const int Aj_begin)
{
    int i, counter = 0;

    /* check whether the interval imposed on each variable is too tight */

    for (i = 0; i < n_alt1; i++)
    {
        if (P_upbo[Ai_begin + i] == P_lobo[Ai_begin + i])
        {
            counter++;
        }
    }

    for (i = 0; i < n_alt2; i++)
    {
        if (P_upbo[Aj_begin + i] == P_lobo[Aj_begin + i])
        {
            counter++;
        }
    }

    return counter;
}

void dataset_constructor(
    struct dataset *p,
    int noc,
    int nov_x,
    int nov_y,
    int CUTS_MAX)
{
    p->valid_A = managed::matrix_calloc(noc + CUTS_MAX, nov_x + 1);
    p->xl = managed::vector_calloc(nov_x);
    p->xu = managed::vector_calloc(nov_x);
    p->Q = managed::matrix_calloc(nov_x, nov_y);
    p->dy = managed::vector_calloc(nov_y);

    gsl_vector_set_all(p->xu, 1);
}

void dataset_destructor(struct dataset *p)
{
    managed::matrix_free(p->valid_A);
    managed::vector_free(p->xl);
    managed::vector_free(p->xu);
    managed::matrix_free(p->Q);
    managed::vector_free(p->dy);
    managed::destroy(p);
    p = nullptr;
    // mem_free(p, "exit_cpm", 4);
}

void global_solution_constructor(
    struct global_solution *p,
    int dim_x,
    int dim_y)
{
    p->alpha = 1e7;
    p->x_solution = managed::vector_calloc(dim_x);
    p->y_solution = managed::vector_calloc(dim_y);
}

void global_solution_destructor(struct global_solution *p)
{
    managed::vector_free(p->x_solution);
    managed::vector_free(p->y_solution);
    managed::destroy(p);
    p = nullptr;
    // mem_free(p, "exit_cpm", 2);
}

void local_solution_constructor(struct local_solution *p,
                                int dim_x,
                                int dim_y)
{
    p->cbov = 1e7;
    p->adj_vertices = managed::matrix_calloc(dim_x, adjacency_capacity(dim_x));
    p->n_adj_vertices = 0;
    p->degenerate_vertex = false;
    p->x_solution = managed::vector_calloc(dim_x);
    p->y_solution = managed::vector_calloc(dim_y);
}

void local_solution_destructor(struct local_solution *p)
{
    managed::matrix_free(p->adj_vertices);
    managed::vector_free(p->x_solution);
    managed::vector_free(p->y_solution);
    managed::destroy(p);
    p = nullptr;
    // mem_free(p, "exit_cpm", 3);
}

void dblp_constructor(struct dblp *p,
                      int *dimension)
{
    p->Ax = managed::matrix_calloc(dimension[2], dimension[3]);
    p->bx = managed::vector_calloc(dimension[2]);
    p->xl = managed::vector_calloc(dimension[3]);
    p->xu = managed::vector_calloc(dimension[3]);
    gsl_vector_set_all(p->xl, 0);
    gsl_vector_set_all(p->xu, 1);

    p->Q = managed::matrix_calloc(dimension[1], dimension[3]);
}

void dblp_destructor(struct dblp *p)
{
    managed::matrix_free(p->Ax);
    managed::vector_free(p->bx);
    managed::vector_free(p->xl);
    managed::vector_free(p->xu);
    managed::matrix_free(p->Q);
    managed::destroy(p);
    p = nullptr;
    // mem_free(p, "exit_cpm", 1);
}

void free_all_fixed_variables(LpModel *lp_Y,
                              struct dblp *instance,
                              gsl_vector *temp_v,
                              gsl_vector *temp_coef,
                              gsl_vector_int *effective_constraints,
                              gsl_vector_int *fixed_variables)
{
    free_lp(&lp_Y);
    dblp_destructor(instance);
    managed::vector_free(temp_v);
    managed::vector_free(temp_coef);
    managed::vector_int_free(effective_constraints);
    managed::vector_int_free(fixed_variables);
}

void solve_all_fixed_variables(
    LpModel *lp_Y,
    struct dblp *instance,
    d_row V_min_point,
    d_row P_min_point,
    int Ai_begin,
    int Ai_end,
    int Aj_begin,
    int Aj_end,
    d_row V_max_point,
    d_row P_max_point,
    int n_variable_y)
{
    int i, j;

    gsl_vector *coef_y = managed::vector_calloc(n_variable_y + 1);
    gsl_vector *solution_y = managed::vector_calloc(n_variable_y);

    gsl_vector_view dest_v;

    dest_v = gsl_vector_subvector(coef_y, 1, n_variable_y);

    set_minim(lp_Y);

    gsl_vector_set_zero(coef_y);
    gsl_blas_dgemv(CblasTrans, 1, instance->Q, instance->xl, 1, &dest_v.vector); // the coefficient of y, Q' * x
    set_obj_fn(lp_Y, coef_y->data);
    solve(lp_Y);

    V_min_point[0] = 0;
    P_min_point[0] = get_objective(lp_Y);
    get_variables(lp_Y, solution_y->data);

    j = 0;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        V_min_point[i] = gsl_vector_get(instance->xl, j);
        P_min_point[i] = gsl_vector_get(solution_y, j++);
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        V_min_point[i] = gsl_vector_get(instance->xl, j);
        P_min_point[i] = gsl_vector_get(solution_y, j++);
    }

    gsl_matrix_scale(instance->Q, -1);
    gsl_vector_set_zero(coef_y);
    gsl_blas_dgemv(CblasTrans, 1, instance->Q, instance->xl, 1, &dest_v.vector); // the coefficient of y, Q' * x
    set_obj_fn(lp_Y, coef_y->data);
    solve(lp_Y);

    V_max_point[0] = 0;
    P_max_point[0] = get_objective(lp_Y);
    get_variables(lp_Y, solution_y->data);

    j = 0;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        V_max_point[i] = gsl_vector_get(instance->xl, j);
        P_max_point[i] = gsl_vector_get(solution_y, j++);
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        V_max_point[i] = gsl_vector_get(instance->xl, j);
        P_max_point[i] = gsl_vector_get(solution_y, j++);
    }

    managed::vector_free(coef_y);
    managed::vector_free(solution_y);
}

void free_local_variables(struct local_solution *local_amc,
                          LpModel *lp_Y,
                          struct dblp *instance,
                          struct global_solution *solution,
                          gsl_vector *temp_v,
                          gsl_vector *temp_coef,
                          gsl_vector_int *effective_constraints,
                          gsl_vector_int *fixed_variables,
                          struct dataset *dataset_min,
                          struct dataset *dataset_max,
                          gsl_vector_int *unfixed_indices,
                          LpModel *lp_X_std_min,
                          LpModel *lp_X_std_max,
                          gsl_vector *coef_x)
{
    local_solution_destructor(local_amc);
    free_lp(&lp_Y);
    dblp_destructor(instance);
    global_solution_destructor(solution);
    managed::vector_free(temp_v);
    managed::vector_free(temp_coef);
    managed::vector_int_free(effective_constraints);
    managed::vector_int_free(fixed_variables);
    dataset_destructor(dataset_min);
    dataset_destructor(dataset_max);
    managed::vector_int_free(unfixed_indices);
    free_lp(&lp_X_std_min);
    free_lp(&lp_X_std_max);
    managed::vector_free(coef_x);
}

void dataset_copy(struct dataset *dataset_max,
                  struct dataset *dataset_min)
{
    if (dataset_max->valid_A->size1 != dataset_min->valid_A->size1)
    {
        managed::matrix_free(dataset_max->valid_A);
        dataset_max->valid_A = managed::matrix_calloc(dataset_min->valid_A->size1, dataset_min->valid_A->size2);
    }

    gsl_matrix_memcpy(dataset_max->valid_A, dataset_min->valid_A);
    gsl_matrix_memcpy(dataset_max->Q, dataset_min->Q);
    gsl_vector_memcpy(dataset_max->xl, dataset_min->xl);
    gsl_vector_memcpy(dataset_max->xu, dataset_min->xu);
    gsl_vector_memcpy(dataset_max->dy, dataset_min->dy);
    dataset_max->ori_num_cons = dataset_min->ori_num_cons;
}

void construct_parameters(
    struct dblp *instance,
    gsl_vector_int *effective_constraints,
    gsl_vector_int *fixed_variables,
    struct dataset *dataset_in,
    int num_alt_i)
{
    const int n_variable = instance->xl->size; // the original # of variables in X
    const int n_fixed_variable = gsl_vector_int_get(fixed_variables, 0);
    int i, j, which_row;
    REAL entry_value;

    gsl_vector_view src_v, dest_v;

    dataset_in->ori_num_cons = gsl_vector_int_get(effective_constraints, 0);

    if (n_fixed_variable == 0) // establish all parameters directly as no variable is fixed
    {
        for (i = 0; i < gsl_vector_int_get(effective_constraints, 0); i++)
        {
            which_row = gsl_vector_int_get(effective_constraints, i + 1);
            src_v = gsl_matrix_row(instance->Ax, which_row);
            dest_v = gsl_matrix_subrow(dataset_in->valid_A, i, 1, n_variable);
            gsl_vector_memcpy(&dest_v.vector, &src_v.vector);
            gsl_matrix_set(dataset_in->valid_A, i, 0, gsl_vector_get(instance->bx, which_row));
        }

        gsl_vector_memcpy(dataset_in->xu, instance->xu);
        gsl_vector_memcpy(dataset_in->xl, instance->xl);
        gsl_matrix_memcpy(dataset_in->Q, instance->Q);
    }
    else // some variables are fixed
    {
        /*
         * Compact each effective row directly into the reduced variable space.
         * Keeping the preallocated matrix is important: its extra CUTS_MAX rows
         * are used later by the cutting-plane controller.
         */
        const int original_effective = gsl_vector_int_get(effective_constraints, 0);
        int reduced_rows = 0;
        gsl_matrix_set_zero(dataset_in->valid_A);

        for (i = 0; i < original_effective; i++)
        {
            which_row = gsl_vector_int_get(effective_constraints, i + 1);
            double rhs = gsl_vector_get(instance->bx, which_row);
            int reduced_col = 1;
            bool has_unfixed_coefficient = false;

            for (j = 1; j <= n_variable; j++)
            {
                entry_value = gsl_matrix_get(instance->Ax, which_row, j - 1);
                if (gsl_vector_int_get(fixed_variables, j) != 0)
                {
                    rhs -= entry_value * gsl_vector_get(instance->xl, j - 1);
                }
                else
                {
                    gsl_matrix_set(dataset_in->valid_A, reduced_rows, reduced_col++, entry_value);
                    has_unfixed_coefficient = has_unfixed_coefficient || std::fabs(entry_value) > DBLP_EPSILON;
                }
            }

            /*
             * A zero row with a nonnegative right-hand side is redundant.
             * A zero row with a negative right-hand side must be retained so
             * that the feasibility LP reports the model as infeasible.
             */
            if (!has_unfixed_coefficient && rhs >= -DBLP_EPSILON)
            {
                dest_v = gsl_matrix_row(dataset_in->valid_A, reduced_rows);
                gsl_vector_set_zero(&dest_v.vector);
                continue;
            }

            gsl_matrix_set(dataset_in->valid_A, reduced_rows, 0, rhs);
            gsl_vector_int_set(effective_constraints, reduced_rows + 1, which_row);
            reduced_rows++;
        }

        gsl_vector_int_set(effective_constraints, 0, reduced_rows);
        for (i = reduced_rows + 1; i < static_cast<int>(effective_constraints->size); i++)
            gsl_vector_int_set(effective_constraints, i, 0);
        dataset_in->ori_num_cons = reduced_rows;

        i = 0;

        /* construct Q, lower and upper bounds, and dy in dataset_in */
        for (j = 1; j <= n_variable; j++)
        {
            if (gsl_vector_int_get(fixed_variables, j) == 0) // the variable is unfixed
            {
                src_v = gsl_matrix_row(instance->Q, j - 1);
                gsl_matrix_set_row(dataset_in->Q, i, &src_v.vector);
                gsl_vector_set(dataset_in->xl, i, gsl_vector_get(instance->xl, j - 1));
                gsl_vector_set(dataset_in->xu, i, gsl_vector_get(instance->xu, j - 1));
                i++;
            }
            else // the variable is fixed
            {
                if (j <= num_alt_i)
                {
                    gsl_vector_set(dataset_in->dy, j - 1, gsl_vector_get(instance->xl, j - 1)); // the entry in Q is 1
                }
                else
                {
                    gsl_vector_set(dataset_in->dy, j - 1, -gsl_vector_get(instance->xl, j - 1)); // the entry in Q is -1
                }
            }
        }
    }

    // show_gsl_matrix(dataset_in->Q);
    // show_gsl_vector(dataset_in->dy);
}

void construct_standard_lp(
    struct dataset *dataset_in,
    LpModel *lp_X_std,
    const int noc)
{
    const int n_variable = get_Ncolumns(lp_X_std);
    int i, j;
    REAL result;

    gsl_vector *coef = managed::vector_calloc(1 + n_variable);
    gsl_vector_view src_v, dest_v;

    dest_v = gsl_vector_subvector(coef, 1, n_variable);

    for (i = 0; i < noc; i++)
    {
        src_v = gsl_matrix_subrow(dataset_in->valid_A, i, 1, n_variable);
        gsl_vector_memcpy(&dest_v.vector, &src_v.vector);

        set_row(lp_X_std, i + 1, coef->data);
        set_constr_type(lp_X_std, i + 1, kLpLe);
        gsl_blas_ddot(&src_v.vector, dataset_in->xl, &result);
        set_rh(lp_X_std, i + 1, gsl_matrix_get(dataset_in->valid_A, i, 0) - result); // rhs = b_i - A_i' * xl
    }

    j = noc + 1;

    for (i = 0; i < n_variable; i++)
    {
        set_lowbo(lp_X_std, i + 1, 0);
        set_upbo(lp_X_std, i + 1, get_infinite(lp_X_std));

        gsl_vector_set_zero(coef);
        gsl_vector_set(coef, i + 1, 1);
        set_row(lp_X_std, j, coef->data);
        set_constr_type(lp_X_std, j, kLpLe);
        set_rh(lp_X_std, j, gsl_vector_get(dataset_in->xu, i) - gsl_vector_get(dataset_in->xl, i));
        j++;
    }

    managed::vector_free(coef);
}

LpModel *construct_Y(
    int n_variable_y,
    int *mY,
    std::vector<LinCons> &LC_P,
    int num_alt_i,
    int num_alt_j,
    const int Ai_begin,
    const int Ai_end,
    const int Aj_begin,
    const int Aj_end,
    const d_row P_lobo,
    const d_row P_upbo)
{
    int i, j;
    LpModel *lp_Y = make_lp(*mY + 2, n_variable_y);

    /* construct Y, which does not change along the entire algorithm */
    set_verbose(lp_Y, kLpNeutral);
    set_minim(lp_Y);

    gsl_vector *temp_coef = managed::vector_calloc(n_variable_y + 1);

    // set linear constraints in P-base for lp_Y
    if (*mY > 0)
    {
        for (i = 0; i < *mY; i++)
        {
            gsl_vector_set_zero(temp_coef);
            for (size_t k = 0; k < LC_P[i].index.size(); ++k)
                gsl_vector_set(temp_coef, LC_P[i].index[k], LC_P[i].coefficient[k]);
            set_row(lp_Y, i + 1, temp_coef->data);
            set_constr_type(lp_Y, i + 1, kLpLe);
            set_rh(lp_Y, i + 1, LC_P[i].b);
        }

        i++;
    }
    else
    {
        i = 1;
    }

    // set two normalization constraints in for lp_Y
    gsl_vector_set_zero(temp_coef);

    for (j = 1; j <= num_alt_i; j++)
    {
        gsl_vector_set(temp_coef, j, 1);
    }

    set_row(lp_Y, i, temp_coef->data);
    set_constr_type(lp_Y, i, kLpEq);
    set_rh(lp_Y, i++, 1);

    gsl_vector_set_zero(temp_coef);

    for (j = 1; j <= num_alt_j; j++)
    {
        gsl_vector_set(temp_coef, num_alt_i + j, 1);
    }

    set_row(lp_Y, i, temp_coef->data);
    set_constr_type(lp_Y, i, kLpEq);
    set_rh(lp_Y, i++, 1);

    // set lower and upper bounds in P-base for lp_Y
    j = 1;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        set_lowbo(lp_Y, j, P_lobo[i]);
        set_upbo(lp_Y, j, P_upbo[i]);
        j++;
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        set_lowbo(lp_Y, j, P_lobo[i]);
        set_upbo(lp_Y, j, P_upbo[i]);
        j++;
    }

    managed::vector_free(temp_coef);
    return lp_Y;
    // gsl_vector_set_all(temp_coef, 1);
    // set_obj_fn(lp_Y, temp_coef->data);

    // return_value = solve(lp_Y);

    // if (return_value == 2) // if solve(lp_Y) != 0, it is always assumed that the feasible region is infeasible
    // {
    //     free_lp(&lp_Y);
    //     dblp_destructor(instance);
    //     managed::vector_free(temp_coef);

    //     return INFEASIBLE;
    // }
}

/*
    cpm_solver() solves a minimizing DBLP with a cutting plane method without degeneracy removal and lower bounding technique

    inputs
        P_lobo, P_upbo, V_lobo, V_upbo: lower and upper bounds for P-/V-bases
        P_links, V_links: linear constraints for P-/V-bases, P_links contains no normalization constraints
        P_min_point, V_min_point: final P and V results with P_min_point[0] containing the objective value of minimizing problem
        P_max_point, V_max_point: final P and V results with P_max_point[0] containing the objective value of maximizing problem
        Ai_begin, Ai_end, Aj_begin, Aj_end: starting and ending indices for alternative i and j in P_lobo, P_upbo, V_lobo, V_upbo

    outputs
        P_min_point, V_min_point: final P and V results with P_min_point[0] containing the objective value of minimizing problem
        P_max_point, V_max_point: final P and V results with P_max_point[0] containing the objective value of maximizing problem
        P_lobo[0]/P_upbo[0]
            00      default value, everything is ok
            -1      memory allocation failed
            -2      input error
            -3      LP subproblem failed
            -4      GSL LU_decomp/det/invert/solve failed
            -5      CUT_LIMIT
	            -6      degeneracy
	            -9      no transition column at the starting degenerate basis
*/

static int cpm_solver_impl(
    d_row P_lobo, d_row P_upbo, const std::vector<LinCons> &P_rows,
    d_row V_lobo, d_row V_upbo, const std::vector<LinCons> &V_rows,
    d_row P_min_point, d_row V_min_point, d_row P_max_point, d_row V_max_point,
    int Ai_begin, int Ai_end, int Aj_begin, int Aj_end,
    int CUTS_MAX,
    std::vector<LinCons> &LC_V, std::vector<LinCons> &LC_P,
    int *mX, int *mY,
    const double *P_min_start,
    const double *P_max_start,
    RltCrossLevelState *min_rlt_state,
    RltCrossLevelState *max_rlt_state)
{
    gsl_set_error_handler_off();
    reset_lp_backend_stats();

    int exitFlag;

    /*
        X contains V-base, and Y constaints P-base.
        dimension[2] and [3] contain the rows and columns of linear constraints in V-base;
        dimension[0] and [1] contain the rows and columns of linear constraints in P-base.
    */

    /* construct an instance to save the original data for the current DBLP */
    int dimension[4], i, j, counter, return_value;

    dimension[3] = Ai_end - Ai_begin + 1 + Aj_end - Aj_begin + 1;
    dimension[1] = dimension[3]; // the # of variables in P- and V-base are identical

    if (*mX == -1)
    {
        dimension[2] = static_cast<int>(V_rows.size());
        dimension[0] = static_cast<int>(P_rows.size()) + 2;
    }
    else
    {
        dimension[0] = *mY + 2;
        dimension[2] = *mX;
    }

    struct dblp *instance = managed::create<dblp>();

    // struct dblp *instance = (struct dblp *)mem_alloc(sizeof(struct dblp), "dblp", "cpm", 1, FALSE);

    if (!instance)
    {
        P_lobo[0] = -1;
        P_upbo[0] = -1;

        return (OUT_OF_MEMORY); // allocation failed
    }

    dblp_constructor(instance, dimension);

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

    // set linear constraints in V-base
#if 0
    if (*mX == -1)
    {
        for (i = 1; i <= V_related_constraints[0]; i++)
        {
            if (V_links[V_related_constraints[i]].cmp == '<')
            {
                gsl_vector_set(instance->bx, i - 1, V_links[V_related_constraints[i]].value);
                temp_sign = 1;
            }
            else
            {
                gsl_vector_set(instance->bx, i - 1, -V_links[V_related_constraints[i]].value);
                temp_sign = -1;
            }

            head = V_links[V_related_constraints[i]].head;

            if ((head >= Ai_begin) && (head <= Ai_end))
            {
                gsl_matrix_set(instance->Ax, i - 1, head - Ai_begin, temp_sign);
            }
            else if ((head >= Aj_begin) && (head <= Aj_end))
            {
                gsl_matrix_set(instance->Ax, i - 1, head - Aj_begin + num_alt_i, temp_sign);
            }

            tail = V_links[V_related_constraints[i]].tail;

            if ((tail >= Ai_begin) && (tail <= Ai_end))
            {
                gsl_matrix_set(instance->Ax, i - 1, tail - Ai_begin, -temp_sign);
            }
            else if ((tail >= Aj_begin) && (tail <= Aj_end))
            {
                gsl_matrix_set(instance->Ax, i - 1, tail - Aj_begin + num_alt_i, -temp_sign);
            }
        }
    }
    else
    {
        for (i = 0; i < *mX; i++)
        {
            gsl_matrix_set(instance->Ax, i, LC_V[i].j1 - 1, 1.0);
            gsl_matrix_set(instance->Ax, i, LC_V[i].j2 - 1, -1.0);
            gsl_vector_set(instance->bx, i, LC_V[i].b);
        }
    }
#endif
    const std::vector<LinCons> &active_v = (*mX == -1) ? V_rows : LC_V;
    for (size_t row = 0; row < active_v.size(); ++row)
    {
        for (size_t k = 0; k < active_v[row].index.size(); ++k)
            gsl_matrix_set(instance->Ax, row, active_v[row].index[k] - 1,
                           active_v[row].coefficient[k]);
        gsl_vector_set(instance->bx, row, active_v[row].b);
    }

    // construct two LPs for X and Y, respectively
    int n_variable_x = instance->Ax->size2;
    const int n_variable_y = n_variable_x;

    LpModel *lp_Y = NULL;

    gsl_vector *temp_coef = managed::vector_calloc(n_variable_y + 1);

    gsl_vector *temp_v = managed::vector_calloc(n_variable_x + 1);
    gsl_vector_int *effective_constraints = managed::vector_int_calloc(instance->Ax->size1 + 1);
    gsl_vector_int *fixed_variables = managed::vector_int_calloc(n_variable_x + 1);

    // Along an inward contraction path the feasible box only shrinks.  A
    // general linear row proved redundant at the preceding level therefore
    // remains redundant.  Reuse that active set instead of resolving one LP
    // per row at every level.  The first level and standalone queries still
    // run the selected row-filter policy.
    const bool reuse_reduction = std::getenv("DECIDEIT_REUSE_REDUCTION") &&
        *mX >= 0 && *mY >= 0;
    if (reuse_reduction)
    {
        gsl_vector_int_set_zero(fixed_variables);
        gsl_vector_int_set(effective_constraints, 0,
                           static_cast<int>(instance->Ax->size1));
        for (size_t row = 0; row < instance->Ax->size1; ++row)
            gsl_vector_int_set(effective_constraints, row + 1,
                               static_cast<int>(row));
    }
    else
    {
        redundancy(instance, effective_constraints, fixed_variables);
    }

    if (gsl_vector_int_get(fixed_variables, 0) == n_variable_x) // solve the dblp problem directly as all variables are fixed
    {
        P_lobo[0] = 0;
        P_upbo[0] = 0;

        if (check_PBase(P_lobo, P_upbo, num_alt_i, num_alt_j, Ai_begin, Aj_begin) == n_variable_y)
        {
            double obj_value = 0.0;
            double temp_v, temp_p;

            j = 1;

            for (i = Ai_begin; i <= Ai_end; i++)
            {
                temp_v = (V_lobo[i] + V_upbo[i]) / 2;
                temp_p = (P_lobo[i] + P_upbo[i]) / 2;
                obj_value = obj_value + temp_v * temp_p;

                V_min_point[j] = temp_v;
                P_min_point[j] = temp_p;
                V_max_point[j] = temp_v;
                P_max_point[j] = temp_p;

                j++;
            }

            for (i = Aj_begin; i <= Aj_end; i++)
            {
                temp_v = (V_lobo[i] + V_upbo[i]) / 2;
                temp_p = (P_lobo[i] + P_upbo[i]) / 2;
                obj_value = obj_value - temp_v * temp_p;

                V_min_point[j] = temp_v;
                P_min_point[j] = temp_p;
                V_max_point[j] = temp_v;
                P_max_point[j] = temp_p;

                j++;
            }

            V_min_point[0] = 0;
            P_min_point[0] = obj_value;
            V_max_point[0] = 0;
            P_max_point[0] = -obj_value;
        }
        else
        {
            if (*mY == -1)
            {
                LC_P = P_rows;
                *mY = static_cast<int>(LC_P.size());
            }

            lp_Y = construct_Y(n_variable_y, mY, LC_P, num_alt_i, num_alt_j, Ai_begin, Ai_end, Aj_begin, Aj_end, P_lobo, P_upbo);
            solve_all_fixed_variables(lp_Y, instance, V_min_point, P_min_point, Ai_begin, Ai_end, Aj_begin, Aj_end, V_max_point, P_max_point, n_variable_y);
        }

        free_all_fixed_variables(lp_Y, instance, temp_v, temp_coef, effective_constraints, fixed_variables);

        return OK;
    }

    if (!reuse_reduction)
        DataPreparation(V_rows, effective_constraints, P_rows, P_lobo, P_upbo,
                        Ai_begin, Ai_end, Aj_begin, Aj_end, LC_V, LC_P, mX, mY);

    lp_Y = construct_Y(n_variable_y, mY, LC_P, num_alt_i, num_alt_j, Ai_begin, Ai_end, Aj_begin, Aj_end, P_lobo, P_upbo);

    LinearizationParams LinPara;
    initialize_LinearizationParameters(&LinPara);
    LinearizationParams LinParaMax;
    initialize_LinearizationParameters(&LinParaMax);

    int noc = gsl_vector_int_get(effective_constraints, 0);

    counter = n_variable_x - gsl_vector_int_get(fixed_variables, 0);

    struct global_solution *solution = managed::create<global_solution>();
    struct local_solution *local_amc = managed::create<local_solution>();
    struct dataset *dataset_min = managed::create<dataset>();
    struct dataset *dataset_max = managed::create<dataset>();

    // struct global_solution *solution = (struct global_solution *)mem_alloc(sizeof(struct global_solution), "global_solution", "cpm", 2, FALSE);
    // struct local_solution *local_amc = (struct local_solution *)mem_alloc(sizeof(struct local_solution), "local_solution", "cpm", 3, FALSE);
    // struct dataset *dataset_min = (struct dataset *)mem_alloc(sizeof(struct dataset), "dataset_min", "cpm", 4, FALSE);
    // struct dataset *dataset_max = (struct dataset *)mem_alloc(sizeof(struct dataset), "dataset_max", "cpm", 4, FALSE);

    if (solution == NULL || local_amc == NULL || dataset_min == NULL || dataset_max == NULL)
    {
        free_all_fixed_variables(lp_Y, instance, temp_v, temp_coef, effective_constraints, fixed_variables);
        return (OUT_OF_MEMORY);
    }

    local_solution_constructor(local_amc, counter, n_variable_y);
    global_solution_constructor(solution, counter, n_variable_y);

    dataset_constructor(dataset_min, noc, counter, n_variable_y, CUTS_MAX);
    dataset_constructor(dataset_max, noc, counter, n_variable_y, CUTS_MAX);

    construct_parameters(instance, effective_constraints, fixed_variables, dataset_min, num_alt_i);
    noc = gsl_vector_int_get(effective_constraints, 0);

    /* construct the dataset for the maximizing problem */
    dataset_copy(dataset_max, dataset_min); // make a copy for the maximizing problem

    gsl_matrix_scale(dataset_max->Q, -1);
    gsl_vector_scale(dataset_max->dy, -1);

    /* establish the LP problem for X */
    const int nov_x = dataset_min->xl->size; // reset the # of variables in X

    auto make_probability_start = [&](const double *source) -> gsl_vector *
    {
        if (!source) return nullptr;
        gsl_vector *result = managed::vector_alloc(n_variable_y);
        int z = 0;
        for (int k = Ai_begin; k <= Ai_end; ++k)
            gsl_vector_set(result, z++, source[k]);
        for (int k = Aj_begin; k <= Aj_end; ++k)
            gsl_vector_set(result, z++, source[k]);
        return result;
    };
    gsl_vector *min_start = make_probability_start(P_min_start);
    gsl_vector *max_start = make_probability_start(P_max_start);

    LpModel *lp_X_std_min = make_lp(noc + nov_x, nov_x);
    LpModel *lp_X_std_max = make_lp(noc + nov_x, nov_x);

    if (lp_X_std_min == NULL || lp_X_std_max == NULL)
    {
        free_all_fixed_variables(lp_Y, instance, temp_v, temp_coef, effective_constraints, fixed_variables);
        global_solution_destructor(solution);
        local_solution_destructor(local_amc);
        dataset_destructor(dataset_min);
        dataset_destructor(dataset_max);

        return (OUT_OF_MEMORY);
    }

    set_verbose(lp_X_std_min, kLpNeutral);
    set_minim(lp_X_std_min);

    construct_standard_lp(dataset_min, lp_X_std_min, noc);

    gsl_vector *coef_x = managed::vector_calloc(nov_x + 1);

    gsl_vector_set_all(coef_x, 0);
    set_obj_fn(lp_X_std_min, coef_x->data);

    return_value = std::getenv("DECIDEIT_SKIP_INITIAL_X_SOLVE")
        ? OPTIMAL : solve(lp_X_std_min);

    if (return_value != OPTIMAL)
    {
        free_all_fixed_variables(lp_Y, instance, temp_v, temp_coef, effective_constraints, fixed_variables);
        global_solution_destructor(solution);
        local_solution_destructor(local_amc);
        dataset_destructor(dataset_min);
        dataset_destructor(dataset_max);
        managed::vector_free(coef_x);
        free_lp(&lp_X_std_min);
        free_lp(&lp_X_std_max);

        return return_value;
    }

    /* if the code runs here, lp_X_std_max should be feasible too */
    set_verbose(lp_X_std_max, kLpNeutral);
    set_minim(lp_X_std_max);

    construct_standard_lp(dataset_max, lp_X_std_max, noc);

    gsl_vector_int *unfixed_indices = managed::vector_int_calloc(1 + counter);

    gsl_vector_int_set(unfixed_indices, 0, counter);
    j = 1;

    for (i = 1; i <= n_variable_x; i++)
    {
        if (gsl_vector_int_get(fixed_variables, i) == 0)
        {
            if (i <= num_alt_i)
            {
                gsl_vector_int_set(unfixed_indices, j, i);
            }
            else
            {
                gsl_vector_int_set(unfixed_indices, j, -i);
            }

            j++;
        }
    }

    /* The two objective senses share only immutable model data.  On machines
       with at least two cores they can therefore be solved concurrently.  It
       is the production elapsed-time default; --sequential-directions remains
       available for single-core algorithm comparisons. */
    const bool parallel_directions =
        std::getenv("DECIDEIT_PARALLEL_DIRECTIONS") != nullptr;
    int min_exit = OK;
    int max_exit = OK;
    struct local_solution *local_amc_max = nullptr;
    struct global_solution *solution_max = nullptr;
    LpModel *lp_Y_max = nullptr;
    int noc_max = gsl_vector_int_get(effective_constraints, 0);
    RltCrossLevelState cross_direction_rlt_state;
    const bool cross_direction_rlt_basis =
        std::getenv("DECIDEIT_CROSS_DIRECTION_RLT_BASIS") != nullptr;
    RltCrossLevelState *effective_min_rlt_state = cross_direction_rlt_basis
        ? &cross_direction_rlt_state : min_rlt_state;
    RltCrossLevelState *effective_max_rlt_state = cross_direction_rlt_basis
        ? &cross_direction_rlt_state : max_rlt_state;

    P_lobo[0] = 0;
    P_upbo[0] = 0;
    if (parallel_directions)
    {
        local_amc_max = managed::create<local_solution>();
        solution_max = managed::create<global_solution>();
        local_solution_constructor(local_amc_max, counter, n_variable_y);
        global_solution_constructor(solution_max, counter, n_variable_y);
        lp_Y_max = construct_Y(n_variable_y, mY, LC_P, num_alt_i, num_alt_j,
                               Ai_begin, Ai_end, Aj_begin, Aj_end,
                               P_lobo, P_upbo);

        struct DirectionResult
        {
            int status;
            DecideItDirectionStats stats;
        };
        std::future<DirectionResult> max_future = std::async(
            std::launch::async,
            [&]() -> DirectionResult
            {
                const int status = cpm_controller(
                    instance, local_amc_max, lp_X_std_max, lp_Y_max,
                    dataset_max, &noc_max,
                    P_lobo, P_upbo, V_lobo, V_upbo,
                    solution_max, CUTS_MAX,
                    P_max_point, V_max_point,
                    Ai_begin, Ai_end, Aj_begin, Aj_end,
                    fixed_variables,
                    LC_V, LC_P, P_rows, *mX, *mY, &LinParaMax,
                    FALSE, max_start, effective_max_rlt_state);
                return DirectionResult{status, last_decideit_stats().maximum};
            });

        min_exit = cpm_controller(
            instance, local_amc, lp_X_std_min, lp_Y,
            dataset_min, &noc,
            P_lobo, P_upbo, V_lobo, V_upbo,
            solution, CUTS_MAX,
            P_min_point, V_min_point,
            Ai_begin, Ai_end, Aj_begin, Aj_end,
            fixed_variables,
            LC_V, LC_P, P_rows, *mX, *mY, &LinPara,
            TRUE, min_start, effective_min_rlt_state);
        const DecideItDirectionStats min_stats = last_decideit_stats().minimum;
        const DirectionResult max_result = max_future.get();
        max_exit = max_result.status;
        DecideItRunStats combined_stats;
        combined_stats.minimum = min_stats;
        combined_stats.maximum = max_result.stats;
        set_decideit_stats(combined_stats);
    }
    else
    {
        min_exit = cpm_controller(
            instance, local_amc, lp_X_std_min, lp_Y,
            dataset_min, &noc,
            P_lobo, P_upbo, V_lobo, V_upbo,
            solution, CUTS_MAX,
            P_min_point, V_min_point,
            Ai_begin, Ai_end, Aj_begin, Aj_end,
            fixed_variables,
            LC_V, LC_P, P_rows, *mX, *mY, &LinPara,
            TRUE, min_start, effective_min_rlt_state);

        if (min_exit == OK)
        {
            local_amc->cbov = 1e7;
            solution->alpha = 1e7;
            noc = gsl_vector_int_get(effective_constraints, 0);
            max_exit = cpm_controller(
                instance, local_amc, lp_X_std_max, lp_Y,
                dataset_max, &noc,
                P_lobo, P_upbo, V_lobo, V_upbo,
                solution, CUTS_MAX,
                P_max_point, V_max_point,
                Ai_begin, Ai_end, Aj_begin, Aj_end,
                fixed_variables,
                LC_V, LC_P, P_rows, *mX, *mY, &LinPara,
                FALSE, max_start, effective_max_rlt_state);
        }
    }

    if (parallel_directions)
    {
        local_solution_destructor(local_amc_max);
        global_solution_destructor(solution_max);
        free_lp(&lp_Y_max);
        destroy_LinearizationParameters(&LinParaMax);
    }

    free_local_variables(local_amc, lp_Y, instance, solution, temp_v, temp_coef, effective_constraints, fixed_variables, dataset_min, dataset_max, unfixed_indices, lp_X_std_min, lp_X_std_max, coef_x);
    destroy_LinearizationParameters(&LinPara);
    managed::vector_free(min_start);
    managed::vector_free(max_start);

    exitFlag = min_exit != OK ? min_exit : max_exit;
    return exitFlag;
}

#if 0
int cpm_solver(
    d_row pl,d_row pu,link_set pc,d_row vl,d_row vu,link_set vc,
    d_row pmin,d_row vmin,d_row pmax,d_row vmax,
    int ai,int ae,int bi,int be,int cuts,
    std::vector<LinCons>& lv,std::vector<LinCons>& lp,int* mx,int* my)
{
    managed::Scope query_resources;
    if(ai<1 || ae<ai || bi<=ae || be<bi || be>VX_MAXDIM || cuts<1 || cuts>10000)
        throw std::invalid_argument("Invalid pair dimensions or cut limit");
    if(pc[0].head<0||pc[0].head>5*VX_MAXDIM||vc[0].head<0||vc[0].head>5*VX_MAXDIM)
        throw std::invalid_argument("Too many statements");
    auto included=[&](int i){return (i>=ai&&i<=ae)||(i>=bi&&i<=be);};
    for(int i=ai;i<=be;i++)if(included(i)) {
        if(!std::isfinite(pl[i])||!std::isfinite(pu[i])||!std::isfinite(vl[i])||!std::isfinite(vu[i])||
           pl[i]<0||pu[i]>1||vl[i]<0||vu[i]>1||pl[i]>pu[i]||vl[i]>vu[i])
            throw std::invalid_argument("Invalid normalized probability/value bounds");
    }
    int result=cpm_solver_impl(pl,pu,pc,vl,vu,vc,pmin,vmin,pmax,vmax,ai,ae,bi,be,cuts,lv,lp,mx,my);
    if(result!=OK) return result;
    auto verify=[&](const double* p,const double* v,double reported){
        double objective=0,sum_a=0,sum_b=0;
        for(int i=ai;i<=be;i++)if(included(i)){
            if(!std::isfinite(p[i])||!std::isfinite(v[i])||p[i]<pl[i]-1e-7||p[i]>pu[i]+1e-7||v[i]<vl[i]-1e-7||v[i]>vu[i]+1e-7)
                throw std::runtime_error("Returned solution violates an original bound");
            objective+=(i<=ae?1.:-1.)*p[i]*v[i];
            if(i<=ae)sum_a+=p[i];else sum_b+=p[i];
        }
        if(std::fabs(sum_a-1)>1e-7||std::fabs(sum_b-1)>1e-7)throw std::runtime_error("Returned probabilities do not sum to one");
        for(int base=0;base<2;base++){
            const link_rec* links=base?vc:pc;const double* value=base?v:p;
            for(int k=1;k<=links[0].head;k++){
                int i=links[k].head,j=links[k].tail;
                if(included(i)&&included(j)){
                    double diff=value[i]-value[j];
                    if((links[k].cmp=='<'&&diff>links[k].value+1e-7)||(links[k].cmp=='>'&&diff<links[k].value-1e-7))
                        throw std::runtime_error("Returned solution violates an original difference statement");
                }
            }
        }
        if(!std::isfinite(reported)||std::fabs(objective-reported)>DBLP_EPSILON)
            throw std::runtime_error("Reported objective disagrees with returned solution");
    };
    verify(pmin,vmin,pmin[0]);verify(pmax,vmax,-pmax[0]);
    return OK;
}
#endif

namespace
{
bool selected_index(int index, int ai, int ae, int bi, int be)
{
    return (index >= ai && index <= ae) || (index >= bi && index <= be);
}

int selected_local_index(int index, int ai, int ae, int bi)
{
    return index <= ae ? index - ai + 1 : ae - ai + 2 + index - bi;
}

std::vector<LinCons> localize_rows(const std::vector<LinCons> &rows,
                                   int ai, int ae, int bi, int be)
{
    std::vector<LinCons> result;
    for (const LinCons &source : rows)
    {
        bool included = true;
        for (int index : source.index)
            included = included && selected_index(index, ai, ae, bi, be);
        if (!included) continue;
        LinCons row = source;
        for (int &index : row.index)
            index = selected_local_index(index, ai, ae, bi);
        result.push_back(row);
    }
    return result;
}

void validate_rows(const std::vector<LinCons> &rows, int limit)
{
    if (rows.size() > static_cast<size_t>(5 * VX_MAXDIM))
        throw std::invalid_argument("Too many linear statements");
    for (const LinCons &row : rows)
    {
        if (row.index.empty() || row.index.size() != row.coefficient.size() ||
            !std::isfinite(row.b))
            throw std::invalid_argument("Invalid sparse linear statement");
        std::vector<bool> seen(limit + 1, false);
        for (size_t k = 0; k < row.index.size(); ++k)
        {
            if (row.index[k] < 1 || row.index[k] > limit || seen[row.index[k]] ||
                !std::isfinite(row.coefficient[k]) || row.coefficient[k] == 0.0)
                throw std::invalid_argument("Invalid sparse linear-statement term");
            seen[row.index[k]] = true;
        }
    }
}
} // namespace

int cpm_solver_rows(
    d_row pl, d_row pu, const std::vector<LinCons> &probability_rows,
    d_row vl, d_row vu, const std::vector<LinCons> &value_rows,
    d_row pmin, d_row vmin, d_row pmax, d_row vmax,
    int ai, int ae, int bi, int be, int cuts,
    std::vector<LinCons> &lc_v, std::vector<LinCons> &lc_p, int *mx, int *my,
    const double *pmin_start, const double *pmax_start,
    RltCrossLevelState *min_rlt_state,
    RltCrossLevelState *max_rlt_state)
{
    managed::Scope query_resources;
    reset_decideit_stats();
    reset_row_reduction_stats();
    if (ai < 1 || ae < ai || bi <= ae || be < bi || be > VX_MAXDIM ||
        cuts < 1 || cuts > 10000)
        throw std::invalid_argument("Invalid pair dimensions or cut limit");
    validate_rows(probability_rows, be);
    validate_rows(value_rows, be);
    for (int i = ai; i <= be; ++i)
    {
        if (!selected_index(i, ai, ae, bi, be)) continue;
        if (!std::isfinite(pl[i]) || !std::isfinite(pu[i]) ||
            !std::isfinite(vl[i]) || !std::isfinite(vu[i]) ||
            pl[i] < 0 || pu[i] > 1 || vl[i] < 0 || vu[i] > 1 ||
            pl[i] > pu[i] || vl[i] > vu[i])
            throw std::invalid_argument("Invalid normalized probability/value bounds");
    }

    const std::vector<LinCons> local_p = localize_rows(probability_rows, ai, ae, bi, be);
    const std::vector<LinCons> local_v = localize_rows(value_rows, ai, ae, bi, be);
    const int result = cpm_solver_impl(
        pl, pu, local_p, vl, vu, local_v, pmin, vmin, pmax, vmax,
        ai, ae, bi, be, cuts, lc_v, lc_p, mx, my,
        pmin_start, pmax_start, min_rlt_state, max_rlt_state);
    if (result != OK) return result;

    auto verify = [&](const char *sense, const double *p, const double *v,
                      double reported)
    {
        double objective = 0.0, sum_a = 0.0, sum_b = 0.0;
        for (int i = ai; i <= be; ++i)
        {
            if (!selected_index(i, ai, ae, bi, be)) continue;
            if (!std::isfinite(p[i]) || !std::isfinite(v[i]) ||
                p[i] < pl[i] - 1e-7 || p[i] > pu[i] + 1e-7 ||
                v[i] < vl[i] - 1e-7 || v[i] > vu[i] + 1e-7)
                throw std::runtime_error("Returned solution violates an original bound");
            objective += (i <= ae ? 1.0 : -1.0) * p[i] * v[i];
            (i <= ae ? sum_a : sum_b) += p[i];
        }
        if (std::fabs(sum_a - 1.0) > 1e-7 || std::fabs(sum_b - 1.0) > 1e-7)
            throw std::runtime_error("Returned probabilities do not sum to one");
        for (int base = 0; base < 2; ++base)
        {
            const std::vector<LinCons> &rows = base ? value_rows : probability_rows;
            const double *point = base ? v : p;
            for (const LinCons &row : rows)
            {
                bool included = true;
                double lhs = 0.0;
                for (size_t k = 0; k < row.index.size(); ++k)
                {
                    included = included && selected_index(row.index[k], ai, ae, bi, be);
                    lhs += row.coefficient[k] * point[row.index[k]];
                }
                if (included && lhs > row.b + 1e-7)
                    throw std::runtime_error("Returned solution violates an original linear statement");
            }
        }
        if (!std::isfinite(reported) || std::fabs(objective - reported) > DBLP_EPSILON)
            throw std::runtime_error("Reported objective disagrees with returned solution");
        if (std::getenv("DECIDEIT_TRACE_POINT"))
        {
            std::fprintf(stderr,
                         "TRACE {\"event\":\"solution_point\",\"sense\":\"%s\",\"objective\":%.17g,\"p\":[",
                         sense, objective);
            bool first = true;
            for (int i = ai; i <= be; ++i)
            {
                if (!selected_index(i, ai, ae, bi, be)) continue;
                std::fprintf(stderr, "%s%.17g", first ? "" : ",", p[i]);
                first = false;
            }
            std::fprintf(stderr, "],\"v\":[");
            first = true;
            for (int i = ai; i <= be; ++i)
            {
                if (!selected_index(i, ai, ae, bi, be)) continue;
                std::fprintf(stderr, "%s%.17g", first ? "" : ",", v[i]);
                first = false;
            }
            std::fprintf(stderr, "]}\n");
        }
    };
    verify("min", pmin, vmin, pmin[0]);
    verify("max", pmax, vmax, -pmax[0]);
    return OK;
}

int cpm_solver(
    d_row pl, d_row pu, link_set pc, d_row vl, d_row vu, link_set vc,
    d_row pmin, d_row vmin, d_row pmax, d_row vmax,
    int ai, int ae, int bi, int be, int cuts,
    std::vector<LinCons> &lc_v, std::vector<LinCons> &lc_p, int *mx, int *my)
{
    if (pc[0].head < 0 || pc[0].head > 5 * VX_MAXDIM ||
        vc[0].head < 0 || vc[0].head > 5 * VX_MAXDIM)
        throw std::invalid_argument("Too many statements");
    std::vector<LinCons> probability_rows, value_rows;
    for (int k = 1; k <= pc[0].head; ++k)
        probability_rows.push_back(pc[k].cmp == '<'
            ? difference_row(pc[k].head, pc[k].tail, pc[k].value)
            : difference_row(pc[k].tail, pc[k].head, -pc[k].value));
    for (int k = 1; k <= vc[0].head; ++k)
        value_rows.push_back(vc[k].cmp == '<'
            ? difference_row(vc[k].head, vc[k].tail, vc[k].value)
            : difference_row(vc[k].tail, vc[k].head, -vc[k].value));
    return cpm_solver_rows(pl, pu, probability_rows, vl, vu, value_rows,
                           pmin, vmin, pmax, vmax, ai, ae, bi, be, cuts,
                           lc_v, lc_p, mx, my);
}
