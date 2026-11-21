#include "dblp.h"

void show_LC(LinCons *LC, int num)
{
    printf("\n");
    for (int i = 0; i < num; i++)
    {
        printf("%d:", i);
        for (size_t k = 0; k < LC[i].index.size(); ++k)
            printf(" %+.12g*z%d", LC[i].coefficient[k], LC[i].index[k]);
        printf(" <= %.12g\n", LC[i].b);
    }
}

#if 0
void DataPreparation(
    link_set V_links,
    gsl_vector_int *effective_constraints,
    link_set P_links,
    d_row P_lobo,
    d_row P_upbo,
    int Ai_begin,
    int Ai_end,
    int Aj_begin,
    int Aj_end,
    std::vector<LinCons> &LC_V,
    std::vector<LinCons> &LC_P,
    int *mX,
    int *mY)
{
    // set lower and upper bounds in V-base
    int i, j, head, tail, counter, pointer;
    const int num_alt_i = Ai_end - Ai_begin + 1, num_alt_j = Aj_end - Aj_begin + 1;
    const int n_variable = num_alt_i + num_alt_j;

    if (*mX == -1)
    {
        counter = 1;
        pointer = 0;

        for (i = 1; i <= V_links[0].head; i++)
        {
            if (gsl_vector_int_get(effective_constraints, counter) + 1 == i)
            {
                switch (V_links[i].cmp)
                {
                case '<':
                    LC_V[pointer].j1 = V_links[i].head;
                    LC_V[pointer].j2 = V_links[i].tail;
                    LC_V[pointer].b = V_links[i].value;

                    break;

                case '>':
                    LC_V[pointer].j1 = V_links[i].tail;
                    LC_V[pointer].j2 = V_links[i].head;
                    LC_V[pointer].b = -V_links[i].value;

                    break;

                default:
                    break;
                }

                pointer++;
                counter++;
            }
        }
    }
    else
    {
        if (*mX > 0)
        {
            j = 0;
            pointer = 1;

            for (i = 0; i < *mX; i++)
            {
                if (gsl_vector_int_get(effective_constraints, pointer) == i)
                {
                    LC_V[j].j1 = LC_V[i].j1;
                    LC_V[j].j2 = LC_V[i].j2;
                    LC_V[j].b = LC_V[i].b;

                    pointer++;
                    j++;
                }
            }
        }
    }

    *mX = gsl_vector_int_get(effective_constraints, 0);

    LpModel *lp = NULL;

    if (*mY == -1)
    {
        lp = make_lp(P_links[0].head + 2, n_variable); // two normalizations
    }
    else
    {
        lp = make_lp(*mY + 2, n_variable); // two normalizations
    }

    /* set up the original LP */
    set_verbose(lp, kLpNeutral);
    set_maxim(lp);

    /* set lower and upper bound for lp */
    j = 1;

    for (i = Ai_begin; i <= Ai_end; i++)
    {
        set_lowbo(lp, j, P_lobo[i]);
        set_upbo(lp, j, P_upbo[i]);

        j++;
    }

    for (i = Aj_begin; i <= Aj_end; i++)
    {
        set_lowbo(lp, j, P_lobo[i]);
        set_upbo(lp, j, P_upbo[i]);

        j++;
    }

    /* set linear constraints for lp */
    int colno[2], row_no = 1;
    REAL row[2];

    if (*mY == -1)
    {
        for (i = 1; i <= P_links[0].head; i++)
        {
            switch (P_links[i].cmp)
            {
            case '<':
                set_rh(lp, i, P_links[i].value);
                row[0] = 1;
                row[1] = -1;

                break;

            case '>':
                set_rh(lp, i, -P_links[i].value);
                row[0] = -1;
                row[1] = 1;

                break;

            default:
                break;
            }

            head = P_links[i].head;
            tail = P_links[i].tail;

            if ((head >= Ai_begin) && (head <= Ai_end))
            {
                colno[0] = head - Ai_begin + 1;
            }
            else
            {
                colno[0] = head - Aj_begin + num_alt_i + 1;
            }

            if ((tail >= Ai_begin) && (tail <= Ai_end))
            {
                colno[1] = tail - Ai_begin + 1;
            }
            else
            {
                colno[1] = tail - Aj_begin + num_alt_i + 1;
            }

            set_rowex(lp, i, 2, row, colno); // the # of nonzero elements is always 2 for comparaive statements
            set_constr_type(lp, i, kLpLe);
        }

        row_no = P_links[0].head + 1;
    }
    else
    {
        if (*mY > 0)
        {
            for (i = 0; i < *mY; i++)
            {
                set_rh(lp, i + 1, LC_P[i].b);

                row[0] = 1;
                row[1] = -1;
                colno[0] = LC_P[i].j1;
                colno[1] = LC_P[i].j2;

                set_rowex(lp, i + 1, 2, row, colno); // the # of nonzero elements is always 2 for comparaive statements
                set_constr_type(lp, i + 1, kLpLe);
            }

            row_no = *mY + 1;
        }
    }

    gsl_vector *normal_vector = managed::vector_calloc(n_variable + 1);

    gsl_vector_view temp = gsl_vector_subvector(normal_vector, 1, num_alt_i);
    gsl_vector_set_all(&temp.vector, 1);
    set_row(lp, row_no, normal_vector->data);
    set_constr_type(lp, row_no, kLpEq);
    set_rh(lp, row_no, 1);
    row_no++;

    gsl_vector_set_zero(normal_vector);

    temp = gsl_vector_subvector(normal_vector, num_alt_i + 1, num_alt_j);
    gsl_vector_set_all(&temp.vector, 1);
    set_row(lp, row_no, normal_vector->data);
    set_constr_type(lp, row_no, kLpEq);
    set_rh(lp, row_no, 1);

    managed::vector_free(normal_vector);

    // print_lp(lp);

    /* check Ax <= b to figure out the indices of effective constraints */
    counter = 0;

    if (*mY == -1)
    {
        for (i = 1; i <= P_links[0].head; i++)
        {
            /* solve the LP for the ith constraint */
            switch (P_links[i].cmp)
            {
            case '<':
                row[0] = 1;
                row[1] = -1;

                break;

            case '>':
                row[0] = -1;
                row[1] = 1;

                break;

            default:
                break;
            }

            head = P_links[i].head;
            tail = P_links[i].tail;

            if ((head >= Ai_begin) && (head <= Ai_end))
            {
                colno[0] = head - Ai_begin + 1;
            }
            else
            {
                colno[0] = head - Aj_begin + num_alt_i + 1;
            }

            if ((tail >= Ai_begin) && (tail <= Ai_end))
            {
                colno[1] = tail - Ai_begin + 1;
            }
            else
            {
                colno[1] = tail - Aj_begin + num_alt_i + 1;
            }

            set_obj_fnex(lp, 2, row, colno);
            // print_lp(lp);
            solve(lp);

            if (get_objective(lp) >= get_rh(lp, i) - DBLP_EPSILON) // a valid constraint
            {
                switch (P_links[i].cmp)
                {
                case '<':
                    LC_P[counter].j1 = P_links[i].head;
                    LC_P[counter].j2 = P_links[i].tail;
                    LC_P[counter].b = P_links[i].value;

                    break;

                case '>':
                    LC_P[counter].j1 = P_links[i].tail;
                    LC_P[counter].j2 = P_links[i].head;
                    LC_P[counter].b = -P_links[i].value;

                    break;

                default:
                    break;
                }

                counter++;
            }
        }
    }
    else
    {
        if (*mY > 0)
        {
            for (i = 0; i < *mY; i++)
            {
                /* solve the LP for the ith constraint */
                row[0] = 1;
                row[1] = -1;
                colno[0] = LC_P[i].j1;
                colno[1] = LC_P[i].j2;

                set_obj_fnex(lp, 2, row, colno);
                // print_lp(lp);
                solve(lp);

                if (get_objective(lp) >= get_rh(lp, i + 1) - DBLP_EPSILON) // a valid constraint
                {
                    LC_P[counter].j1 = LC_P[i].j1;
                    LC_P[counter].j2 = LC_P[i].j2;
                    LC_P[counter].b = LC_P[i].b;

                    counter++;
                }
            }
        }
    }

    *mY = counter;

    delete_lp(lp);
}
#endif

namespace
{
thread_local RowReductionStats row_reduction_stats;

enum class BaseRelaxation
{
    Box,
    TwoSimplexes
};

void set_sparse_row(LpModel *lp, int row_number, const LinCons &row)
{
    std::vector<REAL> coefficient(row.coefficient.begin(), row.coefficient.end());
    std::vector<int> column(row.index.begin(), row.index.end());
    set_rowex(lp, row_number, static_cast<int>(column.size()),
              coefficient.data(), column.data());
    set_constr_type(lp, row_number, kLpLe);
    set_rh(lp, row_number, row.b);
}

void set_sparse_objective(LpModel *lp, const LinCons &row)
{
    std::vector<REAL> coefficient(row.coefficient.begin(), row.coefficient.end());
    std::vector<int> column(row.index.begin(), row.index.end());
    set_obj_fnex(lp, static_cast<int>(column.size()), coefficient.data(), column.data());
}

enum class RowFilterMode
{
    Legacy,
    Hybrid,
    Cheap
};

RowFilterMode row_filter_mode()
{
    const char *value = std::getenv("DECIDEIT_ROW_FILTER");
    const std::string selected = value ? value : "cheap";
    if (selected == "legacy") return RowFilterMode::Legacy;
    if (selected == "hybrid") return RowFilterMode::Hybrid;
    if (selected == "cheap") return RowFilterMode::Cheap;
    throw std::invalid_argument("row filter must be cheap, hybrid, or legacy");
}

std::vector<std::vector<double>> dense_rows(const std::vector<LinCons> &rows,
                                            int n_variable)
{
    std::vector<std::vector<double>> result(
        rows.size(), std::vector<double>(n_variable, 0.0));
    for (size_t i = 0; i < rows.size(); ++i)
        for (size_t k = 0; k < rows[i].index.size(); ++k)
            result[i].at(rows[i].index[k] - 1) += rows[i].coefficient[k];
    return result;
}

double dot_product(const std::vector<double> &a,
                   const std::vector<double> &x)
{
    long double result = 0.0L;
    for (size_t j = 0; j < a.size(); ++j)
        result += static_cast<long double>(a[j]) * x[j];
    return static_cast<double>(result);
}

double row_guard(const LinCons &row,
                 const std::vector<double> &lower,
                 const std::vector<double> &upper)
{
    long double scale = 1.0L + std::fabs(row.b);
    for (size_t k = 0; k < row.index.size(); ++k)
    {
        const int j = row.index[k] - 1;
        const double magnitude = std::fmax(std::fabs(lower[j]),
                                           std::fabs(upper[j]));
        scale += std::fabs(row.coefficient[k]) * magnitude;
    }
    return static_cast<double>(1.0e-10L * scale);
}

bool same_left_hand_side(const std::vector<double> &left,
                         const std::vector<double> &right)
{
    return left == right;
}

bool maximize_simplex_group(const std::vector<double> &coefficient,
                            const std::vector<double> &lower,
                            const std::vector<double> &upper,
                            int begin, int end,
                            std::vector<double> &point)
{
    double remaining = 1.0;
    std::vector<int> order;
    order.reserve(end - begin);
    for (int j = begin; j < end; ++j)
    {
        point[j] = lower[j];
        remaining -= lower[j];
        order.push_back(j);
    }
    if (remaining < -1.0e-10) return false;
    std::sort(order.begin(), order.end(), [&](int left, int right)
    {
        if (coefficient[left] != coefficient[right])
            return coefficient[left] > coefficient[right];
        return left < right;
    });
    for (int j : order)
    {
        const double addition = std::fmin(remaining, upper[j] - lower[j]);
        if (addition > 0.0)
        {
            point[j] += addition;
            remaining -= addition;
        }
    }
    return remaining <= 1.0e-10;
}

bool relaxed_maximizer(const std::vector<double> &coefficient,
                       const std::vector<double> &lower,
                       const std::vector<double> &upper,
                       BaseRelaxation relaxation,
                       int simplex_split,
                       std::vector<double> &point)
{
    point.resize(coefficient.size());
    if (relaxation == BaseRelaxation::Box)
    {
        for (size_t j = 0; j < coefficient.size(); ++j)
            point[j] = coefficient[j] > 0.0 ? upper[j] :
                coefficient[j] < 0.0 ? lower[j] : 0.5 * (lower[j] + upper[j]);
        return true;
    }
    return maximize_simplex_group(coefficient, lower, upper, 0,
                                  simplex_split, point) &&
        maximize_simplex_group(coefficient, lower, upper, simplex_split,
                               static_cast<int>(coefficient.size()), point);
}

bool feasible_witness(const std::vector<double> &point,
                      const std::vector<std::vector<double>> &coefficient,
                      const std::vector<LinCons> &rows,
                      const std::vector<double> &lower,
                      const std::vector<double> &upper,
                      BaseRelaxation relaxation,
                      int simplex_split)
{
    for (size_t j = 0; j < point.size(); ++j)
        if (point[j] < lower[j] - 1.0e-10 || point[j] > upper[j] + 1.0e-10)
            return false;
    if (relaxation == BaseRelaxation::TwoSimplexes)
    {
        double first = 0.0, second = 0.0;
        for (int j = 0; j < simplex_split; ++j) first += point[j];
        for (size_t j = simplex_split; j < point.size(); ++j) second += point[j];
        if (std::fabs(first - 1.0) > 1.0e-9 ||
            std::fabs(second - 1.0) > 1.0e-9)
            return false;
    }
    for (size_t i = 0; i < rows.size(); ++i)
        if (dot_product(coefficient[i], point) >
            rows[i].b + row_guard(rows[i], lower, upper))
            return false;
    return true;
}

std::vector<int> filter_effective_rows(
    const std::vector<LinCons> &rows,
    const std::vector<double> &lower,
    const std::vector<double> &upper,
    BaseRelaxation relaxation,
    int simplex_split,
    const char *failure_message)
{
    const auto started = std::chrono::steady_clock::now();
    const int n_variable = static_cast<int>(lower.size());
    const int n_rows = static_cast<int>(rows.size());
    row_reduction_stats.input_rows += n_rows;
    const std::vector<std::vector<double>> coefficient =
        dense_rows(rows, n_variable);
    std::vector<int> state(n_rows, -1); // -1 unknown, 0 remove, 1 retain

    const RowFilterMode mode = row_filter_mode();
    if (mode != RowFilterMode::Legacy)
    {
        // Exact duplicate/dominance check.  Keeping the smallest RHS preserves
        // the feasible set and avoids generating duplicate RLT products.
        for (int i = 0; i < n_rows; ++i)
        {
            for (int j = 0; j < n_rows; ++j)
            {
                if (i == j || !same_left_hand_side(coefficient[i], coefficient[j]))
                    continue;
                if (rows[j].b < rows[i].b ||
                    (rows[j].b == rows[i].b && j < i))
                {
                    state[i] = 0;
                    row_reduction_stats.cheap_removed++;
                    break;
                }
            }
        }

        // Optimize over only the box, or over the two bounded probability
        // simplexes.  A low upper bound proves redundancy.  If the relaxed
        // maximizer also satisfies all general rows, it proves effectiveness.
        std::vector<double> point;
        for (int i = 0; i < n_rows; ++i)
        {
            if (state[i] != -1) continue;
            if (!relaxed_maximizer(coefficient[i], lower, upper, relaxation,
                                   simplex_split, point))
                throw std::runtime_error(failure_message);
            const double value = dot_product(coefficient[i], point);
            const double guard = row_guard(rows[i], lower, upper);
            if (value < rows[i].b - DBLP_EPSILON - guard)
            {
                state[i] = 0;
                row_reduction_stats.cheap_removed++;
            }
            else if (value >= rows[i].b - DBLP_EPSILON + guard &&
                     feasible_witness(point, coefficient, rows, lower, upper,
                                      relaxation, simplex_split))
            {
                state[i] = 1;
                row_reduction_stats.cheap_kept++;
            }
        }
    }

    if (mode == RowFilterMode::Cheap)
    {
        for (int &value : state)
            if (value == -1)
            {
                value = 1;
                row_reduction_stats.conservative_kept++;
            }
    }

    int unresolved = 0;
    for (int value : state) unresolved += value == -1;
    if (unresolved > 0)
    {
        const int equality_rows =
            relaxation == BaseRelaxation::TwoSimplexes ? 2 : 0;
        LpModel *lp = make_lp(n_rows + equality_rows, n_variable);
        set_verbose(lp, kLpNeutral);
        set_maxim(lp);
        for (int j = 0; j < n_variable; ++j)
        {
            set_lowbo(lp, j + 1, lower[j]);
            set_upbo(lp, j + 1, upper[j]);
        }
        for (int i = 0; i < n_rows; ++i)
            set_sparse_row(lp, i + 1, rows[i]);
        if (equality_rows)
        {
            std::vector<REAL> normal(n_variable + 1, 0.0);
            for (int j = 1; j <= simplex_split; ++j) normal[j] = 1.0;
            set_row(lp, n_rows + 1, normal.data());
            set_constr_type(lp, n_rows + 1, kLpEq);
            set_rh(lp, n_rows + 1, 1.0);
            std::fill(normal.begin(), normal.end(), 0.0);
            for (int j = simplex_split + 1; j <= n_variable; ++j) normal[j] = 1.0;
            set_row(lp, n_rows + 2, normal.data());
            set_constr_type(lp, n_rows + 2, kLpEq);
            set_rh(lp, n_rows + 2, 1.0);
        }

        std::vector<double> solution(n_variable, 0.0);
        for (int i = 0; i < n_rows; ++i)
        {
            if (state[i] != -1) continue;
            set_sparse_objective(lp, rows[i]);
            if (solve(lp) != OPTIMAL)
                throw std::runtime_error(failure_message);
            row_reduction_stats.lp_calls++;
            row_reduction_stats.simplex_iterations +=
                get_simplex_iterations(lp);
            state[i] = get_objective(lp) >= rows[i].b - DBLP_EPSILON ? 1 : 0;

            if (mode != RowFilterMode::Hybrid) continue;
            get_variables(lp, solution.data());
            if (!feasible_witness(solution, coefficient, rows, lower, upper,
                                  relaxation, simplex_split))
                continue;
            // One optimum is often on several facets.  Reuse it as a rigorous
            // primal witness for all remaining rows before changing objective.
            for (int j = 0; j < n_rows; ++j)
            {
                if (state[j] != -1) continue;
                const double guard = row_guard(rows[j], lower, upper);
                if (dot_product(coefficient[j], solution) >=
                    rows[j].b - DBLP_EPSILON + guard)
                {
                    state[j] = 1;
                    row_reduction_stats.witness_kept++;
                }
            }
        }
        delete_lp(lp);
    }

    std::vector<int> retained;
    retained.reserve(n_rows);
    for (int i = 0; i < n_rows; ++i)
        if (state[i] == 1) retained.push_back(i);
    row_reduction_stats.retained_rows += static_cast<int>(retained.size());
    row_reduction_stats.wall_s += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return retained;
}
} // namespace

void reset_row_reduction_stats()
{
    row_reduction_stats = RowReductionStats{};
}

const RowReductionStats &last_row_reduction_stats()
{
    return row_reduction_stats;
}

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
    int *mY)
{
    const int num_alt_i = Ai_end - Ai_begin + 1;
    const int num_alt_j = Aj_end - Aj_begin + 1;
    const int n_variable = num_alt_i + num_alt_j;
    const int effective = gsl_vector_int_get(effective_constraints, 0);

    const std::vector<LinCons> original_v = (*mX == -1) ? V_rows : LC_V;
    std::vector<LinCons> reduced_v;
    reduced_v.reserve(effective);
    for (int k = 1; k <= effective; ++k)
    {
        const int row = gsl_vector_int_get(effective_constraints, k);
        if (row < 0 || row >= static_cast<int>(original_v.size()))
            throw std::runtime_error("Invalid effective-constraint index");
        reduced_v.push_back(original_v[row]);
    }
    LC_V.swap(reduced_v);
    *mX = static_cast<int>(LC_V.size());

    const std::vector<LinCons> original_p = (*mY == -1) ? P_rows : LC_P;
    std::vector<double> lower(n_variable), upper(n_variable);
    int local = 0;
    for (int original = Ai_begin; original <= Ai_end; ++original, ++local)
    {
        lower[local] = P_lobo[original];
        upper[local] = P_upbo[original];
    }
    for (int original = Aj_begin; original <= Aj_end; ++original, ++local)
    {
        lower[local] = P_lobo[original];
        upper[local] = P_upbo[original];
    }
    const std::vector<int> retained = filter_effective_rows(
        original_p, lower, upper, BaseRelaxation::TwoSimplexes, num_alt_i,
        "Probability-base LP is infeasible or failed");
    std::vector<LinCons> reduced_p;
    reduced_p.reserve(retained.size());
    for (int row : retained) reduced_p.push_back(original_p[row]);
    LC_P.swap(reduced_p);
    *mY = static_cast<int>(LC_P.size());
}

/*
    redundancy() removes redundant linear constraints.

    inputs
        instance: the original instance;
        eff_cons: to record the indices of effective constraints;
        fixed_variables: to indicate whether the variable is fixed.

    outputs
        instance: the instance with the original linear constraints, but with fixed, original, or untightened lower and upper bounds;
        eff_cons: the indices of effective constraints with element 0 for the #;
        fixed_variables: the indicators for each variable, 1 for fixed and 0 otherwise, with element 0 for the #.
 */

int redundancy(struct dblp *instance,
               gsl_vector_int *eff_cons,
               gsl_vector_int *fixed_variables)
{
    const int n_variable = instance->Ax->size2, n_constraints = instance->Ax->size1;
    int j, counter = 0;
    REAL objective_value;

    /* check whether the interval imposed on each variable is too tight */
    gsl_vector_int_set_zero(fixed_variables); // each variable is considered free (0) at the beginning

    for (j = 0; j < n_variable; j++)
    {
        if (gsl_vector_get(instance->xu, j) == gsl_vector_get(instance->xl, j)) // only truly fixed variables
        {
            objective_value = (gsl_vector_get(instance->xu, j) + gsl_vector_get(instance->xl, j)) / 2; // to fix the variable with the mean of its interval
            gsl_vector_int_set(fixed_variables, 0, gsl_vector_int_get(fixed_variables, 0) + 1);        // increase the # of fixed variables
            gsl_vector_int_set(fixed_variables, j + 1, 1);                                             // the corresponding indicator is set to 1
            gsl_vector_set(instance->xl, j, objective_value);                                          // fix the variable
            gsl_vector_set(instance->xu, j, objective_value);                                          // fix the variable
        }
    }

    // show_gsl_vector_int(fixed_variables);

    if (gsl_vector_int_get(fixed_variables, 0) == n_variable) // nothing to do if all variables are fixed
    {
        return -1;
    }

    std::vector<LinCons> rows;
    rows.reserve(n_constraints);
    for (int i = 0; i < n_constraints; ++i)
    {
        LinCons row;
        row.b = gsl_vector_get(instance->bx, i);
        for (j = 0; j < n_variable; ++j)
        {
            const double value = gsl_matrix_get(instance->Ax, i, j);
            if (value != 0.0)
            {
                row.index.push_back(j + 1);
                row.coefficient.push_back(value);
            }
        }
        rows.push_back(row);
    }
    std::vector<double> lower(n_variable), upper(n_variable);
    for (j = 0; j < n_variable; ++j)
    {
        lower[j] = gsl_vector_get(instance->xl, j);
        upper[j] = gsl_vector_get(instance->xu, j);
    }
    const std::vector<int> retained = filter_effective_rows(
        rows, lower, upper, BaseRelaxation::Box, 0,
        "Value-base LP is infeasible or failed");
    gsl_vector_int_set(eff_cons, 0, static_cast<int>(retained.size()));
    for (size_t i = 0; i < retained.size(); ++i)
        gsl_vector_int_set(eff_cons, i + 1, retained[i]);

    return counter;
}

int safe_inverse(const gsl_matrix *A, gsl_matrix *A_inv, double tol)
{
    size_t n = A->size1;
    gsl_matrix *U = managed::matrix_alloc(n, n);
    gsl_matrix *V = managed::matrix_alloc(n, n);
    gsl_vector *S = managed::vector_alloc(n);
    gsl_vector *work = managed::vector_alloc(n);

    double s, inv_s, val, vi, ui;
    size_t i;

    gsl_matrix_memcpy(U, A);

    if (gsl_linalg_SV_decomp(U, V, S, work) != GSL_SUCCESS)
    {
        managed::matrix_free(U);
        managed::matrix_free(V);
        managed::vector_free(S);
        managed::vector_free(work);

        return -1;
    }

    for (i = 0; i < n; i++)
    {
        if (fabs(gsl_vector_get(S, i)) < tol)
        {
            managed::matrix_free(U);
            managed::matrix_free(V);
            managed::vector_free(S);
            managed::vector_free(work);

            return -1;
        }
    }

    gsl_matrix_set_zero(A_inv);

    for (i = 0; i < n; i++)
    {
        s = gsl_vector_get(S, i);
        inv_s = 1.0 / s;

        for (size_t j = 0; j < n; j++)
        {
            for (size_t k = 0; k < n; k++)
            {
                val = gsl_matrix_get(A_inv, j, k);
                vi = gsl_matrix_get(V, j, i);
                ui = gsl_matrix_get(U, k, i);
                gsl_matrix_set(A_inv, j, k, val + vi * inv_s * ui);
            }
        }
    }

    managed::matrix_free(U);
    managed::matrix_free(V);
    managed::vector_free(S);
    managed::vector_free(work);

    return 0;
}

int svd_solve_safe(gsl_matrix *A, gsl_vector *b, gsl_vector *x, double tol)
{
    size_t N = A->size1, j;
    int status;
    double sigma, dot, xi;

    if (A->size2 != N || b->size != N || x->size != N)
    {
        return -1;
    }

    gsl_matrix *U = managed::matrix_alloc(N, N);
    gsl_matrix *V = managed::matrix_alloc(N, N);
    gsl_vector *S = managed::vector_alloc(N);
    gsl_vector *tmp = managed::vector_alloc(N);
    gsl_vector *c = managed::vector_alloc(N);

    gsl_matrix_memcpy(U, A);

    status = gsl_linalg_SV_decomp(U, V, S, tmp);

    if (status != GSL_SUCCESS)
    {
        managed::matrix_free(U);
        managed::matrix_free(V);
        managed::vector_free(S);
        managed::vector_free(tmp);
        managed::vector_free(c);

        return -1;
    }

    for (size_t i = 0; i < N; i++)
    {
        sigma = gsl_vector_get(S, i);
        dot = 0.0;

        for (j = 0; j < N; j++)
        {
            dot += gsl_matrix_get(U, j, i) * gsl_vector_get(b, j);
        }

        gsl_vector_set(c, i, (fabs(sigma) > tol) ? dot / sigma : 0.0);
    }

    for (size_t i = 0; i < N; i++)
    {
        xi = 0.0;

        for (size_t j = 0; j < N; j++)
        {
            xi += gsl_matrix_get(V, i, j) * gsl_vector_get(c, j);
        }

        gsl_vector_set(x, i, xi);
    }

    managed::matrix_free(U);
    managed::matrix_free(V);
    managed::vector_free(S);
    managed::vector_free(tmp);
    managed::vector_free(c);

    return 0;
}

int results_BARON(d_row p_point, d_row v_point)
{
    FILE *fp = fopen("data/results.txt", "r");

    if (!fp)
    {
        perror("Cannot open file");
        return 1;
    }

    char var[10];
    double x_val, y_val;

    // Skip header line
    char line[100];
    fgets(line, sizeof(line), fp); // Skip the ": x y :=" line

    while (fgets(line, sizeof(line), fp))
    {
        if (strchr(line, ';'))
            break; // End of data

        int index = 0;
        x_val = y_val = 0.0;
        sscanf(line, "%s", var);

        if (var[0] == 'x')
        {
            sscanf(var + 1, "%d", &index);
            sscanf(line + strlen(var), "%lf", &x_val);
            v_point[index] = x_val;
        }
        else if (var[0] == 'y')
        {
            sscanf(var + 1, "%d", &index);
            char *last = strrchr(line, ' ');
            if (last)
            {
                y_val = std::stod(last);
                p_point[index] = y_val;
            }
        }
    }

    fclose(fp);

    // // Print the arrays
    // printf("x_values:\n");
    // for (int i = 0; i < dimension; i++)
    // {
    //     printf("x%d = %f\n", i + 1, v_point[i]);
    // }

    // printf("\ny_values:\n");
    // for (int i = 0; i < dimension; i++)
    // {
    //     printf("y%d = %f\n", i + 1, p_point[i]);
    // }

    return 0;
}

/*
    check_consistency() checks whether a solution is consistent with the constraints.

    inputs
        A, b: the constraint matrix and rhs;
        solution: the solution to be checked.

    output
        ConsistencyFlag: true for consistency, false otherwise.
*/

bool check_consistency_reverse_search(
    gsl_matrix *A,
    gsl_vector *b,
    gsl_vector *solution)
{
    size_t i;
    bool ConsistencyFlag = true;

    gsl_vector *result = managed::vector_calloc(b->size);
    gsl_vector_memcpy(result, b);

    gsl_blas_dgemv(CblasNoTrans, 1, A, solution, -1, result);
    resolve_zero_vector(result);

    for (i = 0; i < b->size; i++)
    {
        if (gsl_vector_get(result, i) > 0)
        {
            ConsistencyFlag = false;
            break;
        }
    }

    if (ConsistencyFlag)
    {
        for (i = 0; i < solution->size; i++)
        {
            if (gsl_vector_get(solution, i) < 0)
            {
                ConsistencyFlag = false;
                break;
            }
        }
    }

    managed::vector_free(result);
    return ConsistencyFlag;
}

// bool check_consistency_reverse_search(gsl_matrix *A,
//                                       gsl_vector *b,
//                                       gsl_vector *solution)
// {
//     int i, j, counter, temp_1, temp_2;
//     REAL value, sum;
//     const int m = A->size1, n = A->size2;

//     counter = 0;

//     for (i = 0; i < m; i++)
//     {
//         sum = 0;

//         for (j = 0; j < n; j++)
//         {
//             value = gsl_matrix_get(A, i, j);

//             if (fabs(value) > DBLP_EPSILON)
//             {
//                 if (value == 1)
//                 {
//                     sum = sum + gsl_vector_get(solution, j);
//                     temp_1 = j;
//                 }
//                 else
//                 {
//                     sum = sum - gsl_vector_get(solution, j);
//                     temp_2 = j;
//                 }
//             }
//         }

//         if (sum > gsl_vector_get(b, i) + 1e-5)
//         {
//             // if (counter == 0)
//             // {
//             //     printf("\n");
//             // }

//             // printf("the %dth linear constraint in Ax <= b is violated, v[%d] = %f, v[%d] = %f, v[%d] - v[%d] = %f, rhs = %f\n", i, temp_1 + 1, gsl_vector_get(solution, temp_1), temp_2 + 1, gsl_vector_get(solution, temp_2), temp_1 + 1, temp_2 + 1, sum, gsl_vector_get(b, i));

//             counter++;
//         }
//     }

//     if (counter > 0)
//     {
//         return (FALSE);
//     }
//     else
//     {
//         return (TRUE);
//     }
// }

/*
    consistency_check() checks whether a solution is consistent with the constraints.

    inputs
        dataset_in: ;
        vertex_x: the solution to be checked;
        noc: the # of constraints.

    output
        ConsistencyFlag: true for consistency, false otherwise.
*/

bool consistency_check(
    struct dataset *dataset_in,
    gsl_vector *vertex_x,
    int noc)
{
    int i, j, counter;
    REAL value, sum;
    const int m = noc, n = dataset_in->valid_A->size2 - 1;
    bool flag_a = true, flag_b = true;
    bool debug = false;

    counter = 0;

    for (i = 0; i < m; i++)
    {
        sum = 0;

        for (j = 1; j <= n; j++)
        {
            value = gsl_matrix_get(dataset_in->valid_A, i, j);

            if (fabs(value) > DBLP_EPSILON)
            {
                sum = sum + value * gsl_vector_get(vertex_x, j - 1);
            }
        }

        if (sum > gsl_matrix_get(dataset_in->valid_A, i, 0) + DBLP_OBJ_EPSILON)
        {
            if (counter == 0)
            {
                if (debug)
                {
                    printf("\n");
                }
            }

            if (debug)
            {
                printf("the %dth (No. %d) linear constraint in Ax <= b is violated, the sum is %f, rhs = %f\n", i, i + 1, sum, gsl_matrix_get(dataset_in->valid_A, i, 0));
            }

            counter++;
        }
    }

    if (counter > 0)
    {
        flag_a = false;
    }

    counter = 0;

    for (j = 0; j < n; j++)
    {
        value = gsl_vector_get(vertex_x, j);

        if ((value < gsl_vector_get(dataset_in->xl, j) - DBLP_EPSILON) || (value > gsl_vector_get(dataset_in->xu, j) + DBLP_EPSILON))
        {
            if (counter == 0)
            {
                if (debug)
                {
                    printf("\n");
                }
            }

            if (debug)
            {
                printf("the box constraints for variable %d is violated, v[%d]_l = %f, v[%d]_u = %f, v[%d] = %f\n", j + 1, j + 1, gsl_vector_get(dataset_in->xl, j), j + 1, gsl_vector_get(dataset_in->xu, j), j + 1, gsl_vector_get(vertex_x, j));
            }

            counter++;
        }
    }

    if (counter > 0)
    {
        flag_b = false;
    }

    if (flag_a && flag_b)
    {
        return (true);
    }
    else
    {
        return (false);
    }
}

void show_array_float(REAL my_array[])
{
    int i;

    printf("\n");

    for (i = 0; i <= my_array[0]; i++)
    {
        printf("% 7.4f", my_array[i]);
    }

    printf("\n");
}

void show_array_int(int my_array[])
{
    int i;

    printf("\n");

    for (i = 0; i <= my_array[0]; i++)
    {
        printf("% 4d", my_array[i]);
    }

    printf("\n");
}

/*
    construct_compelete_solution() establishes the solution vector in full dimensions.

    inputs
        vector_new: the vector to save the full dimensional solution vector;
        basic_index: the indices of basic variables;
        column_solution: the solution corresponding to basic_index.

    output
        vector_new: the established full dimensional solution vector.
*/

void construct_complete_solution(gsl_vector *vector_new,
                                 gsl_vector_int *basic_index,
                                 gsl_vector *column_solution)
{
    int i;
    const int m = column_solution->size;

    gsl_vector_set_zero(vector_new);

    for (i = 0; i < m; i++)
    {
        gsl_vector_set(vector_new, gsl_vector_int_get(basic_index, i + 1), gsl_vector_get(column_solution, i));
    }

    gsl_vector_set(vector_new, 0, vector_new->size - 1);
}

// /*
//     condition_number() tells whether a square matrix is invertible using its condition number

//     inputs
//         A: the square matrix whose condition number will be evaluated;
//         threshhold: 1e12 by default.

//     output
//         flag: true for an invertible square matrix, false otherwise.
// */

// bool condition_number(gsl_matrix *A, double threshold)
// {
//     int n = A->size1, status;
//     double s_max, s_min, cond;
//     bool flag;

//     gsl_matrix *U = managed::matrix_alloc(n, n);
//     gsl_matrix_memcpy(U, A);
//     gsl_matrix *V = managed::matrix_calloc(n, n);
//     gsl_vector *S = managed::vector_calloc(n);
//     gsl_vector *work = managed::vector_calloc(n);

//     status = gsl_linalg_SV_decomp(U, V, S, work);

//     if (status != GSL_SUCCESS) {
//         flag = false;
//     } else {
//         // 按降序排列
//         s_max = gsl_vector_get(S, 0);
//         s_min = gsl_vector_get(S, n - 1);

//         double eps = 1e-15; // double 精度机器精度
//         if (fabs(s_min) < eps * fabs(s_max)) {
//             flag = false; // 认为不可逆
//         } else {
//             cond = s_max / s_min;
//             flag = cond < threshold; // 阈值判断病态程度
//         }
//     }

//     managed::matrix_free(U);
//     managed::matrix_free(V);
//     managed::vector_free(S);
//     managed::vector_free(work);

//     return flag;
// }

// bool condition_number(gsl_matrix *A, double threshold)
// {
//     int n = A->size1, status, i;
//     bool flag;
//     double s_max, s_min, s, cond;

//     gsl_matrix *U = managed::matrix_alloc(n, n);
//     gsl_matrix_memcpy(U, A);
//     gsl_matrix *V = managed::matrix_calloc(n, n);
//     gsl_vector *S = managed::vector_calloc(n);
//     gsl_vector *work = managed::vector_calloc(n);

//     status = gsl_linalg_SV_decomp(U, V, S, work);

//     if (status != GSL_SUCCESS)
//     {
//         flag = false;
//     }
//     else
//     {
//         s_max = gsl_vector_get(S, 0);
//         s_min = gsl_vector_get(S, 0);

//         for (i = 1; i < n; ++i)
//         {
//             s = gsl_vector_get(S, i);

//             if (s > s_max)
//             {
//                 s_max = s;
//             }

//             if (s < s_min)
//             {
//                 s_min = s;
//             }
//         }

//         if (s_min < 1e-15)
//         {
//             flag = false;
//         }
//         else
//         {
//             cond = s_max / s_min;
//             flag = cond < threshold;
//         }
//     }

//     managed::matrix_free(U);
//     managed::matrix_free(V);
//     managed::vector_free(S);
//     managed::vector_free(work);
//     return flag;
// }

/*
    setDifference() gives set3 = set1 \ set2

    inputs
        set1: the large set containing set2;
        set2: the elements to be omitted from set1;
        set3: to save the identified set1 \ set2.

    output
        set3: set1 \ set2.
*/

void setDifference(std::vector<int> &set1, std::vector<int> &set2, std::vector<int> &set3)
{
    int isInSecondSet, counter = 1, i, j;

    for (i = 1; i <= set1[0]; i++)
    {
        isInSecondSet = 0;

        for (j = 1; j <= set2[0]; j++)
        {
            if (set1[i] == set2[j])
            {
                isInSecondSet = 1;
                break;
            }
        }

        if (!isInSecondSet)
        {
            set3[counter] = set1[i];
            counter++;
        }
    }

    set3[0] = counter - 1;
}

/*
    nchoosek() computes the # of combinations of choosing k from n.

    inputs
        n: the total # of elements;
        k: the # of elements to be selected.

    output
        n! / (k! * (n-k)!).
*/

unsigned long long nchoosek(int n, int k)
{
    int i;

    if (k > n)
    {
        return 0;
    }

    if (k > n - k)
    {
        k = n - k;
    }

    unsigned long long result = 1;

    for (i = 1; i <= k; i++)
    {
        result *= n--;
        result /= i;
    }

    return result;
}

/*
    nchoosek_enum() only gives the ith combination of choosing k from n elements.

    inputs
        n, k: choose k from n elements;
        i: the ith combination.

    outputs
        result: the ith combination.
*/

void nchoosek_enum(int n, int k, int i, std::vector<int> &result)
{
    int index = 0, j, m;
    unsigned long long combinations;

    for (j = 0; j < k; j++)
    {
        for (m = index; m < n; m++)
        {
            combinations = nchoosek(n - m - 1, k - j - 1);

            if (combinations >= (unsigned long long)i)
            {
                result[j] = m + 1;
                index = m + 1;
                break;
            }
            else
            {
                i -= combinations;
            }
        }
    }
}

void resolve_zero_matrix(gsl_matrix *m)
{
    size_t i, j;

    for (i = 0; i < m->size1; i++)
    {
        for (j = 0; j < m->size2; j++)
        {
            if (fabs(gsl_matrix_get(m, i, j)) < DBLP_EPSILON)
            {
                gsl_matrix_set(m, i, j, 0);
            }
        }
    }
}

void resolve_zero_vector(gsl_vector *v)
{
    size_t i;

    for (i = 0; i < v->size; i++)
    {
        if (fabs(gsl_vector_get(v, i)) < DBLP_EPSILON)
        {
            gsl_vector_set(v, i, 0);
        }
    }
}

void show_gsl_matrix_int_rows(gsl_matrix_int *p, int rows)
{
    int i, j, columns = p->size2;

    for (i = 0; i < rows; i++)
    {
        printf("\n");

        for (j = 0; j < columns; j++)
            printf("%2d ", gsl_matrix_int_get(p, i, j));
    }

    printf("\n");
}

void show_gsl_matrix_cuts(gsl_matrix *p, int n_cuts)
{
    int i, j, rows = n_cuts, columns = p->size2;

    for (i = 0; i < rows; i++)
    {
        printf("\n");

        for (j = 0; j < columns; j++)
            printf("% 7.5f ", gsl_matrix_get(p, i, j));
    }

    printf("\n");
}

// void fpe_handler(int signum)
// {
//     printf("Floating point exception detected! (Signal: %d)\n", signum);
//     exit(EXIT_FAILURE);
// }

void show_gsl_matrix_neighbors(gsl_matrix *neighbors, size_t num_neighbors)
{
    size_t i, j;

    for (i = 0; i < neighbors->size1; i++)
    {
        printf("\n");

        for (j = 0; j < num_neighbors; j++)
            printf("% 7.2f ", gsl_matrix_get(neighbors, i, j));
    }

    printf("\n");
}

void show_gsl_matrix(gsl_matrix *p)
{
    int i, j, rows = p->size1, columns = p->size2;

    for (i = 0; i < rows; i++)
    {
        printf("\n");

        for (j = 0; j < columns; j++)
            printf("% 8.4f ", gsl_matrix_get(p, i, j));
    }

    printf("\n");
}

void show_gsl_matrix_columns_int(gsl_matrix_int *p, int columns)
{
    int i, j, rows = p->size1;

    for (i = 0; i < rows; i++)
    {
        printf("\n");

        for (j = 0; j < columns; j++)
        {
            printf("% .2d ", gsl_matrix_int_get(p, i, j));
        }
    }

    printf("\n");
}

void show_gsl_matrix_columns(gsl_matrix *p, int columns)
{
    int i, j, rows = p->size1;

    for (i = 0; i < rows; i++)
    {
        printf("\n");

        for (j = 0; j < columns; j++)
        {
            printf("% 10.4f", gsl_matrix_get(p, i, j));
        }
    }

    fflush(stdout);
}

void show_gsl_matrix_int(gsl_matrix_int *p)
{
    int i, j, rows = p->size1, columns = p->size2;

    for (i = 0; i < rows; i++)
    {
        printf("\n");

        for (j = 0; j < columns; j++)
            printf("% 7.2d ", gsl_matrix_int_get(p, i, j));
    }

    printf("\n");
}

void show_gsl_vector(gsl_vector *p)
{
    size_t i;

    printf("\n");

    for (i = 0; i < p->size; i++)
    {
        printf("% 8.4f ", gsl_vector_get(p, i));
    }

    fflush(stdout);
}

void show_gsl_vector_int_all(gsl_vector_int *p)
{
    int i, ceiling = p->size;

    for (i = 0; i < ceiling; i++)
    {
        printf("%2d ", gsl_vector_int_get(p, i));
    }

    printf("\n");
}

void show_gsl_vector_int_counter(gsl_vector_int *p, int start, int end)
{
    int i;

    for (i = start; i < end; i++)
    {
        printf("%2d ", gsl_vector_int_get(p, i));
    }

    printf("\n");
}

void show_gsl_vector_counter(gsl_vector *p, int start, int end)
{
    int i;

    for (i = start; i < end; i++)
    {
        printf("% 8.4f ", gsl_vector_get(p, i));
    }

    printf("\n");
}

void show_gsl_vector_int(gsl_vector_int *p)
{
    int i, ceiling = gsl_vector_int_get(p, 0) + 1;

    if (gsl_vector_int_get(p, 0) == -1)
    {
        printf("%2d ", gsl_vector_int_get(p, 0));
    }
    else
    {
        for (i = 0; i < ceiling; i++)
        {
            printf("%2d ", gsl_vector_int_get(p, i));
        }
    }

    printf("\n");
}

/*
    show_cutting_values() calculates and shows the objective values at the derived cutting points.
    it is used to judge the correctness of the derived cutting points.
    for a normal cutting point, i.e., "cuts_indicator != 0", the objective value should equal CBOV.

    inputs
        cutting_points: the derived cutting points in columns with their corresponding elements in "cuts_indicator";
        instance: the structure containing the original problem;
        lp_Y: the already constructed LP model for Y space;
        cuts_indicator: 1 for positive, -1 for negative, and 0 for abnormal.

    outputs
        see the information on the screen.
*/

void show_cutting_values(gsl_vector *b, struct dataset *dataset_in, LpModel *lp_Y, int cuts_indicator, int No)
{
    const int n_variable_y = dataset_in->dy->size;
    std::vector<REAL> temp_array_y(n_variable_y + 1);
    REAL obj_value;

    gsl_vector *X = managed::vector_calloc(b->size);
    gsl_vector *temp_y = managed::vector_calloc(n_variable_y);

    temp_array_y[0] = 0;

    // derive the coefficient of y, Q' * b + dy
    gsl_vector_memcpy(X, b);
    gsl_vector_memcpy(temp_y, dataset_in->dy);
    gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, X, 1, temp_y);

    if (cuts_indicator == 1)
    {
        gsl_vector_scale(temp_y, -1); // a max LP along the negative extension of an edge
    }

    gsl_vector_view v = gsl_vector_view_array(temp_array_y.data(), 1 + n_variable_y);
    gsl_vector_view dest_v = gsl_vector_subvector(&v.vector, 1, n_variable_y);
    gsl_vector_memcpy(&dest_v.vector, temp_y);

    set_obj_fn(lp_Y, temp_array_y.data());
    solve(lp_Y); // always solve a min LP

    if (cuts_indicator == 1)
    {
        obj_value = -get_objective(lp_Y); // retrieve the proper objective value for a max LP
    }
    else
    {
        obj_value = get_objective(lp_Y);
    }

    printf("the objective value of the %dth cutting point is %f\n", No, obj_value);

    managed::vector_free(temp_y);
    managed::vector_free(X);
}

/*
    check_objective() computes the objective value based on given vectors.
    it is used to check the consistency of computations from another aspect.

    inputs
        x, y: the solution vectors in X and Y, respectively;
        instance: the original dblp instance.

    outputs
        the information of the objective value on the screen.
*/

void check_objective(gsl_vector *x, gsl_vector *y, struct dblp *instance)
{
    REAL obj_value_y;

    gsl_vector *temp_y = managed::vector_calloc(y->size);

    /* derive (Q' * x + dy)' * y */
    gsl_vector_set_zero(temp_y);
    gsl_blas_dgemv(CblasTrans, 1, instance->Q, x, 1, temp_y);
    gsl_blas_ddot(temp_y, y, &obj_value_y);

    printf("\nthe objective value is: %f\n", obj_value_y);

    managed::vector_free(temp_y);
}
