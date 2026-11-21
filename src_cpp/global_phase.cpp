#include "dblp.h"

double vector_norm(const double *v, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += v[i] * v[i];
    return sqrt(sum);
}

int is_good_vector(const double *v, int n)
{
    double norm = vector_norm(v, n);

    if (norm < DBLP_EPSILON || norm > 1 / DBLP_EPSILON) // the norm of v should not be too small or large
    {
        return 0;
    }

    int nonzero_count = 0;
    double max_val = 0.0, val = 0.0;

    for (int i = 0; i < n; i++)
    {
        val = fabs(v[i]);

        if (val > DBLP_EPSILON)
        {
            nonzero_count++;
        }

        if (val > max_val)
        {
            max_val = val;
        }
    }

    if (nonzero_count < n / 3) // the proportion of nonzero elements should be more than 2 / 3
    {
        return 0;
    }

    double ratio = max_val / norm; // the vector should not be dominated by a single element

    if (ratio > 0.99)
    {
        return 0;
    }

    return 1;
}

/*
    active_constraint() locates an active constraint constructed by benchmark and n - 1 adjacent vertices.

    inputs
        points: vertices (in rows) adjacent to the benchmark vertex, which occupies the last row;
        A, b: the constraint set;
        which_constraint: the starting active constraint used to resume the last search;
        indices_new: to record the row indices of the vertices in "points" lying on the located constraint.

    outputs
        flag: true for an active constraint with at least n vertices on it, false otherwise;
        indices_new: row indices of the vertices in "points" lying on the located constraint;
        which_constraint: the oridinal number of an active constraint where we have located n vertices on it.
*/

bool active_constraint(gsl_matrix *points,
                       gsl_matrix *A,
                       gsl_vector *b,
                       int *which_constraint,
                       gsl_vector_int *indices_new)
{
    const int rows = points->size1, columns = points->size2, n_constraints = A->size1;
    int i, j;
    bool flag = false;

    gsl_vector_int *indices_active = managed::vector_int_calloc(n_constraints);
    gsl_vector *result = managed::vector_calloc(n_constraints), *residual = managed::vector_calloc(rows);
    gsl_vector_view src_v;

    /* result = A * benchmark - b */
    gsl_vector_memcpy(result, b);
    src_v = gsl_matrix_row(points, rows - 1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, &src_v.vector, -1, result);
    resolve_zero_vector(result);

    /* locate an active constraint with enough vertices */
    for (j = *which_constraint; j < n_constraints; j++)
    {
        if (gsl_vector_get(result, j) == 0)
        {
            src_v = gsl_matrix_row(A, j);
            gsl_vector_set_all(residual, gsl_vector_get(b, j));
            gsl_blas_dgemv(CblasNoTrans, 1.0, points, &src_v.vector, -1, residual);
            resolve_zero_vector(residual);

            gsl_vector_int_set(indices_new, 0, 0);

            for (i = 0; i < rows; i++)
            {
                if (gsl_vector_get(residual, i) == 0)
                {
                    gsl_vector_int_set(indices_new, 0, gsl_vector_int_get(indices_new, 0) + 1);
                    gsl_vector_int_set(indices_new, gsl_vector_int_get(indices_new, 0), i);
                }
            }

            if (gsl_vector_int_get(indices_new, 0) >= columns)
            {
                *which_constraint = j;
                flag = true;
                break;
            }
        }
    }

    managed::vector_free(result);
    managed::vector_free(residual);
    managed::vector_int_free(indices_active);

    return flag;
}

/*
    conservativeCut_new() locates a conservative hyperplane by a smart idea.

    inputs
        points: vertices (in rows) adjacent to the benchmark vertex, which occupies the last row;
        indices_new: indices used to set up the current extreme face;
        normal_vector: to save the normal vaector of a identified conservative hyperplane.

    outputs
        flag: GSL_LU, 1 or 0;
        normal_vector: the normal vector of a conservative hyperplane, valid only when flag is true.
*/

int conservativeCut_New(gsl_matrix *points,
                        gsl_vector_int *indices_new,
                        gsl_vector *normal_vector)
{
    const int rows = points->size1, columns = points->size2, n_remainings = rows - columns;
    int flag = 0, i, j, counter;
    std::vector<int> remainings(1 + n_remainings);
    std::vector<int> row_indices(1 + rows);
    std::vector<int> rest_indices(1 + n_remainings);
    double criterion;
    bool indicator;

    for (i = 0; i < rows; i++)
    {
        row_indices[i + 1] = i;
    }

    row_indices[0] = rows;

    gsl_vector *benchmark = managed::vector_calloc(columns), *rhs = managed::vector_calloc(columns), *rhs_rest_vertices = managed::vector_calloc(n_remainings - 1);
    gsl_vector_int *candidate = managed::vector_int_calloc(1 + columns);
    gsl_matrix *chosenPlane = managed::matrix_calloc(columns, columns), *rest_vertices = managed::matrix_calloc(n_remainings - 1, columns);
    gsl_permutation *p = managed::permutation_alloc(columns);
    // gsl_matrix *chosenPlane_copy = managed::matrix_calloc(columns, columns);

    gsl_vector_view row_from_points;
    gsl_vector_int_view v_src, v_dest;

    v_src = gsl_vector_int_subvector(indices_new, 0, columns);
    v_dest = gsl_vector_int_subvector(candidate, 0, columns);
    gsl_vector_int_memcpy(&v_dest.vector, &v_src.vector);
    gsl_vector_int_set(candidate, 0, gsl_vector_int_get(indices_new, 0) - 1);

    gsl_matrix_get_row(benchmark, points, rows - 1);

    std::vector<int> indices_new_vec(indices_new->data, indices_new->data + indices_new->size);
    setDifference(row_indices, indices_new_vec, remainings);

    gsl_vector_set_all(rhs, 1.0);

    for (i = 0; i < remainings[0]; i++)
    {
        for (j = 0; j < columns - 1; j++)
        {
            row_from_points = gsl_matrix_row(points, gsl_vector_int_get(candidate, j + 1));
            gsl_matrix_set_row(chosenPlane, j, &row_from_points.vector);
        }

        row_from_points = gsl_matrix_row(points, remainings[i + 1]);
        gsl_matrix_set_row(chosenPlane, columns - 1, &row_from_points.vector);

        /* to derive the location of the benchmark w.r.t. the constructed hyperplane */
        // gsl_matrix_memcpy(chosenPlane_copy, chosenPlane);

        // if (gsl_linalg_LU_decomp(chosenPlane, p, &s) != GSL_SUCCESS)
        // {
        //     continue;
        // }

        // if (gsl_linalg_LU_solve(chosenPlane, p, rhs, normal_vector) != GSL_SUCCESS)
        // {
        //     continue;
        // }

        if (svd_solve_safe(chosenPlane, rhs, normal_vector, 1e-12) != 0)
        {
            continue;
        }

        gsl_blas_ddot(normal_vector, benchmark, &criterion);
        criterion--; // criterion - 1

        if (fabs(criterion) > DBLP_EPSILON)
        {
            /* to determine whether all remaining vertices are lying on the other side */
            gsl_vector_int_set(candidate, 0, columns);
            gsl_vector_int_set(candidate, columns, remainings[i + 1]);

            std::vector<int> indices_new_vec(candidate->data, candidate->data + candidate->size);
            setDifference(row_indices, indices_new_vec, rest_indices);

            counter = 0;

            for (j = 1; j <= rest_indices[0]; j++)
            {
                if (rest_indices[j] != rows - 1)
                {
                    row_from_points = gsl_matrix_row(points, rest_indices[j]);
                    gsl_matrix_set_row(rest_vertices, counter++, &row_from_points.vector);
                }
            }

            gsl_vector_set_all(rhs_rest_vertices, 1);
            gsl_blas_dgemv(CblasNoTrans, 1, rest_vertices, normal_vector, -1, rhs_rest_vertices);

            indicator = false;

            if (criterion > 0)
            {
                for (j = 0; j < n_remainings - 1; j++)
                {
                    if (gsl_vector_get(rhs_rest_vertices, j) > 0)
                    {
                        break;
                    }
                }

                if (j == n_remainings - 1)
                {
                    indicator = true;
                }
            }

            if (criterion < 0)
            {
                for (j = 0; j < n_remainings - 1; j++)
                {
                    if (gsl_vector_get(rhs_rest_vertices, j) < 0)
                    {
                        break;
                    }
                }

                if (j == n_remainings - 1)
                {
                    indicator = true;
                }
            }

            if (indicator)
            {
                flag = 1;
                break;
            }
        }
    }

    managed::vector_free(benchmark);
    managed::vector_int_free(candidate);
    managed::matrix_free(chosenPlane);
    managed::permutation_free(p);
    managed::vector_free(rhs);
    managed::vector_free(rhs_rest_vertices);
    managed::matrix_free(rest_vertices);

    return flag;
}

/* Find a supporting separator without enumerating candidate facets.  The LP
   maximizes the L1-normalized margin between the current pseudo-minimizer and
   the convex hull of all recovered cutting points. */
bool separating_cut_lp(gsl_matrix *points, gsl_vector *normal_vector,
                       double *cut_rhs)
{
    if (!points || points->size1 < 2 || points->size2 < 1) return false;
    const int n = static_cast<int>(points->size2);
    const int point_count = static_cast<int>(points->size1) - 1;
    const int columns = 2 * n + 1;
    LpModel *lp = make_lp(point_count + 1, columns);
    set_verbose(lp, kLpNeutral);
    set_maxim(lp);
    for (int column = 1; column <= columns; ++column)
        set_upbo(lp, column, 1.0);

    std::vector<REAL> row(columns + 1, 0.0), objective(columns + 1, 0.0);
    gsl_vector_const_view benchmark = gsl_matrix_const_row(points, point_count);
    for (int i = 0; i < point_count; ++i)
    {
        std::fill(row.begin(), row.end(), 0.0);
        for (int j = 0; j < n; ++j)
        {
            const double difference = gsl_vector_get(&benchmark.vector, j) -
                                      gsl_matrix_get(points, i, j);
            row[1 + j] = difference;
            row[1 + n + j] = -difference;
        }
        row[columns] = -1.0;
        set_row(lp, i + 1, row.data());
        set_constr_type(lp, i + 1, kLpGe);
        set_rh(lp, i + 1, 0.0);
    }
    std::fill(row.begin(), row.end(), 0.0);
    for (int j = 0; j < 2 * n; ++j) row[1 + j] = 1.0;
    set_row(lp, point_count + 1, row.data());
    set_constr_type(lp, point_count + 1, kLpLe);
    set_rh(lp, point_count + 1, 1.0);
    objective[columns] = 1.0;
    set_obj_fn(lp, objective.data());

    const int lp_status = solve(lp);
    const double lp_margin = lp_status == OPTIMAL ? get_objective(lp) : 0.0;
    if (std::getenv("DECIDEIT_TRACE"))
        fprintf(stderr, "TRACE {\"event\":\"neighbor_separator_lp\",\"status\":%d,\"margin\":%.12g,\"euclidean_margin_upper_bound\":%.12g,\"dimension\":%d,\"neighbors\":%d}\n",
                lp_status, lp_margin, std::sqrt(static_cast<double>(n)) * lp_margin,
                n, point_count);
    if (lp_status != OPTIMAL || lp_margin <= 1e-10)
    {
        delete_lp(lp);
        return false;
    }
    std::vector<double> solution(columns, 0.0);
    get_variables(lp, solution.data());
    delete_lp(lp);
    for (int j = 0; j < n; ++j)
        gsl_vector_set(normal_vector, j, solution[j] - solution[n + j]);
    // Do not round small separator coefficients away here: a tiny positive
    // margin may rely on their combined contribution.  The final geometric
    // verification and the tableau's own tolerances decide usability.
    if (vector_norm(normal_vector->data, n) <= 1e-10) return false;

    *cut_rhs = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < point_count; ++i)
    {
        double value = 0.0;
        gsl_vector_const_view point = gsl_matrix_const_row(points, i);
        gsl_blas_ddot(normal_vector, &point.vector, &value);
        if (value > *cut_rhs) *cut_rhs = value;
    }
    double benchmark_value = 0.0;
    gsl_blas_ddot(normal_vector, &benchmark.vector, &benchmark_value);
    return std::isfinite(*cut_rhs) && benchmark_value > *cut_rhs + 1e-10;
}

static bool verify_neighbor_separator(gsl_matrix *points, gsl_vector *normal,
                                      double *rhs, double *margin)
{
    const size_t n = points->size2;
    const size_t neighbors = points->size1 - 1;
    const double norm = vector_norm(normal->data, static_cast<int>(n));
    if (!std::isfinite(norm) || norm <= 1e-12) return false;
    gsl_vector_scale(normal, 1.0 / norm);
    double highest_neighbor = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < neighbors; ++i)
    {
        double value = 0.0;
        for (size_t j = 0; j < n; ++j)
            value += gsl_vector_get(normal, j) * gsl_matrix_get(points, i, j);
        if (!std::isfinite(value)) return false;
        if (value > highest_neighbor) highest_neighbor = value;
    }
    double vertex_value = 0.0;
    for (size_t j = 0; j < n; ++j)
        vertex_value += gsl_vector_get(normal, j) *
                        gsl_matrix_get(points, neighbors, j);
    const double scale = std::fabs(vertex_value) > std::fabs(highest_neighbor)
        ? std::fabs(vertex_value) : std::fabs(highest_neighbor);
    // Reverse-search tableaux already round entries below DBLP_EPSILON to
    // zero.  A cut with a smaller geometric gap would create new vertices
    // that the neighborhood code cannot distinguish reliably.
    const double relative_tolerance = 1e-10 * (1.0 + scale);
    const double tolerance = relative_tolerance > DBLP_EPSILON
        ? relative_tolerance : DBLP_EPSILON;
    *margin = vertex_value - highest_neighbor;
    if (!std::isfinite(*margin) || *margin <= 8 * tolerance)
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"neighbor_separator_rejected\",\"margin\":%.12g,\"minimum\":%.12g,\"scale\":%.12g,\"neighbors\":%zu}\n",
                    *margin, 8 * tolerance, scale, neighbors);
        return false;
    }
    // Keep a numerical cushion for all neighbors without losing the strict
    // separation of the current vertex.
    double guard_factor = 2.0;
    if (const char *guard = std::getenv("DECIDEIT_NEIGHBOR_GUARD_FACTOR"))
    {
        char *end = nullptr;
        const double parsed = std::strtod(guard, &end);
        if (end == guard || *end != '\0' || !std::isfinite(parsed) ||
            parsed < 0.0 || parsed > 4.0)
            throw std::invalid_argument("neighbor guard factor must be in [0,4]");
        guard_factor = parsed;
    }
    *rhs = highest_neighbor + guard_factor * tolerance;
    return true;
}

// The F3 affine-hull branch of the Complete-Boundary construction applied to
// the first layer. If all neighbors span one hyperplane, its normal gives a
// cut supported by every neighbor without a separation LP or facet search.
static bool affine_neighbor_normal(gsl_matrix *points, gsl_vector *normal)
{
    const size_t n = points->size2;
    const size_t neighbors = points->size1 - 1;
    if (neighbors < n) return false;
    std::vector<std::vector<double>> basis;
    double largest = 0.0;
    for (size_t i = 1; i < neighbors; ++i)
    {
        double squared = 0.0;
        for (size_t j = 0; j < n; ++j)
        {
            const double difference = gsl_matrix_get(points, i, j) -
                                      gsl_matrix_get(points, 0, j);
            squared += difference * difference;
        }
        const double length = std::sqrt(squared);
        if (length > largest) largest = length;
    }
    if (largest <= 1e-12) return false;
    const double rank_tolerance = 1e-9 * (largest > 1.0 ? largest : 1.0) * n;
    for (size_t i = 1; i < neighbors; ++i)
    {
        std::vector<double> direction(n);
        for (size_t j = 0; j < n; ++j)
            direction[j] = gsl_matrix_get(points, i, j) -
                           gsl_matrix_get(points, 0, j);
        for (int pass = 0; pass < 2; ++pass)
            for (const auto &column : basis)
            {
                double projection = 0.0;
                for (size_t j = 0; j < n; ++j)
                    projection += direction[j] * column[j];
                for (size_t j = 0; j < n; ++j)
                    direction[j] -= projection * column[j];
            }
        double squared = 0.0;
        for (double value : direction) squared += value * value;
        const double length = std::sqrt(squared);
        if (length <= rank_tolerance) continue;
        for (double &value : direction) value /= length;
        basis.push_back(std::move(direction));
        if (basis.size() == n) return false;
    }
    if (basis.size() != n - 1) return false;

    double best_length = 0.0;
    std::vector<double> best(n);
    for (size_t axis = 0; axis < n; ++axis)
    {
        std::vector<double> candidate(n, 0.0);
        candidate[axis] = 1.0;
        for (int pass = 0; pass < 2; ++pass)
            for (const auto &column : basis)
            {
                double projection = 0.0;
                for (size_t j = 0; j < n; ++j)
                    projection += candidate[j] * column[j];
                for (size_t j = 0; j < n; ++j)
                    candidate[j] -= projection * column[j];
            }
        double squared = 0.0;
        for (double value : candidate) squared += value * value;
        const double length = std::sqrt(squared);
        if (length > best_length)
        {
            best_length = length;
            best.swap(candidate);
        }
    }
    if (best_length <= 1e-10) return false;
    double center_gap = 0.0;
    for (size_t j = 0; j < n; ++j)
        center_gap += best[j] * (gsl_matrix_get(points, neighbors, j) -
                                 gsl_matrix_get(points, 0, j));
    const double orientation = center_gap >= 0.0 ? 1.0 : -1.0;
    for (size_t j = 0; j < n; ++j)
        gsl_vector_set(normal, j, orientation * best[j] / best_length);
    return true;
}

// A cheap representation-independent candidate.  It often avoids a
// separator LP when active general-linear normals nearly cancel, but is
// accepted only after checking every enumerated neighbor.
static void centroid_neighbor_normal(gsl_matrix *points, gsl_vector *normal)
{
    const size_t neighbors = points->size1 - 1;
    for (size_t j = 0; j < points->size2; ++j)
    {
        long double mean = 0.0L;
        for (size_t i = 0; i < neighbors; ++i)
            mean += gsl_matrix_get(points, i, j);
        mean /= neighbors;
        gsl_vector_set(normal, j,
                       gsl_matrix_get(points, neighbors, j) -
                           static_cast<double>(mean));
    }
}

bool conservative_neighbor_cut(gsl_matrix *points, gsl_matrix *valid_A,
                               int valid_rows, const gsl_vector *xl,
                               const gsl_vector *xu, gsl_vector *normal,
                               double *rhs)
{
    if (!points || points->size1 < 2 || points->size2 < 1 || !normal ||
        !rhs || points->size2 != normal->size ||
        (valid_A && (valid_rows < 0 ||
                     static_cast<size_t>(valid_rows) > valid_A->size1 ||
                     valid_A->size2 != points->size2 + 1))) return false;

    const size_t n = points->size2;
    const size_t vertex_row = points->size1 - 1;
    gsl_vector_set_zero(normal);
    double margin = 0.0;
    const char *affine_mode = std::getenv("DECIDEIT_NEIGHBOR_AFFINE_HULL");
    if (affine_mode && std::strcmp(affine_mode, "1") == 0 &&
        affine_neighbor_normal(points, normal) &&
        verify_neighbor_separator(points, normal, rhs, &margin))
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"conservative_neighbor_cut\",\"mode\":\"affine_hull\",\"neighbors\":%zu,\"margin\":%.12g}\n",
                    vertex_row, margin);
        return true;
    }
    gsl_vector_set_zero(normal);
    int active_count = 0;
    if (valid_A)
    {
        for (int i = 0; i < valid_rows; ++i)
        {
            double value = 0.0, norm_sq = 0.0;
            for (size_t j = 0; j < n; ++j)
            {
                const double a = gsl_matrix_get(valid_A, i, j + 1);
                value += a * gsl_matrix_get(points, vertex_row, j);
                norm_sq += a * a;
            }
            const double bound = gsl_matrix_get(valid_A, i, 0);
            if (norm_sq <= 1e-24 || !std::isfinite(value) ||
                std::fabs(value - bound) >
                    1e-8 * (1.0 + std::fabs(value) + std::fabs(bound)))
                continue;
            const double inv_norm = 1.0 / std::sqrt(norm_sq);
            for (size_t j = 0; j < n; ++j)
                gsl_vector_set(normal, j, gsl_vector_get(normal, j) +
                               gsl_matrix_get(valid_A, i, j + 1) * inv_norm);
            ++active_count;
        }
    }
    for (size_t j = 0; j < n; ++j)
    {
        const double value = gsl_matrix_get(points, vertex_row, j);
        if (xl && std::fabs(value - gsl_vector_get(xl, j)) <=
                  1e-8 * (1.0 + std::fabs(value) +
                          std::fabs(gsl_vector_get(xl, j))))
        {
            gsl_vector_set(normal, j, gsl_vector_get(normal, j) - 1.0);
            ++active_count;
        }
        if (xu && std::fabs(value - gsl_vector_get(xu, j)) <=
                  1e-8 * (1.0 + std::fabs(value) +
                          std::fabs(gsl_vector_get(xu, j))))
        {
            gsl_vector_set(normal, j, gsl_vector_get(normal, j) + 1.0);
            ++active_count;
        }
    }

    const char *cut_mode = std::getenv("DECIDEIT_CONSERVATIVE_CUT");
    const bool lp_first = cut_mode && std::strcmp(cut_mode, "lp") == 0;
    bool active_valid = false;
    double active_rhs = 0.0, active_margin = 0.0;
    const char *backup_mode = nullptr;
    std::vector<double> active_normal;
    if (!lp_first && active_count &&
        verify_neighbor_separator(points, normal, rhs, &margin))
    {
        if (margin >= 1e-4)
        {
            if (std::getenv("DECIDEIT_TRACE"))
                fprintf(stderr, "TRACE {\"event\":\"conservative_neighbor_cut\",\"mode\":\"active_normal\",\"neighbors\":%zu,\"margin\":%.12g}\n",
                        vertex_row, margin);
            return true;
        }
        active_valid = true;
        active_rhs = *rhs;
        active_margin = margin;
        backup_mode = "active_normal";
        active_normal.assign(normal->data, normal->data + n);
    }

    if (!lp_first)
    {
        centroid_neighbor_normal(points, normal);
        if (verify_neighbor_separator(points, normal, rhs, &margin))
        {
            if (margin >= 1e-4)
            {
                if (std::getenv("DECIDEIT_TRACE"))
                    fprintf(stderr, "TRACE {\"event\":\"conservative_neighbor_cut\",\"mode\":\"centroid_normal\",\"neighbors\":%zu,\"margin\":%.12g}\n",
                            vertex_row, margin);
                return true;
            }
            if (!active_valid || margin > active_margin)
            {
                active_valid = true;
                active_rhs = *rhs;
                active_margin = margin;
                backup_mode = "centroid_normal";
                active_normal.assign(normal->data, normal->data + n);
            }
        }
    }

    // Strict separation of a vertex from all its adjacent vertices implies
    // it uniquely maximizes this linear form over the current polytope.
    // Hence the supporting half-space of their convex hull preserves every
    // other current vertex, not merely the neighbors.  The LP is a numerical
    // fallback when active-row sums cancel or the representation is noisy.
    if (separating_cut_lp(points, normal, rhs) &&
        verify_neighbor_separator(points, normal, rhs, &margin))
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"conservative_neighbor_cut\",\"mode\":\"separator_lp\",\"neighbors\":%zu,\"margin\":%.12g}\n",
                    vertex_row, margin);
        return true;
    }
    if (active_valid)
    {
        for (size_t j = 0; j < n; ++j)
            gsl_vector_set(normal, j, active_normal[j]);
        *rhs = active_rhs;
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"conservative_neighbor_cut\",\"mode\":\"verified_fallback\",\"source\":\"%s\",\"neighbors\":%zu,\"margin\":%.12g}\n",
                    backup_mode, vertex_row, active_margin);
        return true;
    }
    if (!lp_first && affine_neighbor_normal(points, normal) &&
        verify_neighbor_separator(points, normal, rhs, &margin))
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"conservative_neighbor_cut\",\"mode\":\"affine_fallback\",\"neighbors\":%zu,\"margin\":%.12g}\n",
                    vertex_row, margin);
        return true;
    }
    return false;
}

// Experimental true-cone parent cut from Cone_Decomposition, Section 3.3.
// Unlike a neighbor-hull cut, this may remove adjacent vertices, but only
// after independently checking local dominance and every forward radial LP.
// It is used only on explicit opt-in after the conservative cut fails.
bool positive_ray_parent_cut(gsl_matrix *points, LpModel *lp_y,
                             const gsl_matrix *Q, const gsl_vector *dy,
                             gsl_vector *normal, double *rhs,
                             ConeRayGeometry *geometry)
{
    if (geometry) *geometry = ConeRayGeometry{};
    if (!points || !lp_y || !Q || !dy || !normal || !rhs ||
        points->size1 < 2 || points->size2 != Q->size1 ||
        Q->size2 != static_cast<size_t>(get_Ncolumns(lp_y)) ||
        dy->size != Q->size2 || normal->size != points->size2)
        return false;
    const int nx = static_cast<int>(points->size2);
    const int ny = get_Ncolumns(lp_y);
    const int m = get_Nrows(lp_y);
    const int neighbors = static_cast<int>(points->size1) - 1;
    const double infinity = get_infinite(lp_y);
    auto rejected = [&](const char *reason, int ray) {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"positive_ray_rejected\",\"reason\":\"%s\",\"ray\":%d,\"neighbors\":%d}\n",
                    reason, ray, neighbors);
        return false;
    };
    std::vector<REAL> objective(ny + 1, 0.0), row(ny + 2, 0.0);
    std::vector<double> rays(neighbors * nx, 0.0), gauges(neighbors, 0.0),
        ray_lengths(neighbors, 0.0), neighbor_gaps(neighbors, 0.0);
    gsl_vector_const_view vertex = gsl_matrix_const_row(points, neighbors);
    auto set_y_objective = [&](const double *x) {
        objective[0] = 0.0;
        for (int j = 0; j < ny; ++j)
        {
            double value = gsl_vector_get(dy, j);
            for (int i = 0; i < nx; ++i)
                value += gsl_matrix_get(Q, i, j) * x[i];
            objective[j + 1] = value;
        }
        set_obj_fn(lp_y, objective.data());
        return solve(lp_y) == OPTIMAL;
    };
    std::vector<double> x(nx);
    for (int i = 0; i < nx; ++i) x[i] = gsl_vector_get(&vertex.vector, i);
    if (!set_y_objective(x.data())) return rejected("apex_lp", -1);
    const double alpha = get_objective(lp_y);
    if (!std::isfinite(alpha)) return rejected("apex_value", -1);
    const double value_tolerance = 1e-8 * (1.0 + std::fabs(alpha));
    for (int k = 0; k < neighbors; ++k)
    {
        double length_sq = 0.0;
        for (int i = 0; i < nx; ++i)
        {
            const double direction = gsl_matrix_get(points, k, i) -
                                     gsl_vector_get(&vertex.vector, i);
            rays[k * nx + i] = direction;
            length_sq += direction * direction;
        }
        const double length = std::sqrt(length_sq);
        ray_lengths[k] = length;
        if (!std::isfinite(length) || length <= 1e-8)
            return rejected("short_ray", k);
        for (int i = 0; i < nx; ++i) rays[k * nx + i] /= length;
        for (int i = 0; i < nx; ++i) x[i] = gsl_matrix_get(points, k, i);
        if (!set_y_objective(x.data())) return rejected("neighbor_lp", k);
        neighbor_gaps[k] = get_objective(lp_y) - alpha;
        if (get_objective(lp_y) < alpha - value_tolerance)
            return rejected("neighbor_improves", k);
    }
    for (int i = 0; i < nx; ++i) x[i] = gsl_vector_get(&vertex.vector, i);

    // Charnes--Cooper: z=t*y, -delta(z)=1, minimize
    // (Q'v+dy)'z-alpha*t over the homogenized Y polytope.
    std::unique_ptr<LpModel, decltype(&delete_lp)> radial(
        make_lp(m + 2 * ny + 1, ny + 1), &delete_lp);
    set_verbose(radial.get(), kLpNeutral);
    set_minim(radial.get());
    for (int j = 1; j <= ny; ++j)
    {
        set_lowbo(radial.get(), j, -infinity);
        set_upbo(radial.get(), j, infinity);
    }
    for (int k = 0; k < m; ++k)
    {
        std::fill(row.begin(), row.end(), 0.0);
        for (int j = 0; j < ny; ++j) row[j + 1] = lp_y->rows[k][j];
        row[ny + 1] = -lp_y->rhs[k];
        set_row(radial.get(), k + 1, row.data());
        set_constr_type(radial.get(), k + 1, lp_y->types[k]);
        set_rh(radial.get(), k + 1, 0.0);
    }
    for (int j = 0; j < ny; ++j)
    {
        std::fill(row.begin(), row.end(), 0.0);
        if (lp_y->upper[j] < infinity / 2)
        {
            row[j + 1] = 1.0;
            row[ny + 1] = -lp_y->upper[j];
        }
        set_row(radial.get(), m + 2 * j + 1, row.data());
        std::fill(row.begin(), row.end(), 0.0);
        if (lp_y->lower[j] > -infinity / 2)
        {
            row[j + 1] = -1.0;
            row[ny + 1] = lp_y->lower[j];
        }
        set_row(radial.get(), m + 2 * j + 2, row.data());
    }
    std::vector<double> slopes(ny, 0.0), cc_objective(ny + 2, 0.0);
    for (int j = 0; j < ny; ++j)
    {
        cc_objective[j + 1] = gsl_vector_get(dy, j);
        for (int i = 0; i < nx; ++i)
            cc_objective[j + 1] += gsl_matrix_get(Q, i, j) * x[i];
    }
    cc_objective[ny + 1] = -alpha;
    set_obj_fn(radial.get(), cc_objective.data());
    int finite = 0;
    for (int k = 0; k < neighbors; ++k)
    {
        std::fill(objective.begin(), objective.end(), 0.0);
        for (int j = 0; j < ny; ++j)
        {
            double slope = 0.0;
            for (int i = 0; i < nx; ++i)
                slope += gsl_matrix_get(Q, i, j) * rays[k * nx + i];
            slopes[j] = slope;
            objective[j + 1] = slope;
        }
        set_obj_fn(lp_y, objective.data());
        if (solve(lp_y) != OPTIMAL) return rejected("slope_lp", k);
        if (get_objective(lp_y) >= -1e-9) continue; // recession ray
        std::fill(row.begin(), row.end(), 0.0);
        for (int j = 0; j < ny; ++j) row[j + 1] = -slopes[j];
        set_row(radial.get(), m + 2 * ny + 1, row.data());
        set_constr_type(radial.get(), m + 2 * ny + 1, kLpEq);
        set_rh(radial.get(), m + 2 * ny + 1, 1.0);
        if (solve(radial.get()) != OPTIMAL) return rejected("radial_lp", k);
        const double tau = get_objective(radial.get());
        if (!std::isfinite(tau) || tau <= 1e-7)
        {
            if (std::getenv("DECIDEIT_TRACE"))
                fprintf(stderr, "TRACE {\"event\":\"positive_ray_radius\",\"ray\":%d,\"tau\":%.12g,\"min_slope\":%.12g,\"alpha\":%.12g,\"neighbor_distance\":%.12g,\"neighbor_gap\":%.12g}\n",
                        k, tau, get_objective(lp_y), alpha,
                        ray_lengths[k], neighbor_gaps[k]);
            return rejected("small_radius", k);
        }
        // A slightly shorter, directly checked radial length avoids using
        // an optimistic LP endpoint as a lower bound on the true gauge.
        const double safe_tau = tau - 1e-7 * (1.0 + tau);
        if (safe_tau <= 1e-7) return rejected("small_safe_radius", k);
        for (int i = 0; i < nx; ++i)
            x[i] = gsl_vector_get(&vertex.vector, i) +
                   safe_tau * rays[k * nx + i];
        if (!set_y_objective(x.data())) return rejected("contact_lp", k);
        if (get_objective(lp_y) < alpha - value_tolerance)
            return rejected("contact_outside", k);
        gauges[k] = 1.0 / safe_tau;
        ++finite;
    }
    if (!finite) return rejected("all_recession", -1);

    // Everything needed by a child-local polar construction has now been
    // independently certified.  Preserve it even when the single parent
    // separator below is too steep to use numerically: avoiding that common
    // separator is precisely the purpose of the cone-decomposition pilot.
    if (geometry)
    {
        geometry->valid = true;
        geometry->dimension = nx;
        geometry->ray_count = neighbors;
        geometry->alpha = alpha;
        geometry->apex.resize(nx);
        for (int i = 0; i < nx; ++i)
            geometry->apex[i] = gsl_vector_get(&vertex.vector, i);
        geometry->rays = rays;
        geometry->gauges = gauges;
    }

    // Any feasible ray majorant is valid.  Minimize its L1 norm rather than
    // the paper's directional depth objective here: the latter can select
    // extremely large, nearly cancelling coefficients and leave an unusably
    // tiny normalized geometric gap at the apex.
    std::unique_ptr<LpModel, decltype(&delete_lp)> separator(
        make_lp(neighbors, 2 * nx), &delete_lp);
    set_verbose(separator.get(), kLpNeutral);
    set_minim(separator.get());
    std::vector<REAL> selection(2 * nx + 1, 0.0),
        ray_row(2 * nx + 1, 0.0);
    for (int i = 1; i <= 2 * nx; ++i) selection[i] = 1.0;
    for (int k = 0; k < neighbors; ++k)
    {
        for (int i = 0; i < nx; ++i)
        {
            const double value = rays[k * nx + i];
            ray_row[i + 1] = value;
            ray_row[nx + i + 1] = -value;
        }
        set_row(separator.get(), k + 1, ray_row.data());
        set_constr_type(separator.get(), k + 1, kLpGe);
        set_rh(separator.get(), k + 1, gauges[k]);
    }
    set_obj_fn(separator.get(), selection.data());
    if (solve(separator.get()) != OPTIMAL)
        return rejected("separator_lp", -1);
    std::vector<double> a_solution(2 * nx, 0.0);
    get_variables(separator.get(), a_solution.data());
    std::vector<double> a(nx);
    for (int i = 0; i < nx; ++i)
        a[i] = a_solution[i] - a_solution[nx + i];
    const char *repair_setting = std::getenv("DECIDEIT_POSITIVE_RAY_REPAIR");
    const bool repair_enabled = repair_setting &&
        std::strcmp(repair_setting, "1") == 0;
    if (repair_enabled)
    {
        // This is an opt-in numerical experiment.  Every candidate must
        // satisfy the full true-ray system after reconstruction; neither an
        // LP status nor a small negative residual is treated as a certificate.
        std::vector<std::vector<double>> candidates(1, a);
        std::vector<double> iterate = a;
        for (int pass = 0; pass < 2; ++pass)
        {
            double largest = 0.0;
            for (double coefficient : iterate)
                largest = (std::max)(largest, std::fabs(coefficient));
            if (!std::isfinite(largest) || largest <= 0.0) break;
            for (int i = 0; i < nx; ++i)
            {
                const double weight = (std::max)(0.02,
                    std::fabs(iterate[i]) / largest);
                selection[i + 1] = weight;
                selection[nx + i + 1] = weight;
            }
            set_obj_fn(separator.get(), selection.data());
            if (solve(separator.get()) != OPTIMAL) break;
            get_variables(separator.get(), a_solution.data());
            for (int i = 0; i < nx; ++i)
                iterate[i] = a_solution[i] - a_solution[nx + i];
            candidates.push_back(iterate);
        }

        auto ray_dot = [&](const std::vector<double> &candidate, int k) {
            long double value = 0.0L;
            for (int i = 0; i < nx; ++i)
                value += static_cast<long double>(rays[k * nx + i]) *
                         candidate[i];
            return value;
        };
        bool needs_interior = false;
        for (const auto &candidate : candidates)
            for (int k = 0; k < neighbors; ++k)
                if (ray_dot(candidate, k) < 0.0L)
                    needs_interior = true;

        std::vector<double> interior;
        if (needs_interior)
        {
            // A pointed full-dimensional true cone has a strictly positive
            // dual direction.  It can repair a tiny negative recession-ray
            // residual without relaxing any required ray inequality.
            for (int k = 0; k < neighbors; ++k)
                set_rh(separator.get(), k + 1, 1.0);
            for (int i = 1; i <= 2 * nx; ++i) selection[i] = 1.0;
            set_obj_fn(separator.get(), selection.data());
            if (solve(separator.get()) == OPTIMAL)
            {
                get_variables(separator.get(), a_solution.data());
                interior.resize(nx);
                for (int i = 0; i < nx; ++i)
                    interior[i] = a_solution[i] - a_solution[nx + i];
                for (int k = 0; k < neighbors; ++k)
                    if (!(ray_dot(interior, k) >= 0.5L))
                    {
                        interior.clear();
                        break;
                    }
            }
        }

        double best_norm = std::numeric_limits<double>::infinity();
        int best_index = -1;
        double best_inflation = 0.0;
        for (size_t index = 0; index < candidates.size(); ++index)
        {
            std::vector<double> candidate = candidates[index];
            if (!interior.empty())
            {
                long double inflation = 0.0L;
                for (int k = 0; k < neighbors; ++k)
                {
                    const long double target = gauges[k] > 0.0 ?
                        static_cast<long double>(gauges[k]) +
                        1e-7L * (1.0L + gauges[k]) : 0.0L;
                    const long double deficit = target - ray_dot(candidate, k);
                    if (deficit > 0.0L)
                        inflation = (std::max)(inflation,
                            deficit / ray_dot(interior, k));
                }
                if (!std::isfinite(static_cast<double>(inflation))) continue;
                for (int i = 0; i < nx; ++i)
                    candidate[i] += static_cast<double>(inflation) * interior[i];
            }
            long double scale = 1.0L;
            bool feasible = true;
            for (int k = 0; k < neighbors; ++k)
            {
                const long double lhs = ray_dot(candidate, k);
                if (gauges[k] == 0.0)
                {
                    if (lhs < 0.0L) feasible = false;
                }
                else if (lhs <= 0.0L)
                    feasible = false;
                else
                    scale = (std::max)(scale,
                        (static_cast<long double>(gauges[k]) +
                         1e-7L * (1.0L + gauges[k])) / lhs);
            }
            if (!feasible || !std::isfinite(static_cast<double>(scale))) continue;
            for (double &coefficient : candidate)
                coefficient *= static_cast<double>(scale);
            double norm_sq = 0.0;
            for (double coefficient : candidate)
                norm_sq += coefficient * coefficient;
            const double candidate_norm = std::sqrt(norm_sq);
            if (!std::isfinite(candidate_norm)) continue;
            for (int k = 0; k < neighbors; ++k)
            {
                const long double lhs = ray_dot(candidate, k);
                const long double target = gauges[k] > 0.0 ?
                    static_cast<long double>(gauges[k]) +
                    1e-8L * (1.0L + gauges[k]) : 0.0L;
                if (lhs < target)
                {
                    feasible = false;
                    break;
                }
            }
            if (feasible && candidate_norm < best_norm)
            {
                best_norm = candidate_norm;
                best_index = static_cast<int>(index);
                best_inflation = static_cast<double>(scale);
                a.swap(candidate);
            }
        }
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"positive_ray_repair\",\"candidates\":%zu,\"selected\":%d,\"normal_norm\":%.12g,\"scale\":%.12g,\"interior_lp\":%s}\n",
                    candidates.size(), best_index,
                    best_index >= 0 ? best_norm : 0.0, best_inflation,
                    needs_interior ? "true" : "false");
        if (best_index < 0) return rejected("repair_failed", -1);
    }
    double norm_sq = 0.0;
    for (int i = 0; i < nx; ++i) norm_sq += a[i] * a[i];
    const double norm = std::sqrt(norm_sq);
    if (!std::isfinite(norm) || norm <= 1e-12 ||
        norm >= 1.0 / (8 * DBLP_EPSILON))
    {
        if (std::getenv("DECIDEIT_TRACE"))
        {
            double largest_gauge = 0.0;
            for (double gauge : gauges)
                if (gauge > largest_gauge) largest_gauge = gauge;
            double pair_norm_lower_bound = largest_gauge;
            for (int first = 0; first < neighbors; ++first)
                for (int second = first + 1; second < neighbors; ++second)
                {
                    double sum_sq = 0.0;
                    for (int i = 0; i < nx; ++i)
                    {
                        const double sum = rays[first * nx + i] +
                                           rays[second * nx + i];
                        sum_sq += sum * sum;
                    }
                    if (sum_sq > 1e-24)
                    {
                        const double bound = (gauges[first] + gauges[second]) /
                                             std::sqrt(sum_sq);
                        if (bound > pair_norm_lower_bound)
                            pair_norm_lower_bound = bound;
                    }
                }
            fprintf(stderr, "TRACE {\"event\":\"positive_ray_scale\",\"normal_norm\":%.12g,\"max_gauge\":%.12g,\"pair_norm_lower_bound\":%.12g}\n",
                    norm, largest_gauge, pair_norm_lower_bound);
        }
        return rejected("separator_scale", -1);
    }
    for (int k = 0; k < neighbors; ++k)
    {
        double lhs = 0.0;
        for (int i = 0; i < nx; ++i) lhs += rays[k * nx + i] * a[i];
        if (lhs < gauges[k] - 1e-7 * (1.0 + gauges[k]))
            return rejected("separator_residual", k);
    }
    double vertex_value = 0.0;
    for (int i = 0; i < nx; ++i)
    {
        gsl_vector_set(normal, i, -a[i] / norm);
        vertex_value += a[i] * gsl_vector_get(&vertex.vector, i);
    }
    *rhs = (-1.0 - vertex_value) / norm;
    if (std::getenv("DECIDEIT_TRACE"))
        fprintf(stderr, "TRACE {\"event\":\"positive_ray_parent_cut\",\"rays\":%d,\"finite\":%d,\"normal_norm\":%.12g}\n",
                neighbors, finite, norm);
    return true;
}

void free_cp_normal(gsl_vector *temp_y,
                    gsl_vector *temp_v,
                    gsl_vector *X_extension,
                    gsl_vector *coef)
{
    managed::vector_free(temp_y);
    managed::vector_free(temp_v);
    managed::vector_free(X_extension);
    managed::vector_free(coef);
}

/*
    save_cp_normal() saves the results derived in cutting_points_normal().
*/

void save_cp_normal(gsl_vector *X,
                    gsl_vector *edge,
                    REAL alpha,
                    gsl_vector *x_solution,
                    int indicator)
{
    gsl_vector_memcpy(X, edge);

    switch (indicator)
    {
    case -1: // a negative cutting point
        gsl_vector_scale(X, -alpha);
        break;

    case 0: // exceptional case when alpha is infinite along both extensions
        gsl_vector_scale(X, alpha);
        break;

    case 1: // a positive cutting point
        gsl_vector_scale(X, alpha);
        break;

    default:
        break;
    }

    gsl_vector_sub(X, x_solution);
    gsl_vector_scale(X, -1);
}

/*
    cutting_point_normal() derives the cutting points along a given edge.
    it is designed for a NON-DEGENERATE x, thus the word "normal".

    inputs
        x_solution: the current local solution in X;
        edge: one edge emanating from "x_solution";
        instance: original data;
        alpha: the current best objective value;
        lp_Y: established LP for Y, i.e., P-base;
        cuts_indicator: the derived cutting point being along the positive (1), negative (-1), or abnormal (0) extension;
        which_edge: the index of the given edge, used to fill in "cuts_indicator".

    outputs
        1: a negative cutting point;
        0: a positive cutting point;
*/

int cutting_point(
    gsl_vector *x_solution,
    gsl_vector *edge,
    double alpha,
    LpModel *lp_Y,
    gsl_vector *X,
    struct dataset *dataset_in)
{
    const int n_variable_x = dataset_in->xl->size, n_variable_y = get_Ncolumns(lp_Y);
    int counter;
    REAL lambda[100], slope[100], obj, sum_1;
    bool loopFlag;

    gsl_vector *X_extension = managed::vector_calloc(n_variable_x);
    gsl_vector *temp_v = managed::vector_calloc(n_variable_x);   // to save the results as a temporary vector
    gsl_vector *temp_y = managed::vector_calloc(n_variable_y);   // to save the results of y's coefficients as a temporary vector
    gsl_vector *coef = managed::vector_calloc(1 + n_variable_y); // to save the results of y's coefficients as a temporary vector

    gsl_vector_view dest_v;

    lambda[0] = 1e8; // initializing lambda

    // X_p = X - lambda[0] * e, the initial cutting point along e's positive extension
    gsl_vector_memcpy(temp_v, edge);
    gsl_vector_scale(temp_v, lambda[0]);
    gsl_vector_sub(temp_v, x_solution);
    gsl_vector_scale(temp_v, -1);
    gsl_vector_memcpy(X_extension, temp_v);

    loopFlag = TRUE;

    // construct the coefficient of y, i.e., Q' * X_p + dy, send it to lp_Y, and solve the LP
    // gsl_vector_memcpy(temp_y, instance->dy);

    gsl_vector_set_zero(temp_y);
    gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, X_extension, 1, temp_y);
    gsl_vector_add(temp_y, dataset_in->dy);
    dest_v = gsl_vector_subvector(coef, 1, n_variable_y);
    gsl_vector_memcpy(&dest_v.vector, temp_y);
    set_obj_fn(lp_Y, coef->data);
    solve(lp_Y);

    // derive the true objective value, fval + cx' * X_p
    // gsl_blas_ddot(X_p, instance->cx, &obj);

    obj = get_objective(lp_Y);

    /* derive the cutting point along e's positive extension */
    if (obj < alpha - DBLP_OBJ_EPSILON) // obj is really less than alpha
    {
        counter = 0;

        while (fabs(obj - alpha) > DBLP_OBJ_EPSILON) // obj and alpha are not equal yet
        {
            if(counter>=98) { free_cp_normal(temp_y,temp_v,X_extension,coef); return SOLUTION_LP; }
            get_variables(lp_Y, temp_y->data);

            // derive the slope -(e' * Q * y_solution)
            gsl_blas_dgemv(CblasNoTrans, 1, dataset_in->Q, temp_y, 0, temp_v);
            gsl_blas_ddot(edge, temp_v, &sum_1);

            if (fabs(sum_1) <= DBLP_sum_EPSILON)
            {
                loopFlag = FALSE;
                break;
            }

            slope[counter] = -sum_1;
            lambda[counter + 1] = (alpha - obj + slope[counter] * lambda[counter]) / slope[counter]; // the new lambda

            if (fabs(lambda[counter + 1]) > 1e8)
            {
                loopFlag = FALSE;
                break;
            }

            counter = counter + 1;

            // derive the new point, X_p = X - lambda[counter] * e
            gsl_vector_memcpy(X_extension, edge);
            gsl_vector_scale(X_extension, lambda[counter]);
            gsl_vector_sub(X_extension, x_solution);
            gsl_vector_scale(X_extension, -1);

            // derive the new coefficient of y, Q' * X_p + dy, send it to lp_Y, and solve the LP
            gsl_vector_set_zero(temp_y);
            gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, X_extension, 1, temp_y);
            gsl_vector_add(temp_y, dataset_in->dy);
            gsl_vector_memcpy(&dest_v.vector, temp_y);
            set_obj_fn(lp_Y, coef->data);
            solve(lp_Y);

            // derive the true objective value
            obj = get_objective(lp_Y);
        }

        // given alpha, the cutting point is X_p = X - lambda * e along the positive extension of the edge
        if (!loopFlag)
        {
            save_cp_normal(X, edge, 1, x_solution, 0);
        }
        else
        {
            save_cp_normal(X, edge, lambda[counter], x_solution, 1);
        }

        free_cp_normal(temp_y, temp_v, X_extension, coef);

        return (0);
    }

    gsl_vector_memcpy(X_extension, edge);
    gsl_vector_scale(X_extension, lambda[0]);
    gsl_vector_add(X_extension, x_solution); // X_n = X + lambda[0] * e, the initial cutting point along e's negative extension

    loopFlag = TRUE;

    // construct the coefficient of y,  Q' * X_n, send it to lp_Y, and solve the LP
    gsl_vector_set_zero(temp_y);
    gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, X_extension, 1, temp_y);
    gsl_vector_add(temp_y, dataset_in->dy);
    gsl_vector_scale(temp_y, -1); // it becomes a max LP along e's negative extension
    gsl_vector_memcpy(&dest_v.vector, temp_y);
    set_obj_fn(lp_Y, coef->data);
    solve(lp_Y);

    // derive the true objective value
    obj = -get_objective(lp_Y);

    /* derive the cutting point along e's negative extension */
    if (obj < alpha - DBLP_OBJ_EPSILON) // obj is really less than alpha
    {
        counter = 0;

        while (fabs(obj - alpha) > DBLP_OBJ_EPSILON) // obj and alpha are not equal yet
        {
            // construct y_solution
            if(counter>=98) { free_cp_normal(temp_y,temp_v,X_extension,coef); return SOLUTION_LP; }
            get_variables(lp_Y, temp_y->data);

            // derive -(e' * Q * y_solution)
            gsl_blas_dgemv(CblasNoTrans, 1, dataset_in->Q, temp_y, 0, temp_v);
            gsl_blas_ddot(edge, temp_v, &sum_1);

            if (fabs(sum_1) <= DBLP_sum_EPSILON)
            {
                loopFlag = FALSE;
                break;
            }

            slope[counter] = -sum_1;
            lambda[counter + 1] = (alpha - obj - slope[counter] * lambda[counter]) / (-slope[counter]); // derive next lambda

            if (fabs(lambda[counter + 1]) > 1e8)
            {
                loopFlag = FALSE;
                break;
            }

            counter = counter + 1;

            // derive the new X_n = X + lambda(counter) * e
            gsl_vector_memcpy(X_extension, edge);
            gsl_vector_scale(X_extension, lambda[counter]);
            gsl_vector_add(X_extension, x_solution);

            // derive the new coefficient of y, Q' * X_n + dy, send it to lp_Y, and solve the LP
            gsl_vector_set_zero(temp_y);
            gsl_blas_dgemv(CblasTrans, 1, dataset_in->Q, X_extension, 1, temp_y);
            gsl_vector_add(temp_y, dataset_in->dy);
            gsl_vector_scale(temp_y, -1); // it now becomes a max lp
            gsl_vector_memcpy(&dest_v.vector, temp_y);
            set_obj_fn(lp_Y, coef->data);
            solve(lp_Y);

            // derive the true objective value
            obj = -get_objective(lp_Y);
        }

        // given alpha, the cutting point is X_n = X + lambda * e along the negative extension of the edge
        if (!loopFlag)
        {
            save_cp_normal(X, edge, 1, x_solution, 0);
            free_cp_normal(temp_y, temp_v, X_extension, coef);

            return (0);
        }
        else
        {
            save_cp_normal(X, edge, lambda[counter], x_solution, -1);
            free_cp_normal(temp_y, temp_v, X_extension, coef);

            return (1);
        }
    }

    // set lambda to 1e8 along the positive extension when lambdas are infinities along both extensions
    save_cp_normal(X, edge, 1e8, x_solution, 0);
    free_cp_normal(temp_y, temp_v, X_extension, coef);

    return (0);
}

void free_global_phase(gsl_vector *temp_array_x,
                       gsl_matrix *cutting_points,
                       gsl_matrix *points,
                       gsl_vector *a,
                       gsl_vector *b,
                       gsl_vector *c,
                       gsl_vector *temp_rhs,
                       gsl_vector *normal_vector,
                       gsl_permutation *p,
                       gsl_matrix *all_Ax,
                       gsl_vector *all_bx)
{
    managed::vector_free(temp_array_x);
    managed::matrix_free(cutting_points);
    managed::matrix_free(points);
    managed::vector_free(a);
    managed::vector_free(b);
    managed::vector_free(c);
    managed::vector_free(temp_rhs);
    managed::vector_free(normal_vector);
    managed::permutation_free(p);
    managed::matrix_free(all_Ax);
    managed::vector_free(all_bx);
}

/*
    global_phase() generates cutting points along each edge at x

    inputs
        instance: the structure containing the original data;
        local_amc: the structure to save a PGM;
        alpha: the current best objective value;
        lp_Y: the constructed LPs for P-base;
        dataset_in: the structure containing the actual linear constraints, bounds and variables in X;
        noc: the # of linear constraints in X.

    outputs
        1       successful derivations along all edges
        -3      LP subproblem failed
        -4      GSL LU_decomp/det/invert/solve failed
        -7      ALL_NEGATIVE    this is the only good sign for an early termination
*/

// Experimental implementation of the Complete-Boundary idea.  It resumes
// reverse search from the exact pivoted dictionary saved for every first-
// layer neighbor.  An incomplete parent search invalidates the entire B set.
// Diagnostic mode only reports geometry; cut mode returns a candidate only
// after a scale-independent margin check.
static void complete_boundary_state_diagnostic(
    const local_solution *first, const dataset *data,
    gsl_vector *candidate_normal, double *candidate_rhs,
    bool *candidate_found)
{
    *candidate_found = false;
    const char *setting =
        std::getenv("DECIDEIT_COMPLETE_BOUNDARY_DIAGNOSTIC");
    const char *cut_setting =
        std::getenv("DECIDEIT_COMPLETE_BOUNDARY_CUT");
    const bool enabled = (setting && std::strcmp(setting, "1") == 0) ||
        (cut_setting && std::strcmp(cut_setting, "1") == 0);
    if (!enabled) return;
    long long budget_ms = 5000;
    if (const char *configured_budget =
            std::getenv("DECIDEIT_COMPLETE_BOUNDARY_BUDGET_MS"))
    {
        char *end = nullptr;
        errno = 0;
        const long parsed = std::strtol(configured_budget, &end, 10);
        if (errno || end == configured_budget || *end ||
            parsed < 0 || parsed > 60000)
            throw std::invalid_argument(
                "DECIDEIT_COMPLETE_BOUNDARY_BUDGET_MS must be 0..60000");
        budget_ms = parsed;
    }
    const int dimension = static_cast<int>(first->x_solution->size);
    const int parent_count = first->n_adj_vertices;
    const auto started = std::chrono::steady_clock::now();
    const auto elapsed_ms = [&]() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - started).count();
    };
    size_t state_bytes = 0;
    for (const auto &state : first->neighbor_states)
        state_bytes += state.tableau.size() * sizeof(double) +
            (state.basic_index.size() + state.nonbasic_index.size()) *
                sizeof(int) + state.point.size() * sizeof(double);
    std::vector<std::vector<double>> boundary;
    const auto emit = [&](const char *status, int processed, double margin) {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"complete_boundary_state\",\"status\":\"%s\",\"dimension\":%d,\"first_layer\":%d,\"states\":%zu,\"processed\":%d,\"boundary\":%zu,\"margin_l1\":%.12g,\"state_mib\":%.6g,\"elapsed_ms\":%lld}\n",
                    status, dimension, parent_count,
                    first->neighbor_states.size(), processed,
                    boundary.size(), margin,
                    state_bytes / (1024.0 * 1024.0),
                    static_cast<long long>(elapsed_ms()));
    };
    if (dimension < 1 || dimension > 20 || parent_count < 1 ||
        parent_count > 100)
    {
        emit("size_limit", 0, 0.0);
        return;
    }
    if (first->neighbor_states.size() !=
        static_cast<size_t>(parent_count))
    {
        emit("missing_state", 0, 0.0);
        return;
    }
    managed::Scope scope;
    const auto same_point = [&](const std::vector<double> &point,
                                const gsl_vector *other) {
        for (int j = 0; j < dimension; ++j)
            if (std::fabs(point[j] - gsl_vector_get(other, j)) > 1e-7)
                return false;
        return true;
    };
    for (int parent_index = 0; parent_index < parent_count; ++parent_index)
    {
        if (elapsed_ms() >= budget_ms)
        {
            emit("time_limit", parent_index, 0.0);
            return;
        }
        const NeighborDictionaryState &state =
            first->neighbor_states[parent_index];
        if (state.rows < 1 || state.columns < 2 ||
            state.tableau.size() !=
                static_cast<size_t>(state.rows) * state.columns ||
            state.basic_index.size() != static_cast<size_t>(state.rows + 1) ||
            state.point.size() != static_cast<size_t>(dimension))
        {
            emit("invalid_state", parent_index, 0.0);
            return;
        }
        gsl_matrix *tableau = managed::matrix_alloc(state.rows, state.columns);
        gsl_vector_int *basic =
            managed::vector_int_alloc(state.basic_index.size());
        gsl_vector_int *nonbasic =
            managed::vector_int_alloc(state.nonbasic_index.size());
        for (int i = 0; i < state.rows; ++i)
            for (int j = 0; j < state.columns; ++j)
                gsl_matrix_set(tableau, i, j,
                    state.tableau[static_cast<size_t>(i) * state.columns + j]);
        for (size_t i = 0; i < state.basic_index.size(); ++i)
            gsl_vector_int_set(basic, i, state.basic_index[i]);
        for (size_t i = 0; i < state.nonbasic_index.size(); ++i)
            gsl_vector_int_set(nonbasic, i, state.nonbasic_index[i]);

        // A state must reproduce the coordinate that it accompanies before
        // it is trusted as a second-layer root.
        std::vector<double> reconstructed(dimension, 0.0);
        for (int j = 0; j < dimension; ++j)
            reconstructed[j] = gsl_vector_get(data->xl, j);
        for (int row = 0; row < state.rows; ++row)
        {
            const int variable = state.basic_index[row + 1];
            if (variable >= 1 && variable <= dimension)
                reconstructed[variable - 1] = gsl_matrix_get(tableau, row, 0) +
                    gsl_vector_get(data->xl, variable - 1);
        }
        gsl_vector_const_view parent =
            gsl_matrix_const_column(first->adj_vertices, parent_index);
        std::vector<double> stored_point(dimension);
        for (int j = 0; j < dimension; ++j)
            stored_point[j] = state.point[j] + gsl_vector_get(data->xl, j);
        if (!same_point(stored_point, &parent.vector) ||
            !same_point(reconstructed, &parent.vector))
        {
            double largest_error = 0.0, largest_order_error = 0.0;
            int largest_index = -1;
            for (int j = 0; j < dimension; ++j)
            {
                const double dictionary_error =
                    std::fabs(reconstructed[j] -
                              gsl_vector_get(&parent.vector, j));
                if (dictionary_error > largest_error)
                {
                    largest_error = dictionary_error;
                    largest_index = j;
                }
                largest_order_error = (std::max)(largest_order_error,
                    std::fabs(stored_point[j] -
                              gsl_vector_get(&parent.vector, j)));
            }
            if (std::getenv("DECIDEIT_TRACE"))
                fprintf(stderr, "TRACE {\"event\":\"complete_boundary_state_mismatch\",\"parent\":%d,\"coordinate\":%d,\"dictionary_error\":%.12g,\"order_error\":%.12g,\"dictionary_value\":%.12g,\"stored_value\":%.12g,\"lower\":%.12g}\n",
                        parent_index, largest_index, largest_error,
                        largest_order_error, reconstructed[largest_index],
                        stored_point[largest_index],
                        gsl_vector_get(data->xl, largest_index));
            emit("state_mismatch", parent_index, 0.0);
            return;
        }
        gsl_matrix *second = managed::matrix_calloc(
            dimension, adjacency_capacity(dimension));
        const int neighbor_count = reverse_search(
            tableau, basic, nonbasic, second);
        if (neighbor_count < 1)
        {
            emit(neighbor_count == NEIGHBOR_TIME_LIMIT ?
                     "parent_time_limit" : "parent_incomplete",
                 parent_index, 0.0);
            return;
        }
        for (int neighbor_index = 0;
             neighbor_index < neighbor_count; ++neighbor_index)
        {
            std::vector<double> point(dimension);
            for (int j = 0; j < dimension; ++j)
                point[j] = gsl_matrix_get(second, j, neighbor_index) +
                           gsl_vector_get(data->xl, j);
            if (same_point(point, first->x_solution)) continue;
            bool first_layer = false;
            for (int i = 0; i < parent_count && !first_layer; ++i)
            {
                gsl_vector_const_view candidate =
                    gsl_matrix_const_column(first->adj_vertices, i);
                first_layer = same_point(point, &candidate.vector);
            }
            if (first_layer) continue;
            bool duplicate = false;
            for (const auto &old : boundary)
            {
                duplicate = true;
                for (int j = 0; j < dimension; ++j)
                    if (std::fabs(old[j] - point[j]) > 1e-7)
                    {
                        duplicate = false;
                        break;
                    }
                if (duplicate) break;
            }
            if (!duplicate) boundary.push_back(std::move(point));
            if (boundary.size() > 10000)
            {
                emit("boundary_limit", parent_index + 1, 0.0);
                return;
            }
        }
        managed::matrix_free(second);
        managed::matrix_free(tableau);
        managed::vector_int_free(basic);
        managed::vector_int_free(nonbasic);
    }
    if (boundary.empty())
    {
        emit("empty_boundary", parent_count, 0.0);
        return;
    }

    // Variables are [a+,a-,beta,delta].  Retain every B point and maximize
    // one common deletion margin for the center and all first-layer points.
    const int columns = 2 * dimension + 2;
    const int deleted_count = parent_count + 1;
    std::unique_ptr<LpModel, decltype(&delete_lp)> separator(
        make_lp(static_cast<int>(boundary.size()) + deleted_count + 1,
                columns), &delete_lp);
    separator->backend = "highs";
    set_verbose(separator.get(), kLpNeutral);
    set_maxim(separator.get());
    const double infinity = get_infinite(separator.get());
    set_lowbo(separator.get(), 2 * dimension + 1, -infinity);
    set_upbo(separator.get(), 2 * dimension + 1, infinity);
    std::vector<double> row(columns + 1, 0.0),
        objective(columns + 1, 0.0);
    int row_number = 0;
    for (const auto &point : boundary)
    {
        std::fill(row.begin(), row.end(), 0.0);
        for (int j = 0; j < dimension; ++j)
        {
            row[j + 1] = -point[j];
            row[dimension + j + 1] = point[j];
        }
        row[2 * dimension + 1] = 1.0;
        set_row(separator.get(), ++row_number, row.data());
    }
    for (int point_index = 0; point_index < deleted_count; ++point_index)
    {
        std::fill(row.begin(), row.end(), 0.0);
        for (int j = 0; j < dimension; ++j)
        {
            const double value = point_index == parent_count ?
                gsl_vector_get(first->x_solution, j) :
                gsl_matrix_get(first->adj_vertices, j, point_index);
            row[j + 1] = value;
            row[dimension + j + 1] = -value;
        }
        row[2 * dimension + 1] = -1.0;
        row[2 * dimension + 2] = 1.0;
        set_row(separator.get(), ++row_number, row.data());
    }
    std::fill(row.begin(), row.end(), 0.0);
    for (int j = 1; j <= 2 * dimension; ++j) row[j] = 1.0;
    set_row(separator.get(), ++row_number, row.data());
    set_rh(separator.get(), row_number, 1.0);
    objective[2 * dimension + 2] = 1.0;
    set_obj_fn(separator.get(), objective.data());
    if (solve(separator.get()) != OPTIMAL)
    {
        emit("separator_lp", parent_count, 0.0);
        return;
    }
    const double common_margin = get_objective(separator.get());
    if (common_margin > 1e-10)
    {
        std::vector<double> solution(columns, 0.0);
        get_variables(separator.get(), solution.data());
        double norm_squared = 0.0;
        for (int j = 0; j < dimension; ++j)
        {
            // The separator LP uses a*x >= beta for retained boundary
            // points.  X cuts use normal*x <= rhs, hence normal=-a.
            const double value = solution[dimension + j] - solution[j];
            gsl_vector_set(candidate_normal, j, value);
            norm_squared += value * value;
        }
        const double norm = std::sqrt(norm_squared);
        double highest_boundary = -std::numeric_limits<double>::infinity();
        double lowest_deleted = std::numeric_limits<double>::infinity();
        if (norm > 1e-12)
        {
            for (const auto &point : boundary)
            {
                double value = 0.0;
                for (int j = 0; j < dimension; ++j)
                    value += gsl_vector_get(candidate_normal, j) * point[j];
                highest_boundary = (std::max)(highest_boundary, value);
            }
            for (int point_index = 0; point_index < deleted_count;
                 ++point_index)
            {
                double value = 0.0;
                for (int j = 0; j < dimension; ++j)
                {
                    const double coordinate = point_index == parent_count ?
                        gsl_vector_get(first->x_solution, j) :
                        gsl_matrix_get(first->adj_vertices, j, point_index);
                    value += gsl_vector_get(candidate_normal, j) * coordinate;
                }
                lowest_deleted = (std::min)(lowest_deleted, value);
            }
        }
        const double geometric_margin = norm > 1e-12 ?
            (lowest_deleted - highest_boundary) / norm : 0.0;
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"complete_boundary_common\",\"dimension\":%d,\"first_layer\":%d,\"boundary\":%zu,\"margin_l1\":%.12g,\"geometric_margin\":%.12g}\n",
                    dimension, parent_count, boundary.size(), common_margin,
                    geometric_margin);
        if (geometric_margin > 8 * DBLP_EPSILON)
        {
            *candidate_rhs = highest_boundary + 2 * DBLP_EPSILON * norm;
            *candidate_found = lowest_deleted > *candidate_rhs +
                4 * DBLP_EPSILON * norm;
            emit(*candidate_found ? "complete_common_strong" :
                                   "complete_common_guard_failed",
                 parent_count, geometric_margin);
            return;
        }
    }

    // Complete-Boundary fallback: retaining B is the safety requirement;
    // C may lie on either side.  Test whether the center alone is outside
    // conv(B), and report how many first-layer points the resulting cut would
    // also delete.  Cut mode may return this weaker candidate, again only
    // after a scale-independent margin check.
    gsl_matrix *points = managed::matrix_calloc(
        boundary.size() + 1, dimension);
    gsl_vector *normal = managed::vector_calloc(dimension);
    for (size_t i = 0; i < boundary.size(); ++i)
        for (int j = 0; j < dimension; ++j)
            gsl_matrix_set(points, i, j, boundary[i][j]);
    gsl_matrix_set_row(points, boundary.size(), first->x_solution);
    double rhs = 0.0;
    if (!separating_cut_lp(points, normal, &rhs))
    {
        emit("complete_no_separator", parent_count, 0.0);
        return;
    }
    double center_value = 0.0;
    gsl_blas_ddot(normal, first->x_solution, &center_value);
    const double norm = gsl_blas_dnrm2(normal);
    const double geometric_margin = norm > 0.0 ?
        (center_value - rhs) / norm : 0.0;
    int deleted_first = 0;
    for (int i = 0; i < parent_count; ++i)
    {
        gsl_vector_const_view point =
            gsl_matrix_const_column(first->adj_vertices, i);
        double value = 0.0;
        gsl_blas_ddot(normal, &point.vector, &value);
        if (value > rhs + 1e-10) ++deleted_first;
    }
    if (std::getenv("DECIDEIT_TRACE"))
        fprintf(stderr, "TRACE {\"event\":\"complete_boundary_center\",\"dimension\":%d,\"first_layer\":%d,\"boundary\":%zu,\"geometric_margin\":%.12g,\"deleted_first\":%d,\"state_mib\":%.6g,\"elapsed_ms\":%lld}\n",
                dimension, parent_count, boundary.size(), geometric_margin,
                deleted_first, state_bytes / (1024.0 * 1024.0),
                static_cast<long long>(elapsed_ms()));
    const bool strong = geometric_margin > 8 * DBLP_EPSILON;
    emit(strong ? "complete_center_strong" : "complete_center_weak",
         parent_count, geometric_margin);
    if (strong)
    {
        gsl_vector_memcpy(candidate_normal, normal);
        *candidate_rhs = rhs + 2 * DBLP_EPSILON * norm;
        *candidate_found = center_value > *candidate_rhs +
            4 * DBLP_EPSILON * norm;
    }
}

static int append_neighbor_hull_cut(struct local_solution *local_amc,
                                    LpModel *lp_X_std,
                                    LpModel *lp_Y,
                                    struct dataset *dataset_in,
                                    int n_constraint)
{
    const int n_variable = get_Ncolumns(lp_X_std);
    const int n_edge = local_amc->n_adj_vertices;
    if (n_edge < 1) return CUT_SEPARATION_FAILED;
    gsl_matrix *neighbor_points = managed::matrix_calloc(n_edge + 1, n_variable);
    gsl_vector *neighbor_normal = managed::vector_calloc(n_variable);
    gsl_vector *cut_row = managed::vector_calloc(n_variable + 1);
    for (int edge = 0; edge < n_edge; ++edge)
    {
        gsl_vector_const_view neighbor =
            gsl_matrix_const_column(local_amc->adj_vertices, edge);
        gsl_matrix_set_row(neighbor_points, edge, &neighbor.vector);
    }
    gsl_matrix_set_row(neighbor_points, n_edge, local_amc->x_solution);
    double cut_rhs = 0.0;
    bool found = conservative_neighbor_cut(neighbor_points, dataset_in->valid_A,
                                           n_constraint, dataset_in->xl,
                                           dataset_in->xu, neighbor_normal,
                                           &cut_rhs);
    const char *positive_mode = std::getenv("DECIDEIT_DEGENERATE_CUT");
    const bool cone_pilot = std::getenv("DECIDEIT_CONE_PILOT_DIR");
    if (!found && ((positive_mode &&
        std::strcmp(positive_mode, "positive_on_failure") == 0) ||
        cone_pilot))
        found = positive_ray_parent_cut(neighbor_points, lp_Y,
                                        dataset_in->Q, dataset_in->dy,
                                        neighbor_normal, &cut_rhs,
                                        &local_amc->cone_geometry);
    if (!found)
    {
        bool boundary_candidate = false;
        complete_boundary_state_diagnostic(
            local_amc, dataset_in, neighbor_normal, &cut_rhs,
            &boundary_candidate);
        const char *boundary_cut =
            std::getenv("DECIDEIT_COMPLETE_BOUNDARY_CUT");
        found = boundary_candidate && boundary_cut &&
            std::strcmp(boundary_cut, "1") == 0;
    }
    if (!found)
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"cut_separation_failed\",\"dimension\":%d,\"neighbors\":%d}\n",
                    n_variable, n_edge);
        managed::matrix_free(neighbor_points);
        managed::vector_free(neighbor_normal);
        managed::vector_free(cut_row);
        return CUT_SEPARATION_FAILED;
    }
    gsl_vector_view coefficients = gsl_vector_subvector(cut_row, 1, n_variable);
    gsl_vector_memcpy(&coefficients.vector, neighbor_normal);
    double shift = 0.0;
    gsl_blas_ddot(neighbor_normal, dataset_in->xl, &shift);
    add_constraint(lp_X_std, cut_row->data, kLpLe, cut_rhs - shift);
    gsl_matrix_set_row(dataset_in->valid_A, n_constraint, cut_row);
    gsl_matrix_set(dataset_in->valid_A, n_constraint, 0, cut_rhs);
    managed::matrix_free(neighbor_points);
    managed::vector_free(neighbor_normal);
    managed::vector_free(cut_row);
    return 1;
}

int global_phase(
    struct local_solution *local_amc,
    double alpha,
    LpModel *lp_X_std,
    LpModel *lp_Y,
    struct dataset *dataset_in,
    int noc)
{
    int i, j, n_negative = 0, exitFlag, which_constraint, counter;
    const int n_variable = get_Ncolumns(lp_X_std), n_edge = local_amc->n_adj_vertices, n_constraint = noc;
    REAL criterion, result, cut_rhs = 1.0;

    const char *cut_mode = std::getenv("DECIDEIT_CONSERVATIVE_CUT");
    const bool use_neighbor_hull = local_amc->degenerate_vertex &&
        !(cut_mode && std::strcmp(cut_mode, "legacy") == 0);
    if (use_neighbor_hull)
    {
        // The degenerate case needs no facet-combination search.  Work with
        // the actual adjacent vertices, not extrapolated cutting points:
        // a strict separator from their convex hull retains every other
        // vertex of the bounded current X polytope.
        return append_neighbor_hull_cut(local_amc, lp_X_std, lp_Y,
                                        dataset_in, n_constraint);
    }

    gsl_matrix *cutting_points = managed::matrix_calloc(n_variable, n_edge);
    gsl_vector *a = managed::vector_calloc(n_variable), *b = managed::vector_calloc(n_variable), *c = managed::vector_calloc(n_variable);
    gsl_vector *temp_array_x = managed::vector_calloc(n_variable + 1);

    /* construct the edges emanating from "x_solution" and cutting points in an ordinary manner */
    gsl_vector_memcpy(c, local_amc->x_solution);

    for (i = 0; i < n_edge; i++)
    {
        gsl_matrix_get_col(a, local_amc->adj_vertices, i);
        gsl_vector_sub(a, c);
        gsl_vector_scale(a, -1);

        exitFlag = cutting_point(local_amc->x_solution, a, alpha, lp_Y, b, dataset_in);

        if(exitFlag<0) {
            managed::matrix_free(cutting_points);managed::vector_free(a);managed::vector_free(b);managed::vector_free(c);managed::vector_free(temp_array_x);
            if (exitFlag == SOLUTION_LP &&
                !(cut_mode && std::strcmp(cut_mode, "legacy") == 0))
                return append_neighbor_hull_cut(local_amc, lp_X_std, lp_Y,
                                                dataset_in, n_constraint);
            return exitFlag;
        }
        if (exitFlag == 1)
        {
            n_negative++;
        }

        // show_cutting_values(b, dataset_in, lp_Y, exitFlag, i + 1);

        if(std::getenv("DECIDEIT_GEOMETRY_LOG")) gsl_vector_fprintf(stdout,b,"% .12g");

        gsl_matrix_set_col(cutting_points, i, b);
    }

    if (n_negative == n_edge)
    {
        gsl_vector_set_zero(temp_array_x);
        add_constraint(lp_X_std, temp_array_x->data, kLpEq, 0);

        managed::matrix_free(cutting_points);
        managed::vector_free(a);
        managed::vector_free(b);
        managed::vector_free(c);
        managed::vector_free(temp_array_x);

        return (ALL_NEGATIVE);
    }

    // show_gsl_matrix_columns(cutting_points, n_edge);

    const int total_n_constraint = n_constraint + 2 * n_variable;
    bool CutFlag;
    int ceiling, ceiling_num_constraint, n_qualified;
    std::vector<int> index(n_variable - 1);

    gsl_vector *temp_rhs = managed::vector_calloc(n_variable);
    gsl_vector_set_all(temp_rhs, 1);
    gsl_vector *normal_vector = managed::vector_calloc(n_variable);
    gsl_permutation *p = managed::permutation_calloc(n_variable);
    gsl_vector_int *indices_qualified = NULL, *indices_new = NULL;
    gsl_matrix *all_Ax = managed::matrix_calloc(total_n_constraint, n_variable);
    gsl_vector *all_bx = managed::vector_calloc(total_n_constraint);
    gsl_matrix *points = managed::matrix_calloc(n_edge + 1, n_variable);

    gsl_matrix_view sub_dest, sub_src;
    gsl_vector_view destination, source;
    gsl_vector_int_view v_src;

    if (n_edge > n_variable)
    {
        sub_dest = gsl_matrix_submatrix(points, 0, 0, n_edge, n_variable);
        sub_src = gsl_matrix_submatrix(cutting_points, 0, 0, n_variable, n_edge);

        gsl_matrix_transpose_memcpy(&sub_dest.matrix, &sub_src.matrix);
        gsl_matrix_set_row(points, n_edge, local_amc->x_solution);

        // show_gsl_matrix(points);

        sub_dest = gsl_matrix_submatrix(all_Ax, 0, 0, n_constraint, n_variable);
        sub_src = gsl_matrix_submatrix(dataset_in->valid_A, 0, 1, n_constraint, n_variable);
        gsl_matrix_memcpy(&sub_dest.matrix, &sub_src.matrix);

        destination = gsl_vector_subvector(all_bx, 0, n_constraint);
        source = gsl_matrix_subcolumn(dataset_in->valid_A, 0, 0, n_constraint);
        gsl_vector_memcpy(&destination.vector, &source.vector);

        /* incorporate lower and upper bounds for each effective variable */
        counter = 0;
        ceiling = n_constraint + n_variable;

        for (i = n_constraint; i < ceiling; i++)
        {
            gsl_vector_set(all_bx, i, -gsl_vector_get(dataset_in->xl, counter));
            gsl_matrix_set(all_Ax, i, counter++, -1);
        }

        counter = 0;

        for (i = ceiling; i < total_n_constraint; i++)
        {
            gsl_vector_set(all_bx, i, gsl_vector_get(dataset_in->xu, counter));
            gsl_matrix_set(all_Ax, i, counter++, 1);
        }

        which_constraint = -1;
        indices_qualified = managed::vector_int_calloc(1 + points->size1);
        indices_new = managed::vector_int_calloc(1 + n_variable);
        ceiling_num_constraint = all_bx->size;

        while (which_constraint < ceiling_num_constraint)
        {
            which_constraint++;

            if (active_constraint(points, all_Ax, all_bx, &which_constraint, indices_qualified))
            {
                n_qualified = gsl_vector_int_get(indices_qualified, 0);

                if (n_qualified > n_variable)
                {
                    ceiling = nchoosek(n_qualified - 1, n_variable - 1);

                    if (ceiling <= 1e4)
                    {
                        for (i = 1; i <= ceiling; i++)
                        {
                            nchoosek_enum(n_qualified - 1, n_variable - 1, i, index);

                            for (j = 1; j <= n_variable - 1; j++)
                            {
                                gsl_vector_int_set(indices_new, j, gsl_vector_int_get(indices_qualified, index[j - 1]));
                            }

                            gsl_vector_int_set(indices_new, n_variable, points->size1 - 1);
                            gsl_vector_int_set(indices_new, 0, n_variable);

                            CutFlag = conservativeCut_New(points, indices_new, normal_vector);

                            if (CutFlag == 1)
                            {
                                if (is_good_vector(normal_vector->data, normal_vector->size))
                                {
                                    break;
                                }
                                else
                                {
                                    CutFlag = 0;
                                }
                            }
                        }

                        if (CutFlag == 1)
                        {
                            break;
                        }
                    }
                }
                else // n_qualified == n_variable
                {
                    v_src = gsl_vector_int_subvector(indices_qualified, 0, 1 + n_variable);
                    gsl_vector_int_memcpy(indices_new, &v_src.vector);

                    CutFlag = conservativeCut_New(points, indices_new, normal_vector);

                    if (CutFlag == 1)
                    {
                        break;
                    }
                }
            }
        }

        managed::vector_int_free(indices_qualified);
        managed::vector_int_free(indices_new);

        if (which_constraint == ceiling_num_constraint)
        {
            if (!separating_cut_lp(points, normal_vector, &cut_rhs))
            {
                if (std::getenv("DECIDEIT_TRACE"))
                    fprintf(stderr, "TRACE {\"event\":\"cut_separation_failed\",\"dimension\":%d,\"neighbors\":%d}\n",
                            n_variable, n_edge);
                free_global_phase(temp_array_x, cutting_points, points, a, b, c, temp_rhs, normal_vector, p, all_Ax, all_bx);
                return (CUT_SEPARATION_FAILED);
            }
        }
    }
    else // n_edge == n_variable
    {
        /* derive the normal vector and introduce the cut into dataset_in->valid_A as well as lp_X_std */
        gsl_matrix_transpose(cutting_points);

        // if (!condition_number(cutting_points, 1e12)) // the matrix is not invertible
        // {
        //     free_global_phase(temp_array_x, cutting_points, points, a, b, c, temp_rhs, normal_vector, p, all_Ax, all_bx);
        //     return (GSL_LU);
        // }

        if (svd_solve_safe(cutting_points, temp_rhs, normal_vector, 1e-12) != 0)
        {
            free_global_phase(temp_array_x, cutting_points, points, a, b, c, temp_rhs, normal_vector, p, all_Ax, all_bx);
            return (GSL_LU);
        }

        // if (gsl_linalg_LU_decomp(cutting_points, p, &s) != GSL_SUCCESS)
        // {
        //     free_global_phase(temp_array_x, cutting_points, points, a, b, c, temp_rhs, normal_vector, p, all_Ax, all_bx);
        //     return (GSL_LU);
        // }

        // if (gsl_linalg_LU_solve(cutting_points, p, temp_rhs, normal_vector) != GSL_SUCCESS)
        // {
        //     free_global_phase(temp_array_x, cutting_points, points, a, b, c, temp_rhs, normal_vector, p, all_Ax, all_bx);
        //     return (GSL_LU);
        // }
    }

    resolve_zero_vector(normal_vector);

    if(std::getenv("DECIDEIT_GEOMETRY_LOG")) gsl_vector_fprintf(stdout,normal_vector,"% .12g");


    destination = gsl_vector_subvector(temp_array_x, 1, n_variable);
    gsl_vector_memcpy(&destination.vector, normal_vector);
    gsl_blas_ddot(normal_vector, local_amc->x_solution, &criterion);

    if (n_edge > n_variable && std::getenv("DECIDEIT_TRACE"))
    {
        int violations = 0;
        double worst = 0.0;
        double highest_oriented_neighbor = -std::numeric_limits<double>::infinity();
        const double orientation = criterion < cut_rhs ? -1.0 : 1.0;
        for (int edge = 0; edge < n_edge; ++edge)
        {
            gsl_vector_const_view neighbor =
                gsl_matrix_const_column(local_amc->adj_vertices, edge);
            double value = 0.0;
            gsl_blas_ddot(normal_vector, &neighbor.vector, &value);
            if (orientation * value > highest_oriented_neighbor)
                highest_oriented_neighbor = orientation * value;
            const double violation = criterion < cut_rhs
                ? cut_rhs - value : value - cut_rhs;
            if (violation > 1e-8)
            {
                ++violations;
                if (violation > worst) worst = violation;
            }
        }
        fprintf(stderr, "TRACE {\"event\":\"legacy_cut_neighbor_check\",\"neighbors\":%d,\"violations\":%d,\"worst\":%.12g,\"separation_margin\":%.12g}\n",
                n_edge, violations, worst,
                orientation * criterion - highest_oriented_neighbor);
    }

    if (criterion < cut_rhs)
    {
        gsl_vector_scale(temp_array_x, -1);

        gsl_blas_ddot(&destination.vector, dataset_in->xl, &result);
        result = -cut_rhs - result;
        add_constraint(lp_X_std, temp_array_x->data, kLpLe, result);

        gsl_matrix_set_row(dataset_in->valid_A, n_constraint, temp_array_x);
        gsl_matrix_set(dataset_in->valid_A, n_constraint, 0, -cut_rhs);
    }
    else
    {
        gsl_blas_ddot(&destination.vector, dataset_in->xl, &result);
        result = cut_rhs - result;
        add_constraint(lp_X_std, temp_array_x->data, kLpLe, result);

        gsl_matrix_set_row(dataset_in->valid_A, n_constraint, temp_array_x);
        gsl_matrix_set(dataset_in->valid_A, n_constraint, 0, cut_rhs);
    }

    free_global_phase(temp_array_x, cutting_points, points, a, b, c, temp_rhs, normal_vector, p, all_Ax, all_bx);

    return (1);
}
