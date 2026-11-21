#include "dblp.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

double value(const gsl_vector *normal, const std::vector<double> &point)
{
    double result = 0.0;
    for (size_t j = 0; j < point.size(); ++j)
        result += gsl_vector_get(normal, j) * point[j];
    return result;
}

void check_cut(const std::vector<std::vector<double>> &neighbors,
               const std::vector<double> &vertex, bool with_bounds)
{
    const size_t n = vertex.size();
    gsl_matrix *points = gsl_matrix_calloc(neighbors.size() + 1, n);
    gsl_vector *normal = gsl_vector_calloc(n);
    gsl_vector *lower = gsl_vector_calloc(n);
    gsl_vector *upper = gsl_vector_calloc(n);
    gsl_vector_set_all(upper, 1.0);
    for (size_t i = 0; i < neighbors.size(); ++i)
        for (size_t j = 0; j < n; ++j)
            gsl_matrix_set(points, i, j, neighbors[i][j]);
    for (size_t j = 0; j < n; ++j)
        gsl_matrix_set(points, neighbors.size(), j, vertex[j]);
    double rhs = 0.0;
    const bool found = conservative_neighbor_cut(
        points, nullptr, 0, with_bounds ? lower : nullptr,
        with_bounds ? upper : nullptr, normal, &rhs);
    require(found, "failed to separate a vertex from its neighbor hull");
    require(value(normal, vertex) > rhs + 1e-9,
            "cut did not strictly remove the current vertex");
    for (const auto &neighbor : neighbors)
        require(value(normal, neighbor) <= rhs + 1e-9,
                "cut removed an adjacent vertex");
    gsl_matrix_free(points);
    gsl_vector_free(normal);
    gsl_vector_free(lower);
    gsl_vector_free(upper);
}

} // namespace

int main()
{
    gsl_set_error_handler_off();
    setenv("DECIDEIT_NEIGHBOR_AFFINE_HULL", "1", 1);
    setenv("DECIDEIT_NEIGHBOR_GUARD_FACTOR", "0", 1);
    check_cut({{1, 0}, {0, 1}}, {0, 0}, true);
    check_cut({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
              {0, 0, 0}, true);
    check_cut({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0.5}},
              {0, 0, 0}, true); // full-rank neighbor hull: affine branch falls back
    unsetenv("DECIDEIT_NEIGHBOR_AFFINE_HULL");
    unsetenv("DECIDEIT_NEIGHBOR_GUARD_FACTOR");
    check_cut({{1, 0}, {0, 1}}, {0, 0}, true);
    reset_lp_backend_stats();
    check_cut({{1, 0}, {0, 1}}, {0, 0}, false);
    LpBackendStats stats = last_lp_backend_stats();
    require(stats.highs_solves + stats.pivotcopy_solves == 0,
            "centroid candidate unexpectedly needed a separator LP");
    reset_lp_backend_stats();
    check_cut({{1, 0}, {-1, 1}}, {0, 0}, false);
    stats = last_lp_backend_stats();
    require(stats.highs_solves + stats.pivotcopy_solves > 0,
            "failed to use separator LP after centroid rejection");
    check_cut({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
              {0, 0, 0}, true);
    check_cut({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
              {0, 0, 0}, false);

    gsl_matrix *general_points = gsl_matrix_calloc(3, 2);
    gsl_matrix_set(general_points, 0, 0, 2.0 / 3.0);
    gsl_matrix_set(general_points, 0, 1, 1.0 / 3.0);
    gsl_matrix_set(general_points, 1, 0, 1.0 / 3.0);
    gsl_matrix_set(general_points, 1, 1, 2.0 / 3.0);
    gsl_matrix *general_rows = gsl_matrix_calloc(2, 3);
    gsl_matrix_set(general_rows, 0, 1, 1.0);
    gsl_matrix_set(general_rows, 0, 2, -2.0);
    gsl_matrix_set(general_rows, 1, 1, -2.0);
    gsl_matrix_set(general_rows, 1, 2, 1.0);
    gsl_vector *general_normal = gsl_vector_calloc(2);
    double general_rhs = 0.0;
    require(conservative_neighbor_cut(general_points, general_rows, 2,
                                      nullptr, nullptr, general_normal,
                                      &general_rhs),
            "general linear active rows did not yield a hull cut");
    require(value(general_normal, {0, 0}) > general_rhs + 1e-9 &&
            value(general_normal, {2.0 / 3.0, 1.0 / 3.0}) <= general_rhs &&
            value(general_normal, {1.0 / 3.0, 2.0 / 3.0}) <= general_rhs,
            "general-linear hull cut lost a neighbor");
    gsl_matrix_free(general_points);
    gsl_matrix_free(general_rows);
    gsl_vector_free(general_normal);

    // If the benchmark lies in the convex hull, the requested strict cut
    // does not exist; reporting failure is safer than inventing a cut.
    gsl_matrix *impossible = gsl_matrix_calloc(3, 1);
    gsl_vector *normal = gsl_vector_calloc(1);
    gsl_matrix_set(impossible, 0, 0, -1.0);
    gsl_matrix_set(impossible, 1, 0, 1.0);
    double rhs = 0.0;
    require(!conservative_neighbor_cut(impossible, nullptr, 0, nullptr,
                                       nullptr, normal, &rhs),
            "accepted a separator despite hull containment");
    gsl_matrix_free(impossible);
    gsl_vector_free(normal);

    gsl_matrix *tiny = gsl_matrix_calloc(3, 2);
    normal = gsl_vector_calloc(2);
    gsl_matrix_set(tiny, 0, 0, 1e-8);
    gsl_matrix_set(tiny, 1, 1, 1e-8);
    require(!conservative_neighbor_cut(tiny, nullptr, 0, nullptr,
                                       nullptr, normal, &rhs),
            "accepted a cut below the tableau's numerical resolution");
    gsl_matrix_free(tiny);
    gsl_vector_free(normal);

    // g(x)=min(x1+x2, 5-x1-x2) on [0,4]^2.  At the origin both adjacent
    // vertices are non-improving, yet the opposite vertex improves.  The
    // true-cone parent cut may remove those neighbors and must keep every
    // strict improver; the conservative hull cut cannot do that.
    gsl_matrix *radial_points = gsl_matrix_calloc(3, 2);
    gsl_matrix_set(radial_points, 0, 0, 4.0);
    gsl_matrix_set(radial_points, 1, 1, 4.0);
    gsl_matrix *Q = gsl_matrix_calloc(2, 2);
    for (int i = 0; i < 2; ++i)
    {
        gsl_matrix_set(Q, i, 0, 1.0);
        gsl_matrix_set(Q, i, 1, -1.0);
    }
    gsl_vector *dy = gsl_vector_calloc(2);
    gsl_vector_set(dy, 1, 5.0);
    normal = gsl_vector_calloc(2);
    LpModel *lp_y = make_lp(1, 2);
    const double simplex_row[] = {0.0, 1.0, 1.0};
    set_row(lp_y, 1, simplex_row);
    set_constr_type(lp_y, 1, kLpEq);
    set_rh(lp_y, 1, 1.0);
    set_upbo(lp_y, 1, 1.0);
    set_upbo(lp_y, 2, 1.0);
    ConeRayGeometry cone_geometry;
    require(positive_ray_parent_cut(radial_points, lp_y, Q, dy,
                                    normal, &rhs, &cone_geometry),
            "failed to construct a forward positive-ray parent cut");
    require(cone_geometry.valid && cone_geometry.dimension == 2 &&
            cone_geometry.ray_count == 2 &&
            cone_geometry.apex.size() == 2 &&
            cone_geometry.rays.size() == 4 &&
            cone_geometry.gauges.size() == 2 &&
            cone_geometry.gauges[0] > 0.0 &&
            cone_geometry.gauges[1] > 0.0,
            "positive-ray cut did not retain certified cone geometry");
    setenv("DECIDEIT_POSITIVE_RAY_REPAIR", "1", 1);
    require(positive_ray_parent_cut(radial_points, lp_y, Q, dy,
                                    normal, &rhs),
            "repaired forward positive-ray parent cut failed");
    unsetenv("DECIDEIT_POSITIVE_RAY_REPAIR");
    require(value(normal, {0, 0}) > rhs &&
            value(normal, {4, 0}) > rhs &&
            value(normal, {0, 4}) > rhs &&
            value(normal, {4, 4}) <= rhs,
            "parent cut did not strengthen the conservative neighbor cut");
    for (int i = 0; i <= 16; ++i)
        for (int j = 0; j <= 16; ++j)
        {
            const double x1 = i / 4.0, x2 = j / 4.0;
            const double sum = x1 + x2;
            if (std::min(sum, 5.0 - sum) < -1e-7)
                require(value(normal, {x1, x2}) <= rhs + 1e-7,
                        "parent cut removed a strict improver");
        }
    gsl_vector_set(dy, 1, 3.0);
    require(!positive_ray_parent_cut(radial_points, lp_y, Q, dy,
                                     normal, &rhs),
            "parent cut ignored an improving adjacent vertex");
    delete_lp(lp_y);
    gsl_vector_free(normal);
    gsl_vector_free(dy);
    gsl_matrix_free(Q);
    gsl_matrix_free(radial_points);

    std::cout << "Conservative neighbor-hull cut passed\n";
}
