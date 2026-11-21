#include "dblp.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

void chooseCandidateLex(gsl_vector_int *, gsl_matrix *,
                        const std::vector<int> &, double,
                        gsl_matrix_int *);

namespace {

void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

struct SearchResult
{
    int status;
    std::vector<std::vector<double>> neighbors;
    std::vector<NeighborDictionaryState> states;
};

SearchResult search(const std::vector<std::vector<double>> &rows,
                    size_t neighbor_capacity = 8,
                    const char *mode = nullptr,
                    bool capture_states = false)
{
    if (mode) setenv("DECIDEIT_REVERSE_SEARCH", mode, 1);
    else unsetenv("DECIDEIT_REVERSE_SEARCH");
    const size_t m = rows.size();
    const size_t dimension = rows[0].size() - 1;
    gsl_matrix *tableau = gsl_matrix_calloc(m, 1 + dimension + m);
    gsl_vector_int *basic = gsl_vector_int_calloc(m + 1);
    gsl_vector_int *nonbasic = gsl_vector_int_calloc(dimension + 1);
    gsl_matrix *neighbors = gsl_matrix_calloc(dimension, neighbor_capacity);

    gsl_vector_int_set(basic, 0, static_cast<int>(m));
    gsl_vector_int_set(nonbasic, 0, static_cast<int>(dimension));
    for (size_t j = 0; j < dimension; ++j)
        gsl_vector_int_set(nonbasic, j + 1, static_cast<int>(j + 1));
    for (size_t i = 0; i < m; ++i)
    {
        require(rows[i].size() == dimension + 1, "inconsistent row width");
        gsl_matrix_set(tableau, i, 0, rows[i][0]);
        for (size_t j = 0; j < dimension; ++j)
            gsl_matrix_set(tableau, i, j + 1, rows[i][j + 1]);
        gsl_matrix_set(tableau, i, dimension + 1 + i, 1.0);
        gsl_vector_int_set(basic, i + 1, static_cast<int>(dimension + 1 + i));
    }

    SearchResult result;
    result.status = reverse_search(tableau, basic, nonbasic, neighbors,
                                   capture_states ? &result.states : nullptr);
    if (result.status > 0)
    {
        for (int j = 0; j < result.status; ++j)
        {
            std::vector<double> point(dimension);
            for (size_t i = 0; i < dimension; ++i)
                point[i] = gsl_matrix_get(neighbors, i, j);
            result.neighbors.push_back(point);
        }
    }

    gsl_matrix_free(tableau);
    gsl_vector_int_free(basic);
    gsl_vector_int_free(nonbasic);
    gsl_matrix_free(neighbors);
    return result;
}

std::vector<double> state_point(const NeighborDictionaryState &state,
                                size_t dimension)
{
    std::vector<double> point(dimension, 0.0);
    for (int row = 0; row < state.rows; ++row)
    {
        const int variable = state.basic_index[row + 1];
        if (variable >= 1 && variable <= static_cast<int>(dimension))
            point[variable - 1] =
                state.tableau[static_cast<size_t>(row) * state.columns];
    }
    return point;
}

SearchResult resume(const NeighborDictionaryState &state, size_t dimension)
{
    gsl_matrix *tableau = gsl_matrix_alloc(state.rows, state.columns);
    gsl_vector_int *basic = gsl_vector_int_alloc(state.basic_index.size());
    gsl_vector_int *nonbasic =
        gsl_vector_int_alloc(state.nonbasic_index.size());
    gsl_matrix *neighbors = gsl_matrix_calloc(
        dimension, adjacency_capacity(dimension));
    for (int i = 0; i < state.rows; ++i)
        for (int j = 0; j < state.columns; ++j)
            gsl_matrix_set(tableau, i, j,
                state.tableau[static_cast<size_t>(i) * state.columns + j]);
    for (size_t i = 0; i < state.basic_index.size(); ++i)
        gsl_vector_int_set(basic, i, state.basic_index[i]);
    for (size_t i = 0; i < state.nonbasic_index.size(); ++i)
        gsl_vector_int_set(nonbasic, i, state.nonbasic_index[i]);
    SearchResult result;
    result.status = reverse_search(tableau, basic, nonbasic, neighbors);
    if (result.status > 0)
        for (int column = 0; column < result.status; ++column)
        {
            std::vector<double> point(dimension);
            for (size_t row = 0; row < dimension; ++row)
                point[row] = gsl_matrix_get(neighbors, row, column);
            result.neighbors.push_back(std::move(point));
        }
    gsl_matrix_free(tableau);
    gsl_vector_int_free(basic);
    gsl_vector_int_free(nonbasic);
    gsl_matrix_free(neighbors);
    return result;
}

bool contains(const SearchResult &result, double x, double y)
{
    for (const auto &point : result.neighbors)
        if (std::fabs(point[0] - x) < 1e-8 &&
            std::fabs(point[1] - y) < 1e-8)
            return true;
    return false;
}

bool contains_point(const SearchResult &result,
                    const std::vector<double> &expected)
{
    for (const auto &point : result.neighbors)
    {
        double squared_distance = 0.0;
        for (size_t i = 0; i < expected.size(); ++i)
        {
            const double difference = point[i] - expected[i];
            squared_distance += difference * difference;
        }
        if (squared_distance < 1e-16) return true;
    }
    return false;
}

double determinant3(const std::array<std::array<double, 3>, 3> &a)
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

// An independent vertex/edge oracle for small, full-dimensional 3D
// polytopes.  Vertices are intersections of three supporting planes;
// two distinct vertices are adjacent iff their common active rows span
// a rank-two face.  This uses no tableau pivots or Geue traversal.
std::vector<std::vector<double>> oracle_neighbors3(
    const std::vector<std::vector<double>> &rows)
{
    std::vector<std::vector<double>> all = rows;
    for (int j = 0; j < 3; ++j)
    {
        std::vector<double> lower(4, 0.0);
        lower[j + 1] = -1.0;
        all.push_back(lower);
    }
    std::vector<std::vector<double>> vertices;
    for (size_t i = 0; i < all.size(); ++i)
        for (size_t j = i + 1; j < all.size(); ++j)
            for (size_t k = j + 1; k < all.size(); ++k)
            {
                const size_t ids[3] = {i, j, k};
                std::array<std::array<double, 3>, 3> a{};
                std::array<double, 3> b{};
                for (int r = 0; r < 3; ++r)
                {
                    b[r] = all[ids[r]][0];
                    for (int c = 0; c < 3; ++c)
                        a[r][c] = all[ids[r]][c + 1];
                }
                const double det = determinant3(a);
                if (std::fabs(det) < 1e-10) continue;
                std::vector<double> point(3);
                for (int c = 0; c < 3; ++c)
                {
                    auto replaced = a;
                    for (int r = 0; r < 3; ++r) replaced[r][c] = b[r];
                    point[c] = determinant3(replaced) / det;
                }
                bool feasible = true;
                for (const auto &row : all)
                    if (row[1] * point[0] + row[2] * point[1] +
                        row[3] * point[2] > row[0] + 1e-8)
                        feasible = false;
                if (!feasible) continue;
                bool duplicate = false;
                for (const auto &old : vertices)
                    if (std::fabs(old[0] - point[0]) < 1e-8 &&
                        std::fabs(old[1] - point[1]) < 1e-8 &&
                        std::fabs(old[2] - point[2]) < 1e-8)
                        duplicate = true;
                if (!duplicate) vertices.push_back(point);
            }
    std::vector<std::vector<double>> neighbors;
    for (const auto &vertex : vertices)
    {
        if (std::fabs(vertex[0]) + std::fabs(vertex[1]) +
            std::fabs(vertex[2]) < 1e-8) continue;
        std::vector<std::array<double, 3>> common;
        for (const auto &row : all)
        {
            if (std::fabs(row[0]) > 1e-8) continue;
            const double at_vertex = row[1] * vertex[0] +
                                     row[2] * vertex[1] +
                                     row[3] * vertex[2];
            if (std::fabs(at_vertex - row[0]) < 1e-8)
                common.push_back({{row[1], row[2], row[3]}});
        }
        bool rank_two = false;
        for (size_t i = 0; i < common.size(); ++i)
            for (size_t j = i + 1; j < common.size(); ++j)
            {
                const auto &u = common[i], &v = common[j];
                const double cross_sq =
                    std::pow(u[1] * v[2] - u[2] * v[1], 2) +
                    std::pow(u[2] * v[0] - u[0] * v[2], 2) +
                    std::pow(u[0] * v[1] - u[1] * v[0], 2);
                if (cross_sq > 1e-12) rank_two = true;
            }
        if (rank_two) neighbors.push_back(vertex);
    }
    return neighbors;
}

} // namespace

int main()
{
    gsl_set_error_handler_off();

    // Both zero rows have the same first TNP ratio.  The second column of
    // the negative basis must choose row 1, not the first tableau row.
    gsl_matrix *tie_table = gsl_matrix_calloc(2, 4);
    gsl_vector_int *tie_nonbasic = gsl_vector_int_calloc(2);
    gsl_matrix_int *tie_candidates = gsl_matrix_int_calloc(3, 2);
    gsl_matrix_int_set_all(tie_candidates, -1);
    gsl_vector_int_set(tie_nonbasic, 0, 1);
    gsl_vector_int_set(tie_nonbasic, 1, 1);
    gsl_matrix_set(tie_table, 0, 1, 1);
    gsl_matrix_set(tie_table, 1, 1, 1);
    gsl_matrix_set(tie_table, 0, 2, -1);
    gsl_matrix_set(tie_table, 1, 2, -1);
    gsl_matrix_set(tie_table, 0, 3, -2);
    gsl_matrix_set(tie_table, 1, 3, -1);
    chooseCandidateLex(tie_nonbasic, tie_table, {2, 3}, -1.0,
                       tie_candidates);
    require(gsl_matrix_int_get(tie_candidates, 0, 1) == 1,
            "TNP did not use the second perturbation column to break a tie");
    gsl_matrix_int_free(tie_candidates);
    gsl_vector_int_free(tie_nonbasic);
    gsl_matrix_free(tie_table);

    // The origin has three active constraints in two dimensions, and (1,0)
    // is itself degenerate because both x<=1 and 2x+y<=2 are tight there.
    const std::vector<std::vector<double>> degenerate_neighbor = {
        {0, -1, -1}, {1, 1, 0}, {1, 0, 1}, {2, 2, 1}
    };
    const SearchResult found = search(degenerate_neighbor);
    require(found.status == 2, "incorrect degenerate-vertex neighbor count");
    require(contains(found, 1, 0) && contains(found, 0, 1),
            "incorrect degenerate-vertex neighbors");
    const SearchResult captured = search(
        degenerate_neighbor, 8, nullptr, true);
    require(captured.states.size() == captured.neighbors.size(),
            "one dictionary state was not saved per geometric neighbor");
    for (size_t i = 0; i < captured.states.size(); ++i)
    {
        const std::vector<double> reconstructed =
            state_point(captured.states[i], 2);
        require(std::fabs(reconstructed[0] - captured.neighbors[i][0]) < 1e-10 &&
                std::fabs(reconstructed[1] - captured.neighbors[i][1]) < 1e-10,
                "saved dictionary does not reproduce its neighbor");
        const SearchResult second_layer = resume(captured.states[i], 2);
        require(second_layer.status > 0 && contains(second_layer, 0, 0),
                "saved dictionary cannot resume the reverse search");
    }
    const SearchResult tnp_found = search(degenerate_neighbor, 8, "tnp");
    require(tnp_found.status == 2 && contains(tnp_found, 1, 0) &&
            contains(tnp_found, 0, 1),
            "TNP missed a degenerate neighbor from a feasible root basis");
    const SearchResult lex_found = search(degenerate_neighbor, 8, "lex");
    require(lex_found.status == 2 && contains(lex_found, 1, 0) &&
            contains(lex_found, 0, 1),
            "lexicographic search missed a degenerate neighbor");

    // The feasible wedge x<=2y, y<=2x, x+y<=1 has two adjacent vertices
    // (2/3,1/3) and (1/3,2/3), but the initial slack basis has no transition
    // column.  TNP alone cannot certify this vertex, whereas the
    // lexicographic N-tree fallback must enumerate both neighbors.
    const std::vector<std::vector<double>> internal_basis = {
        {0, 1, -2}, {0, -2, 1}, {1, 1, 1}
    };
    require(search(internal_basis, 8, "tnp").status == NO_T_COL,
            "TNP-only internal basis must report missing transition column");
    const SearchResult recovered = search(internal_basis);
    require(recovered.status == 2 &&
            contains(recovered, 2.0 / 3.0, 1.0 / 3.0) &&
            contains(recovered, 1.0 / 3.0, 2.0 / 3.0),
            "lexicographic fallback missed internal-basis neighbors");
    const SearchResult explicit_lex = search(internal_basis, 8, "lex");
    require(explicit_lex.status == 2 &&
            contains(explicit_lex, 2.0 / 3.0, 1.0 / 3.0) &&
            contains(explicit_lex, 1.0 / 3.0, 2.0 / 3.0),
            "explicit lexicographic search missed internal-basis neighbors");

    // A transition column is not by itself a Geue TNP starting certificate.
    // Here z is a transition column, but x+y<=0 fixes x=y=0, and no negative
    // basis can make the slack starting basis TNP-lexico-feasible.  The
    // standard N-tree fallback still finds the sole z-axis neighbor.
    const std::vector<std::vector<double>> infeasible_tnp_root = {
        {0, 1, 1, 0}, {1, 0, 0, 1}
    };
    require(search(infeasible_tnp_root, 8, "tnp").status == NO_T_COL,
            "TNP accepted a transition column without an L-feasible root");
    const SearchResult fallback_neighbor = search(infeasible_tnp_root);
    require(fallback_neighbor.status == 1 &&
            contains_point(fallback_neighbor, {0, 0, 1}),
            "standard N-tree did not recover a non-TNP root");

    // This root needs two ordered perturbation columns: z is negative in
    // x-z<=0 but zero in y-x<=0, so x must resolve the latter row.
    const std::vector<std::vector<double>> two_stage_tnp = {
        {0, 1, 0, -1}, {0, -1, 1, 0}, {1, 0, 0, 1}
    };
    const SearchResult staged_tnp = search(two_stage_tnp, 8, "tnp");
    require(staged_tnp.status == 3 &&
            contains_point(staged_tnp, {0, 0, 1}) &&
            contains_point(staged_tnp, {1, 0, 1}) &&
            contains_point(staged_tnp, {1, 1, 1}),
            "complete TNP perturbation missed a two-stage neighbor");

    // Equal first perturbation ratios are resolved by the later columns of
    // Bhat, not by tableau row order alone.
    const std::vector<std::vector<double>> tied_tnp = {
        {0, 1, 1, -1}, {0, 1, 1, -1}, {1, 0, 0, 1}
    };
    require(search(tied_tnp, 8, "tnp").status == 3 &&
            search(tied_tnp, 8, "lex").status == 3,
            "TNP or standard N-tree missed a duplicate-row neighbor");

    // Geue (1993), Example 1.1, translated from apex (0,0,0,1) to the
    // origin by z=1-x4.  Its other eight vertices are apex neighbors.
    const std::vector<std::vector<double>> geue_example = {
        {0, 3, -4, 0, -5}, {0, 0, 5, 12, -13},
        {0, -6, 0, 8, -10}, {0, 15, 36, 0, -39},
        {1, 0, 0, 0, 1}
    };
    const SearchResult geue_tnp = search(geue_example, 16, "tnp");
    const SearchResult geue_lex = search(geue_example, 16, "lex");
    require(geue_tnp.status == 8 && geue_lex.status == 8,
            "Geue's four-dimensional example has an incorrect neighbor count");
    for (const auto &point : geue_lex.neighbors)
        require(contains_point(geue_tnp, point),
                "TNP and standard N-tree disagree on Geue's example");

    require(search(degenerate_neighbor, 1).status == DEGENERACY,
            "neighbor-capacity exhaustion must remain non-certified");
    require(search(internal_basis, 1, "lex").status == DEGENERACY,
            "lexicographic capacity exhaustion must remain non-certified");

    // The leaving-row search must not silently drop a finite edge merely
    // because its ratio exceeds the old hard-coded 1e7 sentinel.
    const SearchResult long_edge = search({{2e8, 1.0}});
    require(long_edge.status == 1 &&
            std::fabs(long_edge.neighbors[0][0] - 2e8) < 1e-5,
            "finite large-ratio transition edge was dropped");

    // At the origin of x_i <= t (i=1..5), t+sum(u_j)<=1, x,t,u>=0,
    // each subset of the five x_i defines one edge with t>0.  Four
    // independent u axes add four more: 32+4=36 neighbors in dimension 10.
    // The old 3*n limit truncated this case; the new 5*n limit must finish.
    std::vector<std::vector<double>> many_neighbors(
        6, std::vector<double>(11, 0.0));
    for (int i = 0; i < 5; ++i)
    {
        many_neighbors[i][1 + i] = 1.0;
        many_neighbors[i][6] = -1.0;
    }
    many_neighbors[5][0] = 1.0;
    for (int i = 6; i <= 10; ++i) many_neighbors[5][i] = 1.0;
    unsetenv("DECIDEIT_NEIGHBOR_CAP_MULTIPLIER");
    unsetenv("DECIDEIT_REVERSE_NODE_CAP");
    require(adjacency_capacity(10) == 50, "10-dimensional capacity must be 5*n");
    setenv("DECIDEIT_NEIGHBOR_CAP_MULTIPLIER", "10", 1);
    require(adjacency_capacity(10) == 100,
            "runtime neighbor-cap override did not take effect");
    unsetenv("DECIDEIT_NEIGHBOR_CAP_MULTIPLIER");
    setenv("DECIDEIT_REVERSE_NODE_CAP", "4096", 1);
    require(reverse_node_capacity(20) == 4096,
            "runtime visited-basis cap override did not take effect");
    unsetenv("DECIDEIT_REVERSE_NODE_CAP");
    require(search(many_neighbors, 30).status == DEGENERACY,
            "old 3*n capacity should truncate the 36-neighbor test");
    const SearchResult expanded = search(many_neighbors, adjacency_capacity(10));
    const SearchResult lex_expanded = search(
        many_neighbors, adjacency_capacity(10), "lex");
    setenv("DECIDEIT_REVERSE_TRAVERSAL", "recursive", 1);
    const SearchResult recursive_expanded = search(
        many_neighbors, adjacency_capacity(10));
    unsetenv("DECIDEIT_REVERSE_TRAVERSAL");
    require(expanded.status == 36 && lex_expanded.status == 36,
            "5*n capacity should enumerate all 36 neighbors");
    require(recursive_expanded.status == expanded.status,
            "iterative and recursive Geue traversals disagree");
    setenv("DECIDEIT_DYNAMIC_SEARCH", "1", 1);
    setenv("DECIDEIT_DYNAMIC_SEARCH_BUDGET_MS", "1000", 1);
    setenv("DECIDEIT_DYNAMIC_NEIGHBOR_MAX", "50", 1);
    const SearchResult grown_neighbors = search(many_neighbors, 30);
    require(grown_neighbors.status == 36,
            "in-place neighbor growth did not finish the Geue traversal");
    for (const auto &point : expanded.neighbors)
        require(contains_point(grown_neighbors, point),
                "in-place neighbor growth lost an earlier neighbor");
    setenv("DECIDEIT_DYNAMIC_NEIGHBOR_MAX", "32", 1);
    require(search(many_neighbors, 30).status == DEGENERACY,
            "dynamic hard ceiling must not certify an incomplete traversal");
    setenv("DECIDEIT_DYNAMIC_NEIGHBOR_MAX", "50", 1);
    setenv("DECIDEIT_REVERSE_NODE_CAP", "1", 1);
    setenv("DECIDEIT_DYNAMIC_NODE_MAX", "4096", 1);
    require(search(many_neighbors, 50).status == 36,
            "in-place visited-basis growth did not finish the traversal");
    setenv("DECIDEIT_DYNAMIC_SEARCH_BUDGET_MS", "0", 1);
    require(search(many_neighbors, 30).status == NEIGHBOR_TIME_LIMIT,
            "exhausted time budget must not be treated as complete");
    unsetenv("DECIDEIT_DYNAMIC_SEARCH");
    unsetenv("DECIDEIT_DYNAMIC_SEARCH_BUDGET_MS");
    unsetenv("DECIDEIT_DYNAMIC_NEIGHBOR_MAX");
    unsetenv("DECIDEIT_DYNAMIC_NODE_MAX");
    unsetenv("DECIDEIT_REVERSE_NODE_CAP");
    for (int subset = 0; subset < 32; ++subset)
    {
        std::vector<double> point(10, 0.0);
        for (int i = 0; i < 5; ++i)
            point[i] = (subset >> i) & 1;
        point[5] = 1.0;
        require(contains_point(expanded, point) &&
                contains_point(lex_expanded, point),
                "TNP or lexicographic search missed a cube edge");
    }
    for (int axis = 6; axis < 10; ++axis)
    {
        std::vector<double> point(10, 0.0);
        point[axis] = 1.0;
        require(contains_point(expanded, point) &&
                contains_point(lex_expanded, point),
                "TNP or lexicographic search missed an independent edge");
    }
    for (int seed = 0; seed < 40; ++seed)
    {
        std::vector<std::vector<double>> rows = {
            {1, 1, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 1},
            {2.2, 1, 1, 1}
        };
        for (int t = 0; t < 3; ++t)
        {
            const int from = (seed + 2 * t) % 3;
            const int to = (seed / 3 + t + 1) % 3;
            if (from == to) continue;
            std::vector<double> row(4, 0.0);
            row[from + 1] = 1.0;
            row[to + 1] = -(2.0 + (seed + t) % 3);
            rows.push_back(row);
        }
        const auto expected = oracle_neighbors3(rows);
        const SearchResult actual = search(rows, 256);
        const SearchResult lex_actual = search(rows, 256, "lex");
        setenv("DECIDEIT_REVERSE_TRAVERSAL", "recursive", 1);
        const SearchResult recursive_actual = search(rows, 256);
        unsetenv("DECIDEIT_REVERSE_TRAVERSAL");
        require(actual.status == static_cast<int>(expected.size()) &&
                lex_actual.status == static_cast<int>(expected.size()) &&
                recursive_actual.status == static_cast<int>(expected.size()),
                "Geue search disagreed with independent 3D edge count");
        for (const auto &point : expected)
            require(contains_point(actual, point) &&
                    contains_point(lex_actual, point) &&
                    contains_point(recursive_actual, point),
                    "Geue search missed an independent-oracle 3D edge");
    }
    std::cout << "Reverse search: TNP, lexicographic fallback, capacity passed\n";
}
