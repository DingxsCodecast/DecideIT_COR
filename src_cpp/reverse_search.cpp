#include "dblp.h"
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <memory>
#include <unordered_set>

static thread_local bool search_truncated=false;
static thread_local bool search_tnp_invalid=false;
// Recursive Geue traversal retains a tableau at each level.  A much larger
// experimental node cap must not turn a long degenerate path into a C-stack
// overflow; reaching this guard is an incomplete search, not a certificate.
static const size_t kMaximumSearchDepth = 512;
static thread_local size_t search_max_depth = 0;

namespace {
struct BasisHash
{
    size_t operator()(const std::vector<int> &basis) const
    {
        size_t result = 1469598103934665603ULL;
        for (int index : basis)
            result = (result ^ static_cast<size_t>(index)) * 1099511628211ULL;
        return result;
    }
};

size_t dynamic_budget_milliseconds()
{
    const char *value = std::getenv("DECIDEIT_DYNAMIC_SEARCH_BUDGET_MS");
    if (!value) return 1000;
    if (!*value || *value == '-')
        throw std::invalid_argument("DECIDEIT_DYNAMIC_SEARCH_BUDGET_MS must be 0..60000");
    char *end = nullptr;
    errno = 0;
    const unsigned long parsed = std::strtoul(value, &end, 10);
    if (errno || *end || parsed > 60000)
        throw std::invalid_argument("DECIDEIT_DYNAMIC_SEARCH_BUDGET_MS must be 0..60000");
    return static_cast<size_t>(parsed);
}

struct SearchExpansionBudget
{
    bool enabled = false;
    bool timed_out = false;
    size_t node_limit = 0;
    size_t node_max = 0;
    size_t neighbor_max = 0;
    std::chrono::steady_clock::time_point deadline;

    SearchExpansionBudget(size_t dimension, size_t initial_nodes,
                          size_t initial_neighbors)
        : node_limit(initial_nodes), node_max(initial_nodes),
          neighbor_max(initial_neighbors)
    {
        const char *setting = std::getenv("DECIDEIT_DYNAMIC_SEARCH");
        if (setting && std::strcmp(setting, "0") != 0 &&
            std::strcmp(setting, "1") != 0)
            throw std::invalid_argument("DECIDEIT_DYNAMIC_SEARCH must be 0 or 1");
        enabled = setting && std::strcmp(setting, "1") == 0;
        if (!enabled) return;
        const size_t default_neighbors = dimension > VX_MAXVER / 20
            ? VX_MAXVER : 20 * dimension;
        neighbor_max = (std::max)(initial_neighbors,
            reverse_limit_override("DECIDEIT_DYNAMIC_NEIGHBOR_MAX",
                                   default_neighbors));
        const size_t default_nodes = initial_nodes > VX_MAXVER / 4
            ? VX_MAXVER : std::min(static_cast<size_t>(5000), 4 * initial_nodes);
        node_max = (std::max)(initial_nodes,
            reverse_limit_override("DECIDEIT_DYNAMIC_NODE_MAX",
                                   default_nodes));
        deadline = std::chrono::steady_clock::now() +
            std::chrono::milliseconds(dynamic_budget_milliseconds());
    }

    bool expired()
    {
        if (!enabled || timed_out) return timed_out;
        if (std::chrono::steady_clock::now() < deadline) return false;
        timed_out = true;
        search_truncated = true;
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"reverse_search_time_limit\"}\n");
        return true;
    }

    bool grow_nodes()
    {
        if (!enabled || expired() || node_limit >= node_max) return false;
        const size_t previous = node_limit;
        node_limit = std::min(node_max,
            previous > VX_MAXVER / 2 ? node_max : 2 * previous);
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"reverse_search_node_growth\",\"from\":%zu,\"to\":%zu}\n",
                    previous, node_limit);
        return true;
    }

    bool grow_neighbors(gsl_matrix *&neighbors)
    {
        if (!enabled || expired() || neighbors->size2 >= neighbor_max)
            return false;
        const size_t previous = neighbors->size2;
        const size_t expanded = std::min(neighbor_max,
            previous > VX_MAXVER / 2 ? neighbor_max : 2 * previous);
        grow_neighbor_matrix(neighbors, expanded);
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"reverse_search_neighbor_growth\",\"from\":%zu,\"to\":%zu}\n",
                    previous, expanded);
        return true;
    }
};
}

void grow_neighbor_matrix(gsl_matrix *&neighbors, size_t columns)
{
    if (columns <= neighbors->size2) return;
    if (columns > VX_MAXVER)
        throw std::invalid_argument("neighbor matrix exceeds safety capacity");
    gsl_matrix *grown = managed::matrix_calloc(neighbors->size1, columns);
    gsl_matrix_view prefix = gsl_matrix_submatrix(
        grown, 0, 0, neighbors->size1, neighbors->size2);
    gsl_matrix_memcpy(&prefix.matrix, neighbors);
    managed::matrix_free(neighbors);
    neighbors = grown;
}

void free_choose_t_col(gsl_matrix *o_table,
                       gsl_matrix *basis,
                       gsl_matrix *temp_inverse,
                       gsl_vector *vector_new,
                       //    gsl_vector_int *abscent,
                       gsl_vector_int *o_basic_index,
                       gsl_vector_int *o_nonbasic_index,
                       gsl_permutation *p)
{
    managed::matrix_free(o_table);
    managed::matrix_free(basis);
    managed::matrix_free(temp_inverse);
    managed::vector_free(vector_new);
    // managed::vector_int_free(abscent);
    managed::vector_int_free(o_basic_index);
    managed::vector_int_free(o_nonbasic_index);
    managed::permutation_free(p);
}

/*
    check_zeros() identifies which elements in a vector are zeros.

    inputs
        zero_index: to indicate zero or non-zero;
        b_bar: the vector to be checked.

    outputs
        number_zero: the # of zeros in the vector;
        zero_index: for each element, true for zero, false otherwise.
*/

int check_zeros(
    // bool zero_index[],
    std::vector<bool> &zero_index,
    gsl_vector *b_bar)
{
    size_t i, number_zero = 0;

    for (i = 0; i < b_bar->size; i++)
    {
        if (gsl_vector_get(b_bar, i) == 0)
        {
            zero_index[i] = true;
            number_zero++;
        }
        else
        {
            zero_index[i] = false;
        }
    }

    return number_zero;
}

// Compare the complete lexicographic ratio for L=(b,sign*Bhat).  The
// standard N-tree uses sign=+1 and Bhat=B0; Geue's improved N-tree uses
// sign=-1 and Bhat a basis chosen so that the root is L-feasible.  Testing
// only the first column of Bhat is not a lexicographic pivot rule when two
// zero rows have equal first perturbation ratios.
bool tnp_lex_feasible_node(const gsl_matrix *table,
                           const std::vector<int> &perturbation_columns)
{
    for (size_t row = 0; row < table->size1; ++row)
    {
        const double b = gsl_matrix_get(table, row, 0);
        if (b < 0) return false;
        if (b > 0) continue;
        bool positive = false;
        for (int column : perturbation_columns)
        {
            const double entry = gsl_matrix_get(table, row, column);
            if (entry == 0) continue;
            positive = entry < 0; // the perturbation basis is -Bhat
            break;
        }
        if (!positive) return false;
    }
    return true;
}

void chooseCandidateLex(gsl_vector_int *nonbasic_index, gsl_matrix *table,
                        const std::vector<int> &perturbation_columns,
                        double perturbation_sign,
                        gsl_matrix_int *candidates_pivots)
{
    const int m = static_cast<int>(table->size1);
    const int count = gsl_vector_int_get(nonbasic_index, 0);
    for (int j = 1; j <= count; ++j)
    {
        const int column = gsl_vector_int_get(nonbasic_index, j);
        int best_row = -1;
        for (int i = 0; i < m; ++i)
        {
            const double entry = gsl_matrix_get(table, i, column);
            if (entry <= 1e-12) continue;
            if (best_row < 0)
            {
                best_row = i;
                continue;
            }
            const double best_entry = gsl_matrix_get(table, best_row, column);
            // Compare b/a first, then every component of
            // sign*B^{-1}Bhat/a.  The full basis resolves all exact ties.
            for (size_t k = 0; k <= perturbation_columns.size(); ++k)
            {
                const int index = k == 0 ? 0 : perturbation_columns[k - 1];
                const double sign = k == 0 ? 1.0 : perturbation_sign;
                const double value = sign * gsl_matrix_get(table, i, index) /
                                     entry;
                const double best = sign *
                    gsl_matrix_get(table, best_row, index) / best_entry;
                const double tolerance = 1e-12 *
                    (1.0 + std::fabs(value) + std::fabs(best));
                if (value < best - tolerance)
                {
                    best_row = i;
                    break;
                }
                if (value > best + tolerance) break;
            }
        }
        if (best_row < 0) continue; // unbounded ray, not a neighbor
        gsl_matrix_int_set(candidates_pivots, 0, j, best_row);
        gsl_matrix_int_set(candidates_pivots, 1, j, column);
        gsl_matrix_int_set(candidates_pivots, 2, j,
                           gsl_matrix_get(table, best_row, 0) > 1e-10 ? 1 : 0);
    }
}

/*
    check_existance() checks the existance of an adjacent vertex.
    Note that we can only check a specific vector rather than its basis because we cannot ensure the degeneracy of the adjacent vertex.

    inputs
        vector_new: the newly located adjacent vertex to be checked;
        number_adjacent: the # of already located adjacent vertices;
        adj_vertices: the matrix containing already located adjacent vertices (in columns).

    outputs
        flag: true for existance, false otherwise.
*/

bool check_existance(
    gsl_vector *vector_new,
    int number_adjacent,
    gsl_matrix *adj_vertices)
{

    const int m = adj_vertices->size1;
    // The vertices are stored in original-variable coordinates.  Compare
    // directly instead of allocating two GSL vectors for every candidate.
    constexpr double squared_tolerance = 1e-16;
    for (int j = 0; j < number_adjacent; ++j)
    {
        double squared_distance = 0.0;
        for (int i = 0; i < m; ++i)
        {
            const double difference = gsl_matrix_get(adj_vertices, i, j) -
                                      gsl_vector_get(vector_new, i + 1);
            squared_distance += difference * difference;
            if (squared_distance > squared_tolerance) break;
        }
        if (squared_distance <= squared_tolerance) return true;
    }
    return false;
}

/*
    dfs() implements depth first search algorithm on the tableau.

    inputs
        table: the current tableau;
        basic_index, nonbasic_index: the indices of basic and non-basic variables;
        number_adjacent: the # of identified adjacent vertices;
        adj_vertices: the identified adjacent vertices;
        perturbation_columns/sign: the fixed complete lexicographic rule;
        visited: hash set of visited bases during the traversal;
        number_nodes: the # of visited nodes.

    outputs
        table: the tableau to be operated in the next round;
        basis: the basis corresponding to the tableau.
*/

void dfs(
    gsl_matrix *table,
    gsl_vector_int *basic_index,
    gsl_vector_int *nonbasic_index,
    size_t *number_adjacent,
    gsl_matrix *adj_vertices,
    const std::vector<int> &perturbation_columns,
    double perturbation_sign,
    size_t node_capacity,
    std::unordered_set<std::vector<int>, BasisHash> &visited,
    size_t *number_nodes,
    size_t depth)
{
    if (search_truncated) return;
    if (depth + 1 > search_max_depth) search_max_depth = depth + 1;
    if (perturbation_sign < 0 &&
        !tnp_lex_feasible_node(table, perturbation_columns))
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"reverse_search_tnp_lex_infeasible\"}\n");
        search_truncated = true;
        search_tnp_invalid = true;
        return;
    }

    const int m = table->size1, n = table->size2;
    int i, j, temp_swap, candidate_row, candidate_col;
    double theta, pivot_val, factor;

    gsl_vector *a = nullptr;
    gsl_vector *b = nullptr;
    gsl_vector *vector_new = nullptr;
    gsl_matrix *pre_table = nullptr;

    gsl_vector_view src_v, dest_v;

    gsl_matrix_int *candidates = managed::matrix_int_calloc(3, gsl_vector_int_get(nonbasic_index, 0) + 1);
    gsl_matrix_int_set_all(candidates, -1);

    chooseCandidateLex(nonbasic_index, table, perturbation_columns,
                       perturbation_sign, candidates);

    // show_gsl_matrix_int(candidates);

    for (j = 1; j <= gsl_vector_int_get(nonbasic_index, 0); j++)
    {
        if (search_truncated) break;
        if (gsl_matrix_int_get(candidates, 2, j) != -1)
        {
            candidate_row = gsl_matrix_int_get(candidates, 0, j);
            candidate_col = gsl_matrix_int_get(candidates, 1, j);

            temp_swap = gsl_vector_int_get(basic_index, candidate_row + 1);
            gsl_vector_int_set(basic_index, candidate_row + 1, candidate_col);
            // Candidate j corresponds to nonbasic_index[j].  Undoing this
            // pair of assignments is enough after the child search returns.
            gsl_vector_int_set(nonbasic_index, j, temp_swap);

            // show_gsl_vector_int(basic_index);
            // show_gsl_vector_int(nonbasic_index);

            if (gsl_matrix_int_get(candidates, 2, j) == 0) // degenerate basis pivot
            {
                std::vector<int> basis_key(m);
                for (i = 0; i < m; ++i)
                    basis_key[i] = gsl_vector_int_get(basic_index, i + 1);
                std::sort(basis_key.begin(), basis_key.end());

                if (visited.find(basis_key) != visited.end())
                {
                    gsl_vector_int_set(basic_index, candidate_row + 1, temp_swap);
                    gsl_vector_int_set(nonbasic_index, j, candidate_col);
                    continue;
                }
                if (*number_nodes >= node_capacity)
                {
                    if (std::getenv("DECIDEIT_TRACE"))
                        fprintf(stderr, "TRACE {\"event\":\"reverse_search_node_limit\",\"nodes\":%zu,\"capacity\":%zu}\n",
                                *number_nodes, node_capacity);
                    search_truncated = true;
                    gsl_vector_int_set(basic_index, candidate_row + 1, temp_swap);
                    gsl_vector_int_set(nonbasic_index, j, candidate_col);
                    break;
                }
                if (depth + 1 >= kMaximumSearchDepth)
                {
                    if (std::getenv("DECIDEIT_TRACE"))
                        fprintf(stderr, "TRACE {\"event\":\"reverse_search_depth_limit\",\"depth\":%zu,\"capacity\":%zu}\n",
                                depth + 1, kMaximumSearchDepth);
                    search_truncated = true;
                    gsl_vector_int_set(basic_index, candidate_row + 1, temp_swap);
                    gsl_vector_int_set(nonbasic_index, j, candidate_col);
                    break;
                }
                visited.insert(std::move(basis_key));
                ++(*number_nodes);
            }

            if (gsl_matrix_int_get(candidates, 2, j) == 1) // geometric neighbor
            {
                // for a non-degenerate pivot, it is only necessary to update the solution

                if (!a)
                {
                    a = managed::vector_calloc(m);
                    b = managed::vector_calloc(m);
                    vector_new = managed::vector_calloc(n);
                }

                gsl_matrix_get_col(a, table, 0);
                gsl_matrix_get_col(b, table, candidate_col);
                theta = gsl_vector_get(a, candidate_row) / gsl_vector_get(b, candidate_row);
                gsl_vector_set(a, candidate_row, theta);

                for (i = 0; i < m; i++)
                {
                    if (i != candidate_row)
                    {
                        gsl_vector_set(a, i, gsl_vector_get(a, i) - theta * gsl_vector_get(b, i));
                    }
                }

                construct_complete_solution(vector_new, basic_index, a);
                // show_gsl_vector(vector_new);

                if (!check_existance(vector_new, *number_adjacent, adj_vertices))
                {
                    src_v = gsl_vector_subvector(vector_new, 1, n - m - 1);

                    if (*number_adjacent >= adj_vertices->size2)
                    {
                        if (std::getenv("DECIDEIT_TRACE"))
                            fprintf(stderr, "TRACE {\"event\":\"reverse_search_neighbor_limit\",\"neighbors\":%zu,\"capacity\":%zu}\n",
                                    *number_adjacent, adj_vertices->size2);
                        search_truncated=true;
                        gsl_vector_int_set(basic_index, candidate_row + 1, temp_swap);
                        gsl_vector_int_set(nonbasic_index, j, candidate_col);
                        break;
                    }
                    else
                    {
                        gsl_matrix_set_col(adj_vertices, (*number_adjacent)++, &src_v.vector);
                    }
                    // show_gsl_matrix_columns(adj_vertices, *number_adjacent);
                }
            }
            else
            {
                // for a degenerate pivot, it is necessary to update the entire simplex tableau

                if (!pre_table)
                {
                    pre_table = managed::matrix_alloc(table->size1, table->size2);
                    gsl_matrix_memcpy(pre_table, table);
                }

                // show_gsl_matrix(table);
                src_v = gsl_matrix_row(table, candidate_row);
                pivot_val = gsl_matrix_get(table, candidate_row, candidate_col);
                gsl_vector_scale(&src_v.vector, 1.0 / pivot_val);

                for (i = 0; i < m; i++)
                {
                    if (i != candidate_row)
                    {
                        factor = gsl_matrix_get(table, i, candidate_col);
                        src_v = gsl_matrix_row(table, candidate_row);
                        dest_v = gsl_matrix_row(table, i);

                        gsl_blas_daxpy(-factor, &src_v.vector, &dest_v.vector);
                    }
                }

                resolve_zero_matrix(table);
                // show_gsl_matrix(table);

                dfs(table, basic_index, nonbasic_index, number_adjacent,
                    adj_vertices, perturbation_columns, perturbation_sign,
                    node_capacity, visited, number_nodes, depth + 1);

                gsl_matrix_memcpy(table, pre_table);
            }

            gsl_vector_int_set(basic_index, candidate_row + 1, temp_swap);
            gsl_vector_int_set(nonbasic_index, j, candidate_col);
        }
    }

    managed::matrix_int_free(candidates);
    if (a) managed::vector_free(a);
    if (b) managed::vector_free(b);
    if (vector_new) managed::vector_free(vector_new);
    if (pre_table) managed::matrix_free(pre_table);
}

// The same Geue candidate rule and depth-first order as dfs(), but the
// traversal state lives on the heap.  A parent's tableau snapshot is reused
// for all its children, so neither recursion depth nor repeated sibling
// copies depend on the C stack.
void iterative_dfs(
    gsl_matrix *table, gsl_vector_int *basic_index,
    gsl_vector_int *nonbasic_index, size_t *number_adjacent,
    gsl_matrix *&adj_vertices,
    const std::vector<int> &perturbation_columns, double perturbation_sign,
    SearchExpansionBudget &budget,
    std::unordered_set<std::vector<int>, BasisHash> &visited,
    size_t *number_nodes,
    std::vector<NeighborDictionaryState> *neighbor_states)
{
    const auto save_neighbor_state = [&](int pivot_row, int pivot_column,
                                         const gsl_vector *point) {
        if (!neighbor_states) return;
        NeighborDictionaryState state;
        state.rows = static_cast<int>(table->size1);
        state.columns = static_cast<int>(table->size2);
        state.tableau.resize(table->size1 * table->size2);
        gsl_matrix *pivoted = managed::matrix_alloc(table->size1,
                                                    table->size2);
        gsl_matrix_memcpy(pivoted, table);
        const double pivot = gsl_matrix_get(pivoted, pivot_row, pivot_column);
        gsl_vector_view pivot_values = gsl_matrix_row(pivoted, pivot_row);
        gsl_vector_scale(&pivot_values.vector, 1.0 / pivot);
        for (int i = 0; i < state.rows; ++i)
        {
            if (i == pivot_row) continue;
            const double factor = gsl_matrix_get(pivoted, i, pivot_column);
            gsl_vector_view source = gsl_matrix_row(pivoted, pivot_row);
            gsl_vector_view destination = gsl_matrix_row(pivoted, i);
            gsl_blas_daxpy(-factor, &source.vector, &destination.vector);
        }
        resolve_zero_matrix(pivoted);
        for (int i = 0; i < state.rows; ++i)
            for (int j = 0; j < state.columns; ++j)
                state.tableau[static_cast<size_t>(i) * state.columns + j] =
                    gsl_matrix_get(pivoted, i, j);
        state.basic_index.resize(basic_index->size);
        state.nonbasic_index.resize(nonbasic_index->size);
        state.point.resize(point->size);
        for (size_t i = 0; i < basic_index->size; ++i)
            state.basic_index[i] = gsl_vector_int_get(basic_index, i);
        for (size_t i = 0; i < nonbasic_index->size; ++i)
            state.nonbasic_index[i] = gsl_vector_int_get(nonbasic_index, i);
        for (size_t i = 0; i < point->size; ++i)
            state.point[i] = gsl_vector_get(point, i);
        neighbor_states->push_back(std::move(state));
        managed::matrix_free(pivoted);
    };
    struct Frame
    {
        gsl_matrix_int *candidates = nullptr;
        gsl_matrix *snapshot = nullptr;
        gsl_vector *a = nullptr;
        gsl_vector *b = nullptr;
        gsl_vector *vector_new = nullptr;
        int next_j = 1;
        int incoming_row = -1;
        int incoming_j = -1;
        int incoming_col = -1;
        int old_basic = -1;

        Frame(gsl_matrix *table, gsl_vector_int *nonbasic,
              const std::vector<int> &perturbation, double sign,
              int row = -1, int j = -1,
              int col = -1, int previous_basic = -1)
            : incoming_row(row), incoming_j(j), incoming_col(col),
              old_basic(previous_basic)
        {
            candidates = managed::matrix_int_calloc(
                3, gsl_vector_int_get(nonbasic, 0) + 1);
            gsl_matrix_int_set_all(candidates, -1);
            if (sign < 0 && !tnp_lex_feasible_node(table, perturbation))
            {
                if (std::getenv("DECIDEIT_TRACE"))
                    fprintf(stderr, "TRACE {\"event\":\"reverse_search_tnp_lex_infeasible\"}\n");
                search_truncated = true;
                search_tnp_invalid = true;
                return;
            }
            chooseCandidateLex(nonbasic, table, perturbation, sign, candidates);
        }
        ~Frame()
        {
            managed::matrix_int_free(candidates);
            if (snapshot) managed::matrix_free(snapshot);
            if (a) managed::vector_free(a);
            if (b) managed::vector_free(b);
            if (vector_new) managed::vector_free(vector_new);
        }
    };

    const int m = static_cast<int>(table->size1);
    const int n = static_cast<int>(table->size2);
    const int nonbasic_count = gsl_vector_int_get(nonbasic_index, 0);
    std::vector<std::unique_ptr<Frame>> stack;
    stack.emplace_back(new Frame(table, nonbasic_index,
                                 perturbation_columns, perturbation_sign));
    search_max_depth = 1;

    while (!stack.empty())
    {
        if (budget.expired()) search_truncated = true;
        // A capacity/time stop will return an incomplete result, so no
        // tableau restoration is needed. Only a TNP-to-lex fallback needs
        // the root tableau reconstructed while unwinding.
        if (search_truncated && (!search_tnp_invalid || budget.timed_out))
            break;
        Frame &frame = *stack.back();
        if (search_truncated || frame.next_j > nonbasic_count)
        {
            if (frame.incoming_row >= 0)
            {
                Frame &parent = *stack[stack.size() - 2];
                gsl_matrix_memcpy(table, parent.snapshot);
                gsl_vector_int_set(basic_index, frame.incoming_row + 1,
                                   frame.old_basic);
                gsl_vector_int_set(nonbasic_index, frame.incoming_j,
                                   frame.incoming_col);
            }
            stack.pop_back();
            continue;
        }

        const int j = frame.next_j++;
        const int type = gsl_matrix_int_get(frame.candidates, 2, j);
        if (type == -1) continue;
        const int row = gsl_matrix_int_get(frame.candidates, 0, j);
        const int col = gsl_matrix_int_get(frame.candidates, 1, j);
        const int old_basic = gsl_vector_int_get(basic_index, row + 1);
        gsl_vector_int_set(basic_index, row + 1, col);
        gsl_vector_int_set(nonbasic_index, j, old_basic);

        if (type == 1)
        {
            if (!frame.a)
            {
                frame.a = managed::vector_calloc(m);
                frame.b = managed::vector_calloc(m);
                frame.vector_new = managed::vector_calloc(n);
            }
            gsl_matrix_get_col(frame.a, table, 0);
            gsl_matrix_get_col(frame.b, table, col);
            const double theta = gsl_vector_get(frame.a, row) /
                                 gsl_vector_get(frame.b, row);
            gsl_vector_set(frame.a, row, theta);
            for (int i = 0; i < m; ++i)
                if (i != row)
                    gsl_vector_set(frame.a, i, gsl_vector_get(frame.a, i) -
                                   theta * gsl_vector_get(frame.b, i));
            construct_complete_solution(frame.vector_new, basic_index,
                                        frame.a);
            if (!check_existance(frame.vector_new, *number_adjacent,
                                 adj_vertices))
            {
                if (*number_adjacent >= adj_vertices->size2)
                {
                    if (!budget.grow_neighbors(adj_vertices))
                    {
                        if (!budget.timed_out && std::getenv("DECIDEIT_TRACE"))
                            fprintf(stderr, "TRACE {\"event\":\"reverse_search_neighbor_limit\",\"neighbors\":%zu,\"capacity\":%zu}\n",
                                    *number_adjacent, adj_vertices->size2);
                        search_truncated = true;
                    }
                }
                if (!search_truncated)
                {
                    gsl_vector_view point = gsl_vector_subvector(
                        frame.vector_new, 1, n - m - 1);
                    gsl_matrix_set_col(adj_vertices, *number_adjacent,
                                       &point.vector);
                    save_neighbor_state(row, col, &point.vector);
                    ++(*number_adjacent);
                }
            }
            gsl_vector_int_set(basic_index, row + 1, old_basic);
            gsl_vector_int_set(nonbasic_index, j, col);
            continue;
        }

        std::vector<int> basis_key(m);
        for (int i = 0; i < m; ++i)
            basis_key[i] = gsl_vector_int_get(basic_index, i + 1);
        std::sort(basis_key.begin(), basis_key.end());
        if (visited.find(basis_key) != visited.end())
        {
            gsl_vector_int_set(basic_index, row + 1, old_basic);
            gsl_vector_int_set(nonbasic_index, j, col);
            continue;
        }
        if (*number_nodes >= budget.node_limit && !budget.grow_nodes())
        {
            if (!budget.timed_out && std::getenv("DECIDEIT_TRACE"))
                fprintf(stderr, "TRACE {\"event\":\"reverse_search_node_limit\",\"nodes\":%zu,\"capacity\":%zu}\n",
                        *number_nodes, budget.node_limit);
            search_truncated = true;
            gsl_vector_int_set(basic_index, row + 1, old_basic);
            gsl_vector_int_set(nonbasic_index, j, col);
            continue;
        }
        visited.insert(std::move(basis_key));
        ++(*number_nodes);
        if (!frame.snapshot)
        {
            frame.snapshot = managed::matrix_alloc(table->size1,
                                                   table->size2);
            gsl_matrix_memcpy(frame.snapshot, table);
        }
        const double pivot = gsl_matrix_get(table, row, col);
        gsl_vector_view pivot_row = gsl_matrix_row(table, row);
        gsl_vector_scale(&pivot_row.vector, 1.0 / pivot);
        for (int i = 0; i < m; ++i)
        {
            if (i == row) continue;
            const double factor = gsl_matrix_get(table, i, col);
            gsl_vector_view source = gsl_matrix_row(table, row);
            gsl_vector_view destination = gsl_matrix_row(table, i);
            gsl_blas_daxpy(-factor, &source.vector, &destination.vector);
        }
        resolve_zero_matrix(table);
        stack.emplace_back(new Frame(table, nonbasic_index,
                                     perturbation_columns, perturbation_sign,
                                     row, j, col, old_basic));
        if (stack.size() > search_max_depth)
            search_max_depth = stack.size();
    }
}

void union_set(gsl_vector_int *o_basic_index,
               gsl_vector_int *abscent,
               gsl_vector_int *basic_index)
{
    int i, index = 1;

    gsl_vector_int_set(basic_index, 0, 0);

    for (i = 1; i <= gsl_vector_int_get(o_basic_index, 0); i++)
    {
        gsl_vector_int_set(basic_index, index++, gsl_vector_int_get(o_basic_index, i));
    }

    for (i = 1; i <= gsl_vector_int_get(abscent, 0); i++)
    {
        gsl_vector_int_set(basic_index, index++, gsl_vector_int_get(abscent, i));
    }

    gsl_sort_vector_int(basic_index);
    gsl_vector_int_set(basic_index, 0, gsl_vector_int_get(o_basic_index, 0) + gsl_vector_int_get(abscent, 0));
}

/*
    check_transition_column() tells whether a column corresponding to a non-basic column in the tableau is qualified as a transition column.

    inputs
        column: the column corresponding to a non-basic variables;
        b_bar: the solution column corresponding to basic_index.

    output
        true for a qualified transition column, false otherwise.
*/

bool check_transition_column(
    gsl_vector *column,
    // bool zero_index[],
    std::vector<bool> &zero_index,
    int number_zero)
{
    const int m = column->size;
    int i, counter_zero = 0, counter_nonzero = 0;

    for (i = 0; i < m; i++)
    {
        if (gsl_vector_get(column, i) <= 0 && zero_index[i] == true)
        {
            counter_zero++;
        }
        else if (gsl_vector_get(column, i) > 0 && zero_index[i] == false)
        {
            counter_nonzero++;
        }
    }

    if (counter_zero == number_zero && counter_nonzero > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

namespace {

// The columns of Bhat are represented by their indices in the root tableau.
// Modified Gram-Schmidt is used only to reject singular/ill-conditioned
// candidate bases; all search pivots still use the original tableau.
bool independent_column(const gsl_matrix *table, int column,
                        const std::vector<std::vector<double>> &orthogonal,
                        std::vector<double> &normalized)
{
    const size_t m = table->size1;
    normalized.resize(m);
    double original_norm = 0.0;
    for (size_t i = 0; i < m; ++i)
    {
        normalized[i] = gsl_matrix_get(table, i, column);
        original_norm += normalized[i] * normalized[i];
    }
    original_norm = std::sqrt(original_norm);
    if (!std::isfinite(original_norm) || original_norm <= 1e-10)
        return false;
    // Reorthogonalization keeps the rank test stable when many slack columns
    // are needed to complete Bhat.
    for (int pass = 0; pass < 2; ++pass)
        for (const auto &q : orthogonal)
        {
            double projection = 0.0;
            for (size_t i = 0; i < m; ++i)
                projection += normalized[i] * q[i];
            for (size_t i = 0; i < m; ++i)
                normalized[i] -= projection * q[i];
        }
    double residual_norm = 0.0;
    for (double value : normalized)
        residual_norm += value * value;
    residual_norm = std::sqrt(residual_norm);
    if (!std::isfinite(residual_norm) ||
        residual_norm <= 1e-10 * (original_norm > 1.0 ? original_norm : 1.0))
        return false;
    for (double &value : normalized) value /= residual_norm;
    return true;
}

// Geue's TNP theorem requires L=(b,-Bhat) with Bhat a genuine basis and the
// root L-feasible.  A transition column alone is only a necessary consequence
// of that construction, not a sufficient starting condition.  Try transition
// columns with the most strictly negative entries on zero rows first, complete
// an L-feasible Bhat greedily, and let the caller use standard N-tree if none
// can be certified.  This selection is heuristic; correctness comes from the
// explicit checks, not from finding an optimal Bhat.
bool build_tnp_perturbation(const gsl_matrix *table,
                            const gsl_vector_int *basic_index,
                            const gsl_vector_int *nonbasic_index,
                            std::vector<int> &columns)
{
    const int m = static_cast<int>(table->size1);
    const int n = static_cast<int>(table->size2);
    std::vector<bool> zero_rows(m, false);
    int zero_count = 0;
    for (int row = 0; row < m; ++row)
        if (gsl_matrix_get(table, row, 0) == 0)
        {
            zero_rows[row] = true;
            ++zero_count;
        }

    std::vector<std::pair<int, int>> transitions;
    for (int j = 1; j <= gsl_vector_int_get(nonbasic_index, 0); ++j)
    {
        const int column = gsl_vector_int_get(nonbasic_index, j);
        int strict = 0;
        bool valid = true, finite_edge = false;
        for (int row = 0; row < m; ++row)
        {
            const double value = gsl_matrix_get(table, row, column);
            if (zero_rows[row])
            {
                if (value > 0) valid = false;
                if (value < 0) ++strict;
            }
            else if (value > 0)
                finite_edge = true;
        }
        if (valid && finite_edge) transitions.emplace_back(-strict, column);
    }
    std::sort(transitions.begin(), transitions.end());

    for (const auto &transition : transitions)
    {
        const int t_col = transition.second;
        if (-transition.first == zero_count)
        {
            // Strict negativity on every zero row makes root feasibility
            // immediate.  Replacing one root-basic column by t gives a
            // nonsingular Bhat without an O(m^3) rank construction.
            int drop_row = -1;
            double largest = 0.0;
            for (int row = 0; row < m; ++row)
            {
                const double magnitude =
                    std::fabs(gsl_matrix_get(table, row, t_col));
                if (magnitude > largest)
                {
                    largest = magnitude;
                    drop_row = row;
                }
            }
            if (drop_row >= 0 && largest > 1e-10)
            {
                columns.clear();
                columns.push_back(t_col);
                for (int row = 0; row < m; ++row)
                    if (row != drop_row)
                        columns.push_back(
                            gsl_vector_int_get(basic_index, row + 1));
                return true;
            }
        }
        std::vector<int> trial;
        std::vector<bool> selected(n, false), unresolved = zero_rows;
        std::vector<std::vector<double>> orthogonal;
        std::vector<double> normalized;
        if (!independent_column(table, t_col, orthogonal, normalized))
            continue;
        trial.push_back(t_col);
        selected[t_col] = true;
        orthogonal.push_back(normalized);
        int remaining = zero_count;
        for (int row = 0; row < m; ++row)
            if (unresolved[row] && gsl_matrix_get(table, row, t_col) < 0)
            {
                unresolved[row] = false;
                --remaining;
            }

        while (remaining > 0)
        {
            int best_column = -1, best_count = 0;
            std::vector<double> best_normalized;
            for (int column = 1; column < n; ++column)
            {
                if (selected[column]) continue;
                int resolved = 0;
                bool valid = true;
                for (int row = 0; row < m; ++row)
                    if (unresolved[row])
                    {
                        const double value = gsl_matrix_get(table, row, column);
                        if (value > 0) { valid = false; break; }
                        if (value < 0) ++resolved;
                    }
                if (!valid || resolved <= best_count ||
                    !independent_column(table, column, orthogonal, normalized))
                    continue;
                best_column = column;
                best_count = resolved;
                best_normalized = normalized;
            }
            if (best_column < 0) break;
            trial.push_back(best_column);
            selected[best_column] = true;
            orthogonal.push_back(std::move(best_normalized));
            for (int row = 0; row < m; ++row)
                if (unresolved[row] &&
                    gsl_matrix_get(table, row, best_column) < 0)
                {
                    unresolved[row] = false;
                    --remaining;
                }
        }
        if (remaining > 0) continue;

        // Root basic columns are identity columns in this tableau and are
        // inexpensive, well-conditioned choices for completing the basis.
        for (int pass = 0; pass < 2 && trial.size() < table->size1; ++pass)
        {
            const int count = pass == 0 ? m : n - 1;
            for (int j = 1; j <= count && trial.size() < table->size1; ++j)
            {
                const int column = pass == 0 ?
                    gsl_vector_int_get(basic_index, j) : j;
                if (selected[column] ||
                    !independent_column(table, column, orthogonal, normalized))
                    continue;
                trial.push_back(column);
                selected[column] = true;
                orthogonal.push_back(normalized);
            }
        }
        if (trial.size() != table->size1) continue;

        bool feasible = true;
        for (int row = 0; row < m && feasible; ++row)
        {
            if (!zero_rows[row]) continue;
            bool resolved = false;
            for (int column : trial)
            {
                const double value = gsl_matrix_get(table, row, column);
                if (value == 0) continue;
                resolved = value < 0; // row of B_root^{-1}(-Bhat) is positive
                break;
            }
            if (!resolved) feasible = false;
        }
        if (feasible)
        {
            columns = std::move(trial);
            return true;
        }
    }
    return false;
}
} // namespace

/*
    choose_t_col() tests the identified starting basis for a transition
    column.  It does not enumerate alternative degenerate bases; if this
    basis is internal, the caller must use a complete fallback or report
    NO_T_COL, never optimality.

    inputs
        basic_index: the identified basic indices;
        nonbasic_index: the identified non-basic indices;
        o_table: [A_bar, b_bar].

    outputs
        t_col: the index of the selected t column in the table;
        table: the current tableau corresponding to the selected t column;
        basic_index, nonbasic_index: the basic and non-basic indices corresponding to the table.
*/

int choose_t_col(gsl_vector_int *basic_index,
                 gsl_vector_int *nonbasic_index,
                 gsl_matrix *working_table)
{
    const int m = working_table->size1, n = working_table->size2, full_dimension = n - 1, dimension = full_dimension - m;
    int k, number_zero;
    // bool zero_index[m];
    std::vector<bool> zero_index(m);

    gsl_vector_int *o_basic_index = managed::vector_int_calloc(m + 1);
    gsl_vector_int_memcpy(o_basic_index, basic_index);

    gsl_vector_int *o_nonbasic_index = managed::vector_int_calloc(n);
    gsl_vector_int_memcpy(o_nonbasic_index, nonbasic_index);

    gsl_matrix *o_table = managed::matrix_calloc(m, n);
    gsl_matrix_memcpy(o_table, working_table);

    gsl_permutation *p = managed::permutation_alloc(m);
    gsl_vector *vector_new = managed::vector_calloc(full_dimension + 1);
    gsl_matrix *basis = managed::matrix_calloc(m, m);
    gsl_matrix *temp_inverse = managed::matrix_calloc(m, m);

    gsl_matrix_view matrix_consistency = gsl_matrix_submatrix(o_table, 0, 1, m, dimension);
    gsl_vector_view vector_consistency = gsl_matrix_column(o_table, 0);

    // number_absent = m - gsl_vector_int_get(o_basic_index, 0);
    // number_available = gsl_vector_int_get(o_nonbasic_index, 0);
    // number_total = nchoosek(number_available, number_absent);

    // gsl_vector_int *abscent = managed::vector_int_calloc(number_absent + 1);

    gsl_vector_view src_v, col_v_from_matrix, solution_v;

    // for (j = 0; j < number_total; j++)
    // {
    //     nchoosek_enum(number_available, number_absent, j + 1, combination);
    //     // show_array_int_all(combination, number_absent);

    //     gsl_vector_int_set(abscent, 0, number_absent);

    //     for (i = 1; i <= number_absent; i++)
    //     {
    //         gsl_vector_int_set(abscent, i, gsl_vector_int_get(o_nonbasic_index, combination[i - 1]));
    //     }

    //     union_set(o_basic_index, abscent, basic_index);

    //     for (k = 1; k <= m; k++)
    //     {
    //         src_v = gsl_matrix_column(o_table, gsl_vector_int_get(basic_index, k));
    //         gsl_matrix_set_col(basis, k - 1, &src_v.vector);
    //     }

    //     if (condition_number(basis, 1e12))
    //     {
    //         setDifference(o_nonbasic_index->data, abscent->data, nonbasic_index->data);

    //         if (gsl_linalg_LU_decomp(basis, p, &sign_permutation) != 0)
    //         {
    //             free_choose_t_col(o_table, basis, temp_inverse, vector_new, abscent, o_basic_index, o_nonbasic_index, p);
    //             return (GSL_LU);
    //         }

    //         if (gsl_linalg_LU_invert(basis, p, temp_inverse) != 0)
    //         {
    //             free_choose_t_col(o_table, basis, temp_inverse, vector_new, abscent, o_basic_index, o_nonbasic_index, p);
    //             return (GSL_LU);
    //         }

    //         // show_gsl_matrix(temp_inverse);
    //         gsl_matrix_set_zero(working_table);
    //         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, temp_inverse, o_table, 0, working_table);
    //         resolve_zero_matrix(working_table);

    //         src_v = gsl_matrix_column(working_table, 0);
    //         construct_complete_solution(vector_new, basic_index, &src_v.vector);
    //         solution_v = gsl_vector_subvector(vector_new, 1, dimension);

    //         // show_gsl_vector(&solution_v.vector);

    //         if (check_consistency_reverse_search(&matrix_consistency.matrix, &vector_consistency.vector, &solution_v.vector))
    //         {
    //             number_zero = check_zeros(zero_index, &src_v.vector);

    //             for (k = 0; k < gsl_vector_int_get(nonbasic_index, 0); k++)
    //             {
    //                 col_v_from_matrix = gsl_matrix_column(working_table, gsl_vector_int_get(nonbasic_index, k + 1));

    //                 if (check_transition_column(&col_v_from_matrix.vector, zero_index, number_zero))
    //                 {
    //                     free_choose_t_col(o_table, basis, temp_inverse, vector_new, abscent, o_basic_index, o_nonbasic_index, p);

    //                     return gsl_vector_int_get(nonbasic_index, k + 1);
    //                 }
    //             }
    //         }
    //     }
    // }

    for (k = 1; k <= m; k++)
    {
        src_v = gsl_matrix_column(o_table, gsl_vector_int_get(basic_index, k));
        gsl_matrix_set_col(basis, k - 1, &src_v.vector);
    }

    if (safe_inverse(basis, temp_inverse, 1e-12) != 0)
    {
        free_choose_t_col(o_table, basis, temp_inverse, vector_new, o_basic_index, o_nonbasic_index, p);

        return (GSL_LU);
    }

    // if (gsl_linalg_LU_decomp(basis, p, &sign_permutation) != GSL_SUCCESS)
    // {
    //     free_choose_t_col(o_table, basis, temp_inverse, vector_new, o_basic_index, o_nonbasic_index, p);
    //     return (GSL_LU);
    // }

    // if (gsl_linalg_LU_invert(basis, p, temp_inverse) != GSL_SUCCESS)
    // {
    //     free_choose_t_col(o_table, basis, temp_inverse, vector_new, o_basic_index, o_nonbasic_index, p);
    //     return (GSL_LU);
    // }

    // show_gsl_matrix(temp_inverse);
    gsl_matrix_set_zero(working_table);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, temp_inverse, o_table, 0, working_table);
    resolve_zero_matrix(working_table);

    src_v = gsl_matrix_column(working_table, 0);
    construct_complete_solution(vector_new, basic_index, &src_v.vector);
    solution_v = gsl_vector_subvector(vector_new, 1, dimension);

    // show_gsl_vector(&solution_v.vector);

    if (check_consistency_reverse_search(&matrix_consistency.matrix, &vector_consistency.vector, &solution_v.vector))
    {
        number_zero = check_zeros(zero_index, &src_v.vector);

        for (k = 0; k < gsl_vector_int_get(nonbasic_index, 0); k++)
        {
            col_v_from_matrix = gsl_matrix_column(working_table, gsl_vector_int_get(nonbasic_index, k + 1));

            if (check_transition_column(&col_v_from_matrix.vector, zero_index, number_zero))
            {
                free_choose_t_col(o_table, basis, temp_inverse, vector_new, o_basic_index, o_nonbasic_index, p);

                return gsl_vector_int_get(nonbasic_index, k + 1);
            }
        }
    }

    free_choose_t_col(o_table, basis, temp_inverse, vector_new, o_basic_index, o_nonbasic_index, p);

    return (0);
}

/*
    reverse_search() implements the two neighborhood-enumeration variants in
    Geue (1993): the improved N-tree with transition-node pivoting (TNP), and
    the standard lexicographic N-tree.  TNP is used only when a full negative
    basis Bhat makes the root L-feasible; otherwise the standard variant is
    the complete fallback.  Merely finding a transition column is not enough.

    inputs
        tableau: the tableau containing the transformed effective linear constraints, intervals (0 <= x <= ub - lb) and slack variables;
        lb: lower bounds for the effective variables in X.
        basic_index and nonbasic_index: the indices of identified basic and non-basic variables;
        adj_vertices: the matrix containing all vertices adjacent to x.

    outputs
        number_adjacent: the # of identified adjacent vertices;
        adj_vertices: all adjacent vertices in the original space and dimension (in columns).
*/

int reverse_search(gsl_matrix *tableau,
                   gsl_vector_int *basic_index,
                   gsl_vector_int *nonbasic_index,
                   gsl_matrix *&adj_vertices,
                   std::vector<NeighborDictionaryState> *neighbor_states)
{
    int t_col;
    size_t num_neighbors = 0, num_nodes = 0;
    size_t *number_adjacent = &num_neighbors, *number_nodes = &num_nodes;

    search_truncated=false;
    search_tnp_invalid=false;
    search_max_depth=0;
    if (neighbor_states) neighbor_states->clear();
    const char *configured_search = std::getenv("DECIDEIT_REVERSE_SEARCH");
    const bool force_lex = configured_search &&
        std::strcmp(configured_search, "lex") == 0;
    const bool force_tnp = configured_search &&
        std::strcmp(configured_search, "tnp") == 0;
    const size_t dimension = gsl_vector_int_get(nonbasic_index, 0);
    const size_t node_capacity = reverse_node_capacity(dimension);
    SearchExpansionBudget budget(dimension, node_capacity, adj_vertices->size2);
    t_col = choose_t_col(basic_index, nonbasic_index, tableau);
    if (budget.expired()) return NEIGHBOR_TIME_LIMIT;
    if (t_col == GSL_LU) return GSL_LU;
    std::vector<int> tnp_columns;
    const bool tnp_available = !force_lex && t_col > 0 &&
        build_tnp_perturbation(tableau, basic_index, nonbasic_index,
                               tnp_columns);
    if (budget.expired()) return NEIGHBOR_TIME_LIMIT;
    if (force_tnp && !tnp_available)
    {
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"reverse_search\",\"mode\":\"improved_tnp\",\"nodes\":0,\"neighbors\":0,\"status\":%d}\n",
                    NO_T_COL);
        return NO_T_COL;
    }
    bool lex_search = !tnp_available;

    std::vector<int> root_basis(tableau->size1);
    for (size_t i = 0; i < tableau->size1; ++i)
        root_basis[i] = gsl_vector_int_get(basic_index, i + 1);
    const std::vector<int> perturbation_columns =
        lex_search ? root_basis : tnp_columns;
    const double perturbation_sign = lex_search ? 1.0 : -1.0;
    const std::vector<int> standard_columns = root_basis;
    std::sort(root_basis.begin(), root_basis.end());
    std::unordered_set<std::vector<int>, BasisHash> visited;
    visited.insert(root_basis);
    *number_nodes = 1;

    const char *traversal = std::getenv("DECIDEIT_REVERSE_TRAVERSAL");
    if (budget.enabled && traversal &&
        std::strcmp(traversal, "recursive") == 0)
        throw std::invalid_argument(
            "dynamic neighbor search requires iterative traversal");
    if (neighbor_states && traversal &&
        std::strcmp(traversal, "recursive") == 0)
        throw std::invalid_argument(
            "neighbor-state capture requires iterative traversal");
    const auto traverse = [&](const std::vector<int> &columns, double sign)
    {
        if (traversal && std::strcmp(traversal, "recursive") == 0)
            dfs(tableau, basic_index, nonbasic_index, number_adjacent,
                adj_vertices, columns, sign, node_capacity, visited,
                number_nodes, 0);
        else
            iterative_dfs(tableau, basic_index, nonbasic_index,
                          number_adjacent, adj_vertices, columns, sign,
                          budget, visited, number_nodes, neighbor_states);
    };
    if (!budget.expired()) traverse(perturbation_columns, perturbation_sign);
    if (search_tnp_invalid && !force_tnp && !budget.timed_out)
    {
        // A floating-point pivot left the verified TNP subgraph.  The
        // partially enumerated neighbor list is not a certificate; restart
        // from the restored root with the standard positive-basis rule.
        if (std::getenv("DECIDEIT_TRACE"))
            fprintf(stderr, "TRACE {\"event\":\"reverse_search_tnp_fallback\"}\n");
        search_truncated = false;
        search_tnp_invalid = false;
        search_max_depth = 0;
        *number_adjacent = 0;
        if (neighbor_states) neighbor_states->clear();
        *number_nodes = 1;
        visited.clear();
        visited.insert(root_basis);
        lex_search = true;
        if (!budget.expired()) traverse(standard_columns, 1.0);
    }

    // Also check after the final candidate: a long last pivot or duplicate
    // check must not turn an expired search into a completion certificate.
    budget.expired();

    if (std::getenv("DECIDEIT_TRACE"))
        fprintf(stderr, "TRACE {\"event\":\"reverse_search\",\"mode\":\"%s\",\"nodes\":%zu,\"neighbors\":%zu,\"max_depth\":%zu,\"status\":%d}\n",
                lex_search ? "standard_lex" : "improved_tnp",
                *number_nodes, *number_adjacent, search_max_depth,
                search_truncated ? (budget.timed_out ? NEIGHBOR_TIME_LIMIT : DEGENERACY)
                                 : static_cast<int>(*number_adjacent));
    return search_truncated
        ? (budget.timed_out ? NEIGHBOR_TIME_LIMIT : DEGENERACY)
        : static_cast<int>(*number_adjacent);
}
