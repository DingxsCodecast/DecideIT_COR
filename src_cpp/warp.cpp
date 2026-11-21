/* The vertex algorithm is described in the paper
   Sundgren, D., Danielson, M., and Ekenberg, L.,
   Warp effects on calculating interval probabilities,
   International Journal of Approximate Reasoning 50 (2009),
   pp. 1360-1368. Code written by Mats Danielson 2008. */

// #define VX_MAXDIM 16    // max nbr of prob vars
// #define VX_MAXVER 32768 // max nbr of vertices inside hyper-pyramide (=2^(VX_MAXDIM-1))
// #define WARP_EPS 1.0E-6 // minimum (epsilon) mass difference
// #define WARP_EPS10 10.0 * WARP_EPS

// #include "zor.h"

#include "dblp.h"

// generate_mp is called concurrently by the minimum and maximum directions.
// These arrays and counters are per-search scratch state, not shared model
// data. Keep the vector objects thread-local, but their storage on the heap:
// a fixed s_path[VX_MAXDIM+1][VX_MAXVER+1] would put about 120 MB in TLS and
// can make pthread creation fail on Linux. The path table is sized only for
// the current problem dimension and its capacity is reused by that thread.
static thread_local std::vector<double> sigma(VX_MAXVER + 1);
// static double delta[VX_MAXVER + 1];
static thread_local std::vector<double> s_pow(VX_MAXVER + 1);
static thread_local std::vector<int> upnodes(VX_MAXVER + 1);
static thread_local std::vector<int> s_path;
static thread_local int s_count;
static thread_local bool vertex_limit;
static thread_local unsigned long long search_nodes;
static const unsigned long long SEARCH_NODE_LIMIT = 2000000ULL;

static int &path_entry(int dimension, int vertex)
{
    return s_path[static_cast<size_t>(dimension) * (VX_MAXVER + 1) + vertex];
}

static int f1_T(double value, double target, int cur, int stop, int path[], h_vector lobo, h_vector upbo)
{
    int i, upn, ofl;

    if (vertex_limit || ++search_nodes > SEARCH_NODE_LIMIT)
    {
        vertex_limit = true;
        return 1;
    }

    /* Check if still inside hyper-pyramid */
    if (value > target - 1.0E-6)
        return 1;
    /* Check if end of path */
    if (cur == stop)
    {
        if (s_count >= VX_MAXVER) {
            vertex_limit=true;
            return 1; // path exhausted
        }
        ++s_count;
        sigma[s_count] = value;
        s_pow[s_count] = pow(target - value, (double)(stop - 1));
        upn = 0;
        for (i = 1; i <= stop; i++)
        {
            /* Record the up/down path... */
            path_entry(i, s_count) = path[i];
            /* ...and the nbr of up-nodes */
            if (path[i])
                upn++;
        }
        upnodes[s_count] = upn;
        return 0;
    }
    /* Split into two paths and try both */
    path[cur + 1] = 0;
    ofl = f1_T(value + lobo[cur + 1], target, cur + 1, stop, path, lobo, upbo);
    path[cur + 1] = 1;
    if (!ofl)
        f1_T(value + upbo[cur + 1], target, cur + 1, stop, path, lobo, upbo);
    return ofl;
}

/* Vertex-based mass point algorithm */

static double ijar_warp_T(int x, int dim, double target, double sum2, h_vector lobo, h_vector upbo)
{
    int i;
    double cij, sum1, term, mp;

    if (dim <= 1)
        /* Only one dim -> no vertices */
        return 1.0;
    /* Sum weights from all vertices for this dimension x */
    sum1 = 0.0;
    for (i = 1; i <= s_count; i++)
    {
        if (path_entry(x, i))
            cij = upbo[x];
        else
            cij = lobo[x];
        term = s_pow[i] * (dim * cij + target - sigma[i]);
        if (upnodes[i] % 2)
            sum1 -= term;
        else
            sum1 += term;
    }
    mp = sum1 / (dim * sum2);
    return mp;
}

static double ijar_warp_sum2()
{
    int i;
    double sum2;

    /* Sum normalisation weights from all vertices */
    sum2 = 0.0;
    for (i = 1; i <= s_count; i++)
    {
        if (upnodes[i] % 2)
            sum2 -= s_pow[i];
        else
            sum2 += s_pow[i];
    }
    return sum2;
}

/* Generate mass point according to the vertex algorithm.
 * 'lobo' and 'upbo' contain the lower and upper bounds.
 * The mass point is returned in the 'mid' vector.
 * Collapsed dimensions should be excluded before calling.
 * Their total probability is deducted from 1.0 -> 'target'.
 * Only for one-level P/W-trees -> sum to 'target' overall. */

int generate_mp(int n_cons, double target, h_vector lobo, h_vector mid, h_vector upbo)
{
    int j, sd_path[VX_MAXDIM + 1];
    double sum2;

    /* Check if function can start */
    if ((n_cons < 1) || (n_cons > VX_MAXDIM))
    {
        printf("Wrong number of dimensions (%d)\n", n_cons);
        return -1;
    }
    if ((target < WARP_EPS10) || (target > 1.0))
    {
        printf("Illegal mass target (%.3lf)\n", target);
        return -1;
    }
    for (j = 1; j <= n_cons; j++)
        /* Detect collapsed dimensions */
        if (upbo[j] - lobo[j] < WARP_EPS10)
        {
            printf("P%d: too narrow hull (%lg)\n", j, upbo[j] - lobo[j]);
            return -1;
        }
    /* Initialise invariants */
    s_path.resize(static_cast<size_t>(n_cons + 1) * (VX_MAXVER + 1));
    s_count = 0;
    vertex_limit = false;
    search_nodes = 0;
    f1_T(0.0, target, 0, n_cons, sd_path, lobo, upbo); // modifies s_count
    if(vertex_limit) {
        fprintf(stderr,"Exact focal search limit reached; using projected fallback "
                "(vertices=%d, nodes=%llu)\n",
                VX_MAXVER, SEARCH_NODE_LIMIT);
        return -1;
    }
    if (!s_count)
    {
        printf("No path\n");
        return -1;
    }
    sum2 = ijar_warp_sum2();
    if (sum2 < NORM_LIMIT)
    {
        printf("No normaliser\n");
        return -1;
    }
    /* Calculate mass point */
    for (j = 1; j <= n_cons; j++)
    {
        mid[j] = ijar_warp_T(j, n_cons, target, sum2, lobo, upbo);
    }
    return n_cons;
}

static int projected_midpoint_focal(int n, h_vector lo, h_vector mid,
                                    h_vector hi)
{
    double sum_lo = 0.0, sum_hi = 0.0;
    for (int i = 1; i <= n; ++i)
    {
        sum_lo += lo[i];
        sum_hi += hi[i];
    }
    if (sum_lo > 1.0 + 1e-10 || sum_hi < 1.0 - 1e-10)
        return -1;

    // Euclidean projection of the interval midpoint onto
    // {p : lo <= p <= hi, sum(p)=1}.  The scalar multiplier is found by
    // monotone bisection, so this path is O(n) and remains practical at the
    // 60-variable design limit when exact vertex integration is exponential.
    double lower_shift = -1.0, upper_shift = 1.0;
    for (int iteration = 0; iteration < 100; ++iteration)
    {
        const double shift = 0.5 * (lower_shift + upper_shift);
        double sum = 0.0;
        for (int i = 1; i <= n; ++i)
        {
            const double center = 0.5 * (lo[i] + hi[i]);
            sum += std::fmax(lo[i], std::fmin(hi[i], center + shift));
        }
        if (sum < 1.0) lower_shift = shift;
        else upper_shift = shift;
    }
    const double shift = 0.5 * (lower_shift + upper_shift);
    double sum = 0.0;
    for (int i = 1; i <= n; ++i)
    {
        const double center = 0.5 * (lo[i] + hi[i]);
        mid[i] = std::fmax(lo[i], std::fmin(hi[i], center + shift));
        sum += mid[i];
    }

    // Remove the last few ulps of bisection error while preserving bounds.
    double residual = 1.0 - sum;
    for (int i = 1; i <= n && std::fabs(residual) > 1e-14; ++i)
    {
        const double room = residual > 0.0 ? hi[i] - mid[i] : mid[i] - lo[i];
        const double change = std::copysign(std::fmin(std::fabs(residual), room),
                                            residual);
        mid[i] += change;
        residual -= change;
    }
    return std::fabs(residual) <= 1e-10 ? n : -1;
}

int probability_focal(int n,h_vector lo,h_vector mid,h_vector hi) {
    h_vector l={},u={},v={};std::vector<int> free_indices;
    double target=1.;
    for(int i=1;i<=n;i++){
        if(lo[i]==hi[i]){mid[i]=lo[i];target-=lo[i];}
        else {free_indices.push_back(i);l[free_indices.size()]=lo[i];u[free_indices.size()]=hi[i];}
    }
    const int m=free_indices.size();
    if(!m)return std::fabs(target)<=1e-8?n:-1;
    if(m==1){if(target<l[1]-1e-10||target>u[1]+1e-10)return -1;mid[free_indices[0]]=target;return n;}
    // Normalize the remaining mass so generate_mp's one-dimensional and
    // collapsed-coordinate conventions are never applied to the original tree.
    if(target<=0||target>1)return -1;
    for(int i=1;i<=m;i++){l[i]/=target;u[i]/=target;}
    // Exact mass-point integration enumerates vertices.  Preserve it for the
    // small cases for which it is effective; use a deterministic central
    // feasible point at larger dimension or if enumeration reaches its guard.
    if(m>20 || generate_mp(m,1.,l,v,u)!=m)
        return projected_midpoint_focal(n,lo,mid,hi);
    for(int i=0;i<m;i++)mid[free_indices[i]]=target*v[i+1];
    return n;
}
