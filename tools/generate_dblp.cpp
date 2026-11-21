#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

struct Range { double lo, hi; };
struct SparseRow { std::vector<int> index; std::vector<double> coefficient; double rhs; int type; };

struct Options
{
    int n1 = 0, n2 = 0, m_probability = 0, m_value = 0;
    unsigned long long seed = 1;
    std::string output;
    Range probability_slack{0.05, 0.25};
    Range value_slack{0.05, 0.35};
    Range x_bound_slack{0.05, 0.35};
    Range y_bound_slack{0.01, 0.25};
    double type_probability[3] = {0.45, 0.20, 0.35};
    double compound_probability[3] = {0.45, 0.35, 0.20};
    double minimum_bound_width = 1e-3;
    double lower_floor_fraction = 0.2;
    double target_objective_lower_bound = -0.8;
    double anchor_target_safety = 0.05;
    double block2_anchor_floor_fraction = 0.2;
    double simplex_alpha = 2.0;
    int round_digits = 12;
    int max_total_dimension = 60;
    long long max_rlt_rows = 100000;
    bool overwrite = false;
};

double parse_double(const char *text, const char *name)
{
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(text, &end);
    if (errno || !end || *end || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid ") + name);
    return value;
}

long long parse_integer(const char *text, const char *name)
{
    char *end = nullptr;
    errno = 0;
    const long long value = std::strtoll(text, &end, 10);
    if (errno || !end || *end)
        throw std::invalid_argument(std::string("invalid ") + name);
    return value;
}

void require_range(const Range &range, const char *name)
{
    if (range.lo < 0.0 || range.hi < range.lo)
        throw std::invalid_argument(std::string(name) + " must satisfy 0 <= lo <= hi");
}

Options parse_options(int argc, char **argv)
{
    Options result;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument = argv[i];
        auto value = [&](const char *name) -> const char * {
            if (++i >= argc) throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[i];
        };
        auto range = [&](const char *name) -> Range {
            const double lo = parse_double(value(name), name);
            const double hi = parse_double(value(name), name);
            return {lo, hi};
        };
        if (argument == "--help")
        {
            std::cout << "Usage: generate_dblp --output FILE --n1 N --n2 N "
                         "--m-probability N --m-value N [--seed N] [options]\n"
                         "Options: --probability-slack LO HI --value-slack LO HI "
                         "--x-bound-slack LO HI --y-bound-slack LO HI "
                         "--type-probabilities P2 P3 P4 "
                         "--compound-probabilities P2 P3 P4 "
                         "--target-objective-lower-bound X --simplex-alpha X "
                         "--max-total-dimension N --max-rlt-rows N --overwrite\n";
            std::exit(0);
        }
        else if (argument == "--output") result.output = value("--output");
        else if (argument == "--n1") result.n1 = static_cast<int>(parse_integer(value("--n1"), "n1"));
        else if (argument == "--n2") result.n2 = static_cast<int>(parse_integer(value("--n2"), "n2"));
        else if (argument == "--m-probability") result.m_probability = static_cast<int>(parse_integer(value("--m-probability"), "m-probability"));
        else if (argument == "--m-value") result.m_value = static_cast<int>(parse_integer(value("--m-value"), "m-value"));
        else if (argument == "--seed") result.seed = static_cast<unsigned long long>(parse_integer(value("--seed"), "seed"));
        else if (argument == "--probability-slack") result.probability_slack = range("--probability-slack");
        else if (argument == "--value-slack") result.value_slack = range("--value-slack");
        else if (argument == "--x-bound-slack") result.x_bound_slack = range("--x-bound-slack");
        else if (argument == "--y-bound-slack") result.y_bound_slack = range("--y-bound-slack");
        else if (argument == "--type-probabilities")
            for (double &entry : result.type_probability) entry = parse_double(value("--type-probabilities"), "type probability");
        else if (argument == "--compound-probabilities")
            for (double &entry : result.compound_probability) entry = parse_double(value("--compound-probabilities"), "compound probability");
        else if (argument == "--target-objective-lower-bound") result.target_objective_lower_bound = parse_double(value("--target-objective-lower-bound"), "target objective lower bound");
        else if (argument == "--simplex-alpha") result.simplex_alpha = parse_double(value("--simplex-alpha"), "simplex alpha");
        else if (argument == "--min-bound-width") result.minimum_bound_width = parse_double(value("--min-bound-width"), "minimum bound width");
        else if (argument == "--lower-floor-fraction") result.lower_floor_fraction = parse_double(value("--lower-floor-fraction"), "lower floor fraction");
        else if (argument == "--anchor-target-safety") result.anchor_target_safety = parse_double(value("--anchor-target-safety"), "anchor target safety");
        else if (argument == "--block2-anchor-floor-fraction") result.block2_anchor_floor_fraction = parse_double(value("--block2-anchor-floor-fraction"), "block2 anchor floor fraction");
        else if (argument == "--round-digits") result.round_digits = static_cast<int>(parse_integer(value("--round-digits"), "round digits"));
        else if (argument == "--max-total-dimension") result.max_total_dimension = static_cast<int>(parse_integer(value("--max-total-dimension"), "max total dimension"));
        else if (argument == "--max-rlt-rows") result.max_rlt_rows = parse_integer(value("--max-rlt-rows"), "max RLT rows");
        else if (argument == "--overwrite") result.overwrite = true;
        else throw std::invalid_argument("unexpected argument: " + argument);
    }
    const int n = result.n1 + result.n2;
    if (result.max_total_dimension < 1 || result.max_total_dimension > 300)
        throw std::invalid_argument("max total dimension must lie in [1,300]");
    if (result.output.empty() || result.n1 < 1 || result.n2 < 1 ||
        n > result.max_total_dimension)
        throw std::invalid_argument("output is required and dimensions must satisfy n1,n2 >= 1 and n1+n2 <= max total dimension (default 60)");
    if (result.m_probability < 0 || result.m_value < 0 ||
        result.m_probability > 60 || result.m_value > 60)
        throw std::invalid_argument("constraint counts must lie in [0,60]");
    require_range(result.probability_slack, "probability slack");
    require_range(result.value_slack, "value slack");
    require_range(result.x_bound_slack, "x-bound slack");
    require_range(result.y_bound_slack, "y-bound slack");
    if (result.seed > static_cast<unsigned long long>(std::numeric_limits<long long>::max()) ||
        result.simplex_alpha <= 0.0 || result.minimum_bound_width <= 0.0 ||
        result.minimum_bound_width >= 1.0 || result.lower_floor_fraction < 0.0 ||
        result.lower_floor_fraction >= 1.0 || result.target_objective_lower_bound <= -1.0 ||
        result.target_objective_lower_bound >= 1.0 || result.anchor_target_safety < 0.0 ||
        result.block2_anchor_floor_fraction < 0.0 || result.block2_anchor_floor_fraction >= 1.0 ||
        result.round_digits < 0 || result.round_digits > 15 || result.max_rlt_rows < 1)
        throw std::invalid_argument("invalid generator option range");
    const double type_sum = std::accumulate(result.type_probability, result.type_probability + 3, 0.0);
    const double compound_sum = std::accumulate(result.compound_probability, result.compound_probability + 3, 0.0);
    if (type_sum <= 0.0 || compound_sum <= 0.0 ||
        std::any_of(result.type_probability, result.type_probability + 3, [](double x){ return x < 0.0; }) ||
        std::any_of(result.compound_probability, result.compound_probability + 3, [](double x){ return x < 0.0; }))
        throw std::invalid_argument("probability weights must be nonnegative with positive sums");
    const long long rows = result.m_probability + result.m_value +
        4LL * n * n + 2LL * n * (result.m_probability + result.m_value) + 2LL + 2LL * n;
    if (rows > result.max_rlt_rows)
        throw std::invalid_argument("estimated cross-term RLT exceeds --max-rlt-rows");
    return result;
}

class Generator
{
public:
    explicit Generator(const Options &options) : options_(options), engine_(options.seed) {}

    double uniform(const Range &range)
    {
        std::uniform_real_distribution<double> distribution(range.lo, range.hi);
        return distribution(engine_);
    }

    int weighted(const std::vector<double> &weights)
    {
        std::vector<double> normalized = weights;
        if (std::accumulate(normalized.begin(), normalized.end(), 0.0) <= 0.0)
            std::fill(normalized.begin(), normalized.end(), 1.0);
        std::discrete_distribution<int> distribution(normalized.begin(), normalized.end());
        return distribution(engine_);
    }

    std::vector<int> sample_indices(int count)
    {
        std::vector<int> values(options_.n1 + options_.n2);
        std::iota(values.begin(), values.end(), 1);
        std::shuffle(values.begin(), values.end(), engine_);
        values.resize(count);
        return values;
    }

    std::vector<double> simplex(int count)
    {
        if (count == 1) return {1.0};
        std::gamma_distribution<double> gamma(options_.simplex_alpha, 1.0);
        std::vector<double> values(count);
        double sum = 0.0;
        for (double &value : values) { value = gamma(engine_); sum += value; }
        for (double &value : values) value /= sum;
        return values;
    }

    std::vector<double> value_anchor()
    {
        const double target = options_.target_objective_lower_bound;
        double upper = 1.0 - target - options_.anchor_target_safety - 1e-3;
        if (upper > 0.9) upper = 0.9;
        if (upper < 1e-3) upper = 1e-3;
        double lower = 0.35 * upper;
        if (lower < 1e-3) lower = 1e-3;
        if (lower > upper) lower = upper;
        double maximum_block2 = uniform({lower, upper});
        double minimum_block1 = maximum_block2 + target + options_.anchor_target_safety;
        if (minimum_block1 < 0.0) minimum_block1 = 0.0;
        if (minimum_block1 >= 1.0)
        {
            maximum_block2 = 1.0 - target - options_.anchor_target_safety - 1e-3;
            if (maximum_block2 < 0.0) maximum_block2 = 0.0;
            minimum_block1 = maximum_block2 + target + options_.anchor_target_safety;
        }
        if (minimum_block1 >= 1.0) throw std::runtime_error("target cannot be reached");
        std::vector<double> result(options_.n1 + options_.n2);
        for (int i = 0; i < options_.n1; ++i)
            result[i] = uniform({minimum_block1, 1.0});
        const double minimum_block2 = options_.block2_anchor_floor_fraction * maximum_block2;
        for (int i = options_.n1; i < options_.n1 + options_.n2; ++i)
            result[i] = uniform({minimum_block2, maximum_block2});
        return result;
    }

    void bounds(const std::vector<double> &anchor, const Range &slack,
                std::vector<double> &lower, std::vector<double> &upper)
    {
        lower.resize(anchor.size()); upper.resize(anchor.size());
        for (size_t i = 0; i < anchor.size(); ++i)
        {
            lower[i] = anchor[i] - uniform(slack);
            upper[i] = anchor[i] + uniform(slack);
            if (lower[i] < 0.0) lower[i] = 0.0;
            if (upper[i] > 1.0) upper[i] = 1.0;
            const double floor = options_.lower_floor_fraction * anchor[i];
            if (lower[i] < floor) lower[i] = floor;
            if (upper[i] - lower[i] < options_.minimum_bound_width)
            {
                lower[i] = std::max(0.0,
                    std::min(lower[i], anchor[i] - 0.5 * options_.minimum_bound_width));
                upper[i] = std::min(1.0,
                    std::max(upper[i], anchor[i] + 0.5 * options_.minimum_bound_width));
            }
            if (upper[i] - lower[i] < options_.minimum_bound_width)
            {
                if (lower[i] <= 0.0) upper[i] = std::min(1.0, options_.minimum_bound_width);
                else if (upper[i] >= 1.0) lower[i] = std::max(0.0, 1.0 - options_.minimum_bound_width);
                else
                {
                    lower[i] = std::max(0.0, anchor[i] - options_.minimum_bound_width);
                    upper[i] = std::min(1.0, anchor[i] + options_.minimum_bound_width);
                }
            }
            if (anchor[i] < lower[i]) lower[i] = anchor[i];
            if (anchor[i] > upper[i]) upper[i] = anchor[i];
        }
    }

    void enforce_target(const std::vector<double> &anchor,
                        std::vector<double> &lower, std::vector<double> &upper)
    {
        double maximum_block2 = 0.0;
        for (int i = options_.n1; i < options_.n1 + options_.n2; ++i)
            if (anchor[i] > maximum_block2) maximum_block2 = anchor[i];
        const double required_block1 = options_.target_objective_lower_bound + maximum_block2;
        for (int i = 0; i < options_.n1; ++i)
        {
            if (lower[i] < required_block1) lower[i] = required_block1;
            if (lower[i] > anchor[i]) lower[i] = anchor[i];
            if (upper[i] < anchor[i]) upper[i] = anchor[i];
            if (upper[i] < lower[i] + options_.minimum_bound_width)
                upper[i] = std::min(1.0, std::max(anchor[i], lower[i] + options_.minimum_bound_width));
        }
        for (int i = options_.n1; i < options_.n1 + options_.n2; ++i)
        {
            const double cap = uniform({anchor[i], maximum_block2});
            if (upper[i] > cap) upper[i] = cap;
            if (upper[i] < anchor[i]) upper[i] = anchor[i];
            if (lower[i] > upper[i] - options_.minimum_bound_width)
                lower[i] = std::max(0.0, std::min(anchor[i], upper[i] - options_.minimum_bound_width));
        }
        double certified_left = lower[0], certified_right = upper[options_.n1];
        for (int i = 1; i < options_.n1; ++i) certified_left = std::min(certified_left, lower[i]);
        for (int i = options_.n1 + 1; i < options_.n1 + options_.n2; ++i)
            certified_right = std::max(certified_right, upper[i]);
        if (certified_left - certified_right < options_.target_objective_lower_bound - 1e-10)
            throw std::runtime_error("target bound enforcement failed");
    }

    std::vector<SparseRow> rows(int count, const std::vector<double> &anchor,
                                const Range &slack)
    {
        std::vector<SparseRow> result;
        result.reserve(count);
        const int n = options_.n1 + options_.n2;
        const double scale = std::pow(10.0, options_.round_digits);
        for (int row_number = 0; row_number < count; ++row_number)
        {
            std::vector<double> type_weights{options_.type_probability[0],
                                             options_.type_probability[1],
                                             options_.type_probability[2]};
            if (n < 4) type_weights[1] = 0.0;
            const int choice = weighted(type_weights);
            SparseRow row;
            row.type = choice == 0 ? 2 : choice == 1 ? 3 : 4;
            int cardinality = 2;
            if (row.type == 3) cardinality = 4;
            if (row.type == 4)
            {
                std::vector<double> card_weights;
                for (int card = 2; card <= std::min(4, n); ++card)
                    card_weights.push_back(options_.compound_probability[card - 2]);
                cardinality = 2 + weighted(card_weights);
            }
            row.index = sample_indices(cardinality);
            row.coefficient.assign(cardinality, 1.0);
            if (row.type == 2) row.coefficient[1] = -1.0;
            if (row.type == 3) row.coefficient = {1.0, -1.0, -1.0, 1.0};
            double lhs = 0.0;
            for (int k = 0; k < cardinality; ++k)
                lhs += row.coefficient[k] * anchor[row.index[k] - 1];
            row.rhs = std::round((lhs + uniform(slack)) * scale) / scale;
            if (row.rhs < lhs) row.rhs = lhs + 1.0 / scale;
            result.push_back(row);
        }
        return result;
    }

private:
    const Options &options_;
    std::mt19937_64 engine_;
};

std::string stem(const std::string &path)
{
    const size_t slash = path.find_last_of("/\\");
    const size_t begin = slash == std::string::npos ? 0 : slash + 1;
    const size_t dot = path.find_last_of('.');
    return path.substr(begin, dot == std::string::npos || dot < begin ? std::string::npos : dot - begin);
}

void emit_row(std::ostream &out, const SparseRow &row, int n1)
{
    out << "3 " << row.index.size();
    for (size_t k = 0; k < row.index.size(); ++k)
    {
        const int index = row.index[k];
        out << ' ' << (index <= n1 ? 1 : 2)
            << ' ' << (index <= n1 ? index : index - n1)
            << ' ' << row.coefficient[k];
    }
    out << ' ' << row.rhs << '\n';
}

void write_instance(const Options &options,
                    const std::vector<double> &probability_lower,
                    const std::vector<double> &probability_upper,
                    const std::vector<double> &value_lower,
                    const std::vector<double> &value_upper,
                    const std::vector<SparseRow> &probability_rows,
                    const std::vector<SparseRow> &value_rows)
{
    if (!options.overwrite && std::ifstream(options.output.c_str()).good())
        throw std::runtime_error("output already exists (use --overwrite)");
    const std::string temporary = options.output + ".tmp." + std::to_string(static_cast<long long>(getpid()));
    std::ofstream out(temporary.c_str());
    if (!out) throw std::runtime_error("cannot create temporary output");
    out << std::fixed << std::setprecision(options.round_digits);
    const int n = options.n1 + options.n2;
    out << "8.8.4\n" << stem(options.output) << "\n1\n0\n2 "
        << options.n1 << ' ' << options.n2 << "\nalternative_A\nalternative_B\n"
        << n + probability_rows.size() << '\n';
    for (int i = 0; i < n; ++i)
        out << "1 " << (i < options.n1 ? 1 : 2) << ' '
            << (i < options.n1 ? i + 1 : i - options.n1 + 1)
            << " 1 " << probability_lower[i] << ' ' << probability_upper[i] << '\n';
    for (const SparseRow &row : probability_rows) emit_row(out, row, options.n1);
    out << "0 1\n" << n + value_rows.size() << '\n';
    for (int i = 0; i < n; ++i)
        out << "1 " << (i < options.n1 ? 1 : 2) << ' '
            << (i < options.n1 ? i + 1 : i - options.n1 + 1)
            << " 1 " << value_lower[i] << ' ' << value_upper[i] << '\n';
    for (const SparseRow &row : value_rows) emit_row(out, row, options.n1);
    out << "0\n0\n";
    out.close();
    if (!out) { std::remove(temporary.c_str()); throw std::runtime_error("failed writing output"); }
    if (std::rename(temporary.c_str(), options.output.c_str()) != 0)
    {
        std::remove(temporary.c_str());
        throw std::runtime_error("cannot install generated output");
    }
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        Generator generator(options);
        const std::vector<double> value_anchor = generator.value_anchor();
        std::vector<double> probability_anchor = generator.simplex(options.n1);
        const std::vector<double> second = generator.simplex(options.n2);
        probability_anchor.insert(probability_anchor.end(), second.begin(), second.end());
        std::vector<double> xl, xu, yl, yu;
        generator.bounds(value_anchor, options.x_bound_slack, xl, xu);
        generator.enforce_target(value_anchor, xl, xu);
        generator.bounds(probability_anchor, options.y_bound_slack, yl, yu);
        const std::vector<SparseRow> probability_rows = generator.rows(
            options.m_probability, probability_anchor, options.probability_slack);
        const std::vector<SparseRow> value_rows = generator.rows(
            options.m_value, value_anchor, options.value_slack);
        write_instance(options, yl, yu, xl, xu, probability_rows, value_rows);
        const int n = options.n1 + options.n2;
        const long long rlt_rows = options.m_probability + options.m_value +
            4LL * n * n + 2LL * n * (options.m_probability + options.m_value) + 2LL + 2LL * n;
        std::cout << "GENERATED {\"file\":\"" << options.output
                  << "\",\"seed\":" << options.seed
                  << ",\"n1\":" << options.n1 << ",\"n2\":" << options.n2
                  << ",\"m_probability\":" << options.m_probability
                  << ",\"m_value\":" << options.m_value
                  << ",\"max_total_dimension\":" << options.max_total_dimension
                  << ",\"estimated_rlt_columns\":" << 2 * n + n * n
                  << ",\"estimated_rlt_rows\":" << rlt_rows << "}\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 2;
    }
}
