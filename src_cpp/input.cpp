#include "dblp.h"

#include <fstream>
#include <set>
#include <sstream>
#include <string>

namespace
{
std::string next_line(std::ifstream &in)
{
    std::string text;
    if (!std::getline(in, text))
        throw std::invalid_argument("Truncated DDT file");
    return text;
}

void require_complete(std::istringstream &in)
{
    std::string extra;
    if (!in || (in >> extra))
        throw std::invalid_argument("Invalid DDT row");
}

int variable(int alternative, int consequence, const std::vector<int> &dimension)
{
    if (alternative < 1 || alternative > 2 || consequence < 1 ||
        consequence > dimension[alternative])
        throw std::invalid_argument("DDT variable index out of range");
    return consequence + (alternative == 2 ? dimension[1] : 0);
}

void read_base(std::ifstream &in, const std::vector<int> &dimension,
               d_row lower, d_row upper, link_set links,
               std::vector<LinCons> &rows)
{
    std::istringstream header(next_line(in));
    int statement_count = 0;
    header >> statement_count;
    require_complete(header);
    const int variables = dimension[1] + dimension[2];
    if (statement_count < 0 || statement_count > 5 * VX_MAXDIM)
        throw std::invalid_argument("DDT statement count exceeds capacity");

    std::vector<bool> seen(variables + 1, false);
    rows.clear();
    int link_count = 0;
    for (int statement = 0; statement < statement_count; ++statement)
    {
        std::istringstream input(next_line(in));
        int type = 0;
        input >> type;
        if (!input)
            throw std::invalid_argument("Invalid DDT statement");

        if (type == 1)
        {
            int alternative = 0, consequence = 0, sign = 0;
            double lo = 0.0, hi = 0.0;
            input >> alternative >> consequence >> sign >> lo >> hi;
            require_complete(input);
            const int index = variable(alternative, consequence, dimension);
            if (sign != 1 || seen[index])
                throw std::invalid_argument("Invalid or duplicate interval statement");
            if (!std::isfinite(lo) || !std::isfinite(hi) || lo < 0.0 ||
                hi > 1.0 || lo > hi)
                throw std::invalid_argument(
                    "Intervals must satisfy 0 <= lower <= upper <= 1");
            lower[index] = lo;
            upper[index] = hi;
            seen[index] = true;
        }
        else if (type == 2)
        {
            int alternative_1 = 0, consequence_1 = 0, sign_1 = 0;
            int alternative_2 = 0, consequence_2 = 0, sign_2 = 0;
            double lo = 0.0, hi = 0.0;
            input >> alternative_1 >> consequence_1 >> sign_1
                  >> alternative_2 >> consequence_2 >> sign_2 >> lo >> hi;
            require_complete(input);
            int first = variable(alternative_1, consequence_1, dimension);
            int second = variable(alternative_2, consequence_2, dimension);
            if ((sign_1 != 1 && sign_1 != -1) || sign_2 != -sign_1 ||
                first == second || !std::isfinite(lo) || !std::isfinite(hi) ||
                lo > hi)
                throw std::invalid_argument("Invalid difference statement");
            if (sign_1 == -1)
                std::swap(first, second);
            if (link_count + 2 > 5 * VX_MAXDIM)
                throw std::invalid_argument("Too many difference bounds");
            links[++link_count] = {first, second, '>', lo};
            links[++link_count] = {first, second, '<', hi};
            rows.push_back(difference_row(second, first, -lo));
            rows.push_back(difference_row(first, second, hi));
        }
        else if (type == 3)
        {
            /* 3 nnz (alternative consequence coefficient)* rhs
               represents sum(coefficient_k*z_k) <= rhs. */
            int nonzeros = 0;
            input >> nonzeros;
            if (nonzeros < 1 || nonzeros > variables)
                throw std::invalid_argument("Invalid sparse DDT row size");
            LinCons row;
            std::set<int> used;
            for (int k = 0; k < nonzeros; ++k)
            {
                int alternative = 0, consequence = 0;
                double coefficient = 0.0;
                input >> alternative >> consequence >> coefficient;
                const int index = variable(alternative, consequence, dimension);
                if (!input || !std::isfinite(coefficient) || coefficient == 0.0 ||
                    !used.insert(index).second)
                    throw std::invalid_argument("Invalid sparse DDT row term");
                row.index.push_back(index);
                row.coefficient.push_back(coefficient);
            }
            input >> row.b;
            require_complete(input);
            if (!std::isfinite(row.b))
                throw std::invalid_argument("Sparse DDT row has a non-finite rhs");
            rows.push_back(row);
        }
        else
        {
            throw std::invalid_argument("Unsupported DDT statement type");
        }
    }

    if (std::find(seen.begin() + 1, seen.end(), false) != seen.end())
        throw std::invalid_argument("Missing variable intervals");
    links[0].head = link_count;
}
} // namespace

LinCons difference_row(int head, int tail, double rhs)
{
    LinCons row;
    row.index = {head, tail};
    row.coefficient = {1.0, -1.0};
    row.b = rhs;
    return row;
}

void read_delta(const char *filename, std::vector<int> &dimension,
                d_row probability_lower, d_row probability_upper,
                link_set probability_links, std::vector<LinCons> &probability_rows,
                d_row value_lower, d_row value_upper,
                link_set value_links, std::vector<LinCons> &value_rows)
{
    std::ifstream input(filename);
    if (!input)
        throw std::invalid_argument("Cannot open DDT input");
    for (int i = 0; i < 4; ++i)
        next_line(input);

    std::istringstream dimensions(next_line(input));
    int alternatives = 0, first = 0, second = 0;
    dimensions >> alternatives >> first >> second;
    require_complete(dimensions);
    if (alternatives != 2)
        throw std::invalid_argument("DDT input must contain exactly two alternatives");
    if (first < 1 || second < 1 || first > VX_MAXDIM || second > VX_MAXDIM ||
        first + second > VX_MAXDIM)
        throw std::invalid_argument("DDT dimension exceeds VX_MAXDIM");
    dimension = {2, first, second};

    next_line(input);
    next_line(input);
    read_base(input, dimension, probability_lower, probability_upper,
              probability_links, probability_rows);
    next_line(input);
    read_base(input, dimension, value_lower, value_upper, value_links, value_rows);

    for (int alternative = 1; alternative <= 2; ++alternative)
    {
        double lo = 0.0, hi = 0.0;
        for (int consequence = 1; consequence <= dimension[alternative]; ++consequence)
        {
            const int index = variable(alternative, consequence, dimension);
            lo += probability_lower[index];
            hi += probability_upper[index];
        }
        if (lo > 1.0 + 1e-10 || hi < 1.0 - 1e-10)
            throw std::invalid_argument("Probability intervals cannot sum to one");
    }
}
