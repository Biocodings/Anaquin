#include <map>
#include "data/tokens.hpp"
#include "data/experiment.hpp"

using namespace Anaquin;

// Eg: A1,A2,A3,B1,B2,B3
void Experiment::addNames(const std::string &str)
{
    Tokens::split(str, ",", _names);
    assert(_samples.empty() || _names.size() == _samples.size());
}

void Experiment::addFactors(const std::string &str)
{
    std::vector<std::string> toks;
    Tokens::split(str, ",", toks);

    // Factor mapping
    std::map<std::string, Factor> m;
    
    for (const auto &tok : toks)
    {
        if (!m.count(tok))
        {
            m[tok] = static_cast<unsigned>(m.size());
        }

        _samples.push_back(m[tok]);
        _factors.insert(m[tok]);
    }
    
    assert(_samples.size() == toks.size());
    assert(_names.empty() || _samples.size() == _names.size());
}

std::vector<std::size_t> Experiment::cond(Factor f) const
{
    std::vector<std::size_t> r;
    
    for (auto i = 0; i < _samples.size(); i++)
    {
        if (_samples[i] == f)
        {
            // It's the index we're looking for
            r.push_back(i);
        }
    }

    return r;
}
