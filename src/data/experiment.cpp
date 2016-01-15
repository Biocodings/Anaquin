#include <map>
#include "data/tokens.hpp"
#include "data/experiment.hpp"

using namespace Anaquin;

Experiment::Experiment(const std::string &str)
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

        _reps.push_back(m[tok]);
        _factors.insert(m[tok]);
    }
    
    assert(_reps.size() == toks.size());
}

std::vector<std::size_t> Experiment::factor(Factor f) const
{
    std::vector<std::size_t> r;
    
    for (auto i = 0; i < _reps.size(); i++)
    {
        if (_reps[i] == f)
        {
            // It's the index we're looking for
            r.push_back(i);
        }
    }

    return r;
}
