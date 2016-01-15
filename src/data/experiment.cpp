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
        
        _factors.push_back(m[tok]);
    }
    
    assert(_factors.size() == toks.size());
}