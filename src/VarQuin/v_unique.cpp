#include "VarQuin/v_unique.hpp"
#include "parsers/parser_fa.hpp"

using namespace Anaquin;

void VUnique::unique(const FileName &file)
{
    std::map<SequinID, Sequence> s2s;

    // Unique sequin names;
    std::set<SequinID> uniqs;
    
    ParserFA::parse(file, [&](const ParserFA::Data &x, const ParserProgress &)
    {
        if (isSubstr(x.id, "_R") || isSubstr(x.id, "_V"))
        {
            uniqs.insert(noLast(x.id, "_"));
            s2s[noLast(x.id, "_")] = x.seq;
        }
    });
    
    
    
    
}
