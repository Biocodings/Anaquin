#include "VarQuin/v_unique.hpp"
#include "parsers/parser_fa.hpp"

using namespace Anaquin;

void VUnique::unique(const FileName &file)
{
//    std::map<SequinID, Sequence> s2s;
//
//    // Unique sequin names;
//    std::set<SequinID> uniqs;
//    
//    ParserFA::parse(file, [&](const ParserFA::Data &x, const ParserProgress &)
//    {
//        if (isSubstr(x.id, "_R") || isSubstr(x.id, "_V"))
//        {
//            uniqs.insert(noLast(x.id, "_"));
//            s2s[noLast(x.id, "_")] = x.seq;
//        }
//    });
//
//    A_ASSERT(!uniqs.empty());
//    
//    for (const auto &uniq : uniqs)
//    {
//        const auto R  = uniq + "_R";
//        const auto V  = uniq + "_R";
//        auto &R_ = s2s[R];
//        auto &V_ = s2s[V];
//
//        A_ASSERT(R_.size() != V_.size());
//        
//        // Breakpoint relative to the reference
//        Base b = 0;
//        
//        Base i1 = 0;
//        Base i2 = 0;
//        char c1, c2;
//        
//        /*
//         * Let's find the breakpoint
//         */
//        
//        for (; R_[b] != V_[b]; b++) {}
//
//        /*
//         * What's this variant? SNP? Indel?
//         */
//        
//        if R_
//        
//        
//    }
}
