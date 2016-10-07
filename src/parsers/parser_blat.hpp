#ifndef PARSER_BLAT_HPP
#define PARSER_BLAT_HPP

#include <functional>
#include "data/data.hpp"
#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserBlat
    {
        enum Field
        {
            PSL_Matches,
            PSL_MisMatches,
            PSL_RepMatches,
            PSL_NCount,
            PSL_QGap_Count,
            PSL_QGap_Bases,
            PSL_TGap_Count,
            PSL_TGap_Bases,
            PSL_Strand,
            PSL_QName,
            PSL_QSize,
            PSL_QStart,
            PSL_QEnd,
            PSL_TName,
            PSL_TSize,
            PSL_TStart,
            PSL_TEnd,
            PSL_Block_Count,
            PSL_Block_Sizes,
            PSL_Q_Starts,
            PSL_T_Starts
        };
        
        struct Data
        {
            // Target sequence name
            std::string tName;
            
            // Query sequence name
            std::string qName;
            
            // Alignment start position in target
            Base tStart;
            
            // Alignment end position in target
            Base tEnd;
            
            // Target sequence size
            Base tSize;
            
            // Alignment start position in query
            Base qStart;
            
            // Alignment end position in query
            Base qEnd;
            
            // Query sequence size
            Base qSize;
            
            // Number of inserts in query
            Counts qGapCount;
            
            // Number of bases inserted into query
            Base qGap;
            
            // Number of inserts in target
            Counts tGapCount;
            
            // Number of bases inserted into target
            Base tGap;
            
            // Number of matching bases
            Base match;
            
            // Number of mismatching bases
            Base mismatch;
        };
        
        template <typename F> static void parse(const Reader &r, F f)
        {
            protectParse("PSL format", [&]()
                         {
                             ParserProgress p;
                             
                             std::string line;
                             std::vector<std::string> toks;
                             
                             Data l;
                             
                             while (r.nextLine(line))
                             {
                                 if (p.i++ <= 3)
                                 {
                                     continue;
                                 }
                                 
                                 Tokens::split(line, "\t", toks);
                                 
                                 if (toks.size() != 21)
                                 {
                                     throw std::runtime_error("Invalid line: " + line);
                                 }
                                 
                                 l.qName     = toks[PSL_QName];
                                 l.tName     = toks[PSL_TName];
                                 
                                 l.tStart    = stoi(toks[PSL_TStart]);
                                 l.tEnd      = stoi(toks[PSL_TEnd]);
                                 l.tSize     = stoi(toks[PSL_TSize]);
                                 
                                 l.qStart    = stoi(toks[PSL_QStart]);
                                 l.qEnd      = stoi(toks[PSL_QEnd]);
                                 l.qSize     = stoi(toks[PSL_QSize]);
                                 
                                 l.match     = stoi(toks[PSL_Matches]);
                                 l.mismatch  = stoi(toks[PSL_MisMatches]);
                                 
                                 l.qGap      = stoi(toks[PSL_QGap_Bases]);
                                 l.tGap      = stoi(toks[PSL_TGap_Bases]);
                                 
                                 l.qGapCount = stoi(toks[PSL_QGap_Count]);
                                 l.tGapCount = stoi(toks[PSL_TGap_Count]);
                                 
                                 f(l, p);
                             }
                         });
        }
    };
}

#endif
