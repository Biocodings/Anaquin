/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef PARSER_DIFF_HPP
#define PARSER_DIFF_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "data/convert.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserDiff
    {
        typedef enum
        {
            ChrID,
            GeneID,
            IsoformID,
            Sample1,
            Sample2,
            LogFold,
            LogFoldSE,
            PValue,
            QValue,
            Mean
        } Field;

        typedef DiffTest Data;
        
        static bool isDiff(const Reader &r)
        {
            std::string line;
            std::vector<Tokens::Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                
                if (toks.size() == 10      &&
                    toks[0] == "ChrID"     &&
                    toks[1] == "GeneID"    &&
                    toks[2] == "IsoformID" &&
                    toks[3] == "Sample1"   &&
                    toks[4] == "Sample2"   &&
                    toks[5] == "LogFold"   &&
                    toks[6] == "LogFoldSE" &&
                    toks[7] == "PValue"    &&
                    toks[8] == "QValue"    &&
                    toks[9] == "Average")
                {
                    return true;
                }
            }
            
            return false;
        }
        
        template <typename F> static void parse(const FileName &file, F f)
        {
            const auto &r = Standard::instance().r_trans;
            
            Reader rr(file);
            ParserProgress p;
            std::vector<Tokens::Token> toks;

            std::string line;
            
            while (rr.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                Data x;

                if (p.i)
                {
                    x.gID    = toks[Field::GeneID]    != "-" ? toks[Field::GeneID]    : "";
                    x.iID    = toks[Field::IsoformID] != "-" ? toks[Field::IsoformID] : "";
                    x.p      = ss2ld(toks[Field::PValue]);
                    x.q      = ss2ld(toks[Field::QValue]);
                    x.mean   = s2d(toks[Field::Mean]);
                    x.logF_   = s2d(toks[Field::LogFold]);
                    x.logFSE = s2d(toks[Field::LogFoldSE]);
                    x.samp1  = s2d(toks[Field::Sample1]);
                    x.samp2  = s2d(toks[Field::Sample2]);

                    // Eg: DESeq2 wouldn't give the chromoname name
                    x.cID = toks[Field::ChrID];
                    
                    if (x.cID == "-")
                    {
                        if (r.findGene(ChrIS, x.gID) || r.findTrans(ChrIS, x.iID))
                        {
                            x.cID = ChrIS;
                        }
                        else
                        {
                            x.cID = Geno;
                        }
                    }

                    f(x, p);
                }

                p.i++;
            }
        }
        
        static bool isIsoform(const Reader &r)
        {
            std::vector<Tokens::Token> toks;
            
            if (r.nextTokens(toks, "\t"))
            {
                return toks[1] == "IsoID";
            }
            
            return true;
        }
    };
}

#endif
