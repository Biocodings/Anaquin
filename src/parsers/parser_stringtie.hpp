#ifndef PARSER_STRINGTIE_HPP
#define PARSER_STRINGTIE_HPP

#include "data/tokens.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserStringTie
    {
        struct Data
        {
            ChrID cID;
            
            IsoformID id;
            
            Locus l;
            
            // The expression for the abundance
            FPKM abund;
        };
        
        template <typename F> static void parse(const FileName &file, F f)
        {
            Reader r(file);
            
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                f(toks, p);
                p.i++;
            }
        }

        enum GeneField
        {
            GeneID,
            GeneName,
            GeneStrand,
            GeneStart,
            GeneEnd,
            GeneLength,
            GeneCoverage,
            GeneFPKM,
            GeneTPM
        };
        
        enum CTabField
        {
            TabID,
            TabChrID,
            TabStrand,
            TabStart,
            TabEnd,
            TabTName,
            TabNumExons,
            TabLength,
            TabGeneID,
            TabGeneName,
            TabCov,
            TabFPKM
        };

        static void parseGenes(const FileName &file, std::function<void (const Data &, const ParserProgress &)> f)
        {
            Data t;
            
            ParserStringTie::parse(file, [&](const std::vector<Tokens::Token> &toks, const ParserProgress &p)
            {
                if (p.i)
                {
                    t.id = toks[GeneField::GeneID];
                    
                    if (Standard::instance().r_trans.findGene(ChrT, t.id))
                    {
                        t.cID = ChrT;
                    }
                    else
                    {
                        t.cID = Endo;
                    }
                    
                    t.l     = Locus(stod(toks[GeneField::GeneStart]), stod(toks[GeneField::GeneEnd]));
                    t.abund = stod(toks[GeneField::GeneFPKM]);

                    f(t, p);
                }
            });
        }
            
        static void parseCTab(const FileName &file, std::function<void (const Data &, const ParserProgress &)> f)
        {
            Data t;
                      
            ParserStringTie::parse(file, [&](const std::vector<Tokens::Token> &toks, const ParserProgress &p)
            {
                if (p.i)
                {
                    t.id    = toks[CTabField::TabTName];
                    t.cID   = toks[CTabField::TabChrID];
                    t.l     = Locus(stod(toks[CTabField::TabStart]), stod(toks[CTabField::TabEnd]));
                    t.abund = stod(toks[CTabField::TabFPKM]);

                    f(t, p);
                }
            });
        }
    };
}

#endif