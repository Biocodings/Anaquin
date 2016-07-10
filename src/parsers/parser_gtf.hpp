#ifndef PARSER_GTF_HPP
#define PARSER_GTF_HPP

#include <assert.h>
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "stats/analyzer.hpp"

extern bool __hack__;

namespace Anaquin
{
    struct ParserGTF
    {
        struct Data
        {
            inline bool overlap(const Locus &l) const
            {
                return this->l.overlap(l);
            }
            
            operator Locus() const { return l; }
            
            void operator=(const Data &f)
            {
                l    = f.l;
                type = f.type;
                tID  = f.tID;
                gID  = f.gID;
            }
            
            ChrID cID;
            
            // Forward or reverse strand?
            Strand str;
            
            // The location of the feature relative to the chromosome
            Locus l;
            
            RNAFeature type;
            
            // Empty if the information is unavailable
            GeneID gID;
            
            // Empty if the information is unavailable
            TransID tID;
            
            // Available only for Cufflink
            double fpkm = NAN;
        };
        
        template <typename F> static void parse(const Reader &r, F f)
        {
            std::map<std::string, RNAFeature> mapper =
            {
                { "exon",       Exon },
                { "gene",       Gene },
                { "transcript", Transcript }
            };
            
            std::string line;
            Data x;
            
            /*
             * Fields must be tab-separated. Also, all but the final field in each feature line must contain a value;
             * "empty" columns should be denoted with a '.'. Please refer to the online documentation for more details.
             *
             *    1. seqname
             *    2. sourcre
             *    3. feature
             *    4. start
             *    5. end
             *    6. score
             *    7. strand
             *    8. frame
             *    9. attribute
             */
            
            ParserProgress p;
            
            std::vector<std::string> opts;
            std::vector<std::string> toks;
            std::vector<std::string> nameVal;
            
            while (r.nextLine(line))
            {
                if (__hack__)
                {
                    std::cout << line << std::endl;                    
                }
                
                p.i++;
                boost::split(toks, line, boost::is_any_of("\t"));
                
                // Empty line? Unknown feature such as mRNA?
                if (toks.size() == 1 || !mapper.count(toks[2]))
                {
                    continue;
                }
                
                x.cID  = toks[0];
                x.type = mapper[toks[2]];

                x.l.start = stoi(toks[3]);
                x.l.end   = stoi(toks[4]);

                if (toks[6] != "+" || toks[6] != "-")
                {
                    throw std::runtime_error("File: " + r.src() + ". Invalid strand: [" + toks[6] + "]. Line: " + line );
                }

                x.str = toks[6] == "+" ? Strand::Forward : Strand::Backward;
                
                /*
                 * Eg: "gene_id "R_5_3"; transcript_id "R_5_3_R";"
                 */
                
                boost::split(opts, toks[8], boost::is_any_of(";"));
                
                for (auto option : opts)
                {
                    if (!option.empty())
                    {
                        boost::trim(option);
                        boost::split(nameVal, option, boost::is_any_of(" "));
                        
                        if (nameVal.size() == 2)
                        {
                            // Make sure the silly characters are removed
                            nameVal[1].erase(std::remove(nameVal[1].begin(), nameVal[1].end(), '\"'), nameVal[1].end());
                            
                            if (nameVal[0] == "gene_id")
                            {
                                x.gID = nameVal[1];
                            }
                            else if (nameVal[0] == "transcript_id")
                            {
                                x.tID = nameVal[1];
                            }
                            else if (nameVal[0] == "FPKM")
                            {
                                x.fpkm = stod(nameVal[1]);
                            }
                        }
                    }
                }

                f(x, line, p);
            }            
        }
    };
}

#endif