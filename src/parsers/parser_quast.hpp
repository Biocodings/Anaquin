#ifndef PARSER_QUAST_HPP
#define PARSER_QUAST_HPP

#include "data/types.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    class ParserQuast
    {
        public:
        
            struct GenomeData
            {
                // Eg: MG_29
                SequinID id;
            
                // Total length
                Base total;
            
                // Covered length
                Base covered;
            };
        
            struct ContigData
            {
                // Eg: MG_29
                SequinID id;
            
                std::vector<ContigID> contigs;
            };
        
        /*
         * contigs_reports/alignments_Contigs.tsv
         *
         *   CME003.v013_M3_G   contig-2042000000_1793_nucleotides
         *   CME003.v013_MG_28  contig-1936000000_2957_nucleotides  contig-21984000000_278_nucleotides
         */
        
        static void parseContigs(const Reader &r, std::function<void(const ContigData &, const ParserProgress &)> f)
        {
            ContigData x;
            ParserProgress p;
            
            Line line;
            std::vector<Tokens::Token> toks;
            
            while (r.nextLine(line))
            {
                p.i++;
                boost::trim(line);
                
                /*
                 * Eg: CME003.v013_MG_28  contig-1936000000_2957_nucleotides  contig-21984000000_278_nucleotides
                 */
                
                Tokens::split(line, "\t", toks);
                assert(toks.size() >= 2);
                
                // Eg: MG_28
                x.id = parseSequin(toks[0]);

                x.contigs.clear();

                for (auto i = 1; i < toks.size(); i++)
                {
                    x.contigs.push_back(toks[i]);
                }

                f(x, p);
            }
        }
        
        /*
         *  genome_stats/genome_info.txt
         *
         *  reference chromosomes:
         *      CME003.v013_MG_29 (total length: 2974 bp, maximal covered length: 0 bp)
         *      CME003.v013_M3_G (total length: 1824 bp, maximal covered length: 1793 bp)
         */

        static void parseGenome(const Reader &r, std::function<void(const GenomeData &, const ParserProgress &)> f)
        {
            protectParse("genome_info.txt", [&]()
            {
                GenomeData x;
                ParserProgress p;
                
                Line line;
                std::vector<Tokens::Token> toks;
                
                while (r.nextLine(line))
                {
                    // Skip: "referebce chromosome:"
                    if (p.i++ == 0)
                    {
                        continue;
                    }
                    else if (line.find("total length") == std::string::npos)
                    {
                        continue;
                    }
                    
                    boost::trim(line);

                    /*
                     * Eg: CME003.v013_MG_29 (total length: 2974 bp, maximal covered length: 100 bp)
                     */
                    
                    Tokens::split(line, " ", toks);
                    
                    // Eg: 2974
                    x.total = stod(toks[3]);
                    
                    // Eg: 100
                    x.covered = stod(toks[8]);
                    
                    // Eg: CME003.v013_MG_29
                    Tokens::split(toks[0], "_", toks);
                    
                    // Eg: MG_29
                    x.id = toks[1] + "_" + toks[2];
                    
                    f(x, p);
                }
            });
        }
        
        private:
        
            // Eg: CME003.v013_MG_29
            static std::string parseSequin(const std::string &x)
            {
                std::vector<Tokens::Token> toks;
                Tokens::split(x, "_", toks);
            
                // Eg: MG_29
                return toks[1] + "_" + toks[2];
            }
    };
}

#endif