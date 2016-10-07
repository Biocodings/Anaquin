#ifndef PARSER_TSV_HPP
#define PARSER_TSV_HPP

#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserTSV
    {
        enum TSVField
        {
            Contig,
            KMerLength,
            ContigLength,
            ColoredKMers,
            PropColoredKMers,
            ModeKMer,
            KMerObs,
            Total,
            ProportionKMers
        };
        
        struct TSV
        {
            ContigID id;
            
            Depth dep;
            
            // Unnormalized k-mer observations
            KMers kmer;
            
            // K-mer length
            Base klen;
        };
        
        typedef std::function<void(const TSV &, const ParserProgress &)> Functor;
        
        static void parse(const Reader &r, Functor f)
        {
            /*
             * 1. Contig name
             * 2. K-mer length
             * 3. Contig length in k-mers
             * 4. Colored k-mers
             * 5. Proportion of colored k-mers in contig
             * 6. Mode k-mer coverage depth
             * 7. K-mer observations
             * 8. Total
             * 9. Proportion of k-mer observations in sample
             */
            
            TSV t;
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            while (r.nextLine(line))
            {
                if (p.i++ == 0)
                {
                    continue;
                }
                
                Tokens::split(line, "\t", toks);
                
                t.id   = toks[TSVField::Contig];
                t.dep  = stoi(toks[TSVField::ModeKMer]);
                t.kmer = stoi(toks[TSVField::KMerObs]);
                t.klen = stoi(toks[TSVField::ContigLength]);
                
                f(t, p);
            }
        }
    };
}

#endif
