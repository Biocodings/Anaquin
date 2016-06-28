#ifndef PARSER_SAM_HPP
#define PARSER_SAM_HPP

#include "data/alignment.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserSAM
    {
        struct Info
        {
            ParserProgress p;

            // Whether this is a multi-alignment
            bool multi;
            
            // Whether there is insertion
            bool ins;
            
            // Whether there is deletion
            bool del;
            
            // Whether there is skipped region
            bool skip;
            
            // Size of the chromosome of the alignment
            Base length;
            
            // Internal data representation
            void *data;
            
            // Internal data representation
            void *header;
        };

        class Data : public Alignment
        {
            friend class ParserSAM;
            
            public:
            
                bool nextCigar(Locus &l, bool &spliced);

            private:
            
                mutable int _i, _n;
            
                void *data;
                void *head;
        };
        
        typedef std::function<void (Data &, const Info &)> Functor;
        static void parse(const FileName &, Functor);
    };
}

#endif