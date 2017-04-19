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
            
            // Whether there is clipping
            bool clip;
            
            // Size of the chromosome of the alignment
            Base length;
            
            void *b;
            void *h;
        };

        class Data : public Alignment
        {
            friend class ParserSAM;
            
            public:
            
                bool nextCigar(Locus &l, bool &spliced);

                inline void *b() const { return _b; }
                inline void *h() const { return _h; }

            private:
            
                mutable int _i, _n;

                void *_b;
                void *_h;
        };
        
        typedef std::function<void (Data &, const Info &)> Functor;
        
        /*
         * In order to improve the efficiency, not everything is computed. Set the last
         * argument to true will force it to happen.
         */

        static void parse(const FileName &, Functor, bool details = false);
    };
}

#endif
