#ifndef PARSER_BAM_HPP
#define PARSER_BAM_HPP

#include "data/alignment.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserBAM
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
            friend class ParserBAM;
            
            public:
            
                bool nextCigar(Locus &l, bool &spliced);

                /*
                 * Optional fields
                 */
            
                // Eg: B7_591:6:155:12:674
                void lName();

                // Segment sequence (optional)
                void lSeq();

                // ASCII of base QUALity (optional)
                void lQual();

                inline void *b() const { return _b; }
                inline void *h() const { return _h; }

            private:
            
                mutable int _i, _n;

                void *_b;
                void *_h;
        };
        
        typedef std::function<void (Data &, const Info &)> Functor;
        
        static std::map<ChrID, Base> header(const FileName &);
        
        /*
         * In order to improve the efficiency, not everything is computed. Set the last
         * argument to true will force it to happen.
         */

        static void parse(const FileName &, Functor, bool details = false);
    };
}

#endif
