#ifndef F_EXPRESS_HPP
#define F_EXPRESS_HPP

#include "fusion/f_discover.hpp"

namespace Anaquin
{
    struct FExpress
    {
        typedef FDiscover::Options Options;

        struct Stats : public FDiscover::Stats
        {
            Limit ss;
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif