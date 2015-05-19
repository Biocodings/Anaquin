#ifndef SS_PROB_HPP
#define SS_PROB_HPP

namespace SS
{
    struct P
    {
        P(double p = 0.0) : p(p) {}

        inline double comp() const { return (1 - p); }

        inline void operator=(const P &p) { this->p = p.p; }

        operator double() const { return p; }

        double p;
    };
}

#endif
