#ifndef SS_DATA_HPP
#define SS_DATA_HPP

#include <ss/c.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> struct Data
        {
            static Data frame(const C<T> &c1, const C<T> &c2)
            {
                Data d;
                
                d.cs.push_back(c1);
                d.cs.push_back(c2);

                return d;
            }

            std::vector<C<T>> cs;
        };
    }

    static Internal::Data<Real> data;
}

#endif