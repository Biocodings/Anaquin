#ifndef SS_C_HPP
#define SS_C_HPP

#include <vector>
#include <algorithm>
#include <ss/type.hpp>
#include <ss/types.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> class C
        {
            public:
                typedef Real value_type;

                typename std::vector<T>::const_iterator begin() const { return _data.begin(); }
                typename std::vector<T>::const_iterator  end()  const { return _data.end();   }

                /*
                 * Simulates giving a list of values to a column vector in R.
                 * For example, c(1,2,3,4,5,6) in R can be simulated by
                 *
                 *     c({1,2,3,4,5,6}) or c("1,2,3,4,5,6")
                 */

                template <typename Iter> C operator()(const Iter &x)
                {
                    _data.clear();
                    _data.resize(x.size());

                    std::copy(x.begin(), x.end(), _data.begin());

                    return *this;
                }

                inline std::size_t size() const { return _data.size(); }

            //private:
                //std::vector<Type<T>> _data;
                std::vector<T> _data;
        };
    }

    static Internal::C<Real> c;
}

#endif