#ifndef SS_TYPE_HPP
#define SS_TYPE_HPP

#include <string>

namespace SS
{
    namespace Internal
    {
        template <typename T> class Type
        {
            public:
                operator T() const { return _t; }
            
                void operator=(const std::string &s)
                {
                    _s = s;
                    _type = Numeric;
                }

                void operator=(T t)
                {
                    _t = t;
                    _type = String;
                }
            
            private:

                enum ActiveType
                {
                    Numeric,
                    String,
                };

                ActiveType _type;

                // Numerical data-type in R
                T _t;

                // Another common data-type in R
                std::string _s;
        };
    }
}

#endif