#ifndef SS_SUMMARY_HPP
#define SS_SUMMARY_HPP

#include <string>

namespace SS
{
    template <typename T> void summary(const T &t)
    {
        t.summary();
    }
}

#endif