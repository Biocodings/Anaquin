#ifndef SS_INTERNAL_RANK_HPP
#define SS_INTERNAL_RANK_HPP

#include <ss/internal/sort.hpp>

namespace SS
{
    namespace Internal
    {
        enum class TieMethod
        {
            TieAverage,
            TieMax,
            TieMin,
        };
        
        template <typename T> inline void cummin(const T &x, T &r)
        {
            auto min = x.front();
            auto xi  = x.begin();
            auto ri  = r.begin();
            
            while (xi != x.end())
            {
                if (*xi < min)
                {
                    min = *xi;
                }

                *ri = min;

                xi++;
                ri++;
            }
        }

        template <typename T1, typename T2> inline void order(const T1 &x, T2 &r, bool increase=true)
        {
            Index i = 0;
            
            for (auto &j : r)
            {
                j = ++i;
            }
            
            permSort(x, r);
            
            if (!increase)
            {
                std::reverse(r.begin(), r.end());
            }
        }
   
        template <typename T1, typename T2> inline void rank(const T1 &t, T2 &r, TieMethod meth)
        {
            auto x1 = t;
            auto x2 = t;

            Counts i = 0;
            
            for (auto iter = x2.begin(); iter != x2.end(); iter++)
            {
                *iter = i++;
            }
            
            Internal::permSort(x1, x2);
            
            const auto n = t.size();

            for (int i = 0, j = 0; i < n; i = j+1)
            {
                j = i;
                
                while ((j < n-1) && (x1[j] == x1[j+1]))
                {
                    j++;
                }

                switch (meth)
                {
                    case TieMethod::TieAverage:
                    {
                        for (int k = i; k <= j; k++) { r[(Counts)x2[k]] = (i+j+2) / 2.; }
                        break;
                    }

                    case TieMethod::TieMax:
                    {
                        for (int k = i; k <= j; k++) { r[(Counts)x2[k]] = j+1; }
                        break;
                    }

                    case TieMethod::TieMin:
                    {
                        for (int k = i; k <= j; k++) { r[(Counts)x2[k]] = i+1; }
                        break;
                    }
                }
            }
        }
    }
}

#endif
