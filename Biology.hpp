#ifndef AS_BIOLOGY_HPP
#define AS_BIOLOGY_HPP

template <typename Iter, typename F> void exonsToIntrons(const Iter &exons, F f)
{
    for (auto i = 0; i < exons.size(); i++)
    {
        if (i)
        {
            f(exons[i - 1].end, exons[i].start);
        }
    }
}

#endif