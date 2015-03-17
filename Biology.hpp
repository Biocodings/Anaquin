#ifndef AS_BIOLOGY_HPP
#define AS_BIOLOGY_HPP

/*
 * Intron is a sequence within a gene that is removed by RNA splicing. This function calculates
 * the gaps between the given exons and assume them be introns (also known as spliced junctions).
 */

template <typename Iter, typename F> void extractIntrons(const Iter &exons, F f)
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