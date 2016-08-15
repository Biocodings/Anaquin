#include <klib/khash.h>
#include "tools/random.hpp"
#include "VarQuin/v_sample2.hpp"
#include "writers/writer_sam.hpp"
#include "readers/reader_bam.hpp"

// Defined in main.cpp
extern Anaquin::FileName BedRef();

using namespace Anaquin;

typedef std::map<ChrID, std::map<Locus, Proportion>> NormFactors;

static ReaderBam::Stats sample(const FileName &file,
                               const std::map<ChrID, std::set<Locus>> &c2l,
                               const NormFactors &norms,
                               const VSample2::Options &o)
{
    typedef std::map<ChrID, std::map<Locus, std::shared_ptr<RandomSelection>>> Selection;
    
    Selection select;
    
    for (const auto &i : norms)
    {
        for (const auto &j : i.second)
        {
            assert(j.second > 0 && j.second <= 1.0);
            select[i.first][j.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - j.second));
        }
    }

    assert(select.size() == norms.size());

    /*
     * Subsampling alignments. It's expected that coverage would roughly match between
     * the genome and synthetic chromosome.
     */
    
    o.info("Sampling the alignments");
    
    WriterSAM writer;
    writer.openTerm();
    
    // Statistics within the sampling region (ReaderBam gives everything)
    auto sampled = ReaderBam::reg2Inters(c2l);

    const auto r = ReaderBam::stats(file, c2l, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logInfo(std::to_string(info.p.i));
        }
        
        // Must write to the new file
        bool mustWrite = false;
        
        // Write to the new file because it's been sampled
        bool sampledWrite = false;
        
        // Always write anything other than sequins
        mustWrite = !x.mapped || !Standard::isSynthetic(x.cID);
        
        if (x.mapped && !mustWrite)
        {
            // Within sampling regions?
            const auto inRegion = sampled.count(x.cID) ? (bool) sampled.at(x.cID).contains(x.l) : false;
            
            if (inRegion && select[x.cID][x.l]->select(x.name))
            {
                sampledWrite = true;
                sampled.at(x.cID).overlap(x.l)->map(x.l);
            }
        }
        
        if (mustWrite || sampledWrite)
        {
            if (sampledWrite)
            {
                o.logInfo("Sampled " + x.name);
            }
            
            writer.write(x);
        }

        return (mustWrite || sampledWrite);
    });
    
    return r;
}

VSample2::Stats VSample2::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    // Intervals required for subsampling
    std::map<ChrID, std::set<Locus>> c2v;
    
    const auto bStats = ReaderBam::stats(file, c2v, [&](const ParserSAM::Data &, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        return true;
    });
    
    VSample2::Stats stats;

    // Normalization for each region
    NormFactors norms;
    
    // For each chromosome...
    for (auto &chr : c2v)
    {
        // For each region...
        for (auto &reg : chr.second)
        {
            // Genomic statistics within the region
            const auto gStats = bStats.gen.at(chr.first).find(reg.key())->stats();
            
            // Synthetic statistics within the region
            const auto sStats = bStats.gen.at(chr.first).find(reg.key())->stats();
            
            o.info("Calculating coverage for the synthetic and genome");
            
            /*
             * Now we have the data, we'll need to compare the coverage and determine the fraction that
             * the synthetic alignments needs to be sampled.
             */

            Coverage synC = NAN;
            Coverage genC = NAN;

            switch (o.meth)
            {
                case VSample2::Method::Mean:
                {
                    synC = sStats.mean;
                    genC = gStats.mean;
                    break;
                }

                case VSample2::Method::Median:
                {
                    synC = sStats.p50;
                    genC = gStats.p50;
                    break;
                }

                case VSample2::Method::Prop:
                {
                    synC = stats.samp.syn;
                    genC = stats.samp.gen;
                    break;
                }

                case VSample2::Method::Reads:
                {
                    throw std::runtime_error("Not implemented");
                    break;
                }
            }
            
            const auto norm = std::min(genC / synC, 1.0);
            
            if (norm == 1.0)
            {
                o.warn((boost::format("Normalization factor is 1 for %1%:%2%-%3%") % chr.first
                                                                                   % reg.start
                                                                                   % reg.end).str());
            }
            
            norms[chr.first][reg] = norm;
        }
    }
    
    /*
     * Now, we have the normalization factors. We can proceed with subsampling
     */
    
    sample(file, c2v, norms, o);
    
    return stats;
}

void VSample2::report(const FileName &file, const Options &o)
{
    stats(file, o);
}