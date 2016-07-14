#include <klib/khash.h>
#include "VarQuin/v_sample.hpp"
#include "writers/writer_sam.hpp"

#define DEBUG_VSAMPLE

// Defined in main.cpp
extern Anaquin::FileName BedRef();

using namespace Anaquin;

class SamplingTool
{
    public:
    
        SamplingTool(double prob) : _prob(prob)
        {
            assert(prob >= 0.0);
            _seed = rand();
        }
    
        inline bool select(const std::string &hash) const
        {
            const uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(hash.c_str()) ^ _seed);
            return ((double)(k&0xffffff) / 0x1000000 >= _prob);
        }
    
    private:
    
        // Random seed
        int _seed;
    
        // The probability of selection
        const double _prob;
};

struct Subsampler
{
    template <typename Options> static VSample::Stats stats(const FileName &file, const Options &o)
    {
        const auto &r = Standard::instance().r_var;
        
        o.analyze(file);
        
        VSample::Stats stats;
        
        // Reference regions for synthetic and genome
        auto inters = r.dInters();
        
        // There must be at least synthetic and genome
        assert(inters.size() >= 2);
        
        // Statistics for alignments mapped to the selected regions
        const auto rr = CoverageTool::stats(file, inters);
        
        // Number of alignments within the sampling regions for synthetic
        auto s_syn = 0;
        
        // Number of alignments within the sampling regions for genome
        auto s_gen = 0;
        
        // For each chromsome within the sampling regions...
        for (const auto &i : rr.hist)
        {
            if (Standard::isSynthetic(i.first))
            {
                s_syn += i.second;
            }
            else
            {
                s_gen += i.second;
            }
        }
        
        stats.tot.syn  = rr.n_syn;
        stats.tot.gen  = rr.n_gen;
        stats.samp.syn = s_syn;
        stats.samp.gen = s_gen;
        
        assert(stats.samp.syn <= stats.tot.syn);
        assert(stats.samp.gen <= stats.tot.gen);
        
        if (!stats.samp.syn)
        {
            throw std::runtime_error("No alignment found within synthetic regions");
        }
        
        if (!stats.samp.gen)
        {
            throw std::runtime_error("No alignment found within genomic regions");
        }
        
        o.logInfo(toString(stats.tot.syn) + " total alignments to synthetic");
        o.logInfo(toString(stats.tot.gen) + " total alignments to genome");
        o.logInfo(toString(stats.tot.syn + stats.tot.gen) + " total alignments in total");

        o.logInfo(toString(stats.samp.syn) + " alignments (within sampling regions) to synthetic");
        o.logInfo(toString(stats.samp.gen) + " alignments (within sampling regions)to genome");
        o.logInfo(toString(stats.samp.syn + stats.samp.gen) + " alignments (within sampling regions) in total");
        
        /*
         * Calculate statistics for both synthetic and genome
         */
        
        for (const auto &i : inters)
        {
            const auto &cID = i.first;
            
            if (Standard::isSynthetic(cID))
            {
                stats.syn.add(cID, i.second);
            }
            else if (Standard::isGenomic(cID))
            {
                stats.gen.add(cID, i.second);
            }
        }
        
        assert(!stats.syn.empty());
        assert(!stats.gen.empty());
        
        o.info(std::to_string(stats.syn.size()) + " reference synthetic regions");
        o.info(std::to_string(stats.gen.size()) + " reference genomic regions");
        
        const auto ss = stats.syn.stats();
        const auto gs = stats.gen.stats();
        
        assert(ss.mean && gs.mean);
        
        o.info("Calculating coverage for the synthetic and genome");
        
        /*
         * Now we have the data, we'll need to compare the coverage and determine the fraction that
         * the synthetic alignments needs to be sampled.
         */
        
        switch (o.meth)
        {
            case VSample::Method::Mean:
            {
                stats.synC = ss.mean;
                stats.genC = gs.mean;
                break;
            }
                
            case VSample::Method::Median:
            {
                stats.synC = ss.p50;
                stats.genC = gs.p50;
                break;
            }
                
            case VSample::Method::Prop:
            {
                stats.synC = stats.samp.syn;
                stats.genC = stats.samp.gen;
                break;
            }

            case VSample::Method::Reads:
            {
                stats.synC = s_syn;
                stats.genC = s_gen;
                break;
            }
        }
        
        assert(stats.synC && stats.genC);
        
        o.info("Synthetic coverage: " + std::to_string(stats.synC));
        o.info("Genomic coverage: " + std::to_string(stats.genC));
        
        return stats;
    }
    
    struct SampleStats
    {
        // Sequence coverage after subsampling (synthetic only)
        Coverage cov = NAN;
        
        SynGenAligns tot;
        SynGenAligns samp;
    };
    
    template <typename Options> static SampleStats sample(const FileName &src,
                                                          const FileName &dst,
                                                          Proportion prop,
                                                          ID2Intervals &inters,
                                                          const Options &o)
    {
        assert(prop > 0 && prop <= 1.0);
        assert(!src.empty() && !dst.empty());
        
        o.info("Subsampling: " + std::to_string(prop));
        
        /*
         * Subsampling alignments. It's expected that coverage would roughly match between
         * the genome and synthetic chromosome.
         */
        
        o.info("Sampling the alignments");
        
        WriterSAM writ;
        writ.openTerm();
        
        if (prop == 0.0)
        {
            o.warn("Sampling proportion is zero. This could be an error in the inputs.");
        }
        else if (prop == 1.0)
        {
            o.warn("Sampling proportion is one. This could be an error in the inputs.");
        }
        
        SamplingTool sampler(1.0 - prop);
        SampleStats r;
        
        /*
         * Subsample the input alignment file
         */
        
        ParserSAM::parse(src, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
        {
            if (info.p.i && !(info.p.i % 1000000))
            {
                o.logInfo(std::to_string(info.p.i));
            }
            
            const auto isSyn = Standard::isSynthetic(x.cID);
            const auto isGen = !isSyn;
            
            const auto shouldWrite = !x.mapped || !isSyn;
            
            // This is the key, randomly write the reads with certain probability
            if (shouldWrite || sampler.select(x.name))
            {
                // Within sampling regions?
                const auto inRegion = inters.count(x.cID) ? (bool) inters.at(x.cID).contains(x.l) : false;

                if (x.mapped && isSyn)
                {
                    assert(Standard::isSynthetic(x.cID));
                    const auto m = inters.at(x.cID).contains(x.l);
                    
                    if (m)
                    {
                        m->map(x.l);
                    }
                }
                
                if (x.cID != "*")
                {
                    if (isSyn) { r.tot.syn++; }
                    if (isGen) { r.tot.gen++; }
                    
                    if (inRegion)
                    {
                        if (isSyn) { r.samp.syn++; }
                        if (isGen) { r.samp.gen++; }
                    }
                }
                
                if (!shouldWrite)
                {
                    o.logInfo("Sampled " + x.name);
                }
                
                // Print the SAM line
                writ.write(x);
            }
        }, true);
        
        writ.close();
        
        /*
         * Calculating sequence coverage on the in-silico chromosome after subsampling
         */
        
        switch (o.meth)
        {
            case VSample::Method::Prop:   { r.cov = r.tot.syn;           break; }
            case VSample::Method::Mean:   { r.cov = inters.stats().mean; break; }
            case VSample::Method::Median: { r.cov = inters.stats().p50;  break; }
                
            /*
             * Total synthetic alignments within the sampling regions after subsampling.
             * For example, we might subsample from 100 million alignments down to 2 million
             * alignments (this implies the genomic regions have 2 million alignments).
             */

            case VSample::Method::Reads:  { r.cov = r.samp.syn; break; }
        }
        
        assert(!isnan(r.cov));
        return r;
    }

    template <typename Options> static void report(const FileName &file, const Options &o)
    {
        auto meth2Str = [&]()
        {
            switch (o.meth)
            {
                case VSample::Method::Mean:   { return "Mean";       }
                case VSample::Method::Median: { return "Median";     }
                case VSample::Method::Reads:  { return "Reads";      }
                case VSample::Method::Prop:   { return "Proportion"; }
            }
        };
        
        o.info(meth2Str());
        
        const auto sampled = "VarSubsample_sampled.sam";
        
        // Statistics before sampling
        const auto before = Subsampler::stats(file, o);

        switch (o.meth)
        {
            case VSample::Method::Mean:
            case VSample::Method::Median:
            {
                if (before.genC > before.synC)
                {
                    throw std::runtime_error("Coverage for the genome is higher than the synthetic chromosome. Unexpected because the genome should be much wider.");
                }
                
                break;
            }
                
            case VSample::Method::Reads:
            {
                if (before.genC > before.synC)
                {
                    throw std::runtime_error("Anaquin is not able to subsample because there are more alignments in the genomic regions than the in-silico regions. Genomic alignments: " + std::to_string(before.genC) + ". Synthetic alignments: " + std::to_string(before.synC));
                }
                
                break;
            }

            default: { break; }
        }
        
        const auto &r = Standard::instance().r_var;
        
        // Reference regions for synthetic chromosomes
        auto inters = r.dIntersSyn();
        
        // Proportion of reads be sampled
        Proportion norm;
        
        switch (o.meth)
        {
            case VSample::Method::Prop:
            {
                norm = o.p;
                break;
            }

            default:
            {
                norm = before.genC / before.synC;
                break;
            }
        }

        assert(norm > 0 && norm < 1.0);
        
        o.info("Normalization: " + std::to_string(norm));

        o.info("Coverage (before): " + std::to_string(before.synC));
        o.info("Coverage (before): " + std::to_string(before.genC));
        
        // Subsample the alignments
        const auto after = Subsampler::sample(file, o.work + "/" + sampled, norm, inters, o);
        
        o.info("Coverage (after): " + std::to_string(after.cov));
        o.info("Coverage (after): " + std::to_string(before.genC));
        
#ifdef DEBUG_VSAMPLE
        /*
         * Generating bedgraph before subsampling (only synthetic)
         */
        
        auto pre = CoverageTool::CoverageBedGraphOptions();
        
        pre.writer = o.writer;
        pre.file   = "VarSubsample_before.bedgraph";
        
        CoverageTool::bedGraph(before.syn, pre);
        
        /*
         * Generating bedgraph after subsampling (only synthetic)
         */
        
        auto post = CoverageTool::CoverageBedGraphOptions();
        
        post.writer = o.writer;
        post.file   = "VarSubsample_after.bedgraph";
        
        CoverageTool::bedGraph(inters, pre);
#endif
        
        /*
         * Reproduce %6%: samtools view sampled.bam | cut -f3,6 | grep chrT | grep -v '*' | wc
         * Reprodcue %7%: samtools view sampled.bam | cut -f3,6 | grep -v chrT | grep -v '*' | wc
         */
        
        /*
         * Generating VarSubsample_summary.stats
         */
        
        const auto summary = "-------VarSubsample Summary Statistics\n\n"
                             "       Reference annotation file: %1%\n"
                             "       User alignment file: %2%\n\n"
                             "-------Reference regions\n\n"
                             "       Synthetic regions: %3%\n"
                             "       Genomic regions:   %4%\n\n"
                             "       Method: %5%\n\n"
                             "-------Total alignments (before subsampling)\n\n"
                             "       Synthetic: %6%\n"
                             "       Genome:    %7%\n\n"
                             "-------Total alignments (after subsampling)\n\n"
                             "       Synthetic: %8%\n"
                             "       Genome:    %9%\n\n"
                             "-------Alignments within sampling regions (before subsampling)\n\n"
                             "       Synthetic: %10%\n"
                             "       Genome:    %11%\n\n"
                             "-------Alignments within sampling regions (after subsampling)\n\n"
                             "       Synthetic: %12%\n"
                             "       Genome:    %13%\n\n"
                             "       Normalization: %14%\n\n"
                             "-------Before subsampling (within sampling regions)\n\n"
                             "       Synthetic coverage: %15%\n"
                             "       Genome coverage:    %16%\n\n"
                             "-------After subsampling (within sampling regions)\n\n"
                             "       Synthetic coverage: %17%\n"
                             "       Genome coverage:    %18%\n";
        
        o.generate("VarSubsample_summary.stats");
        o.writer->open("VarSubsample_summary.stats");
        o.writer->write((boost::format(summary) % BedRef()                 // 1
                                                % file                     // 2
                                                % before.syn.countInters() // 3
                                                % before.gen.countInters() // 4
                                                % meth2Str()               // 5
                                                % before.tot.syn           // 6
                                                % before.tot.gen           // 7
                                                % after.tot.syn            // 8
                                                % after.tot.gen            // 9
                                                % before.samp.syn          // 10
                                                % before.samp.gen          // 11
                                                % after.samp.syn           // 12
                                                % before.samp.gen          // 13
                                                % norm                     // 14
                                                % before.synC              // 15
                                                % before.genC              // 16
                                                % after.cov                // 17
                                                % before.genC              // 18
                         ).str());
        o.writer->close();
    }
};

VSample::Stats VSample::stats(const FileName &file, const Options &o)
{
    return Subsampler::stats(file, o);
}

void VSample::report(const FileName &file, const Options &o)
{
    Subsampler::report(file, o);
}