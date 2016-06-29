#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <klib/khash.h>
#include "data/standard.hpp"
#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/writer_sam.hpp"

extern Anaquin::FileName BedRef();

namespace Anaquin
{
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
        enum CoverageMethod
        {
            Mean,
            Median,
            ReadCount,
        };

        struct Stats
        {
            // Intervals for synthetic
            ID2Intervals syn;
            
            // Intervals for genome
            ID2Intervals gen;
            
            // Raw coverage
            CoverageTool::Stats cov;
            
            // Calculated coverage for synthetic
            Coverage synC;
            
            // Calculated coverage for the query (eg: chr21)
            Coverage genC;
            
            Counts n_syn = 0;
            Counts n_gen = 0;

            /*
             * Fraction required to subsample in chrT. This works because chrT is a short
             * chromosome and almost certianly will have the highest coverage.
             */
            
            inline Proportion sample() const { return genC / synC; }
        };

        template <typename Options> static Stats stats(const FileName &file, const Options &o)
        {
            const auto &r = Standard::instance().r_var;
            
            o.analyze(file);
            
            Stats stats;
            
            // Reference regions for synthetic and genome
            auto inters = r.dInters();
            
            // There must be at least synthetic and genome
            assert(inters.size() >= 2);
            
            // Statistics for alignments mapped to the selected regions
            const auto rr = CoverageTool::stats(file, inters);

            for (const auto &i : rr.hist)
            {
                if (Standard::isSynthetic(i.first))
                {
                    stats.n_syn += i.second;
                }
                else if (Standard::isGenomic(i.first))
                {
                    stats.n_gen += i.second;
                }
            }
            
            if (!stats.n_syn)
            {
                throw std::runtime_error("Failed to find alignments for synthetic");
            }
            else if (!stats.n_gen)
            {
                throw std::runtime_error("Failed to find alignments for the genome");
            }
            
            o.info(toString(stats.n_syn) + " alignments to synthetic");
            o.info(toString(stats.n_gen) + " alignments to genome");
            o.info(toString(stats.n_syn + stats.n_gen) + " alignments in total");

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
            
            o.info(toString(stats.syn.size()) + " reference synthetic regions");
            o.info(toString(stats.gen.size()) + " reference genomic regions");
            
            const auto ss = stats.syn.stats();
            const auto gs = stats.gen.stats();

            assert(ss.mean && gs.mean);

            o.info("Calculating coverage for the synthetic and genome");
            
            /*
             * Now we have the data, we'll need to compare the coverage and determine the fraction that
             * the synthetic alignments needs to be sampled.
             */
            
            switch (o.method)
            {
                case Mean:
                {
                    stats.synC = ss.mean;
                    stats.genC = gs.mean;
                    break;
                }

                case Median:
                {
                    stats.synC = ss.p50;
                    stats.genC = gs.p50;
                    break;
                }

                case ReadCount:
                {
                    stats.synC = stats.n_syn;
                    stats.genC = stats.n_gen;
                    break;
                }
            }
            
            assert(stats.synC && stats.genC);
            
            o.info("Synthetic coverage: " + toString(stats.synC));
            o.info("Genomic coverage: " + toString(stats.genC));

            return stats;
        }

        template <typename Stats, typename Options> static Coverage meth2cov(const Stats &stats, const Options &o)
        {
            switch (o.method)
            {
                case Mean:      { return stats.mean; }
                case Median:    { return stats.p50;  }
                case ReadCount: { throw "????"; }
            }
        }
        
        struct SampleStats
        {
            // Sequence coverage after subsampling
            Coverage cov;
            
            // Number of reads after subsampling
            Counts reads = 0;
        };
        
        template <typename Options> static SampleStats sample(const FileName &src,
                                                              const FileName &dst,
                                                              Proportion prop,
                                                              ID2Intervals &inters,
                                                              const Options &o)
        {
            assert(prop >= 0 && prop <= 1.0);
            assert(!src.empty() && !dst.empty());
         
            o.info("Subsampling: " + toString(prop));
            
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
            
            ParserSAM::parse(src, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
            {
                if (info.p.i && !(info.p.i % 1000000))
                {
                    o.wait(std::to_string(info.p.i));
                }
                
                const auto isGen = !Standard::isSynthetic(x.cID);
                
                // This is the key, randomly write the reads with certain probability
                if (isGen || sampler.select(x.name))
                {
                    if (x.mapped && !isGen)
                    {
                        const auto m = inters.at(x.cID).contains(x.l);
                        
                        if (m)
                        {
                            m->map(x.l);
                        }
                    }

                    if (!isGen)
                    {
                        r.reads++;
                    }
                    
                    // Print the SAM line
                    writ.write(x);
                }
            }, true);

            writ.close();

            // Calculate the coverage after subsampling
            r.cov = meth2cov(inters.stats(), o);
            
            return r;
        }
        
        template <typename Options> static void report(const FileName &file, const Options &o)
        {
            auto meth2Str = [&]()
            {
                switch (o.method)
                {
                    case Mean:      { return "Mean";   }
                    case Median:    { return "Median"; }
                    case ReadCount: { return "Reads";  }
                }
            };
            
            o.info(meth2Str());
            
            const auto sampled = "VarSubsample_sampled.sam";
            
            // Statistics before sampling
            const auto before = Subsampler::stats(file, o);
            
            if (before.genC > before.synC)
            {
                throw std::runtime_error("Coverage for the genome is higher than the synthetic chromosome. Unexpected because the genome should be much wider.");
            }

            // Genomic coverage remains unchanged
            const auto g_after = before.genC;
            
            const auto &r = Standard::instance().r_var;

            // Reference regions for synthetic chromosomes
            auto inters = r.dIntersSyn();
            
            // Subsample the alignments
            const auto samp = Subsampler::sample(file, o.work + "/" + sampled, before.sample(), inters, o);
            
            o.info("Coverage (after): " + toString(samp.cov));
            o.info("Coverage (after): " + toString(g_after));

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

            /*
             * Generating VarSubsample_summary.stats
             */
            
            const auto summary = "VarSubsample Output Results\n\n"
                                 "-------VarSubsample Output\n\n"
                                 "       Reference regions: %1%\n"
                                 "       User generated alignment: %2%\n\n"
                                 "-------Reference regions\n\n"
                                 "       Synthetic regions: %3%\n"
                                 "       Genomic regions:   %4%\n\n"
                                 "       Method: %5%\n\n"
                                 "-------User alignments (before subsampling)\n\n"
                                 "       Synthetic: %6%\n"
                                 "       Genome:    %7%\n\n"
                                 "-------User alignments (after subsampling)\n\n"
                                 "       Synthetic: %8%\n"
                                 "       Genome:    %9%\n\n"
                                 "-------Before subsampling\n\n"
                                 "       Synthetic coverage: %10%\n"
                                 "       Genome coverage:    %11%\n\n"
                                 "-------After subsampling\n\n"
                                 "       Synthetic coverage: %12%\n"
                                 "       Genome coverage:    %13%\n";

            o.generate("VarSubsample_summary.stats");
            o.writer->open("VarSubsample_summary.stats");
            o.writer->write((boost::format(summary) % BedRef()
                                                    % file
                                                    % before.syn.countInters()
                                                    % before.gen.countInters()
                                                    % meth2Str()
                                                    % before.n_syn
                                                    % before.n_gen
                                                    % samp.reads
                                                    % before.n_gen
                                                    % before.genC
                                                    % before.synC
                                                    % before.genC
                                                    % samp.cov).str());
            o.writer->close();
        }
    };
}

#endif