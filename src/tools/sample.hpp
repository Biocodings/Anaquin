#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <klib/khash.h>
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
            ArithAverage,
            Maximum,
            Median,
            Percentile75,
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
            
            /*
             * Fraction required to subsample in chrT. This works because chrT is a short
             * chromosome and almost certianly will have the highest coverage.
             */
            
            inline Proportion sample() const { return genC / synC; }
        };
        
        /*
         * Histogram manipulations and operations
         */
        
        template <typename T> static Counts sums(const std::map<T, Counts> &m)
        {
            Counts c = 0;
            
            for (const auto &i : m)
            {
                if (i.second == 0)
                {
                    c++;
                }
                else
                {
                    c += i.second;
                }
            }
            
            return c;
        }
        
        template <typename Options> static Stats stats(const FileName &file, const Options &o)
        {
            const auto &r = Standard::instance().r_var;
            
            o.analyze(file);
            
            Stats stats;
            
            // Intervals for synthetic and genome
            auto inters = r.inters();
            
            assert(inters.size() >= 2);
            
            // Statistics for all reads mapped to the selected regions
            const auto rr = CoverageTool::stats__(file, inters);

            // Number of alignments to the synthetic
            Counts n_syn = 0;
            
            // Number of alignments to the genome
            Counts n_gen = 0;
            
            for (const auto &i : rr.hist)
            {
                if (Standard::isSynthetic(i.first))
                {
                    n_syn += i.second;
                }
                else if (Standard::isGenomic(i.first))
                {
                    n_gen += i.second;
                }
            }
            
            if (!n_syn)
            {
                throw std::runtime_error("Failed to find any alignment for synthetic");
            }
            else if (!n_gen)
            {
                throw std::runtime_error("Failed to find any alignment for the genome");
            }
            
            o.info(toString(n_syn) + " alignments to synthetic");
            o.info(toString(n_gen) + " alignments to genome");
            o.info(toString(n_syn + n_gen) + " alignments in total");
            o.info(toString(stats.cov.inters.size()) + " intervals generated");

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
            
            const auto ss = stats.syn.stats();
            const auto gs = stats.gen.stats();

            assert(ss.mean && gs.mean);

            o.info("Calculating coverage for the synthetic and genome");
            
            /*
             * Now we have the data, we'll need to compare the coverage and determine the fraction that
             * the synthetic genome needs to be sampled.
             */
            
            switch (o.method)
            {
                case ArithAverage:
                {
                    stats.synC = ss.mean;
                    stats.genC = gs.mean;
                    break;
                }

                case Maximum:
                {
                    stats.synC = ss.max;
                    stats.genC = gs.max;
                    break;
                }

                case Median:
                {
                    stats.synC = ss.p50;
                    stats.genC = gs.p50;
                    break;
                }

                case Percentile75:
                {
                    stats.synC = ss.p75;
                    stats.genC = gs.p75;
                    break;
                }
            }
            
            assert(stats.synC && stats.genC);
            
            if (stats.genC > stats.synC)
            {
                throw std::runtime_error("Coverage for the genome is higher than the synthetic chromosome. Unexpected because the genome should be much wider.");
            }
            
            return stats;
        }

        /*
         * Perform subsampling to a proportion of input synthetic reads
         */

        template <typename Options> static void sample(const FileName &src,
                                                       const FileName &dst,
                                                       Proportion prop,
                                                       const Options &o)
        {
            assert(prop >= 0 && prop <= 1.0);
            assert(!src.empty() && !dst.empty());
            
            /*
             * Subsampling alignments. It's expected that coverage would roughly match between
             * the genome and synthetic chromosome.
             */
            
            o.info("Sampling the alignment");
            
            WriterSAM writer;
            writer.open(dst);
            
            if (prop == 0.0)
            {
                o.warn("Sampling fraction is zero. This could be an error in the inputs.");
            }
            else if (prop == 1.0)
            {
                o.warn("Sampling fraction is one. This could be an error in the inputs.");
            }
            
            SamplingTool sampler(1 - prop);
            
            ParserSAM::parse(src, [&](const Alignment &align, const ParserSAM::Info &info)
            {
                if (!align.i && info.p.i && !(info.p.i % 1000000))
                {
                    o.wait(std::to_string(info.p.i));
                }
                
                const auto *b = reinterpret_cast<bam1_t *>(info.data);
                const auto *h = reinterpret_cast<bam_hdr_t *>(info.header);
                
                if (!align.i)
                {
                    // This is the key, randomly write the reads with certain probability
                    if (!Standard::isSynthetic(align.cID) || sampler.select(bam_get_qname(b)))
                    {
                        writer.write(h, b);
                    }
                }
            });
            
            writer.close();
        }
        
        template <typename Options> static void report(const FileName &file, const Options &o)
        {
            const auto &r = Standard::instance().r_var;
            
            auto meth2Str = [&]()
            {
                switch (o.method)
                {
                    case Percentile75: { return "75th";    }
                    case ArithAverage: { return "Mean";    }
                    case Maximum:      { return "Maximum"; }
                    case Median:       { return "Median";  }
                }
            };
            
            o.info(meth2Str());
            
            const auto sampled = "VarSubsample_sampled.sam";
            
            // Statistics before alignment
            const auto before = Subsampler::stats(file, o);
            
            // Subsample the alignments
            Subsampler::sample(file, o.work + "/" + sampled, before.sample(), o);
            
            // Statistics after alignment
            const auto after = Subsampler::stats(o.work + "/" + sampled, o);

            o.info("Proportion (after): " + toString(after.sample()));

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
            
            CoverageTool::bedGraph(after.syn, pre);

            /*
             * Generating VarSubsample_summary.stats
             */
            
            const auto summary = "VarSubsample Output Results\n\n"
                                 "-------VarSubsample Output\n\n"
                                 "       Reference sequin regions: %1%\n"
                                 "       User generated alignment: %2%\n\n"
                                 "-------Reference regions\n\n"
                                 "       Synthetic regions: %3%\n"
                                 "       Genomic regions:   %4%\n\n"
                                 "       Method: %5%\n\n"
                                 "-------User alignments (before subsampling)\n\n"
                                 "       Genome:    %6%\n\n"
                                 "       Synthetic: %7%\n\n"
                                 "-------User alignments (after subsampling)\n\n"
                                 "       Genome:    %8%\n\n"
                                 "       Synthetic: %9\n\n"
                                 "-------Before subsampling\n\n"
                                 "       Genome coverage:    %10%\n"
                                 "       Synthetic coverage: %11%\n\n"
                                 "-------After subsampling\n\n"
                                 "       Genome coverage:    %12%\n"
                                 "       Synthetic coverage: %13%\n\n";

            o.generate("VarSubsample_summary.stats");
            o.writer->open("VarSubsample_summary.stats");
            o.writer->write((boost::format(summary) % BedRef()
                                                    % file
                                                    % r.sInters().size()
                                                    % r.gInters().size()
                                                    % meth2Str()
                                                    % before.n_gen
                                                    % before.n_syn
                                                    % after.n_gen
                                                    % after.n_syn
                                                    % before.synC
                                                    % before.genC
                                                    % after.synC
                                                    % after.genC).str());
            o.writer->close();
        }
    };
}

#endif