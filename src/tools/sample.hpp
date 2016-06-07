#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <klib/khash.h>
#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/writer_sam.hpp"

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
            // Coverage statistics for chrT
            Interval::Stats chrT;
            
            // Coverage statistics for endogenous (eg: chr21)
            Interval::Stats endo;
            
            // Raw coverage
            CoverageTool::Stats cov;
            
            // Calculated coverage for chrT
            Coverage chrTC;
            
            // Calculated coverage for the query (eg: chr21)
            Coverage endoC;
            
            /*
             * Fraction required to subsample in chrT. This works because chrT is a short
             * chromosome and almost certianly will have the highest coverage.
             */
            
            inline Proportion sample() const { return endoC / chrTC; }
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
        
        static bool checkAlign(const ChrID &genoID, const ChrID &id, const Locus &l)
        {
            const auto &r = Standard::instance().r_var;
            
            if (id == ChrT)
            {
                return r.match(l, MatchRule::Contains);
            }
            else if (id == genoID)
            {
                return r.findGeno(genoID, l);
            }
            
            return false;
        }
        
        struct StatsImpl
        {
            // Returns the chromosome for the genomic region
            virtual ChrID genoID() const = 0;

            // Whether the genomic region should be selected
            virtual bool shouldGenomic(const ChrID &, const Locus &) const = 0;

            // Whether the synthetic region should be selected
            virtual bool shouldSynthetic(const ChrID &, const Locus &) const = 0;
        };
        
        struct ReportImpl
        {
            // File name for the summary statistics
            virtual FileName summary() const = 0;
            
            // File name for the bedgraph before subsampling
            virtual FileName beforeBG() const = 0;
            
            // File Name for the bedgraph after subsampling
            virtual FileName afterBG() const = 0;
          
            // File name for the subsampled alignment
            virtual FileName sampled() const = 0;
            
            // Number of synthetic sequins
            virtual Counts countSeqs() const = 0;
            
            // Number of genomic intervals
            virtual Counts countInters() const = 0;
        };

        template <typename Options> static Stats stats(const FileName &file,
                                                       const Options &o,
                                                       const StatsImpl &impl)
        {
            const auto genoID = impl.genoID();
            
            o.info("Genome: " + genoID);
            o.analyze(file);
            
            Stats stats;
            
            stats.cov = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &p)
            {
                if (!align.i && !(p.i % 1000000))
                {
                    o.wait(std::to_string(p.i));
                }
                
                return checkAlign(genoID, align.cID, align.l);
            });
            
            if (!stats.cov.hist.count(ChrT))
            {
                throw std::runtime_error("Failed to find any alignment for " + ChrT);
            }
            else if (!stats.cov.hist.count(genoID))
            {
                throw std::runtime_error("Failed to find any alignment for: " + genoID);
            }
            
            o.info(toString(sums(stats.cov.hist)) + " alignments in total");
            o.info(toString(stats.cov.hist.at(ChrT)) + " alignments to chrT");
            o.info(toString(stats.cov.hist.at(genoID)) + " alignments to " + genoID);
            o.info(toString(stats.cov.inters.size()) + " intervals generated");
            
            o.info("Generating statistics for: " + ChrT);
            stats.chrT = stats.cov.inters.find(ChrT)->stats([&](const ChrID &id, Base i, Base j, Coverage cov)
            {
                return impl.shouldSynthetic(id, Locus(i,j));
            });
            
            o.info("Generating statistics for: " + genoID);
            stats.endo = stats.cov.inters.find(genoID)->stats([&](const ChrID &id, Base i, Base j, Coverage cov)
            {
                return impl.shouldGenomic(id, Locus(i,j));
            });
            
            assert(stats.chrT.mean && stats.endo.mean);
            
            o.info("Calculating coverage for " + ChrT + " and " + genoID);
            
            /*
             * Now we have the data, we'll need to compare the coverage and determine what fraction that
             * the synthetic chromosome needs to be subsampled.
             */
            
            switch (o.method)
            {
                case ArithAverage:
                {
                    stats.chrTC = stats.chrT.mean;
                    stats.endoC = stats.endo.mean;
                    break;
                }
                    
                case Maximum:
                {
                    stats.chrTC = stats.chrT.max;
                    stats.endoC = stats.endo.max;
                    break;
                }
                    
                case Median:
                {
                    stats.chrTC = stats.chrT.p50;
                    stats.endoC = stats.endo.p50;
                    break;
                }
                    
                case Percentile75:
                {
                    stats.chrTC = stats.chrT.p75;
                    stats.endoC = stats.endo.p75;
                    break;
                }
            }
            
            assert(stats.chrTC && stats.endoC);
            
            if (stats.endoC >= stats.chrTC)
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
                if (!align.i && !(info.p.i % 1000000))
                {
                    o.wait(std::to_string(info.p.i));
                }
                
                const auto *b = reinterpret_cast<bam1_t *>(info.data);
                const auto *h = reinterpret_cast<bam_hdr_t *>(info.header);
                
                if (!align.i)
                {
                    /*
                     * This is the key, randomly write the reads with certain probability
                     */
                    
                    if (align.cID != ChrT || sampler.select(bam_get_qname(b)))
                    {
                        writer.write(h, b);
                    }
                }
            });
            
            writer.close();
        }
        
        template <typename Options> static void report(const FileName &file,
                                                       const Options &o,
                                                       const StatsImpl &si,
                                                       const ReportImpl &ri)
        {
            const auto genoID = si.genoID();
            
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
            
            const auto sampled = ri.sampled();
            
            // Statistics before alignment
            const auto before = Subsampler::stats(file, o, si);
            
            // Subsample the alignment
            Subsampler::sample(file, o.work + "/" + sampled, before.sample(), o);
            
            // Statistics after alignment
            const auto after = Subsampler::stats(o.work + "/" + sampled, o, si);

            o.info("Proportion (after): " + toString(after.sample()));

            /*
             * Generating bedgraph before subsampling
             */
            
            auto pre = CoverageTool::CoverageBedGraphOptions();
            
            pre.writer = o.writer;
            pre.file   = ri.beforeBG();
            
            CoverageTool::bedGraph(before.cov, pre, [&](const ChrID &id, Base i, Base j, Coverage)
            {
                return checkAlign(genoID, id, Locus(i, j));
            });
            
            /*
             * Generating statistics after subsampling
             */
            
            auto post = CoverageTool::CoverageBedGraphOptions();
            
            post.writer = o.writer;
            post.file   = ri.afterBG();
            
            CoverageTool::bedGraph(after.cov, post, [&](const ChrID &id, Base i, Base j, Coverage)
            {
                return checkAlign(genoID, id, Locus(i, j));
            });
            
            /*
             * Generating summary statistics
             */
            
            const auto summary = "Summary for input: %1%\n\n"
                                 "   ***\n"
                                 "   *** Proportion of alignments mapped to the synthetic and genome\n"
                                 "   ***\n\n"
                                 "   Unmapped:  %2% aligns\n"
                                 "   Synthetic: %3% aligns\n"
                                 "   Genome:    %4% aligns\n\n"
                                 "   ***\n"
                                 "   *** Reference annotation\n"
                                 "   ***\n\n"
                                 "   File: %5%\n\n"
                                 "   Synehtic: %6% sequins\n\n"
                                 "   Genome:   %7% intervals\n\n"
                                 "   ***                                                      \n"
                                 "   *** How the coverage is calculated? Possibilities:       \n"
                                 "   ***                                                      \n"
                                 "   ***    - Mean                                            \n"
                                 "   ***    - Quartile                                        \n"
                                 "   ***    - Median                                          \n"
                                 "   ***    - Maximum                                         \n"
                                 "   ***                                                      \n"
                                 "   *** Please refer to the online documentation for details \n"
                                 "   ***                                                      \n\n"
                                 "   Method: %8%\n\n"
                                 "   ***                               \n"
                                 "   *** Statistics before subsampling \n"
                                 "   ***                               \n\n"
                                 "   Coverage (Synthetic): %9%\n"
                                 "   Coverage (Genome):    %10%\n\n"
                                 "   ***                              \n"
                                 "   *** Statistics after subsampling \n"
                                 "   ***                              \n\n"
                                 "   Coverage (Synthetic): %11%\n"
                                 "   Coverage (Genome):    %12%\n";
            
            o.generate(ri.summary());
            o.writer->open(ri.summary());
            o.writer->write((boost::format(summary) % file
                                                    % before.cov.n_unmap
                                                    % before.cov.n_syn
                                                    % before.cov.n_gen
                                                    % o.rAnnot
                                                    % ri.countSeqs()
                                                    % ri.countInters()
                                                    % meth2Str()
                                                    //% sums(before.cov.hist)
                                                    % before.chrTC
                                                    % before.endoC
                                                    //% sums(after.cov.hist)
                                                    % after.chrTC
                                                    % after.endoC).str());
            o.writer->close();
        }
    };
}

#endif