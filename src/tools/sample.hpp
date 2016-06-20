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
            // Coverage statistics for synthetic
            Intervals<> syn;
            
            // Coverage statistics for genome
            Intervals<> gen;
            
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
        
        struct StatsImpl
        {
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
        };

        template <typename Options> static Stats stats(const FileName &file,
                                                       const Options &o,
                                                       const StatsImpl &impl)
        {
            const auto &r = Standard::instance().r_var;
            
            o.analyze(file);
            
            Stats stats;
            
            // Intervals for synthetic and genome
            auto inters = r.inters();
            
            assert(inters.size() >= 2);
            
            // Does the read aligned within a region?
            auto inside = [&](const Alignment &x)
            {
                if (!inters.count(x.cID))
                {
                    return false;
                }
                
                // Does the read aligned within the reference regions?
                return (bool)inters[x.cID].contains(x.l);
            };
            
            /*
             * Collecting statistics for all alignments to the synthetic and genome
             */
            
            stats.cov = CoverageTool::stats(file, [&](const Alignment &x, const ParserProgress &p)
            {
                if (!x.i && p.i && !(p.i % 1000000))
                {
                    o.wait(std::to_string(p.i));
                }
                
                return inside(x);
            });
            
            // Number of alignments to the synthetic
            Counts n_syn = 0;
            
            // Number of alignments to the genome
            Counts n_gen = 0;
            
            for (const auto &i : stats.cov.hist)
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
            
            // There should be at least synthetic and genome
            assert(stats.cov.inters.size() >= 2);

            o.info(toString(n_syn) + " alignments to synthetic");
            o.info(toString(n_gen) + " alignments to genome");
            o.info(toString(n_syn + n_gen) + " alignments in total");
            o.info(toString(stats.cov.inters.size()) + " intervals generated");

            /*
             * Filter out to the genomic regions
             */
            
            o.info("Generating statistics for genome");

            for (const auto &i : stats.cov.inters.data())
            {
                const auto &cID = i.first;
                
                if (Standard::isSynthetic(cID))
                {
                    const auto inter = stats.cov.inters.find(i.first)->stats
                            ([&](const ChrID &id, Base i, Base j, Coverage cov)
                    {
                         return inters[cID].contains(Locus(i,j));
                    });
                    
                    //stats.syn.add(inter);
                }
                
                if (Standard::isGenomic(cID))
                {
                    const auto inter = stats.cov.inters.find(i.first)->stats
                            ([&](const ChrID &id, Base i, Base j, Coverage cov)
                    {
                            return inters[cID].contains(Locus(i,j));
                    });

                    //stats.gen.add(inter);
                }
            }
            
            const auto sstats = stats.syn.stats();
            const auto gstats = stats.gen.stats();

            assert(sstats.mean && gstats.mean);

            o.info("Calculating coverage for the synthetic and genome");
            
            /*
             * Now we have the data, we'll need to compare the coverage and determine what fraction that
             * the synthetic chromosome needs to be subsampled.
             */
            
            switch (o.method)
            {
                case ArithAverage:
                {
                    stats.synC = sstats.mean;
                    stats.genC = gstats.mean;
                    break;
                }
                    
                case Maximum:
                {
                    stats.synC = sstats.max;
                    stats.genC = gstats.max;
                    break;
                }
                    
                case Median:
                {
                    stats.synC = sstats.p50;
                    stats.genC = gstats.p50;
                    break;
                }
                    
                case Percentile75:
                {
                    stats.synC = sstats.p75;
                    stats.genC = gstats.p75;
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
            
            const auto sampled = ri.sampled();
            
            // Statistics before alignment
            const auto before = Subsampler::stats(file, o, si);
            
            // Subsample the alignment
            Subsampler::sample(file, o.work + "/" + sampled, before.sample(), o);
            
            // Statistics after alignment
            const auto after = Subsampler::stats(o.work + "/" + sampled, o, si);

            o.info("Proportion (after): " + toString(after.sample()));

            // Intervals for synthetic and genome
            auto inters = r.inters();
            
            // Does the read aligned within a region?
            auto inside = [&](const ChrID &cID, const Locus &l)
            {
                if (!inters.count(cID))
                {
                    return false;
                }
                
                // Does the read aligned within the reference regions?
                return (bool)inters[cID].contains(l);
            };

            /*
             * Generating bedgraph before subsampling
             */
            
            auto pre = CoverageTool::CoverageBedGraphOptions();
            
            pre.writer = o.writer;
            pre.file   = ri.beforeBG();
            
            CoverageTool::bedGraph(before.cov, pre, [&](const ChrID &cID, Base i, Base j, Coverage)
            {
                return inside(cID, Locus(i, j));
            });
            
            /*
             * Generating statistics after subsampling
             */
            
            auto post = CoverageTool::CoverageBedGraphOptions();
            
            post.writer = o.writer;
            post.file   = ri.afterBG();
            
            CoverageTool::bedGraph(after.cov, post, [&](const ChrID &cID, Base i, Base j, Coverage)
            {
                return inside(cID, Locus(i, j));
            });
            
            /*
             * Generating summary statistics
             */
            
            const auto summary = "VarSubsample Output Results\n\n"
                                 "-------VarSubsample Output\n\n"
                                 "       Reference sequin regions: %1%\n"
                                 "       User generated alignment: %2%\n\n"
                                 "-------Reference regions\n\n"
                                 "       Genome regions:   %3%\n"
                                 "       Synthetic regions: %4%\n\n"
                                 "       Method: %8%\n\n"            
                                 "-------User alignments (before subsampling)\n\n"
                                 "       Unmapped:  %5%\n"
                                 "       Genome:    %6%\n\n"
                                 "       Synthetic: %7%\n\n"
                                 "-------User alignments (after subsampling)\n\n"
                                 "       Unmapped:  %5%\n"
                                 "       Genome:    %6%\n\n"
                                 "       Synthetic: %7%\n\n"
                                 "-------Before subsampling\n\n"
                                 "       Genome coverage:    %9%\n"
                                 "       Synthetic coverage: %10%\n\n"
                                 "-------After subsampling\n\n"
                                 "       Genome coverage:    %11%\n"
                                 "       Synthetic coverage: %12%\n\n";
            
            o.generate(ri.summary());
            o.writer->open(ri.summary());
            o.writer->write((boost::format(summary) % o.rAnnot
                                                    % file
                                                    % "????" //ri.countInters()
                                                    % "????" //ri.countSeqs()
                                                    % before.cov.n_unmap
                                                    % before.cov.n_gen
                                                    % before.cov.n_syn
                                                    % meth2Str()
                                                    % before.synC
                                                    % before.genC
                                                    % after.synC
                                                    % after.genC).str());
            o.writer->close();
        }
    };
}

#endif