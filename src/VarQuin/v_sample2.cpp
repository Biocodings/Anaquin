#include "tools/random.hpp"
#include "VarQuin/v_sample2.hpp"
#include "writers/writer_sam.hpp"
#include "readers/reader_bam.hpp"

// Defined in main.cpp
extern Anaquin::FileName BedRef();

using namespace Anaquin;

typedef std::map<ChrID, std::map<Locus, Proportion>> NormFactors;

static ReaderBam::Stats sample(const FileName &file,
                               const NormFactors &norms,
                               const VSample2::Options &o)
{
    typedef std::map<ChrID, std::map<Locus, std::shared_ptr<RandomSelection>>> Selection;
    
    /*
     * Initalize independnet random generators for each of the sampling region
     */
    
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

    o.info("Sampling: " + file);
    
    WriterSAM writer;
    writer.openTerm();

    // Subsampling regions
    const auto sampled = Standard::instance().r_var.dInters();
    
    return ReaderBam::stats(file, sampled, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.logInfo(std::to_string(info.p.i));
        }
        
        auto shouldSampled = !x.mapped;
        
        if (x.mapped && !shouldSampled)
        {
            Interval *inter;
            
            // Within sampling regions?
            if (sampled.count(x.cID) && (inter = sampled.at(x.cID).contains(x.l)))
            {
                if (select[x.cID][inter->l()]->select(x.name))
                {
                    shouldSampled = true;
                }
            }
        }
        
        if (shouldSampled)
        {
            //writer.write(x);
            return true;
        }

        return false;
    });
}

template <typename Stats> Coverage stats2cov(const VSample2::Method meth, const Stats &stats)
{
    switch (meth)
    {
        case VSample2::Method::Mean:   { return stats.mean; }
        case VSample2::Method::Median: { return stats.p50;  }

        case VSample2::Method::Prop:
        {
            throw std::runtime_error("Not implemented");
            break;
        }

        case VSample2::Method::Reads:
        {
            throw std::runtime_error("Not implemented");
            break;
        }
    }
}

VSample2::Stats VSample2::analyze(const FileName &gen, const FileName &seq, const Options &o)
{
    o.analyze(gen);
    o.analyze(seq);

    const auto &r = Standard::instance().r_var;
    
    VSample2::Stats stats;
    
    // Regions to subsample
    const auto refs = r.dInters();
    
    A_ASSERT(!refs.empty(), "Empty reference sampling regions");
    
    // Checking genomic alignments before subsampling
    const auto gStats = ReaderBam::stats(gen, refs, [&](const ParserSAM::Data &, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        return true;
    });
    
    // Checking synthetic alignments before subsampling
    const auto sStats = ReaderBam::stats(seq, refs, [&](const ParserSAM::Data &, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        return true;
    });
    
    // Normalization for each region
    NormFactors norms;
    
    // For each chromosome...
    for (auto &i : refs)
    {
        const auto &cID = i.first;
        
        // For each region...
        for (auto &j : i.second.data())
        {
            const auto &l = j.second.l();
            
            // Genomic statistics within the region
            const auto gs = gStats.inters.at(cID).find(l.key())->stats();
            
            // Synthetic statistics within the region
            const auto ss = sStats.inters.at(cID).find(l.key())->stats();
            
            o.info("Calculating coverage for the synthetic and genome");
            
            /*
             * Now we have the data, we'll need to compare the coverage and determine the fraction that
             * the synthetic alignments needs to be sampled.
             */
            
            const auto synC = stats2cov(o.meth, ss);
            const auto genC = stats2cov(o.meth, gs);
            
            const auto norm = std::min(genC / synC, 1.0);
            
            if (norm == 1.0)
            {
                o.warn((boost::format("Normalization factor is 1 for %1%:%2%-%3%") % i.first
                                                                                   % l.start
                                                                                   % l.end).str());
            }
            else
            {
                o.logInfo((boost::format("Normalization for %1%-%2% - %3%") % i.first
                                                                            % l.start
                                                                            % l.end).str());
            }
            
            if (!genC) { stats.noGAlign++; }
            if (!synC) { stats.noSAlign++; }
            
            stats.c2v[cID][l].gen    = genC;
            stats.c2v[cID][l].before = synC;
            stats.c2v[cID][l].norm   = norms[i.first][l] = norm;
            
            stats.totLen += l.length();
        }
    }
    
    // Now, we have the normalization factors. We can proceed with subsampling.
    const auto ss = sample(seq, norms, o);
    
    auto n = 0;
    
    // For each chromosome...
    for (auto &i : refs)
    {
        // For each region...
        for (auto &j : i.second.data())
        {
            n++;
            const auto &l = j.second.l();
            stats.c2v[i.first][l].after = stats2cov(o.meth, j.second.stats());
        }
    }
    
    // Average sampling length
    stats.averLen = static_cast<Base>(stats.totLen / (double) n);
    
    return stats;
}

static void generateCSV(const FileName &file, const VSample2::Stats &stats, const VSample2::Options &o)
{
    o.generate(file);

    const auto format = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%");
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ChrID"
                                           % "Position"
                                           % "Start"
                                           % "End"
                                           % "Ref"
                                           % "Alt"
                                           % "Genome"
                                           % "Before"
                                           % "After"
                                           % "Norm"
                                           % "Type").str());

    for (const auto &i : stats.c2v)
    {
        
    }
    
    o.writer->close();
}

static void generateSummary(const FileName &file,
                            const FileName &gen,
                            const FileName &seq,
                            const VSample2::Stats &stats,
                            const VSample2::Options &o)
{
    o.generate(file);
    
    auto meth2Str = [&]()
    {
        switch (o.meth)
        {
            case VSample2::Method::Mean:   { return "Mean";       }
            case VSample2::Method::Median: { return "Median";     }
            case VSample2::Method::Reads:  { return "Reads";      }
            case VSample2::Method::Prop:   { return "Proportion"; }
        }
    };

    const auto summary = "-------VarSubsample Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Alignment file (genome):  %2%\n"
                         "       Alignment file (sequins): %3%\n\n"
                         "-------Reference regions\n\n"
                         "       Variant regions: %4% regions\n\n"
                         "       Method: %5%\n\n"
                         "-------Total alignments (before subsampling)\n\n"
                         "       Synthetic: %6%\n"
                         "       Genome:    %7%\n\n"
                         "-------Total alignments (after subsampling)\n\n"
                         "       Genome:    %8%\n"
                         "       Synthetic: %9%\n\n"
                         "-------Alignments within sampling regions (before subsampling)\n\n"
                         "       Genome:    %10%\n"
                         "       Synthetic: %11%\n\n"
                         "-------Alignments within sampling regions (after subsampling)\n\n"
                         "       Genome:    %12%\n"
                         "       Synthetic: %13%\n\n"
                         "       Normalization: %14% \u00B1 %15%\n\n"
                         "-------Before subsampling (within sampling regions)\n\n"
                         "       Synthetic coverage (average): %16%\n"
                         "       Genome coverage (average):    %17%\n\n"
                         "-------After subsampling (within sampling regions)\n\n"
                         "       Synthetic coverage (average): %18%\n"
                         "       Genome coverage (average):    %19%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()    // 1
                                            % gen         // 2
                                            % seq         // 3
                                            % stats.count // 4
                                            % meth2Str()  // 5
                                            % "????"      // 6
                                            % "????"      // 7
                                            % "????"      // 8
                                            % "????"      // 9
                                            % "????"      // 10
                                            % "????"      // 11
                                            % "????"      // 12
                                            % "????"      // 13
                                            % "????"      // 14
                                            % "????"      // 15
                                            % "????"      // 16
                                            % "????"      // 17
                                            % "????"      // 18
                                            % "????"      // 19
                     ).str());
    o.writer->close();
}

void VSample2::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = analyze(files[0], files[1], o);
    
    /*
     * Generating VarSubsample_summary.stats
     */
    
    generateSummary("VarSubsample_summary.stats", files[0], files[1], stats, o);

    /*
     * Generating VarSubsample_sequins.csv
     */

    generateCSV("VarSubsample_sequins.csv", stats, o);
}