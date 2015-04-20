#include "assembly.hpp"
#include "classify.hpp"
#include "standard.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "parsers/parser_gtf.hpp"

using namespace Spike;

AssemblyStats Assembly::analyze(const std::string &file, const Assembly::Options &options)
{
    AssemblyStats stats;
    const auto &r = Standard::instance();

    /*
     * Create counters for various levels. The number of counts will give out the
     * distribution.
     */

    INIT_COUNTER(c);
    INIT_COUNTER(c_trans);
    INIT_COUNTER(c_exons);
    INIT_COUNTER(c_introns);
    
    std::vector<Feature> exons;
    
	ParserGTF::parse(file, [&](const Feature &f)
	{
        classify(stats, f,
                 [&](const Feature &) // Positive
                 {
                     switch (f.type)
                     {
                         case Transcript:
                         {
                             assert(r.seqs_iA.count(f.iID));
                             const auto &seq = r.seqs_iA.at(f.iID);
                             
                             assert(c.count(f.iID) && c_trans.count(f.iID));
                             
                             // True if the locus match
                             if (tfp(seq.l == f.l, &stats.m_trans))
                             {
                                 c[f.iID]++;
                                 c_trans[f.iID]++;
                                 return true;
                             }
                             else
                             {
                                 return false;
                             }
                             
                             break;
                         }

                         case Exon:
                         {
                             exons.push_back(f);
                             assert(c.count(f.iID) && c_exons.count(f.iID));

                             if (tfp(find(r.exons, f), &stats.m_exon))
                             {
                                 c[f.iID]++;
                                 c_exons[f.iID]++;
                                 return true;
                             }
                             else
                             {
                                 return false;
                             }
                             
                             break;
                         }

                         default:
                         {
                             throw std::runtime_error("Unknown assembly type!");
                         }
                     }
                 });
	});

    assert(!r.introns.empty());
    
    extractIntrons(exons, [&](const Feature &, const Feature &, Feature &f)
    {
        f.id = "chrT"; // TODO: Fix it later...
        
        classify(stats, f,
                 [&](const Feature &)
                 {
                     if (tfp(find(r.introns, f), &stats.m_intron))
                     {
                         c_introns[f.iID]++;
                         return true;
                     }
                     else
                     {
                         return false;
                     }
                 });
    });

    ANALYZE_COUNTS(c, base_r);
    ANALYZE_COUNTS(c_trans, trans_r);
    ANALYZE_COUNTS(c_exons, exon_r);
    ANALYZE_COUNTS(c_introns, intron_r);

    stats.s        = base_r.sens();
    stats.s_trans  = trans_r.sens();
    stats.s_exon   = exon_r.sens();
    stats.s_intron = intron_r.sens();
        
    assert(stats.n && stats.nr + stats.nq == stats.n);
    
    /*
     * Reports various statistics related to the "accuracy" of the transcripts in each sample when
     * compared to the reference silico-data. The typical gene finding measures of "sensitivity" and
     * "specificity" are calculated at various levels (nucleotide, exon, intron, transcript, gene) for
     * for the input file. The Sn and Sp columns show specificity and sensitivity values at each level.
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%";

    options.writer->open("base.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m.sp()
                                                 % (stats.m.sn() == NAN ? 99 : 1)
                                                 % stats.s.exp).str());
    options.writer->close();

    options.writer->open("exons.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m_exon.sp()
                                                 % stats.m_exon.sn()
                                                 % stats.s.exp).str());
    options.writer->close();

    options.writer->open("intron.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m_intron.sp()
                                                 % stats.m_intron.sn()
                                                 % stats.s.exp).str());
    options.writer->close();
    
    options.writer->open("transcripts.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m_trans.sp()
                                                 % stats.m_trans.sn()
                                                 % stats.s.exp).str());
    options.writer->close();

	return stats;
}