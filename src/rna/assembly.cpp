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

    INIT_COUNTER(c_base);
    INIT_COUNTER(c_trans);
    INIT_COUNTER(c_exons);
    INIT_COUNTER(c_introns);
    
    std::vector<Feature> exons;
    
	ParserGTF::parse(file, [&](const Feature &f, ParserProgress &p)
	{
        classify(r, stats, f,
                 [&]()
                 {
                     switch (f.type)
                     {
                         case Transcript:
                         {
                             assert(r.seqs_iA.count(f.iID));
                             const auto &seq = r.seqs_iA.at(f.iID);

                             assert(c_base.count(f.iID));
                             assert(c_trans.count(f.iID));

                             // It's a true-positive if the locus match
                             if (tfp(seq.l == f.l, &stats.m_base, &stats.m_trans))
                             {
                                 c_base[f.iID]++;
                                 c_trans[f.iID]++;
                             }

                             break;
                         }

                         case Exon:
                         {
                             exons.push_back(f);
                             tfp(find(r.exons, f), &stats.m_base, &stats.m_exon);
                             break;
                         }

                         default:
                         {
                             throw std::runtime_error("Unknown assembly type!");
                         }
                     }
                 },

                 [&](bool mapped)
                 {
                     switch (f.type)
                     {
                         case Transcript:
                         {
                             tfn(mapped, &stats.m_base, &stats.m_trans);
                             break;
                         }

                         case Exon:
                         {
                             tfn(mapped, &stats.m_base, &stats.m_exon);
                             break;
                         }

                         default: { throw std::runtime_error("Unknown assembly type!"); };
                     }
                 });
	});

    assert(!r.introns.empty());

    
    extractIntrons(exons, [&](const Feature &, const Feature &, Feature &f)
    {
        f.id = "chrT"; // TODO: Fix it later...
        
        classify(r, stats, f,
                 [&]()
                 {
                     tfp(find(r.introns, f), &stats.m_base, &stats.m_intron);
                 },
                 [&](bool mapped)
                 {
                     tfn(mapped, &stats.m_base, &stats.m_intron);
                 });
    });
    
    ANALYZE_COUNTS(c_base,base_r);
    ANALYZE_COUNTS(c_trans, trans_r);
    ANALYZE_COUNTS(c_exons, exon_r);
    ANALYZE_COUNTS(c_introns, intron_r);

    stats.sens_base.id     = base_r.limit_key;
    stats.sens_base.counts = base_r.limit_count;
    stats.sens_base.exp    = base_r.limit_count ? r.seqs_iA.at(base_r.limit_key).exp +
                                                  r.seqs_iA.at(base_r.limit_key).exp: NAN;
    stats.sens_trans.id     = trans_r.limit_key;
    stats.sens_trans.counts = trans_r.limit_count;
    stats.sens_trans.exp    = trans_r.limit_count ? r.seqs_iA.at(trans_r.limit_key).exp +
                                                    r.seqs_iA.at(trans_r.limit_key).exp: NAN;
    stats.sens_exon.id      = exon_r.limit_key;
    stats.sens_exon.counts  = exon_r.limit_count;
    stats.sens_exon.exp     = exon_r.limit_count ? r.seqs_iA.at(exon_r.limit_key).exp +
                                                   r.seqs_iA.at(exon_r.limit_key).exp: NAN;
    stats.sens_intron.id      = intron_r.limit_key;
    stats.sens_intron.counts  = intron_r.limit_count;
    stats.sens_intron.exp     = intron_r.limit_count ? r.seqs_iA.at(intron_r.limit_key).exp +
                                                       r.seqs_iA.at(intron_r.limit_key).exp: NAN;
    
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
                                                 % stats.m_base.sp()
                           % (stats.m_base.sn() == NAN ? 99 : 1)
                                                 % stats.sens_base.exp).str());
    options.writer->close();

    options.writer->open("exons.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m_exon.sp()
                                                 % stats.m_exon.sn()
                                                 % stats.sens_exon.exp).str());
    options.writer->close();

    options.writer->open("intron.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m_intron.sp()
                                                 % stats.m_intron.sn()
                                                 % stats.sens_intron.exp).str());
    options.writer->close();
    
    options.writer->open("transcripts.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m_trans.sp()
                                                 % stats.m_trans.sn()
                                                 % stats.sens_trans.exp).str());
    options.writer->close();
    
	return stats;
}
