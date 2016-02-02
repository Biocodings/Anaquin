#include "VARQuin/v_discover.hpp"

using namespace Anaquin;

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    VDiscover::Stats stats;
    
    stats.data[ChrT];
    
    parseVariant(file, o.caller, [&](const VariantMatch &m)
    {
        if (m.query->chrID == ChrT)
        {
            Stats::Classifed cls;

            cls.seq   = m.seq;
            cls.query = *(m.query);
            cls.eAllFreq = m.eAllFreq;
            
            if (m.match && m.ref && m.alt)
            {
                stats.data.at(ChrT).tps.push_back(cls);
            }
            else
            {
                /*
                 * By definition, this is a false-positive because there is no match by position.
                 */

                if (!m.seq)
                {
                    assert(false);
                }
                
                stats.data.at(ChrT).fps.push_back(cls);
            }
        }
    });
    
    return stats;
}

static void writeClass(const FileName &file,
                       const std::vector<VDiscover::Stats::Classifed> &data,
                       const VDiscover::Options &o)
{
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";

    o.writer->open(file);
    o.writer->write((boost::format(format) % "Sequin"
                                           % "Position"
                                           % "PValue"
                                           % "RefRead"
                                           % "VarRead"
                                           % "EAlleleF"
                                           % "Type").str());

    for (const auto &i : data)
    {
        std::string type;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { type = "SNP";   break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { type = "Indel"; break; }
        }

        o.writer->write((boost::format(format) % i.seq->id
                                               % i.query.l.start
                                               % i.query.pval
                                               % i.query.readR
                                               % i.query.readV
                                               % i.eAllFreq
                                               % type).str());
    }

    o.writer->close();
}

void VDiscover::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    const auto &data = stats.data.at(ChrT);
    
    o.logInfo("Number of false positives: " + std::to_string(data.fps.size()));
    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */
    
    /*
     * Generating true positives
     */

    writeClass("VarDiscover_TP.csv", stats.data.at(ChrT).tps, o);

    /*
     * Generating false positives
     */
    
    writeClass("VarDiscover_FP.csv", stats.data.at(ChrT).fps, o);
    
    /*
     * Generating ROC curve
     */
    
    o.writer->open("VarDiscover_ROC.R");
    o.writer->write(RWriter::createROC_V(o.working, "VarDiscover_TP.csv", "VarDiscover_FP.csv"));
    o.writer->close();

    /*
     * Generating LODR
     */

    o.writer->open("VarDiscover_LODR.R");
    o.writer->close();
    

//    const auto summary = "Summary for dataset: %1%\n\n"
//                         "   Experiment:  %2% variants\n"
//                         "   Synthetic:   %3% variants\n"
//                         "   Reference:   %4% variants\n"
//                         "   Detected:    %5% variants\n"
//                         "   False-Pos:   %6% variants\n\n"
//                         "   Sensitivity: %7%\n"
//                         "   Specificity: %8%";
//
//    o.writer->open("VarDiscover_summary.stats");
//    o.writer->write((boost::format(summary) % file
//                                            % stats.chrT->n_endo
//                                            % stats.chrT->n_chrT
//                                            % r.countVars()
//                                            % stats.chrT->m.tp()
//                                            % (stats.chrT->n_chrT - stats.chrT->m.tp())
//                                            % stats.chrT->m.sn()
//                                            % stats.chrT->m.pc()).str());
//    o.writer->close();
//
//    /*
//     * Generate statistics for each sequin
//     */
//
//    o.writer->open("VarDiscover_quins.stats");
//    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
//    
//    const auto format = "%1%\t%2%";
//    o.writer->write((boost::format(format) % "id" % "detected").str());
//
//    for (const auto &i : stats.chrT->h)
//    {
//        o.writer->write((boost::format(format) % i.first % stats.chrT->h.at(i.first)).str());
//    }
//    
//    o.writer->close();
}