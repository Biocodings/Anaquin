#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats MAbundance::report(const std::string &file, const MAbundance::Options &o)
{
    const auto stats = MAssembly::analyze(file, o);
    
    o.info("Generating linaer model");

    o.writer->open("MetaAssembly_quins.stats");
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";
    
    o.writer->write((boost::format(format) % "ID"
                                           % "con"
                                           % "status"
                                           % "avg_align"
                                           % "avg_sequin"
                                           % "covered").str());
    
    for (const auto &meta : stats.blat.metas)
    {
        const auto &align = meta.second;
        const auto detect = align->contigs.size() != 0;
        const auto status = detect ? std::to_string(align->covered) : "-";
    
        o.writer->write((boost::format(format) % align->seq->id
                                               % align->seq->mixes.at(Mix_1)
                                               % status
                                               % (detect ? std::to_string(align->depthAlign)  : "-")
                                               % (detect ? std::to_string(align->depthSequin) : "-")
                                               % (detect ? std::to_string(align->covered)     : "-")).str());
    }

    AnalyzeReporter::scatter(stats.lm, "MetaAbundance", "Expected abudnance (attomol/ul)", "K-Mer average", o.writer);

    return stats;
}