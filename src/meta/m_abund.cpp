#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats MAbundance::report(const std::string &file, const MAbundance::Options &o)
{
    const auto stats = MAssembly::analyze(file, o);
    
    o.info("Generating linaer model");

    AnalyzeReporter::scatter(stats.lm, "MetaAbundance", "Expected abudnance (attomol/ul)", "K-Mer average", o.writer);

    return stats;
}