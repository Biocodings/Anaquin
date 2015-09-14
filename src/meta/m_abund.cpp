#include "meta/m_abund.hpp"

using namespace Anaquin;

MAbundance::Stats report(const std::string &file, const MAbundance::Options &o)
{
    const auto stats = MAssembly::analyze(file, o);
    
    o.info("Generating linaer model");
    AnalyzeReporter::linear(stats.lm, "MetaAssembly", "k-mer average", o.writer);

    return stats;
}