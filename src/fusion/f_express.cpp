#include "fusion/f_express.hpp"
#include "fusion/f_analyzer.hpp"

using namespace Anaquin;

ModelStats FExpress::analyze(const std::string &file, const Options &options)
{
    const auto stats = FAnalyzer::analyze(file, options);

    {
        AnalyzeReporter::linear(stats, "FusionExpress", "FPKM", options.writer);
    }

    {
        AnalyzeReporter::missing("FusionExpress_miss.csv", stats, options.writer);
    }

    return ModelStats(stats);
}