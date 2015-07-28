#include "fusion/f_express.hpp"
#include "fusion/f_analyzer.hpp"

using namespace Anaquin;

ModelStats FExpress::analyze(const std::string &file, const Options &options)
{
    const auto stats = FAnalyzer::analyze(file, options);

    /*
     * Generate summary statistics
     */

    {
        AnalyzeReporter::linear(stats, "FExpress", "FPKM", options.writer);
    }
 
    // Conver to the statistics expected
    return ModelStats(stats);
}