#include "tools/script.hpp"
#include "VarQuin/v_allele.hpp"
#include "VarQuin/v_kreport.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

// Shared with other modules
extern Path __output__;

// Defined in resources.cpp
extern Scripts KReportScript();

static Scripts reportScript()
{
    return Script::trim(KReportScript());
}

VKReport::Stats VKReport::analyze(const FileName &data, const Options &o)
{
    VKReport::Stats stats;

    // Where the analysis should be saved
    const std::string output = tmpFile();
    
    std::cout << "QVarQuin " + o.index + " " + output + data << std::endl;
    
    Script::run(reportScript(), "python", "QVarQuin " + o.index + " " + output + " " + data);

    stats.output = output;

    {
        VAllele::Options ro;
        
        ro.writer = o.writer;
        ro.format = VAllele::Format::Salmon;

        // Statistics for allele frequency detection
        stats.allele = VAllele::analyze(output + "/quant.sf", ro);
    }

    return stats;
}

template <typename T> void setWorkDir(const Path &path, T &o)
{
    // Make sure to reuse the temporary directory
    o.writer = std::shared_ptr<FileWriter>(new FileWriter(path));
    
    __output__ = path;
}

void VKReport::report(const FileName &file, const Options &o)
{
    const auto stats = VKReport::analyze(file, o);
    
    o.info("Temporary: " + stats.output);
    
    /*
     * Allele frequency analysis
     */
    
    {
        VAllele::Options ro;
        
        ro.logger = o.logger;
        ro.format = VAllele::Format::Salmon;
        
        setWorkDir(stats.output + "/Allele", ro);
        
        VAllele::generateCSV("VarAllele_sequins.csv", stats.allele, ro);
        VAllele::generateSummary("VarAllele_summary.stats", stats.allele, ro);
        VAllele::generateR("VarAllele_linear.R", "VarAllele_sequins.csv", stats.allele, ro);
    }
    
    /*
     * Create a PDF report based on the generated files
     */

    Script::run(reportScript(), "python", "RVarQuin " + o.work + "/VarKReport_report.pdf " + stats.output);
}