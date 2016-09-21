#include "tools/script.hpp"
#include "RnaQuin/r_kreport.hpp"
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

RKReport::Stats RKReport::analyze(const FileName &data, const Options &o)
{
    RKReport::Stats stats;

    stats.exp = ParserExp::parse(data);
    
    // Where the analysis should be saved
    const std::string output = tmpFile();
    
    std::cout << "QRnaQuin " + o.index + " " + output + data << std::endl;
    
    /*
     * Generate Kallisto and Sleuth analysis.
     *
     *   Eg: python kexpress.py RNAQuin ARN004.v032.index 1,0 /tmp/kallisto experiment.txt
     */
    
    Script::run(reportScript(), "python", "RnaQuin " + o.index + " " + output + " " + data);

    stats.output = output;

    stats.kFiles.push_back(output + "/A1/abundance.tsv");
    stats.kFiles.push_back(output + "/A2/abundance.tsv");
    stats.kFiles.push_back(output + "/A3/abundance.tsv");
    stats.kFiles.push_back(output + "/B1/abundance.tsv");
    stats.kFiles.push_back(output + "/B2/abundance.tsv");
    stats.kFiles.push_back(output + "/B3/abundance.tsv");
    
    {
        RExpress::Options ro;
        
        ro.writer = o.writer;
        ro.format = RExpress::Format::Kallisto;
        
        {
            ro.metrs = RExpress::Metrics::Isoform;
            stats.iExpress = RExpress::analyze(stats.kFiles, ro);
        }
        
        {
            ro.metrs = RExpress::Metrics::Gene;
            stats.gExpress = RExpress::analyze(stats.kFiles, ro);
        }
    }
    
    {
        RFold::Options ro;
        
        ro.writer = o.writer;
        ro.format = RFold::Format::Sleuth;
        
        {
            ro.metrs = RFold::Metrics::Isoform;
            stats.iFold = RFold::analyze(output + "/sleuth.csv", ro);
        }
        
        {
            ro.metrs = RFold::Metrics::Gene;
            stats.gFold = RFold::analyze(output + "/sleuth.csv", ro);
        }
    }

    return stats;
}

template <typename T> void setWorkDir(const Path &path, T &o)
{
    // Make sure to reuse the temporary directory
    o.writer = std::shared_ptr<FileWriter>(new FileWriter(path));
    
    __output__ = path;
}

void RKReport::report(const FileName &file, const Options &o)
{
    const auto stats = RKReport::analyze(file, o);
    
    o.info("Temporary: " + stats.output);
    
    /*
     * Expression analysis at the isoform and gene level
     */
    
    {
        RExpress::Options ro;
        
        ro.logger = o.logger;
        ro.format = RExpress::Format::Kallisto;
        
        {
            setWorkDir(stats.output + "/ExpressG", ro);
            ro.metrs = RExpress::Metrics::Gene;
            
            RExpress::generateCSV("RnaExpress_sequins.csv", stats.gExpress, ro);
            RExpress::generateSummary("RnaExpress_summary.stats", stats.kFiles, stats.gExpress, ro, "genes");
            RExpress::generateR("RnaExpress_linear.R", "RnaExpress_sequins.csv", stats.gExpress, ro);
        }
        
        {
            setWorkDir(stats.output + "/ExpressI", ro);
            ro.metrs = RExpress::Metrics::Isoform;

            RExpress::generateCSV("RnaExpress_sequins.csv", stats.iExpress, ro);
            RExpress::generateSummary("RnaExpress_summary.stats", stats.kFiles, stats.iExpress, ro, "genes");
            RExpress::generateR("RnaExpress_linear.R", "RnaExpress_sequins.csv", stats.iExpress, ro);
        }
    }
    
    /*
     * Differential analysis at the isoform and gene level
     */
    
    {
        RFold::Options ro;
        
        ro.logger = o.logger;
        ro.format = RFold::Format::Sleuth;
        
        const auto src = stats.output + "/sleuth.csv";
        
        {
            setWorkDir(stats.output + "/FoldG", ro);
            ro.metrs = RFold::Metrics::Gene;
            
            RFold::generateCSV("RnaFoldChange_sequins.csv", stats.gFold, ro);
            RFold::generateSummary("RnaFoldChange_summary.stats", src, stats.gFold, ro, "genes");
            RFold::generateR(stats.gFold, ro);
        }
        
        {
            setWorkDir(stats.output + "/FoldI", ro);
            ro.metrs = RFold::Metrics::Isoform;
            
            RFold::generateCSV("RnaFoldChange_sequins.csv", stats.iFold, ro);
            RFold::generateSummary("RnaFoldChange_summary.stats", src, stats.iFold, ro, "isoforms");
            RFold::generateR(stats.iFold, ro);
        }
    }
    
    /*
     * Create a PDF report based on the generated files
     */

    Script::run(reportScript(), "python", "RRnaQuin " + o.work + "/RnaKReport_report.pdf " + stats.output);
}