#include <thread>
#include <boost/format.hpp>
#include "tools/script.hpp"
#include "tools/markdown.hpp"
#include "RnaQuin/r_kreport.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

RKReport::Stats RKReport::analyze(const FileName &data, const Options &o)
{
    if (!System::checkConsole("kallisto"))
    {
        throw MissingDependencyException("Kallisto is not installed. Please consult the user guide on www.sequin.xyz and try again.");
    }
    else if (!System::checkConsole("R"))
    {
        throw MissingDependencyException("R is not installed. Please consult the user guide on www.sequin.xyz and try again.");
    }
    else if (!System::checkRPack("sleuth"))
    {
        throw MissingDependencyException("Sleuth is not installed. Please consult the user guide on www.sequin.xyz and try again.");
    }

    // Where the analysis files should be saved
    const auto output = System::tmpFile();

    // Create the directory structure
    System::runCmd("mkdir -p " + output);
    
    RKReport::Stats stats;

    // Parse the metadata
    stats.exp = ParserExp::parse(data);

    if (stats.exp.samps.empty())
    {
        throw InvalidInputError("Empty metadata: " + data);
    }
    
    struct KalSample
    {
        Path path;
        FileName p1, p2;
    };
    
    // Eg: First and second paired
    std::vector<KalSample> samps;

    /*
     * Generate Kallisto commands for each sample
     */
    
    // Eg: A1, A2 A3
    Counts i = 0;
    
    if (!stats.exp.samps.at(Mix_1).empty())
    {
        for (const auto &info : stats.exp.samps[Mix_1])
        {
            KalSample samp;

            samp.p1   = info.p1;
            samp.p2   = info.p2;
            samp.path = output + "/A" + std::to_string(++i);

            samps.push_back(samp);
            stats.abunds.push_back(samp.path + "/abundance.tsv");
        }
    }

    // Eg: B1, B2 B3
    i = 0;

    if (!stats.exp.samps.at(Mix_2).empty())
    {
        for (const auto &info : stats.exp.samps[Mix_2])
        {
            KalSample samp;
            
            samp.p1   = info.p1;
            samp.p2   = info.p2;
            samp.path = output + "/B" + std::to_string(++i);
            
            samps.push_back(samp);
            stats.abunds.push_back(samp.path + "/abundance.tsv");
        }
    }
    
    auto runKallisto = [&](const KalSample &samp)
    {
        const auto format = "kallisto quant -b 500 -i %1% -o %2% %3% %4%";
        System::runCmd((boost::format(format) % o.index % samp.path % samp.p1 % samp.p2).str());
    };
    
    // Multi-threaded instances
    std::vector<std::thread> kals;
    
    std::for_each(samps.begin(), samps.end(), [&](const KalSample &samp)
    {
        kals.push_back(std::thread(runKallisto, samp));
    });

    /*
     * Wait until Kallisto quantification completes
     */
    
    std::for_each(kals.begin(), kals.end(), [&](std::thread &tID)
    {
        tID.join();
    });
    
    /*
     * Kallisto output files are in the output directory. We should analyze those files.
     */
    
    {
        RExpress::Options ro;
        
        ro.writer = o.writer;
        ro.format = RExpress::Format::Kallisto;
        
        {
            ro.metrs = RExpress::Metrics::Isoform;
            stats.iExpress = RExpress::analyze(stats.abunds, ro);
        }

        {
            ro.metrs = RExpress::Metrics::Gene;
            stats.gExpress = RExpress::analyze(stats.abunds, ro);
        }
    }
    
//    {
//        RFold::Options ro;
//        
//        ro.writer = o.writer;
//        ro.format = RFold::Format::Sleuth;
//        
//        {
//            ro.metrs = RFold::Metrics::Isoform;
//            stats.iFold = RFold::analyze(output + "/sleuth.csv", ro);
//        }
//        
//        {
//            ro.metrs = RFold::Metrics::Gene;
//            stats.gFold = RFold::analyze(output + "/sleuth.csv", ro);
//        }
//    }

    return stats;
}

void RKReport::report(const FileName &file, const Options &o)
{
    const auto stats = RKReport::analyze(file, o);

    // Directory where the temporary files should be saved
    const auto tmp = System::tmpFile();
    
    MarkDown mark;
    
    /*
     * Expression analysis at the isoform and gene level
     */
    
    {
        RExpress::Options ro;
        
        ro.logger = o.logger;
        ro.format = RExpress::Format::Kallisto;
        
        {
            ro.metrs = RExpress::Metrics::Gene;
            
            const auto x = RExpress::generateSummary(stats.abunds, stats.gExpress, ro, "genes");
            const auto y = RExpress::generateCSV(stats.gExpress, ro);
            
            mark.start("Gene Analysis");
            mark.addText("Summary Statistics", x);
            mark.addText("Sequin Statistics", y);
            mark.end();
        }

        {
            ro.metrs = RExpress::Metrics::Isoform;

            const auto x = RExpress::generateSummary(stats.abunds, stats.iExpress, ro, "isoforms");
            const auto y = RExpress::generateCSV(stats.iExpress, ro);

            mark.start("Isoform Analysis");
            mark.addText("Summary Statistics", x);
            mark.addText("Sequin Statistics", y);
            mark.end();
        }
    }

    /*
     * C++ doesn't have the functionality to create a PDF report. Generate an Rmarkdown document and use it to create
     * a PDF document (other document types are also possible).
     */

    std::cout << tmp + "/report.Rmd" << std::endl;
    
    FileWriter fw("/Users/tedwong/Sources/QA");
    fw.open("report.Rmd");
    fw.write(mark.generate("RnaQuin Report"));
    fw.close();

    /*
     * Differential analysis at the isoform and gene level
     */
    
//    {
//        RFold::Options ro;
//        
//        ro.logger = o.logger;
//        ro.format = RFold::Format::Sleuth;
//        
//        const auto src = stats.output + "/sleuth.csv";
//        
//        {
//            setWorkDir(stats.output + "/FoldG", ro);
//            ro.metrs = RFold::Metrics::Gene;
//            
//            RFold::generateCSV("RnaFoldChange_sequins.csv", stats.gFold, ro);
//            RFold::generateSummary("RnaFoldChange_summary.stats", src, stats.gFold, ro, "genes");
//            RFold::generateR(stats.gFold, ro);
//        }
//        
//        {
//            setWorkDir(stats.output + "/FoldI", ro);
//            ro.metrs = RFold::Metrics::Isoform;
//            
//            RFold::generateCSV("RnaFoldChange_sequins.csv", stats.iFold, ro);
//            RFold::generateSummary("RnaFoldChange_summary.stats", src, stats.iFold, ro, "isoforms");
//            RFold::generateR(stats.iFold, ro);
//        }
//    }
}
