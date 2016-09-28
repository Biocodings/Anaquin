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

    std::cout << output << std::endl;
    
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
            stats.tsvs[Mix_1].push_back(samp.path + "/abundance.tsv");
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
            stats.tsvs[Mix_2].push_back(samp.path + "/abundance.tsv");
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
     * Perform sleuth differential analysis (recommended by the Kallisto team)
     */
    
    // Do we have at least two replicates for each mixture?
    if (stats.exp.samps.size() == 2 && stats.exp.samps[Mix_1].size() >= 2 && stats.exp.samps[Mix_2].size() >= 2)
    {
        const auto script = "library(sleuth)\n\n"
                            "# Where the Kallisto files are\n"
                            "path <- '%1%'\n\n"
                            "# Samples (eg. A1, A2...)\n"
                            "samps <- dir(file.path(path))\n\n"
                            "samps <- samps[samps!='sleuth.R']\n\n"
                            "# Factors for the samples\n"
                            "conds <- c(rep(1,sum(startsWith(samps, 'A'))), rep(0,sum(startsWith(samps, 'B'))))\n\n"
                            "# Construct full path for the samples\n"
                            "path <- paste(path, samps, sep='/')\n\n"
                            "s2c <- data.frame(sample=samps, condition=conds)\n"
                            "s2c <- dplyr::mutate(s2c, path=path)\n\n"
                            "so <- sleuth_prep(s2c, ~condition)\n"
                            "so <- sleuth_fit(so)\n"
                            "so <- sleuth_wt(so, 'condition1')\n\n"
                            "results <- sleuth_results(so, 'condition1')\n"
                            "write.csv(results, file='sleuth.csv', row.names=FALSE, quote=FALSE)";

        FileWriter fw(output);
        fw.open("sleuth.R");
        fw.write(((boost::format(script)) % output).str());
        fw.close();

        // Run differential analysis and save results to sleuth.csv
        System::runCmd("Rscript " + output + "/sleuth.R");
    }

    /*
     * Kallisto output files are in the output directory. We should analyze those files.
     */
    
    {
        RExpress::Options ro;
        
        ro.writer = o.writer;
        ro.format = RExpress::Format::Kallisto;
        
        for (auto &mix : stats.exp.samps)
        {
            {
                ro.metrs = RExpress::Metrics::Isoform;
                stats.iExpress[mix.first] = RExpress::analyze(stats.tsvs.at(mix.first), ro);
            }

            {
                ro.metrs = RExpress::Metrics::Gene;
                stats.gExpress[mix.first] = RExpress::analyze(stats.tsvs.at(mix.first), ro);
            }
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
        
        ro.work   = tmp;
        ro.logger = o.logger;
        ro.format = RExpress::Format::Kallisto;
        
        // Gene expression analysis for a mixture
        auto geneExpress = [&](const Title &title,
                               const std::vector<FileName> &files,
                               const std::vector<RExpress::Stats> &stats)
        {
            ro.metrs = RExpress::Metrics::Gene;
            
            const auto x = RExpress::generateSummary(files, stats, ro, "genes");
            const auto y = RExpress::generateCSV(stats, ro);
            const auto z = RExpress::generateRLinear("/RnaKReportGeneExpress.csv", stats, ro);

            FileWriter fw(tmp);
            fw.open("RnaKReportGeneExpress.csv");
            fw.write(y);
            fw.close();
            
            mark.start(title);
            mark.addText("Summary Statistics", x);
            mark.addText("Sequin Statistics",  y);
            mark.addRCode("Plot for Gene Expression", z);
            mark.end();
        };
        
        // Isoform expression analysis for a mixture
        auto isoExpress = [&](const Title &title,
                              const std::vector<FileName> &files,
                              const std::vector<RExpress::Stats> &stats)
        {
            ro.metrs = RExpress::Metrics::Isoform;
            
            const auto x = RExpress::generateSummary(files, stats, ro, "genes");
            const auto y = RExpress::generateCSV(stats, ro);
            const auto z = RExpress::generateRLinear("/RnaKReportIsoformExpress.csv", stats, ro);
            
            FileWriter fw(tmp);
            fw.open("RnaKReportIsoformExpress.csv");
            fw.write(y);
            fw.close();
            
            mark.start(title);
            mark.addText("Summary Statistics", x);
            mark.addText("Sequin Statistics",  y);
            mark.addRCode("Plot for Isoform Expression", z);
            mark.end();
        };
        
        for (auto &mix : stats.exp.samps)
        {
            const auto mixStr = mix.first == Mix_1 ? "Mixture A" : "Mixture B";
            
            {
                const auto format = "Gene Expression (%1%)";
                geneExpress(((boost::format(format) % mixStr).str()), stats.tsvs.at(mix.first), stats.gExpress.at(mix.first));
            }
            
            {
                const auto format = "Isoform Expression (%1%)";
                isoExpress(((boost::format(format) % mixStr).str()), stats.tsvs.at(mix.first), stats.iExpress.at(mix.first));
            }
        }
    }

    std::cout << tmp << std::endl;
    
    /*
     * C++ doesn't have the functionality to create a PDF report. Generate an Rmarkdown document and use it to create
     * a PDF document (other document types are also possible).
     */

    std::cout << tmp + "/report.Rmd" << std::endl;

    {
        FileWriter fw(tmp);
        fw.open("report.Rmd");
        fw.write(mark.generate("RnaQuin Report"));
        fw.close();
    }

    /*
     * Convert the markup document to PDF
     */

    {
        FileWriter fw(tmp);
        fw.open("r2pdf.R");
        fw.write("library(Anaquin)\nlibrary(rmarkdown)\nrender('report.Rmd', 'pdf_document')\n");
        fw.close();

        System::runCmd("cd " + tmp + "; Rscript " + tmp + "/r2pdf.R");
        System::runCmd("mv " + tmp + "/report.pdf " + o.work + "/RnaKReport_report.pdf");
    }
    
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
