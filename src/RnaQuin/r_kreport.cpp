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
    
    for (const auto &info : stats.exp.samps[Mix_1])
    {
        KalSample samp;
        
        samp.p1   = info.p1;
        samp.p2   = info.p2;
        samp.path = output + "/A" + std::to_string(++i);
        
        samps.push_back(samp);
        stats.tsvs[Mix_1].push_back(samp.path + "/abundance.tsv");
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
                            "path <- list.dirs(path, recursive=FALSE)\n\n"
                            "samps <- basename(path)\n\n"
                            "conds <- c(rep('0',sum(startsWith(samps, 'A'))), rep('1',sum(startsWith(samps, 'B'))))\n\n"
                            "s2c <- data.frame(sample=samps, condition=conds)\n"
                            "s2c <- dplyr::mutate(s2c, path=path)\n\n"
                            "so <- sleuth_prep(s2c, ~condition)\n"
                            "so <- sleuth_fit(so)\n"
                            "so <- sleuth_wt(so, 'condition1')\n\n"
                            "results <- sleuth_results(so, 'condition1')\n"
                            "write.csv(results, file='%1%/sleuth.csv', row.names=FALSE, quote=FALSE)";
        
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
            const auto z = RExpress::generateRLinear("RnaKReportGeneExpress.csv", stats, ro);
            
            FileWriter::create(tmp, "RnaKReportGeneExpress.csv", y);
            
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
            const auto z = RExpress::generateRLinear("RnaKReportIsoformExpress.csv", stats, ro);

            // Required for R
            FileWriter::create(tmp, "RnaKReportIsoformExpress.csv", y);
            
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
                // Gene expression for each sample
                for (auto i = 0; i < mix.second.size(); i++)
                {
                    const auto format = "Gene Expression (%1%) for %2%";
                    geneExpress(((boost::format(format) % mixStr % mix.second[i].p1).str()),
                                  std::vector<FileName> { stats.tsvs.at(mix.first)[i] },
                                  std::vector<RExpress::Stats> { stats.gExpress.at(mix.first)[i] });
                }
            }
            
            {
                // Isoform expression for each sample
                for (auto i = 0; i < mix.second.size(); i++)
                {
                    const auto format = "Isoform Expression (%1%) for %2%";
                    isoExpress(((boost::format(format) % mixStr % mix.second[i].p1).str()),
                                 std::vector<FileName> { stats.tsvs.at(mix.first)[i] },
                                 std::vector<RExpress::Stats> { stats.iExpress.at(mix.first)[i] });
                }
            }

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
    
    /*
     * Differential analysis at the isoform and gene level
     */
    
    {
        RFold::Options ro;
        
        ro.work   = tmp;
        ro.logger = o.logger;
        ro.format = RFold::Format::Sleuth;
        
        const auto src = tmp + "/sleuth.csv";
        
        auto geneFold = [&](const Title &title,
                            const RFold::Stats &stats)
        {
            ro.metrs = RFold::Metrics::Gene;
            
            const auto x = RFold::generateSummary(src, stats, ro, "genes");
            const auto y = RFold::generateCSV(stats, ro);
            const auto z = RFold::generateRFold(stats, ro);
            
            // Required for R
            FileWriter::create(tmp, "RnaFoldChange_sequins.csv", y);

            mark.start(title);
            mark.addText("Summary Statistics", x);
            mark.addText("Sequin Statistics",  y);
            mark.addRCode("Plot for Fold Change", z);            
            mark.end();
        };
        
        geneFold("Differential Gene", stats.gFold);
    }
    
    /*
     * C++ doesn't have the functionality to create a PDF report. Generate an Rmarkdown document and use it to create
     * a PDF document (other document types are also possible).
     */
    
    std::cout << tmp + "/report.Rmd" << std::endl;
    
    FileWriter::create(tmp, "report.Rmd", mark.generate("RnaQuin Report"));
    FileWriter::create(tmp, "r2pdf.R", "library(Anaquin)\nlibrary(rmarkdown)\nrender('report.Rmd', 'pdf_document')\n");
    
    System::runCmd("cd " + tmp + "; Rscript " + tmp + "/r2pdf.R");
    System::runCmd("mv " + tmp + "/report.pdf " + o.work + "/RnaKReport_report.pdf");
}
