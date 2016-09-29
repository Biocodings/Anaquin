#include <thread>
#include "tools/script.hpp"
#include "tools/markdown.hpp"
#include "VarQuin/v_allele.hpp"
#include "VarQuin/v_kreport.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

VKReport::Stats VKReport::analyze(const FileName &data, const Options &o)
{
    if (!System::checkConsole("salmon"))
    {
        throw MissingDependencyException("Salmon is not installed. Please consult the user guide on www.sequin.xyz and try again.");
    }
    else if (!System::checkConsole("R"))
    {
        throw MissingDependencyException("R is not installed. Please consult the user guide on www.sequin.xyz and try again.");
    }
    
    // Where the analysis files should be saved
    const auto output = System::tmpFile();
    
    std::cout << output << std::endl;
    
    // Create the directory structure
    System::runCmd("mkdir -p " + output);
    
    VKReport::Stats stats;
    
    // Parse the metadata
    stats.exp = ParserExp::parse(data);
    
    if (stats.exp.samps.empty())
    {
        throw InvalidInputError("Empty metadata: " + data);
    }
    
    std::vector<Sample> samps;
    
    /*
     * Generate Salmon command for the sample
     */
    
    auto i = 0;
    
    for (const auto &info : stats.exp.samps[Mix_1])
    {
        Sample samp;
        
        samp.p1   = info.p1;
        samp.p2   = info.p2;
        samp.path = output + "/A" + std::to_string(++i);
        
        samps.push_back(samp);
        //stats.tsvs[Mix_1].push_back(samp.path + "/abundance.tsv");
    }

    auto runSalmon = [&](const Sample &samp)
    {
        const auto format = "salmon quant -i %1% -l A -1 %2% -2 %3% -o %4%/salmon";
        System::runCmd((boost::format(format) % o.index % samp.p1 % samp.p2 % output).str());
    };
    
    // Multi-threaded instances
    std::vector<std::thread> kals;
    
    std::for_each(samps.begin(), samps.end(), [&](const Sample &samp)
    {
        kals.push_back(std::thread(runSalmon, samp));
    });

    std::for_each(kals.begin(), kals.end(), [&](std::thread &tID)
    {
        tID.join();
    });
    
    /*
     * Analyze Salmon output files
     */
    
    {
        VAllele::Options ro;
        
        ro.writer = o.writer;
        ro.format = VAllele::Format::Salmon;
        
        stats.allele = VAllele::analyze(output + "/salmon/quant.sf", ro);
    }

    return stats;
}

void VKReport::report(const FileName &file, const Options &o)
{
    const auto stats = VKReport::analyze(file, o);
    
    // Directory where the temporary files should be saved
    const auto tmp = System::tmpFile();
    
    MarkDown mark;
    
    /*
     * Allele frequency analysis
     */
    
    {
        VAllele::Options ro;
        
        ro.work   = tmp;
        ro.logger = o.logger;
        ro.format = VAllele::Format::Salmon;
        
        auto alleleFreq = [&](const Title &title,
                              const VAllele::Stats &stats)
        {
            const auto x = VAllele::generateSummary(stats, ro);
            const auto y = VAllele::generateCSV(stats, ro);
            const auto z = VAllele::generateRLinear("/VarKReportAllelle.csv", stats, ro);
            
            FileWriter::create(tmp, "VarKReportAllelle.csv", y);
            
            mark.start(title);
            mark.addText("Summary Statistics", x);
            mark.addText("Sequin Statistics",  y);
            mark.addRCode("Plot for Allele Frequency", z);
            mark.end();
        };
        
        alleleFreq("Allele Frequency", stats.allele);
    }
    
    /*
     * C++ doesn't have the functionality to create a PDF report. Generate an Rmarkdown document and use it to create
     * a PDF document (other document types are also possible).
     */
    
    std::cout << tmp + "/report.Rmd" << std::endl;
    
    FileWriter::create(tmp, "report.Rmd", mark.generate("VarQuin Report"));
    FileWriter::create(tmp, "r2pdf.R", "library(Anaquin)\nlibrary(rmarkdown)\nrender('report.Rmd', 'pdf_document')\n");
    
    System::runCmd("cd " + tmp + "; Rscript " + tmp + "/r2pdf.R");
    System::runCmd("mv " + tmp + "/report.pdf " + o.work + "/VarKReport_report.pdf");
}
