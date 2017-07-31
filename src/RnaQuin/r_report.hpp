#ifndef R_REPORT_HPP
#define R_REPORT_HPP

#include "tools/system.hpp"
#include "stats/analyzer.hpp"

extern std::string PythonReport();

namespace Anaquin
{
    struct RReport
    {
        typedef AnalyzerOptions Options;
        
        static void report(const FileName &index, const FileName &rMix, const FileName &meta, const Options &o)
        {
            const auto format = "RnaReport %1% %2% %3% %4%/RnaReport.pdf";
            const auto cmd = ((boost::format(format) % index
                                                     % meta
                                                     % rMix
                                                     % o.work)).str();
            
            // Everything is done in the Python script
            System::runScript(PythonReport(), cmd);
        }
    };
}

#endif
