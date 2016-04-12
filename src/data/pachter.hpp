#ifndef PACHTER_HPP
#define PACHTER_HPP

#include <cstring>
#include <unistd.h>
#include "data/types.hpp"
#include "data/script.hpp"
#include "stats/stats.hpp"

// Defined by Kallisto
extern int __main__(int argc, char *argv[]);

// Defined by resources.cpp
extern Anaquin::Scripts SleuthR();

namespace Anaquin
{
    struct Pachter
    {
        static FileName sleuth(const std::vector<FileName> &files,
                               const std::vector<std::string> &names,
                               const std::vector<std::string> &facts)
        {
            assert(!files.empty() && files.size() == facts.size());
            
            // Where the differential analysis should be saved
            const auto output = tmpFile();

            // The R-script to run
            const auto script = (boost::format(SleuthR()) % concat(files)
                                                          % concat(names)
                                                          % concat(facts)
                                                          % output).str();
            Script::run(script, "R CMD BATCH", "");
            
            return output;
        }
        
        static FileName externalQuant(const FileName &index,
                                      const FileName &file1,
                                      const FileName &file2,
                                      bool  bootstrap = false)
        {
            const auto output = tmpFile();
            
            const auto cmd = (boost::format("kallisto quant -i %1% -o %2% -b %3% %4% %5%")
                                                    % index
                                                    % output
                                                    % (bootstrap ? 200 : 10)
                                                    % file1
                                                    % file2).str();
            const int status = system(cmd.c_str());
            
            if (status)
            {
                throw std::runtime_error("Failed to run kallisto. Status code: " + std::to_string(status));
            }
            
            return output + "/abundance.tsv";
        }
        
        static FileName internalQuant(const FileName &index,
                                      const FileName &file1,
                                      const FileName &file2,
                                      bool  bootstrap = false)
        {
            char *argv[10];
            
            /*
             * Eg: kallisto quant -i /data/index/VarQuin/AVA010.v032.index -o output LVA086.1_val_1.fq LVA086.2_val_2.fq
             */
            
            const auto output = tmpFile();
            const auto abund  = output + "/abundance.tsv";
            
            argv[0] = new char[9];
            argv[1] = new char[6];
            argv[2] = new char[3];
            argv[3] = new char[index.size()+1];
            argv[4] = new char[3];
            argv[5] = new char[output.size()+1];
            argv[6] = new char[3];
            argv[7] = new char[5];
            argv[8] = new char[file1.size()+1];
            argv[9] = new char[file2.size()+1];
            
            strcpy(argv[0], "kallisto");
            strcpy(argv[1], "quant");
            strcpy(argv[2], "-i");
            strcpy(argv[3], index.c_str());
            strcpy(argv[4], "-o");
            strcpy(argv[5], output.c_str());
            strcpy(argv[6], "-b");
            strcpy(argv[7], (bootstrap ? "200" : "10"));
            strcpy(argv[8], file1.c_str());
            strcpy(argv[9], file2.c_str());
            
            optind = 1;
            
            // Execute Kallisto as if it were a standalone program
            __main__(10, argv);
            
            return abund;
        }
    };
}

#endif