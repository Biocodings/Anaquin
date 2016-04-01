#ifndef KALLISTO_HPP
#define KALLISTO_HPP

#include <string>
#include <cstring>
#include <unistd.h>

// Defined by Kallisto
extern int __main__(int argc, char *argv[]);

namespace Anaquin
{
    struct Kallisto
    {
        // Where the output should be saved to
        static constexpr auto output = "/tmp/output";

        // Where the abundance file generated
        static constexpr auto abundFile = "/tmp/output/abundance.tsv";
        
        static void quant(const std::string &index, const std::string &file1, const std::string &file2)
        {            
            char *argv[8];
            
            /*
             * Eg: kallisto quant -i /data/index/VarQuin/AVA010.v032.index -o output LVA086.1_val_1.fq LVA086.2_val_2.fq
             */
            
            const auto output = std::string(Kallisto::output);
            
            argv[0] = new char[9];
            argv[1] = new char[6];
            argv[2] = new char[3];
            argv[3] = new char[index.size()+1];
            argv[4] = new char[3];
            argv[5] = new char[output.size()+1];
            argv[6] = new char[file1.size()+1];
            argv[7] = new char[file2.size()+1];
            
            strcpy(argv[0], "kallisto");
            strcpy(argv[1], "quant");
            strcpy(argv[2], "-i");
            strcpy(argv[3], index.c_str());
            strcpy(argv[4], "-o");
            strcpy(argv[5], output.c_str());
            strcpy(argv[6], file1.c_str());
            strcpy(argv[7], file2.c_str());
            
            optind = 1;
            
            /*
             * Execute Kallisto as if it were a standalone program
             */
            
            __main__(8, argv);
        }
        
    };
}

#endif
