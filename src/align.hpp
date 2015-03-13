#ifndef ALIGN_ANALYZER_HPP
#define ALIGN_ANALYZER_HPP

#include <string>

struct AlignAnalyzer
{
    /*
     * Analyze sequence alignments for sequins for a given SAM file.
     */
    
    static void analyze(const std::string &file);
};

#endif