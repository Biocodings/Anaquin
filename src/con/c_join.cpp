#include "con/c_join.hpp"
#include "stats/expression.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_sam.hpp"

using namespace Spike;

CJoin::Stats CJoin::analyze(const std::string &file, const Options &options)
{
    CJoin::Stats stats;

    // Construct a histogram of the aligned sequins
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        // Don't repeat the same read if it's spliced
        if (align.i == 0)
        {
            stats.hist[align.id]++;
        }
    });

    const auto &s = Standard::instance();

    
    
    
    
    
    
    
    
//    // Library size for expected
//    const auto libSizeSequin = 31.4572;
//    
//    const auto libSizeActual = 38.4;
//
//    for (const auto &i : s.c_seqs_A)
//    {
//        const std::string base = "GA116"; // i.first;
//        
//        /*
//         * A -> 2x
//         * B -> 4x
//         * C -> 6x
//         * D -> 8x
//         */
//        
//        const auto base_A = base + "_A";
//        const auto base_B = base + "_B";
//        const auto base_C = base + "_C";
//        const auto base_D = base + "_D";
//        
//        const auto countA = (double)stats.hist.at(base_A) / libSizeActual;
//        const auto countB = (double)stats.hist.at(base_B) / libSizeActual;
//        const auto countC = (double)stats.hist.at(base_C) / libSizeActual;
//        const auto countD = (double)stats.hist.at(base_D) / libSizeActual;
//        
//        const auto sequin = s.c_seqs_A.at(base);
//        const auto foldC  = sequin.abund() / 10;
//        
//        const auto expA = 1.0 * foldC / libSizeSequin;
//        const auto expB = 2.0 * foldC / libSizeSequin;
//        const auto expC = 4.0 * foldC / libSizeSequin;
//        const auto expD = 8.0 * foldC / libSizeSequin;
//
//        const auto excessA = countA / expA;
//        const auto excessB = countB / expB;
//        const auto excessC = countC / expC;
//        const auto excessD = countD / expD;
//        
//        std::vector<double> x = { expA, expB, expC, expD };
//        std::vector<double> y = { countA, countB, countC, countD };
//        
//        const auto m = SS::lm("y ~ x", SS::data.frame(SS::c(y), SS::c(x)));
//
//        std::cout << m.coeffs[1].v << std::endl;
//        
//        const auto correctA = countA / m.coeffs[1].v;
//        const auto correctB = countB / m.coeffs[1].v;
//        const auto correctC = countC / m.coeffs[1].v;
//        const auto correctD = countD / m.coeffs[1].v;
//
//        
//        
//    }
//    
//    
    
    
    //c_seqs_A
    
    
    /*
     * Write out histogram
     */

    options.writer->open("conjoin_histogram.stats");
    
    for (const auto &i : stats.hist)
    {
        options.writer->write((boost::format("%1%\t%2%") % i.first % i.second).str());
    }

    options.writer->close();
    
	return stats;
}