#include "assembly.hpp"
#include "classify.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Spike;

AssemblyStats Assembly::analyze(const std::string &file, const AssemblyOptions &options)
{
    AssemblyStats stats;
    const auto r = StandardFactory::reference();

	ParserGTF::parse(file, [&](const Feature &f, ParserProgress &p)
	{
        classify(r, f,
                 [&]() // Positive
                 {
                     switch (f.type)
                     {
                         case Exon:
                         {
                             verifytPositive(find(r.exons, f), &stats.base, &stats.exon);
                             break;
                         }
                             
                         case Transcript:
                         {
                             // TODO: Need information for the transcript...
                             break;
                         }

                         default:
                         {
                             throw std::runtime_error("Unknown assembly type!");
                         }
                     }
                 },
                 [&](bool mapped) // Negative
                 {
                     verifyNegative(mapped, &stats.base, &stats.exon);
                 });
	});

    /*
     * Reports various statistics related to the "accuracy" of the transcripts in each sample
     * when compared to the reference silico data. The typical gene finding measures of "sensitivity"
     * and "specificity" are calculated at various levels (nucleotide, exon, intron, transcript,
     * gene) for the input file. The Sn and Sp columns show specificity and sensitivity values at each
     * level.
     */
    
	return stats;
}
