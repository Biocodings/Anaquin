#include "variant/v_allele.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

static double alleleFreq(const BaseSeq &m)
{
    assert(m.sequins.size() == 2);
    
    const auto ref = m.sequins.begin()->first;
    const auto var = m.sequins.rbegin()->first;
    
    // Abundance for the reference
    const auto r = m.sequins.at(ref).abund();
    
    // Abundance for the variant
    const auto v = m.sequins.at(var).abund();

    // Abundance ratio of reference to variant DNA standard
    return v / (r + v);
}

VAllele::Stats VAllele::analyze(const std::string &file, const Options &o)
{
    VAllele::Stats stats;
    const auto &s = Standard::instance();

    long n = 0;

    o.info("Parsing VCF file");

    ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
    {
        n++;
        
        if (v.id != s.id)
        {
            return;
        }
        
        Variation match;
        
        Confusion m;

        if (classify(m, v, [&](const VCFVariant &)
        {
            // Can we find this variant?
            if (!s.v_vars.count(v.l))
            {
                return Negative;
            }

            match = s.v_vars.at(v.l);

            // Does the variant match with the meta?
            if (match.type != v.type || match.alt != v.alt || match.ref != v.ref)
            {
                return Negative;
            }

            assert(s.bases_1.count(match.id));
            
            const auto &base = s.bases_1.at(match.id);
            
            /*
             * Plotting the relative allele frequency that is established by differences
             * in the concentration of reference and variant DNA standards.
             */
            
            // The measured coverage is the number of base calls aligned and used in variant calling
            const auto measured = (double) v.dp_a / (v.dp_r + v.dp_a);

            // The known coverage for allele frequnece
            const auto known = alleleFreq(base);

            stats.add(match.id, known, measured);
  
            return Positive;
        }))
        {
            //stats.h.at(match)++;
        }
    });

    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Ignored: %2% variants not in chrT\n"
                         "   Found: %3% variants in chrT\n"
                         "   Reference: %4% variants\n\n"
                         "Fuzzy: %5%\n\n"
                         "Correlation:     %6%\n"
                         "Slope:     %7%\n"
                         "R2:     %8%\n"
    ;
    
    const auto lm = stats.linear();
    
    o.writer->open("VarAllele_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % (n - stats.size())
                                            % stats.size()
                                            % s.v_vars.size()
                                            % o.fuzzy
                                            % lm.r
                                            % lm.m
                                            % lm.r2).str());
    
    o.writer->close();

    AnalyzeReporter::scatter(stats, "VarAllele", "", o.writer);

    return stats;
}