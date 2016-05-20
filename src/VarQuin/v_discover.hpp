#ifndef V_DISCOVER_HPP
#define V_DISCOVER_HPP

#include <vector>
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VDiscover
    {
        enum class Software
        {
            GATK,
            VarScan,
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            Software soft;
            
            // Significance level
            Probability sign = 0.1;
        };

        struct Stats : public MappingStats, public VariantStats
        {
            typedef VariantMatch ChrTData;
            
            struct ChrTStats
            {
                inline Counts count(const std::vector<ChrTData> &data, Mutation type) const
                {
                    return std::count_if(data.begin(), data.end(), [&](const ChrTData &d)
                    {
                        return (d.query.type() == type) ? 1 : 0;
                    });
                }

                inline Counts tpTot() const { return tps.size(); }
                inline Counts tpSNP() const { return count(tps, Mutation::SNP); }
                inline Counts tpInd() const { return tpTot() - tpSNP(); }

                inline Counts fpTot() const { return fps.size(); }
                inline Counts fpSNP() const { return count(fps, Mutation::SNP); }
                inline Counts fpInd() const { return fpTot() - fpSNP(); }

                inline Counts tnTot() const { return tns.size(); }
                inline Counts tnSNP() const { return count(tns, Mutation::SNP); }
                inline Counts tnInd() const { return tnTot() - tnSNP(); }

                inline Counts fnTot() const { return fns.size(); }
                inline Counts fnSNP() const { return count(fns, Mutation::SNP); }
                inline Counts fnInd() const { return fnTot() - fnSNP(); }
                
                inline Counts sTot() const { return tpTot() + fpTot(); }
                inline Counts sSNP() const { return tpSNP() + fpSNP(); }
                inline Counts sInd() const { return tpInd() + fpInd(); }

                inline Counts dTot() const { return tpTot() + fpTot() + tnTot() + fnTot(); }
                inline Counts dSNP() const { return tpSNP() + fpSNP() + tnSNP() + fnSNP(); }
                inline Counts dInd() const { return tpInd() + fpInd() + tnInd() + fnInd(); }
                
                std::vector<ChrTData> fps, tps, tns, fns;

                // Performance metrics
                Confusion m, m_snp, m_ind;
            };

            typedef std::vector<CalledVariant> GenomeStats;
            
            // Statistics for synthetic variants
            ChrTStats chrT;

            // Statistics for genomic variants
            GenomeStats geno;
            
            // Distribution for variants
            HashHist hist;            
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
