#ifndef V_DISCOVER_HPP
#define V_DISCOVER_HPP

#include <vector>
#include "stats/analyzer.hpp"
#include "tools/vcf_data.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VDiscover
    {
        typedef VarInput Input;
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Input input;
            
            // Significance level
            Probability sign = 0.1;
        };

        struct Stats : public MappingStats, public VariantStats
        {
            struct Data
            {
                inline Counts count(const std::vector<VariantMatch> &data, Mutation type) const
                {
                    return std::count_if(data.begin(), data.end(), [&](const VariantMatch &d)
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
                
                std::vector<VariantMatch> fps, tps, tns, fns;

                std::map<long, VariantMatch *> fns_, tps_;
                
                // Performance metrics
                Confusion m, m_snp, m_ind;
            };
            
            VCFData vData;
            
            // Distribution for the variants
            std::map<ChrID, HashHist> hist;

            std::map<ChrID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
