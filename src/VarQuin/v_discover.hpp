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

                inline Counts fnTot() const { return m.fn(); }
                inline Counts fnSNP() const { return m_snp.fn(); }
                inline Counts fnInd() const { return m_ind.fn(); }
                
                inline Counts sTot() const { return tpTot() + fpTot(); }
                inline Counts sSNP() const { return tpSNP() + fpSNP(); }
                inline Counts sInd() const { return tpInd() + fpInd(); }

                inline Counts dTot() const { return tpTot() + fpTot() /*+ tnTot() */+ fnTot(); }
                inline Counts dSNP() const { return tpSNP() + fpSNP() /*+ tnSNP() */+ fnSNP(); }
                inline Counts dInd() const { return tpInd() + fpInd() /*+ tnInd() */+ fnInd(); }
                
                std::vector<VariantMatch> fps, tps;

                //std::map<long, VariantMatch *> fns_, tps_;
                std::map<long, VariantMatch> fns_, tps_;
                
                // Performance metrics
                Confusion m, m_snp, m_ind;
            };

            inline Counts countTP(const ChrID& cID) const { return data.at(cID).tps.size(); }
            inline Counts countFP(const ChrID& cID) const { return data.at(cID).fps.size(); }
            inline Counts countFN(const ChrID& cID) const { return data.at(cID).m.fn(); }

            inline Counts countSNP_TP(const ChrID& cID) const { return data.at(cID).tpSNP(); }
            inline Counts countInd_TP(const ChrID& cID) const { return data.at(cID).tpInd(); }
            inline Counts countVar_TP(const ChrID& cID) const { return data.at(cID).tpTot(); }
            inline Counts countSNP_FP(const ChrID& cID) const { return data.at(cID).fpSNP(); }
            inline Counts countInd_FP(const ChrID& cID) const { return data.at(cID).fpInd(); }
            inline Counts countVar_FP(const ChrID& cID) const { return data.at(cID).fpTot(); }
            inline Counts countSNP_FN(const ChrID& cID) const { return data.at(cID).fnSNP(); }
            inline Counts countInd_FN(const ChrID& cID) const { return data.at(cID).fnInd(); }
            inline Counts countVar_FN(const ChrID& cID) const { return data.at(cID).fnTot(); }

            inline Counts countSNP_TP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countSNP_TP(cID) : 0;
                });
            }

            inline Counts countSNP_TP_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countSNP_TP(cID) : 0;
                });
            }
            
            inline Counts countInd_TP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countInd_TP(cID) : 0;
                });
            }
            
            inline Counts countInd_TP_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countInd_TP(cID) : 0;
                });
            }

            inline Counts countVar_TP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countVar_TP(cID) : 0;
                });
            }
            
            inline Counts countVar_TP_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countVar_TP(cID) : 0;
                });
            }
            
            inline Counts countSNP_FP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countSNP_FP(cID) : 0;
                });
            }
            
            inline Counts countSNP_FP_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countSNP_FP(cID) : 0;
                });
            }
            
            inline Counts countInd_FP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countInd_FP(cID) : 0;
                });
            }
            
            inline Counts countInd_FP_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countInd_FP(cID) : 0;
                });
            }
            
            inline Counts countVar_FP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countVar_FP(cID) : 0;
                });
            }
            
            inline Counts countVar_FP_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countVar_FP(cID) : 0;
                });
            }
            
            inline Counts countSNP_FN_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countSNP_FN(cID) : 0;
                });
            }
            
            inline Counts countSNP_FN_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countSNP_FN(cID) : 0;
                });
            }
            
            inline Counts countInd_FN_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countInd_FN(cID) : 0;
                });
            }
            
            inline Counts countInd_FN_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countInd_FN(cID) : 0;
                });
            }
            
            inline Counts countVar_FN_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return Standard::isSynthetic(cID) ? countVar_FN(cID) : 0;
                });
            }
            
            inline Counts countVar_FN_Gen() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return !Standard::isSynthetic(cID) ? countVar_FN(cID) : 0;
                });
            }

            inline Proportion countVarSN_Syn() const
            {
                return (Proportion)countVar_TP_Syn() / (countVar_TP_Syn() + countVar_FN_Syn());
            }
            
            inline Proportion countVarSN_Gen() const
            {
                return (Proportion)countVar_TP_Gen() / (countVar_TP_Gen() + countVar_FN_Gen());
            }

            inline Proportion countSNPSN_Syn() const
            {
                return (Proportion)countSNP_TP_Syn() / (countSNP_TP_Syn() + countSNP_FN_Syn());
            }
            
            inline Proportion countSNPSN_Gen() const
            {
                return (Proportion)countSNP_TP_Gen() / (countSNP_TP_Gen() + countSNP_FN_Gen());
            }
            
            inline Proportion countIndSN_Syn() const
            {
                return (Proportion)countInd_TP_Syn() / (countInd_TP_Syn() + countInd_FN_Syn());
            }
            
            inline Proportion countIndSN_Gen() const
            {
                return (Proportion)countInd_TP_Gen() / (countInd_TP_Gen() + countInd_FN_Gen());
            }
            
            inline Proportion countVarPC_Syn() const
            {
                return (Proportion)countVar_TP_Syn() / (countVar_TP_Syn() + countVar_FP_Syn());
            }
            
            inline Proportion countVarPC_Gen() const
            {
                return (Proportion)countVar_TP_Gen() / (countVar_TP_Gen() + countVar_FP_Gen());
            }

            inline Proportion countSNPPC_Syn() const
            {
                return (Proportion)countSNP_TP_Syn() / (countSNP_TP_Syn() + countSNP_FP_Syn());
            }
            
            inline Proportion countSNPPC_Gen() const
            {
                return (Proportion)countSNP_TP_Gen() / (countSNP_TP_Gen() + countSNP_FP_Gen());
            }
            
            inline Proportion countIndPC_Syn() const
            {
                return (Proportion)countInd_TP_Syn() / (countInd_TP_Syn() + countInd_FP_Syn());
            }
            
            inline Proportion countIndPC_Gen() const
            {
                return (Proportion)countInd_TP_Gen() / (countInd_TP_Gen() + countInd_FP_Gen());
            }

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
