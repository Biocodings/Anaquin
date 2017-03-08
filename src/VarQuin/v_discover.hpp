#ifndef V_DISCOVER_HPP
#define V_DISCOVER_HPP

#include <set>
#include <vector>
#include "stats/analyzer.hpp"
#include "tools/vcf_data.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VDiscover
    {
        struct Options : public AnalyzerOptions
        {
            Options() {}

            // Matching alleles?
            bool matchAllele = true;
            
            VarFormat format;
        };

        struct Stats : public MappingStats, public VariantStats
        {
            struct Data
            {
                // Measured minor allele frequency
                Proportion af;
                
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
                    return countSNP_TP(cID);
                });
            }

            inline Counts countInd_TP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countInd_TP(cID);
                });
            }
            
            inline Counts countVar_TP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countVar_TP(cID);
                });
            }
            
            inline Counts countSNP_FP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countSNP_FP(cID);
                });
            }
            
            inline Counts countInd_FP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countInd_FP(cID);
                });
            }
            
            inline Counts countVar_FP_Syn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countVar_FP(cID);
                });
            }
            
            inline Counts countSNP_FnSyn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countSNP_FN(cID);
                });
            }
            
            inline Counts countInd_FnSyn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countInd_FN(cID);
                });
            }
            
            inline Counts countVar_FnSyn() const
            {
                return ::Anaquin::count(data, [&](const ChrID &cID, const Data &x)
                {
                    return countVar_FN(cID);
                });
            }

            inline Proportion varF1() const
            {
                const auto sn = countVarSnSyn();
                const auto pc = countVarPC_Syn();
                
                return 2.0 * ((pc * sn) / (pc + sn));
            }

            inline Proportion countVarSnSyn() const
            {
                return (Proportion)countVar_TP_Syn() / (countVar_TP_Syn() + countVar_FnSyn());
            }
            
            inline Proportion countSNPSnSyn() const
            {
                return (Proportion)countSNP_TP_Syn() / (countSNP_TP_Syn() + countSNP_FnSyn());
            }
            
            inline Proportion indelF1() const
            {
                const auto sn = countIndSnSyn();
                const auto pc = countIndPC_Syn();
                
                return 2.0 * ((pc * sn) / (pc + sn));
            }
            
            inline Proportion countIndSnSyn() const
            {
                return (Proportion)countInd_TP_Syn() / (countInd_TP_Syn() + countInd_FnSyn());
            }
            
            inline Proportion countVarPC_Syn() const
            {
                return (Proportion)countVar_TP_Syn() / (countVar_TP_Syn() + countVar_FP_Syn());
            }
            
            inline Proportion countSNPPC_Syn() const
            {
                return (Proportion)countSNP_TP_Syn() / (countSNP_TP_Syn() + countSNP_FP_Syn());
            }

            inline Proportion SNPF1() const
            {
                const auto sn = countSNPSnSyn();
                const auto pc = countSNPPC_Syn();
                
                return 2.0 * ((pc * sn) / (pc + sn));
            }
            
            inline Proportion countIndPC_Syn() const
            {
                return (Proportion)countInd_TP_Syn() / (countInd_TP_Syn() + countInd_FP_Syn());
            }
            
            VCFData vData;
            
            // Distribution for the variants
            std::map<ChrID, HashHist> hist;

            std::map<ChrID, Data> data;

            /*
             * Structure for the query data. This is important because we might not be given
             * an annoation for the genome but we still might need to calculate number of variants
             * with allele frequency below and above the LOQ.
             */
            
            struct QueryData
            {
                // Measured allele frequency
                std::set<Proportion> af;
            };

            struct VarStats : public SequinStats, public LimitStats {};
            
            std::map<ChrID, QueryData> query;

            /*
             * Statistics for allele frequency
             */
            
            // Statistics for all variants
            VarStats vars;
            
            // Statistics for SNPs
            VarStats snp;
            
            // Statistics for indels
            VarStats ind;
            
            std::map<long, Counts> readR;
            std::map<long, Counts> readV;
            std::map<long, Counts> depth;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
