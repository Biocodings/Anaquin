#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include "data/hist.hpp"
#include "data/bData.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "data/minters.hpp"
#include "data/dinters.hpp"
#include "tools/gtf_data.hpp"
//#include "RnaQuin/RnaQuin.hpp"
//#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct SequinVariant
    {
        enum class Context
        {
            Common,
            VeryLowGC,
            LowGC,
            HighGC,
            VeryHighGC,
            ShortDinRep,  // Dinucleotide repeats
            LongDinRep,   // Dinucleotide repeats
            ShortHompo,
            LongHompo,
            ShortQuadRep, // Quad-nucleotide repeats
            LongQuadRep,  // Quad-nucleotide repeats
            ShortTrinRep, // Trinucleotide repeats
            LongTrinRep,  // Trinucleotide repeats
            Cancer,
        } ctx;
        
        Genotype gt;
        
        // Copy number
        unsigned copy = 1;
    };

    struct VCFLadder
    {
        VCFData data;
        
        Ladder af;
        std::set<SequinID> vIDs;
        std::map<VarKey, SequinVariant> sVars;
    };

    enum class Tool
    {
        Test,
        Help,
        
        RnaAlign,
        RnaAssembly,
        RnaExpress,
        RnaFoldChange,
        RnaSubsample,
        RnaReport,
        
        VarCopy,
        VarAlign,
        VarDetect,
        VarSomatic,
        VarSample,
        VarTrim,
        VarFlip,
        VarKmer,
        VarStructure,
        VarConjoint,
        
        MetaCoverage,
        MetaAssembly,
        MetaSubsample,
    };
    
    /*
     * Different rules how two positions can be compared
     */

//    enum MatchRule
//    {
//        Exact,
//        Overlap,
//        Contains,
//    };
    
    class  Ladder;
    class  GTFData;
    struct VCFLadder;

    struct UserReference
    {
        std::shared_ptr<Ladder> l1, l2, l3, l4, l5, l6;

        // Ladder for VCF references and allele frequency
        std::shared_ptr<VCFLadder> vcf;
        
        // Translation
        std::shared_ptr<Translate> t1, t2;

        // First bed regions (not trimmed)
        std::shared_ptr<BedData> r1;

        // Second bed regions (trimmed)
        std::shared_ptr<BedData> r2;
        
        // GTF annoation
        std::shared_ptr<GTFData> g1;
    };
    
    class Reference
    {
        public:

            typedef std::set<SequinID> SequinIDs;
        
            inline SequinIDs seqs()   const { return _seqs; }
            inline SequinIDs seqsL1() const { return _l1->seqs; }
            inline SequinIDs seqsL2() const { return _l2->seqs; }
            inline SequinIDs seqsL3() const { return _l3->seqs; }

            inline Name t1(const Name &x) const { return _t1->translate(x); }
            inline Name t2(const Name &x) const { return _t2->translate(x); }

            inline Concent input1(const SequinID &x, Mixture m = Mix_1) const { return _l1->input(x, m); }
            inline Concent input2(const SequinID &x, Mixture m = Mix_1) const { return _l2->input(x, m); }
            inline Concent input3(const SequinID &x, Mixture m = Mix_1) const { return _l3->input(x, m); }
            inline Concent input4(const SequinID &x, Mixture m = Mix_1) const { return _l4->input(x, m); }
            inline Concent input5(const SequinID &x, Mixture m = Mix_1) const { return _l5->input(x, m); }
            inline Concent input6(const SequinID &x, Mixture m = Mix_1) const { return _l6->input(x, m); }

            // Allele frequency from VCF reference (not from mixture file)
            inline Proportion af(const SequinID &x) const { return _vcf->af.input(x, Mix_1); }

            // Position in the reference annoation
            inline Locus locus(const SequinID &id) const
            {
                for (const auto &i : *(_r1))
                {
                    if (i.second.r2d.count(id))
                    {
                        return i.second.r2d.at(id).l;
                    }
                }
                
                throw std::runtime_error("Region not found for " + id);
            }

            inline std::shared_ptr<GTFData> gtf() const { return _g1; }

            inline Chr2DInters regs1()  const { return _r1->inters();  }
            inline Chr2DInters regs2()  const { return _r2->inters();  }
            inline Chr2MInters mRegs1() const { return _r1->minters(); }
            inline Chr2MInters mRegs2() const { return _r2->minters(); }
        
            inline Counts nRegs() const { return _r1->count();  }
            inline Counts lRegs() const { return _r1->length(); }
        
            inline MergedIntervals<> mInters(const ChrID &cID) const { return _r1->minters(cID); }
            inline Chr2MInters mInters() const { return _r1->minters(); }

            inline void finalize(Tool x, const UserReference &r)
            {
                validate(x, r);
            }

        protected:

            inline void build(std::shared_ptr<BedData> r1)
            {
                _r1 = r1;
            }

            inline void build(std::shared_ptr<BedData> r1,
                              std::shared_ptr<BedData> r2,
                              std::shared_ptr<VCFLadder> vcf)
            {
                _r1  = r1;
                _r2  = r2;
                _vcf = vcf;
            }
        
            inline void build(std::shared_ptr<BedData> r1, std::shared_ptr<BedData> r2)
            {
                _r1 = r1;
                _r2 = r2;
            }

            inline void build(std::shared_ptr<Ladder> l1)
            {
                _l1   = l1;
                _seqs = l1->seqs;
            }

            inline void build(std::shared_ptr<Ladder> l1, std::shared_ptr<Ladder> l2)
            {
                _l1 = l1;
                _l2 = l2;
            }

            inline void build(std::shared_ptr<Ladder>    l1,
                              std::shared_ptr<Ladder>    l2,
                              std::shared_ptr<Translate> t1,
                              std::shared_ptr<Translate> t2)
            {
                _l1 = l1;
                _l2 = l2;
                _t1 = t1;
                _t2 = t2;
            }
        
            inline void build(std::shared_ptr<Ladder> l1,
                              std::shared_ptr<Ladder> l2,
                              std::shared_ptr<Ladder> l3,
                              std::shared_ptr<Ladder> l4,
                              std::shared_ptr<Ladder> l5,
                              std::shared_ptr<Ladder> l6)
            {
                _l1 = l1;
                _l2 = l2;
                _l3 = l3;
                _l4 = l4;
                _l5 = l5;
                _l6 = l6;
            }

            inline void build(std::shared_ptr<GTFData> g1)
            {
                _g1 = g1;
            }
        
            inline void build(std::shared_ptr<Ladder>  l1,
                              std::shared_ptr<Ladder>  l2,
                              std::shared_ptr<Ladder>  l3,
                              std::shared_ptr<Ladder>  l4,
                              std::shared_ptr<Ladder>  l5,
                              std::shared_ptr<Ladder>  l6,
                              std::shared_ptr<GTFData> g1)
            {
                _g1 = g1;
                _l1 = l1;
                _l2 = l2;
                _l3 = l3;
                _l4 = l4;
                _l5 = l5;
                _l6 = l6;
            }
        
            inline void build(std::shared_ptr<Ladder> l1,
                              std::shared_ptr<Ladder> l2,
                              std::shared_ptr<Ladder> l3,
                              std::shared_ptr<Ladder> l4,
                              std::shared_ptr<Ladder> l5)
            {
                _l1 = l1;
                _l2 = l2;
                _l3 = l3;
                _l4 = l4;
                _l5 = l5;
            }
        
            inline void build(std::shared_ptr<Ladder> l1, std::shared_ptr<BedData> r1)
            {
                _l1   = l1;
                _r1   = r1;
                _seqs = l1->seqs;
            }
        
            inline void build(std::shared_ptr<Ladder>  l1,
                              std::shared_ptr<BedData> r1,
                              std::shared_ptr<BedData> r2)
            {
                _l1   = l1;
                _r1   = r1;
                _r2   = r2;
                _seqs = l1->seqs;
            }

            virtual void validate(Tool, const UserReference &) = 0;

            template <typename T> std::vector<SequinID> merge(const std::set<T> &t1, const std::set<T> &t2)
            {
                std::set<SequinID> x, y;
                
                for (const auto &i : t1) { x.insert(static_cast<SequinID>(i)); }
                for (const auto &i : t2) { y.insert(static_cast<SequinID>(i)); }

                assert(!x.empty() && !y.empty());

                std::vector<SequinID> diffs, inters;
            
                /*
                 * Check for any sequin defined in x but not in y
                 */
            
                std::set_difference(x.begin(),
                                    x.end(),
                                    y.begin(),
                                    y.end(),
                                    std::back_inserter(diffs));

                /*
                 * Check for any sequin defined in both sets
                 */
            
                std::set_intersection(x.begin(),
                                      x.end(),
                                      y.begin(),
                                      y.end(),
                                      std::back_inserter(inters));

                return diffs;
            }
        
            template <typename T> std::vector<SequinID> merge(const std::set<T> &x)
            {
                return merge(x, x);
            }

            // Sequins
            std::set<SequinID> _seqs;
        
            // Sequin regions
            std::shared_ptr<BedData> _r1, _r2;

            // VCF references and allele frequency ladder
            std::shared_ptr<VCFLadder> _vcf;
        
            // Sequin regions
            std::shared_ptr<GTFData> _g1;
        
            // Sequin ladders
            std::shared_ptr<Ladder> _l1, _l2, _l3, _l4, _l5, _l6;

            // Translation
            std::shared_ptr<Translate> _t1, _t2;
    };

    /*
     * -------------------- Metagenomic Reference --------------------
     */
    
    class MetaRef : public Reference
    {
        public:
            MetaRef();

        protected:
        
            void validate(Tool, const UserReference &) override;
        
        private:

            struct MetaRefImpl;
            std::shared_ptr<MetaRefImpl> _impl;
    };

    /*
     * -------------------- Variant Reference --------------------
     */
    
    class VarRef : public Reference
    {
        public:

            VarRef();

            Counts nCNV(int)  const;
            Counts nGeno(Genotype)  const;
            Counts nType(Variation) const;
            Counts nContext(SequinVariant::Context) const;

            // Returns all reference variants
            std::set<Variant> vars() const;

            const SequinVariant &findSeqVar(long) const;
        
            const Variant *findVar(const ChrID &, const Locus &) const;

        protected:

            void validate(Tool, const UserReference &) override;

        private:

            struct VarRefImpl;

            std::shared_ptr<VarRefImpl> _impl;
    };
    
    /*
     * -------------------- Transcriptome Referenceb --------------------
     */
    
    struct GeneData;
    struct TransData;
    
    class RnaRef : public Reference
    {
        public:

            RnaRef();

            void readRef(const Reader &);

        protected:
        
            void validate(Tool, const UserReference &) override;
        
        private:

            struct RnaRefImpl;

            std::shared_ptr<RnaRefImpl> _impl;        
    };
}

#endif
