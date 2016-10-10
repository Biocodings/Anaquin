#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_flip.hpp"

using namespace Anaquin;

TEST_CASE("VFlip_1")
{
    Test::clear();
    
    struct Impl : public VFlip::Impl
    {
        bool isReverse(const ChrID &)
        {
            return true;
        }
        
        void paired(const ParserSAM::Data &x, const ParserSAM::Data& y)
        {
            REQUIRE(x.name == y.name);
            pairs1[x.name] = x;
            pairs2[x.name] = y;
        }
        
        void nonPaired(const ParserSAM::Data &x)
        {
            unpairs[x.name] = x;
        }
        
        void unknownPaired(const ParserSAM::Data &x)
        {
            unknowns[x.name] = x;
        }

        std::map<ReadName, ParserSAM::Data> unpairs;
        std::map<ReadName, ParserSAM::Data> unknowns;
        std::map<ReadName, ParserSAM::Data> pairs1;
        std::map<ReadName, ParserSAM::Data> pairs2;
    };

    Impl impl;
    
    const auto r = VFlip::analyze("tests/data/genome.bam", VFlip::Options(), impl);

    REQUIRE(impl.pairs1.size()   == 112);
    REQUIRE(impl.pairs2.size()   == 112);
    REQUIRE(impl.unpairs.size()  == 0);
    REQUIRE(impl.unknowns.size() == 42);

    REQUIRE(impl.unknowns.count("1-hg38.fwd.NA12878_hets.sim_reads11906977"));
    REQUIRE(impl.unknowns["1-hg38.fwd.NA12878_hets.sim_reads11906977"].flag == 163);
    REQUIRE(impl.unknowns["1-hg38.fwd.NA12878_hets.sim_reads11906977"].cID  == "chr1");
    REQUIRE(impl.unknowns["1-hg38.fwd.NA12878_hets.sim_reads11906977"].seq  == "TCCTACACTACCGGTATCGATGGTATCGATAAAACGATAACATCTTACAGGTTCATCAACTCCCTGATGGTATATCTCGGACTTCTACTACGGTTGTCCTCCATTCCGTATAGTCTCTCGGTCTT");

    REQUIRE(impl.pairs1.count("1-hg38.fwd.NA12878_homos.sim_reads12533664"));
    REQUIRE(impl.pairs1["1-hg38.fwd.NA12878_homos.sim_reads12533664"].flag == 99);
    REQUIRE(impl.pairs1["1-hg38.fwd.NA12878_homos.sim_reads12533664"].cID  == "chr1");
    REQUIRE(impl.pairs1["1-hg38.fwd.NA12878_homos.sim_reads12533664"].seq  == "TTAACCTCCTATTTGTGGTGAAGACAGATCAAATATCAAAGAGTTCTTTAATCTATGACTGTCGGTCTACGCCACCGAGTGTGGACATTAGGGTCGTGAAACCCTCCGGTCCCTCCCGCCTAGTG");
    
    REQUIRE(impl.pairs2.count("1-hg38.fwd.NA12878_homos.sim_reads12533664"));
    REQUIRE(impl.pairs2["1-hg38.fwd.NA12878_homos.sim_reads12533664"].flag == 147);
    REQUIRE(impl.pairs2["1-hg38.fwd.NA12878_homos.sim_reads12533664"].cID  == "chr1");
    REQUIRE(impl.pairs2["1-hg38.fwd.NA12878_homos.sim_reads12533664"].seq  == "TACGTTACCGCACTAGAGCCGAGTGACGTTGGAGACGGAGGACCCAAGTTCGCTAAGAACACGGAGTCGGAGGACTCATCGACCCTGATGTCAACGGGTGGTGATGCGGGTCGATTAAAAACATA");
}
