#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_flip.hpp"

using namespace Anaquin;

TEST_CASE("VFlip_1")
{
    clrTest();
    
    struct Impl : public VFlip::Impl
    {
        typedef VFlip::Status Status;
        
        bool isReverse(const ChrID &)
        {
            return true;
        }
        
        void process(const ParserBAM::Data &x, const ParserBAM::Data &y, Status status)
        {
            switch (status)
            {
                case Status::ReverseReverse:
                {
                    REQUIRE(x.name == y.name);
                    rever1[x.name] = x;
                    rever2[x.name] = y;
                    break;
                }
                    
                case Status::ReverseNotMapped:
                {
                    break;
                }
                    
                case Status::ForwardForward:
                {
                    break;
                }
                    
                case Status::ForwardReverse:
                {
                    break;
                }
                    
                case Status::ForwardNotMapped:
                {
                    break;
                }
                    
                case Status::NotMappedNotMapped:
                {
                    break;
                }
                    
                case Status::RevHang:
                {
                    hangs[x.name] = x;
                    break;
                }
                    
                case Status::ForHang:
                {
                    hangs[x.name] = x;
                    break;
                }                    
            }
        }

        std::map<ReadName, ParserBAM::Data> hangs;
        std::map<ReadName, ParserBAM::Data> rever1;
        std::map<ReadName, ParserBAM::Data> rever2;
    };

    Impl impl;
    
    const auto r = VFlip::analyze("tests/data/genome.bam", VFlip::Options(), impl);

    REQUIRE(impl.rever1.size()  == 112);
    REQUIRE(impl.rever2.size()  == 112);
    REQUIRE(impl.hangs.size() == 42);

    REQUIRE(impl.hangs.count("1-hg38.fwd.NA12878_hets.sim_reads11906977"));
    REQUIRE(impl.hangs["1-hg38.fwd.NA12878_hets.sim_reads11906977"].flag == 163);
    REQUIRE(impl.hangs["1-hg38.fwd.NA12878_hets.sim_reads11906977"].cID  == "chr1");
    REQUIRE(impl.hangs["1-hg38.fwd.NA12878_hets.sim_reads11906977"].seq  == "TCCTACACTACCGGTATCGATGGTATCGATAAAACGATAACATCTTACAGGTTCATCAACTCCCTGATGGTATATCTCGGACTTCTACTACGGTTGTCCTCCATTCCGTATAGTCTCTCGGTCTT");

    REQUIRE(impl.rever1.count("1-hg38.fwd.NA12878_homos.sim_reads12533664"));
    REQUIRE(impl.rever1["1-hg38.fwd.NA12878_homos.sim_reads12533664"].flag == 99);
    REQUIRE(impl.rever1["1-hg38.fwd.NA12878_homos.sim_reads12533664"].cID  == "chr1");
    REQUIRE(impl.rever1["1-hg38.fwd.NA12878_homos.sim_reads12533664"].seq  == "TTAACCTCCTATTTGTGGTGAAGACAGATCAAATATCAAAGAGTTCTTTAATCTATGACTGTCGGTCTACGCCACCGAGTGTGGACATTAGGGTCGTGAAACCCTCCGGTCCCTCCCGCCTAGTG");
    
    REQUIRE(impl.rever2.count("1-hg38.fwd.NA12878_homos.sim_reads12533664"));
    REQUIRE(impl.rever2["1-hg38.fwd.NA12878_homos.sim_reads12533664"].flag == 147);
    REQUIRE(impl.rever2["1-hg38.fwd.NA12878_homos.sim_reads12533664"].cID  == "chr1");
    REQUIRE(impl.rever2["1-hg38.fwd.NA12878_homos.sim_reads12533664"].seq  == "TACGTTACCGCACTAGAGCCGAGTGACGTTGGAGACGGAGGACCCAAGTTCGCTAAGAACACGGAGTCGGAGGACTCATCGACCCTGATGTCAACGGGTGGTGATGCGGGTCGATTAAAAACATA");
}
