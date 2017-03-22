#include <catch.hpp>
#include "readers/reader_bam.hpp"

using namespace Anaquin;

TEST_CASE("ReaderBam_1")
{
    Intervals<> chr10;
    chr10.add(Interval("Everything", Locus(10020, 133777009)));
    
    ID2Intervals c2i;
    c2i["chr10"] = chr10;

    const auto r = ReaderBam::stats("tests/data/sampled.bam", c2i, [&](const ParserSAM::Data &x, const ParserSAM::Info &info, const Interval *inter)
    {
        return ReaderBam::Response::OK;
    }).inters.stats();
    
    /*
     * samtools view -F 4 sampled.bam | cut -f1,2,3 | grep chr10 | grep -v alt | cut -f1 | sort > A.txt
     */
    
    REQUIRE(r.aligns == 18480);
}
