#include <catch.hpp>
#include "parsers/parser_bambed.hpp"

using namespace Anaquin;

typedef ParserBAMBED::Response Response;

TEST_CASE("ParserBAMBED_1")
{
    DIntervals<> chr10;
    chr10.add(DInter("Everything", Locus(10020, 133777009)));
    
    ID2Intervals c2i;
    c2i["chr10"] = chr10;

    const auto r = ParserBAMBED::parse("tests/data/sampled.bam", c2i,
                    [&](const ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *inter)
    {
        return Response::OK;
    }).inters.stats();
    
    /*
     * samtools view -F 4 sampled.bam | cut -f1,2,3 | grep chr10 | grep -v alt | cut -f1 | sort > A.txt
     */
    
    REQUIRE(r.aligns == 18480);
}
