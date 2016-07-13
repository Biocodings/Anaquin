#include <catch.hpp>
#include "unit/test.hpp"
#include "VarQuin/v_sample.hpp"

using namespace Anaquin;

extern Scripts AVA019_Bed();

TEST_CASE("VSubsample_ZeroProp")
{
    Test::clear();

    Standard::instance().addVStd(Reader(AVA019_Bed(),  DataMode::String));
    Standard::instance().r_var.finalize();
    
    VSample::Options o;
    o.p = 0.0;
    o.meth = VSample::Method::Prop;
    
    const auto r = VSample::stats("tests/data/sampled.bam");
    
    REQUIRE(r.tot.syn == 32323);
    REQUIRE(r.tot.gen == 422673);
    REQUIRE(r.syn.size() == 1);
    REQUIRE(r.gen.size() == 1);
    REQUIRE(r.synC == Approx(78.3321393416));
    REQUIRE(r.genC == Approx(0.0112222222));
}