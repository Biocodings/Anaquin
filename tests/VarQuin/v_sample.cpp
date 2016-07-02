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
    
    REQUIRE(r.n_gen == 421855);
    REQUIRE(r.n_syn == 32313);
    REQUIRE(r.syn.size() == 36);
    REQUIRE(r.gen.size() == 36);
    REQUIRE(r.synC == Approx(36));
    REQUIRE(r.genC == Approx(36));
}