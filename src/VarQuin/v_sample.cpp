#include "VarQuin/v_sample.hpp"

using namespace Anaquin;

typedef Subsampler::Stats Stats;

Stats VSample::stats(const FileName &file, const Options &o)
{
    return Subsampler::stats(file, o);
}

void VSample::report(const FileName &file, const Options &o)
{
    Subsampler::report(file, o);
}
