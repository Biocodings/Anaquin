#include "tools/system.hpp"
#include "VarQuin/v_germ.hpp"
#include "VarQuin/v_somatic.hpp"
#include "writers/vcf_writer.hpp"
#include "VarQuin/v_mutation.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

// Defined main.cpp
extern FileName Bed1Ref();

static void writeRegions(const FileName &file, const VMutation::Options &o)
{
    extern FileName Bed2Ref();
    System::copy(Bed2Ref(), o.work + "/" + file);
}

void VMutation::report(const FileName &endo, const FileName &seqs, const VMutation::Options &o)
{
    o.info("Edge: " + toString(o.edge));

    if (o.isGerm)
    {
        o.info("Germline analysis");
        
        VGerm::Options o_;
        
        o_.work = o.work;
        o_.logger = o.logger;
        o_.filter = o.filter;
        o_.writer = o.writer;
        o_.output = o.output;

        VGerm::report(endo, seqs, o_);
    }
    else
    {
        o.info("Somatic analysis");

        VSomatic::Options o_;
        
        o_.work = o.work;
        o_.logger = o.logger;
        o_.filter = o.filter;
        o_.writer = o.writer;
        o_.output = o.output;

        VSomatic::report(endo, seqs, o_);
    }

    /*
     * Generating VarMuation_sequins.bed
     */
    
    writeRegions("VarMuation_sequins.bed", o);
}
