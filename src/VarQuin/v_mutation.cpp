#include "tools/system.hpp"
#include "VarQuin/v_germ.hpp"
#include "VarQuin/v_somatic.hpp"
#include "writers/vcf_writer.hpp"
#include "VarQuin/v_mutation.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

// Defined main.cpp
extern FileName Bed1Ref();

template <typename T, typename O> void writeSamples(const FileName &file, const T &stats, const O &o)
{
    const auto r2 = Standard::instance().r_var.regs2();
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "Qual").str());
    
    for (const auto &i : stats.es.vs)
    {
        A_ASSERT(!i.name.empty());
        
        o.writer->write((boost::format(format) % i.name
                                               % i.cID
                                               % i.l.start
                                               % i.readR
                                               % i.readV
                                               % i.depth
                                               % toString(i.qual)).str());
    }
    
    o.writer->close();
}

void VMutation::report(const FileName &endo, const FileName &seqs, const VMutation::Options &o)
{
    o.info("Edge: " + toString(o.edge));

    BaseCallerStats stats;
    
    if (o.isGerm)
    {
        o.info("Germline analysis");
        
        VGerm::Options o_;
        
        o_.work = o.work;
        o_.logger = o.logger;
        o_.filter = o.filter;
        o_.writer = o.writer;
        o_.output = o.output;

        stats = VGerm::report(endo, seqs, o_);
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

        stats = VSomatic::report(endo, seqs, o_);
    }

    /*
     * Generating VarMutation_sample.vcf (done in VGerm or VSomatic)
     */
    
    /*
     * Generating VarMutation_sample.tsv
     */
    
    writeSamples("VarMutation_sample.tsv", stats, o);
    
    /*
     * Generating VarMuation_sequins.bed
     */
    
    System::copy(Bed1Ref(), o.work + "/VarMuation_sequins.bed");
}
