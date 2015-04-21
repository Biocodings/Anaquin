#include <iostream>
#include <assert.h>
#include "dalign.hpp"
#include "biology.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

DAlignStats DAlign::analyze(const std::string &file, const DAlign::Options &options)
{
    DAlignStats stats;
    const auto &r = Standard::instance();

    //Feature f1, f2;

    /*
     * At the alignment, we're only interested in counting for the gene-level.
     * Isoform-level is another possibiltiy but there wouldn't be information
     * to distinguish ambiguous reads from alternative splicing.
     */
    
    auto c = countsForGenes();
    
    ParserSAM::parse(file, [&](const Alignment &align)
    {
        classify(stats, align, [&](const Alignment &)
        {
            return true;
        });
 
        stats.n++;        
    });
//
//    assert(stats.nr + stats.nq == stats.n);
//    
//    const auto cr = Expression::analyze(c);
//
//    // Either the samples are independent or the least detectable-abundant sequin is known
////    assert(!cr.counts || r.r_seqs_gA.count(cr.limit_key));
//
///*
//    stats.sens.id = cr.limit_key;
//    stats.sens.counts = cr.limit_count;
//    stats.sens.exp = cr.limit_count ? r.r_seqs_gA.at(cr.limit_key).r.raw +
//                                      r.r_seqs_gA.at(cr.limit_key).v.raw: NAN;
//*/
//    
//    const std::string format = "%1%\t%2%\t%3%\t%4%";
//    
//    options.writer->open("align.stats");
//    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
//    options.writer->write((boost::format(format) % stats.dilution()
//                                                 % stats.m.sp()
//                                                 % stats.m.sn()
//                                                 % stats.s.exp).str());
//    options.writer->close();

	return stats;
}