#include <iostream>
#include <assert.h>
#include "d_align.hpp"
#include "biology.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

DAlignStats DAlign::analyze(const std::string &file, const Options &options)
{
    DAlignStats stats;
    //const auto &r = Standard::instance();

    //Feature f1, f2;

    //auto c = countsForGenes();

//    ParserSAM::parse(file, [&](const Alignment &align)
//    {
        // There shouldn't be any splicing for DNA
//        assert(!align.spliced);
        
//        classify(stats, align, [&](const Alignment &)
//        {
//            /*
//            
//            if ((!align.spliced && find(r.fs.begin(), r.fs.end(), align, f1)) ||
//                (align.spliced && checkSplice(r, align, f1, f2)))
//            {
//                assert(r.iso2Gene.count(f1.iID));
//                
//                ce[r.iso2Gene.at(f1.iID)]++;
//                cb[r.iso2Gene.at(f1.iID)]++;
//                
//                return true;
//            }
//            else
//            {
//                return false;
//            }
//            */
//            
//            return true;
//        });
//    });
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

	return DAlignStats();
}