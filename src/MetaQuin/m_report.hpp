//#ifndef M_REPORT_HPP
//#define M_REPORT_HPP
//
//#include "stats/analyzer.hpp"
//#include "RnaQuin/r_fold.hpp"
//#include "RnaQuin/r_express.hpp"
//#include "parsers/parser_exp.hpp"
//
//namespace Anaquin
//{
//    struct MReport : public Analyzer
//    {
//        struct Options : public AnalyzerOptions
//        {
//            FileName index;
//        };
//
//        struct Stats
//        {
//            // Kallsito quantified files in TSV format
//            std::map<Mixture, std::vector<FileName>> tsvs;
//            
//            // Statistics for expression at the isoform level
//            std::map<Mixture, std::vector<RExpress::Stats>> iExpress;
//            
//            // Statistics for expression at the gene level
//            std::map<Mixture, std::vector<RExpress::Stats>> gExpress;
//            
//            // Statistics for differential at the isoform level
//            RFold::Stats iFold;
//            
//            // Statistics for differential at the gene level
//            RFold::Stats gFold;
//            
//            ParserExp::Experiment exp;
//        };
//        
//        static Stats analyze(const FileName &, const Options &);
//        static void  report (const FileName &, const Options &o = Options());
//    };
//}
//
//#endif
