//#ifndef L_DIFFS_HPP
//#define L_DIFFS_HPP
//
//#include "stats/analyzer.hpp"
//
//namespace Anaquin
//{
//    struct LDiffs
//    {
//        typedef AnalyzerOptions Options;
//
//        struct Stats : public AlignmentStats
//        {
//            struct Data : public LinearStats
//            {
//                // Empty Implementation
//            };
//            
//            std::map<ChrID, Data> data;
//            
//            // Sensitivity at the joined level
//            Limit s_joined;
//
//            // Sensitivity at the unjoined level
//            Limit ss;
//
//            // Histogram at the unjoined level
//            SequinHist h = Standard::instance().r_lad.hist();
//            
//            // Histogram at the joined level
//            LadderRef::JoinHist h_joined  = Standard::instance().r_lad.joinHist();
//        };
//
//        static Stats report(const FileName &, const FileName &, const Options &o = Options());
//    };
//}
//
//#endif