//#ifndef L_COPY_HPP
//#define L_COPY_HPP
//
//#include "stats/analyzer.hpp"
//
//namespace Anaquin
//{
//    struct LCopy
//    {
//        typedef SingleMixtureOption Options;
//
//        struct Stats : public AlignmentStats
//        {
//            std::vector<std::string> ids;
//            std::vector<double> expected;
//            std::vector<double> measured;
//
//            // Detected duplicates
//            std::set<Locus> detected;
//        };
//
//        static Stats analyze(const FileName &, const Options &o = Options());
//        static void  report (const FileName &, const Options &o = Options());
//    };
//}
//
//#endif