#include <catch.hpp>
#include "writers/r_writer.hpp"

// Defined in resources.cpp
extern std::string PlotTROC();

using namespace Anaquin;

//TEST_CASE("LinearInflection_Summary")
//{
//    const auto r1 = StatsWriter::inflectSummary();
//    const auto r2 = StatsWriter::inflectSummary("chrT.gtf", "endo.gtf", SInflectStats(), "genes");
//
//    REQUIRE(!r1.empty());
//    REQUIRE(!r2.empty());
//}

//TEST_CASE("R_ROC_Plot")
//{
//    const auto seqs = std::vector<std::string> { "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31" };
//    
//    const auto pvals = std::vector<double> { 0.88150509946551,0.790579826554692,0.703769035448768,0.327173030430958,0.675541902369305,0.27463411000412,0.128385466840651,0.89663342888243,0.0113084185194,0.747229352440727,0.399276537215399,0.292443920793433,0.0112254429527438,0.954407036854045,0.455910328939659,0.0962627125313444,0.0368961848301403,0.00993122073795335,0.839274592316032,0.816494961078089,0.0808554362362422,0.870104248850572,0.0639398485934761,0.0503499345800991,0.634624821872723,0.00767564931491453,0.658344706952809,0.0569646760950342,0.857124116442357,0.0468269349419728,0.966861022197228 };
//    
//    const auto labels = std::vector<std::string> { "FP","FP","FP","FP","FP","TP","TP","FP","TP","FP","FP","TP","TP","FP","FP","TP","TP","TP","TP","FP","TP","FP","TP","TP","TP","TP","FP","TP","FP","TP","FP" };
//
//    REQUIRE_NOTHROW(PlotTROC(seqs, pvals, "gene"));
//}

TEST_CASE("R_LODR_Plot")
{
    const auto seqs = std::vector<std::string> { "R1_43","R1_52","R1_91","R2_1","R2_26","R2_32","R2_60","R1_102","R1_11","R1_22","R1_24","R1_32","R1_51","R1_72","R2_67","R2_68","R1_31","R1_93","R2_105","R2_143","R2_47","R2_6","R2_76","R1_101","R1_42","R1_61","R1_73","R2_150","R2_19","R2_27","R2_57","R1_63","R1_83","R2_115","R2_28","R2_38","R1_62","R2_37","R2_41","R2_55","R2_66","R2_72","R2_14","R1_12","R1_53","R2_42","R2_46","R2_73","R1_103","R1_14","R1_33","R1_41","R1_92","R2_153","R2_154","R1_13","R1_21","R1_23","R1_71","R1_81","R1_82","R2_116","R2_117","R2_140","R2_151","R2_152","R2_18","R2_20","R2_24","R2_45","R2_53","R2_54","R2_59","R2_63","R2_65","R2_7","R2_71" };
    
    const auto avgs = std::vector<double> { 2705.086868,7.49993785,1.592708121,1.230225393,30245.16372,5.815765526,57.31422832,7.861308236,49.19302776,50.41496679,651.9783904,7.985523155,133.7062386,1.915751089,1.073442671,29.01734957,1137.120055,119.4224592,0.18917613,1.097868629,885.0918089,56.41690991,1.682529281,3.924674515,4603.840008,1.063762344,1009.656343,799.0494006,2483.998333,4.952578082,4.375024791,15084.0071,72.92686167,371.1565224,2.953588894,1.68953155,1.01873706,1.777813952,212.6304813,21327.89169,14764.17376,2.492149405,19030.73113,8.859330357,14.00238258,2.767311968,1.93160761,7.175418831,691.7270976,197.5698542,2.636168967,2381.400168,88.00397943,1.751980522,1690.384465,5658.731489,20826.87313,13.07524022,8130.48873,66.37341747,3234.696422,2.420004433,49.36082846,70.83321355,0.715429424,43.70764867,5519.955865,15.24277683,19.99246051,7.024112983,5.195766319,2147.18809,1.256549113,966.0105647,1.362640676,858.8064128,1.881553695 };
    
    const auto pvals = std::vector<double> { 1e-300,0.001169934,0.139450731,0.073978738,9.99999999999997e-311,0.092481602,1.13e-20,0.002833409,2.4e-19,6.48e-21,1.04e-204,0.00057502,1.12e-50,0.69161262,0.45297249,2.21e-13,1.99e-246,2.62e-35,0.581316493,0.358522136,2.7e-247,3.87e-18,0.209233963,0.105514837,1e-305,0.441015205,2.38e-209,2.8e-147,1.09e-295,0.025384555,0.02590173,4.41e-226,7.42e-15,3.42e-83,0.978385972,0.667864534,0.930023203,0.580193256,5e-48,3.07e-122,6.09e-197,0.473146607,7.42e-49,0.307018313,0.475080251,0.922872647,0.714491618,0.475080251,2.29e-26,2.35e-09,0.967411086,3.24e-168,0.012707059,0.575509154,7.88e-31,0.69161262,0.427165219,0.954603427,0.954603427,0.94632463,0.010126637,0.682926366,0.564752562,0.226542104,0.572511389,0.667864534,0.241245248,0.922872647,0.486450354,0.466319074,0.842866512,0.166400526,0.358522136,3.74e-39,0.332977339,0.550500262,0.493058443 };
    
    const auto logFCs = std::vector<double> { 4,4,4,4,4,4,4,-4,-4,-4,-4,-4,-4,-4,-4,-4,3,3,3,3,3,3,3,-3,-3,-3,-3,-3,-3,-3,-3,2,2,2,2,2,-2,-2,-2,-2,-2,-2,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

    //REQUIRE_NOTHROW(RWriter::createLODR(seqs, avgs, pvals, logFCs));
}