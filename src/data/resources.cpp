#include <string>

/*
 * Scripts
 */

#include "resources/sleuth.R"
#include "resources/viewer.py"
#include "resources/manual.txt"
#include "resources/reports.py"

/*
 * Fusion Resources
 */

#include "resources/plotROC_F.R"
#include "resources/plotExpress_F.R"

#include "resources/AFU004.v032.ref"
#include "resources/AFU005.v032.bed"
#include "resources/MFU007.v013.csv"

/*
 * Ladder Resources
 */

#include "resources/plotLadderAbund.R"

#include "resources/MLA014.v013.csv"
#include "resources/MLA016.v013.csv"
#include "resources/MLA020.v013.csv"

/*
 * Transcriptome Resources
 */

#include "resources/plotMA.R"
#include "resources/plotFold.R"
#include "resources/plotTROC.R"
#include "resources/plotMajor.R"
#include "resources/plotLODR.R"
#include "resources/plotRAbundAbund.R"
#include "resources/plotTAbundAbund.R"

#include "resources/ATR001.v032.gtf"
#include "resources/MTR002.v013.csv"
#include "resources/MTR003.v013.csv"
#include "resources/MTR004.v013.csv"
#include "resources/MTR005.v013.csv"

/*
 * Metagenomics Resources
 */

#include "resources/plotMFold.R"

#include "resources/MME023.v013.csv"
#include "resources/AME015.v032.bed"

/*
 * Variant Resources
 */

#include "resources/plotROC_V.R"
#include "resources/plotLODR_V.R"
#include "resources/plotDensity.R"
#include "resources/plotSubsample.R"
#include "resources/plotVAbundAbund.R"
#include "resources/plotAlleleReads.R"
#include "resources/plotAlleleAllele.R"

#include "resources/AVA009.v032.vcf"
#include "resources/MVA011.v013.csv"
#include "resources/MVA012.v013.csv"
#include "resources/AVA017.v032.bed"

#define ToString(x) std::string(reinterpret_cast<char*>(x))

typedef std::string Scripts;

Scripts ReportScript()
{
    return ToString(scripts_reports_py);
}

Scripts ViewerScript()
{
    return ToString(scripts_viewer_py);
}

Scripts Manual()
{
    return ToString(data_manual_txt);
}

/*
 * Fusion Resources
 */

Scripts PlotROC_F()
{
    return ToString(src_r_plotROC_F_R);
}

Scripts FusionDataMixA()
{
    return ToString(data_FusQuin_MFU007_v013_csv);
}

Scripts FusionDataRef()
{
    return ToString(data_FusQuin_AFU004_v032_ref);
}

Scripts PlotExpress_F()
{
    return ToString(src_r_plotExpress_F_R);
}

/*
 * Ladder Resources
 */

Scripts PlotLadderAbund()
{
    return ToString(src_r_plotLadderAbund_R);
}

Scripts LadderDataMixA()
{
    return ToString(data_LadQuin_MLA014_v013_csv);
}

Scripts LadderDataMixB()
{
    return ToString(data_LadQuin_MLA016_v013_csv);
}

Scripts LadderDataMixAB()
{
    return ToString(data_LadQuin_MLA020_v013_csv);
}

/*
 * Metagenomics Resources
 */

Scripts PlotMFlod()
{
    return ToString(src_r_plotMFold_R);
}

Scripts MetaDataMix()
{
    return ToString(data_MetaQuin_MME023_v013_csv);
}

Scripts MetaDataBed()
{
    return ToString(data_MetaQuin_AME015_v032_bed);
}

/*
 * Transcriptome Resources
 */

Scripts SleuthR()
{
    return ToString(scripts_sleuth_R);
}

Scripts PlotFold()
{
    return ToString(src_r_plotFold_R);
}

Scripts PlotRAbundAbund()
{
    return ToString(src_r_plotRAbundAbund_R);
}

Scripts PlotMajor()
{
    return ToString(src_r_plotMajor_R);
}

Scripts PlotMA()
{
    return ToString(src_r_plotMA_R);
}

Scripts PlotTAbundAbund()
{
    return ToString(src_r_plotTAbundAbund_R);
}

Scripts PlotTROC()
{
    return ToString(src_r_plotTROC_R);
}

Scripts PlotLODR()
{
    return ToString(src_r_plotLODR_R);
}

Scripts TransStandGTF()
{
    return ToString(data_TransQuin_ATR001_v032_gtf);
}

Scripts TransDataMixA()
{
    return ToString(data_TransQuin_MTR002_v013_csv);
}

Scripts TransDataMixB()
{
    return ToString(data_TransQuin_MTR003_v013_csv);
}

Scripts TransDataMixF()
{
    return ToString(data_TransQuin_MTR005_v013_csv);
}

Scripts TransDataMixAB()
{
    return ToString(data_TransQuin_MTR004_v013_csv);
}

/*
 * Variant Resources
 */

Scripts PlotVAbundAbund()
{
    return ToString(src_r_plotVAbundAbund_R);
}

Scripts PlotAlleleReads()
{
    return ToString(src_r_plotAlleleReads_R);
}

Scripts PlotAlleleAllele()
{
    return ToString(src_r_plotAlleleAllele_R);
}

Scripts PlotROC_V()
{
    return ToString(src_r_plotROC_V_R);
}

Scripts PlotLODR_V()
{
    return ToString(src_r_plotLODR_V_R);
}

Scripts PlotDensity()
{
    return ToString(src_r_plotDensity_R);
}

Scripts PlotSubsample()
{
    return ToString(src_r_plotSubsample_R);
}

Scripts VarDataMixA()
{
    return ToString(data_VarQuin_MVA011_v013_csv);
}

Scripts VarDataMixF()
{
    return ToString(data_VarQuin_MVA012_v013_csv);
}

Scripts VarDataVCF()
{
    return ToString(data_VarQuin_AVA009_v032_vcf);
}

Scripts VarDataBed()
{
    return ToString(data_VarQuin_AVA017_v032_bed);
}
