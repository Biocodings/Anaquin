#include <string>
#include <algorithm>

/*
 * Scripts
 */

#include "resources/sleuth.R"
#include "resources/viewer.py"
#include "resources/manual.txt"
#include "resources/reports.py"

#include "resources/plotScatter.R"
#include "resources/plotSensitivity.R"
#include "resources/plotMultScatter.R"

/*
 * Fusion Resources
 */

#include "resources/plotFROC.R"
#include "resources/plotFFold.R"
#include "resources/plotFNormal.R"
#include "resources/plotFFusion.R"

#include "resources/AFU004.v032.bed"
#include "resources/AFU005.v032.bed"
#include "resources/MFU007.v013.csv"

/*
 * Ladder Resources
 */

#include "resources/plotLNorm.R"

#include "resources/MLA014.v013.csv"
#include "resources/MLA016.v013.csv"
#include "resources/MLA020.v013.csv"

/*
 * Transcriptome Resources
 */

#include "resources/plotTMA.R"
#include "resources/plotTROC.R"
#include "resources/plotTLODR.R"
#include "resources/plotTMinor.R"

#include "resources/ATR001.v032.gtf"
#include "resources/MTR002.v013.csv"
#include "resources/MTR003.v013.csv"
#include "resources/MTR004.v013.csv"
#include "resources/MTR005.v013.csv"

/*
 * Metagenomics Resources
 */

#include "resources/plotMKMer.R"
#include "resources/plotMFold.R"
#include "resources/plotMReads.R"
#include "resources/plotMKAbund.R"
#include "resources/plotMAssembly.R"

#include "resources/MME023.v013.csv"
#include "resources/AME015.v032.bed"

/*
 * Variant Resources
 */

#include "resources/plotVLOD.R"
#include "resources/plotVROC1.R"
#include "resources/plotVROC2.R"
#include "resources/plotVDensity.R"

#include "resources/AVA009.v032.vcf"
#include "resources/MVA011.v013.csv"
#include "resources/MVA012.v013.csv"
#include "resources/AVA017.v032.bed"

typedef std::string Scripts;

#define ToString(x) std::string(reinterpret_cast<char*>(x))

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

Scripts PlotScatter()
{
    return ToString(src_r_plotScatter_R);
}

Scripts plotMultScatter()
{
    return ToString(src_r_plotMultScatter_R);
}

Scripts PlotSensitivity()
{
    return ToString(src_r_plotSensitivity_R);
}

/*
 * Fusion Resources
 */

Scripts PlotFROC()
{
    return ToString(src_r_plotFROC_R);
}

Scripts FusionDataMixA()
{
    return ToString(data_FusQuin_MFU007_v013_csv);
}

Scripts FusionDataRef()
{
    return ToString(data_FusQuin_AFU004_v032_bed);
}

Scripts PlotFNormal()
{
    return ToString(src_r_plotFNormal_R);
}

Scripts PlotFFusion()
{
    return ToString(src_r_plotFFusion_R);
}

Scripts PlotFFold()
{
    return ToString(src_r_plotFFold_R);
}

/*
 * Ladder Resources
 */

Scripts PlotLNorm()
{
    return ToString(src_r_plotLNorm_R);
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

Scripts PlotMReads()
{
 	return ToString(src_r_plotMReads_R);
}

Scripts PlotMAssembly()
{
	return ToString(src_r_plotMAssembly_R);
}

Scripts PlotMKMer()
{
    return ToString(src_r_plotMKMer_R);
}

Scripts PlotMKAbund()
{
    return ToString(src_r_plotMKAbund_R);
}

Scripts PlotMFold()
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

Scripts PlotTMinor()
{
    return ToString(src_r_plotTMinor_R);
}

Scripts PlotTMA()
{
    return ToString(src_r_plotTMA_R);
}

Scripts PlotTROC()
{
    return ToString(src_r_plotTROC_R);
}

Scripts PlotTLODR()
{
    return ToString(src_r_plotTLODR_R);
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

Scripts PlotVLOD()
{
    return ToString(src_r_plotVLOD_R);
}

Scripts PlotVROC1()
{
    return ToString(src_r_plotVROC1_R);
}

Scripts PlotVROC2()
{
    return ToString(src_r_plotVROC2_R);
}

Scripts PlotVDensity()
{
    return ToString(src_r_plotVDensity_R);
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
