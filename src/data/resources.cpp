#include <string>
#include <algorithm>

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

#include "resources/plotFROC.R"
#include "resources/plotFExpress.R"

#include "resources/AFU004.v032.ref"
#include "resources/AFU005.v032.bed"
#include "resources/MFU007.v013.csv"

/*
 * Ladder Resources
 */

#include "resources/plotLExpress.R"

#include "resources/MLA014.v013.csv"
#include "resources/MLA016.v013.csv"
#include "resources/MLA020.v013.csv"

/*
 * Transcriptome Resources
 */

#include "resources/plotTMA.R"
#include "resources/plotTFold.R"
#include "resources/plotTROC.R"
#include "resources/plotTMinor.R"
#include "resources/plotTLODR.R"
#include "resources/plotTMultiple.R"
#include "resources/plotTExpress.R"

#include "resources/ATR001.v032.gtf"
#include "resources/MTR002.v013.csv"
#include "resources/MTR003.v013.csv"
#include "resources/MTR004.v013.csv"
#include "resources/MTR005.v013.csv"

/*
 * Metagenomics Resources
 */

#include "resources/plotMFold.R"
#include "resources/plotMExpress.R"
#include "resources/plotMFraction.R"

#include "resources/MME023.v013.csv"
#include "resources/AME015.v032.bed"

/*
 * Variant Resources
 */

#include "resources/plotVROC.R"
#include "resources/plotVProb.R"
#include "resources/plotVAllele.R"
#include "resources/plotVDensity.R"
#include "resources/plotVExpress.R"
#include "resources/plotVSubsample.R"
#include "resources/plotVAlleleReads.R"

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
    return ToString(data_FusQuin_AFU004_v032_ref);
}

Scripts PlotFExpress()
{
    return ToString(src_r_plotFExpress_R);
}

/*
 * Ladder Resources
 */

Scripts PlotLExpress()
{
    return ToString(src_r_plotLExpress_R);
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

Scripts PlotMFraction()
{
    return ToString(src_r_plotMFraction_R);
}

Scripts PlotMExpress()
{
    return ToString(src_r_plotMExpress_R);
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

Scripts PlotTFold()
{
    return ToString(src_r_plotTFold_R);
}

Scripts PlotTMultiple()
{
    return ToString(src_r_plotTMultiple_R);
}

Scripts PlotTMinor()
{
    return ToString(src_r_plotTMinor_R);
}

Scripts PlotTMA()
{
    return ToString(src_r_plotTMA_R);
}

Scripts PlotTExpress()
{
    return ToString(src_r_plotTExpress_R);
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

Scripts PlotVExpress()
{
    return ToString(src_r_plotVExpress_R);
}

Scripts PlotVAlleleReads()
{
    return ToString(src_r_plotVAlleleReads_R);
}

Scripts PlotVAllele()
{
    return ToString(src_r_plotVAllele_R);
}

Scripts PlotVROC()
{
    return ToString(src_r_plotVROC_R);
}

Scripts PlotVProb()
{
    return ToString(src_r_plotVProb_R);
}

Scripts PlotVDensity()
{
    return ToString(src_r_plotVDensity_R);
}

Scripts PlotVSubsample()
{
    return ToString(src_r_plotVSubsample_R);
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
