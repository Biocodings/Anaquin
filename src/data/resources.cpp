#include <string>

/*
 * Anaquin R-bioconductor
 */

#include "resources/plotMA.R"
#include "resources/plotSplice.R"
#include "resources/plotScatterPool.R"

/*
 * Scripts
 */

#include "resources/viewer.py"
#include "resources/manual.txt"

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

#include "resources/plotROC_T.R"
#include "resources/plotLODR_T.R"
#include "resources/plotScatter_T.R"

#include "resources/ATR001.v032.gtf"
#include "resources/MTR002.v013.csv"
#include "resources/MTR003.v013.csv"
#include "resources/MTR004.v013.csv"
#include "resources/MTR005.v013.csv"

/*
 * Metagenomics Resources
 */

#include "resources/MME023.v013.csv"
#include "resources/AME015.v032.bed"

/*
 * Variant Resources
 */

#include "resources/plotROC_V.R"
#include "resources/plotAllCov.R"
#include "resources/plotLODR_V.R"
#include "resources/plotAlleleF.R"
#include "resources/plotDensity.R"
#include "resources/plotSubsample.R"

#include "resources/AVA009.v032.vcf"
#include "resources/MVA011.v013.csv"
#include "resources/MVA012.v013.csv"
#include "resources/AVA017.v032.bed"

#define ToString(x) std::string(reinterpret_cast<char*>(x))

typedef std::string Scripts;

Scripts ViewerScript()
{
    return ToString(scripts_viewer_py);
}

Scripts Manual()
{
    return ToString(data_manual_txt);
}

/*
 * Scripts for Anaquin Bioconductor
 */

Scripts PlotScatterPool()
{
    return ToString(src_r_plotScatterPool_R);
}

Scripts PlotSplice()
{
    return ToString(src_r_plotSplice_R);
}

Scripts PlotMA()
{
    return ToString(src_r_plotMA_R);
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
    return ToString(data_fusion_MFU007_v013_csv);
}

Scripts FusionDataRef()
{
    return ToString(data_fusion_AFU004_v032_ref);
}

Scripts PlotExpress_F()
{
    return ToString(src_r_plotExpress_F_R);
}

Scripts PlotAllCov()
{
    return ToString(src_r_plotAllCov_R);
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
    return ToString(data_ladder_MLA014_v013_csv);
}

Scripts LadderDataMixB()
{
    return ToString(data_ladder_MLA016_v013_csv);
}

Scripts LadderDataMixAB()
{
    return ToString(data_ladder_MLA020_v013_csv);
}

/*
 * Metagenomics Resources
 */

Scripts MetaDataMix()
{
    return ToString(data_meta_MME023_v013_csv);
}

Scripts MetaDataBed()
{
    return ToString(data_meta_AME015_v032_bed);
}

/*
 * Transcriptome Resources
 */

Scripts PlotScatter_T()
{
    return ToString(src_r_plotScatter_T_R);
}

Scripts PlotROC_T()
{
    return ToString(src_r_plotROC_T_R);
}

Scripts PlotLODR_T()
{
    return ToString(src_r_plotLODR_T_R);
}

Scripts TransStandGTF()
{
    return ToString(data_trans_ATR001_v032_gtf);
}

Scripts TransDataMixA()
{
    return ToString(data_trans_MTR002_v013_csv);
}

Scripts TransDataMixB()
{
    return ToString(data_trans_MTR003_v013_csv);
}

Scripts TransDataMixF()
{
    return ToString(data_trans_MTR005_v013_csv);
}

Scripts TransDataMixAB()
{
    return ToString(data_trans_MTR004_v013_csv);
}

/*
 * Variant Resources
 */

Scripts PlotAlleleF()
{
    return ToString(src_r_plotAlleleF_R);
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
