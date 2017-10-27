#include <string>
#include <algorithm>

#include "resources/RKmersForVarKStats.tsv"
#include "resources/structural_trimming.bed"

#include "resources/anaquin.txt"
#include "resources/VarCopy.txt"
#include "resources/VarTrim.txt"
#include "resources/VarKmer.txt"
#include "resources/VarAlign.txt"
#include "resources/VarKStats.txt"
#include "resources/VarMutation.txt"
#include "resources/VarPartition.txt"
#include "resources/VarStructure.txt"

#include "resources/RnaAlign.txt"
#include "resources/RnaReport.txt"
#include "resources/RnaAssembly.txt"
#include "resources/RnaSubsample.txt"
#include "resources/RnaExpression.txt"
#include "resources/RnaFoldChange.txt"

#include "resources/MetaAssembly.txt"
#include "resources/MetaCoverage.txt"
#include "resources/MetaSubsample.txt"

#include "resources/plotCNV.R"
#include "resources/plotFold.R"
#include "resources/plotTROC.R"
#include "resources/plotVGROC.R"
#include "resources/plotVCROC.R"
#include "resources/plotTLODR.R"
#include "resources/plotAllele.R"
#include "resources/plotLinear.R"
#include "resources/plotKAllele.R"
#include "resources/plotConjoint.R"
#include "resources/plotLogistic.R"

typedef std::string Scripts;

#define ToString(x) std::string(reinterpret_cast<char*>(x))

Scripts Manual() { return ToString(data_manuals_anaquin_txt); }

Scripts StructBED()    { return ToString(data_structural_trimming_bed); }
Scripts KMVarKStats()  { return ToString(data_RKmersForVarKStats_tsv);  }

Scripts PlotFold()     { return ToString(src_r_plotFold_R);     }
Scripts PlotCNV()      { return ToString(src_r_plotCNV_R);      }
Scripts PlotLinear()   { return ToString(src_r_plotLinear_R);   }
Scripts PlotConjoint() { return ToString(src_r_plotConjoint_R); }
Scripts PlotAllele()   { return ToString(src_r_plotAllele_R);   }
Scripts PlotKAllele()  { return ToString(src_r_plotKAllele_R);  }
Scripts PlotLogistic() { return ToString(src_r_plotLogistic_R); }

Scripts RnaAlign()      { return ToString(data_manuals_RnaAlign_txt);      }
Scripts RnaReport()     { return ToString(data_manuals_RnaReport_txt);     }
Scripts RnaSubsample()  { return ToString(data_manuals_RnaSubsample_txt);  }
Scripts RnaAssembly()   { return ToString(data_manuals_RnaAssembly_txt);   }
Scripts RnaExpression() { return ToString(data_manuals_RnaExpression_txt); }
Scripts RnaFoldChange() { return ToString(data_manuals_RnaFoldChange_txt); }

Scripts VarTrim()       { return ToString(data_manuals_VarTrim_txt);      }
Scripts VarCopy()       { return ToString(data_manuals_VarCopy_txt);      }
Scripts VarAlign()      { return ToString(data_manuals_VarAlign_txt);     }
Scripts VarKmer()       { return ToString(data_manuals_VarKmer_txt);      }
Scripts VarKStats()     { return ToString(data_manuals_VarKStats_txt);    }
Scripts VarMutation()   { return ToString(data_manuals_VarMutation_txt);  }
Scripts VarPartition()  { return ToString(data_manuals_VarPartition_txt); }
Scripts VarStructure()  { return ToString(data_manuals_VarStructure_txt); }

Scripts MetaCoverage()  { return ToString(data_manuals_MetaCoverage_txt);   }
Scripts MetaSubsample() { return ToString(data_manuals_MetaSubsample_txt);  }
Scripts MetaAssembly()  { return ToString(data_manuals_MetaAssembly_txt);   }

Scripts PlotTROC()  { return ToString(src_r_plotTROC_R);  }
Scripts PlotTLODR() { return ToString(src_r_plotTLODR_R); }

Scripts PlotVGROC() { return ToString(src_r_plotVGROC_R); }
Scripts PlotVCROC() { return ToString(src_r_plotVCROC_R); }

