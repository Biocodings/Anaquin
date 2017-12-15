#include <string>
#include <algorithm>

#include "resources/anaquin.txt"

#include "resources/RnaAlign.txt"
#include "resources/RnaAssembly.txt"
#include "resources/RnaSubsample.txt"
#include "resources/RnaExpression.txt"
#include "resources/RnaFoldChange.txt"

#include "resources/plotFold.R"
#include "resources/plotTROC.R"
#include "resources/plotTLODR.R"
#include "resources/plotLinear.R"
#include "resources/plotLogistic.R"

typedef std::string Scripts;

#define ToString(x) std::string(reinterpret_cast<char*>(x))

Scripts Manual() { return ToString(data_manuals_anaquin_txt); }

Scripts PlotFold()       { return ToString(src_r_plotFold_R);       }
Scripts PlotLinear()     { return ToString(src_r_plotLinear_R);     }
Scripts PlotLogistic()   { return ToString(src_r_plotLogistic_R);   }

Scripts RnaAlign()      { return ToString(data_manuals_RnaAlign_txt);      }
Scripts RnaSubsample()  { return ToString(data_manuals_RnaSubsample_txt);  }
Scripts RnaAssembly()   { return ToString(data_manuals_RnaAssembly_txt);   }
Scripts RnaExpression() { return ToString(data_manuals_RnaExpression_txt); }
Scripts RnaFoldChange() { return ToString(data_manuals_RnaFoldChange_txt); }

Scripts PlotTROC()  { return ToString(src_r_plotTROC_R);  }
Scripts PlotTLODR() { return ToString(src_r_plotTLODR_R); }
