#include <limits>
#include <iostream>
#include <tclap/CmdLine.h>
#include "gtest/gtest.h"
#include "aligner.hpp"
#include "assembly.hpp"

int main(int argc, char ** argv)
{
    typedef TCLAP::SwitchArg SArg;
    typedef TCLAP::ValueArg<std::string> VArg;
    
    try
    {
        TCLAP::CmdLine cmd("", ' ', "1.0");
        
        SArg t("t", "test", "Internal testing", cmd, false);
        VArg a("a", "align", "Assess alignment", false, "", "string", cmd);
        
        cmd.parse(argc, argv);

        if (t.getValue())
        {
            ::testing::InitGoogleTest(&argc, argv);
            return RUN_ALL_TESTS();
        }
        else if (!a.getValue().empty())
        {
            Aligner::base(a.getValue(), Sequins(), 1000);
        }
    }
    catch (TCLAP::ArgException &e)
    {
        throw;
    }
    
    return 0;
}
