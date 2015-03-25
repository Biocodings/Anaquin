#include <limits>
#include <iostream>
#include "gtest/gtest.h"
#include "AlignerAnalyst.hpp"
#include "AssemblyAnalyst.hpp"
#include <tclap/CmdLine.h>

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
            AlignerAnalyst::base(a.getValue(), Sequins(), 1000);
        }
    }
    catch (TCLAP::ArgException &e)
    {
        throw;
    }
    
    return 0;
}
