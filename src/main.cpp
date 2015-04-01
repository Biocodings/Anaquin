#include <limits>
#include <iostream>
#include <tclap/CmdLine.h>
#include "aligner.hpp"
#include "assembly.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

int main(int argc, char ** argv)
{
    typedef TCLAP::SwitchArg SArg;
    typedef TCLAP::ValueArg<std::string> VArg;

	std::cout << "Anquin by Garvan Institute - Chromosome data-analysis tool\n" << std::endl;

    try
    {
        TCLAP::CmdLine cmd("", ' ', "1.0");
        
        SArg t("t", "test", "Internal testing", cmd, false);
        VArg a("a", "align", "Assess alignment", false, "", "string", cmd);
        
        cmd.parse(argc, argv);

        if (t.getValue())
        {
            return Catch::Session().run(1, argv);
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
