#include <iostream>
#include <tclap/CmdLine.h>
#include <unistd.h>

#include "rna/aligner.hpp"
#include "rna/assembly.hpp"
#include "rna/abundance.hpp"
#include "rna/differential.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

using namespace Spike;

static void printManual()
{
    std::cout << std::ifstream("docs/manual.txt").rdbuf() << std::endl;
}

int main(int argc, char ** argv)
{
    typedef TCLAP::SwitchArg SArg;
    typedef TCLAP::ValueArg<std::string> VArg;

    printManual();
    return 0;
    
	std::cout << "Anquin by Garvan Institute - sequin data-analysis tool\n" << std::endl;

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
            Aligner::analyze(a.getValue());
        }
    }
    catch (TCLAP::ArgException &e)
    {
        throw;
    }
    
    return 0;
}