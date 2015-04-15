#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "rna/aligner.hpp"
#include "rna/assembly.hpp"
#include "rna/abundance.hpp"
#include "rna/differential.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#define OPT_ALIGN        265
#define OPT_ASSEMBLY     266
#define OPT_ABUNDANCE    267
#define OPT_DIFFERENTIAL 268

using namespace Spike;

/*
 * Variables used in argument parsing
 */

// The path that output files are written
static std::string output;

static const char *short_options = "o:";

static const struct option long_options[] =
{
    {"test",     no_argument, 0, 't' },
    {"version",  no_argument, 0, 'v' },

    {"al",           required_argument, 0, OPT_ALIGN        },
    {"align",        required_argument, 0, OPT_ALIGN        },
    {"as",           required_argument, 0, OPT_ASSEMBLY     },
    {"assembly",     required_argument, 0, OPT_ASSEMBLY     },
    {"ab",           required_argument, 0, OPT_ABUNDANCE    },
    {"abundance",    required_argument, 0, OPT_ABUNDANCE    },
    {"df",           required_argument, 0, OPT_DIFFERENTIAL },
    {"differential", required_argument, 0, OPT_DIFFERENTIAL },

    {0,  0, 0,  0 }
};

static void print_usage()
{
    std::cout << std::ifstream("docs/manual.txt").rdbuf() << std::endl;
}

static void print_version()
{
    std::cout << "Version 1.0. Garvan Institute, copyright 2015" << std::endl;
}

template <typename Analyzer> void analyze(const std::string &file)
{
    typename Analyzer::Options o;

    // If nothing is provided, output will simply be an empty string
    o.output = output;

    Analyzer::analyze(file);
}

static int parse_options(int argc, char ** argv)
{
    int next, index;

    do
    {
        switch (next = getopt_long_only(argc, argv, short_options, long_options, &index))
        {
            case 'o': { output = optarg; break; }

            case 'v': { print_version(); break; }

            case 't':
            {
                return Catch::Session().run(1, argv);
            }

            /*
             * RNA options
             */

            case OPT_ALIGN:
            {
                analyze<Aligner>(optarg);
                break;
            }

            case OPT_ASSEMBLY:
            {
                analyze<Assembly>(optarg);
                break;
            }

            case OPT_ABUNDANCE:
            {
                analyze<Abundance>(optarg);
                break;
            }

            case OPT_DIFFERENTIAL:
            {
                analyze<Differential>(optarg);
                break;
            }

            /*
             * DNA options
             */

            /*
             * Metagenomics options
             */
                
            default:
            {
                break;
            }
        }
    } while (next != -1);

    return 0;
}

int main(int argc, char ** argv)
{
    if (argc == 1)
    {
        print_usage();
        return 0;
    }
    else
    {
        return parse_options(argc, argv);
    }
}