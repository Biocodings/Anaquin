#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "rna/aligner.hpp"
#include "rna/assembly.hpp"
#include "rna/abundance.hpp"
#include "rna/differential.hpp"

#include "dna/structural.hpp"

#include "writers/path_writer.hpp"
#include "parsers/parser_sequins.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#define O_RNA          265
#define O_DNA          266
#define O_META         267
#define O_SEQUENCING   268
#define O_ALIGN        269
#define O_ASSEMBLY     270
#define O_ABUNDANCE    271
#define O_DIFFERENTIAL 272
#define O_VARIATION    273

using namespace Spike;

/*
 * Variables used in argument parsing
 */

// The path that output files are written
static std::string output;

// The sequins that have been restricted
static std::vector<SequinID> filtered;

// O_RNA, O_DNA or O_META
static int mode = 0;

/*
 * Argument options
 */

static const char *short_options = "o:";

static const struct option long_options[] =
{
    { "test",    no_argument, 0, 't' },
    { "version", no_argument, 0, 'v' },

    { "restrict", required_argument, 0, 'r' },

    { "rna",          required_argument, 0, O_RNA  },
    { "dna",          required_argument, 0, O_DNA  },
    { "meta",         required_argument, 0, O_META },

    { "al",           required_argument, 0, O_ALIGN        },
    { "align",        required_argument, 0, O_ALIGN        },
    { "as",           required_argument, 0, O_ASSEMBLY     },
    { "assembly",     required_argument, 0, O_ASSEMBLY     },
    { "ab",           required_argument, 0, O_ABUNDANCE    },
    { "abundance",    required_argument, 0, O_ABUNDANCE    },
    { "df",           required_argument, 0, O_DIFFERENTIAL },
    { "differential", required_argument, 0, O_DIFFERENTIAL },
    { "var",          required_argument, 0, O_VARIATION    },
    { "variation",    required_argument, 0, O_VARIATION    },

    /*
     * Metagenomics options
     */

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

template <typename Analyzer, typename Level> void analyze(const std::string &file, Level lv)
{
    typename Analyzer::Options o;

    // If nothing is provided, default directory (empty string) will be used
    o.writer = std::shared_ptr<PathWriter>(new PathWriter(output));

    o.level = lv;
    
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Analyze " << Analyzer::name() << " data-analyzer..." << std::endl;

    Analyzer::analyze(file, o);

    std::cout << "Completed." << std::endl;
}

/*
 * Some commands such as RNA-abdunance takes an input file at the gene-level or isoform-level.
 * If the purpose isn't explicity specified, maybe we can automatically detect its usage.
 */

static RNALevel detect(const std::string &file)
{
    const bool found_gene = file.find("gene") != std::string::npos;
    const bool found_isoform = file.find("isoform") != std::string::npos;
    
    if (found_gene && !found_isoform)
    {
        std::cout << "Gene input detected" << std::endl;
        return LevelGene;
    }
    else if (!found_gene && found_isoform)
    {
        std::cout << "Isoform input detected" << std::endl;
        return LevelIsoform;
    }

    throw std::runtime_error("Unknown file usage");
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

            case 'r':
            {
                filtered = ParserSequins::parse(optarg);
                break;
            }
                
            case O_RNA:
            {
                if (mode != 0)
                {
                    throw std::runtime_error("Unknown command");
                }

                mode = next;
                break;
            }

            /*
             * RNA options
             */

            case O_ALIGN:
            {
                analyze<Aligner>(optarg, Aligner::LevelBase);
                break;
            }

            case O_ASSEMBLY:
            {
                analyze<Assembly>(optarg, Assembly::Mode::Assembly_Base);
                break;
            }

            case O_ABUNDANCE:
            {
                analyze<Abundance>(optarg, detect(optarg));
                break;
            }

            case O_DIFFERENTIAL:
            {
                analyze<Differential>(optarg, detect(optarg));
                break;
            }

            /*
             * DNA options
             */

            case O_VARIATION:
            {
                analyze<Structural>(optarg, Structural::LevelBase);
                break;
            }

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