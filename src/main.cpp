#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "rna/ralign.hpp"
#include "rna/assembly.hpp"
#include "rna/abundance.hpp"
#include "rna/differential.hpp"

#include "dna/dalign.hpp"
#include "dna/structural.hpp"

#include "parsers/parser_csv.hpp"
#include "writers/path_writer.hpp"
#include "parsers/parser_sequins.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#define RNA_MIX_PATH "data/RNA"
#define DNA_MIX_PATH "data/DNA"

#define MODE_VER   'v'
#define MODE_TEST  265
#define MODE_RNA   266
#define MODE_DNA   267
#define MODE_META  268

#define CMD_RESTRICTS    280
#define CMD_SEQS         281
#define CMD_SEQUENCING   282
#define CMD_ALIGN        283
#define CMD_ASSEMBLY     284
#define CMD_ABUNDANCE    285
#define CMD_DIFFERENTIAL 286
#define CMD_VARIATION    287

using namespace Spike;

/*
 * Variables used in argument parsing
 */

// The path that output files are written
static std::string _output;

// The sequins that have been restricted
static std::vector<SequinID> filtered;

// The mode for the command such as RNA and DNA
static int _mode = 0;

// The command for the specified mode
static int _cmd = 0;

/*
 * Argument options
 */

static const char *short_options = "o:";

static const struct option long_options[] =
{
    { "test",         no_argument,       0, MODE_TEST },
    { "version",      no_argument,       0, MODE_VER  },
    { "rna",          required_argument, 0, MODE_DNA  },
    { "dna",          required_argument, 0, MODE_RNA  },
    { "meta",         required_argument, 0, MODE_META },

    { "restrict",     required_argument, 0, CMD_RESTRICTS    },
    { "l",            no_argument,       0, CMD_SEQS         },
    { "seqs",         no_argument,       0, CMD_SEQS         },
    { "al",           required_argument, 0, CMD_ALIGN        },
    { "align",        required_argument, 0, CMD_ALIGN        },
    { "as",           required_argument, 0, CMD_ASSEMBLY     },
    { "assembly",     required_argument, 0, CMD_ASSEMBLY     },
    { "ab",           required_argument, 0, CMD_ABUNDANCE    },
    { "abundance",    required_argument, 0, CMD_ABUNDANCE    },
    { "df",           required_argument, 0, CMD_DIFFERENTIAL },
    { "differential", required_argument, 0, CMD_DIFFERENTIAL },
    { "var",          required_argument, 0, CMD_VARIATION    },
    { "variation",    required_argument, 0, CMD_VARIATION    },

    {0,  0, 0,  0 }
};

static void print_usage()
{
    std::cout << std::ifstream("docs/manual.txt").rdbuf() << std::endl;
}

static void print_version()
{
    std::cout << "Version 1.0. Garvan Institute, copyright 2015." << std::endl;
}

static void print_sequins(const std::string &file)
{
    const std::string format = "%1%  %2%  %3%  %4%  %5%  %6%  %7%  %8%";

    std::cout << (boost::format(format) % "r_name"
                                        % "v_name"
                                        % "r_ratio"
                                        % "v_ratio"
                                        % "r_con"
                                        % "v_con").str() << std::endl;

    ParserCSV::parse(file, [&](const Fields &fields)
    {
        std::cout << (boost::format(format) % fields[0]
                                            % fields[3]
                                            % fields[8]
                                            % fields[9]
                                            % fields[11]
                                            % fields[12]).str() << std::endl;
    });
}

template <typename Analyzer, typename Level> void analyze(const std::string &file, Level lv)
{
    typename Analyzer::Options o;

    // If nothing is provided, default directory (empty string) will be used
    o.writer = std::shared_ptr<PathWriter>(new PathWriter(_output));

    o.level = lv;

    std::cout << "-----------------------------------------" << std::endl;
    //std::cout << "Analyze " << Analyzer::name() << " data-analyzer..." << std::endl;

    Analyzer::analyze(file, o);

    std::cout << "Completed." << std::endl;
}

/*
 * Some commands such as RNA-abdunance takes an input file at the gene-level or isoform-level.
 * If the purpose isn't explicity specified, maybe we can automatically detect its usage.
 */

template <typename Level> static Level detect(const std::string &file)
{
    const bool found_gene = file.find("gene") != std::string::npos;
    const bool found_isoform = file.find("isoform") != std::string::npos;

    if (found_gene && !found_isoform)
    {
        return Level::Gene;
    }
    else if (!found_gene && found_isoform)
    {
        return Level::Isoform;
    }

    throw std::runtime_error("Unknown type. Have you specified the level?");
}

static int parse_options(int argc, char ** argv)
{
    int next, index;

    do
    {
        switch (next = getopt_long_only(argc, argv, short_options, long_options, &index))
        {
            case 'o': { _output = optarg; break; }

            case MODE_VER:
            case MODE_DNA:
            case MODE_RNA:
            case MODE_META:
            case MODE_TEST:
            {
                _mode = next;
                break;
            }

            default:
            {
                _cmd = next;
                break;
            }
        }
    } while (next != -1);

    /*
     * The arguments have been parsed. It's time to check if it's a valid command.
     */
    
    // It's an error to provide anything for a void mode
    if ((_mode == MODE_TEST || _mode == MODE_VER) && (_output.empty()))
    {
        print_usage();
    }
    else
    {
        switch (_mode)
        {
            case MODE_TEST:
            {
                return Catch::Session().run(1, argv);
            }
                
            case MODE_VER:
            {
                print_version();
                break;
            }

            case MODE_RNA:
            {
                if (_cmd != CMD_SEQUENCING &&
                    _cmd != CMD_ALIGN      &&
                    _cmd != CMD_ASSEMBLY   &&
                    _cmd != CMD_ABUNDANCE  &&
                    _cmd != CMD_DIFFERENTIAL)
                {
                    print_usage();
                }
                else
                {
                    switch (_cmd)
                    {
                        case CMD_SEQS:         { print_sequins(RNA_MIX_PATH);                   break; }
                        case CMD_ALIGN:        { analyze<RAlign>(optarg, RAlign::Base);     break; }
                        case CMD_ASSEMBLY:     { analyze<Assembly>(optarg, Assembly::Base);     break; }
                        case CMD_ABUNDANCE:
                        {
                            analyze<Abundance>(optarg, detect<Abundance::Level>(optarg));
                            break;
                        }

                        case CMD_DIFFERENTIAL:
                        {
                            analyze<Differential>(optarg, detect<Differential::Level>(optarg));
                            break;
                        }
                    }
                }
                
                break;
            }
                
            case MODE_DNA:
            {
                if (_cmd != CMD_SEQUENCING &&
                    _cmd != CMD_ALIGN      &&
                    _cmd != CMD_VARIATION)
                {
                    print_usage();
                }
                else
                {
                    switch (_cmd)
                    {
                        case CMD_SEQS:      { print_sequins(RNA_MIX_PATH);                   break; }
                        case CMD_ALIGN:     { analyze<DAlign>(optarg, DAlign::Base);         break; }
                        case CMD_VARIATION: { analyze<Structural>(optarg, Structural::Base); break; }
                    }
                }

                break;
            }

            default:
            {
                assert(false);
            }
        }
    }

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