#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "rna/r_align.hpp"
#include "rna/assembly.hpp"
#include "rna/abundance.hpp"
#include "rna/differential.hpp"

#include "dna/d_align.hpp"
#include "dna/d_sequence.hpp"
#include "dna/structural.hpp"

#include "meta/denovo.hpp"
#include "meta/m_align.hpp"
#include "meta/m_sequence.hpp"

#include "parsers/parser_csv.hpp"
#include "writers/path_writer.hpp"
#include "parsers/parser_sequins.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#define RNA_MIX_PATH "data/RNA"
#define DNA_MIX_PATH "data/DNA"

#define CMD_VER   'v'
#define CMD_TEST  't'
#define CMD_RNA   265
#define CMD_DNA   266
#define CMD_META  267

#define MODE_RESTRICTS    280
#define MODE_SEQS         281
#define MODE_SEQUENCING   282
#define MODE_ALIGN        283
#define MODE_ASSEMBLY     284
#define MODE_ABUNDANCE    285
#define MODE_DIFFERENTIAL 286
#define MODE_VARIATION    287
#define MODE_DE_NOVO      288

using namespace Spike;

/*
 * Variables used in argument parsing
 */

// The path that output files are written
static std::string _output;

// The sequins that have been restricted
static std::vector<SequinID> filtered;

static int _cmd  = 0;
static int _mode = 0;

// The operand for _cmd
static std::string _opt;

/*
 * Argument options
 */

static const char *short_options = "o";

static const struct option long_options[] =
{
    { "t",    no_argument,       0, CMD_TEST },
    { "v",    no_argument,       0, CMD_VER  },
    { "rna",  required_argument, 0, CMD_RNA  },
    { "dna",  required_argument, 0, CMD_DNA  },
    { "meta", required_argument, 0, CMD_META },

    { "f",    required_argument, 0, MODE_RESTRICTS    },
    { "l",    no_argument,       0, MODE_SEQS         },
    { "se",   required_argument, 0, MODE_SEQS         },
    { "al",   required_argument, 0, MODE_ALIGN        },
    { "as",   required_argument, 0, MODE_ASSEMBLY     },
    { "ab",   required_argument, 0, MODE_ABUNDANCE    },
    { "df",   required_argument, 0, MODE_DIFFERENTIAL },
    { "va",   required_argument, 0, MODE_VARIATION    },
    { "de",   required_argument, 0, MODE_DE_NOVO      },

    {0,  0, 0,  0 }
};

typedef int Command;

static int invalid_cmd(Command cmd = 0)
{
    const std::map<Command, std::string> mapper =
    {
        { CMD_TEST,  "-t"    },
        { CMD_VER,   "-v"    },
        { CMD_RNA,   "-rna"  },
        { CMD_DNA,   "-dna"  },
        { CMD_META,  "-meta" },
    };

    if (cmd)
    {
        std::cerr << "Invalid command: " << mapper.at(cmd) << std::endl;
    }
    
    return 1;
}

static void print_usage()
{
    std::cout << std::ifstream("docs/manual.txt").rdbuf() << std::endl;
}

static void print_version()
{
    std::cout << "Version 1.0. Garvan Institute, copyright 2015." << std::endl;
}

static void print_rna()
{
    const auto &r = Standard::instance();
    const std::string format = "%1%  %2%  %3%  %4%  %5%  %6%  %7%";

    auto f = [&](const SequinsMap &seqs)
    {
        std::cout << (boost::format(format) % "r"
                                            % "v"
                                            % "fold"
                                            % "r_con"
                                            % "v_con"
                                            % "r_norm"
                                            % "v_norm").str() << std::endl;

        for (std::size_t i = A; i <= D; i++)
        {
            for (auto j : seqs)
            {
                const auto &p = j.second;
                
                if (p.grp == i)
                {
                    std::cout << (boost::format(format) % p.r.id
                                                        % p.v.id
                                                        % p.fold
                                                        % p.r.raw
                                                        % p.v.raw
                                                        % p.r.fpkm
                                                        % p.r.fpkm).str() << std::endl;
                }
            }
        }
    };

    std::cout << "\n----------------------------- A -----------------------------" << std::endl;
    f(r.r_seqs_gA);
    std::cout << "\n----------------------------- B -----------------------------" << std::endl;
    f(r.r_seqs_gB);
}

void print_dna()
{
    // Empty Implementation
}

void print_meta()
{
    // Empty Implementation    
}

template <typename Analyzer, typename Level> void analyze(const std::string &file, Level lv)
{
    typename Analyzer::Options o;

    // If nothing is provided, default directory (empty string) will be used
    o.writer = std::shared_ptr<PathWriter>(new PathWriter(_output));

    o.level = lv;
    
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "------------- Sequin Analysis -----------" << std::endl;
    std::cout << "-----------------------------------------" << std::endl << std::endl;

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
    /*
     * It's quite tricky and complicated to parse all commands by getopt_long_only().
     * Check for the commands that don't require '-' prefix.
     */

    if (argc <= 1)
    {
        print_usage();
    }
    else
    {
        const auto tmp = std::string(argv[1]);

        /*
         * The following commands will be ignored by getopt_long_only()
         */
        
        if (tmp == "rna")  { _cmd = CMD_RNA;  }
        if (tmp == "dna")  { _cmd = CMD_DNA;  }
        if (tmp == "meta") { _cmd = CMD_META; }
    }
    
    int next, index;

    while ((next = getopt_long_only(argc, argv, short_options, long_options, &index)) != -1)
    {
        switch (next)
        {
            case 'o': { _output = optarg; break; }

            case CMD_VER:
            case CMD_DNA:
            case CMD_RNA:
            case CMD_META:
            case CMD_TEST:
            {
                if (_cmd != 0)
                {
                    return invalid_cmd(_cmd);
                }
                
                _cmd = next;
                break;
            }
                
            default:
            {
                if (_mode != 0 || !_opt.empty())
                {
                    return invalid_cmd();
                }
                
                _mode = next;
                _opt  = optarg ? optarg : _opt;

                break;
            }
        }
    }

    if (_cmd == 0)
    {
        print_usage();
    }
    else if ((_cmd == CMD_TEST || _cmd == CMD_VER) && (!_output.empty() || _mode != 0 || !_opt.empty()))
    {
        invalid_cmd(_cmd);
    }
    else
    {
        switch (_cmd)
        {
            case CMD_TEST:
            {
                return Catch::Session().run(1, argv);
            }

            case CMD_VER:
            {
                print_version();
                break;
            }
                
            case CMD_RNA:
            {
                if (_mode != MODE_SEQS       &&
                    _mode != MODE_SEQUENCING &&
                    _mode != MODE_ALIGN      &&
                    _mode != MODE_ASSEMBLY   &&
                    _mode != MODE_ABUNDANCE  &&
                    _mode != MODE_DIFFERENTIAL)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_SEQS:      { print_rna();                             break; }
                        case MODE_ALIGN:     { analyze<RAlign>(_opt, RAlign::Base);     break; }
                        case MODE_ASSEMBLY:  { analyze<Assembly>(_opt, Assembly::Base); break; }

                        case MODE_ABUNDANCE:
                        {
                            analyze<Abundance>(optarg, detect<Abundance::Level>(_opt));
                            break;
                        }

                        case MODE_DIFFERENTIAL:
                        {
                            analyze<Differential>(optarg, detect<Differential::Level>(_opt));
                            break;
                        }
                    }
                }
                
                break;
            }

            case CMD_DNA:
            {
                if (_mode != MODE_SEQS       &&
                    _mode != MODE_SEQUENCING &&
                    _mode != MODE_ALIGN      &&
                    _mode != MODE_VARIATION)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_SEQS:      { print_dna();                                   break; }
                        case MODE_ALIGN:     { analyze<DAlign>(optarg, DAlign::Base);         break; }
                        case MODE_VARIATION: { analyze<Structural>(optarg, Structural::Base); break; }
                    }
                }
                
                break;
            }
                
            case CMD_META:
            {
                if (_mode != MODE_SEQS  &&
                    _mode != MODE_ALIGN &&
                    _mode != MODE_DE_NOVO)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_SEQS:    { print_meta();                          break; }
                        case MODE_ALIGN:   { analyze<MAlign>(optarg, MAlign::Base); break; }
                        case MODE_DE_NOVO: { analyze<Denovo>(optarg, Denovo::Base); break; }
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