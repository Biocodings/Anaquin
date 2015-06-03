#include <ctime>
#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "resources.hpp"
#include "rna/r_align.hpp"
#include "rna/r_assembly.hpp"
#include "rna/r_abundance.hpp"
#include "rna/r_differential.hpp"

#include "dna/d_align.hpp"
#include "dna/d_variant.hpp"
#include "dna/d_sequence.hpp"

#include "meta/m_assembly.hpp"

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
#define MODE_SEQUINS      288

#define OPT_THREAD        320
#define OPT_MIN           321
#define OPT_MAX           322
#define OPT_LIMIT         323
#define OPT_OUTPUT        324

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

// The operands for the command
std::vector<std::string> _opts;

/*
 * Argument options
 */

static const char *short_options = "o";

static const struct option long_options[] =
{
    { "t",    no_argument,       0, CMD_TEST },
    { "v",    no_argument,       0, CMD_VER  },

    { "o",     required_argument, 0, OPT_OUTPUT },
    { "p",     required_argument, 0, OPT_THREAD },
    { "min",   required_argument, 0, OPT_MIN    },
    { "max",   required_argument, 0, OPT_MAX    },
    { "limit", required_argument, 0, OPT_LIMIT  },

    { "rna",  required_argument, 0, CMD_RNA  },
    { "dna",  required_argument, 0, CMD_DNA  },
    { "meta", required_argument, 0, CMD_META },

    { "l",    no_argument,       0, MODE_SEQS         },
    { "seqs", no_argument,       0, MODE_SEQUINS      },
    { "f",    required_argument, 0, MODE_RESTRICTS    },
    { "al",   required_argument, 0, MODE_ALIGN        },
    { "as",   required_argument, 0, MODE_ASSEMBLY     },
    { "ab",   required_argument, 0, MODE_ABUNDANCE    },
    { "df",   required_argument, 0, MODE_DIFFERENTIAL },
    { "va",   required_argument, 0, MODE_VARIATION    },
    { "de",   required_argument, 0, MODE_ASSEMBLY     },

    { "sequins",   no_argument,       0, MODE_SEQUINS      },
    { "align",     required_argument, 0, MODE_ALIGN        },
    { "variant",   required_argument, 0, MODE_VARIATION    },
    { "assembly",  required_argument, 0, MODE_ASSEMBLY     },
    { "abundance", required_argument, 0, MODE_ABUNDANCE    },
    { "diffs",     required_argument, 0, MODE_DIFFERENTIAL },

    {0,  0, 0, 0 }
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
    extern std::string manual();
    std::cout << manual() << std::endl;
}

static void print_version()
{
    const auto c = Resources::chromo();
    const auto m = Resources::mixture();

    std::cout << "Version 1.0. Garvan Institute of Medical Research, 2015." << std::endl;
    std::cout << std::endl;
    std::cout << "Chromosome: " << c.id << " version " << c.v << std::endl;
    std::cout << "Mixture A: version " << m.va << std::endl;
    std::cout << "Mixture B: version " << m.vb << std::endl;
}

static void print_rna_sequins()
{
    for (const auto &i : Standard::instance().r_sequins)
    {
        std::cout << i.id << std::endl;
    }
}

static void print_rna()
{
    const auto &r = Standard::instance();
    const std::string format = "%1%  %2%  %3%  %4%";

    auto f = [&](const Standard::PairMap &seqs)
    {
        std::cout << (boost::format(format) % "r"
                                            % "v"
                                            % "r_con"
                                            % "v_con").str() << std::endl;

        for (std::size_t i = A; i <= D; i++)
        {
            for (auto j : seqs)
            {
                const auto &p = j.second;
                
                if (p.grp == i)
                {
                    std::cout << (boost::format(format) % p.r.id
                                                        % p.v.id
                                                        % p.r.raw
                                                        % p.v.raw).str() << std::endl;
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
    const auto &s = Standard::instance();
    const std::string format = "%1%  %2%  %3%  %4%";
    
    auto f = [&](const Standard::PairMap &seqs)
    {
        std::cout << (boost::format(format) % "r"
                                            % "v"
                                            % "r_con"
                                            % "v_con").str() << std::endl;

        for (std::size_t i = A; i <= D; i++)
        {
            for (auto j : seqs)
            {
                const auto &p = j.second;
                
                if (p.grp == i)
                {
                    std::cout << (boost::format(format) % p.r.id
                                  % p.v.id
                                  % p.r.raw
                                  % p.v.raw).str() << std::endl;
                }
            }
        }
    };

    std::cout << "\n----------------------------- A -----------------------------" << std::endl;
    f(s.d_pair_A);
    std::cout << "\n----------------------------- B -----------------------------" << std::endl;
    f(s.d_pair_B);
}

void print_meta()
{
    // Empty Implementation
}

template <typename Analyzer> void analyze(const std::string &file, typename Analyzer::Options o = typename Analyzer::Options())
{
    o.writer = std::shared_ptr<PathWriter>(new PathWriter(_output.empty() ? "spike_out" : _output));

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "------------- Sequin Analysis -----------" << std::endl;
    std::cout << "-----------------------------------------" << std::endl << std::endl;

    std::clock_t begin = std::clock();
    
    Analyzer::analyze(file, o);

    std::clock_t end = std::clock();
    const double elapsed = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Completed. Elpased: " << elapsed << " seconds" << std::endl;
}

template <typename Options> static Options detect(const std::string &file)
{
    const bool found_gene = file.find("gene") != std::string::npos;
    const bool found_isoform = file.find("isoform") != std::string::npos;

    Options o;
    
    if (found_gene && !found_isoform)
    {
        std::cout << "Calculating for the genes" << std::endl;
        o.level = RNALevel::Gene;
    }
    else if (!found_gene && found_isoform)
    {
        std::cout << "Calcualting for the isoforms" << std::endl;
        o.level = RNALevel::Isoform;
    }
    else
    {
        throw std::runtime_error("Unknown type. Have you specified the level?");
    }

    return o;
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

        if (tmp == "rna")  { _cmd = CMD_RNA;  }
        if (tmp == "dna")  { _cmd = CMD_DNA;  }
        if (tmp == "meta") { _cmd = CMD_META; }
    }
    
    int next, index;

    while ((next = getopt_long_only(argc, argv, short_options, long_options, &index)) != -1)
    {
        switch (next)
        {
            case OPT_OUTPUT: { _output = optarg; break; }

            case OPT_MAX:
            case OPT_MIN:
            case OPT_THREAD:
            case OPT_LIMIT:
            {
                break;
            }

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
                if (_mode != 0)
                {
                    return invalid_cmd();
                }
                
                _mode = next;
                
                if (optarg)
                {
                    _opts.push_back(optarg);
                }

                break;
            }
        }
    }

    if (_cmd == 0)
    {
        print_usage();
    }
    else if ((_cmd == CMD_TEST || _cmd == CMD_VER) && (!_output.empty() || _mode != 0 || !_opts.empty()))
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
                if (_mode != MODE_SEQUINS    &&
                    _mode != MODE_SEQS       &&
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
                        case MODE_SEQUINS:  { print_rna_sequins();          break; }
                        case MODE_SEQS:     { print_rna();                  break; }
                        case MODE_ALIGN:    { analyze<RAlign>(_opts[0]);    break; }
                        case MODE_ASSEMBLY: { analyze<RAssembly>(_opts[0]); break; }

                        case MODE_ABUNDANCE:
                        {
                            analyze<RAbundance>(_opts[0], detect<RAbundance::Options>(_opts[0]));
                            break;
                        }

                        case MODE_DIFFERENTIAL:
                        {
                            analyze<RDifferential>(_opts[0], detect<RDifferential::Options>(_opts[0]));
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
                        case MODE_SEQS:      { print_dna();                break; }
                        case MODE_ALIGN:     { analyze<DAlign>(_opts[0]);   break; }
                        case MODE_VARIATION: { analyze<DVariant>(_opts[0]); break; }
                    }
                }

                break;
            }

            case CMD_META:
            {
                if (_mode != MODE_SEQS && _mode != MODE_ASSEMBLY)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_ASSEMBLY: { analyze<MAssembly>(_opts[0]); break; }
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