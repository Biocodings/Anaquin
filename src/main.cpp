#include <ctime>
#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/resources.hpp"

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

#define CMD_VER    'v'
#define CMD_TEST   't'
#define CMD_RNA   265
#define CMD_DNA   266
#define CMD_META  267

#define MODE_FILTER       280
#define MODE_SEQUENCE     282
#define MODE_ALIGN        283
#define MODE_ASSEMBLY     284
#define MODE_ABUNDANCE    285
#define MODE_DIFFERENTIAL 286
#define MODE_VARIATION    287
#define MODE_SEQUINS      288

#define OPT_MIN           321
#define OPT_MAX           322
#define OPT_LIMIT         323
#define OPT_OUTPUT        324
#define OPT_PSL           325

using namespace Spike;

/*
 * Variables used in argument parsing
 */

// An alignment file by BLAST used for metagenomics
static std::string _psl;

// The path that output files are written
static std::string _output;

// The sequins that have been restricted
static std::vector<SequinID> _filters;

static int _cmd  = 0;
static int _mode = 0;

// The operands for the command
std::vector<std::string> _opts;

/*
 * Argument options
 */

static const char *short_options = "";

static const struct option long_options[] =
{
    { "t",    no_argument,       0, CMD_TEST },
    { "v",    no_argument,       0, CMD_VER  },

    { "o",     required_argument, 0, OPT_OUTPUT },
    { "min",   required_argument, 0, OPT_MIN    },
    { "max",   required_argument, 0, OPT_MAX    },
    { "limit", required_argument, 0, OPT_LIMIT  },

    { "rna",  required_argument, 0, CMD_RNA  },
    { "dna",  required_argument, 0, CMD_DNA  },
    { "meta", required_argument, 0, CMD_META },

    { "l",       no_argument, 0, MODE_SEQUINS },
    { "sequins", no_argument, 0, MODE_SEQUINS },

    { "r",        required_argument, 0, MODE_FILTER },
    { "restrict", required_argument, 0, MODE_FILTER },
    
    { "al",   required_argument, 0, MODE_ALIGN        },
    { "as",   required_argument, 0, MODE_ASSEMBLY     },
    { "ab",   required_argument, 0, MODE_ABUNDANCE    },
    { "df",   required_argument, 0, MODE_DIFFERENTIAL },
    { "va",   required_argument, 0, MODE_VARIATION    },
    { "de",   required_argument, 0, MODE_ASSEMBLY     },

    { "psl",  required_argument, 0, OPT_PSL },
    
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

// Print a file of mixture A and B
void print(Reader &r)
{
    /*
     * Format: <ID, Mix A, Mix B>
     */

    std::string l;
    
    // Skip the first line
    r.nextLine(l);

    std::cout << "ID\tMix A\tMix B" << std::endl;
    
    while (r.nextLine(l))
    {
        if (l == "\r" || l == "\n" || l == "\r\n")
        {
            continue;
        }

        std::vector<std::string> tokens;
        Tokens::split(l, ",", tokens);

        if (tokens.size() != 3)
        {
            std::cerr << "Malformed file: [" << l << "]" << std::endl;
        }

        std::cout << tokens[0] << "\t" << tokens[1] << "\t" << tokens[2] << std::endl;
    }
}

void print_meta()
{
    extern std::string m_mix_f();
    Reader r(m_mix_f(), String);
    print(r);
}

static void print_rna()
{
    // Empty Implementation
}

void print_dna()
{
    // Empty Implementation
}

// Read sequins from a file, one per line. The identifiers must match.
static bool readFilters(const std::string &file)
{
    Reader r(file);
    std::string line;
    
    // We'll use it to compare the sequins
    const auto &s = Standard::instance();

    while (r.nextLine(line))
    {
        switch (_cmd)
        {
            case CMD_RNA:
            {
                break;
            }

            case CMD_DNA:
            {
                break;
            }

            case CMD_META:
            {
                if (!s.m_seq_A.count(line))
                {
                    std::cerr << "Unknown sequin for metagenomics: " << line << std::endl;
                    return false;
                }

                _filters.push_back(line);
                break;
            }

            default: { print_usage(); break; }
        }
    }

    if (_filters.empty())
    {
        std::cerr << "No sequin found in: " << line << std::endl;
        return false;
    }
    
    return true;
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

int parse_options(int argc, char ** argv)
{
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
            case OPT_PSL:    { _psl = optarg;    break; }
            case OPT_OUTPUT: { _output = optarg; break; }

            case MODE_FILTER:
            {
                if (!readFilters(optarg))
                {
                    return 1;
                }

                break;
            }

            case OPT_MAX:
            case OPT_MIN:
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
                std::cout << "RNA command detected" << std::endl;
                
                if (_mode != MODE_SEQUINS   &&
                    _mode != MODE_SEQUENCE  &&
                    _mode != MODE_ALIGN     &&
                    _mode != MODE_ASSEMBLY  &&
                    _mode != MODE_ABUNDANCE &&
                    _mode != MODE_DIFFERENTIAL)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_SEQUENCE: { break; }
                        case MODE_SEQUINS:  { print_rna();                  break; }
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
                std::cout << "DNA command detected" << std::endl;
                
                if (_mode != MODE_SEQUINS  &&
                    _mode != MODE_SEQUENCE &&
                    _mode != MODE_ALIGN    &&
                    _mode != MODE_VARIATION)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_SEQUENCE:  { break; }
                        case MODE_SEQUINS:   { print_dna();                 break; }
                        case MODE_ALIGN:     { analyze<DAlign>(_opts[0]);   break; }
                        case MODE_VARIATION: { analyze<DVariant>(_opts[0]); break; }
                    }
                }

                break;
            }

            case CMD_META:
            {
                std::cout << "Metagenomics command detected" << std::endl;

                if (_mode != MODE_SEQUINS  &&
                    _mode != MODE_SEQUENCE &&
                    _mode != MODE_ASSEMBLY)
                {
                    print_usage();
                }
                else
                {
                    switch (_mode)
                    {
                        case MODE_SEQUENCE: { break; }
                        case MODE_SEQUINS:  { print_meta(); break; }

                        case MODE_ASSEMBLY:
                        {
                            MAssembly::Options o;

                            // We'd also take an alignment file from a user
                            o.psl = _psl;
                            
                            analyze<MAssembly>(_opts[0], o);
                            break;
                        }
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