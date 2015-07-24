#include <map>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <execinfo.h>

#include "trans/t_diffs.hpp"
#include "trans/t_align.hpp"
#include "trans/t_express.hpp"
#include "trans/t_assembly.hpp"

#include "var/v_freq.hpp"
#include "var/v_align.hpp"
#include "var/v_discover.hpp"

#include "meta/m_blast.hpp"
#include "meta/m_diffs.hpp"
#include "meta/m_assembly.hpp"

#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"

#include "fusion/f_discover.hpp"

#include "parsers/parser_csv.hpp"
#include "parsers/parser_sequins.hpp"

#include "writers/file_writer.hpp"
#include "writers/terminal_writer.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

typedef int Tool;
typedef int Input;
typedef int Option;

typedef std::string Value;
typedef std::set<Value> Range;

#define TOOL_VERSION    'v'
#define TOOL_TEST       264
#define TOOL_T_SEQUIN   265
#define TOOL_T_ALIGN    266
#define TOOL_T_ASSEMBLY 267
#define TOOL_T_EXPRESS  268
#define TOOL_T_DIFF     269
#define TOOL_T_NORM     270
#define TOOL_T_IGV      271
#define TOOL_V_ALIGN    272
#define TOOL_V_DISCOVER 273
#define TOOL_V_FREQ     274
#define TOOL_V_DIFF     275
#define TOOL_V_IGV      276
#define TOOL_M_PSL      277
#define TOOL_M_ALIGN    278
#define TOOL_M_ABUND    279
#define TOOL_M_ASSEMBLY 280
#define TOOL_M_DIFF     281
#define TOOL_M_IGV      282
#define TOOL_L_ABUND    283
#define TOOL_L_DIFF     284
#define TOOL_L_IGV      285
#define TOOL_F_DISCOVER 286
#define TOOL_F_EXPRESS  287
#define TOOL_F_IGV      288

/*
 * Options specified in the command line
 */

#define OPT_TOOL    320
#define OPT_MIN     321
#define OPT_MAX     322
#define OPT_LOS     323
#define OPT_PATH    324
#define OPT_MIXTURE 326
#define OPT_FILTER  327
#define OPT_THREAD  328
#define OPT_VERSION 332
#define OPT_TEST    333
#define OPT_R_BED   380
#define OPT_R_GTF   381
#define OPT_U_VCF   384
#define OPT_GTRACK  385
#define OPT_ITRACK  386
#define OPT_U_GTF   387
#define OPT_GDIFF   388
#define OPT_IDIFF   389
#define OPT_BAM_1   392
#define OPT_BAM_2   393
#define OPT_PSL_1   394
#define OPT_PSL_2   395
#define OPT_FA_1    396
#define OPT_FA_2    397
#define OPT_FUS     398
#define OPT_OUT     399
#define OPT_VCF     400

using namespace Anaquin;

// Shared with other modules
std::string __full_command__;

// Shared with other modules
std::string date()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    
    strftime(buffer, 80, "%d-%m-%Y %I:%M:%S", timeinfo);
    std::string str(buffer);
    
    return str;
}

/*
 * Defines the possible tools and representations
 */

static std::map<Value, Tool> _tools =
{
    { "test",             TOOL_TEST       },

    { "TransSequin",      TOOL_T_SEQUIN   },
    { "TransAlign",       TOOL_T_ALIGN    },
    { "TransAssembly",    TOOL_T_ASSEMBLY },
    { "TransExpress",     TOOL_T_EXPRESS  },
    { "TransExpression",  TOOL_T_EXPRESS  },
    { "TransDiff",        TOOL_T_DIFF     },
    { "TransDifferent",   TOOL_T_DIFF     },
    { "TransNorm",        TOOL_T_NORM     },
    { "TransIGV",         TOOL_T_IGV      },

    { "VarAlign",         TOOL_V_ALIGN    },
    { "VarDiscover",      TOOL_V_DISCOVER },
    { "VarFrequency",     TOOL_V_FREQ     },
    { "VarIGV",           TOOL_V_IGV      },

    { "MetaPSL",          TOOL_M_PSL      },
    { "MetaAlign",        TOOL_M_ALIGN    },
    { "MetaAssembly",     TOOL_M_ASSEMBLY },
    { "MetaAbund",        TOOL_M_ABUND    },
    { "MetaAbundance",    TOOL_M_ABUND    },
    { "MetaDiff",         TOOL_M_DIFF     },
    { "MetaDifferent",    TOOL_M_DIFF     },
    { "MetaIGV",          TOOL_V_IGV      },

    { "LadderAbund",      TOOL_L_ABUND    },
    { "LadderAbundance",  TOOL_L_ABUND    },
    { "LadderDiff",       TOOL_L_DIFF     },
    { "LadderDifferent",  TOOL_L_DIFF     },
    { "LadderIGV",        TOOL_L_IGV      },

    { "FusionDiscover",   TOOL_F_DISCOVER },
    { "FusionExpress",    TOOL_F_EXPRESS  },
    { "FusionExpression", TOOL_F_EXPRESS  },
    { "FusionIGV",        TOOL_F_IGV      },
};

/*
 * Defines the tools that require a mixture file
 */

static std::set<Tool> _mixes =
{
    TOOL_T_ASSEMBLY, TOOL_T_EXPRESS, TOOL_T_DIFF, TOOL_T_NORM
};

/*
 * Defines the inputs that each tool expects
 */

static std::map<Tool, std::set<Input>> _inputs =
{
    { TOOL_T_ALIGN, { OPT_R_GTF, OPT_BAM_1 } }
};

/*
 * Variables used in argument parsing
 */

struct Parsing
{
    // The path that output files are written
    std::string path = "output";

    // Context specific reference files
    std::string ref_1, ref_2, ref_3;
    
    // Context specific input files
    std::string input_1, input_2, input_3;

    int name_1, name_2, name_3;
    
    // Context specific input files
    std::map<Input, std::string> inputs;
    
    std::string mix;
    
    // Number of threads
    unsigned threads = 1;
    
    // Custom minmium concentration
    double min = 0;
    
    // Custom maximum concentration
    double max = std::numeric_limits<double>::max();
    
    // Custom sensivitiy
    double los;
    
    // The sequins that have been filtered
    std::set<SequinID> filters;
    
    // How the software is invoked
    std::string command;

    Tool tool = 0;
};

// Wrap the variables so that it'll be easier to reset them
static Parsing _p;

template<typename T> std::string concat(const std::map<Value, T> &m)
{
    std::string str;
    
    for (const auto i : m)
    {
        str += i.first + "|";
    }
    
    return str.substr(0, str.length() - 2);
}

static std::string toolRange()
{
    return "TransAlign,TransAssembly,TransExpression,TransDifferent,TransNorm,TransIGV";
}

struct InvalidCommandException : public std::exception
{
    InvalidCommandException(const std::string &data) : data(data) {}

    // The exact meaning is context-specific
    std::string data;
};

struct InvalidOptionException : public std::exception
{
    InvalidOptionException(const std::string &opt) : opt(opt) {}
    
    std::string opt;
};

/*
 * The type of the argument is invalid, the expected type is integer. For example, giving
 * "ABCD" as the number of threads.
 */

struct InvalidIntegerError : public InvalidCommandException
{
    InvalidIntegerError(const std::string &arg) : InvalidCommandException(arg) {}
};

// An option is being given more than once
struct RepeatOptionError : public InvalidCommandException
{
    RepeatOptionError(const std::string &opt) : InvalidCommandException(opt) {}
};

// A mandatory option is missing, for instance, failing to specify the command
struct MissingOptionError : public std::exception
{
    MissingOptionError(const std::string &opt) : opt(opt) {}
    MissingOptionError(const std::string &opt, const std::string &range) : opt(opt), range(range) {}

    // Option that is missing
    const std::string opt;
    
    // Possible values for the missing option
    const std::string range;
};

struct MissingMixtureError   : public std::exception {};
struct MissingReferenceError : public std::exception {};
struct MissingInputError     : public std::exception {};

struct InvalidToolError : public std::exception
{
    InvalidToolError(const std::string &value) : value(value) {}

    std::string value;
};

struct InvalidInputCountError : std::exception
{
    InvalidInputCountError(std::size_t expected, std::size_t actual) : expected(expected), actual(actual) {}
    
    // Number of inputs detected
    std::size_t actual;
    
    // Number of inputs expected
    std::size_t expected;
};

struct TooLessInputError : std::exception
{
    TooLessInputError(std::size_t n) : n(n) {}
    
    // Number of inputs expected
    std::size_t n;
};

struct TooManyOptionsError : public std::runtime_error
{
    TooManyOptionsError(const std::string &msg) : std::runtime_error(msg) {}
};

struct InvalidFilterError : public std::runtime_error
{
    InvalidFilterError(const std::string &msg) : std::runtime_error(msg) {}
};

/*
 * Argument options
 */

static const char *short_options = "";

static const struct option long_options[] =
{
    { "v", no_argument, 0, OPT_VERSION },

    { "t",    required_argument, 0, OPT_TOOL },
    { "tool", required_argument, 0, OPT_TOOL },

    { "u_sam",    required_argument, 0, OPT_BAM_1 },
    { "u_bam",    required_argument, 0, OPT_BAM_1 },
    { "u_sam_1",  required_argument, 0, OPT_BAM_1 },
    { "u_bam_1",  required_argument, 0, OPT_BAM_1 },
    { "u_sam_2",  required_argument, 0, OPT_BAM_2 },
    { "u_bam_2",  required_argument, 0, OPT_BAM_2 },

    { "r_fus",    required_argument, 0, OPT_FUS },
    { "u_out",    required_argument, 0, OPT_OUT },

    { "r_bed",    required_argument, 0, OPT_R_BED },
    { "r_gtf",    required_argument, 0, OPT_R_GTF },

    { "u_vcf",    required_argument, 0, OPT_VCF    },
    { "u_fa",     required_argument, 0, OPT_FA_1   },
    { "u_fa_1",   required_argument, 0, OPT_FA_1   },
    { "u_fa_2",   required_argument, 0, OPT_FA_2   },
    { "u_gtf",    required_argument, 0, OPT_U_GTF  },
    { "u_gtrack", required_argument, 0, OPT_GTRACK },
    { "u_itrack", required_argument, 0, OPT_ITRACK },
    { "u_gdiff",  required_argument, 0, OPT_GDIFF  },
    { "u_idiff",  required_argument, 0, OPT_IDIFF  },
    { "u_psl",    required_argument, 0, OPT_PSL_1  },
    { "u_psl_1",  required_argument, 0, OPT_PSL_1  },
    { "u_psl_2",  required_argument, 0, OPT_PSL_2  },

    { "min", required_argument, 0, OPT_MIN },
    { "max", required_argument, 0, OPT_MAX },

    { "m",   required_argument, 0, OPT_MIXTURE },
    { "mix", required_argument, 0, OPT_MIXTURE },

    { "o",      required_argument, 0, OPT_PATH  },
    { "output", required_argument, 0, OPT_PATH  },

    { "f",      required_argument, 0, OPT_FILTER },
    { "filter", required_argument, 0, OPT_FILTER },

    {0, 0, 0, 0 }
};

static void printUsage()
{
    extern std::string Manual();
    std::cout << Manual() << std::endl;
}

static void printVersion()
{
    std::cout << "Anaquin v1.1.01" << std::endl;
}

// Print a file of mixture A and B
static void print(Reader &r)
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
        Tokens::split(l, "\t", tokens);

        std::cout << tokens[0] << "\t" << tokens[2] << "\t" << tokens[3] << std::endl;
    }
}

static void printMixture()
{
    Reader r(_p.mix);
    print(r);
}

template <typename Mixture> void applyMix(Mixture mix)
{
    if (_p.mix.empty())
    {
        //if (_mixes.at(_p.mode))
        {
            //throw MissingMixtureError();
        }
    }
    else
    {
        std::cout << "[INFO]: Mixture: " << _p.mix << std::endl;
        mix(Reader(_p.mix));
    }
}

template <typename Reference> void applyRef(Reference ref)
{
    if (_p.ref_1.empty())
    {
        //if (_needRef.at(_p.mode))
        {
            //throw MissingReferenceError();
        }
    }
    else
    {
        std::cout << "[INFO]: Reference: " << _p.ref_1 << std::endl;
        ref(Reader(_p.ref_1));
    }
}

// Read sequins from a file, one per line. The identifiers must match.
static void readFilters(const std::string &file)
{
//    Reader r(file);
//    std::string line;
//    
//    // We'll use it to compare the sequins
//    const auto &s = Standard::instance();
//
//    while (r.nextLine(line))
//    {
//        switch (_p.cmd)
//        {
//            case CMD_FUSION: { break; }
//            case CMD_LADDER: { break; }
//
//            case CMD_RNA:
//            {
//                assert(s.r_seqs_A.size() == s.r_seqs_B.size());
//
//                if (!s.r_seqs_A.count(line))
//                {
//                    throw InvalidFilterError("Unknown sequin for RNA: " + line);
//                }
//                
//                _p.filters.insert(line);
//                break;
//            }
//
//            case CMD_VAR:
//            {
//                assert(s.v_seqs_A.size() == s.v_seqs_B.size());
//                
//                if (!s.v_seqs_A.count(line))
//                {
//                    throw InvalidFilterError("Unknown sequin for DNA: " + line);
//                }
//                
//                _p.filters.insert(line);
//                break;
//            }
//
//            case CMD_META:
//            {
//                assert(s.m_seqs_A.size() == s.m_seqs_B.size());
//
//                if (!s.m_seqs_A.count(line))
//                {
//                    throw InvalidFilterError("Unknown sequin for metagenomics: " + line);
//                }
//
//                _p.filters.insert(line);
//                break;
//            }
//
//            default: { assert(false); }
//        }
//    }
//
//    if (_p.filters.empty())
//    {
//        throw InvalidFilterError("No sequin found in: " + file);
//    }
}

template <typename Analyzer, typename F> void analyzeF(F f, typename Analyzer::Options o)
{
    const auto path = _p.path;

    // This might be needed while scripting
    __full_command__ = _p.command;

#ifndef DEBUG
    o.writer = std::shared_ptr<FileWriter>(new FileWriter(path));
    o.logger = std::shared_ptr<FileWriter>(new FileWriter(path));
    o.output = std::shared_ptr<TerminalWriter>(new TerminalWriter());
    o.logger->open("anaquin.log");
#endif

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    o.info(_p.command);
    o.info(date());
    o.info("Path: " + path);
    o.info("Threads: " + std::to_string(_p.threads) + "\n");

    for (const auto &filter : (o.filters = _p.filters))
    {
        std::cout << "Filter: " << filter << std::endl;
    }

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "------------- Sequin Analysis -----------" << std::endl;
    std::cout << "-----------------------------------------" << std::endl << std::endl;
    
    std::clock_t begin = std::clock();
    
    f(o);
    
    std::clock_t end = std::clock();

    const auto elapsed = (boost::format("Completed. %1% seconds.") % (double(end - begin) / CLOCKS_PER_SEC)).str();
    o.info(elapsed);

#ifndef DEBUG
    o.logger->close();
#endif
}

// Analyze for a single sample
template <typename Analyzer> void analyze_1(typename Analyzer::Options o = typename Analyzer::Options())
{
    if (_p.input_1.empty())
    {
        throw InvalidInputCountError(1, 1);
    }

    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::analyze(_p.input_1, o);
    }, o);
}

// Analyze for a single sample
template <typename Analyzer> void analyze_1(Input x, typename Analyzer::Options o = typename Analyzer::Options())
{
    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::analyze(_p.inputs.at(x), o);
    }, o);
}

// Analyze for two samples
template < typename Analyzer> void analyze_2(Input x, Input y, typename Analyzer::Options o = typename Analyzer::Options())
{
    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::analyze(_p.inputs.at(x), _p.inputs.at(y), o);
    }, o);
}

void parse(int argc, char ** argv)
{
    auto &tool = _p.tool;
    
    _p = Parsing();

    if (argc <= 1)
    {
        printUsage();
    }

    int next, index;

#ifdef UNIT_TESTING
    optind = optreset = 1;
#endif

    /*
     * Reconstruct the overall command
     */
    
    for (int i = 0; i < argc; i++)
    {
        _p.command += std::string(argv[i]) + " ";
    }

    assert(!_p.command.empty());

    // Attempt to parse and store a floating point from string
    auto parseDouble = [&](const std::string &str, double &r)
    {
        assert(next);
        
        try
        {
            r = stof(str);
        }
        catch (...)
        {
            throw;
        }
    };
    
    // Attempt to parse and store an integer from string
    auto parseInt = [&](const std::string &str, unsigned &r)
    {
        assert(next);
        
        try
        {
            r = stoi(str);
        }
        catch (...)
        {
            throw std::runtime_error("eeee");
        }
    };
    
    auto checkFile = [&](const std::string &file)
    {
        if (!std::ifstream(file).good())
        {
            throw InvalidFileError(file);
        }
    };

    /*
     * Pre-process arguments. This way, we can examine the options in whatever order we'd like to
     */

    std::vector<Option> opts;
    std::vector<Value>  vals;

    while ((next = getopt_long_only(argc, argv, short_options, long_options, &index)) != -1)
    {
        opts.push_back(next);

        // Whether this option has an value
        const bool hasValue = optarg;
        
        vals.push_back(hasValue ? std::string(optarg) : "");
    }

    /*
     * Here, we move the command option to the front. Therefore, we also check
     * if we've at least specified the command.
     */
    
    // Find the index for the tool
    auto iter = std::find(opts.begin(), opts.end(), OPT_TOOL);

    if (iter == opts.end() && (iter  = std::find(opts.begin(), opts.end(), OPT_VERSION))  == opts.end())
    {
        throw MissingOptionError("-t", toolRange());
    }

    // This is the index that we'll need to swap
    const auto i = std::distance(opts.begin(), iter);

    std::swap(opts[0], opts[i]);
    std::swap(vals[0], vals[i]);

    /*
     * Now, the first option is also the tool option.
     */

    for (std::size_t i = 0; i < opts.size(); i++)
    {
        const auto opt = opts[i];
        const auto val = vals[i];

        switch (opt)
        {
            case OPT_VERSION:
            {
                _p.tool = opt;
                
                if (argc != 2)
                {
                    throw TooManyOptionsError("Too many options given for -v");
                }

                break;
            }

            case OPT_TOOL:
            {
                if (!_tools.count(val))
                {
                    throw InvalidToolError(val);
                }

                // We'll work with it's internal representation
                _p.tool = _tools.at(val);

                break;
            }
                
            case OPT_MIXTURE: { checkFile(_p.mix = val);   break; }

            case OPT_FUS:
            case OPT_R_BED:
            case OPT_R_GTF:   { checkFile(_p.ref_1 = val); break; }

            /*
             * Options that take a generated input file for the second sample
             */

            case OPT_FA_2:
            case OPT_PSL_2:
            case OPT_BAM_2:
            {
                _p.inputs[opt] = val;
                _p.name_2 = opt;
                checkFile(_p.input_2 = val);
                break;
            }

            case OPT_VCF:
            case OPT_OUT:
            case OPT_FA_1:
            case OPT_PSL_1: { checkFile(_p.inputs[opt] = val); break; }
                
            case OPT_BAM_1:
            case OPT_U_GTF:
            case OPT_IDIFF:
            case OPT_GDIFF:
            case OPT_GTRACK:
            case OPT_ITRACK:
            {
                _p.inputs[opt] = val;
                _p.name_1 = opt;
                checkFile(_p.input_1 = val);
                break;
            }

            case OPT_PATH:    { _p.path = val;             break; }
            case OPT_FILTER:  { readFilters(val);          break; }
            case OPT_MAX:     { parseDouble(val, _p.max);  break; }
            case OPT_MIN:     { parseDouble(val, _p.min);  break; }
            case OPT_LOS:     { parseDouble(val, _p.los);  break; }
            case OPT_THREAD:  { parseInt(val, _p.threads); break; }

            default:
            {
                throw InvalidOptionException(argv[index]);
            }
        }
    }

    // Exception should've already been thrown if tool is not specified
    assert(_p.tool);

    auto &s = Standard::instance();
    
    switch (_p.tool)
    {
        case TOOL_VERSION: { printVersion();                break; }
        case TOOL_TEST:    { Catch::Session().run(1, argv); break; }

        case TOOL_T_IGV:
        case TOOL_T_NORM:
        case TOOL_T_DIFF:
        case TOOL_T_ALIGN:
        case TOOL_T_SEQUIN:
        case TOOL_T_EXPRESS:
        case TOOL_T_ASSEMBLY:
        {
            std::cout << "[INFO]: Transcriptome Analysis" << std::endl;
            
            applyRef(std::bind(&Standard::r_ref, &s, std::placeholders::_1));
            applyMix(std::bind(&Standard::r_mix, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_T_SEQUIN:   { printMixture();         break; }
                case TOOL_T_ALIGN:    { analyze_1<TAlign>();    break; }
                case TOOL_T_ASSEMBLY: { analyze_1<TAssembly>(); break; }

                case TOOL_T_EXPRESS:
                {
                    TExpress::Options o;
                    
                    // The interpreation crticially depends on genes or isoforms
                    o.level = _p.inputs.count(OPT_GTRACK) ? TExpress::Gene : TExpress::Isoform;
                    
                    analyze_1<TExpress>(o);
                    break;
                }

                case TOOL_T_DIFF:
                {
                    TDiffs::Options o;

                    // The interpreation crticially depends on genes or isoforms
                    o.level = _inputs.count(OPT_GTRACK) ? TDiffs::Gene : TDiffs::Isoform;

                    analyze_1<TDiffs>(o);
                    break;
                }
            }
            
            break;
        }

/*
        case CMD_CANCER:
        {
            std::cout << "[INFO]: Cancer Analysis" << std::endl;
            break;
        }
            
        case CMD_CLINIC:
        {
            std::cout << "[INFO]: Clinic Analysis" << std::endl;
            break;
        }
*/
            
        case TOOL_F_IGV:
        case TOOL_F_EXPRESS:
        case TOOL_F_DISCOVER:
        {
            std::cout << "[INFO]: Fusion Analysis" << std::endl;

            applyRef(std::bind(&Standard::f_ref, &s, std::placeholders::_1));
            applyMix(std::bind(&Standard::f_mix, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_F_IGV:      { break; }
                case TOOL_F_EXPRESS:  { break; }
                case TOOL_F_DISCOVER: { analyze_1<FDiscover>(OPT_OUT); break; }
            }

            break;
        }

        case TOOL_L_IGV:
        case TOOL_L_DIFF:
        case TOOL_L_ABUND:
        {
            std::cout << "[INFO]: Ladder Analysis" << std::endl;

            applyMix(std::bind(&Standard::l_mix, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_L_ABUND:  { analyze_1<LAbund>();                     break; }
                case TOOL_L_DIFF:   { analyze_2<LDiffs>(OPT_BAM_1, OPT_BAM_2); break; }
            }

            break;
        }

        case TOOL_V_IGV:
        case TOOL_V_FREQ:
        case TOOL_V_DIFF:
        case TOOL_V_ALIGN:
        case TOOL_V_DISCOVER:
        {
            std::cout << "[INFO]: Variant Analysis" << std::endl;
            
            applyRef(std::bind(&Standard::v_ref, &s, std::placeholders::_1));
            applyMix(std::bind(&Standard::v_mix, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_V_ALIGN:    { analyze_1<VAlign>(OPT_BAM_1);  break; }
                case TOOL_V_DISCOVER: { analyze_1<VDiscover>(OPT_VCF); break; }
                case TOOL_V_FREQ:     { analyze_1<VFreq>();     break; }
                //case TOOL_V_DIFF:     { analyze_1<VVariant>(); break; }
                //case TOOL_V_IGV:      { analyze_1<VVariant>(); break; }
            }

            break;
        }
            
        case TOOL_M_IGV:
        case TOOL_M_PSL:
        case TOOL_M_DIFF:
        case TOOL_M_ALIGN:
        case TOOL_M_ABUND:
        case TOOL_M_ASSEMBLY:
        {
            std::cout << "[INFO]: Metagenomics Analysis" << std::endl;
            
            applyRef(std::bind(&Standard::m_ref, &s, std::placeholders::_1));
            applyMix(std::bind(&Standard::m_mix, &s, std::placeholders::_1));
            
            switch (_p.tool)
            {
                case TOOL_M_PSL: { analyze_1<MBlast>(); break; }

                case TOOL_M_DIFF:
                {
                    MDiffs::Options o;
                    
                    o.pA = _p.inputs.at(OPT_PSL_1);
                    o.pB = _p.inputs.at(OPT_PSL_2);

                    analyze_2<MDiffs>(OPT_FA_1, OPT_FA_2, o);
                    break;
                }

                case TOOL_M_ASSEMBLY:
                {
                    MAssembly::Options o;
                    
                    o.psl = _p.inputs.at(OPT_PSL_1);

                    analyze_1<MAssembly>(o);
                    break;
                }
            }

            break;
        }

        default: { assert(false); }
    }
}

int parse_options(int argc, char ** argv)
{
    auto printError = [&](const std::string &file)
    {
        std::cerr << std::endl;
        std::cerr << "*********************************************************************" << std::endl;
        std::cerr << file << std::endl;
        std::cerr << "*********************************************************************" << std::endl;
    };
    
    try
    {
        parse(argc, argv);
        return 0;
    }
    catch (const InvalidToolError &ex)
    {
        printError("Invalid tool: " + ex.value);
    }
    catch (const InvalidOptionException &ex)
    {
        const auto format = "Unknown option: %1%";
        printError((boost::format(format) % ex.opt).str());
    }
    catch (const MissingOptionError &ex)
    {
        if (!ex.range.empty())
        {
            const auto format = "A mandatory option is missing. Please specify %1%. Possibilities are %2%";
            printError((boost::format(format) % ex.opt % ex.range).str());
        }
        else
        {
            const auto format = "A mandatory option is missing. Please specify %1%.";
            printError((boost::format(format) % ex.opt).str());
        }
    }
    catch (const RepeatOptionError &ex)
    {
        printError((boost::format("The option %1% has been repeated. Please check and try again.") % ex.data).str());
    }
    catch (const MissingInputError &ex)
    {
        printError("No input file given. Please give an input file and try again.");
    }
    catch (const MissingMixtureError &ex)
    {
        printError("Mixture file is missing. Please specify it with -m.");
    }
    catch (const MissingReferenceError &ex)
    {
        printError("Reference file is missing. Please specify it with -r.");
    }
    catch (const InvalidFileError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid file: " % ex.file).str());
    }
    catch (const InvalidFilterError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid filter: " % ex.what()).str());
    }
    catch (const std::runtime_error &ex)
    {
        printError(ex.what());
    }

    return 1;
}

int main(int argc, char ** argv)
{
    return parse_options(argc, argv);
}