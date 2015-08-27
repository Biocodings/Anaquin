#include <map>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <execinfo.h>

#include "trans/t_diffs.hpp"
#include "trans/t_align.hpp"
#include "trans/t_viewer.hpp"
#include "trans/t_express.hpp"
#include "trans/t_assembly.hpp"

#include "var/v_align.hpp"
#include "var/v_viewer.hpp"
#include "var/v_discover.hpp"

#include "meta/m_blast.hpp"
#include "meta/m_diffs.hpp"
#include "meta/m_assembly.hpp"

#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"

#include "fusion/f_viewer.hpp"
#include "fusion/f_express.hpp"
#include "fusion/f_discover.hpp"

#include "parsers/parser_csv.hpp"
#include "parsers/parser_sequins.hpp"

#include "writers/file_writer.hpp"
#include "writers/terminal_writer.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

typedef int Tool;
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

#define OPT_TEST    333
#define OPT_TOOL    320
#define OPT_MIN     321
#define OPT_MAX     322
#define OPT_LOS     323
#define OPT_PATH    324
#define OPT_FILTER  327
#define OPT_THREAD  328
#define OPT_VERSION 332

#define OPT_R_BASE  800
#define OPT_R_BED   801
#define OPT_R_GTF   802
#define OPT_R_FUS   803
#define OPT_R_VCF   804
#define OPT_MIXTURE 805

#define OPT_U_BASE  900
#define OPT_U_VCF   901
#define OPT_U_GTF   902
#define OPT_GTRACK  903
#define OPT_ITRACK  904
#define OPT_GDIFF   905
#define OPT_IDIFF   906
#define OPT_BAM_1   907
#define OPT_BAM_2   908
#define OPT_PSL_1   909
#define OPT_PSL_2   910
#define OPT_FA_1    911
#define OPT_FA_2    912
#define OPT_U_OUT   913

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
 * Defines the options that are expected
 */

static std::map<Tool, std::set<Option>> _required =
{
    { TOOL_T_ALIGN,    { OPT_R_GTF, OPT_MIXTURE, OPT_BAM_1 } },
    { TOOL_T_ASSEMBLY, { OPT_R_GTF, OPT_MIXTURE, OPT_U_GTF } },
    { TOOL_T_EXPRESS,  { OPT_R_GTF, OPT_MIXTURE } },
    { TOOL_T_DIFF,     { OPT_R_GTF, OPT_MIXTURE } },
};

/*
 * Defines the options that one of the possibilites must be defined
 */

static std::map<Tool, std::set<Option>> _pick =
{
    { TOOL_T_EXPRESS, { OPT_GTRACK, OPT_ITRACK } },
    { TOOL_T_DIFF,    { OPT_GDIFF,  OPT_IDIFF  } },
};

/*
 * Variables used in argument parsing
 */

struct Parsing
{
    // The path that output files are written
    std::string path = "output";

    // Context specific options
    std::map<Option, std::string> opts;
    
    // Number of threads
    unsigned threads = 1;
    
    // Minmium concentration
    double min = 0;
    
    // Maximum concentration
    double max = std::numeric_limits<double>::max();

    // Sensivitiy
    double los;
    
    // The sequins that have been filtered
    std::set<SequinID> filters;
    
    // How the software is invoked
    std::string command;

    Tool tool = 0;
};

// Wrap the variables so that it'll be easier to reset them
static Parsing _p;

template<typename T> std::string concat(const std::set<Value, T> &m)
{
    std::string str;
    
    for (const auto i : m)
    {
        str += i.first + "|";
    }
    
    return str.substr(0, str.length() - 2);
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
struct InvalidNegativ30Error : public std::exception {};

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

    { "usam",    required_argument, 0, OPT_BAM_1 },
    { "ubam",    required_argument, 0, OPT_BAM_1 },
    { "usam1",   required_argument, 0, OPT_BAM_1 },
    { "ubam1",   required_argument, 0, OPT_BAM_1 },
    { "usam2",   required_argument, 0, OPT_BAM_2 },
    { "ubam2",   required_argument, 0, OPT_BAM_2 },

    { "rfus",    required_argument, 0, OPT_R_FUS },
    { "uout",    required_argument, 0, OPT_U_OUT },

    { "rbed",    required_argument, 0, OPT_R_BED },
    { "rgtf",    required_argument, 0, OPT_R_GTF },

    { "uvcf",    required_argument, 0, OPT_U_VCF  },
    { "ufa",     required_argument, 0, OPT_FA_1   },
    { "ufa1",    required_argument, 0, OPT_FA_1   },
    { "ufa2",    required_argument, 0, OPT_FA_2   },
    { "ugtf",    required_argument, 0, OPT_U_GTF  },
    { "ugtrack", required_argument, 0, OPT_GTRACK },
    { "uitrack", required_argument, 0, OPT_ITRACK },
    { "ugdiff",  required_argument, 0, OPT_GDIFF  },
    { "uidiff",  required_argument, 0, OPT_IDIFF  },
    { "upsl",    required_argument, 0, OPT_PSL_1  },
    { "upsl1",   required_argument, 0, OPT_PSL_1  },
    { "upsl2",   required_argument, 0, OPT_PSL_2  },

    { "min", required_argument, 0, OPT_MIN },
    { "max", required_argument, 0, OPT_MAX },

    { "m",   required_argument, 0, OPT_MIXTURE },
    { "mix", required_argument, 0, OPT_MIXTURE },

    { "o",      required_argument, 0, OPT_PATH },
    { "output", required_argument, 0, OPT_PATH },

    { "f",      required_argument, 0, OPT_FILTER },
    { "filter", required_argument, 0, OPT_FILTER },

    {0, 0, 0, 0 }
};

static std::string optToStr(int opt)
{
    for (const auto o : long_options)
    {
        if (o.val == opt)
        {
            return o.name;
        }
    }
    
    throw std::runtime_error("Invalid option: " + std::to_string(opt));
}

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

static std::string mixture()
{
    return _p.opts[OPT_MIXTURE];
}

static void printMixture()
{
    Reader r(mixture());
    print(r);
}

static void printError(const std::string &msg)
{
    std::cerr << std::endl;
    std::cerr << "*********************************************************************" << std::endl;
    std::cout << msg << std::endl;
    std::cerr << "*********************************************************************" << std::endl;
}

template <typename Mixture> void applyMix(Mixture mix)
{
    if (mixture().empty())
    {
        return;
    }
    
    std::cout << "[INFO]: Mixture: " << mixture() << std::endl;
    mix(Reader(mixture()));
}

#define CHECK_REF(x) (x != OPT_MIXTURE && x > OPT_R_BASE && x < OPT_U_BASE)

/*
 * Apply reference resource assuming there is only a single reference source
 */

template <typename Reference> void applyRef(Reference ref)
{
    for (const auto &i : _p.opts)
    {
        const auto opt = i.first;
        
        if (CHECK_REF(opt))
        {
            std::cout << "[INFO]: Reference: " << _p.opts[opt] << std::endl;
            ref(Reader(_p.opts[opt]));
        }
    }
}

/*
 * Apply two reference resources for two options
 */

template <typename Reference> void applyRef(Reference ref, Option o1, Option o2)
{
    assert(o1 != o2 && CHECK_REF(o1) && CHECK_REF(o2));
    
    std::cout << "[INFO]: Reference: " << _p.opts[o1] << std::endl;
    std::cout << "[INFO]: Reference: " << _p.opts[o2] << std::endl;

    ref(Reader(_p.opts[o1]), Reader(_p.opts[o2]));
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
    
    std::clock_t begin = std::clock();
    
    f(o);
    
    std::clock_t end = std::clock();

    const auto elapsed = (boost::format("Completed. %1% seconds.") % (double(end - begin) / CLOCKS_PER_SEC)).str();
    o.info(elapsed);

#ifndef DEBUG
    o.logger->close();
#endif
}

template <typename Viewer> void viewer(typename Viewer::Options o = typename Viewer::Options())
{
    // Where the session files are generated
    o.path = _p.path;

    Viewer::generate(_p.opts.at(OPT_BAM_1), o);
}

// Analyze for a single sample
template <typename Analyzer> void analyze_1(Option x, typename Analyzer::Options o = typename Analyzer::Options())
{
    Standard::instance().r_meta.validate();

    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::analyze(_p.opts.at(x), o);
    }, o);
}

// Analyze for two samples
template < typename Analyzer> void analyze_2(Option x, Option y, typename Analyzer::Options o = typename Analyzer::Options())
{
    Standard::instance().r_meta.validate();

    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::analyze(_p.opts.at(x), _p.opts.at(y), o);
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

    unsigned n = 0;

    /*
     * Detect for weird inputs, such as "–" (invalid ASCII)
     */

    for (auto i = 0; i < argc; i++)
    {
        if (argv[i])
        {
            const auto str = std::string(argv[i]);
            
            if (str.size() >= 2)
            {
                const int key = (int) str[0];
                
                if (key == -30) // –
                {
                    throw std::runtime_error("Invalid " + std::string(argv[i]) + ". Please note '–' is NOT the character '-' that you see on your keyboard. The given option is therefore invalid, please type the character manually.");
                }
            }
        }
    }

    while ((next = getopt_long_only(argc, argv, short_options, long_options, &index)) != -1)
    {
        if (next < OPT_TOOL)
        {
            throw InvalidOptionException(argv[n+1]);
        }
        
        opts.push_back(next);

        // Whether this option has an value
        const bool hasValue = optarg;
        
        n += hasValue ? 2 : 1;
        
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
        throw MissingOptionError("-t");
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
                _p.tool = TOOL_VERSION;

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

            case OPT_FA_1:
            case OPT_FA_2:
            case OPT_R_FUS:
            case OPT_U_VCF:
            case OPT_U_OUT:
            case OPT_BAM_1:
            case OPT_U_GTF:
            case OPT_IDIFF:
            case OPT_GDIFF:
            case OPT_PSL_2:
            case OPT_BAM_2:
            case OPT_R_BED:
            case OPT_R_GTF:
            case OPT_PSL_1:
            case OPT_GTRACK:
            case OPT_ITRACK:
            case OPT_MIXTURE: { checkFile(_p.opts[opt] = val); break; }

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
    
    /*
     * Have all the required options given?
     */
    
    if (_required.count(_p.tool))
    {
        auto required = _required[_p.tool];
        
        for (const auto i : _p.opts)
        {
            if (required.count(i.first))
            {
                required.erase(i.first);
            }
        }

        if (!required.empty())
        {
            throw MissingOptionError(optToStr(*required.begin()));
        }
    }

    /*
     * Have all the pick options given?
     */
    
    if (_pick.count(_p.tool))
    {
        for (const auto i : _p.opts)
        {
            if (_pick[_p.tool].count(i.first))
            {
                goto mainSwitch;
            }
        }
        
        throw MissingInputError();
    }

    mainSwitch:

    if (_p.tool != TOOL_VERSION)
    {
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "------------- Sequin Analysis -----------" << std::endl;
        std::cout << "-----------------------------------------" << std::endl << std::endl;        
    }

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
                case TOOL_T_SEQUIN: { printMixture();               break; }
                case TOOL_T_ALIGN:  { analyze_1<TAlign>(OPT_BAM_1); break; }
                
                case TOOL_T_ASSEMBLY:
                {
                    TAssembly::Options o;

                    o.ref   = _p.opts[OPT_R_GTF];
                    o.query = _p.opts[OPT_R_GTF];
                    
                    analyze_1<TAssembly>(OPT_U_GTF, o);
                    break;
                }

                case TOOL_T_EXPRESS:
                {
                    TExpress::Options o;
                    
                    if (_p.opts.count(OPT_GTRACK))
                    {
                        o.level = TExpress::Isoform;
                        analyze_1<TExpress>(OPT_GTRACK);
                    }
                    else
                    {
                        o.level = TExpress::Isoform;
                        analyze_1<TExpress>(OPT_ITRACK);
                    }

                    break;
                }

                case TOOL_T_DIFF:
                {
                    TDiffs::Options o;

                    if (_p.opts.count(OPT_GDIFF))
                    {
                        std::cout << "[INFO]: Gene Analysis" << std::endl;
                        
                        o.level = TDiffs::Gene;
                        analyze_1<TDiffs>(OPT_GDIFF, o);
                    }
                    else
                    {
                        std::cout << "[INFO]: Isoform Analysis" << std::endl;

                        o.level = TDiffs::Isoform;
                        analyze_1<TDiffs>(OPT_IDIFF, o);
                    }

                    break;
                }

                case TOOL_T_IGV: { viewer<TViewer>(); break; }
            }

            break;
        }

        case TOOL_F_IGV:
        case TOOL_F_EXPRESS:
        case TOOL_F_DISCOVER:
        {
            std::cout << "[INFO]: Fusion Analysis" << std::endl;

            applyRef(std::bind(&Standard::f_ref, &s, std::placeholders::_1));
            applyMix(std::bind(&Standard::f_mix, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_F_IGV:      { analyze_1<FViewer>(OPT_U_OUT);   break; }
                case TOOL_F_EXPRESS:  { analyze_1<FExpress>(OPT_GTRACK); break; }
                case TOOL_F_DISCOVER: { analyze_1<FDiscover>(OPT_U_OUT); break; }
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
                case TOOL_L_ABUND:  { analyze_1<LAbund>(OPT_BAM_1);            break; }
                case TOOL_L_DIFF:   { analyze_2<LDiffs>(OPT_BAM_1, OPT_BAM_2); break; }
            }

            break;
        }

        case TOOL_V_IGV:
        case TOOL_V_DIFF:
        case TOOL_V_ALIGN:
        case TOOL_V_DISCOVER:
        {
            std::cout << "[INFO]: Variant Analysis" << std::endl;

            switch (_p.tool)
            {
                case TOOL_V_ALIGN:
                {
                    applyRef(std::bind(&Standard::v_std, &s, std::placeholders::_1)); break;
                }

                case TOOL_V_DIFF:
                case TOOL_V_DISCOVER:
                {
                    applyRef(std::bind(&Standard::v_var, &s, std::placeholders::_1)); break;
                }
            }

            applyMix(std::bind(&Standard::v_mix, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_V_ALIGN:    { analyze_1<VAlign>(OPT_BAM_1);    break; }
                case TOOL_V_DISCOVER: { analyze_1<VDiscover>(OPT_U_VCF); break; }
                //case TOOL_V_DIFF:     { analyze_1<VVariant>(); break; }
                case TOOL_V_IGV:      { analyze_1<VViewer>(OPT_PATH);    break; }
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
            
            applyMix(std::bind(&Standard::m_mix_1, &s, std::placeholders::_1));

            switch (_p.tool)
            {
                case TOOL_M_PSL: { analyze_1<MBlast>(OPT_PSL_1); break; }

                case TOOL_M_DIFF:
                {
                    applyMix(std::bind(&Standard::m_mix_2, &s, std::placeholders::_1));
                    
                    MDiffs::Options o;
                    
                    o.pA = _p.opts.at(OPT_PSL_1);
                    o.pB = _p.opts.at(OPT_PSL_2);

                    analyze_2<MDiffs>(OPT_FA_1, OPT_FA_2, o);
                    break;
                }

                case TOOL_M_ASSEMBLY:
                {
                    MAssembly::Options o;
                    
                    // An alignment file is needed to identify contigs
                    o.psl = _p.opts.at(OPT_PSL_1);

                    analyze_1<MAssembly>(OPT_FA_1, o);
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
    try
    {
        parse(argc, argv);
        return 0;
    }
    catch (const InvalidToolError &ex)
    {
        printError("Invalid tool: " + ex.value + ". Please refer to the documenation for correct usage.");
    }
    catch (const InvalidOptionException &ex)
    {
        const auto format = "Unknown option: %1%";
        printError((boost::format(format) % ex.opt).str());
    }
    catch (const MissingOptionError &ex)
    {
        const auto format = "A mandatory option is missing. Please specify -%1%.";
        printError((boost::format(format) % ex.opt).str());
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
    catch (const std::exception &ex)
    {
        printError(ex.what());
    }

    return 1;
}

int main(int argc, char ** argv)
{
    return parse_options(argc, argv);
}