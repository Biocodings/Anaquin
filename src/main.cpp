#include <map>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <strings.h>
#include <execinfo.h>

#include "data/experiment.hpp"

#include "TransQuin/t_diffs.hpp"
#include "TransQuin/t_align.hpp"
#include "TransQuin/t_viewer.hpp"
#include "TransQuin/t_express.hpp"
#include "TransQuin/t_assembly.hpp"
#include "TransQuin/t_coverage.hpp"

#include "VarQuin/v_align.hpp"
#include "VarQuin/v_allele.hpp"
#include "VarQuin/v_viewer.hpp"
#include "VarQuin/v_sample.hpp"
#include "VarQuin/v_express.hpp"
#include "VarQuin/v_discover.hpp"
#include "VarQuin/v_coverage.hpp"
#include "VarQuin/v_kexpress.hpp"

#include "MetaQuin/m_blat.hpp"
#include "MetaQuin/m_diffs.hpp"
#include "MetaQuin/m_abund.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "MetaQuin/m_coverage.hpp"

#include "LadQuin/l_copy.hpp"
#include "LadQuin/l_diffs.hpp"
#include "LadQuin/l_abund.hpp"
#include "LadQuin/l_coverage.hpp"

#include "FusQuin/f_diff.hpp"
#include "FusQuin/f_viewer.hpp"
#include "FusQuin/f_express.hpp"
#include "FusQuin/f_discover.hpp"
#include "FusQuin/f_coverage.hpp"

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

#define TOOL_VERSION     'v'
#define TOOL_TEST        264
#define TOOL_T_SEQUIN    265
#define TOOL_T_ALIGN     266
#define TOOL_T_ASSEMBLY  267
#define TOOL_T_EXPRESS   268
#define TOOL_T_COUNT     269
#define TOOL_T_DIFF      270
#define TOOL_T_NORM      271
#define TOOL_T_IGV       272
#define TOOL_T_COVERAGE  273
#define TOOL_V_ALIGN     274
#define TOOL_V_DISCOVER  275
#define TOOL_V_IGV       277
#define TOOL_V_ALLELE    278
#define TOOL_V_COVERAGE  279
#define TOOL_V_SUBSAMPLE 280
#define TOOL_M_ABUND     282
#define TOOL_M_ASSEMBLY  283
#define TOOL_M_DIFF      284
#define TOOL_M_IGV       285
#define TOOL_M_COVERAGE  286
#define TOOL_L_ABUND     287
#define TOOL_L_DIFF      288
#define TOOL_L_COVERAGE  289
#define TOOL_F_DISCOVER  290
#define TOOL_F_EXPRESS   291
#define TOOL_F_IGV       292
#define TOOL_F_COVERAGE  293
#define TOOL_F_DIFF      295
#define TOOL_L_COPY      296
#define TOOL_V_EXPRESS   297
#define TOOL_V_KEXPRESS  298
#define TOOL_T_KEXPRESS  299
#define TOOL_V_REPORT    300
#define TOOL_M_REPORT    301
#define TOOL_L_REPORT    302
#define TOOL_F_REPORT    303
#define TOOL_T_KDIFF     305
#define TOOL_T_REPORT    306

/*
 * Options specified in the command line
 */

#define OPT_TEST     320
#define OPT_TOOL     321
#define OPT_MIN      322
#define OPT_MAX      323
#define OPT_LOS      324
#define OPT_PATH     325
#define OPT_FILTER   326
#define OPT_VERSION  338
#define OPT_SOFT     339
#define OPT_C_SOFT   340

/*
 * References - OPT_R_BASE to OPT_U_BASE
 */

#define OPT_R_BASE  800
#define OPT_R_BED   801
#define OPT_R_GTF 803
#define OPT_R_FUS   804
#define OPT_R_VCF   805
#define OPT_MIXTURE 806
#define OPT_FUZZY   807
#define OPT_R_ENDO  808
#define OPT_R_IND   809
#define OPT_U_BASE  900

#define OPT_U_GTF   902
#define OPT_BAM_1   903
#define OPT_BAM_2   904
#define OPT_PSL_1   905
#define OPT_PSL_2   906
#define OPT_FA_1    907
#define OPT_FA_2    908
#define OPT_U_COV   911
#define OPT_U_FACTS 912
#define OPT_LEVEL   913
#define OPT_U_FILES 914
#define OPT_U_NAMES 915
#define OPT_C_FILES 916
#define OPT_SIGN    917

using namespace Anaquin;


// Shared with other modules
std::string __full_command__;

// Shared with other modules
Path __working__;

// Shared with other modules
Path __output__;


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
    { "test",             TOOL_TEST        },

    { "TransSequin",      TOOL_T_SEQUIN    },
    { "TransAlign",       TOOL_T_ALIGN     },
    { "TransAssembly",    TOOL_T_ASSEMBLY  },
    { "TransExpress",     TOOL_T_EXPRESS   },
    { "TransKExpress",    TOOL_T_KEXPRESS  },
    { "TransDiff",        TOOL_T_DIFF      },
    { "TransKDiff",       TOOL_T_KDIFF     },
    { "TransCount",       TOOL_T_COUNT     },
    { "TransNorm",        TOOL_T_NORM      },
    { "TransIGV",         TOOL_T_IGV       },
    { "TransCoverage",    TOOL_T_COVERAGE  },
    { "TransReport",      TOOL_T_REPORT    },

    { "VarAlign",         TOOL_V_ALIGN     },
    { "VarDiscover",      TOOL_V_DISCOVER  },
    { "VarIGV",           TOOL_V_IGV       },
    { "VarAllele",        TOOL_V_ALLELE    },
    { "VarCoverage",      TOOL_V_COVERAGE  },
    { "VarSubsample",     TOOL_V_SUBSAMPLE },
    { "VarExpress",       TOOL_V_EXPRESS   },
    { "VarKExpress",      TOOL_V_KEXPRESS  },
    { "VarReport",        TOOL_V_REPORT    },

    { "MetaAssembly",     TOOL_M_ASSEMBLY  },
    { "MetaAbund",        TOOL_M_ABUND     },
    { "MetaAbundance",    TOOL_M_ABUND     },
    { "MetaDiff",         TOOL_M_DIFF      },
    { "MetaDifferent",    TOOL_M_DIFF      },
    { "MetaIGV",          TOOL_M_IGV       },
    { "MetaCoverage",     TOOL_M_COVERAGE  },
    { "MetaReport",       TOOL_M_REPORT    },

    { "LadderCopy",       TOOL_L_COPY     },
    { "LadderAbund",      TOOL_L_ABUND     },
    { "LadderAbundance",  TOOL_L_ABUND     },
    { "LadderDiff",       TOOL_L_DIFF      },
    { "LadderDifferent",  TOOL_L_DIFF      },
    { "LadderCoverage",   TOOL_L_COVERAGE  },
    { "LadderReport",     TOOL_L_REPORT    },

    { "FusionDiscover",   TOOL_F_DISCOVER  },
    { "FusionExpress",    TOOL_F_EXPRESS   },
    { "FusionExpression", TOOL_F_EXPRESS   },
    { "FusionIGV",        TOOL_F_IGV       },
    { "FusionCoverage",   TOOL_F_COVERAGE  },
    { "FusionDiff",       TOOL_F_DIFF      },
    { "FusionReport",     TOOL_F_REPORT    },
};

static std::map<Tool, std::set<Option>> _required =
{
    /*
     * Transcriptome Analysis
     */
    
    { TOOL_T_IGV,      { OPT_U_FILES                                                             } },
    { TOOL_T_COVERAGE, { OPT_R_GTF, OPT_U_FILES                                                  } },
    { TOOL_T_ALIGN,    { OPT_R_GTF, OPT_MIXTURE, OPT_U_FACTS, OPT_U_NAMES, OPT_U_FILES           } },
    { TOOL_T_ASSEMBLY, { OPT_R_GTF, OPT_MIXTURE, OPT_U_FACTS, OPT_U_NAMES, OPT_U_FILES           } },
    { TOOL_T_EXPRESS,  { OPT_R_GTF, OPT_MIXTURE, OPT_SOFT, OPT_U_FACTS, OPT_U_NAMES, OPT_U_FILES } },
    { TOOL_T_DIFF,     { OPT_R_GTF, OPT_MIXTURE, OPT_SOFT, OPT_U_FACTS, OPT_U_NAMES, OPT_U_FILES } },
    { TOOL_T_COUNT,    { OPT_SOFT, OPT_U_FACTS, OPT_U_NAMES, OPT_U_FILES                         } },

    /*
     * Ladder Analysis
     */

    { TOOL_L_ABUND,    { OPT_U_FILES, OPT_MIXTURE } },

    /*
     * Metagenomics Analysis
     */

    { TOOL_M_IGV,      { OPT_U_FILES                                                     } },
    { TOOL_M_ASSEMBLY, { OPT_R_BED, OPT_PSL_1, OPT_U_FILES, OPT_SOFT                     } },
    { TOOL_M_ABUND,    { OPT_MIXTURE, OPT_PSL_1, OPT_U_FILES, OPT_SOFT                   } },
    { TOOL_M_COVERAGE, { OPT_R_BED, OPT_U_FILES                                          } },
    { TOOL_M_DIFF,     { OPT_MIXTURE, OPT_PSL_1, OPT_PSL_2, OPT_FA_1, OPT_FA_2, OPT_SOFT } },

    /*
     * Fusion Analysis
     */

    { TOOL_F_DISCOVER, { OPT_R_FUS, OPT_SOFT, OPT_U_FILES               } },
    { TOOL_F_EXPRESS,  { OPT_R_BED, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES  } },
    { TOOL_F_COVERAGE, { OPT_R_BED, OPT_U_FILES                         } },
    { TOOL_F_DIFF,     { OPT_R_BED, OPT_R_FUS, OPT_MIXTURE, OPT_U_FILES } },

    /*
     * Variant Analysis
     */
    
    { TOOL_V_REPORT,    { OPT_U_FILES                                   } },
    { TOOL_V_ALIGN,     { OPT_R_BED, OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_V_ALLELE,    { OPT_R_VCF, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES } },
    { TOOL_V_EXPRESS,   { OPT_R_VCF, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES } },
    { TOOL_V_KEXPRESS,  { OPT_R_IND, OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_V_COVERAGE,  { OPT_R_BED, OPT_U_FILES                        } },
    { TOOL_V_DISCOVER,  { OPT_R_VCF, OPT_R_BED, OPT_SOFT, OPT_U_FILES   } },
    { TOOL_V_IGV,       { OPT_BAM_1                                     } },
    { TOOL_V_SUBSAMPLE, { OPT_R_BED, OPT_R_ENDO, OPT_U_FILES            } },
};

/*
 * Variables used in argument parsing
 */

struct Parsing
{
    // Reference annotation file for synthetic
    FileName rChrT;

    // Reference annotation file for endogenous
    FileName rEndo;
    
    // The path that outputs are written
    std::string path = "output";

    // Input files
    std::vector<FileName> inputs;
    
    // Optional input files
    std::vector<FileName> oInputs;
    
    // Context specific options
    std::map<Option, std::string> opts;
    
    // Signifcance level
    Probability sign;
    
    // Minmium concentration
    double min = 0;
    
    // Maximum concentration
    double max = std::numeric_limits<double>::max();

    // Limit of detection
    double limit;
    
    // How Anaquin is invoked
    std::string command;

    unsigned fuzzy = 0;
    
    Tool tool = 0;

    // Experimental meta-data
    std::shared_ptr<Experiment> exp = std::shared_ptr<Experiment>(new Experiment());
};

// Wrap the variables so that it'll be easier to reset them
static Parsing _p;

struct InvalidOptionException : public std::exception
{
    InvalidOptionException(const std::string &opt) : opt(opt) {}
    
    const std::string opt;
};

struct InvalidValueException : public std::exception
{
    InvalidValueException(const std::string &opt, const std::string &val) : opt(opt), val(val) {}

    const std::string opt, val;
};

struct NoValueError: public InvalidValueException
{
    NoValueError(const std::string &opt) : InvalidValueException(opt, "") {}
};

struct InvalidToolError : public InvalidValueException
{
    InvalidToolError(const std::string &val) : InvalidValueException("-t", val) {}
};

struct MissingOptionError : public std::exception
{
    MissingOptionError(const std::string &opt) : opt(opt) {}
    MissingOptionError(const std::string &opt, const std::string &range) : opt(opt), range(range) {}

    // Option that is missing
    const std::string opt;
    
    // Possible values for the missing option
    const std::string range;
};

struct TooManyOptionsError : public std::runtime_error
{
    TooManyOptionsError(const std::string &msg) : std::runtime_error(msg) {}
};

/*
 * Argument options
 */

static const char *short_options = ":";

static const struct option long_options[] =
{
    { "v", no_argument, 0, OPT_VERSION },

    { "t",    required_argument, 0, OPT_TOOL },

    { "sign",    required_argument, 0, OPT_SIGN    },
    { "ufiles",  required_argument, 0, OPT_U_FILES },
    { "cfiles",  required_argument, 0, OPT_C_FILES },

    { "usams",   required_argument, 0, OPT_BAM_1 },
    { "usam1",   required_argument, 0, OPT_BAM_1 },
    { "usam2",   required_argument, 0, OPT_BAM_2 },
    { "ubams",   required_argument, 0, OPT_BAM_1 },
    { "ubam1",   required_argument, 0, OPT_BAM_1 },
    { "ubam2",   required_argument, 0, OPT_BAM_2 },

    { "factors", required_argument, 0, OPT_U_FACTS },
    { "level",   required_argument, 0, OPT_LEVEL   },
    { "names",   required_argument, 0, OPT_U_NAMES },

    { "rexp",    required_argument, 0, OPT_R_ENDO  },
    { "rbed",    required_argument, 0, OPT_R_BED   },
    { "rgtf",    required_argument, 0, OPT_R_GTF   },
    { "rvcf",    required_argument, 0, OPT_R_VCF   },
    { "rfus",    required_argument, 0, OPT_R_FUS   },
    { "rind",    required_argument, 0, OPT_R_IND   },

    { "ufa",     required_argument, 0, OPT_FA_1   },
    { "ufa1",    required_argument, 0, OPT_FA_1   },
    { "ufa2",    required_argument, 0, OPT_FA_2   },
    { "ugtf",    required_argument, 0, OPT_U_GTF  },
    { "upsl",    required_argument, 0, OPT_PSL_1  },
    { "upsl1",   required_argument, 0, OPT_PSL_1  },
    { "upsl2",   required_argument, 0, OPT_PSL_2  },
    { "ucov",    required_argument, 0, OPT_U_COV  },

    { "min", required_argument, 0, OPT_MIN },
    { "max", required_argument, 0, OPT_MAX },

    { "m",   required_argument, 0, OPT_MIXTURE },
    { "mix", required_argument, 0, OPT_MIXTURE },

    { "fuzzy", required_argument, 0, OPT_FUZZY },
    
    { "o",      required_argument, 0, OPT_PATH },
    { "output", required_argument, 0, OPT_PATH },

    { "soft",   required_argument, 0, OPT_SOFT   },
    { "csoft",  required_argument, 0, OPT_C_SOFT },

    { "f",      required_argument, 0, OPT_FILTER },
    { "filter", required_argument, 0, OPT_FILTER },

    {0, 0, 0, 0 }
};

static std::string optToStr(int opt)
{
    for (const auto &o : long_options)
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
    std::cout << "Anaquin v1.1.1" << std::endl;
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

std::string mixture()
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
    std::cerr << msg << std::endl;
    std::cerr << "*********************************************************************" << std::endl << std::endl;
}

static void readFilters(const FileName &file)
{
    // Empty Implementation
}

template <typename Mixture> void addMix(Mixture mix)
{
    if (mixture().empty())
    {
        return;
    }
    
    std::cout << "[INFO]: Mixture: " << mixture() << std::endl;
    mix(Reader(mixture()));
}

#define CHECK_REF(x) (x != OPT_MIXTURE && x > OPT_R_BASE && x < OPT_U_BASE)

template <typename Reference> void addRef(const ChrID &cID, Reference ref, const FileName &file)
{
    if (cID == ChrT)
    {
        std::cout << "[INFO]: Found synthetic reference"  << std::endl;
    }
    else
    {
        std::cout << "[INFO]: Found experimental reference" << std::endl;
    }

    std::cout << "[INFO]: Reference: " << file << std::endl;
    ref(Reader(file));
}

template <typename Reference> void addRef(const ChrID &cID, Reference ref)
{
    for (const auto &i : _p.opts)
    {
        if (CHECK_REF(i.first))
        {
            addRef(cID, ref, _p.opts[i.first]);
            break;
        }
    }
}

/*
 * This provideds a convenient function for adding reference annotations. Experimental annotations are also added if found.
 */

template <typename Reference> void addRef(Reference ref)
{
    for (const auto &i : _p.opts)
    {
        const auto opt = i.first;
        
        if (CHECK_REF(opt))
        {
            switch (opt)
            {
                case OPT_R_ENDO:
                {
                    addRef(Endo, ref, _p.opts[opt]);
                    break;
                }

                default:
                {
                    addRef(ChrT, ref, _p.opts[opt]);
                    break;
                }
            }
        }
    }
}

// Apply a reference source given where it comes from
template <typename Reference> void applyRef(Reference ref, Option opt)
{
    std::cout << "[INFO]: Reference: " << _p.opts[opt] << std::endl;
    ref(Reader(_p.opts[opt]));
}

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
            applyRef(ref, opt);
            break;
        }
    }
}

template <typename Analyzer, typename F> void analyzeF(F f, typename Analyzer::Options o)
{
    const auto path = _p.path;

    // This might be needed for scripting
    __full_command__ = _p.command;

#ifndef DEBUG
    o.writer = std::shared_ptr<FileWriter>(new FileWriter(path));
    o.logger = std::shared_ptr<FileWriter>(new FileWriter(path));
    o.output = std::shared_ptr<TerminalWriter>(new TerminalWriter());
    o.logger->open("anaquin.log");
#endif

    char cwd[1024];

    if (getcwd(cwd, sizeof(cwd)))
    {
        o.working = path; //std::string(cwd) + "/" + path;
    }
    else
    {
        o.working = path;
    }
    
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    o.info(_p.command);
    o.info(date());
    o.info("Path: " + path);

    o.rChrT = _p.rChrT;
    o.rEndo = _p.rEndo.empty() ? "NA" : _p.rEndo;

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

/*
 * Template functions for analyzing
 */

template <typename Analyzer, typename Files> void analyze(const Files &files, typename Analyzer::Options o = typename Analyzer::Options())
{
    // Copying over for the experiment
    o.exp = _p.exp;

    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(files, o);
    }, o);
}

// Analyze for a single sample
template <typename Analyzer> void analyze_1(Option x, typename Analyzer::Options o = typename Analyzer::Options())
{
    return analyze<Analyzer>(_p.opts.at(x), o);
}

// Analyze for a single sample with fuzzy matching
template <typename Analyzer> void analyzeFuzzy(typename Analyzer::Options o = typename Analyzer::Options())
{
    o.fuzzy = _p.fuzzy;
    return analyze<Analyzer>(_p.inputs[0], o);
}

// Analyze for two samples
template < typename Analyzer> void analyze_2(typename Analyzer::Options o = typename Analyzer::Options())
{
    return analyzeF<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.inputs[0], _p.inputs[1], o);
    }, o);
}

/*
 * Functions for parsing string to enums
 */

template <typename T> T parseEnum(const std::string &key, const std::string &str, const std::map<std::string, T> &m)
{
    for (const auto &i : m)
    {
        if (strcasecmp(i.first.c_str(), str.c_str()) == 0)
        {
            return i.second;
        }
    }

    throw InvalidValueException(str, key);
};

template <typename T> T parseCSoft(const std::string &str, const std::string &key)
{
    const static std::map<std::string, T> m =
    {
        { "HTSeqCount", T::HTSeqCount },
    };

    return parseEnum(key, str, m);
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
            throw std::runtime_error("Failed to parse: " + std::to_string(r));
        }
    };
    
    auto checkPath = [&](const Path &path)
    {
        if (path[0] == '/')
        {
            return path;
        }
        else
        {
            return __working__ + "/" + path;
        }
    };
    
    auto checkFile = [&](const FileName &file)
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

    // Prevent error message to stderr
    opterr = 0;
    
    while ((next = getopt_long_only(argc, argv, short_options, long_options, &index)) != -1)
    {
        if (next == ':')
        {
            throw NoValueError(argv[n+1]);
        }
        else if (next < OPT_TOOL)
        {
            throw InvalidOptionException(argv[n+1]);
        }

        opts.push_back(next);

        // Whether this option has an value
        const auto hasValue = optarg;
        
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

    for (auto i = 0; i < opts.size(); i++)
    {
        auto opt = opts[i];
        auto val = vals[i];

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

            case OPT_FUZZY: { parseInt(val, _p.fuzzy); break; }

            /*
             * The following options can only be validated by the tool
             */

            case OPT_SOFT:
            case OPT_LEVEL:
            case OPT_C_SOFT: { _p.opts[opt] = val; break; }

            case OPT_U_FILES:
            {
                std::vector<FileName> temp;
                Tokens::split(val, ",", temp);

                for (auto i = 0; i < temp.size(); i++)
                {
                    checkFile(_p.opts[opt] = temp[i]);
                    _p.exp->addFile(_p.opts[opt]);
                    _p.inputs.push_back(temp[i]);
                }
                
                break;
            }
             
            /*
             * Parse for the optional input files
             */
                
            case OPT_C_FILES:
            {
                Tokens::split(val, ",", _p.oInputs);
                
                for (auto i = _p.oInputs.size(); i-- > 0;)
                {
                    checkFile(_p.opts[opt] = _p.oInputs[i]);
                }

                break;
            }
                
            case OPT_BAM_1:
            {
                Tokens::split(val, ",", _p.inputs);

                for (auto i = _p.inputs.size(); i-- > 0;)
                {
                    checkFile(_p.opts[opt] = _p.inputs[i]);
                }
                
                break;
            }

            case OPT_U_FACTS:
            {
                _p.exp->addFactors(_p.opts[opt] = val);
                break;
            }

            case OPT_U_NAMES:
            {
                _p.exp->addNames(_p.opts[opt] = val);
                break;
            }

            case OPT_FA_1:
            case OPT_FA_2:
            case OPT_U_COV:
            case OPT_U_GTF:
            case OPT_PSL_2:
            case OPT_BAM_2:
            case OPT_PSL_1:
            case OPT_MIXTURE: { checkFile(_p.opts[opt] = val); break; }

            case OPT_R_IND:
            case OPT_R_FUS:
            case OPT_R_VCF:
            case OPT_R_BED:
            case OPT_R_GTF:
            {
                checkFile(_p.opts[opt] = _p.rChrT = val);
                break;
            }
                
            case OPT_R_ENDO:
            {
                checkFile(_p.opts[opt] = _p.rEndo = val);
                break;
            }
                
            case OPT_PATH:    { _p.path = val;              break; }
            case OPT_FILTER:  { readFilters(val);           break; }
            case OPT_MAX:     { parseDouble(val, _p.max);   break; }
            case OPT_MIN:     { parseDouble(val, _p.min);   break; }
            case OPT_LOS:     { parseDouble(val, _p.limit); break; }
            case OPT_SIGN:    { parseDouble(val, _p.sign);  break; }

            default:
            {
                throw InvalidOptionException(argv[index]);
            }
        }
    }

    __output__ = _p.path = checkPath(_p.path);

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
            throw MissingOptionError("-" + optToStr(*required.begin()));
        }
    }

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
        case TOOL_T_COVERAGE:
        {
            std::cout << "[INFO]: Transcriptome Analysis" << std::endl;

            if (_p.tool != TOOL_T_IGV)
            {
                addRef(std::bind(&Standard::addTRef, &s, std::placeholders::_1));
                
                switch (_p.tool)
                {
                    case TOOL_T_DIFF:
                    {
                        addMix(std::bind(&Standard::addTDMix, &s, std::placeholders::_1));
                        break;
                    }

                    default:
                    {
                        addMix(std::bind(&Standard::addTSMix, &s, std::placeholders::_1));
                        break;
                    }
                }

                s.r_trans.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_T_SEQUIN:   { printMixture();                  break; }
                case TOOL_T_COVERAGE: { analyze_1<TCoverage>(OPT_BAM_1); break; }

                case TOOL_T_ALIGN:
                {
                    analyze<TAlign>(_p.inputs);
                    break;
                }

                case TOOL_T_ASSEMBLY:
                {
                    analyze<TAssembly>(_p.inputs);
                    break;
                }

                case TOOL_T_EXPRESS:
                {
                    auto parseLevel = [&](const std::string &key, const std::string &str)
                    {
                        const static std::map<std::string, TExpress::Metrics> m =
                        {
                            { "gene",    TExpress::Metrics::Gene    },
                            { "isoform", TExpress::Metrics::Isoform },
                        };
                        
                        return parseEnum(key, str, m);
                    };

                    auto parseSoft = [&](const std::string &key, const std::string &str)
                    {
                        const static std::map<std::string, TExpress::Software> m =
                        {
                            { "kallisto",  TExpress::Software::Kallisto  },
                            { "cufflinks", TExpress::Software::Cufflinks },
                            { "stringtie", TExpress::Software::StringTie },
                        };
                        
                        return parseEnum(key, str, m);
                    };

                    TExpress::Options o;
                    
                    o.soft = parseSoft("soft", _p.opts.at(OPT_SOFT));
                    
                    if (_p.opts.count(OPT_LEVEL))
                    {
                        o.metrs = parseLevel("level", _p.opts[OPT_LEVEL]);
                    }
                    
                    if (o.soft == TExpress::Software::Kallisto)
                    {
                        o.metrs = TExpress::Metrics::Isoform;
                    }

                    analyze<TExpress>(_p.inputs, o);
                    break;
                }

                case TOOL_T_COUNT:
                {                    
                    break;
                }
                    
                case TOOL_T_DIFF:
                {
                    auto parseMetrs = [&](const std::string &str)
                    {
                        const static std::map<std::string, TDiffs::Metrics> m =
                        {
                            { "gene",    TDiffs::Metrics::Gene    },
                            { "isoform", TDiffs::Metrics::Isoform },
                        };
                        
                        return parseEnum("levels", str, m);
                    };

                    auto parseSoft = [&](const std::string &str)
                    {
                        const static std::map<std::string, TDiffs::DiffSoft> m =
                        {
                            { "edgeR",    TDiffs::DiffSoft::edgeR    },
                            { "deseq2",   TDiffs::DiffSoft::DESeq2   },
                            { "cuffdiff", TDiffs::DiffSoft::Cuffdiff },
                            { "sleuth",   TDiffs::DiffSoft::ParserSleuth },
                        };
                        
                        return parseEnum("soft", str, m);
                    };
                    
                    TDiffs::Options o;

                    o.dSoft = parseSoft(_p.opts[OPT_SOFT]);

                    if (_p.opts.count(OPT_LEVEL))
                    {
                        o.metrs = parseMetrs(_p.opts[OPT_LEVEL]);
                    }
                    
                    if (o.dSoft != TDiffs::DiffSoft::Cuffdiff)
                    {
                        o.metrs = TDiffs::Metrics::Isoform;
                    }

                    /*
                     * Optional count tables (eg: HTSeqCount)
                     */
                    
                    if (_p.opts.count(OPT_C_FILES))
                    {
                        o.cSoft  = parseCSoft<TDiffs::CountSoft>(_p.opts[OPT_C_SOFT], "csoft");
                        o.counts = _p.oInputs;
                    }

                    analyze_1<TDiffs>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_T_IGV: { viewer<TViewer>(); break; }
            }

            break;
        }

        case TOOL_F_IGV:
        case TOOL_F_DIFF:
        case TOOL_F_EXPRESS:
        case TOOL_F_DISCOVER:
        case TOOL_F_COVERAGE:
        {
            auto parseAligner = [&](const std::string &str)
            {
                const static std::map<std::string, FusionCaller> m =
                {
                    { "Star"        , FusionCaller::StarFusion   },
                    { "StarFusion"  , FusionCaller::StarFusion   },
                    { "TopHat"      , FusionCaller::TopHatFusion },
                    { "TopHatFusion", FusionCaller::TopHatFusion },
                    { "Kallisto"    , FusionCaller::Kallisto     },
                };

                return parseEnum("soft", str, m);
            };
            
            std::cout << "[INFO]: Fusion Analysis" << std::endl;

            switch (_p.tool)
            {
                case TOOL_F_EXPRESS:
                {
                    addMix(std::bind(&Standard::addFMix,      &s, std::placeholders::_1));
                    applyRef(std::bind(&Standard::addFSplice, &s, std::placeholders::_1), OPT_R_BED);
                    break;
                }

                case TOOL_F_DIFF:
                {
                    addMix(std::bind(&Standard::addFMix,    &s, std::placeholders::_1));
                    applyRef(std::bind(&Standard::addFRef,    &s, std::placeholders::_1), OPT_R_FUS);
                    applyRef(std::bind(&Standard::addFSplice, &s, std::placeholders::_1), OPT_R_BED);
                    break;
                }

                case TOOL_F_COVERAGE:
                {
                    applyRef(std::bind(&Standard::addFStd, &s, std::placeholders::_1), OPT_R_BED);
                    break;
                }

                case TOOL_F_DISCOVER:
                {
                    applyRef(std::bind(&Standard::addFRef, &s, std::placeholders::_1));
                    addMix(std::bind(&Standard::addFMix,   &s, std::placeholders::_1));
                    break;
                }
                    
                default: { break; }
            }
            
            if (_p.tool != TOOL_F_IGV)
            {
                Standard::instance().r_fus.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_F_IGV:      { viewer<FViewer>();               break; }
                case TOOL_F_COVERAGE: { analyze_1<FCoverage>(OPT_BAM_1); break; }

                case TOOL_F_EXPRESS:
                {
                    FExpress::Options o;
                    o.caller = parseAligner(_p.opts.at(OPT_SOFT));

                    analyze_1<FExpress>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_F_DIFF:
                {
                    analyze_2<FDiff>();
                    break;
                }

                case TOOL_F_DISCOVER:
                {
                    FDiscover::Options o;
                    o.caller = parseAligner(_p.opts.at(OPT_SOFT));

                    analyzeFuzzy<FDiscover>(o);
                    break;
                }
            }

            break;
        }

        case TOOL_L_COPY:
        case TOOL_L_DIFF:
        case TOOL_L_ABUND:
        case TOOL_L_COVERAGE:
        {
            std::cout << "[INFO]: Ladder Analysis" << std::endl;

            addMix(std::bind(&Standard::addLMix, &s, std::placeholders::_1));
            Standard::instance().r_lad.finalize();

            switch (_p.tool)
            {
                case TOOL_L_COPY:     { analyze_1<LCopy>(OPT_U_FILES);             break; }
                case TOOL_L_ABUND:    { analyze_1<LAbund>(OPT_U_FILES);            break; }
                case TOOL_L_DIFF:     { /*analyze_2<LDiffs>(OPT_BAM_1, OPT_BAM_2); break;*/ }
                case TOOL_L_COVERAGE: { analyze_1<LCoverage>(OPT_U_FILES);         break; }
            }

            break;
        }

        case TOOL_V_IGV:
        case TOOL_V_ALIGN:
        case TOOL_V_ALLELE:
        case TOOL_V_REPORT:
        case TOOL_V_EXPRESS:
        case TOOL_V_DISCOVER:
        case TOOL_V_COVERAGE:
        case TOOL_V_KEXPRESS:
        case TOOL_V_SUBSAMPLE:
        {
            auto parseCaller = [&](const std::string &str)
            {
                const static std::map<std::string, Caller> m =
                {
                    { "gatk"   ,  Caller::GATK     },
                    { "VarScan",  Caller::VarScan  },
                    { "VarScan2", Caller::VarScan  },
                };

                return parseEnum("soft", str, m);
            };
            
            auto parseExpress = [&](const std::string &str)
            {
                const static std::map<std::string, VExpress::Software> m =
                {
                    { "kallisto", VExpress::Software::Kallisto },
                };

                return parseEnum("soft", str, m);
            };
            
            std::cout << "[INFO]: Variant Analysis" << std::endl;

            if (_p.tool != TOOL_V_IGV && _p.tool != TOOL_V_REPORT)
            {
                switch (_p.tool)
                {
                    case TOOL_V_SUBSAMPLE:
                    {
                        applyRef(std::bind(&Standard::addStd,    &s, std::placeholders::_1));
                        applyRef(std::bind(&Standard::addInters, &s, std::placeholders::_1), OPT_R_ENDO);
                        break;
                    }

                    case TOOL_V_ALIGN:
                    case TOOL_V_COVERAGE:
                    {
                        applyRef(std::bind(&Standard::addStd, &s, std::placeholders::_1),    OPT_R_BED);
                        applyRef(std::bind(&Standard::addInters, &s, std::placeholders::_1), OPT_R_ENDO);
                        break;
                    }

                    case TOOL_V_ALLELE:
                    case TOOL_V_EXPRESS:
                    {
                        applyRef(std::bind(&Standard::addVar, &s, std::placeholders::_1));
                        break;
                    }

                    case TOOL_V_DISCOVER:
                    {
                        applyRef(std::bind(&Standard::addStd, &s, std::placeholders::_1), OPT_R_BED);
                        applyRef(std::bind(&Standard::addVar, &s, std::placeholders::_1), OPT_R_VCF);
                        break;
                    }

                    default: { break; }
                }

                addMix(std::bind(&Standard::addVMix, &s, std::placeholders::_1));
                Standard::instance().r_var.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_V_IGV:       { viewer<VViewer>();                 break; }
                case TOOL_V_ALIGN:     { analyze_1<VAlign>(OPT_U_FILES);    break; }
                case TOOL_V_COVERAGE:  { analyze_1<VCoverage>(OPT_U_FILES); break; }

                case TOOL_V_KEXPRESS:
                {
                    VKExpress::Options o;
                    o.file = _p.opts[OPT_R_IND];
                    analyze_2<VKExpress>(o);
                    break;
                }

                case TOOL_V_EXPRESS:
                {
                    VExpress::Options o;
                    o.soft = parseExpress(_p.opts.at(OPT_SOFT));
                    
                    analyze_1<VExpress>(OPT_U_FILES, o);
                    break;
                }
                    
                case TOOL_V_ALLELE:
                {
                    VAllele::Options o;
                    o.caller = parseCaller(_p.opts.at(OPT_SOFT));

                    analyze_1<VAllele>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_V_DISCOVER:
                {
                    VDiscover::Options o;
                    o.caller = parseCaller(_p.opts.at(OPT_SOFT));

                    analyze_1<VDiscover>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_V_SUBSAMPLE: { analyze_1<VSample>(OPT_U_FILES); break; }
            }

            break;
        }

        case TOOL_M_IGV:
        case TOOL_M_DIFF:
        case TOOL_M_ABUND:
        case TOOL_M_ASSEMBLY:
        case TOOL_M_COVERAGE:
        {
            auto parseAssembler = [&](const std::string &str)
            {
                const static std::map<std::string, MetaAssembler> m =
                {
                    { "velvet"  , MetaAssembler::Velvet  },
                    { "ray",      MetaAssembler::RayMeta },
                    { "raymeta",  MetaAssembler::RayMeta },
                    { "ray-meta", MetaAssembler::RayMeta },
                };
                
                return parseEnum("soft", str, m);
            };

            std::cout << "[INFO]: Metagenomics Analysis" << std::endl;

            if (_p.tool != TOOL_M_IGV)
            {
                switch (_p.tool)
                {
                    case TOOL_M_ASSEMBLY:
                    case TOOL_M_COVERAGE:
                    {
                        applyRef(std::bind(&Standard::addMRef, &s, std::placeholders::_1));
                        break;
                    }

                    default: { break; }
                }

                if (_p.tool == TOOL_M_ABUND || _p.tool == TOOL_M_DIFF)
                {
                    addMix(std::bind(&Standard::addMMix, &s, std::placeholders::_1));
                }
                
                Standard::instance().r_meta.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_M_IGV:      { viewer<FViewer>();               break; }
                case TOOL_M_COVERAGE: { analyze_1<MCoverage>(OPT_BAM_1); break; }

                case TOOL_M_DIFF:
                case TOOL_M_ABUND:
                case TOOL_M_ASSEMBLY:
                {
                    // Only defined for certain assemblers
                    FileName conts;
                    
                    const auto tool = parseAssembler(_p.opts.at(OPT_SOFT));
                    
                    if (tool == MetaAssembler::RayMeta)
                    {
                        if (!_p.opts.count(OPT_U_COV))
                        {
                            throw MissingOptionError("ucov");
                        }
                        
                        conts = _p.opts.at(OPT_U_COV);
                    }
                    
                    switch (_p.tool)
                    {
                        case TOOL_M_DIFF:
                        {
                            MDiffs::Options o;
                            
                            o.pA = _p.opts.at(OPT_PSL_1);
                            o.pB = _p.opts.at(OPT_PSL_2);
                            
                            //analyze_2<MDiffs>(OPT_FA_1, OPT_FA_2, o);
                            break;
                        }
                            
                        case TOOL_M_ASSEMBLY:
                        {
                            MAssembly::Options o;
                            
                            o.tool    = tool;
                            o.contigs = conts;
                            
                            // An alignment file is needed to identify contigs
                            o.psl = _p.opts.at(OPT_PSL_1);
                            
                            analyze_1<MAssembly>(OPT_FA_1, o);
                            break;
                        }
                            
                        case TOOL_M_ABUND:
                        {
                            MAbundance::Options o;
                            
                            o.tool    = tool;
                            o.contigs = conts;
                            
                            // An alignment file is needed to identify contigs
                            o.psl = _p.opts.at(OPT_PSL_1);
                            
                            analyze_1<MAbundance>(OPT_FA_1, o);
                            break;
                        }

                        default: { break; }
                    }
                }
            }

            break;
        }

        default: { assert(false); }
    }
}

int parse_options(int argc, char ** argv)
{
    char cwd[1024];
    
    if (getcwd(cwd, sizeof(cwd)))
    {
        __working__ = cwd;
    }
    
    try
    {
        parse(argc, argv);
        return 0;
    }
    catch (const NoValueError &ex)
    {
        printError("Invalid command. Need to specify " + ex.opt + ".");
    }
    catch (const InvalidToolError &ex)
    {
        printError("Invalid tool: " + ex.val + ". Please check the user manual and try again.");
    }
    catch (const InvalidOptionException &ex)
    {
        printError((boost::format("Unknown option: %1%") % ex.opt).str());
    }
    catch (const InvalidValueException &ex)
    {
        printError((boost::format("%1% not expected for option -%2%. Please check and try again.") % ex.opt % ex.val).str());
    }
    catch (const MissingOptionError &ex)
    {
        const auto format = "Mandatory option is missing. Please specify %1%.";
        printError((boost::format(format) % ex.opt).str());
    }
    catch (const InvalidFileError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid file: " % ex.file).str());
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