#include <map>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <strings.h>
#include <execinfo.h>

#include "TransQuin/t_diff.hpp"
#include "TransQuin/t_kdiff.hpp"
#include "TransQuin/t_align.hpp"
#include "TransQuin/t_viewer.hpp"
#include "TransQuin/t_report.hpp"
#include "TransQuin/t_sample.hpp"
#include "TransQuin/t_express.hpp"
#include "TransQuin/t_kexpress.hpp"
#include "TransQuin/t_assembly.hpp"
#include "TransQuin/t_coverage.hpp"

#include "VarQuin/v_align.hpp"
#include "VarQuin/v_allele.hpp"
#include "VarQuin/v_viewer.hpp"
#include "VarQuin/v_report.hpp"
#include "VarQuin/v_sample.hpp"
#include "VarQuin/v_express.hpp"
#include "VarQuin/v_kallele.hpp"
#include "VarQuin/v_discover.hpp"
#include "VarQuin/v_coverage.hpp"
#include "VarQuin/v_kexpress.hpp"

#include "MetaQuin/m_diff.hpp"
#include "MetaQuin/m_express.hpp"
#include "MetaQuin/m_assembly.hpp"
#include "MetaQuin/m_coverage.hpp"

#include "StructQuin/s_discover.hpp"

#include "LadQuin/l_norm.hpp"
#include "LadQuin/l_copy.hpp"
#include "LadQuin/l_diff.hpp"
#include "LadQuin/l_coverage.hpp"

#include "FusQuin/f_diff.hpp"
#include "FusQuin/f_viewer.hpp"
#include "FusQuin/f_normal.hpp"
#include "FusQuin/f_fusion.hpp"
#include "FusQuin/f_discover.hpp"
#include "FusQuin/f_coverage.hpp"

#include "parsers/parser_cdiff.hpp"
#include "parsers/parser_quast.hpp"
#include "parsers/parser_cufflink.hpp"

#include "writers/pdf_writer.hpp"
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
#define TOOL_T_ALIGN     266
#define TOOL_T_ASSEMBLY  267
#define TOOL_T_EXPRESS   268
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
#define TOOL_T_SUBSAMPLE 281
#define TOOL_M_EXPRESS   282
#define TOOL_M_ASSEMBLY  283
#define TOOL_M_DIFF      284
#define TOOL_M_IGV       285
#define TOOL_M_COVERAGE  286
#define TOOL_L_NORM       287
#define TOOL_L_DIFF      288
#define TOOL_L_COVERAGE  289
#define TOOL_F_DISCOVER  290
#define TOOL_F_NORMAL    291
#define TOOL_F_FUSION    292
#define TOOL_F_IGV       293
#define TOOL_F_COVERAGE  294
#define TOOL_F_DIFF      295
#define TOOL_L_COPY      296
#define TOOL_V_EXPRESS   297
#define TOOL_V_KEXPRESS  298
#define TOOL_T_KEXPRESS  299
#define TOOL_T_KDIFF     304
#define TOOL_V_KALLELE   306
#define TOOL_S_DISCOVER  307

/*
 * Options specified in the command line
 */

#define OPT_TEST     320
#define OPT_TOOL     321
#define OPT_PATH     325
#define OPT_VERSION  338
#define OPT_SOFT     339
#define OPT_C_SOFT   340

/*
 * References - OPT_R_BASE to OPT_U_BASE
 */

#define OPT_R_BASE  800
#define OPT_R_BED   801
#define OPT_METHOD  802
#define OPT_R_GTF   803
#define OPT_R_FUS   804
#define OPT_R_VCF   805
#define OPT_MIXTURE 806
#define OPT_FUZZY   807
#define OPT_R_GENO  808
#define OPT_R_IND   809
#define OPT_U_BASE  900

#define OPT_U_FILES 909
#define OPT_C_FILES 910
#define OPT_REPORT  911

using namespace Anaquin;


// Shared with other modules
std::string __full_command__;

// Shared with other modules
Path __working__;

// Shared with other modules
Path __output__;

// Full path where Anaquin is
Path __anaquin__;



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
    { "test",           TOOL_TEST        },

    { "TransAlign",     TOOL_T_ALIGN     },
    { "TransAssembly",  TOOL_T_ASSEMBLY  },
    { "TransExpress",   TOOL_T_EXPRESS   },
    { "TransKExpress",  TOOL_T_KEXPRESS  },
    { "TransDiff",      TOOL_T_DIFF      },
    { "TransKDiff",     TOOL_T_KDIFF     },
    { "TransNorm",      TOOL_T_NORM      },
    { "TransIGV",       TOOL_T_IGV       },
    { "TransCoverage",  TOOL_T_COVERAGE  },
    { "TransSubsample", TOOL_T_SUBSAMPLE },

    { "VarAlign",       TOOL_V_ALIGN     },
    { "VarDiscover",    TOOL_V_DISCOVER  },
    { "VarIGV",         TOOL_V_IGV       },
    { "VarAllele",      TOOL_V_ALLELE    },
    { "VarCoverage",    TOOL_V_COVERAGE  },
    { "VarSubsample",   TOOL_V_SUBSAMPLE },
    { "VarExpress",     TOOL_V_EXPRESS   },
    { "VarKExpress",    TOOL_V_KEXPRESS  },
    { "VarKAllele",     TOOL_V_KALLELE   },

    { "MetaAssembly",   TOOL_M_ASSEMBLY  },
    { "MetaExpress",    TOOL_M_EXPRESS   },
    { "MetaDiff",       TOOL_M_DIFF      },
    { "MetaIGV",        TOOL_M_IGV       },
    { "MetaCoverage",   TOOL_M_COVERAGE  },

    { "StructDiscover", TOOL_S_DISCOVER },
    
    { "LadCopy",     TOOL_L_COPY     },
    { "LadNorm",     TOOL_L_NORM     },
    { "LadDiff",     TOOL_L_DIFF     },
    { "LadCoverage", TOOL_L_COVERAGE },

    { "FusDiscover", TOOL_F_DISCOVER },
    { "FusNormal",   TOOL_F_NORMAL   },
    { "FusFusion",   TOOL_F_FUSION   },
    { "FusIGV",      TOOL_F_IGV      },
    { "FusCoverage", TOOL_F_COVERAGE },
    { "FusDiff",     TOOL_F_DIFF     },
};

static std::map<Tool, std::set<Option>> _required =
{
    /*
     * Transcriptome Analysis
     */
    
    { TOOL_T_IGV,       { OPT_U_FILES                                   } },
    { TOOL_T_DIFF,      { OPT_R_GTF, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES } },
    { TOOL_T_ASSEMBLY,  { OPT_R_GTF, OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_T_COVERAGE,  { OPT_R_GTF, OPT_U_FILES                        } },
    { TOOL_T_KEXPRESS,  { OPT_R_IND, OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_T_KDIFF,     { OPT_R_IND, OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_T_ALIGN,     { OPT_R_GTF, OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_T_EXPRESS,   { OPT_R_GTF, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES } },
    { TOOL_T_SUBSAMPLE, { OPT_MIXTURE, OPT_U_FILES } },
    
    /*
     * Ladder Analysis
     */

    { TOOL_L_NORM, { OPT_U_FILES, OPT_MIXTURE } },

    /*
     * Structural analysis
     */
    
    { TOOL_S_DISCOVER, { OPT_R_BED, OPT_U_FILES } },

    /*
     * Metagenomics Analysis
     */

    { TOOL_M_IGV,      { OPT_U_FILES } },
    { TOOL_M_ASSEMBLY, { OPT_R_BED,   OPT_U_FILES, OPT_SOFT } },
    { TOOL_M_EXPRESS,  { OPT_MIXTURE, OPT_U_FILES, OPT_SOFT } },
    { TOOL_M_COVERAGE, { OPT_R_BED, OPT_U_FILES             } },
    { TOOL_M_DIFF,     { OPT_MIXTURE, OPT_U_FILES, OPT_SOFT } },

    /*
     * Fusion Analysis
     */

    { TOOL_F_DISCOVER, { OPT_R_BED, OPT_SOFT, OPT_U_FILES               } },
    { TOOL_F_NORMAL,   { OPT_R_BED, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES  } },
    { TOOL_F_FUSION,   { OPT_R_BED, OPT_MIXTURE, OPT_SOFT, OPT_U_FILES  } },
    { TOOL_F_COVERAGE, { OPT_R_BED, OPT_U_FILES                         } },
    { TOOL_F_DIFF,     { OPT_R_BED, OPT_R_FUS, OPT_MIXTURE, OPT_U_FILES } },

    /*
     * Variant Analysis
     */
    
    { TOOL_V_IGV,       { OPT_U_FILES                                     } },
    { TOOL_V_COVERAGE,  { OPT_R_BED,   OPT_U_FILES                        } },
    { TOOL_V_EXPRESS,   { OPT_MIXTURE, OPT_SOFT, OPT_R_VCF, OPT_U_FILES   } },
    { TOOL_V_SUBSAMPLE, { OPT_R_BED,   OPT_R_GENO,  OPT_U_FILES           } },
    { TOOL_V_ALIGN,     { OPT_R_BED,   OPT_MIXTURE, OPT_U_FILES           } },
    { TOOL_V_ALLELE,    { OPT_R_VCF,   OPT_MIXTURE, OPT_SOFT, OPT_U_FILES } },
    { TOOL_V_KEXPRESS,  { OPT_SOFT,    OPT_R_IND,   OPT_MIXTURE, OPT_U_FILES } },
    { TOOL_V_KALLELE,   { OPT_SOFT,    OPT_R_IND,   OPT_MIXTURE, OPT_U_FILES } },
    { TOOL_V_DISCOVER,  { OPT_R_VCF,   OPT_R_BED, OPT_SOFT, OPT_U_FILES, OPT_MIXTURE } },
};

/*
 * Variables used in argument parsing
 */

struct Parsing
{
    // Reference annotation file for synthetic
    FileName rChrT;

    // Reference annotation file for endogenous
    FileName rGeno;
    
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
    
    // How Anaquin is invoked
    std::string command;

    unsigned fuzzy = 0;
    
    Tool tool = 0;
    
    // Generate a PDF report?
    bool pdf = false;
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

struct NotSingleInputError : public std::exception
{
    // Empty Implementation
};

struct NotDoubleInputError : public std::exception
{
    // Empty Implementation
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

    { "t",       required_argument, 0, OPT_TOOL },
    { "tool",    required_argument, 0, OPT_TOOL },

    { "ufiles",  required_argument, 0, OPT_U_FILES },
    { "cfiles",  required_argument, 0, OPT_C_FILES },

    { "m",       required_argument, 0, OPT_MIXTURE },
    { "mix",     required_argument, 0, OPT_MIXTURE },
    { "meth",    required_argument, 0, OPT_METHOD  },

    { "rgen",    required_argument, 0, OPT_R_GENO  },
    { "rbed",    required_argument, 0, OPT_R_BED   },
    { "rgtf",    required_argument, 0, OPT_R_GTF   },
    { "rvcf",    required_argument, 0, OPT_R_VCF   },
    { "rfus",    required_argument, 0, OPT_R_FUS   },
    { "rind",    required_argument, 0, OPT_R_IND   },

    { "fuzzy",   required_argument, 0, OPT_FUZZY },
    
    { "o",       required_argument, 0, OPT_PATH },
    { "output",  required_argument, 0, OPT_PATH },

    { "soft",    required_argument, 0, OPT_SOFT   },
    { "csoft",   required_argument, 0, OPT_C_SOFT },

    { "report",  required_argument, 0, OPT_REPORT },
    
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
    extern Scripts Manual();
    std::cout << Manual() << std::endl;
}

static void printVersion()
{
    std::cout << "v1.1.1" << std::endl;
}

template <typename F> bool testFile(const FileName &x, F f)
{
    try
    {
        // Anything found?
        if (!f(x))
        {
            return false;
        }
    }
    catch (...)
    {
        return false;
    }
    
    return true;
}

template <typename F1, typename F2> std::vector<FileName> sortInputs(const FileName &x,
                                                                     const FileName &y,
                                                                     const FileName &z,
                                                                     F1 f1,
                                                                     F2 f2)
{
    #define ORDER(a,b,c) std::vector<FileName> { a,b,c };

    if (testFile(x, f1))
    {
        /*
         * [x,?,?]
         */
        
        if (testFile(y, f2))
        {
            return ORDER(x,y,z);
        }
        else
        {
            return ORDER(x,z,y);
        }
    }
    else if (testFile(y, f1))
    {
        /*
         * [y,?,?]
         */
        
        if (testFile(x, f2))
        {
            return ORDER(y,x,z);
        }
        else
        {
            return ORDER(y,z,x);
        }
    }
    else
    {
        /*
         * [z,?,?]
         */
        
        if (testFile(x, f2))
        {
            return ORDER(z,x,y);
        }
        else
        {
            return ORDER(z,y,x);
        }
    }
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

FileName mixture()
{
    return _p.opts[OPT_MIXTURE];
}

#define CHECK_REF(x) (x != OPT_MIXTURE && x > OPT_R_BASE && x < OPT_U_BASE)

FileName refFile()
{
    for (const auto &i : _p.opts)
    {
        const auto opt = i.first;
        
        if (CHECK_REF(opt))
        {
            return _p.opts[opt];
        }
    }
    
    throw "No reference file found";
}

static void printError(const std::string &msg)
{
    std::cerr << std::endl;
    std::cerr << "*********************************************************************" << std::endl;
    std::cerr << msg << std::endl;
    std::cerr << "*********************************************************************" << std::endl << std::endl;
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

template <typename Reference> void addRef(const ChrID &cID, Reference ref, const FileName &file)
{
    if (cID == ChrT)
    {
        std::cout << "[INFO]: Found synthetic reference"  << std::endl;
    }
    else
    {
        std::cout << "[INFO]: Found genome reference" << std::endl;
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
 * This provideds a convenient function for adding reference annotations. Genomic annotations are also
 * added if found.
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
                case OPT_R_IND:
                {
                    continue;
                }
                    
                case OPT_R_GENO:
                {
                    addRef(Geno, ref, _p.opts[opt]);
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

static void saveRef()
{
    for (const auto &i : _p.opts)
    {
        if (CHECK_REF(i.first) && !i.second.empty())
        {
            system(("cp " + i.second + " " + __output__).c_str());
        }
    }
}

// Apply a reference source given where it comes from
template <typename Reference> void applyRef(Reference ref, Option opt)
{
    std::cout << "[INFO]: Reference: " << _p.opts[opt] << std::endl;
    
    if (!_p.opts[opt].empty())
    {
        ref(Reader(_p.opts[opt]));
    }
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

template <typename Analyzer, typename F> void startAnalysis(F f, typename Analyzer::Options o)
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

    o.report = std::shared_ptr<PDFWriter>(new PDFWriter());
    o.work = path;
    
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);

    o.info(_p.command);
    o.info(date());
    o.info("Path: " + path);

    o.rChrT = _p.rChrT;
    o.rGeno = _p.rGeno;

    std::clock_t begin = std::clock();

    f(o);
    
    std::clock_t end = std::clock();

    const auto elapsed = (boost::format("Completed. %1% seconds.") % (double(end - begin) / CLOCKS_PER_SEC)).str();
    o.info(elapsed);

#ifndef DEBUG
    o.logger->close();
#endif

    if (_p.pdf)
    {
        o.report->create(o.work);
    }

    // Always save the reference files
    saveRef();
}

template <typename Report> void report_1(typename Report::Options o = typename Report::Options())
{
    if (_p.inputs.size() != 1)
    {
        throw NotSingleInputError();
    }
    
    o.mix = mixture();
    o.index = _p.opts[OPT_R_IND];
    
    return startAnalysis<Report>([&](const typename Report::Options &o)
    {
        Report::generate(_p.inputs[0], o);
    }, o);
}

//template <typename Report> void report_2(typename Report::Options o = typename Report::Options())
//{
//    if (_p.inputs.size() != 2)
//    {
//        throw NotDoubleInputError();
//    }
//    
//    o.mix = mixture();
//    o.index = _p.opts[OPT_R_IND];
//
//    return startAnalysis<Report>([&](const typename Report::Options &o)
//    {
//        Report::generate(_p.inputs[0], _p.inputs[1], o);
//    }, o);
//}

template <typename Viewer> void viewer(typename Viewer::Options o = typename Viewer::Options())
{
    // Where the session files are generated
    o.path = _p.path;

    Viewer::generate(_p.opts.at(OPT_U_FILES), o);
}

/*
 * Template functions for analyzing
 */

template <typename Analyzer, typename Files> void analyze(const Files &files, typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(files, o);
    }, o);
}

// Analyze for a single sample with fuzzy matching
template <typename Analyzer> void analyzeFuzzy(typename Analyzer::Options o = typename Analyzer::Options())
{
    o.fuzzy = _p.fuzzy;
    return analyze<Analyzer>(_p.inputs[0], o);
}

// Analyze for a single sample
template <typename Analyzer> void analyze_1(Option x, typename Analyzer::Options o = typename Analyzer::Options())
{
    return analyze<Analyzer>(_p.opts.at(x), o);
}

// Analyze for two samples
template < typename Analyzer> void analyze_2(typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        if (_p.inputs.size() != 2)
        {
            throw NotDoubleInputError();
        }

        Analyzer::report(_p.inputs[0], _p.inputs[1], o);
    }, o);
}

// Analyze for n samples
template < typename Analyzer> void analyze_n(typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.inputs, o);
    }, o);
}

/*
 * Functions for parsing string to enums
 */

template <typename T> T parseEnum(const std::string &key, const std::string &str, const std::map<Value, T> &m)
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

template <typename T> T parseCSoft(const Value &str, const std::string &key)
{
    const static std::map<Value, T> m =
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

    // Required for unit-testing
    optind = 1;

    /*
     * Reconstruct the overall command
     */
    
    for (int i = 0; i < argc; i++)
    {
        _p.command += std::string(argv[i]) + " ";
    }

    assert(!_p.command.empty());

    int next, index;

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
    //opterr = 0;
    
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
            
            case OPT_REPORT:
            {
                if (val == "pdf")
                {
                    _p.pdf = true;
                }
                else
                {
                    throw InvalidValueException("-report", val);
                }

                break;
            }

            case OPT_FUZZY: { parseInt(val, _p.fuzzy); break; }

            case OPT_METHOD:
            {
                break;
            }

            /*
             * The following options can only be validated by the tool
             */

            case OPT_SOFT:
            case OPT_C_SOFT: { _p.opts[opt] = val; break; }

            case OPT_U_FILES:
            {
                std::vector<FileName> temp;
                Tokens::split(val, ",", temp);

                for (auto i = 0; i < temp.size(); i++)
                {
                    checkFile(_p.opts[opt] = temp[i]);
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
                
            case OPT_R_GENO:
            {
                checkFile(_p.opts[opt] = _p.rGeno = val);
                break;
            }

            case OPT_PATH: { _p.path = val; break; }

            default:
            {
                throw InvalidOptionException(argv[index]);
            }
        }
    }

    __anaquin__ = argv[0];
    __output__  = _p.path = checkPath(_p.path);

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
        case TOOL_T_KDIFF:
        case TOOL_T_ALIGN:
        case TOOL_T_EXPRESS:
        case TOOL_T_KEXPRESS:
        case TOOL_T_ASSEMBLY:
        case TOOL_T_COVERAGE:
        case TOOL_T_SUBSAMPLE:
        {
            std::cout << "[INFO]: Transcriptome Analysis" << std::endl;

            if (_p.tool != TOOL_T_IGV)
            {
                addRef(std::bind(&Standard::addTRef, &s, std::placeholders::_1));
                
                switch (_p.tool)
                {
                    case TOOL_T_DIFF:
                    case TOOL_T_KDIFF:
                    {
                        addMix(std::bind(&Standard::addTDMix, &s, std::placeholders::_1));
                        break;
                    }
                        
//                    case TOOL_T_REPORT:
//                    {
//                        if (_p.inputs.size() == 2)
//                        {
//                            addMix(std::bind(&Standard::addTMix, &s, std::placeholders::_1));
//                        }
//                        else
//                        {
//                            addMix(std::bind(&Standard::addTDMix, &s, std::placeholders::_1));
//                        }
//                        
//                        break;
//                    }

                    default:
                    {
                        addMix(std::bind(&Standard::addTMix, &s, std::placeholders::_1));
                        break;
                    }
                }

                s.r_trans.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_T_ALIGN:     { analyze_1<TAlign>(OPT_U_FILES);    break; }
                case TOOL_T_COVERAGE:  { analyze_1<TCoverage>(OPT_U_FILES); break; }
                case TOOL_T_ASSEMBLY:  { analyze_1<TAssembly>(OPT_U_FILES); break; }
                case TOOL_T_SUBSAMPLE: { analyze_1<TSample>(OPT_U_FILES);   break; }

                case TOOL_T_KDIFF:
                {
                    TKDiff::Options o;
                    o.index = _p.opts[OPT_R_IND];
                    analyze_1<TKDiff>(OPT_U_FILES, o);

                    break;
                }

                case TOOL_T_KEXPRESS:
                {
                    TKExpress::Options o;
                    o.index = _p.opts[OPT_R_IND];
                    analyze_2<TKExpress>(o);

                    break;
                }

                case TOOL_T_EXPRESS:
                {
                    auto checkCufflink = [&](const FileName &file)
                    {
                        const auto &r = Standard::instance().r_trans;
                        
                        Counts gs = 0;
                        Counts is = 0;
                        
                        ParserCufflink::parse(file, [&](const ParserCufflink::Data &data, const ParserProgress &p)
                        {
                            if (data.cID == ChrT)
                            {
                                try
                                {
                                    if (r.match(data.tID))
                                    {
                                        is++;
                                        
                                        // Important, if there's match for isoform, don't match for it's gene
                                        return;
                                    }
                                }
                                catch (...) {}
                                
                                try
                                {
                                    if (r.findGene(data.cID, data.id)) { gs++; }
                                }
                                catch (...) {}
                            }
                        });
                        
                        return gs > is;
                    };
                    
                    auto parseSoft = [&](const std::string &key, const std::string &str)
                    {
                        const static std::map<Value, TExpress::Software> m =
                        {
                            { "kallisto",  TExpress::Software::Kallisto  },
                            { "cufflink",  TExpress::Software::Cufflinks },
                            { "stringtie", TExpress::Software::StringTie },
                        };
                        
                        return parseEnum(key, str, m);
                    };
                    
                    TExpress::Options o;
                    
                    o.soft  = parseSoft("soft", _p.opts.at(OPT_SOFT));
                    o.metrs = TExpress::Metrics::Gene;
                    
                    if (o.soft == TExpress::Software::Cufflinks)
                    {
                        o.metrs = checkCufflink(_p.inputs[0]) ? TExpress::Metrics::Gene : TExpress::Metrics::Isoform;
                    }
                    
                    if (o.soft == TExpress::Software::Kallisto)
                    {
                        o.metrs = TExpress::Metrics::Isoform;
                    }

                    analyze_n<TExpress>(o);
                    break;
                }

                case TOOL_T_DIFF:
                {
                    auto parseSoft = [&](const std::string &str)
                    {
                        const static std::map<Value, TDiff::Software> m =
                        {
                            { "sleuth",   TDiff::Software::Sleuth   },
                            { "edgeR",    TDiff::Software::edgeR    },
                            { "deseq2",   TDiff::Software::DESeq2   },
                            { "cuffdiff", TDiff::Software::Cuffdiff },
                        };
                        
                        return parseEnum("soft", str, m);
                    };
                    
                    auto checkCuffdiff = [&](const FileName &file)
                    {
                        const auto &r = Standard::instance().r_trans;
                        
                        Counts gs = 0;
                        Counts is = 0;
                        
                        ParserCDiff::parse(file, [&](const ParserCDiff::Data &data, const ParserProgress &p)
                        {
                            if (data.cID == ChrT)
                            {
                                try
                                {
                                    if (r.match(data.id))
                                    {
                                        is++;
                                        
                                        // Important, if there's match for isoform, don't match for it's gene
                                        return;
                                    }
                                }
                                catch (...) {}
                                
                                try
                                {
                                    if (r.findGene(data.cID, data.id)) { gs++; }
                                }
                                catch (...) {}
                            }
                        });
                        
                        return gs > is;
                    };
                    
                    TDiff::Options o;

                    o.dSoft = parseSoft(_p.opts[OPT_SOFT]);
                    o.metrs = TDiff::Metrics::Gene;
                    
                    if (o.dSoft == TDiff::Software::Cuffdiff && !checkCuffdiff(_p.inputs[0]))
                    {
                        o.metrs = TDiff::Metrics::Isoform;
                    }
                    
                    if (o.dSoft == TDiff::Software::Sleuth)
                    {
                        o.metrs = TDiff::Metrics::Isoform;
                    }
                    
                    /*
                     * Optional count tables (eg: HTSeqCount)
                     */
                    
                    if (_p.opts.count(OPT_C_FILES))
                    {
                        o.cSoft  = parseCSoft<TDiff::Counting>(_p.opts[OPT_C_SOFT], "csoft");
                        o.counts = _p.oInputs;
                    }

                    analyze_1<TDiff>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_T_IGV: { viewer<TViewer>(); break; }
            }

            break;
        }

        case TOOL_S_DISCOVER:
        {
            std::cout << "[INFO]: Structual Analysis" << std::endl;
            
            switch (_p.tool)
            {
                case TOOL_S_DISCOVER:
                {
                    applyRef(std::bind(&Standard::addSStruct, &s, std::placeholders::_1), OPT_R_BED);
                    break;
                }

                default: { break; }
            }
            
            switch (_p.tool)
            {
                case TOOL_S_DISCOVER:
                {
                    break;
                }
                    
                default: { break; }
            }

            break;
        }

        case TOOL_F_IGV:
        case TOOL_F_DIFF:
        case TOOL_F_NORMAL:
        case TOOL_F_FUSION:
        case TOOL_F_DISCOVER:
        case TOOL_F_COVERAGE:
        {
            auto parseAligner = [&](const std::string &str)
            {
                const static std::map<Value, FusionCaller> m =
                {
                    { "Star"        , FusionCaller::StarFusion   },
                    { "StarFusion"  , FusionCaller::StarFusion   },
                    { "TopHat"      , FusionCaller::TopHatFusion },
                    { "TopHatFusion", FusionCaller::TopHatFusion },
                };

                return parseEnum("soft", str, m);
            };
            
            std::cout << "[INFO]: Fusion Analysis" << std::endl;

            switch (_p.tool)
            {
                case TOOL_F_NORMAL:
                {
                    addMix(std::bind(&Standard::addFMix,      &s, std::placeholders::_1));
                    applyRef(std::bind(&Standard::addFJunct, &s, std::placeholders::_1), OPT_R_BED);
                    break;
                }

                case TOOL_F_DIFF:
                {
                    addMix(std::bind(&Standard::addFMix,     &s, std::placeholders::_1));
                    applyRef(std::bind(&Standard::addFRef,   &s, std::placeholders::_1), OPT_R_BED);
                    applyRef(std::bind(&Standard::addFJunct, &s, std::placeholders::_1), OPT_R_FUS);
                    break;
                }

                case TOOL_F_COVERAGE:
                {
                    applyRef(std::bind(&Standard::addFStd, &s, std::placeholders::_1), OPT_R_BED);
                    break;
                }

                case TOOL_F_FUSION:
                {
                    applyRef(std::bind(&Standard::addFRef, &s, std::placeholders::_1));
                    addMix(std::bind(&Standard::addFMix,   &s, std::placeholders::_1));
                    break;
                }
                    
                case TOOL_F_DISCOVER:
                {
                    applyRef(std::bind(&Standard::addFRef, &s, std::placeholders::_1));
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
                case TOOL_F_IGV:      { viewer<FViewer>();                 break; }
                case TOOL_F_COVERAGE: { analyze_1<FCoverage>(OPT_U_FILES); break; }

                case TOOL_F_NORMAL:
                {
                    FNormal::Options o;
                    o.soft = parseAligner(_p.opts.at(OPT_SOFT));

                    analyze_1<FNormal>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_F_FUSION:
                {
                    FFusion::Options o;
                    o.soft = parseAligner(_p.opts.at(OPT_SOFT));
                    
                    analyze_1<FFusion>(OPT_U_FILES, o);
                    break;
                }
                    
                case TOOL_F_DIFF:
                {
                    FDiff::Options o;
                    o.soft = parseAligner(_p.opts.at(OPT_SOFT));

                    analyze_2<FDiff>();
                    break;
                }

                case TOOL_F_DISCOVER:
                {
                    FDiscover::Options o;
                    o.soft = parseAligner(_p.opts.at(OPT_SOFT));

                    analyzeFuzzy<FDiscover>(o);
                    break;
                }
            }

            break;
        }

        case TOOL_L_NORM:
        {
            std::cout << "[INFO]: Ladder Analysis" << std::endl;

            addMix(std::bind(&Standard::addLMix, &s, std::placeholders::_1));
            Standard::instance().r_lad.finalize();

            switch (_p.tool)
            {
                case TOOL_L_NORM:  { analyze_1<LNorm>(OPT_U_FILES); break; }
                case TOOL_L_DIFF:  { analyze_2<LDiff>(); break; }
            }

            break;
        }

        case TOOL_V_IGV:
        case TOOL_V_ALIGN:
        case TOOL_V_ALLELE:
        case TOOL_V_EXPRESS:
        case TOOL_V_KALLELE:
        case TOOL_V_DISCOVER:
        case TOOL_V_COVERAGE:
        case TOOL_V_KEXPRESS:
        case TOOL_V_SUBSAMPLE:
        {
            auto parseExpress = [&](const std::string &str)
            {
                const static std::map<Value, VExpress::Software> m =
                {
                    { "kallisto", VExpress::Software::Kallisto },
                };

                return parseEnum("soft", str, m);
            };
            
            std::cout << "[INFO]: Variant Analysis" << std::endl;

            if (_p.tool != TOOL_V_IGV)
            {
                switch (_p.tool)
                {
                    case TOOL_V_SUBSAMPLE:
                    {
                        applyRef(std::bind(&Standard::addStd,    &s, std::placeholders::_1), OPT_R_BED);
                        applyRef(std::bind(&Standard::addInters, &s, std::placeholders::_1), OPT_R_GENO);
                        break;
                    }

                    case TOOL_V_ALIGN:
                    case TOOL_V_COVERAGE:
                    {
                        applyRef(std::bind(&Standard::addStd, &s, std::placeholders::_1),    OPT_R_BED);
                        applyRef(std::bind(&Standard::addInters, &s, std::placeholders::_1), OPT_R_GENO);
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

                case TOOL_V_KALLELE:
                {
                    VKAllele::Options o;
                    o.index = _p.opts[OPT_R_IND];
                    analyze_2<VKAllele>(o);
                    break;
                }
                    
                case TOOL_V_KEXPRESS:
                {
                    VKExpress::Options o;
                    o.index = _p.opts[OPT_R_IND];
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
                    auto parse = [&](const std::string &str)
                    {
                        const static std::map<Value, VAllele::Software> m =
                        {
                            { "gatk"   ,  VAllele::Software::GATK     },
                            { "VarScan",  VAllele::Software::VarScan  },
                            { "VarScan2", VAllele::Software::VarScan  },
                            { "kallisto", VAllele::Software::Kallisto  },
                        };
                        
                        return parseEnum("soft", str, m);
                    };

                    VAllele::Options o;
                    o.soft = parse(_p.opts.at(OPT_SOFT));

                    analyze_1<VAllele>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_V_DISCOVER:
                {
                    auto parse = [&](const std::string &str)
                    {
                        const static std::map<Value, VDiscover::Software> m =
                        {
                            { "gatk"   ,  VDiscover::Software::GATK     },
                            { "VarScan",  VDiscover::Software::VarScan  },
                            { "VarScan2", VDiscover::Software::VarScan  },
                        };
                        
                        return parseEnum("soft", str, m);
                    };

                    VDiscover::Options o;
                    o.soft = parse(_p.opts.at(OPT_SOFT));

                    analyze_1<VDiscover>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_V_SUBSAMPLE: { analyze_1<VSample>(OPT_U_FILES); break; }
            }

            break;
        }

        case TOOL_M_IGV:
        case TOOL_M_DIFF:
        case TOOL_M_EXPRESS:
        case TOOL_M_ASSEMBLY:
        case TOOL_M_COVERAGE:
        {
            std::cout << "[INFO]: Metagenomics Analysis" << std::endl;

            if (_p.tool != TOOL_M_IGV)
            {
                switch (_p.tool)
                {
                    case TOOL_M_EXPRESS:                        
                    case TOOL_M_ASSEMBLY:
                    case TOOL_M_COVERAGE:
                    {
                        applyRef(std::bind(&Standard::addMRef, &s, std::placeholders::_1));
                        break;
                    }

                    default: { break; }
                }
                
                switch (_p.tool)
                {
                    case TOOL_M_ASSEMBLY:
                    case TOOL_M_EXPRESS:
                    case TOOL_M_DIFF:
                    {
                        addMix(std::bind(&Standard::addMMix, &s, std::placeholders::_1));
                        break;
                    }
                }

                Standard::instance().r_meta.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_M_IGV:      { viewer<FViewer>();                 break; }
                case TOOL_M_COVERAGE: { analyze_1<MCoverage>(OPT_U_FILES); break; }

                case TOOL_M_EXPRESS:
                {
                    auto parse = [&](const std::string &str)
                    {
                        const static std::map<Value, MExpress::Software> m =
                        {
                            { "bwa",     MExpress::BWA    },
                            { "bowtie",  MExpress::Bowtie },
                            { "bowtie2", MExpress::Bowtie },
                            { "velvet",  MExpress::Velvet  },
                            { "raymeta", MExpress::RayMeta },
                        };
                        
                        return parseEnum("soft", str, m);
                    };
                    
                    // Only defined for certain assemblers
                    FileName conts;
                    
                    const auto soft = parse(_p.opts.at(OPT_SOFT));

                    MExpress::Options o;
                    
                    o.soft    = soft;
                    o.contigs = conts;
                    
                    analyze_n<MExpress>(o);
                    break;
                }

                case TOOL_M_DIFF:
                {
                    auto parse = [&](const std::string &str)
                    {
                        const static std::map<Value, MDiff::Software> m =
                        {
                            { "stamp", MDiff::Software::STAMP },
                        };
                        
                        return parseEnum("soft", str, m);
                    };

                    MDiff::Options o;

                    o.soft = parse(_p.opts.at(OPT_SOFT));
                    
                    analyze_n<MDiff>(o);
                    break;
                }

                case TOOL_M_ASSEMBLY:
                {
                    auto parse = [&](const std::string &str)
                    {
                        const static std::map<Value, MAssembly::Software> m =
                        {
                            { "quast", MAssembly::MetaQuast },
                        };
                        
                        return parseEnum("soft", str, m);
                    };

                    const auto soft = parse(_p.opts.at(OPT_SOFT));
                    
                    MAssembly::Options o;
                    
                    o.soft = soft;
                    
                    switch (o.soft)
                    {
                        case MAssembly::MetaQuast:
                        {
                            const auto checkGenome = [&](const FileName &x)
                            {
                                Counts n = 0;
                                
                                ParserQuast::parseGenome(Reader(x),
                                                         [&](const ParserQuast::GenomeData &,
                                                             const ParserProgress &)
                                {
                                    n++;
                                });
                                
                                return n;
                            };
                            
                            const auto checkContigs = [&](const FileName &x)
                            {
                                Counts n = 0;
                                
                                ParserQuast::parseContigs(Reader(x),
                                                          [&](const ParserQuast::ContigData &,
                                                              const ParserProgress &)
                                {
                                    n++;
                                });
                                
                                return n;
                            };
                            
                            const auto sorted = sortInputs(_p.inputs[0],
                                                           _p.inputs[1],
                                                           _p.inputs[2],
                                                           checkGenome,
                                                           checkContigs);
                            o.genome  = sorted[0];
                            o.contigs = sorted[1];
                            
                            analyze<MAssembly>(sorted[2], o);
                            break;
                        }
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
    catch (const NotSingleInputError &ex)
    {
        printError("Invalid command. An input file is required.");
    }
    catch (const NotDoubleInputError &ex)
    {
        printError("Invalid command. Two input files are required.");
    }
    catch (const NoValueError &ex)
    {
        printError("Invalid command. Need to specify " + ex.opt + ".");
    }
    catch (const InvalidToolError &ex)
    {
        printError("Invalid command. Unknown tool: " + ex.val + ". Please check the user manual and try again.");
    }
    catch (const InvalidOptionException &ex)
    {
        printError((boost::format("Invalid command. Unknown option: %1%") % ex.opt).str());
    }
    catch (const InvalidValueException &ex)
    {
        printError((boost::format("Invalid command. %1% not expected for option -%2%. Please check and try again.") % ex.opt % ex.val).str());
    }
    catch (const MissingOptionError &ex)
    {
        const auto format = "Invalid command. Mandatory option is missing. Please specify %1%.";
        printError((boost::format(format) % ex.opt).str());
    }
    catch (const InvalidFileError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid command. File is invalid: " % ex.file).str());
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
