#include <ctime>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <strings.h>
#include <execinfo.h>
#include <sys/stat.h>

#include "RnaQuin/r_gene.hpp"
#include "RnaQuin/r_fold.hpp"
#include "RnaQuin/r_align.hpp"
#include "RnaQuin/r_sample.hpp"
#include "RnaQuin/r_report.hpp"
#include "RnaQuin/r_express.hpp"
#include "RnaQuin/r_assembly.hpp"

#include "VarQuin/v_wgs.hpp"
#include "VarQuin/v_trim.hpp"
#include "VarQuin/v_flip.hpp"
#include "VarQuin/v_align.hpp"
#include "VarQuin/v_sample.hpp"
#include "VarQuin/v_report.hpp"
#include "VarQuin/v_cancer.hpp"
#include "VarQuin/v_kexpress.hpp"

#include "MetaQuin/m_diff.hpp"
#include "MetaQuin/m_align.hpp"
#include "MetaQuin/m_abund.hpp"
#include "MetaQuin/m_sample.hpp"
#include "MetaQuin/m_report.hpp"
#include "MetaQuin/m_assembly.hpp"

#include "parsers/parser_gtf.hpp"
#include "parsers/parser_blat.hpp"
#include "parsers/parser_fold.hpp"
#include "parsers/parser_cdiff.hpp"
#include "parsers/parser_edgeR.hpp"
#include "parsers/parser_DESeq2.hpp"
#include "parsers/parser_sleuth.hpp"
#include "parsers/parser_express.hpp"
#include "parsers/parser_cufflink.hpp"
#include "parsers/parser_kallisto.hpp"

#include "writers/file_writer.hpp"
#include "writers/terminal_writer.hpp"

#ifdef UNIT_TEST
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#endif

typedef int Tool;
typedef int Option;

typedef std::string Value;
typedef std::set<Value> Range;

#define TOOL_TEST       264
#define TOOL_HELP       265
#define TOOL_R_ALIGN    266
#define TOOL_R_ASSEMBLY 267
#define TOOL_R_EXPRESS  268
#define TOOL_R_FOLD     269
#define TOOL_R_GENE     270
#define TOOL_R_REPORT   271
#define TOOL_R_SAMPLE   272
#define TOOL_V_KEXPRESS 273
#define TOOL_V_VREPORT  274
#define TOOL_V_ALIGN    275
#define TOOL_V_FLIP     276
#define TOOL_V_WGS      277
#define TOOL_V_CANCER   278
#define TOOL_V_TRIM     279
#define TOOL_V_SAMPLE   280
#define TOOL_M_ALIGN    281
#define TOOL_M_ABUND    282
#define TOOL_M_FOLD     283
#define TOOL_M_SAMPLE   284
#define TOOL_M_ASSEMBLY 285
#define TOOL_M_REPORT   286

/*
 * Options specified in the command line
 */

#define OPT_TEST    320
#define OPT_TOOL    321
#define OPT_PATH    325
#define OPT_VERSION 338

#define OPT_R_BASE  800
#define OPT_R_BED   801
#define OPT_METHOD  802
#define OPT_R_GTF   803
#define OPT_R_VCF   804
#define OPT_TRIM    805
#define OPT_MIXTURE 806
#define OPT_FUZZY   807
#define OPT_REPORT  808
#define OPT_R_IND   809
#define OPT_U_HG    810
#define OPT_U_SEQS  811
#define OPT_U_BAM   812
#define OPT_U_FILES 813
#define OPT_EDGE    814
#define OPT_U_BASE  817

using namespace Anaquin;

// Shared with other modules
bool __showInfo__ = true;

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
    
    strftime(buffer, 80, "%d-%m-%Y %H:%M:%S", timeinfo);
    std::string str(buffer);
    
    return str;
}

static std::map<Value, Tool> _tools =
{
    { "Test",           TOOL_TEST        },
    { "Help",           TOOL_HELP        },

    { "RnaAlign",       TOOL_R_ALIGN     },
    { "RnaAssembly",    TOOL_R_ASSEMBLY  },
    { "RnaExpress",     TOOL_R_EXPRESS   },
    { "RnaExpression",  TOOL_R_EXPRESS   },
    { "RnaReport",      TOOL_R_REPORT    },
    { "RnaFoldChange",  TOOL_R_FOLD      },
    { "RnaSubsample",   TOOL_R_SAMPLE },
    { "RnaGene",        TOOL_R_GENE      },

    { "VarAlign",       TOOL_V_ALIGN     },
    { "VarWGS",         TOOL_V_WGS       },
    { "VarCancer",      TOOL_V_CANCER    },
    { "VarVReport",     TOOL_V_VREPORT   },
    { "VarSubsample",   TOOL_V_SAMPLE    },
    { "VarTrim",        TOOL_V_TRIM      },
    { "VarFlip",        TOOL_V_FLIP      },
    { "VarKExpress",    TOOL_V_KEXPRESS  },

    { "MetaFoldChange", TOOL_M_FOLD      },
    { "MetaAlign",      TOOL_M_ALIGN     },
    { "MetaAbund",      TOOL_M_ABUND     },
    { "MetaAssembly",   TOOL_M_ASSEMBLY  },
    { "MetaSubsample",  TOOL_M_SAMPLE },
};

static std::map<Tool, std::set<Option>> _options =
{
    /*
     * RnaQuin Analysis
     */
    
    { TOOL_R_SAMPLE,   { OPT_U_FILES, OPT_METHOD } },
    { TOOL_R_ASSEMBLY, { OPT_R_GTF, OPT_MIXTURE, OPT_U_FILES } },
    { TOOL_R_FOLD,     { OPT_MIXTURE, OPT_U_FILES, OPT_METHOD } },
    { TOOL_R_EXPRESS,  { OPT_MIXTURE, OPT_U_FILES, OPT_METHOD } },
    { TOOL_R_REPORT,   { OPT_MIXTURE, OPT_R_IND, OPT_U_FILES  } },
    { TOOL_R_ALIGN,    { OPT_R_GTF, OPT_U_FILES } },

    /*
     * VarQuin Analysis
     */

    { TOOL_V_FLIP,     { OPT_U_BAM  } },
    { TOOL_V_TRIM,     { OPT_R_BED,   OPT_U_SEQS } },
    { TOOL_V_ALIGN,    { OPT_R_BED,   OPT_U_HG,   OPT_U_SEQS } },
    { TOOL_V_SAMPLE,   { OPT_R_BED,   OPT_U_HG,   OPT_U_SEQS, OPT_METHOD } },
    { TOOL_V_WGS, { OPT_R_VCF,   OPT_U_SEQS } },
    { TOOL_V_VREPORT,  { OPT_MIXTURE, OPT_U_SEQS } },
    { TOOL_V_KEXPRESS, { OPT_MIXTURE, OPT_U_SEQS } },

    /*
     * MetaQuin Analysis
     */

    { TOOL_M_FOLD,     { OPT_U_FILES } },
    { TOOL_M_ALIGN,    { OPT_R_BED, OPT_U_FILES } },
    { TOOL_M_ASSEMBLY, { OPT_R_BED, OPT_U_FILES } },
    { TOOL_M_SAMPLE,   { OPT_U_FILES, OPT_R_BED, OPT_METHOD } },
};

/*
 * Variables used in argument parsing
 */

struct Parsing
{
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

    // Generating a report?
    bool report = false;
    
    Proportion sampled = NAN;
    
    Tool tool = 0;
};

// Wrap the variables so that it'll be easier to reset them
static Parsing _p;

static FileName __mockGTFRef__;

void SetGTFRef(const FileName &file)
{
    __mockGTFRef__ = file;
}

FileName GTFRef()
{
    return !__mockGTFRef__.empty() ? __mockGTFRef__ : _p.opts.at(OPT_R_GTF);
}

FileName MixRef()
{
    return _p.opts.at(OPT_MIXTURE);
}

FileName VCFRef()
{
    return _p.opts.at(OPT_R_VCF);
}

FileName BedRef()
{
    return _p.opts.at(OPT_R_BED);
}

static Scripts fixManual(const Scripts &str)
{
    auto x = str;
    
    boost::replace_all(x, "<b>", "\e[1m");
    boost::replace_all(x, "</b>", "\e[0m");
    boost::replace_all(x, "<i>", "\e[3m");
    boost::replace_all(x, "</i>", "\e[0m");
    
    return x;
}

typedef std::string ErrorMsg;

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

struct UnknownFormatError : public std::runtime_error
{
    UnknownFormatError() : std::runtime_error("Unknown format") {}
};

/*
 * Argument options
 */

static const char *short_options = ":";

static const struct option long_options[] =
{
    { "v",       no_argument, 0, OPT_VERSION },
    { "version", no_argument, 0, OPT_VERSION },

    { "uhuman",  required_argument, 0, OPT_U_HG    },
    { "useqs",   required_argument, 0, OPT_U_SEQS  },
    { "ufiles",  required_argument, 0, OPT_U_FILES },

    { "m",       required_argument, 0, OPT_MIXTURE },
    { "mix",     required_argument, 0, OPT_MIXTURE },
    { "method",  required_argument, 0, OPT_METHOD  },
    { "trim",    required_argument, 0, OPT_TRIM    },

    { "ubam",    required_argument, 0, OPT_U_BAM  },
    { "rbed",    required_argument, 0, OPT_R_BED  },
    { "rgtf",    required_argument, 0, OPT_R_GTF  },
    { "rvcf",    required_argument, 0, OPT_R_VCF  },
    { "rind",    required_argument, 0, OPT_R_IND  },

    { "edge",    required_argument, 0, OPT_EDGE   },
    { "fuzzy",   required_argument, 0, OPT_FUZZY  },
    { "report",  required_argument, 0, OPT_REPORT },
    
    { "o",       required_argument, 0, OPT_PATH },

    {0, 0, 0, 0 }
};

// Used by modules without std::iostream defined
void printWarning(const std::string &msg)
{
    std::cout << "[Warn]: " << msg << std::endl;
}

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
    std::cout << fixManual(Manual()) << std::endl;
}

static Scripts manual(Tool tool)
{
    extern Scripts VarFlip();
    extern Scripts VarTrim();
    extern Scripts VarAlign();
    extern Scripts VarVReport();
    extern Scripts VarDiscover();
    extern Scripts VarSubsample();
    extern Scripts RnaAlign();
    extern Scripts RnaReport();
    extern Scripts RnaAssembly();
    extern Scripts RnaSubsample();
    extern Scripts RnaExpression();
    extern Scripts RnaFoldChange();
    extern Scripts MetaAlign();
    extern Scripts MetaAbund();
    extern Scripts MetaAssembly();
    extern Scripts MetaSubsample();
    extern Scripts MetaFoldChange();
    
    switch (tool)
    {
        case TOOL_R_ALIGN:     { return RnaAlign();       }
        case TOOL_R_ASSEMBLY:  { return RnaAssembly();    }
        case TOOL_R_EXPRESS:   { return RnaExpression();  }
        case TOOL_R_REPORT:    { return RnaReport();      }
        case TOOL_R_FOLD:      { return RnaFoldChange();  }
        case TOOL_R_SAMPLE:    { return RnaSubsample();   }
        case TOOL_V_FLIP:      { return VarFlip();        }
        case TOOL_V_TRIM:      { return VarTrim();        }
        case TOOL_V_ALIGN:     { return VarAlign();       }
        case TOOL_V_SAMPLE:    { return VarSubsample();   }
        case TOOL_V_WGS:  { return VarDiscover();    }
        case TOOL_V_VREPORT:   { return VarVReport();     }
        case TOOL_M_ALIGN:     { return MetaAlign();      }
        case TOOL_M_SAMPLE:    { return MetaAbund();      }
        case TOOL_M_ASSEMBLY:  { return MetaAssembly();   }
        case TOOL_M_ABUND:     { return MetaSubsample();  }
        case TOOL_M_FOLD:      { return MetaFoldChange(); }
    }

    throw std::runtime_error("Manual not found");
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

        std::vector<std::string> toks;
        Tokens::split(l, "\t", toks);

        std::cout << toks[0] << "\t" << toks[2] << "\t" << toks[3] << std::endl;
    }
}

FileName mixture()
{
    return _p.opts[OPT_MIXTURE];
}

#define CHECK_REF(x) (x != OPT_MIXTURE && x > OPT_R_BASE && x < OPT_U_BASE)

static void printError(const std::string &msg)
{
    std::cerr << "***********************" << std::endl;
    std::cerr << "[ERRO]: " << msg << std::endl;
    std::cerr << "***********************" << std::endl;
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

template <typename Reference> void addRef(const ChrID &cID, Reference ref, const FileName &file)
{
    if (__showInfo__)
    {
        std::cout << "[INFO]: Loading: " << file << std::endl;
    }
    
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

                default:
                {
                    addRef(ChrIS, ref, _p.opts[opt]);
                    break;
                }
            }
        }
    }
}

// Apply a reference source given where it comes from
template <typename Reference> void applyRef(Reference ref, Option opt)
{
    if (__showInfo__)
    {
        std::cout << "[INFO]: Loading: " << _p.opts[opt] << std::endl;
    }
    
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
    
    system(("mkdir -p " + path).c_str());
    
    o.work  = path;
    
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);

    o.info(_p.command);
    o.info(date());
    o.info("Path: " + path);

    std::clock_t begin = std::clock();

    f(o);
    
    std::clock_t end = std::clock();

    const auto elapsed = (boost::format("Completed. %1% seconds.") % (double(end - begin) / CLOCKS_PER_SEC)).str();
    o.info(elapsed);

#ifndef DEBUG
    o.logger->close();
#endif
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

template <typename Analyzer> void analyze_0(typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(o);
    }, o);
}

template <typename Analyzer> void analyze_1(Option x, typename Analyzer::Options o = typename Analyzer::Options())
{
    return analyze<Analyzer>(_p.opts.at(x), o);
}

template <typename Analyzer> void analyze_2(Option x1, Option x2, typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.opts[x1], _p.opts[x2], o);
    }, o);
}

template <typename Analyzer> void analyze_3(typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        A_ASSERT(_p.inputs.size() == 3);
        Analyzer::report(_p.inputs[0], _p.inputs[1], _p.inputs[2], o);
    }, o);
}

template < typename Analyzer> void analyze_n(typename Analyzer::Options o = typename Analyzer::Options())
{
    return startAnalysis<Analyzer>([&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.inputs, o);
    }, o);
}

static void fixInputs(int argc, char ** argv)
{
    for (auto i = 0; i < argc; i++)
    {
        if (argv[i] && strlen(argv[i]) && argv[i][0] == (char)-30)
        {
            auto tmp = std::string(argv[i]);
            
            auto invalid = [&](char c)
            {
                return !(c>=0 && c<128);
            };
            
            tmp.erase(remove_if(tmp.begin(), tmp.end(), invalid), tmp.end());
            tmp = "-" + tmp;

            strcpy(argv[i], tmp.c_str());
        }
    }
}

void parse(int argc, char ** argv)
{
    auto tmp = new char*[argc+1];
    
    for (auto i = 0; i < argc; i++)
    {
        tmp[i] = (char *) malloc((strlen(argv[i]) + 1) * sizeof(char));
        strcpy(tmp[i], argv[i]);
    }
    
    fixInputs(argc, argv=tmp);

    auto &tool = _p.tool;
    
    _p = Parsing();

    if (argc <= 1)
    {
        printUsage();
        return;
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

    A_ASSERT(!_p.command.empty());

    int next, index;

    // Attempt to parse and store a floating point from string
    auto parseDouble = [&](const std::string &str, double &r)
    {
        A_ASSERT(next);
        
        try
        {
            r = stof(str);
        }
        catch (...)
        {
            throw std::runtime_error(str + " is not a valid floating number. Please check and try again.");
        }
    };

    // Attempt to parse and store an integer from string
    auto parseInt = [&](const std::string &str, unsigned &r)
    {
        A_ASSERT(next);
        
        try
        {
            r = stoi(str);
        }
        catch (...)
        {
            throw std::runtime_error(str + " is not a valid integer. Please check and try again.");
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

    if (!_tools.count(argv[1]) && strcmp(argv[1], "-v"))
    {
        throw InvalidToolError(argv[1]);
    }
    else
    {
        _p.tool = _tools[argv[1]];
    }
    
    if (_p.tool == TOOL_V_SAMPLE || _p.tool == TOOL_R_SAMPLE || _p.tool == TOOL_V_TRIM)
    {
        __showInfo__ = false;
    }
    
    const auto isHelp = argc >= 3 && (!strcmp(argv[2], "-h") || !strcmp(argv[2], "--help"));

    if (isHelp)
    {
        if (argc != 3)
        {
            throw std::runtime_error("Too many arguments for help usage. Usage: anaquin <tool> -h or anaquin <tool> --help");
        }

        std::cout << fixManual(manual(_p.tool)) << std::endl << std::endl;
        return;
    }
    
    A_ASSERT(_p.tool);

    unsigned n = 2;

    while ((next = getopt_long_only(argc, argv, short_options, long_options, &index)) != -1)
    {
        if (next == ':')
        {
            throw NoValueError(argv[n]);
        }
        else if (next < OPT_TOOL)
        {
            throw InvalidOptionException(argv[n]);
        }
        
        opts.push_back(next);
        
        // Whether this option has an value
        const auto hasValue = optarg;
        
        n += hasValue ? 2 : 1;
        
        vals.push_back(hasValue ? std::string(optarg) : "");
    }

    for (auto i = 0; i < opts.size(); i++)
    {
        auto opt = opts[i];
        auto val = vals[i];

        switch (opt)
        {
            case OPT_REPORT: { _p.report = true; break; }

            case OPT_EDGE:
            case OPT_FUZZY:
            {
                try
                {
                    stoi(val);
                    _p.opts[opt] = val;
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not an integer. Please check and try again.");
                }

                break;
            }

            case OPT_METHOD:
            {
                switch (_p.tool)
                {
                    case TOOL_R_FOLD:
                    case TOOL_R_EXPRESS:
                    case TOOL_V_SAMPLE: { _p.opts[opt] = val; break; }

                    case TOOL_M_SAMPLE:
                    case TOOL_R_SAMPLE:
                    {
                        parseDouble(_p.opts[opt] = val, _p.sampled);
                        
                        if (_p.sampled <= 0.0)
                        {
                            throw std::runtime_error("Invalid value for -method. Sampling fraction must be greater than zero.");
                        }
                        else if (_p.sampled >= 1.0)
                        {
                            throw std::runtime_error("Invalid value for -method. Sampling fraction must be less than one.");
                        }
                        
                        break;
                    }
                }
                
                break;
            }

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

            case OPT_TRIM: { _p.opts[opt] = val; break; }

            case OPT_U_HG:
            case OPT_U_BAM:
            case OPT_R_IND:
            case OPT_R_VCF:
            case OPT_R_BED:
            case OPT_R_GTF:
            case OPT_U_SEQS:
            case OPT_MIXTURE:
            {
                checkFile(_p.opts[opt] = val); break;
            }

            case OPT_PATH: { _p.path = val; break; }

            default:
            {
                throw InvalidOptionException(argv[n]);
            }
        }
    }

    __output__ = _p.path = checkPath(_p.path);
    
    /*
     * Have all the required options given?
     */
    
    std::set<Option> required;
    
    std::copy_if(_options[_p.tool].begin(), _options[_p.tool].end(), std::inserter(required, required.end()), [&](const Option &x)
    {
        for (auto &o : long_options)
        {
            if (o.val == x)
            {
                return o.has_arg == required_argument;
            }
        }

        A_THROW("Option " + std::to_string(x) + " no found");
    });

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

    if (__showInfo__)
    {
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "------------- Sequin Analysis -----------" << std::endl;
        std::cout << "-----------------------------------------" << std::endl << std::endl;        
    }

    auto &s = Standard::instance();
    
    switch (_p.tool)
    {
        case TOOL_TEST:
        {
#ifdef UNIT_TEST
            Catch::Session().run(1, argv);
#else
            A_THROW("UNIT_TEST is not defined");
#endif
            break;
        }

        case TOOL_R_FOLD:
        case TOOL_R_GENE:
        case TOOL_R_ALIGN:
        case TOOL_R_REPORT:
        case TOOL_R_EXPRESS:
        case TOOL_R_ASSEMBLY:
        case TOOL_R_SAMPLE:
        {
            if (__showInfo__)
            {
                std::cout << "[INFO]: RNA-Seq Analysis" << std::endl;
            }

            if (_p.tool != TOOL_R_SAMPLE)
            {
                switch (_p.tool)
                {
                    case TOOL_R_GENE:
                    {
                        try
                        {
                            applyMix(std::bind(&Standard::addRMix, &s, std::placeholders::_1));
                        }
                        catch (...)
                        {
                            applyMix(std::bind(&Standard::addRDMix, &s, std::placeholders::_1));
                        }

                        break;
                    }

                    case TOOL_R_ALIGN:
                    {
                        addRef(std::bind(&Standard::addRRef, &s, std::placeholders::_1));
                        break;
                    }

                    case TOOL_R_FOLD:
                    {
                        applyMix(std::bind(&Standard::addRDMix, &s, std::placeholders::_1));
                        break;
                    }

                    case TOOL_R_EXPRESS:
                    case TOOL_R_ASSEMBLY:
                    {
                        applyMix(std::bind(&Standard::addRMix, &s, std::placeholders::_1));
                        break;
                    }

                    case TOOL_R_REPORT:
                    {
                        applyMix(std::bind(&Standard::addRDMix, &s, std::placeholders::_1));
                        break;
                    }

                    default:
                    {
                        addRef(std::bind(&Standard::addRRef, &s, std::placeholders::_1));
                        applyMix(std::bind(&Standard::addRMix, &s, std::placeholders::_1));
                        break;
                    }
                }

                s.r_rna.finalize();
            }

            switch (_p.tool)
            {
                case TOOL_R_GENE:     { analyze_0<RGene>();                         break; }
                case TOOL_R_ALIGN:    { analyze_1<RAlign>(OPT_U_FILES);             break; }
                case TOOL_R_ASSEMBLY: { analyze_1<RAssembly>(OPT_U_FILES);          break; }
                case TOOL_R_REPORT:
                {
                    RReport::Options o;
                    o.index = _p.opts[OPT_R_IND];
                    analyze_1<RReport>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_R_SAMPLE:
                {
                    RSample::Options o;
                    o.p = _p.sampled;
                    analyze_1<RSample>(OPT_U_FILES, o);
                    break;
                }
                    
                case TOOL_R_EXPRESS:
                {
                    RExpress::Options o;
                    
                    if (_p.opts[OPT_METHOD] != "gene" && _p.opts[OPT_METHOD] != "isoform")
                    {
                        throw InvalidValueException("-method", _p.opts[OPT_METHOD]);
                    }

                    o.metrs = _p.opts[OPT_METHOD] == "gene" ? RExpress::Metrics::Gene : RExpress::Metrics::Isoform;
                    
                    const auto &file = _p.inputs[0];
                    
                    // Is this a GTF by extension?
                    const auto isGTF = file.find(".gtf") != std::string::npos;
                    
                    if (isGTF)
                    {
                        o.format = RExpress::Format::GTF;
                    }
                    else if (ParserExpress::isExpress(file))
                    {
                        o.format = RExpress::Format::Text;
                    }
                    else if (ParserKallisto::isKallisto(file))
                    {
                        o.format = RExpress::Format::Kallisto;
                    }
                    else
                    {
                        A_THROW("Unknown file type: " + file + ". Input file should be a GTF file or Anaquin format. Please check our user guide for details.");
                    }

                    analyze_n<RExpress>(o);
                    break;
                }

                case TOOL_R_FOLD:
                {
                    RFold::Options o;

                    if (_p.opts[OPT_METHOD] != "gene" && _p.opts[OPT_METHOD] != "isoform")
                    {
                        throw InvalidValueException("-method", _p.opts[OPT_METHOD]);
                    }
                    
                    o.metrs = _p.opts[OPT_METHOD] == "gene" ? RFold::Metrics::Gene : RFold::Metrics::Isoform;
                    
                    const auto &file = _p.inputs[0];
                    
                    if (ParserCDiff::isTracking(Reader(file)))
                    {
                        o.format = RFold::Format::Cuffdiff;
                        std::cout << "[INFO]: Cuffdiff format" << std::endl;
                    }
                    else if (ParserSleuth::isSleuth(Reader(file)))
                    {
                        o.format = RFold::Format::Sleuth;
                        std::cout << "[INFO]: Sleuth format" << std::endl;
                    }
                    else if (ParserDESeq2::isDESeq2(Reader(file)))
                    {
                        o.format = RFold::Format::DESeq2;
                        std::cout << "[INFO]: DESeq2 format" << std::endl;
                    }
                    else if (ParserEdgeR::isEdgeR(Reader(file)))
                    {
                        o.format = RFold::Format::edgeR;
                        std::cout << "[INFO]: edgeR format" << std::endl;
                    }
                    else if (ParserDiff::isDiff(Reader(file)))
                    {
                        o.format = RFold::Format::Anaquin;
                        std::cout << "[INFO]: Anaquin format" << std::endl;
                    }
                    else
                    {
                        throw std::runtime_error("Unknown file format: " + file + ". Anaquin supports Cuffdiff, DESeq2, edgeR and RnaQuin FoldChange format. Please note the input file requires a header.");
                    }

                    analyze_1<RFold>(OPT_U_FILES, o);
                    break;
                }
            }

            break;
        }

        case TOOL_M_FOLD:
        case TOOL_M_ABUND:
        case TOOL_M_ALIGN:
        case TOOL_M_ASSEMBLY:
        case TOOL_M_SAMPLE:
        {
            switch (_p.tool)
            {
                case TOOL_M_FOLD:   { applyMix(std::bind(&Standard::addMDMix, &s, std::placeholders::_1)); break; }
                case TOOL_M_ABUND:  { applyMix(std::bind(&Standard::addMMix,  &s, std::placeholders::_1)); break; }
                case TOOL_M_ALIGN:  { applyRef(std::bind(&Standard::addMBed,  &s, std::placeholders::_1)); break; }
                case TOOL_M_SAMPLE: { applyRef(std::bind(&Standard::addMBed,  &s, std::placeholders::_1)); break; }

                case TOOL_M_ASSEMBLY:
                {
                    applyMix(std::bind(&Standard::addMMix, &s, std::placeholders::_1));
                    applyRef(std::bind(&Standard::addMBed, &s, std::placeholders::_1));
                    break;
                }

                default: { break; }
            }
            
            Standard::instance().r_meta.finalize();
            
            switch (_p.tool)
            {
                case TOOL_M_FOLD:
                {
//                    MDiff::Options o;
//                    
//                    if (_p.inputs.size() == 2 && ParserSAM::isBAM(Reader(_p.inputs[0]))
//                                              && ParserSAM::isBAM(Reader(_p.inputs[1])))
//                    {
//                        o.format = MDiff::Format::BAM;
//                    }
//                    else if (_p.inputs.size() == 4)
//                    {
//                        // TODO: Please fix this
//                        o.format = MDiff::Format::RayMeta;
//                    }
//                    else
//                    {
//                        throw UnknownFormatError();
//                    }
//
//                    analyze_n<MDiff>(o);
                    break;
                }

                case TOOL_M_ABUND:
                {
//                    MAbund::Options o;
//                    
//                    if (_p.inputs.size() == 1 && ParserSAM::isBAM(Reader(_p.inputs[0])))
//                    {
//                        o.format = MAbund::Format::BAM;
//                    }
//                    else if (_p.inputs.size() == 2 && ParserBlat::isBlat(Reader(_p.inputs[1])))
//                    {
//                        o.format = MAbund::Format::RayMeta;
//                    }
//                    else
//                    {
//                        throw UnknownFormatError();
//                    }
//                    
//                    analyze_n<MAbund>(o);
                    break;
                }

                case TOOL_M_SAMPLE:
                {
                    MSample::Options o;
                    o.p = _p.sampled;
                    analyze_1<MSample>(OPT_U_FILES, o);
                    break;
                }

                case TOOL_M_ASSEMBLY:
                {
                    MAssembly::Options o;

                    // This is the only supporting format
                    o.format = MAssembly::Format::Blat;

                    analyze_n<MAssembly>(o);
                    break;
                }

                case TOOL_M_ALIGN: { analyze_n<MAlign>(); break; }

                default: { break; }
            }

            break;
        }

        case TOOL_V_WGS:
        case TOOL_V_FLIP:
        case TOOL_V_TRIM:
        case TOOL_V_ALIGN:
        case TOOL_V_CANCER:
        case TOOL_V_SAMPLE:
        case TOOL_V_VREPORT:
        case TOOL_V_KEXPRESS:
        {
            if (__showInfo__)
            {
                std::cout << "[INFO]: Variant Analysis" << std::endl;
            }

            if (_p.tool != TOOL_V_FLIP)
            {
                switch (_p.tool)
                {
                    case TOOL_V_ALIGN:
                    {
                        applyRef(std::bind(&Standard::addVRef, &s, std::placeholders::_1, 0), OPT_R_BED);
                        break;
                    }
                        
                    case TOOL_V_SAMPLE:
                    {
                        applyRef(std::bind(&Standard::addVRef, &s, std::placeholders::_1,
                                            _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : 0), OPT_R_BED);
                        break;
                    }

                    case TOOL_V_TRIM:
                    {
                        applyRef(std::bind(&Standard::addVRef, &s, std::placeholders::_1, 0), OPT_R_BED);
                        break;
                    }

                    case TOOL_V_WGS:
                    case TOOL_V_CANCER:
                    {
                        applyRef(std::bind(&Standard::addVRef, &s, std::placeholders::_1, 0), OPT_R_BED);
                        applyRef(std::bind(&Standard::addVVar, &s, std::placeholders::_1), OPT_R_VCF);
                        break;
                    }
                        
                    case TOOL_V_KEXPRESS:
                    {
                        applyMix(std::bind(&Standard::addVMix, &s, std::placeholders::_1));
                        break;
                    }

                    default: { break; }
                }

                if (_p.tool != TOOL_V_VREPORT)
                {
                    Standard::instance().r_var.finalize();
                }
            }

            switch (_p.tool)
            {
                case TOOL_V_VREPORT:
                {
                    VVReport::Options o;
                    
                    o.rMix = mixture();
                    o.rBED = _p.opts.at(OPT_R_BED);
                    o.rVCF = _p.opts.at(OPT_R_VCF);
                    
                    analyze_3<VVReport>(o);
                    break;
                }
                    
                case TOOL_V_KEXPRESS: { analyze_1<VKExpress>(OPT_U_SEQS); break; }

                case TOOL_V_FLIP: { analyze_1<VFlip>(OPT_U_BAM); break; }
                
                case TOOL_V_ALIGN:
                {
                    analyze_2<VAlign>(OPT_U_HG, OPT_U_SEQS);
                    break;
                }

                case TOOL_V_TRIM:
                {
                    VTrim::Options o;
                    
                    if (_p.opts.count(OPT_METHOD))
                    {
                        const auto &x = _p.opts.at(OPT_METHOD);
                        
                        if (x == "leftRight")  { o.meth = VTrim::Method::LeftRight; }
                        else if (x == "left")  { o.meth = VTrim::Method::Left;      }
                        else if (x == "right") { o.meth = VTrim::Method::Right;     }
                        else
                        {
                            throw InvalidValueException("-method", x);
                        }
                    }

                    analyze_1<VTrim>(OPT_U_SEQS, o);
                    break;
                }

                case TOOL_V_WGS:
                {
                    VWGS::Options o;
                    analyze_1<VWGS>(OPT_U_SEQS, o);
                    break;
                }

                case TOOL_V_CANCER:
                {
                    VCancer::Options o;
                    analyze_1<VCancer>(OPT_U_SEQS, o);
                    break;
                }

                case TOOL_V_SAMPLE:
                {
                    VSample::Options o;
                    
                    // Eg: "mean", "median", "reads", "0.75"
                    const auto meth = _p.opts[OPT_METHOD];
                    
                    auto isFloat = [&]()
                    {
                        std::istringstream iss(meth);
                        float f;
                        iss >> std::noskipws >> f;
                        return iss.eof() && !iss.fail();
                    };
                    
                    if (meth == "mean")
                    {
                        o.meth = VSample::Method::Mean;
                    }
                    else if (meth == "median")
                    {
                        o.meth = VSample::Method::Median;
                    }
                    else if (meth == "reads")
                    {
                        o.meth = VSample::Method::Reads;
                    }
                    else if (isFloat())
                    {
                        o.p = stod(meth);
                        o.meth = VSample::Method::Prop;
                        
                        if (o.p <= 0.0)
                        {
                            throw std::runtime_error("Normalization factor must be greater than zero");
                        }
                        else if (o.p >= 1.0)
                        {
                            throw std::runtime_error("Normalization factor must be less than one");
                        }
                    }
                    else
                    {
                        throw std::runtime_error("Unknown method: " + meth);
                    }

                    if (_p.opts.count(OPT_EDGE))
                    {
                        o.edge = stoi(_p.opts[OPT_EDGE]);
                    }

                    analyze_2<VSample>(OPT_U_HG, OPT_U_SEQS, o);
                    break;
                }
            }

            break;
        }
    }
}

extern int parse_options(int argc, char ** argv)
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
    catch (const FailedCommandException &ex)
    {
        printError("Invalid command: " + std::string(ex.what()));
    }
    catch (const InvalidFormatException &ex)
    {
        printError("Invalid file format: " + std::string(ex.what()));
    }
    catch (const NoValueError &ex)
    {
        printError("Invalid command. Need to specify " + ex.opt + ".");
    }
    catch (const UnknownFormatError &ex)
    {
        printError("Unknown format for the input file(s)");
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
        printError((boost::format("Invalid command. %1% not expected for -%2%.") % ex.val % ex.opt).str());
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
