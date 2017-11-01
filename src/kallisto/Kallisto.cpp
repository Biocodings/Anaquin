#include <map>
#include <fstream>
#include <assert.h>
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"
#include "tools/errors.hpp"
#include "tools/system.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern std::string RefKKmers();

// Index for all sequin k-mers
static std::shared_ptr<KmerIndex> __allIndex__;

// Running statistics for k-mers
static KMStats __kStats__;

static void LoadRefKKmers()
{
    std::istringstream r(RefKKmers());
    std::string line;
    
    while (std::getline(r, line))
    {
        std::vector<std::string> toks;
        split(line, "\t", toks);
        
        if (toks.size() != 4)
        {
            continue;
        }

        A_ASSERT(toks[3] == "Span");
        
        KMPair p;
        p.norm = toks[1]; // Normal
        p.rcom = toks[2]; // Reverse complement
        
        // Eg: CS_011_R
        const auto sID = noLast(toks[0], "_");
        
        if (isSubstr(toks[0], "_R"))
        {
            __kStats__.vars[sID].R.push_back(p);
        }
        else
        {
            __kStats__.vars[sID].V.push_back(p);
        }
        
        if (toks[3] == "Span")
        {
            __kStats__.spans[p.norm] = 0;
            __kStats__.spans[p.rcom] = 0;
        }
    }
    
    A_ASSERT(!__kStats__.vars.empty());
    A_ASSERT(!__kStats__.spans.empty());
}

void KMInit(const std::string &aIndex, int k)
{
    ProgramOptions o;
    o.k = k;
    o.index = aIndex;
    
    /*
     * Initalize for reference on all sequin k-mers (e.g. sequins.fa.index)
     */
    
    __allIndex__ = std::shared_ptr<KmerIndex>(new KmerIndex(o));
    __allIndex__->load(o);
    
    /*
     * Reference k-mers spanning variants
     */
    
    LoadRefKKmers();
}

static void KMCount(const char *s)
{
    KmerIterator iter(s), end;
    
    // Number of k-mers that are sequins
    unsigned isSeq = 0;
    
    // Number of k-mers that are genome
    unsigned isGen = 0;
    
    for (int i = 0; iter != end; ++i, ++iter)
    {
        const auto k = iter->first.rep().toString();
        
        /*
         * Does the k-mer span sequin variants?
         */
        
        if (__kStats__.spans.count(k))
        {
            __kStats__.spans[k]++;
        }
        
        /*
         * Is this one of the reference k-mers?
         */

        std::vector<std::pair<KmerEntry, int> > v;
        __allIndex__->match(k.c_str(), 31, v);

        if (!v.empty())
        {
            // There's a match, let's update the counting table
            __kStats__.k2c[k]++;
            
            isSeq++;
        }
        else
        {
            isGen++;
        }
    }
    
    if (isSeq > isGen)
    {
        __kStats__.nSeq++;
    }
    else
    {
        __kStats__.nGen++;
    }
}

void KMCount(const char *s1, const char *s2)
{
    KMCount(s1);
    KMCount(s2);
}

SequinID Anaquin::KKM2Sequin(const Kmer &s, unsigned k)
{
    std::vector<std::pair<KmerEntry, int>> v;

    __allIndex__->match(s.c_str(), k, v);
    A_ASSERT(!v.empty());
    
    const Contig &c = __allIndex__->dbGraph.contigs[v[0].first.contig];
    A_ASSERT(!c.transcripts.empty());

    return __allIndex__->target_names_[c.transcripts[0].trid];
}

FileName Anaquin::KHumanFASTA(const FileName &file)
{
    const auto rfa = System::tmpFile();
    
    std::fstream r, w;
    r.open(file, std::fstream::in);
    w.open(rfa,  std::fstream::out);

    int i = -1;
    char buf[30000];
    
    // Reading "> ..."?
    bool isName = false;
    
    auto reset = [&]() { i = -1; };
    auto write = [&]() { w << buf[i];};

    while (((buf[++i] = r.get()) != EOF))
    {
        switch (buf[i])
        {
            case '\n': { if (isName) { write(); reset(); } isName = false; break; }
            case '>':
            {
                if (i)
                {
                    buf[i+1] = NULL;
                    w << (char *) buf;
                }

                isName = true;
                write();
                reset();
                
                break;
            }

            default: { if (isName) { write(); } break; }
        }
    }
    
    r.close();
    w.close();
    
    return rfa;
}

FileName Anaquin::KBuildIndex(const FileName &file, unsigned k)
{
    ProgramOptions opt;

    opt.k = ::Kmer::k = k;
    opt.index = System::tmpFile();
    opt.transfasta.push_back(file);
    
    KmerIndex index(opt);
    index.BuildTranscripts(opt);
    index.write(opt.index);
    
    if (index.dbGraph.contigs.empty() || !index.kmap.size())
    {
        std::runtime_error("Failed to build index for " + file);
    }
    
    return opt.index;
}

bool Anaquin::KQuerySeqs(const Sequence &s, unsigned k)
{
    std::vector<std::pair<KmerEntry, int>> v;
    __allIndex__->match(s.c_str(), k, v);
    return !v.empty();
}

bool Anaquin::KQuery___(const FileName &file, const Sequence &s)
{
    ProgramOptions opt;
    
    opt.k = ::Kmer::k = s.size();
    opt.index = file;

    KmerIndex index(opt);
    index.load(opt);
    
    std::vector<std::pair<KmerEntry, int>> v;
    index.match(s.c_str(), opt.k, v);

    return !v.empty();
}

KMStats Anaquin::KCount(const FileName &i1, const FileName &i2, const FileName &p1, const FileName &p2, unsigned k)
{
    // Initalize for reference k-mers
    KMInit(i1, k);

    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = i1;
    opt.files.push_back(p1);
    opt.files.push_back(p2);

    KmerIndex index(opt);
    MinCollector collection(index, opt);
    
    // Process k-mers
    ProcessReads(index, opt, collection);
    
    // Have we read anything?
    A_ASSERT((__kStats__.nGen + __kStats__.nSeq) > 0);
    
    return __kStats__;
}
