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
        p.normal  = toks[1]; // Normal
        p.revComp = toks[2]; // Reverse complement
        
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
            __kStats__.spans[p.normal]  = 0;
            __kStats__.spans[p.revComp] = 0;
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
#ifdef DEBUG
            __kStats__.all[k]++;
#endif
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

FileName Anaquin::KBuildIndex(const FileName &file, unsigned k)
{
    ProgramOptions opt;

    opt.k = ::Kmer::k = k;
    opt.index = System::tmpFile();
    opt.transfasta.push_back(file);
    
    KmerIndex index(opt);
    index.BuildTranscripts(opt);
    index.write(opt.index);
    
    return opt.index;
}

bool Anaquin::KQuery(const FileName &file, const Sequence &s)
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

KMStats Anaquin::KCount(const FileName &aIndex, const FileName &p1, const FileName &p2, unsigned k)
{
    // Initalize for reference k-mers
    KMInit(aIndex, k);

    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = aIndex;
    opt.files.push_back(p1);
    opt.files.push_back(p2);

    KmerIndex index(opt);
    MinCollector collection(index, opt);
    
    // Process k-mers
    ProcessReads(index, opt, collection);
    
    // Have we read anything?
    A_ASSERT((__kStats__.nGen + __kStats__.nSeq) > 0);
    
#ifdef DEBUG
    std::ofstream w1("KMAll_1.txt");
    
    for (const auto &i : __kStats__.all)
    {
        w1 << i.first << "\t" << i.second << std::endl;
    }
    
    w1.close();
    
    std::ofstream w2("KMAll_2.txt");
    
    for (const auto &i : __kStats__.spans)
    {
        w2 << i.first << "\t" << i.second << std::endl;
    }
    
    w2.close();
#endif
    
    return __kStats__;
}
