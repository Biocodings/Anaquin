#include <map>
#include <fstream>
#include <assert.h>
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"

using namespace Anaquin;

// Index for all sequin k-mers
static std::shared_ptr<KmerIndex> __allIndex__;

// Index for spanning variants
static std::map<std::string, unsigned> __span__;

static KMStats __kStats__;

/*
 * Initlaize k-mers spanning variants, they are useful for estimating allele frequency
 */

static void KMCancerSpan()
{
    std::ifstream r("CancerKMSpan.txt");
    
    if (!r.good())
    {
        throw std::runtime_error("Invalid CancerKMSpan.txt");
    }
    
    std::string line;
    while (std::getline(r, line))
    {
        std::vector<std::string> toks;
        split(line, "t", toks);
        assert(!toks.empty());
        
        if (toks[0] == "Name") { continue; }
        
        __span__[toks[1]] = 0; // Normal
        __span__[toks[2]] = 0; // Reverse complement
    }
    
    r.close();
    assert(!__span__.empty());
}

void KMInit(const std::string &aIndex, int k = 31)
{
    ProgramOptions o;
    o.k = k;
    o.index = aIndex;
    
    /*
     * Initalize for reference on all sequin k-mers (e.g. sequins.fa.index)
     */
    
    __allIndex__ = std::shared_ptr<KmerIndex>(new KmerIndex(o));
    __allIndex__->load(o);
    
    KMCancerSpan();
}

static void KMCount(const char *s)
{
    Kmer::k = 31;
    KmerIterator iter(s), end;
    
    // Number of k-mers that are sequins
    unsigned isSeq = 0;
    
    // Number of k-mers that are genome
    unsigned isGen = 0;
    
    for (int i = 0; iter != end; ++i, ++iter)
    {
        auto k = iter->first.rep().toString();
        
        /*
         * Does the k-mer span sequin variants?
         */
        
        if (__span__.count(k))
        {
            __span__[k]++;
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

void KMDebug()
{
#ifdef DEBUG
//    std::ofstream w("KMAll.txt");
//
//    for (const auto &i : __debug__)
//    {
//        w << i.first << "\t" << i.second << std::endl;
//    }
//
//    w.close();
#endif
}

/*
 * Heavily modified Kallisto k-mer counting (no EM algorithm, bootstrapping etc)
 */

KMStats Kallisto(const std::string &aIndex, const std::string &p1, const std::string &p2, unsigned k)
{
    // Initalize for reference k-mers
    KMInit(aIndex, k);

    ProgramOptions opt;
    
    opt.index = aIndex;
    opt.files.push_back(p1);
    opt.files.push_back(p2);

    KmerIndex index(opt);
    MinCollector collection(index, opt);
    ProcessReads(index, opt, collection);
    
#ifdef DEBUG
    KMDebug();
#endif

    assert((__kStats__.nGen + __kStats__.nSeq) > 0);
    
//    std::cout << __kStats__.nGen << std::endl;
//    std::cout << __kStats__.nSeq << std::endl;
//    std::cout << (float)__kStats__.nSeq / (nSeq + nGen) << std::endl;
    
    return __kStats__;
}
