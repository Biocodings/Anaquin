#include <map>
#include <fstream>
#include <assert.h>
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"

using namespace Anaquin;

extern std::string KMVarKStats();

// Index for all sequin k-mers
static std::shared_ptr<KmerIndex> __allIndex__;

static KMStats __kStats__;

/*
 * Initlaize k-mers spanning variants, they are useful for estimating allele frequency
 */

static void LoadKMSpan()
{
    std::istringstream r(KMVarKStats());
    std::string line;
    
    while (std::getline(r, line))
    {
        std::vector<std::string> toks;
        split(line, "\t", toks);
        
        if (toks.size() == 4 && toks[3] == "Span")
        {
            if (isSubstr(toks[0], "_R"))
            {
                __kStats__.vars[noLast(toks[0], "_")].R.insert(toks[1]);
            }
            else
            {
                __kStats__.vars[noLast(toks[0], "_")].V.insert(toks[2]);
            }
            
            __kStats__.spans[toks[1]] = 0; // Normal
            __kStats__.spans[toks[2]] = 0; // Reverse complement
        }
    }
    
    assert(!__kStats__.vars.empty());
    assert(!__kStats__.spans.empty());
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
    
    LoadKMSpan();
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

/*
 * Heavily modified Kallisto k-mer counting (no EM algorithm, bootstrapping etc)
 */

KMStats Kallisto(const std::string &aIndex, const std::string &p1, const std::string &p2, unsigned k)
{
    Kmer::k = k;
    
    // Initalize for reference k-mers
    KMInit(aIndex, k);

    ProgramOptions opt;
    
    opt.index = aIndex;
    opt.files.push_back(p1);
    opt.files.push_back(p2);

    KmerIndex index(opt);
    MinCollector collection(index, opt);
    
    // Process k-mers
    ProcessReads(index, opt, collection);
    
    // Have we read anything?
    assert((__kStats__.nGen + __kStats__.nSeq) > 0);
    
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
