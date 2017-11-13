#include <mutex>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"
#include "tools/errors.hpp"
#include "tools/system.hpp"

using namespace Anaquin;

typedef std::vector<std::pair<KmerEntry, int>> KmerEntries;

// For writing out to files
static std::mutex __lock__;

// Files for genome and sequin reads
std::shared_ptr<std::ofstream> __gens__, __seqs__, __fors__;

// Kallisto indexes
static std::shared_ptr<KmerIndex> __i1__, __i2__;

// Running statistics for k-mers
static KStats __kStats__;

void Anaquin::KInit(const FileName &i1, const FileName &i2, unsigned k)
{
    auto index = [&](const FileName &i)
    {
        ProgramOptions o;
        o.k = k;
        o.index = i;
        
        auto x = std::shared_ptr<KmerIndex>(new KmerIndex(o));
        x->load(o);
        return x;
    };
    
    __i1__ = index(i1);
    __i2__ = index(i2);
    
    for (const auto &name : __i1__->target_names_)
    {
        __kStats__.seqs.insert(isSubstr(name, "_R") || isSubstr(name, "_V") ? noLast(name, "_") : name);
    }

    A_ASSERT(!__kStats__.seqs.empty());
}

enum class ReadStatus
{
    ReverseSequin,
    ForwardSequin,
    Genome,
};

static ReadStatus KCount(const char *s)
{
    auto __match__ = [&](std::shared_ptr<KmerIndex> index, KStats::KAbund &k)
    {
        KmerIterator ki(s), ke;
        
        Counts nMatch  = 0;
        Counts nNMatch = 0;
        
        for (int i = 0;  ki != ke; ++i, ++ki)
        {
            auto search = index->kmap.find(ki->first.rep());

            // Can we find this k-mer?
            if (search != index->kmap.end())
            {
                nMatch++;
                const auto km = search->first.toString();

                // Matching sequins
                const auto &trans = index->dbGraph.contigs[search->second.contig].transcripts;
                
                std::set<SequinID> seqs;                
                for (const auto &tran : trans)
                {
                    const auto seq = index->target_names_[tran.trid];                    
                    seqs.insert(seq);
                }

                // Unique matching?
                const auto isUniq = seqs.size() == 1;
                
                for (const auto &seq : seqs)
                {
                    if (!isUniq)
                    {
                        k.shared[seq][km]++;
                    }
                    else
                    {
                        k.uniqs[seq][km]++;
                    }

                    k.m[seq]++;
                }
            }
            else
            {
                nNMatch++;
            }
        }
        
        // Match this index
        if (nMatch > nNMatch)
        {
            k.nMatch++;
            
            // How many of the k-mers matching the reference?
            k.nMKMatch += nMatch;
            
            // How many of the k-emrs not matching the reference?
            k.nNMKMatch += nNMatch;
        }
        
        // Doesn't match this index
        else
        {
            k.nNMatch++;
        }

        return nMatch;
    };

    if (__match__(__i1__, __kStats__.R))
    {
        return ReadStatus::ReverseSequin;
    }
    else if (__match__(__i2__, __kStats__.F))
    {
        return ReadStatus::ForwardSequin;
    }
    else
    {
        return ReadStatus::Genome;
    }
}

void KCount(const char *r1, const char *s1, const char *r2, const char *s2)
{
    auto __write__ = [&](const char *r, const char *s)
    {
        const auto x = KCount(s);
        
        if (__seqs__)
        {
            __lock__.lock();
            
            switch (x)
            {
                case ReadStatus::ReverseSequin: { if (__seqs__) { *(__seqs__) << r << std::endl; } break; }
                case ReadStatus::ForwardSequin: { if (__fors__) { *(__fors__) << r << std::endl; } break; }
                case ReadStatus::Genome:        { if (__gens__) { *(__gens__) << r << std::endl; } break; }
            }
            
            __lock__.unlock();
        }
    };
    
    __write__(r1, s1);
    __write__(r2, s2);
}

FileName Anaquin::KHumanFA(const FileName &file, std::map<SequinID, Base> &s2l)
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

    // Current sequin
    SequinID sID;
    
    while (((buf[++i] = r.get()) != EOF))
    {
        switch (buf[i])
        {
            case '\n':
            {
                if (isName)
                {
                    buf[i] = '\0';
                    
                    // We need the sequin name later after it's sequence read
                    w << (sID = std::string(buf));
                    
                    reset();
                }
                else
                {
                    write();
                }
                
                isName = false;
                break;
            }

            case '>':
            {
                if (i)
                {
                    // By definition, the last character before '>' here is '\n'
                    buf[i-1] = '\0';
                    
                    std::string x(buf);
                    s2l[sID] = x.size();
                    
                    // Reverse for the human region
                    std::reverse(x.begin(), x.end());
                    
                    w << x << '\n';
                    buf[i] = '>';
                }

                isName = true;
                write();
                reset();
                
                break;
            }

            default: { break; }
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
        throw std::runtime_error("Failed to build index for " + file);
    }
    
    return opt.index;
}

static Matches KQuery(KmerIndex &i, const KmerEntries &v)
{
    Matches r;
    
    for (auto &e : v)
    {
        const auto &c = i.dbGraph.contigs[e.first.contig];
        A_ASSERT(!c.transcripts.empty());
        
        for (auto &j : c.transcripts)
        {
            r.insert(i.target_names_[j.trid]);
        }
    }
    
    return r;
}

Matches Anaquin::KQuery(const FileName &file, const Sequence &s, unsigned k)
{
    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = file;

    KmerIndex index(opt);
    index.load(opt);

    KmerEntries v;
    index.match(s.c_str(), opt.k, v);

    return KQuery(index, v);
}

KStats Anaquin::KCount(const FileName &i1, const FileName &i2, const FileName &p1, const FileName &p2,                   const FileName &gReads, const FileName &fReads, const FileName &sReads, Counts threads, unsigned k)
{
    KInit(i1, i2, k);

    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = i1;
    opt.files.push_back(p1);
    opt.files.push_back(p2);
    opt.threads = threads;
    opt.fusion = !gReads.empty() && !sReads.empty();
    
    KmerIndex index(opt);
    MinCollector collection(index, opt);
    
    if (!gReads.empty() && !sReads.empty())
    {
        __gens__ = std::shared_ptr<std::ofstream>(new std::ofstream(gReads));
        __fors__ = std::shared_ptr<std::ofstream>(new std::ofstream(fReads));
        __seqs__ = std::shared_ptr<std::ofstream>(new std::ofstream(sReads));
    }
    
    // Process k-mers
    ProcessReads(index, opt, collection);
    
    if (__gens__) { __gens__->close(); }
    if (__fors__) { __fors__->close(); }
    if (__seqs__) { __seqs__->close(); }

    return __kStats__;
}
