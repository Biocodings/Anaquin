#include <fstream>
#include <assert.h>
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"
#include "tools/errors.hpp"
#include "tools/system.hpp"
#include "VarQuin/VarQuin.hpp"

using namespace Anaquin;

typedef std::vector<std::pair<KmerEntry, int>> KmerEntries;

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

static bool KCount(const char *s)
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
        
        if (nMatch > nNMatch) { k.nMatch++;  }
        else                  { k.nNMatch++; }
        
        return nMatch;
    };

    return __match__(__i1__, __kStats__.R) || __match__(__i2__, __kStats__.F);
}

void KCount(const char *s1, const char *s2)
{
    KCount(s1);
    KCount(s2);
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

KStats Anaquin::KCount(const FileName &i1, const FileName &i2, const FileName &p1, const FileName &p2, Counts threads, unsigned k)
{
    KInit(i1, i2, k);

    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = i1;
    opt.files.push_back(p1);
    opt.files.push_back(p2);
    opt.threads = threads;
    
    KmerIndex index(opt);
    MinCollector collection(index, opt);
    
    // Process k-mers
    ProcessReads(index, opt, collection);
    
    return __kStats__;
}
