#include <sstream>
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

struct Node
{
    Node(char base) : base(base) {}
    
    char base;
    int  pos;

    Node *next = nullptr;
    Node *prev = nullptr;
};

struct ModifedVCF : public ParserVCF::Data
{
    // TODO: Why do we need a default constructor?
    ModifedVCF() {}
    
    ModifedVCF(const ParserVCF::Data &d)
    {
        this->cID = d.cID;
        this->id  = d.id;
        this->l   = d.l;
        this->ref = d.ref;
        this->alt = d.alt;
        orgPos = this->l.start;
    }
    
    Base orgPos = 0;
    Base newPos = 0;
};

struct LinkedSequence
{
    // Hashing for each base
    std::map<Base, const Node *> hash;

    Node *start;
    Node *end;
};

struct Flanks
{
    // Sequence to the left
    std::string lSeq;
    
    // Sequence to the right
    std::string rSeq;
};

LinkedSequence createList(const std::string &str)
{
    Node *start = new Node(str[0]);
    Node *prev  = start;
    
    LinkedSequence l;

    for (auto i = 1; i < str.length(); i++)
    {
        Node *node = new Node(str[i]);
        l.hash[i] = node;

        node->pos  = i+1; // This is only used for testing
        
        prev->next = node;
        node->prev = prev;
        
        prev = node;
    }

    l.start = start;
    l.end = prev;
    
    return l;
    
    //return std::pair<Node *, Node *>(start, prev);
}

static Base print(const std::pair<Node *, Node *> &startEnd)
{
    Base n = 0;
    auto node = startEnd.first;
    
    while (node)
    {
        n++;
        std::cout << node->base;
        node = node->next;
    }
    
    std::cout << std::endl;
    return n;
}

static void exec(const std::string &cmd)
{
    system(cmd.c_str());
}

static void writeBegin(const std::string &file, const ChrID &chrID, const LinkedSequence &l)
{
    std::ofstream f;
    f.open(file);

    f << ">" << chrID << std::endl;
    
    auto node = l.start;
    int i = 0;

    while (node)
    {
        f << node->base;
        node = node->next;
        
        if (i++ == 60)
        {
            i = 0;
            f << std::endl;
        }
    }
    
    f.close();
}

static void writeEnd(const std::string &file, const ChrID &chrID, const LinkedSequence &l)
{
    std::ofstream f;
    f.open(file);
    
    f << ">" << chrID << std::endl;

    auto node = l.end;
    
    while (node)
    {
        f << node->base;
        node = node->prev;
    }
    
    f.close();
}

static void checkSNP(const std::pair<Node *, Node *> &startEnd, const std::string &str, const std::map<Base, ModifedVCF> &vars, bool allowAlt=false)
{
    auto i = str.length();
    auto node = startEnd.second;
    
    while (node)
    {
        if (vars.count(i))
        {
            const auto &var = vars.at(i);
            
            switch (var.type())
            {
                case Mutation::SNP:
                {
                    assert(node->pos == i);
                    assert(var.ref.size() == 1);
                    
                    if (!allowAlt)
                    {
                        assert(var.alt.size() == 1 && std::toupper(node->base) == std::toupper(var.ref[0]));
                    }
                    else
                    {
                        assert(var.alt.size() == 1 && std::toupper(node->base) == std::toupper(var.alt[0]));
                    }
                    
                    break;
                }
                
                default : { break; }
            }
        }
        
        i--;
        node = node->prev;
    }
}

static Flanks readFlank(const LinkedSequence &l, Base p, Base lFlank, Base rFlank)
{
    std::stringstream ss;
    
    auto node = l.hash.at(p - lFlank);

    for (auto i = 0; i < lFlank; i++)
    {
        ss << node->base;
        node = node->prev;
    }
    
    Flanks f;
    f.lSeq = ss.str();

    ss.clear();
    node = l.hash.at(p);
    
    for (auto i = 0; i < rFlank; i++)
    {
        ss << node->base;
        node = node->next;
    }

    f.rSeq = ss.str();

    return f;
}

static Base count(const LinkedSequence &l)
{
    auto n = 1;
    auto node = l.start;
    
    while (node)
    {
        n++;
        node = node->next;
    }

    return n;
}

static char fetchBase(const LinkedSequence &l, int n)
{
    auto node = l.start;

    for (auto i = 1; i < n; i++)
    {
        node = node->next;
    }
    
    return node->base;
}

void FastaAlternateReferenceMaker(LinkedSequence &l, const std::string &str, std::map<Base, ModifedVCF> &vars)
{
    auto n = 0;
    auto i = str.length();
    auto node = l.end;
    
    while (node)
    {
        if (!((n++) % 10000000))
        {
            std::cout << n << "/" << str.length() << std::endl;
        }
        
        if (vars.count(i))
        {
            const auto &var = vars.at(i);
            const auto diff = var.diff();
            
            switch (var.type())
            {
                case Mutation::SNP:
                {
                    // Let's check the base before substitution
                    assert(std::toupper(fetchBase(l, (int)i)) == std::toupper(var.ref[0]));
                    
                    assert(var.alt.size() == 1 && std::toupper(node->base) == std::toupper(var.ref[0]));
                    node->base = var.alt[0];

                    // Let's check the base after substitution
                    assert(std::toupper(fetchBase(l, (int)i)) == std::toupper(var.alt[0]));
                    
                    break;
                }
                    
                case Mutation::Insertion:
                {
                    auto temp = node;

                    // Subtitute for the first character
                    temp->base = var.alt[0];
                    
                    for (auto j = 1; j < var.alt.size(); j++)
                    {
                        Node *insert = new Node(var.alt[j]);
                        
                        insert->next = temp->next;
                        insert->prev = temp;
                        
                        temp->next = insert;
                        insert->next->prev = insert;

                        temp = insert;
                    }
                    
                    for (auto &j : vars)
                    {
                        if (j.first > i)
                        {
                            j.second.l += diff;
                        }
                    }

                    break;
                }
                    
                case Mutation::Deletion:
                {
                    auto temp = node;

                    // Subtitute for the first character
                    temp->base = var.alt[0];
                    
                    // We'll start deletion from here...
                    temp = temp->next;
                    
                    for (auto j = 0; j < diff; j++)
                    {
                        temp->prev->next = temp->next;
                        temp->next->prev = temp->prev;
                        
                        auto del = temp;
                        temp = temp->next;
                        delete del;
                    }
                    
                    for (auto &j : vars)
                    {
                        if (j.first > i)
                        {
                            j.second.l -= diff;
                        }
                    }
                    
                   // node = startEnd.first;
                    
                    break;
                }
            }
        }
        
        i--;
        node = node->prev;
    }
    
    std::map<Base, ModifedVCF> x = vars;
    
    // We'll construct the variants again, since the position have changed...
    vars.clear();
    
    for (const auto &i : x)
    {
        vars[i.second.l.start] = i.second;
    }
}

static std::string readFA(const std::string &file, const ChrID &cID)
{
    std::string seq;
    
    ParserFA::parse(file, [&](const ParserFA::Data &chr, const ParserProgress &)
    {
        seq = chr.seq;
    }, cID);

    std::cout << "Number of bases in " + cID + " is: " + std::to_string(seq.size()) << std::endl;
    assert(!seq.empty());
    
    return seq;
}

typedef std::vector<ParserBed::Data> BedData;

static std::vector<ParserBed::Data> readBed(const FileName &file, const ChrID &cID)
{
    std::vector<ParserBed::Data> x;
    
    ParserBed::parse(file, [&](const ParserBed::Data &d, const ParserProgress &)
    {
        if (d.cID == cID)
        {
            x.push_back(d);
        }
    });
    
    return x;
}

static void readVCF(const std::string &file, const ChrID &cID, std::map<Base, ModifedVCF> &vars)
{
    ParserVCF::parse(file, [&](const ParserVCF::Data &d, const ParserProgress &)
    {
        if (d.cID == cID)
        {
            vars[d.l.start] = ModifedVCF(d);

            assert(!vars[d.l.start].ref.empty());
            assert(!vars[d.l.start].alt.empty());
        }
    });
    
    std::cout << "Number of variants: " + std::to_string(vars.size()) << std::endl;
    assert(!vars.empty());
}

static void writeBed(const std::string &file, const std::map<Base, ModifedVCF> &vars, const std::string &chrID)
{
    std::ofstream f;
    f.open(file);

    std::vector<ModifedVCF> x;
    
    for (const auto &var : vars)
    {
        x.push_back(var.second);
    }
    
    std::sort(x.begin(), x.end());
    
    for (const auto &var : x)
    {
        f << (boost::format("%1%\t%2%\t%3%\t%4%\n") % chrID
                                                    % var.newPos
                                                    % var.newPos
                                                    % (std::to_string(var.orgPos) + "_" + var.ref + "_" + var.alt)).str();
    }
    
    f.close();
}

static void writeOldBed(const std::string &file, const std::map<Base, ModifedVCF> &vars, const std::string &chrID)
{
    std::ofstream f;
    f.open(file);
    
    std::vector<ModifedVCF> x;
    
    for (const auto &var : vars)
    {
        x.push_back(var.second);
    }
    
    std::sort(x.begin(), x.end());
    
    for (const auto &var : x)
    {
        f << (boost::format("%1%\t%2%\t%3%\t%4%\n") % chrID
                                                    % var.l.start
                                                    % var.l.start
                                                    % (std::to_string(var.orgPos) + "_" + var.ref + "_" + var.alt)).str();
    }
    
    f.close();
}

static std::pair<std::string, std::map<Base, ModifedVCF>> createTestData()
{
    /*
     * Say we have the following sequence: 0123456789ABCDEFGHIJKLMNOPQ and:
     *
     #    - SNP:       pos 5, from '5' to '#'
     #    - Insertion: pos 9, from '9' to @@@@
     #    - Deletion:  pos 13, DEFGH to &
     #    - SNP:       pos 26, from 'Q' to '!'
     *
     *      01234#678@@@@ABC&IJKLMNOPQ
     *      !PONMLKJI&CBA@@@@876#43210
     */

    ParserVCF::Data v1;
    v1.l = Locus(5, 5);
    v1.ref = "5";
    v1.alt = "#";
    assert(v1.type() == Mutation::SNP);
    
    ParserVCF::Data v2;
    v2.l = Locus(9, 9);
    v2.ref = "9";
    v2.alt = "@@@@";
    assert(v2.type() == Mutation::Insertion);

    ParserVCF::Data v3;
    v3.l = Locus(13, 13);
    v3.ref = "DEFGH";
    v3.alt = "&";
    assert(v3.type() == Mutation::Deletion);

    ParserVCF::Data v4;
    v4.l = Locus(26, 26);
    v4.ref = "Q";
    v4.alt = "!";
    assert(v4.type() == Mutation::SNP);

    std::map<Base, ModifedVCF> m;
    
    m[v1.l.start] = ModifedVCF(v1);
    m[v2.l.start] = ModifedVCF(v2);
    m[v3.l.start] = ModifedVCF(v3);
    m[v4.l.start] = ModifedVCF(v4);

    return std::pair<std::string, std::map<Base, ModifedVCF>>("0123456789ABCDEFGHIJKLMNOPQ", m);
}

static void reverse(const FileName &gFile, const FileName &vFile, const ChrID &chrID)
{
    std::map<Base, ModifedVCF> vars;
    
    const auto begin = clock();
    
    std::cout << "Chromosome: " << chrID << std::endl;
    
    std::cout << "Reading genome..." << std::endl;
    const auto seq = readFA(gFile, chrID);
    
    std::cout << "Reading variants..." << std::endl;
    readVCF(vFile, chrID, vars);
    
    //auto test = createTestData();
    //seq  = test.first;
    //vars = test.second;
    
    std::cout << "Creating linked-list..." << std::endl;
    
    // Create a linked-list representation for the sequence
    auto l = createList(seq);
    
    std::cout << "Linked-list created" << std::endl;
    
    //const auto n_vars = vars.size();
    
    writeBegin("origGenome.fa", chrID, l);
    writeOldBed("origGenome.bed", vars,  chrID);
    
    std::cout << "Generated: origGenome.fa"  << std::endl;
    std::cout << "Generated: origGenome.bed" << std::endl;
    
    //checkSNP(nodes, seq, vars);
    
    //FastaAlternateReferenceMaker(nodes, seq, vars);
    
    // Length of the flipped sequence
    const auto n = count(l);
    
    std::cout << "Number of bases in the flipped sequence: " << n << std::endl;
    
    /*
     * Generating a FASTA file for the variant genome.
     */
    
    writeBegin("/Users/tedwong/Desktop/Ira/varGenome.fa", chrID, l);
    writeOldBed("/Users/tedwong/Desktop/Ira/varGenome.bed", vars, chrID);
    
    std::cout << "Generated: varGenome.fa"  << std::endl;
    std::cout << "Generated: varGenome.bed" << std::endl;
    
    /*
     * Generating the variants after flipping. This is only possible because we know the size of the chromosome.
     */
    
    for (auto &var : vars)
    {
        const auto pos = var.second.l.start;
        
        switch (var.second.type())
        {
            case Mutation::SNP:
            {
                // What's the position of the SNP after reversing? It's n-(i+1)
                var.second.newPos = n-pos;
                
                break;
            }
                
            case Mutation::Insertion:
            {
                // Number of characters inserted
                //const auto offset = var.second.alt.size() - var.second.ref.size();
                
                var.second.newPos = n-pos;
                
                std::cout << var.second.newPos << std::endl;
                
                break;
            }
                
            case Mutation::Deletion:
            {
                // Number of characters inserted
                const auto offset = var.second.ref.size() - var.second.alt.size();
                
                // This'll be our new position
                var.second.newPos = n-pos-offset;
                
                break;
            }
        }
    }
    
    /*
     * Now, flip the genome... This is easy with a linked-list implementation... In fact we don't have to do anything...
     */
    
    /*
     * Generating a FASTA file for the flipped genome
     */
    
    writeEnd("/Users/tedwong/Desktop/Ira/flipGenome.fa", chrID, l);
    std::cout << "Generated: flipGenome.fa" << std::endl;
    
    /*
     * Generating VCF file for the flipped genome
     */
    
    writeBed("/Users/tedwong/Desktop/Ira/flipGenome.bed", vars, chrID);
    std::cout << "Generated: flipGenome.bed" << std::endl;
    
    const auto end = clock();
    std::cout << "Completed in: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;
}

static void writeFA(const FileName &file, const std::map<std::string, std::string> &m)
{
    std::ofstream f;
    f.open(file);
    
    for (const auto &i : m)
    {
        f << ">" << i.first << std::endl;
        f << i.second << std::endl;
    }

    f.close();
}

typedef std::map<std::string, std::string> KeyToSeq;

static void readFlanks(const LinkedSequence &l, const BedData &bs, const BedData &be, KeyToSeq &A, KeyToSeq &B, KeyToSeq &C, KeyToSeq &D, Base flank)
{
    for (const auto &i : bs)
    {
        const auto fls = readFlank(l, i.l.start, flank, flank);
        
        /*
         * Divided the sequence into A and B. We'll merge them later.
         */
        
        A[i.id] = fls.lSeq;
        B[i.id] = fls.rSeq;
    }
    
    for (const auto &i : be)
    {
        const auto fls = readFlank(l, i.l.start, flank, flank);
        
        /*
         * Divided the sequence into C and D. We'll merge them later.
         */
        
        C[i.id] = fls.lSeq;
        D[i.id] = fls.rSeq;
    }
    
    assert(A.size() == B.size());
    assert(B.size() == C.size());
    assert(C.size() == D.size());
}

static void reverseComp(KeyToSeq &x)
{
    for (auto &i : x)
    {
        i.second = i.second;
    }
}

static void invert(KeyToSeq &x)
{
    for (auto &i : x)
    {
        std::reverse(i.second.begin(), i.second.end());
    }
}

static void merge(const KeyToSeq &src1, const KeyToSeq &src2, KeyToSeq &dst)
{
    for (const auto &i : src1)
    {
        const auto x = src1.at(i.first);
        const auto y = src2.at(i.first);
        dst[i.first] = x + y;
    }
}

/*
 * Create sequins for strucutal inversions.
 */

static void inversion(const FileName &gFile,
                      const FileName &bStart,
                      const FileName &bEnd,
                      const ChrID &cID,
                      Base flank = 600)
{
    std::cout << "Reading genome..." << std::endl;
    const auto seq = readFA(gFile, cID);
    
    std::cout << "Reading BED for starting..." << std::endl;
    const auto bs = readBed(bStart, cID);
    
    std::cout << "Reading BED for stopping..." << std::endl;
    const auto be = readBed(bEnd, cID);
    
    std::cout << "Reading linked representation..." << std::endl;
    const auto l = createList(seq);
    
    KeyToSeq A, B, C, D;
    readFlanks(l, bs, be, A, B, C, D, flank);

    KeyToSeq A_, B_, C_, D_;
    readFlanks(l, bs, be, A_, B_, C_, D_, flank);
    
    /*
     * Perform a reverse complement on B and C
     */
    
    reverseComp(B_);
    std::cout << "Reverse complementing B" << std::endl;
    
    reverseComp(C_);
    std::cout << "Reverse complementing C" << std::endl;

    KeyToSeq AB, CD, AC_, B_D;
    
    /*
     * Merging and inverting A and B
     */
    
    merge(A, B, AB);
    invert(AB);
    
    /*
     * Merging and inverting C and D
     */
    
    merge(C, D, CD);
    invert(CD);

    /*
     * Merging and inverting A and C_
     */
    
    merge(A, C_, AC_);
    invert(AC_);
    
    /*
     * Merging and inverting B_ and D
     */
    
    merge(B_, D, B_D);
    invert(B_D);
    
    writeFA("sequins.invertion.AB.rev.fa", AB);
    std::cout << "Generated: sequins.invertion.AB.rev.fa" << std::endl;
    
    writeFA("sequins.invertion.CD.rev.fa", CD);
    std::cout << "Generated: sequins.invertion.CD.rev.fa" << std::endl;
    
    writeFA("sequins.invertion.AC_.rev.fa", AC_);
    std::cout << "Generated: sequins.invertion.AC_.rev.fa" << std::endl;

    writeFA("sequins.invertion.B_D.rev.fa", B_D);
    std::cout << "Generated: sequins.invertion.B_D.rev.fa" << std::endl;
}

// Eg: seqTools invert hg38.fa inversions.large.bed 300
static void invertByScript(const FileName &genome, const FileName &bed, Base flank)
{
    const auto fs = std::to_string(flank);
    
    exec("returnEnds.pl -end=tss  " + bed + " > /tmp/inversions.start.bed");
    exec("returnEnds.pl -end=3utr " + bed + " > /tmp/inversions.stop.bed");

    std::cout << "Generating flanking regions" << std::endl;
    
    auto x = (boost::format("bedEnds2.sh -%1% %2% /tmp/inversions.start.bed /tmp/inversions.AB.bed")
              % flank
              % (flank - 1)).str();
    
    auto y = (boost::format("bedEnds2.sh -%1% %2% /tmp/inversions.stop.bed /tmp/inversions.CD.bed")
              % flank
              % (flank - 1)).str();
    exec(x);
    exec(y);
    
    /*
     * Extract the flanking sequence from the genome
     */
    
    std::cout << "Extracting from the genome..." << std::endl;
    
    exec("sequenceForBed.pl -case=upper " + genome + " /tmp/inversions.AB.bed > /tmp/inversions.AB.tab.fa");
    exec("sequenceForBed.pl -case=upper " + genome + " /tmp/inversions.CD.bed > /tmp/inversions.CD.tab.fa");
    
    /*
     * Break the sequence from AB to A
     */
    
    std::cout << "Generating A..." << std::endl;
    exec("fetchSubsequence.pl -end=5 -length=" + fs + " /tmp/inversions.AB.tab.fa > /tmp/inversions.A.tab.fa");
    
    /*
     * Break the sequence from AB to B
     */
    
    std::cout << "Generating B..." << std::endl;
    exec("fetchSubsequence.pl -end=3 -length=" + fs + " /tmp/inversions.AB.tab.fa > /tmp/inversions.B.tab.fa");

    /*
     * Break the sequence from CD to C
     */
    
    std::cout << "Generating C..." << std::endl;
    exec("fetchSubsequence.pl -end=5 -length=" + fs + " /tmp/inversions.CD.tab.fa > /tmp/inversions.C.tab.fa");
    
    /*
     * Break the sequence from CD to D
     */
    
    std::cout << "Generating C..." << std::endl;
    exec("fetchSubsequence.pl -end=3 -length=" + fs + " /tmp/inversions.CD.tab.fa > /tmp/inversions.D.tab.fa");
    
    exec("flattenFasta.pl -fa /tmp/inversions.B.tab.fa > /tmp/revtmp.fa; revcomp /tmp/revtmp.fa > /tmp/revtmp2.fa; flattenFasta.pl -tab /tmp/revtmp2.fa | sed 's/ //' > /tmp/inversions.B.revcomp.tab.fa");
    exec("flattenFasta.pl -fa /tmp/inversions.C.tab.fa > /tmp/revtmp.fa; revcomp /tmp/revtmp.fa > /tmp/revtmp2.fa; flattenFasta.pl -tab /tmp/revtmp2.fa | sed 's/ //' > /tmp/inversions.C.revcomp.tab.fa");
    
    exec("mergeMultipleTab.pl /tmp/inversions.A.tab.fa /tmp/inversions.B.tab.fa | sed 's/\t//2' > /tmp/inversions.AB.tab.fa");
    exec("mergeMultipleTab.pl /tmp/inversions.C.tab.fa /tmp/inversions.D.tab.fa | sed 's/\t//2' > /tmp/inversion.CD.tab.fa");
    exec("mergeMultipleTab.pl /tmp/inversions.A.tab.fa /tmp/inversions.C.revcomp.tab.fa | sed 's/\t//2' > /tmp/inversions.ACrc.tab.fa");
    exec("mergeMultipleTab.pl /tmp/inversions.B.revcomp.tab.fa /tmp/inversions.D.tab.fa | sed 's/\t//2' > /tmp/inversions.BrcD.tab.fa");
    
    exec("cut -f 2 /tmp/inversions.AB.tab.fa | rev > /tmp/revseq.txt");
    exec("cut -f 1 /tmp/inversions.AB.tab.fa > /tmp/ids.txt");
    exec("paste /tmp/ids.txt /tmp/revseq.txt > sequins.inversions.AB.rev.tab.fa");
    
    exec("cut -f 2 /tmp/inversions.AB.tab.fa | rev > /tmp/revseq.txt");
    exec("cut -f 1 /tmp/inversions.AB.tab.fa > /tmp/ids.txt");
    exec("paste /tmp/ids.txt /tmp/revseq.txt > sequins.inversions.AB.rev.tab.fa");

    exec("cut -f 2 /tmp/inversions.CD.tab.fa | rev > /tmp/revseq.txt");
    exec("cut -f 1 /tmp/inversions.CD.tab.fa > /tmp/ids.txt");
    exec("paste /tmp/ids.txt /tmp/revseq.txt > sequins.inversions.CD.rev.tab.fa");

    exec("cut -f 2 /tmp/inversions.ACrc.tab.fa | rev > /tmp/revseq.txt");
    exec("cut -f 1 /tmp/inversions.ACrc.tab.fa > /tmp/ids.txt");
    exec("paste /tmp/ids.txt /tmp/revseq.txt > sequins.inversions.ACrc.rev.tab.fa");

    exec("cut -f 2 /tmp/inversions.BrcD.tab.fa | rev > /tmp/revseq.txt");
    exec("cut -f 1 /tmp/inversions.BrcD.tab.fa > /tmp/ids.txt");
    exec("paste /tmp/ids.txt /tmp/revseq.txt > sequins.inversions.BrcD.rev.tab.fa");
    
    std::cout << "Please check: sequins.inversions.AB.rev.tab.fa"   << std::endl;
    std::cout << "Please check: sequins.inversions.CD.rev.tab.fa"   << std::endl;
    std::cout << "Please check: sequins.inversions.ACrc.rev.tab.fa" << std::endl;
    std::cout << "Please check: sequins.inversions.BrcD.rev.tab.fa" << std::endl;
}

// Eg: seqTools deletion hg38.fa Personalis.deletions.large.bed 300
static void deletionByScript(const FileName &genome, const FileName &bed, Base flank)
{
    const auto fs = std::to_string(flank);
    
    /*
     * Get the starting base for each variant (tss) or ending base (3utr)
     */
    
    exec("returnEnds.pl -end=tss  " + bed + " > /tmp/deletions.start.bed");
    exec("returnEnds.pl -end=3utr " + bed + " > /tmp/deletions.stop.bed");
    
    /*
     * Generate flanking regions
     */
    
    std::cout << "Generating flanking regions" << std::endl;

    auto x = (boost::format("bedEnds2.sh -%1% %2% /tmp/deletions.start.bed /tmp/deletions.AB.bed")
                                    % flank
                                    % (flank - 1)).str();

    auto y = (boost::format("bedEnds2.sh -%1% %2% /tmp/deletions.stop.bed /tmp/deletions.CD.bed")
                                    % flank
                                    % (flank - 1)).str();
    exec(x);
    exec(y);

    /*
     * Extract the flanking sequence from the genome
     */

    std::cout << "Extracting from the genome..." << std::endl;
    
    exec("sequenceForBed.pl -case=upper " + genome + " /tmp/deletions.AB.bed > /tmp/deletions.AB.tab.fa");
    exec("sequenceForBed.pl -case=upper " + genome + " /tmp/deletions.CD.bed > /tmp/deletions.CD.tab.fa");

    /*
     * Break the sequence from AB to A
     */
    
    std::cout << "Generating A..." << std::endl;
    exec("fetchSubsequence.pl -end=5 -length=" + fs + " /tmp/deletions.AB.tab.fa > /tmp/deletions.A.tab.fa");

    /*
     * Break the sequence from AB to B
     */
    
    std::cout << "Generating B..." << std::endl;
    exec("fetchSubsequence.pl -end=3 -length=" + fs + " /tmp/deletions.AB.tab.fa > /tmp/deletions.B.tab.fa");
    
    /*
     * Break the sequence from CD to C
     */
    
    std::cout << "Generating C..." << std::endl;
    exec("fetchSubsequence.pl -end=5 -length=" + fs + " /tmp/deletions.CD.tab.fa > /tmp/deletions.C.tab.fa");
    
    /*
     * Break the sequence from CD to D
     */
    
    std::cout << "Generating D..." << std::endl;
    exec("fetchSubsequence.pl -end=3 -length=" + fs + " /tmp/deletions.CD.tab.fa > /tmp/deletions.D.tab.fa");
    
    /*
     * Merging A and B
     */
    
    std::cout << "Merging..." << std::endl;
    exec("mergeMultipleTab.pl /tmp/deletions.A.tab.fa /tmp/deletions.B.tab.fa | sed s'/\t//2' > deletions.AB.tab.fa");
    
    /*
     * Merging A and D
     */
    
    std::cout << "Merging..." << std::endl;
    exec("mergeMultipleTab.pl /tmp/deletions.A.tab.fa /tmp/deletions.D.tab.fa | sed s'/\t//2' > deletions.AD.tab.fa");
    
    /*
     * Merging C and D
     */
    
    std::cout << "Merging..." << std::endl;
    exec("mergeMultipleTab.pl /tmp/deletions.C.tab.fa /tmp/deletions.D.tab.fa | sed s'/\t//2' > deletions.CD.tab.fa");

    std::cout << "Please check: deletions.AB.tab.fa" << std::endl;
    std::cout << "Please check: deletions.AD.tab.fa" << std::endl;
    std::cout << "Please check: deletions.CD.tab.fa" << std::endl;
}

/*
 * Create sequins for strucutal deletions.
 */

static void deletion(const FileName &gFile,
                     const FileName &bStart,
                     const FileName &bEnd,
                     const ChrID &cID,
                     Base flank = 600)
{
    std::cout << "Reading genome..." << std::endl;
    const auto seq = readFA(gFile, cID);
    
    std::cout << "Reading BED for starting..." << std::endl;
    const auto bs = readBed(bStart, cID);

    std::cout << "Reading BED for stopping..." << std::endl;
    const auto be = readBed(bEnd, cID);
    
    std::cout << "Reading linked representation..." << std::endl;
    const auto l = createList(seq);

    std::map<std::string, std::string> A, B, C, D;
    readFlanks(l, bs, be, A, B, C, D, flank);

    std::map<std::string, std::string> AB, AD, CB;
    
    /*
     * Merging A and B
     */
    
    merge(A, B, AB);
    std::cout << "Merged A and B. Check out: deletions.AB.fa" << std::endl;

    /*
     * Merging A and D
     */
    
    merge(A, D, AD);
    std::cout << "Merged A and D. Check out: deletions.AD.fa" << std::endl;

    /*
     * Merging C and B
     */
    
    merge(C, B, CB);
    std::cout << "Merged C and D. Check out: deletions.CD.fa" << std::endl;
    
    /*
     * Invert them to minmise cross-alignment errors
     */

    invert(AB);
    invert(AD);
    invert(CB);
    
    writeFA("sequins.deletions.AB.rev.fa", AB);
    std::cout << "Generated: sequins.deletions.AB.rev.fa" << std::endl;
    
    writeFA("sequins.deletions.AD.rev.fa", AD);
    std::cout << "Generated: sequins.deletions.AD.rev.fa" << std::endl;
    
    writeFA("sequins.deletions.CB.rev.fa", CB);
    std::cout << "Generated: sequins.deletions.CB.rev.fa" << std::endl;
}

/*
 *  Usage: seqTools <genomeFile> <variantFile> <Chromosome>
 *
 *    Eg: seqTools reverse chr21.fa chr21.vcf chr21
 *        seqTools
 *
 *        g++ -std=c++11 -c -I /usr/include/boost -I ~/Sources/QA/src -I ~/Sources/SS seqTools.cpp
 *        g++ -std=c++11 -c -I /usr/include/boost -I ~/Sources/QA/src ~/Sources/QA/src/data/reader.cpp
 *        g++ *.o -o seqTools
 */

int main(int argc, const char * argv[])
{
    std::cout << "Sequins Tool, v2.1." << std::endl << std::endl;
    
    if (argc == 1)
    {
        std::cout << "seqTools deletion <genome (eg: hg38.fa)>  <deletion bed (eg: deletion.bed)> <flanking (eg: 300)" << std::endl;
        std::cout << "seqTools invert <genome (eg: hg38.fa)>  <deletion bed (eg: deletion.bed)> <flanking (eg: 300)" << std::endl;
    }
    else
    {
        const auto mode = std::string(argv[1]);
        
        if (mode == "reverse")
        {
            reverse(argv[2], argv[3], argv[4]);
        }
        else if (mode == "deletion")
        {
            deletionByScript(argv[2], argv[3], std::stoi(argv[4]));
        }
        else if (mode == "invert")
        {
            invertByScript(argv[2], argv[3], std::stoi(argv[4]));
        }
    }
    
    return 0;
}