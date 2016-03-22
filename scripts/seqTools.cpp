#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include "parsers/parser_fa.hpp"
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
        this->chrID = d.chrID;
        this->id = d.id;
        this->l = d.l;
        this->ref = d.ref;
        this->alt = d.alt;
        orgPos = this->l.start;
    }
    
    Base orgPos = 0;
    Base newPos = 0;
};

std::pair<Node *, Node *> createList(const std::string &str)
{
    Node *start = new Node(str[0]);
    Node *prev  = start;
    
    for (auto i = 1; i < str.length(); i++)
    {
        Node *node = new Node(str[i]);
        node->pos  = i+1; // This is only used for testing
        
        prev->next = node;
        node->prev = prev;
        
        prev = node;
    }

    return std::pair<Node *, Node *>(start, prev);
}

Base print(const std::pair<Node *, Node *> &startEnd)
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

static void writeBegin(const std::string &file, const ChrID &chrID, const std::pair<Node *, Node *> &startEnd)
{
    std::ofstream f;
    f.open(file);

    f << ">" << chrID << std::endl;
    
    auto node = startEnd.first;
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

static void writeEnd(const std::string &file, const ChrID &chrID, const std::pair<Node *, Node *> &startEnd)
{
    std::ofstream f;
    f.open(file);
    
    f << ">" << chrID << std::endl;

    auto node = startEnd.second;
    
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
                case VarType::SNP:
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

static Base count(const std::pair<Node *, Node *> &startEnd)
{
    auto n = 1;
    auto node = startEnd.first;
    
    while (node)
    {
        n++;
        node = node->next;
    }

    return n;
}

static char fetchBase(const std::pair<Node *, Node *> &startEnd, int n)
{
    auto node = startEnd.first;

    for (auto i = 1; i < n; i++)
    {
        node = node->next;
    }
    
    return node->base;
}

void FastaAlternateReferenceMaker(const std::pair<Node *, Node *> &startEnd, const std::string &str, std::map<Base, ModifedVCF> &vars)
{
    auto n = 0;
    auto i = str.length();
    auto node = startEnd.second;
    
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
                case VarType::SNP:
                {
                    // Let's check the base before substitution
                    assert(std::toupper(fetchBase(startEnd, (int)i)) == std::toupper(var.ref[0]));
                    
                    assert(var.alt.size() == 1 && std::toupper(node->base) == std::toupper(var.ref[0]));
                    node->base = var.alt[0];

                    // Let's check the base after substitution
                    assert(std::toupper(fetchBase(startEnd, (int)i)) == std::toupper(var.alt[0]));
                    
                    break;
                }
                    
                case VarType::Insertion:
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
                    
                case VarType::Deletion:
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

void readFA(const std::string &file, const ChrID &chrID, std::string &seq)
{
    ParserFA::parse(file, [&](const ParserFA::Data &chr, const ParserProgress &)
    {
        seq = chr.seq;
    }, chrID);

    std::cout << "Number of bases in " + chrID + " is: " + std::to_string(seq.size()) << std::endl;
    assert(!seq.empty());
}

void readVCF(const std::string &file, const ChrID &chrID, std::map<Base, ModifedVCF> &vars)
{
    ParserVCF::parse(file, [&](const ParserVCF::Data &d, const ParserProgress &)
    {
        if (d.chrID == chrID)
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
    assert(v1.type() == VarType::SNP);
    
    ParserVCF::Data v2;
    v2.l = Locus(9, 9);
    v2.ref = "9";
    v2.alt = "@@@@";
    assert(v2.type() == VarType::Insertion);

    ParserVCF::Data v3;
    v3.l = Locus(13, 13);
    v3.ref = "DEFGH";
    v3.alt = "&";
    assert(v3.type() == VarType::Deletion);

    ParserVCF::Data v4;
    v4.l = Locus(26, 26);
    v4.ref = "Q";
    v4.alt = "!";
    assert(v4.type() == VarType::SNP);

    std::map<Base, ModifedVCF> m;
    
    m[v1.l.start] = ModifedVCF(v1);
    m[v2.l.start] = ModifedVCF(v2);
    m[v3.l.start] = ModifedVCF(v3);
    m[v4.l.start] = ModifedVCF(v4);

    return std::pair<std::string, std::map<Base, ModifedVCF>>("0123456789ABCDEFGHIJKLMNOPQ", m);
}

/*
 *  Usage: seqTools <genomeFile> <variantFile> <Chromosome>
 *
 *    Eg: seqTools chr21.fa chr21.vcf chr21
 *
 *        g++ -std=c++11 -c -I /usr/include/boost -I ~/Sources/QA/src -I ~/Sources/SS seqTools.cpp
 *        g++ -std=c++11 -c -I /usr/include/boost -I ~/Sources/QA/src ~/Sources/QA/src/data/reader.cpp
 *        g++ *.o -o seqTools
 */

int main(int argc, const char * argv[])
{
    std::string seq;
    std::map<Base, ModifedVCF> vars;

    const auto gFile = argv[1];
    const auto vFile = argv[2];
    const auto chrID = argv[3];

    const auto begin = clock();

    std::cout << "Chromosome: " << chrID << std::endl;
    
    std::cout << "Reading genome..." << std::endl;
    readFA(gFile, chrID, seq);
    std::cout << "Reading genome completed" << std::endl;

    std::cout << "Reading variants..." << std::endl;
    readVCF(vFile, chrID, vars);
    std::cout << "Reading variants completed" << std::endl;

    //auto test = createTestData();
    //seq  = test.first;
    //vars = test.second;

    std::cout << "Creating linked-list..." << std::endl;
    
    // Create a linked-list representation for the sequence
    auto nodes = createList(seq);

    std::cout << "Linked-list created" << std::endl;
    
    //const auto n_vars = vars.size();

    writeBegin("origGenome.fa", chrID, nodes);
    writeOldBed("origGenome.bed", vars,  chrID);
    
    std::cout << "Generated: origGenome.fa"  << std::endl;
    std::cout << "Generated: origGenome.bed" << std::endl;
    
    //checkSNP(nodes, seq, vars);

    //FastaAlternateReferenceMaker(nodes, seq, vars);

    // Length of the flipped sequence
    const auto n = count(nodes);

    std::cout << "Number of bases in the flipped sequence: " << n << std::endl;

    /*
     * Generating a FASTA file for the variant genome.
     */
    
    writeBegin("/Users/tedwong/Desktop/Ira/varGenome.fa", chrID, nodes);
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
            case VarType::SNP:
            {
                // What's the position of the SNP after reversing? It's n-(i+1)
                var.second.newPos = n-pos;
                
                break;
            }

            case VarType::Insertion:
            {
                // Number of characters inserted
                //const auto offset = var.second.alt.size() - var.second.ref.size();
                
                var.second.newPos = n-pos;
                
                std::cout << var.second.newPos << std::endl;
                
                break;
            }

            case VarType::Deletion:
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
    
    writeEnd("/Users/tedwong/Desktop/Ira/flipGenome.fa", chrID, nodes);
    std::cout << "Generated: flipGenome.fa" << std::endl;

    /*
     * Generating VCF file for the flipped genome
     */
    
    writeBed("/Users/tedwong/Desktop/Ira/flipGenome.bed", vars, chrID);
    std::cout << "Generated: flipGenome.bed" << std::endl;
    
    const auto end = clock();
    std::cout << "Completed in: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;
    
    return 0;
}