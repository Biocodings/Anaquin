#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include "parsers/parser_fa.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

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
    }
    
    Base newPos = 0;
};

struct Node
{
    Node(char base) : base(base) {}
    
    char base;

    Node *next = nullptr;
    Node *prev = nullptr;
};

std::pair<Node *, Node *> createList(const std::string &str)
{
    Node *start = new Node(str[0]);
    Node *prev  = start;
    
    for (auto i = 1; i < str.length(); i++)
    {
        Node *node = new Node(str[i]);
     
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
    
    while (node)
    {
        f << node->base;
        node = node->next;
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

void FastaAlternateReferenceMaker(const std::pair<Node *, Node *> &startEnd, const std::string &str, std::map<Base, ModifedVCF> &vars)
{
    auto i = str.length() - 1;
    auto node = startEnd.second;
    
    while (node)
    {
        if (vars.count(i))
        {
            const auto &var = vars.at(i);
            const auto diff = var.diff();
            
            switch (var.type())
            {
                case VarType::SNP:
                {
                    assert(var.alt.size() == 1 && node->base == var.ref[0]);
                    node->base = var.alt[0];
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
                    
                    break;
                }
            }
        }
        
        i--;
        node = node->prev;
    }
    
    std::map<Base, ModifedVCF> x = vars;
    
    // We'll construct the variants again...
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
        }
    });
    
    std::cout << "Number of variants: " + std::to_string(vars.size()) << std::endl;
    assert(!vars.empty());
}

static void writeVCF(const std::string &file, const std::map<Base, ModifedVCF> &vars)
{
    std::ofstream f;
    f.open(file);

    f << "##CHROM POS     ID        REF ALT\n" << std::endl;

    for (const auto &var : vars)
    {
        f << (boost::format("%1%\t%2%\t%3%\t%4%\t%5%\n") % "chrT"
                                                         % var.second.newPos
                                                         % var.second.l.start
                                                         % var.second.ref
                                                         % var.second.alt).str();
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
 *
 *    g++ -std=c++11 -c -I /usr/include/boost -I ~/Sources/QA/src -I ~/Sources/SS main.cpp
 *    g++ -std=c++11 -c -I /usr/include/boost -I ~/Sources/QA/src ~/Sources/QA/src/data/reader.cpp
 *    g++ *.o -o seqTool
 */

int main(int argc, const char * argv[])
{
    std::string seq;
    std::map<Base, ModifedVCF> vars;

    clock_t begin = clock();
    
    const auto gFile = argv[1];
    const auto vFile = argv[2];
    const auto chrID = argv[3];

    std::cout << "Chromosome: " << chrID << std::endl;
    
    std::cout << "Reading genome..." << std::endl;
    readFA(gFile, chrID, seq);
    std::cout << "Reading genome completed" << std::endl;

    std::cout << "Reading variants..." << std::endl;
    readVCF(vFile, chrID, vars);
    std::cout << "Reading variants completed" << std::endl;

    auto test = createTestData();
    seq  = test.first;
    vars = test.second;

    // Create a linked-list representation for the sequence
    auto nodes = createList(seq);
    
    const auto n_vars = vars.size();
    
    FastaAlternateReferenceMaker(nodes, seq, vars);

    if (vars.size() != n_vars)
    {
        throw std::runtime_error("vars.size() != n_vars");
    }
    
    // Length of the reference sequence
    const auto n = seq.length();
    
    std::cout << "Number of bases: " << n << std::endl;

    /*
     * Generating list of variants after flipping. This is only possible because we know the size of the chromosome.
     */

    for (auto &var : vars)
    {
        const auto pos = var.second.l.start;
        
        switch (var.second.type())
        {
            case VarType::SNP:
            {
                // What's the position of the SNP after reversing? It's n-(i+1)
                var.second.newPos = n - (pos+1);
                
                break;
            }
                
            case VarType::Insertion:
            {
                // Number of characters inserted
                const auto offset = var.second.alt.size() - var.second.ref.size();
                
                var.second.newPos = n-(pos+1)-offset;
                break;
            }

            case VarType::Deletion:
            {
                var.second.newPos = n-(pos+1);
                break;
            }
        }
    }
 
    /*
     * Now, flip the genome... This is easy with a linked-list implementation... In fact we don't have to do anything...
     */

    /*
     * Generating a FASTA file for the variant genome.
     */
    
    writeBegin("genome_with_vars.fa", chrID, nodes);
    std::cout << "Generated FASTA for modified genome: genome_with_vars.fa" << std::endl;
    
    /*
     * Generating a FASTA file for the flipped genome
     */
    
    writeEnd("flipped_genome_with_vars.fa", chrID, nodes);
    std::cout << "Generated FASTA for flipped genome: flipped_genome_with_vars.fa" << std::endl;

    /*
     * Generating VCF file for the flipped genome
     */
    
    writeVCF("flipped_vars.vcf", vars);
    std::cout << "Generated VCF for flipped genome: flipped_vars.vcf" << std::endl;
    
    clock_t end = clock();
    std::cout << "Completed in: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;
    
    return 0;
}