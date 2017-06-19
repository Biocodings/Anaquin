#ifndef PARSER_VCF2_HPP
#define PARSER_VCF2_HPP

#include "data/vData.hpp"
#include "tools/tools.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "data/biology.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVCF2
    {
        enum Field
        {
            Chrom   = 0,
            Pos     = 1,
            ID      = 2,
            Ref     = 3,
            Alt     = 4,
            Qual    = 5,
            Filt    = 6,
            Info    = 7,
            Format  = 8,
            Sample1 = 9,
            Sample2 = 10,
        };
        
        typedef std::function<void (Variant &)> Functor;

        static void parse(const Reader &r, Functor f);
//        {
//            std::string line;
//            
//            std::vector<std::string> tmp;
//            std::vector<std::string> fields;
//            std::vector<std::string> formats;
//            
//            while (r.nextLine(line))
//            {
//                Variant x;
//                
//                if (line[0] == '#')
//                {
//                    continue;
//                }
//                
//                Tokens::split(line, "\t", fields);
//                
//                x.opts["QUAL"]   = fields[Field::Qual];
//                x.opts["INFO"]   = fields[Field::Info];
//                x.opts["FILTER"] = fields[Field::Filt];
//                
//                x.cID = standChr(fields[Field::Chrom]);
//                
//                // Eg: D_1_3_R
//                x.name = fields[Field::ID];
//                
//                // VCF has 1-based position
//                x.l.start = x.l.end = stod(fields[Field::Pos]);
//                
//                // Reference allele
//                x.ref = fields[Field::Ref];
//                
//                // Anything but not filtering is "PASS"
//                x.filter = fields[Field::Filt] == "PASS" ? Filter::Pass : Filter::NotFilted;
//                
//                if (fields[Field::Info] != ".")
//                {
//                    static std::vector<std::string> toks;
//                    Tokens::split(fields[Field::Info], ";", toks);
//                    std::map<std::string, std::string> infos;
//                    
//                    for (const auto &tok : toks)
//                    {
//                        Tokens::split(tok, "=", tmp);
//                        
//                        if (tmp.size() == 2)
//                        {
//                            infos[tmp[0]] = tmp[1];
//                            
//                            // Measured allele frequency
//                            if (tmp[0] == "AF") { x.allF = stof(tmp[1]); }
//                            
//                            else if (tmp.size() == 2)
//                            {
//                                x.opts[tmp[0]] = tmp[1];
//                            }
//                        }
//                    }
//
//                    if (infos.count("SVTYPE") && infos.count("END"))
//                    {
//                        x.l.end = stol(infos["END"]);
//                    }
//                }
//
//                std::vector<Sequence> alts;
//                Tokens::split(fields[Field::Alt], ",", alts);
//                
//                // Ignore anything that is not really a variant
//                if ((x.alt = alts[0]) == ".")
//                {
//                    continue;
//                }
//                
//                x.qual = fields[Field::Qual] != "." ? s2d(fields[Field::Qual]) : NAN;
//
//                if (fields.size() > Field::Format)
//                {
//                    Tokens::split(fields[Field::Format], ":", formats);
//                    
//                    auto hack = fields.size() > Field::Sample2 ? Field::Sample2 : Field::Sample1;
//                    
//                    // Eg: 1/2:0,11,5:16:99:694,166,119,378,0,331
//                    Tokens::split(fields[hack], ":", tmp);
//                    
//                    // Check all the format data...
//                    for (auto j = 0; j < tmp.size(); j++)
//                    {
//                        if (formats[j] == "GT")
//                        {
//                            x.gt = (tmp[j] == "1/1" || tmp[j] == "0/0") ? Genotype::Homozygous : Genotype::Heterzygous;
//                        }
//                        else if (formats[j] == "AF")
//                        {
//                            x.allF = stof(tmp[j]);
//                        }
//                        else if (formats[j] == "AD")
//                        {
//                            std::vector<std::string> toks;
//                            Tokens::split(tmp[j], ",", toks);
//                            
//                            if (toks.size() == 1)
//                            {
//                                x.readR = s2d(toks[0]);
//                                x.readV = s2d(toks[0]);
//                            }
//                            else
//                            {
//                                x.readR = s2d(toks[0]);
//                                x.readV = s2d(toks[1]);
//                            }
//                        }
//                        else if (formats[j] == "DP")
//                        {
//                            x.depth = s2d(tmp[j]);
//                        }
//                    }
//                }
//                
//                f(x);
//            }
//        }
    };

    template <typename F, typename T = VData> T parseVCF2(const Reader &r, F f)
    {
        T t;
        
        ParserVCF2::parse(r, [&](const Variant &x)
        {
            t[x.cID].b2v[x.l.start] = x;
            t[x.cID].m2v[x.type()].insert(x);
            f(x);
        });
        
        return t;
    }
}

#endif
