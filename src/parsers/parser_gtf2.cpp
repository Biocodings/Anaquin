//#include <cstdio>
//#include <fstream>
//#include "parsers/parser_gtf2.hpp"
//
//using namespace Anaquin;
//
//////chrIS	Sequin	transcript	6955487	6962816	.	+	.	gene_id "R1_101"; transcript_id "R1_101_1"; gene_type "synthetic"; transcript_type "synthetic"
//void ParserGTF2::parse(const Reader &r)
//{
//    FILE *file = fopen(r.src().c_str(), "r");
//
//    Action action;
//    
//    // This should be large enough...
//    char buf[4096];
//    
//    while ((action = readNext(file, buf)) != Action::EndOfFile)
//    {
//        
//    }
//    
//    
////    char *code = malloc(1000 * sizeof(char));
////    do
////    {
////        *code++ = (char)fgetc(file);
////        
////    } while(*code != EOF);
////    return code;
//}
//
//
////{
////    std::map<std::string, RNAFeature> mapper =
////    {
////        { "exon",       Exon },
////        { "gene",       Gene },
////        { "transcript", Transcript }
////    };
////    
////    std::string line;
////    Data x;
////    
////    ParserProgress p;
////    
////    std::vector<std::string> opts;
////    std::vector<std::string> toks;
////    std::vector<std::string> nameVal;
////    
////    while (r.nextLine(line))
////    {
////        p.i++;
////        boost::split(toks, line, boost::is_any_of("\t"));
////        
////        // Empty line? Unknown feature such as mRNA?
////        if (toks.size() == 1 || !mapper.count(toks[2]))
////        {
////            continue;
////        }
////        
////        x.cID  = toks[0];
////        x.type = mapper[toks[2]];
////        
////        x.l.start = stoi(toks[3]);
////        x.l.end   = stoi(toks[4]);
////        
////        if (toks[6] != "+" && toks[6] != "-" && toks[6] != ".")
////        {
////            throw std::runtime_error("File: " + r.src() + ". Invalid strand: [" + toks[6] + "]. Line: " + line );
////        }
////        
////        if (toks[6] == ".")
////        {
////            x.str = Strand::Unknown;
////        }
////        else
////        {
////            x.str = toks[6] == "+" ? Strand::Forward : Strand::Backward;
////        }
////        
////        /*
////         * Eg: "gene_id "R_5_3"; transcript_id "R_5_3_R";"
////         */
////        
////        boost::split(opts, toks[8], boost::is_any_of(";"));
////        
////        for (auto option : opts)
////        {
////            if (!option.empty())
////            {
////                boost::trim(option);
////                boost::split(nameVal, option, boost::is_any_of(" "));
////                
////                if (nameVal.size() == 2)
////                {
////                    // Make sure the silly characters are removed
////                    nameVal[1].erase(std::remove(nameVal[1].begin(), nameVal[1].end(), '\"'), nameVal[1].end());
////                    
////                    if (nameVal[0] == "gene_id")
////                    {
////                        x.gID = nameVal[1];
////                    }
////                    else if (nameVal[0] == "transcript_id")
////                    {
////                        x.tID = nameVal[1];
////                    }
////                    else if (nameVal[0] == "FPKM")
////                    {
////                        x.fpkm = s2d(nameVal[1]);
////                    }
////                }
////            }
////        }
////        
////        f(x, line, p);
////    }
