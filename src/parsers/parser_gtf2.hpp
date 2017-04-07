#ifndef PARSER_GTF2_HPP
#define PARSER_GTF2_HPP

#include <map>
#include <cstring>
#include "data/locus.hpp"
#include "data/reader.hpp"
#include "data/convert.hpp"
#include "data/biology.hpp"
#include "tools/errors.hpp"

namespace Anaquin
{
    class ParserGTF2
    {
        public:
        
            struct Data
            {
                ChrID cID;
            
                // Forward or reverse strand?
                Strand str;
            
                // The location of the feature relative to the chromosome
                Locus l;
            
                RNAFeature type;
            
                // Empty if the information is unavailable
                GeneID gID;
            
                // Empty if the information is unavailable
                TransID tID;
            
                // Available only for Cufflink
                double fpkm = NAN;
            };
        
            template <typename F> static void parse(const Reader &r, F f)
            {
                std::map<std::string, RNAFeature> mapper =
                {
                    { "exon",       RNAFeature::Exon },
                    { "gene",       RNAFeature::Gene },
                    { "transcript", RNAFeature::Transcript }
                };
                
                FILE *file = fopen(r.src().c_str(), "r");
                
                Data x;
                Action act;
                ParserProgress pg;
                Phase p = Phase::SeqName;
                
                // This should be large enough...
                char mem[1024];

                // We'll need mem[0] for decremental
                auto buf = &(mem[1]);
                
                double *attrD = nullptr;
                std::string *attrS = nullptr;

                Counts i = 0;
                
                auto prepareAttr = [&]()
                {
                    if (!strcmp(buf, "gene_id"))
                    {
                        attrS = &x.gID;
                    }
                    else if (!strcmp(buf, "transcript_id"))
                    {
                        attrS = &x.tID;
                    }
                    else if (!strcmp(buf, "FPKM"))
                    {
                        attrD = &x.fpkm;
                    }
                    else
                    {
                        attrS = nullptr;
                        attrD = nullptr;
                    }
                };

                while ((act = readNext(file, buf)) != Action::EndOfFile)
                {
                    switch (act)
                    {
                        case Action::NewLine:
                        {
                            assert(!x.cID.empty());
                            
                            f(x, pg);
                            pg.i++;
                            p = Phase::SeqName;
                            break;
                        }

                        case Action::NewTab:
                        {
                            switch (p)
                            {
                                case Phase::Attr:
                                case Phase::Frame:
                                case Phase::Score:
                                case Phase::Source:  { break; }
                                case Phase::SeqName: { x.cID = buf; break; }
                                case Phase::Feature: { x.type = mapper[buf];  break; }
                                case Phase::Start:   { x.l.start = std::stoi(buf); break; }
                                case Phase::End:     { x.l.end   = std::stoi(buf); break; }

                                case Phase::Strand:
                                {
                                    switch (buf[0])
                                    {
                                        case '+': { x.str = Forward;  break; }
                                        case '-': { x.str = Backward; break; }
                                        case '.': { x.str = Either;   break; }
                                        default:  { A_THROW("Invalid strand: " + std::string(buf)); }
                                    }
                                    
                                    break;
                                }
                            }

                            p = static_cast<Phase>(p+1);
                            break;
                        }

                        case Action::NewSpace:
                        {
                            // Skip "; "
                            if (strlen(buf))
                            {
                                prepareAttr();
                            }
                            
                            break;
                        }
                            
                        case Action::NewColon:
                        {
                            A_CHECK(p == Phase::Attr, "Invalid line: " + std::to_string(i));
                            
                            if (attrS)
                            {
                                *attrS = buf;
                            }
                            else if (attrD)
                            {
                                *attrD = s2d(buf);
                            }

                            attrS = nullptr;
                            attrD = nullptr;

                            break;
                        }

                        case Action::EndOfFile: { return; }
                    }
                    
                    i++;
                }
            }

        private:
        
            enum Phase
            {
                SeqName,
                Source,
                Feature,
                Start,
                End,
                Score,
                Strand,
                Frame,
                Attr
            };
        
            enum class Action
            {
                NewLine,
                NewTab,
                NewSpace,
                NewColon,
                EndOfFile,
            };
        
            static Action readNext(FILE *file, char *buf)
            {
                static const auto c2a = std::map<char, Action>
                {
                    { '\t', Action::NewTab    },
                    { '\n', Action::NewLine   },
                    { ' ',  Action::NewSpace  },
                    { ';',  Action::NewColon  },
                    { EOF,  Action::EndOfFile },
                };

                // Needed for comments
                auto old = buf;

                auto f = [&]()
                {
                    auto c = *buf;
                    *(buf) = '\0';
                    return c2a.at(c);
                };
                
                for (; (*buf = (char)fgetc(file)); buf++)
                {
                    switch (*buf)
                    {
                        case '#':
                        {
                            while ((*buf = (char)fgetc(file)))
                            {
                                if (*buf == EOF || *buf == '\n')
                                {
                                    break;
                                }
                            }
                            
                            // Throw out everything we've read (includes anything before the comment)
                            buf = old-1;

                            break;
                        }

                        case EOF:
                        case ' ':
                        case ';':
                        case '\t':
                        case '\n':
                        {
                            return f();
                        }

                        case '\"': { buf--; break; }
                    }
                }

                return Action::EndOfFile;
            }
    };
}

#endif
