#include <set>
#include <iostream>
#include "ParserGTF.hpp"
#include "AssemblyAnalyst.hpp"
#include "StandardFactory.hpp"

using namespace std;

AssemblyStats AssemblyAnalyst::analyze(const std::string &file)
{
	AssemblyStats stats;
//    std::set<Position> pos;
//    
//	/*
//	 * Read for the assembled transcript
//	 */
//
//	struct TranscriptReader : public FeatureReader
//	{
//        void exon(const Feature &f)
//        {
//            pos->insert(f.pos);
//        }
//        
//        std::set<Position> *pos;
//	};
//
//	TranscriptReader tReader;
//    tReader.pos = &pos;
//	ParserGTF::parse(file, tReader);
//
//    Reads tp = 0;
//    Reads tn = 0;
//    Reads fp = 0;
//    Reads fn = 0;
//    
//    std::set<Position> pos_2;
//    
//	/*
//	 * Check for the features in the in-sillico chromosome, have they been assembled?
//	 */
//
//	struct SillicoReader : public FeatureReader
//	{
//		void exon(const Feature &f)
//		{
//            if (pos->count(f.pos))
//            {
//                // It's true-positive because the feature has been assembled
//                (*tp)++;
//            }
//            else
//            {
//                // It's false-negative because the feature has not been assembled
//                (*fn)++;
//            }
//            
//            (*pos_2).insert(f.pos);
//		}
//
//        Reads *tp;
//        Reads *fn;
//        std::set<Position> *pos_2;
//        
//        std::set<Position> *pos;
//		AssemblyStats *stats;
//	};
//
//	SillicoReader sReader;
//	sReader.stats = &stats;
//    sReader.tp = &tp;
//    sReader.fn = &fn;
//    sReader.pos = &pos;
//    sReader.pos_2 = &pos_2;
//    
//	// Check the features listed in the sillico transcripts
//	ParserGTF::parse(StandardFactory::transGTF(), sReader);
//
//
//    
//    struct TranscriptReader_2 : public FeatureReader
//    {
//        void exon(const Feature &f)
//        {
//            if (f.id == "chrT")
//            {
//                if (!pos_2->count(f.pos))
//                {
//                    // It's false-positive because the feature has been assembled for the chromosome but the information is incorrect
//                    (*fp)++;
//                }
//            }
//            else
//            {
//                (*tn)++; // ???
//            }
//        }
//        
//        Reads *tn;
//        Reads *fp;
//        std::set<Position> *pos_2;
//    };
//    
//    TranscriptReader_2 tReader_2;
//    
//    tReader_2.tn = &tn;
//    tReader_2.fp = &fp;
//    tReader_2.pos_2 = &pos_2;
//    
//    ParserGTF::parse(file, tReader_2);
//
//    
//    
//    std::cout << tn << std::endl;
//    std::cout << fp << std::endl;
//    
//    
    
return AssemblyStats();
}