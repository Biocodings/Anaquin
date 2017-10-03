#include "VarQuin/v_process.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

template <typename T, typename F> void parse(const FileName &file, T o, F f)
{
    // Required for pooling paired-end reads
    std::map<ReadName, ParserBAM::Data> seenMates;

    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p.i && !(i.p.i % 1000000))
        {
            o.wait(std::to_string(i.p.i));
        }
        
        if (!x.mapped)
        {
//            __stats__.nNA++;
        }
        else if (true /*__impl__->isReverse(x.cID)*/)
        {
//            __stats__.nSeqs++; // Reverse genome
        }
        else
        {
//            __stats__.nEndo++; // Forward genome
        }
        
        if (!x.isPassed || x.isSecondary || x.isSupplement)
        {
            return;
        }
        
        A_CHECK(x.isPaired, x.name + " is not pair-ended. Singled-ended not supported.");
        
        if (!seenMates.count(x.name))
        {
            seenMates[x.name] = x;
        }
        else
        {
//            auto &seen  = seenMates[x.name];
//            auto first  = seen.isFirstPair ? &seen : &x;
//            auto second = seen.isFirstPair ? &x : &seen;
//
//            /*
//             * Only complement reads aligned to the reverse genome
//             */
//
//            if (__impl__->isReverse(first->cID))
//            {
//                if (first->isForward)
//                {
//                    complement(first->seq);
//                }
//                else
//                {
//                    std::reverse(first->seq.begin(), first->seq.end());
//                }
//            }
//
//            if (__impl__->isReverse(second->cID))
//            {
//                if (second->isForward)
//                {
//                    complement(second->seq);
//                }
//                else
//                {
//                    std::reverse(second->seq.begin(), second->seq.end());
//                }
//            }
//
//            const auto bothRev =  __impl__->isReverse(first->cID) && __impl__->isReverse(second->cID);
//            const auto bothFor = !__impl__->isReverse(first->cID) && !__impl__->isReverse(second->cID);
//            const auto anyRev  =  __impl__->isReverse(first->cID) || __impl__->isReverse(second->cID);
//            const auto anyFor  = !__impl__->isReverse(first->cID) || !__impl__->isReverse(second->cID);
//            const auto anyMap  =  first->mapped ||  second->mapped;
//            const auto anyNMap = !first->mapped || !second->mapped;
//
//            VFlip::Status status;
//
//            if (bothRev && !anyNMap)
//            {
//                status = Status::ReverseReverse;
//            }
//            else if (bothFor && !anyNMap)
//            {
//                status = Status::ForwardForward;
//            }
//            else if (anyRev && anyFor && !anyNMap)
//            {
//                status = Status::ForwardReverse;
//            }
//            else if (anyNMap && anyMap && anyRev)
//            {
//                status = Status::ReverseNotMapped;
//            }
//            else if (anyNMap && anyMap && anyFor)
//            {
//                status = Status::ForwardNotMapped;
//            }
//            else
//            {
//                status = Status::NotMappedNotMapped;
//            }
//
//            __stats__.counts[status]++;
//            __impl__->process(*first, *second, status);
//            seenMates.erase(x.name);
        }
    }, true);

    o.logInfo("Found: " + std::to_string(seenMates.size()) + " unpaired mates.");
}

void VProcess::report(const FileName &file, const Options &o)
{
    parse(file, o, [&]()
    {
              
    });
}
