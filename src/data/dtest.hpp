#ifndef DIFF_TEST_HPP
#define DIFF_TEST_HPP

#include "data/types.hpp"

namespace Anaquin
{
    struct DiffTest
    {
        enum class Status
        {
            Tested,
            NotTested,
        };

        ChromoID cID;
        
        // Subject, eg: GeneID, SequinID or ExonID
        GenericID id;
        
        // Log-fold ratio
        LogFold logF;
        
        // The p-value and q-value under null-hypothesis
        double p, q;

        // Not always available
        //FPKM fpkm_1 = 0;

        // Not always available
        //FPKM fpkm_2 = 0;

        Status status;
    };
}

#endif
