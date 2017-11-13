#include "VarQuin/v_germ.hpp"
#include "VarQuin/v_somatic.hpp"
#include "VarQuin/v_mutation.hpp"

using namespace Anaquin;

void VMutation::report(const FileName &endo, const FileName &seqs, const VMutation::Options &o)
{
    if (o.isGerm)
    {
        VGerm::Options o_;
        
        o_.work = o.work;
        o_.logger = o.logger;
        o_.filter = o.filter;
        o_.writer = o.writer;
        o_.output = o.output;

        // Germline mutation caller
        VGerm::report(endo, seqs, o_);
    }
    else
    {
        VSomatic::Options o_;
        
        o_.work = o.work;
        o_.logger = o.logger;
        o_.filter = o.filter;
        o_.writer = o.writer;
        o_.output = o.output;

        // Somatic mutation caller
        VSomatic::report(endo, seqs, o_);
    }
}
