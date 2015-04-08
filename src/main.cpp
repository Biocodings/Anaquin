#include <iostream>
#include <tclap/CmdLine.h>
#include "aligner.hpp"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

using namespace Spike;

/*
 bam_hdr_t *bam_hdr_init(void);
 bam_hdr_t *bam_hdr_read(BGZF *fp);
 int bam_hdr_write(BGZF *fp, const bam_hdr_t *h);
 void bam_hdr_destroy(bam_hdr_t *h);
 int bam_name2id(bam_hdr_t *h, const char *ref);
 bam_hdr_t* bam_hdr_dup(const bam_hdr_t *h0);
 
 bam1_t *bam_init1(void);
 void bam_destroy1(bam1_t *b);
 int bam_read1(BGZF *fp, bam1_t *b);
 int bam_write1(BGZF *fp, const bam1_t *b);
 bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc);
 bam1_t *bam_dup1(const bam1_t *bsrc);
 
 int bam_cigar2qlen(int n_cigar, const uint32_t *cigar);
 int bam_cigar2rlen(int n_cigar, const uint32_t *cigar);
*/

int main(int argc, char ** argv)
{
    typedef TCLAP::SwitchArg SArg;
    typedef TCLAP::ValueArg<std::string> VArg;

	std::cout << "Anquin by Garvan Institute - sequin data-analysis tool\n" << std::endl;

    try
    {
        TCLAP::CmdLine cmd("", ' ', "1.0");
        
        SArg t("t", "test", "Internal testing", cmd, false);
        VArg a("a", "align", "Assess alignment", false, "", "string", cmd);

        cmd.parse(argc, argv);

        if (t.getValue())
        {
            return Catch::Session().run(1, argv);
        }
        else if (!a.getValue().empty())
        {
            Aligner::analyze(a.getValue());
        }
    }
    catch (TCLAP::ArgException &e)
    {
        throw;
    }
    
    return 0;
}