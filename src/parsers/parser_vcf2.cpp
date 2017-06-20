#include <stdio.h>
#include <stdlib.h>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "parsers/parser_vcf2.hpp"

using namespace Anaquin;

void ParserVCF2::parse(const Reader &r, Functor f)
{
    htsFile *fp = bcf_open(r.src().c_str(), "r");
    
    if (!fp)
    {
        throw std::runtime_error("Failed to open: " + r.src());
    }
    
    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    if (!hdr)
    {
        throw std::runtime_error("Failed to open: " + r.src());
    }

    int c;
    char pass[] = "PASS";

    char  *ps = (char  *) malloc(1024 * sizeof(int));
    float *pf = (float *) malloc(1024 * sizeof(float));
    
    bcf1_t *line = bcf_init();
    
    while (bcf_read(fp, hdr, line) == 0)
    {
        Variant x;
        
        bcf_unpack(line, BCF_UN_ALL);
        
        x.cID  = chrom(std::string(bcf_seqname(hdr, line)));
        x.name = line->d.id;
        
        x.l.start = x.l.end = line->pos;
        
        x.ref    =  std::string(line->d.allele[0]);
        x.alt    =  std::string(line->d.allele[1]);
        x.filter = (bcf_has_filter(hdr, line, pass) == 1) ? Filter::Pass : Filter::NotFilted;

        auto info1 = [&](const std::string &key)
        {
            if (bcf_get_info_string(hdr, line, key.c_str(), &ps, &c) > 0)
            {
                x.opts[key] = ps;
            }
        };
        
        auto info2 = [&](const std::string &key)
        {
            if (bcf_get_info_float(hdr, line, key.c_str(), &pf, &c) > 0)
            {
                x.opts[key] = std::to_string(*pf);
            }
        };
        
        info1("CS");
        info1("CX");
        info1("GN");
        info1("GT");
        info2("DP");
        info2("AF");
        info2("CP");
        
        if (x.opts.count("AF")) { x.allF = stod(x.opts.at("AF")); }
        
        int32_t *gt = NULL, ngt_arr = 0;
        if (bcf_get_genotypes(hdr, line, &gt, &ngt_arr) == 2)
        {
            x.gt = gt[0] == gt[1] ? Genotype::Homozygous : Genotype::Heterzygous;
            free(gt);
        }
        
        x.qual = line->qual;
        
        f(x);
    }
    
    hts_close(fp);
}
