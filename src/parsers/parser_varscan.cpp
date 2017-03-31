#include <map>
#include "tools/errors.hpp"
#include "parsers/parser_varscan.hpp"
#include <boost/algorithm/string.hpp>

using namespace Anaquin;

namespace Pileup
{
    enum Field
    {
        Chrom,
        Position,
        Ref,
        Cons,
        Reads1,
        Reads2,
        VarFreq,
        Strands1,
        Strands2,
        Qual1,
        Qual2,       // 10
        PValue,
        MapQual1,
        MapQual2,
        Reads1Plus,
        Reads1Minus,
        Reads2Plus,
        Reads2Minus,
        VarAllele    // 18
    };
}

namespace Somatic
{
    enum Field
    {
        Chrom,
        Position,
        Ref,
        Var,
        NormalReads1,
        NormalReads2,
        NormalVarFreq,
        NormalGT,
        TumorReads1,
        TumorReads2,
        TumorVarFreq,
        TumorGT,
        SomaticStatus,
        VariantPValue,
        SomaticPValue,
        TumorReads1Plus,
        TumorReads1Minus,
        TumorReads2Plus,
        TumorReads2Minus,
        NormalReads1Plus,
        NormalReads1Minus,
        NormalReads2Plus,
        NormalReads2Minus
    };
}

bool ParserVarScan::isPileup(const Reader &r)
{
    std::string line;
    std::vector<Token> toks;
    
    if (r.nextLine(line))
    {
        Tokens::split(line, "\t", toks);
        
        if (toks.size() == 19         &&
            toks[0]  == "Chrom"       &&
            toks[1]  == "Position"    &&
            toks[2]  == "Ref"         &&
            toks[3]  == "Cons"        &&
            toks[4]  == "Reads1"      &&
            toks[5]  == "Reads2"      &&
            toks[6]  == "VarFreq"     &&
            toks[7]  == "Strands1"    &&
            toks[8]  == "Strands2"    &&
            toks[9]  == "Qual1"       &&
            toks[10] == "Qual2"       &&
            toks[11] == "Pvalue"      &&
            toks[12] == "MapQual1"    &&
            toks[13] == "MapQual2"    &&
            toks[14] == "Reads1Plus"  &&
            toks[15] == "Reads1Minus" &&
            toks[16] == "Reads2Plus"  &&
            toks[17] == "Reads2Minus" &&
            toks[18] == "VarAllele")
        {
            return true;
        }
    }
    
    return false;
}

bool ParserVarScan::isSomatic(const Reader &r)
{
    std::string line;
    std::vector<Token> toks;
    
    if (r.nextLine(line))
    {
        Tokens::split(line, "\t", toks);
        
        if (toks.size() == 23                 &&
            toks[0]  == "chrom"               &&
            toks[1]  == "position"            &&
            toks[2]  == "ref"                 &&
            toks[3]  == "var"                 &&
            toks[4]  == "normal_reads1"       &&
            toks[5]  == "normal_reads2"       &&
            toks[6]  == "normal_var_freq"     &&
            toks[7]  == "normal_gt"           &&
            toks[8]  == "tumor_reads1"        &&
            toks[9]  == "tumor_reads2"        &&
            toks[10] == "tumor_var_freq"      &&
            toks[11] == "tumor_gt"            &&
            toks[12] == "somatic_status"      &&
            toks[13] == "variant_p_value"     &&
            toks[14] == "somatic_p_value"     &&
            toks[15] == "tumor_reads1_plus"   &&
            toks[16] == "tumor_reads1_minus"  &&
            toks[17] == "tumor_reads2_plus"   &&
            toks[18] == "tumor_reads2_minus"  &&
            toks[19] == "normal_reads1_plus"  &&
            toks[20] == "normal_reads1_minus" &&
            toks[21] == "normal_reads2_plus"  &&
            toks[22] == "normal_reads2_minus")
        {
            return true;
        }
    }
    
    return false;
}

void ParserVarScan::parseSomatic(const Reader &r, Functor f)
{
    using Somatic::Field;
    
    Data d;
    ParserProgress p;
    
    Line line;
    std::vector<Token> toks;
    
    while (r.nextLine(line))
    {
        if (p.i++ == 0)
        {
            continue;
        }
        
        Tokens::split(line, "\t", toks);
        
        if (toks[Field::SomaticStatus] == "Unknown")
        {
            continue;
        }
        
        d.cID = toks[Field::Chrom];
        
        // Always start and end at the same position
        d.l = Locus(stod(toks[Field::Position]), stod(toks[Field::Position]));
        
        d.allF = s2d(toks[Field::TumorVarFreq]);
        
        const std::map<std::string, Variant::Status> s2s =
        {
            { "LOH",      Variant::Status::LOH },
            { "Somatic",  Variant::Status::Somatic },
            { "Germline", Variant::Status::Germline },
        };
        
        if (!s2s.count(toks[Field::SomaticStatus]))
        {
            A_THROW("Unknown variant status: " + toks[Field::SomaticStatus]);
        }
        
        d.status = s2s.at(toks[Field::SomaticStatus]);

        const auto readR = s2d(toks[Field::TumorReads1]);
        const auto readV = s2d(toks[Field::TumorReads2]);
        
        d.readR = readR;
        d.readV = readV;
        d.depth = d.readR + d.readV;
        
        if (!readV)
        {
            continue;
        }
        
        d.ref = toks[Field::Ref];
        d.alt = toks[Field::Var];
        
        try
        {
            d.p = stold(toks[Field::SomaticPValue]);
        }
        catch (...)
        {
            d.p = 0.0;
        }
        
        A_ASSERT(d.p >= 0 && d.p <= 1.0);
        
        f(d, p);
    }
}

void ParserVarScan::parsePile(const Reader &r, Functor f)
{
    using Pileup::Field;
    
    Data d;
    ParserProgress p;
    
    Line line;
    std::vector<Token> toks;
    
    while (r.nextLine(line))
    {
        if (p.i++ == 0)
        {
            continue;
        }
        
        Tokens::split(line, "\t", toks);
        
        d.cID = toks[Field::Chrom];
        
        // Always start and end at the same position
        d.l = Locus(stod(toks[Field::Position]), stod(toks[Field::Position]));
        
        d.allF = s2d(toks[Field::VarFreq]);
        
        d.status = Variant::Status::Germline;
        
        /*
         * TODO: VarScan does give quality but it's not the same as quality in VCF.
         *       VarScan gives quality score for the reference and alternative allele.
         */
        
        const auto readR = s2d(toks[Field::Reads1]);
        const auto readV = s2d(toks[Field::Reads2]);
        
        d.readR = readR;
        d.readV = readV;
        d.depth = d.readR + d.readV;
        
        if (!readV)
        {
            continue;
        }
        
        /*
         * Is this an insertion?
         */
        
        const auto isInsert = toks[Field::Cons].find("*/+") != std::string::npos;
        const auto isDelete = toks[Field::Cons].find("*/-") != std::string::npos;
        
        auto clean = [&](const std::string &s)
        {
            auto x = s;
            
            boost::replace_all(x, "*", "");
            boost::replace_all(x, "/", "");
            boost::replace_all(x, "-", "");
            boost::replace_all(x, "+", "");
            
            return x;
        };
        
        if (isInsert)
        {
            /*
             * Eg: A * /+CT	+CT
             */
            
            d.ref = clean(toks[Field::Ref]);
            d.alt = clean(toks[Field::Ref] + toks[Field::VarAllele]);
        }
        else if (isDelete)
        {
            /*
             * Eg: G * /-CTTCCTCTTTC CTTCCTCTTTC
             */
            
            d.ref = clean(toks[Field::Ref] + toks[Field::VarAllele]);
            d.alt = clean(toks[Field::Ref]);
        }
        else
        {
            d.ref = clean(toks[Field::Ref]);
            d.alt = clean(toks[Field::VarAllele]);
        }
        
        A_ASSERT(d.ref.find('+') == std::string::npos);
        A_ASSERT(d.alt.find('+') == std::string::npos);
        
        try
        {
            d.p = stold(toks[Field::PValue]);
        }
        catch (...)
        {
            d.p = 0.0;
        }
        
        A_ASSERT(d.p >= 0 && d.p <= 1.0);
        
        f(d, p);
    }
}

void ParserVarScan::parse(const Reader &r, Functor f)
{
    if (ParserVarScan::isPileup(r))
    {
        ParserVarScan::parsePile(Reader(r), f);
    }
    else if (ParserVarScan::isSomatic(Reader(r)))
    {
        ParserVarScan::parseSomatic(Reader(r), f);
    }
    else
    {
        throw InvalidFileError("Unknown file type for ParserVarScan");
    }
}
