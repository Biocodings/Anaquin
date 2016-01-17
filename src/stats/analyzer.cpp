#include "stats/analyzer.hpp"
#include <ss/regression/segmented.hpp>

using namespace Anaquin;

void LinearStats::data(std::vector<double> &x, std::vector<double> &y, bool shouldLog, std::vector<std::string> *ids) const
{
    x.clear();
    y.clear();
    
    if (ids)
    {
        ids->clear();
    }
    
    auto f = [&](double v)
    {
        // Don't log zero...
        assert(!shouldLog || v);
        
        return shouldLog ? (v ? log2(v) : 0) : v;
    };
    
    for (const auto &p : *this)
    {
        if (!isnan(p.second.x) && !isnan(p.second.y))
        {
            x.push_back(f(p.second.x));
            y.push_back(f(p.second.y));
            
            if (ids)
            {
                ids->push_back(p.first);
            }
        }
    }
}

InflectionLimit LinearStats::inflect(bool shouldLog) const
{
    std::vector<std::string> ids;
    std::vector<double> x, y;
    
    data(x, y, shouldLog, &ids);

    const auto r = SS::segmentPieceWise(x, y);
    
    InflectionLimit l;
    
    // The break we're looking for
    l.b = r.b;
    
    l.lR2  = r.bkLR2();
    l.rR2  = r.bkRR2();
    l.lSl  = r.bkLSl();
    l.rSl  = r.bkRSl();
    l.lInt = r.bkLInt();
    l.rInt = r.bkRInt();

    for (auto i = 0; i < ids.size(); i++)
    {
        if (x[i] == l.b)
        {
            l.id = ids[i];
            break;
        }
    }
    
    assert(!l.id.empty());
    
    return l;
}

LinearModel LinearStats::linear(bool shouldLog) const
{
    std::vector<double> x, y;
    data(x, y, shouldLog);
    
    LinearModel lm;
    
    try
    {
        if (std::adjacent_find(x.begin(), x.end(), std::not_equal_to<double>()) == x.end())
        {
            throw std::runtime_error("Failed to perform linear regression. Flat mixture.");
        }
        
        const auto m = SS::lm(SS::R::data.frame(SS::R::c(y), SS::R::c(x)));
        
        lm.f      = m.f;
        lm.p      = m.p;
        lm.r      = SS::cor(x, y);
        lm.c      = m.coeffs[0].value;
        lm.m      = m.coeffs[1].value;
        lm.r2     = m.r2;
        lm.ar2    = m.ar2;
        lm.sst    = m.total.ss;
        lm.ssm    = m.model.ss;
        lm.sse    = m.error.ss;
        lm.sst_df = m.total.df;
        lm.ssm_df = m.model.df;
        lm.sse_df = m.error.df;
    }
    catch(...)
    {
        lm.f      = NAN;
        lm.p      = NAN;
        lm.r      = NAN;
        lm.c      = NAN;
        lm.m      = NAN;
        lm.r2     = NAN;
        lm.ar2    = NAN;
        lm.sst    = NAN;
        lm.ssm    = NAN;
        lm.sse    = NAN;
        lm.sst_df = NAN;
        lm.ssm_df = NAN;
        lm.sse_df = NAN;
    }
    
    return lm;
}
