#include "stats/linear.hpp"
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
        
        lm.F     = m.f;
        lm.p     = m.p;
        lm.r     = SS::cor(x, y);
        lm.c     = m.coeffs[0].value;
        lm.m     = m.coeffs[1].value;
        lm.R2    = m.r2;
        lm.ar2   = m.ar2;
        lm.SST   = m.total.ss;
        lm.SSM   = m.model.ss;
        lm.SSE   = m.error.ss;
        lm.SST_D = m.total.df;
        lm.SSM_D = m.model.df;
        lm.SSE_D = m.error.df;
    }
    catch(...)
    {
        lm.F     = NAN;
        lm.p     = NAN;
        lm.r     = NAN;
        lm.c     = NAN;
        lm.m     = NAN;
        lm.R2    = NAN;
        lm.ar2   = NAN;
        lm.SST   = NAN;
        lm.SSM   = NAN;
        lm.SSE   = NAN;
        lm.SST_D = NAN;
        lm.SSM_D = NAN;
        lm.SSE_D = NAN;
    }

    return lm;
}
