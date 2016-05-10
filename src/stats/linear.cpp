#include "stats/linear.hpp"
#include <ss/regression/segmented.hpp>

using namespace Anaquin;

LinearStats::Data LinearStats::data(bool shouldLog) const
{
    Data r;
    
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
            const auto x = f(p.second.x);
            const auto y = f(p.second.y);
            
            r.x.push_back(x);
            r.y.push_back(y);
            r.ids.push_back(p.first);
            
            r.id2x[p.first] = x;
            r.id2y[p.first] = y;
        }
    }
    
    return r;
}

LOQModel LinearStats::limitQuant(bool shouldLog) const
{
    const auto d = data(shouldLog);
    
    const auto r = SS::segmentPieceWise(d.x, d.y);

    LOQModel l;
    
    // The break we're looking for
    l.b = r.b;
    
    l.lR2  = r.blR2();
    l.rR2  = r.brR2();
    l.lSl  = r.blSl();
    l.rSl  = r.brSl();
    l.lInt = r.blInt();
    l.rInt = r.brInt();

    for (auto i = 0; i < d.ids.size(); i++)
    {
        if (d.x[i] == l.b)
        {
            l.id = d.ids[i];
            break;
        }
    }

    assert(!l.id.empty());
    
    return l;
}

LinearModel LinearStats::linear(bool shouldLog) const
{
    const auto d = data(shouldLog);
    
    LinearModel lm;
    
    try
    {
        if (std::adjacent_find(d.x.begin(), d.x.end(), std::not_equal_to<double>()) == d.x.end())
        {
            throw std::runtime_error("Failed to perform linear regression. Flat mixture.");
        }
        
        const auto m = SS::linearModel(d.y, d.x);

        lm.F     = m.f;
        lm.p     = m.p;
        lm.r     = SS::cor(d.x, d.y);
        lm.c     = m.coeffs[0].est;
        lm.m     = m.coeffs[1].est;
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
