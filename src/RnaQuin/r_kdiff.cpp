/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include <thread>
#include "data/pachter.hpp"
#include "data/experiment.hpp"
#include "RnaQuin/r_kdiff.hpp"

using namespace Anaquin;

// Used in the worker
static std::vector<FileName> outputs;

// Used for the sleuth output
static FileName __sleuth__;

static void worker(unsigned i, const FileName &index, const FileName &file1, const FileName &file2)
{
    // Run quantification in Kallisto
    outputs[i] = Pachter::externalQuant(index, file1, file2, true);
}

TKDiff::Stats TKDiff::analyze(const std::vector<FileName> &a1,
                              const std::vector<FileName> &a2,
                              const std::vector<FileName> &b1,
                              const std::vector<FileName> &b2,
                              const Options &o)
{
    assert(a1.size() == a2.size());
    assert(b1.size() == b2.size());
    
    std::vector<std::thread> thds;
    outputs.resize(a1.size() + b1.size());
    
    auto wait = [&]()
    {
        for (auto &thd : thds)
        {
            if (thd.joinable())
            {
                thd.join();
            }
        }
    };

    std::vector<std::string> names, facts;
    
    for (int i = 0; i < a1.size(); ++i)
    {
        facts.push_back("A");
        names.push_back("A" + std::to_string(i));
        thds.push_back(std::thread(worker, thds.size(), o.index, a1[i], a2[i]));
    }

    for (int i = 0; i < b1.size(); ++i)
    {
        facts.push_back("B");
        names.push_back("B" + std::to_string(i));
        thds.push_back(std::thread(worker, thds.size(), o.index, b1[i], b2[i]));
    }

    wait();
    
    /*
     * We've run Kallisto for all the samples. Next, we'll need to give them to sleuth.
     */

    RDiff::Options o_;
    
    o_.metrs = RDiff::Metrics::Isoform;
    //o_.dSoft = RDiff::Software::Sleuth;

    return RDiff::analyze(__sleuth__ = Pachter::sleuth(outputs, names, facts), o_);
}

void TKDiff::report(const FileName &file, const Options &o)
{
    std::vector<FileName> a1, a2, b1, b2;
    Experiment::readMeta(file, a1, a2, b1, b2);

    const auto stats = analyze(a1, a2, b1, b2, o);
    const auto units = "isoform";
    
    assert(!__sleuth__.empty());
    o.info("Generating statistics");
    
    /*
     * Generating RnaKDiff_summary.stats
     */
    
    RDiff::generateSummary("RnaKDiff_summary.stats", stats, o, units);
    
    /*
     * Generating RnaKDiff_quins.csv
     */
    
    RDiff::generateCSV("RnaKDiff_quins.csv", stats, o);
    
    /*
     * Generating RnaKDiff_fold.R
     */
    
    //RDiff::generateFoldR("RnaKDiff_fold.R", "RnaKDiff_quins.csv", o);
    
    /*
     * Generating RnaKDiff_ROC.R
     */
    
    //RDiff::generateROC("RnaKDiff_ROC.R", "RnaKDiff_quins.csv", o);
    
    /*
     * Generating RnaKDiff_LODR.R
     */
    
    RDiff::generateLODR("RnaKDiff_LODR.R", "RnaKDiff_quins.csv", o);
}