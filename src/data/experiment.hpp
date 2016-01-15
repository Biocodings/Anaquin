#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <set>
#include "data/locus.hpp"

namespace Anaquin
{
    /*
     * Meta-data for the experiment.
     */
    
    class Experiment
    {
        public:
        
            typedef unsigned Factor;
            typedef std::vector<Factor> Replicates;

            Experiment(const std::string &);

            // Return the meta-data for the replicates
            inline const Replicates & reps() const { return _reps; }
        
            // Number of factors, typically control vs treated (ie: 2)
            inline std::size_t countFactors() const { return _factors.size(); }

            // Return indexes for a particular factor
            std::vector<std::size_t> factor(Factor) const;

        private:
            // For all replicates
            Replicates _reps;

            // Unique factors
            std::set<Factor> _factors;
    };
}

#endif