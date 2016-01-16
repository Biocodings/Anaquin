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
            typedef std::vector<std::string> Names;

            Experiment() {}

            Experiment(const std::string &factors, const std::string &names)
            {
                addNames(names);
                addFactors(factors);
            }
        
            // Eg: A1,A2,A3,B1,B2,B3
            void addNames(const std::string &);
        
            // Eg: 1,1,1,2,2,2
            void addFactors(const std::string &);
        
            // Return the meta-data for the replicates
            inline const Replicates & reps() const { return _reps; }

            // Return the meta-names for the replicates
            inline const Names & names() const { return _names; }

            // Number of factors, typically control vs treated (ie: 2)
            inline std::size_t countFactors() const { return _factors.size(); }

            // Return indexes for a particular factor
            std::vector<std::size_t> factor(Factor) const;

        private:
            // For all replicates
            Replicates _reps;
        
            // Eg: A1,A2,A3,B1,B2,B3
            Names _names;

            // Unique factors
            std::set<Factor> _factors;
    };
}

#endif