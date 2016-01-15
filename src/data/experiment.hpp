#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

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
            typedef std::vector<Factor> Factors;

            Experiment(const std::string &);
        
            inline const Factors & factors() const { return _factors; }
        
        private:
            Factors _factors;
    };
}

#endif