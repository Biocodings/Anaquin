#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <set>
#include <map>
#include <vector>

namespace Anaquin
{
    typedef std::string SampleName;
    typedef std::vector<std::string> SampleNames;
    
    class CountTable
    {
        typedef unsigned Count;
        typedef std::vector<Count> Counts;

        typedef std::string FeatureID;
        typedef std::vector<std::string> FeatureIDs;

        public:

            CountTable(const SampleNames &names) : _names(names)
            {
                for (const auto &name : _names)
                {
                    _samples[name];
                }
            }

            // Feature IDs
            inline const FeatureIDs & ids() const { return _ids; };

            // Sample names
            inline const SampleNames & names() const { return _names; }

            // Counts for a sample
            inline const Counts & counts(const SampleName &name) const { return _samples.at(name); }

            // Sort the table given a permutation
            void sort(const std::vector<std::size_t> &);

            /*
             * The following functions should be used cautiously as the table can be out-of-synched if isn't
             */
        
            inline void addCount(const SampleName &name, Count c)
            {
                _samples.at(name).push_back(c);
            }

            inline void addFeature(const FeatureID &id)
            {
                _ids.push_back(id);
            }

        private:

            // Eg: A1
            SampleNames _names;

            // Eg: R1_40
            FeatureIDs _ids;

            // Eg: R1_40, 10, 90, 45
            std::map<SampleName, std::vector<Count>> _samples;
    };

    class Experiment
    {
        public:
        
            typedef unsigned Factor;
            typedef std::vector<Factor> Samples;

            Experiment() {}

            Experiment(const std::string &factors, const std::string &names)
            {
                addNames(names);
                addFactors(factors);
            }

            // Return an empty count table suitable for the experiment
            inline CountTable countTable() const
            {
                return CountTable(_names);
            }

            // Eg: A1,A2,A3,B1,B2,B3
            void addNames(const std::string &);
        
            // Eg: 1,1,1,2,2,2
            void addFactors(const std::string &);
        
            // Return the meta-data for the samples
            inline const Samples & samples() const { return _samples; }

            // Return the meta-names for the replicates
            inline const SampleNames & names() const { return _names; }

            // Number of conditions, typically control vs treated (ie: 2)
            inline std::size_t countConds() const { return _factors.size(); }

            // Return the positions for a condition
            std::vector<std::size_t> cond(Factor) const;

        private:
            Samples _samples;
        
            // Eg: A1,A2,A3,B1,B2,B3
            SampleNames _names;

            // Unique factors
            std::set<Factor> _factors;
    };
}

#endif