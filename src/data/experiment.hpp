/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include <set>
#include <map>
#include <vector>
#include <assert.h>
#include "data/types.hpp"

namespace Anaquin
{
    typedef std::vector<std::string> SampleNames;
    
    class Experiment
    {
        public:
        
            typedef unsigned Factor;
            typedef std::vector<Factor> Samples;

            struct Sample
            {
                // Paired-end files
                FileName file1, file2;
            
                Factor fact;
            };

            Experiment(const std::vector<std::string> &file1,
                       const std::vector<std::string> &file2,
                       const std::vector<std::string> &facts) : _file1(file1), _file2(file2)
            {
                if (_file1.empty())
                {
                    throw std::runtime_error("Empty experimental files");
                }
                else if (_file1.size() != _file2.size())
                {
                    throw std::runtime_error("Paired-end files assumed");
                }
                
                std::map<std::string, Factor> m;
                
                for (const auto &fact : facts)
                {
                    if (!m.count(fact))
                    {
                        const int tmp = m.size();
                        m[fact] = tmp;
                    }
                    
                    _samples.push_back(m[fact]);
                    _factors.insert(m[fact]);
                }

                assert(_samples.size() == _file1.size());
            }

            inline std::size_t countFacts() const { return _factors.size(); }

            inline Sample sample(unsigned i)
            {
                Sample samp;
            
                samp.file1 = _file1[i];
                samp.file2 = _file2[i];
                samp.fact  = _samples[i];

                return samp;
            }
        
        private:
        
            // Factors for the files
            Samples _samples;
        
            // Eg: Paired-end files
            SampleNames _file1, _file2;

            // Unique factors
            std::set<Factor> _factors;
    };
}

#endif