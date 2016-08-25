#ifndef SS_Variable_HPP
#define SS_Variable_HPP

#include <map>
#include <set>
#include <ss/matrix.hpp>
#include <ss/data/data.hpp>
#include <ss/internal/model.hpp>

namespace SS
{
    /*
     * This class represents a variable governed by a random process, including numerical and categorical
     * data type. The class is only designed for internally.
     */
    
    class Variable
    {
        public:
        
            typedef std::set<Level> Levels;

            /*
             * Create a new numerical variable from an iterator.
             */
        
            template <typename Iter> static Variable createNumeric(const Iter &x)
            {
                Variable v;
                
                v._type = Type::Numeric;
                v._data = std::vector<Real> { std::begin(x), std::end(x) };

                return v;
            }

            /*
             * Create a new factor variable.
             */

            static Variable createFactor(const std::vector<Factor> &x)
            {
                Variable v;
            
                v._type   = Type::Categorical;
                v._levels = Levels(x.begin(), x.end());
                
                /*
                 * Convert the factors into simpler representation. For example, we'd convert
                 *
                 *      { "A", "B", "C" } to { 0, 1, 2 }
                 */
            
                v._data.resize(x.size());
                v._factors.resize(x.size());
                
                std::map<Factor, Level> mapping;

                Factor k = 1;
                
                for (auto &level : v._levels)
                {
                    mapping[level] = k++;
                }
                
                for (auto i = 0; i < x.size(); i++)
                {
                    v._data[i] = v._factors[i] = mapping[x[i]];
                }
                
                return v;
            }
        
            // Returns whether this is a numerical variable
            inline bool isNumeric() const { return _type == Type::Numeric; }

            // Returns whether this is a factor variable
            inline bool isFactor() const { return _type == Type::Categorical; }
        
            // Returns the factor levels. Nothing is returned unless this is a categorical variable.
            inline Levels levels() const { return _levels; }
        
            // Number of observations
            inline Counts size() const { return _data.size(); }
        
            // Returns the samples for this variable
            inline const Real * data() const { return ((Real *) &_data[0]); }

            inline double * data__() const
            {
                Matrix m;
                Internal::contrTreatment(m, _levels.size());
                
                double *a = new double[1234];
                Eigen::Map<Eigen::MatrixXd>(a, m.rows(), m.cols()) = m;
                
                return a;
            }
        
            // Returns the samples for this variable
            inline Real * data() { return ((Real *) &_data[0]); }

            // Returns the samples for this variable
            inline const Factor * factors() const { return ((Factor *) &_factors[0]); }
        
            // Returns the samples for this variable
            inline Factor * factors() { return ((Factor *) &_factors[0]); }
        
        private:
        
            enum class Type
            {
                Numeric,
                Categorical,
            };
        
            Type _type;
        
            std::vector<Real> _data;
        
            // Only defined for categorical variable
            Levels _levels;

            // Only defined for a categorical variable
            std::vector<Factor> _factors;
    };
}

#endif