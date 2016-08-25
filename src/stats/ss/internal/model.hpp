#ifndef SS_INTERNAL_MODEL_HPP
#define SS_INTERNAL_MODEL_HPP

#include <numeric>
#include <ss/matrix.hpp>
#include <ss/data/errors.hpp>
#include <ss/internal/vargs.hpp>

namespace SS
{
    namespace Internal
    {
        inline void contrTreatment(Matrix &m, std::size_t n)
        {
            Eigen::VectorXd v(n);
            
            for (auto i = 0; i < n; i++)
            {
                v(i) = 1;
            }
            
            // Diagonse the matrix and remove the first column
            m = Matrix(v.asDiagonal()).block(0, 1, n, n-1);
        }
        
        /*
         * Adopted from modelMatrix.cpp in R
         */
        
        inline void firstFactor(Real *x, long nrx, long ncx, Real *c, long nrc, long ncc, unsigned *v)
        {
            Real *cj, *xj;
            
            for (int j = 0; j < ncc; j++)
            {
                xj = &x[j * nrx];
                cj = &c[j * nrc];
                
                for (int i = 0; i < nrx; i++)
                {
                    xj[i] = cj[v[i]-1];
                }
            }
        }
        
        /*
         * Adopted from modelMatrix.cpp in R
         */
        
        inline void addFactor(Real *x, long nrx, long ncx, Real *c, long nrc, long ncc, unsigned *v)
        {
            Real *xj, *yj, *ck;

            for (long k = ncc - 1; k >= 0; k--)
            {
                for (int j = 0; j < ncx; j++)
                {
                    xj = &x[j * nrx];
                    yj = &x[(k * ncx + j)*nrx];
                    ck = &c[k * nrc];
                    
                    for (int i = 0; i < nrx; i++)
                    {
                        yj[i] = ck[v[i]-1] * xj[i];
                    }
                }
            }
        }
        
        /*
         * Create a model matrix for the input variables. For instance, in a two-way ANOVA:
         *
         *     (nvars) fert  seed fert:seed (nterms)
         *      corn     0     0      0
         *      fert     1     0      1
         *      seed     0     1      1
         */
        
        template <typename T> Matrix modelMatrix(const std::vector<const T*> &vars)
        {
            int i, k;
            
            // Number of variables
            const auto nVars = vars.size();

            const auto nf = std::count_if(vars.begin(), vars.end(), [&](const T *t)
            {
                return (*t).isFactor();
            });
            
            SS_ASSERT(nf, "Factor variable required for model matrix");
            
            /*
             * Build a pattern matrix for the model
             */
            
            Vector p;
            
            if (nf == 1)
            {
                if (nVars == 1)
                {
                    p = Vector(1);
                    p << 1;
                }
                else
                {
                    p = Vector(2);
                    p << 0, 1;
                }
            }
            else
            {
                if (nVars == 2)
                {
                    p = Vector(6);
                    p << 1, 0, 0, 1, 1, 1;
                }
                else
                {
                    p = Vector(9);
                    p << 0, 1, 0, 0, 0, 1, 0, 1, 1;
                }
            }
            
            // Number of terms
            const auto nTerms = nf == 1 ? 1 : 3;
            
            // Number of intercepts
            const auto intrcept = 1;
            
            // Number of samples
            const auto n = (*vars[0]).size();
            
            /*
             * This section of the code checks the types of the variables in the model frame. Note
             * that it should really only check the variables if they appear in a term in the model.
             * Because it does not, we need to allow other types here, as they might well occur
             * on the LHS. The code converts all character variables in the model frame to factors,
             * so the only types that ought to be here are logical, integer (including factor), numeric
             * and complex.
             */
            
            std::vector<T> variable(nVars);
            std::vector<Real> nlevs(nVars);
            std::vector<Real> ordered(nVars);
            std::vector<Real> columns(nVars);
            
            for (auto i = 0; i < nVars; i++)
            {
                variable[i] = (*vars[i]);
                
                if (variable[i].isNumeric())
                {
                    ordered[i] = 0;
                    nlevs[i]   = 0;
                    columns[i] = variable[i].size();
                }
                else if (variable[i].isFactor())
                {
                    ordered[i] = 0;
                    
                    if ((nlevs[i] = variable[i].levels().size()) < 1)
                    {
                        throw std::runtime_error("Variable " + std::to_string(i+1) + " has no levels");
                    }

                    // Will get updated later when contrasts are set
                    columns[i] = 1; //ncols(var_i);
                }
                
                ordered[i] = 0;
                columns[i] = 1;
            }
            
            /*
             * Compute the dummy variable matrices for the each of the categorical variable.
             */
            
            std::vector<Matrix> contr1(nVars);
            
            for (auto i = 0; i < nVars; i++)
            {
                // Don't code unless it's a categorical variable
                if (nlevs[i])
                {
                    Internal::contrTreatment(contr1[i], nlevs[i]);
                }
            }
            
            /*
             * We now have everything needed to build the design matrix. The first step is to
             * compute the matrix size and to allocate it. Note that "count" holds a count of
             * how many columns there are for each term in the model and "nc" gives the total
             * column count.
             */
            
            std::vector<Real> count(nTerms);
            
            // Total number of coefficients required
            auto nc = intrcept ? 1 : 0;
            
            /*
             * Calculating number of coefficents required. For example in an ANOVA 5x3, there'd be
             * 4 + 2 + (4 * 2) = 14. We do that by iterating vertically in the factor matrix.
             */
            
            for (auto j = 0; j < nTerms; j++)
            {
                // Number of coefficent required for the iteration
                Real dk = 1;
                
                for (auto i = 0; i < nVars; i++)
                {
                    // Is this coefficient also a factor?
                    if (p[i + j * nVars])
                    {
                        // Is this a categorical variable?
                        if (nlevs[i])
                        {
                            dk *= contr1[i].cols();
                        }
                    }
                }
                
                count[j] = (int) dk;
                nc += dk;
            }
            
            SS_ASSERT(nc <= INT_MAX, "Matrix requires too many columns");
            
            /*
             * Record which columns of the design matrix are associated with which model terms.
             *
             *     Eg: 0 1 1 1 1 2 2 3 3 3 3 3 3 3 3
             */
            
            auto assign = std::vector<unsigned>();
            assign.resize(nc);
            
            k = 0;
            
            if (intrcept)
            {
                assign[k++] = 0;
            }
            
            for (auto j = 0; j < nTerms; j++)
            {
                for (auto i = 0; i < count[j]; i++)
                {
                    assign[k++] = j+1;
                }
            }
            
            /*
             * Allocate and compute the design matrix
             */
            
            auto rx = new Real[n * nc];
            
            for (int i = 0; i < (n * nc); i++) { rx[i] = 0.0; }
            
            std::size_t jstart, jnext;
            
            /*
             * Begin with a column of 1s for the intercept
             */
            
            if ((jnext = jstart = intrcept) != 0)
            {
                for (i = 0; i < n; i++)
                {
                    rx[i] = 1.0;
                }
            }
            
            /*
             * Now loop over the model terms
             */
            
            for (auto k = 0; k < nTerms; k++)
            {
                if (k == -1)
                {
                    continue;
                }
                
                for (auto i = 0; i < nVars; i++)
                {
                    if (columns[i] == 0)
                    {
                        continue;
                    }
                    
                    auto var_i = variable[i]; //static_cast<Categorical *>(variable[i].get());
                    const auto f = p[i + k * nVars];
                    
                    if (f)
                    {
                        /*
                         * Convert the matrix into a vector. This is easily done in R with:
                         *
                         *      contrast = coerceVector(VECTOR_ELT(contr1, i), REALSXP);
                         */
                        
                        auto m = contr1[i];
                        double *contrast = new double[1234];
                        Eigen::Map<Eigen::MatrixXd>(contrast, m.rows(), m.cols()) = m;
                        
                        if (jnext == jstart)
                        {
                            if (nlevs[i] > 0)
                            {
                                // Add main effects
                                Internal::firstFactor(&rx[jstart * n],
                                                      n,
                                                      jnext - jstart,
                                                      contrast,
                                                      m.rows(),
                                                      m.cols(),
                                                      var_i.factors());
                                jnext += m.cols();
                            }
                        }
                        else
                        {
                            if (nlevs[i] > 0)
                            {
                                const auto nrc = var_i.levels().size();
                                const auto ncc = var_i.levels().size() - 1;
                                
                                // Add interaction effects
                                Internal::addFactor(&rx[jstart * n],
                                                    n,
                                                    jnext - jstart,
                                                    var_i.data__(),
                                                    nrc,
                                                    ncc,
                                                    var_i.factors());
                                
                                jnext += (jnext - jstart)*(var_i.levels().size() - 2);
                            }
                        }
                        
                        // Don't forget to set it back!
                        contr1[i] = Eigen::Map<Eigen::MatrixXd>(contrast, m.rows(), m.cols());
                        
                        delete []contrast;
                    }
                }
                
                jstart = jnext;
            }
            
            Matrix m(n, nc);
            
            for (auto j = 0; j < nc; j++)
            {
                for (auto i = 0; i < n; i++)
                {
                    m(i,j) = rx[j*n + i];
                }
            }
            
            return m;
        }
    }
    
    /*
     * Create a model matrix for the input variables
     */
    
    template <typename T, typename... Args> Matrix modelMatrix(const T &t, Args... args)
    {
        Matrix r;
        
        Internal::vArgs([&](const std::vector<const T *> &p)
        {
            r = Internal::modelMatrix(p);
        }, t, args...);
        
        return r;
    }
}

#endif