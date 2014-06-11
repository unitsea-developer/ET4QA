/* 
 * File:   Regression.hpp
 * Author: matthewsupernaw
 *
 * Created on June 10, 2014, 7:57 AM
 */

#ifndef REGRESSION_HPP
#define	REGRESSION_HPP
#include <vector>

#include "FunctionMinimizer.hpp"
#include "Statistics.hpp"


namespace ad {

    template<class REAL_T>
    class RegressionObject : public ad::FunctionMinimizer<REAL_T> {
    protected:
        std::vector<REAL_T> x_m;
        std::vector<REAL_T> y_m;
        std::vector<REAL_T> residuals_m;
        uint32_t degrees_of_freedom;
        REAL_T residualStandardError;
        REAL_T RSquared;
        REAL_T AdjustedRSquared;
        REAL_T PValue;
    public:

        RegressionObject(const std::vector<REAL_T> &x
                , const std::vector<REAL_T> &y)
        : x_m(x), y_m(y) {

        }

        RegressionObject() {

        }

        std::vector<REAL_T> GetX() const {
            return x_m;
        }

        void SetX(std::vector<REAL_T> x) {
            this->x_m = x;
        }

        std::vector<REAL_T> GetY() const {
            return y_m;
        }

        void SetY(std::vector<REAL_T> y) {
            this->y_m = y;
        }

        std::vector<REAL_T> GetResiduals() const {
            return residuals_m;
        }

        virtual const REAL_T Evaluate(const REAL_T &x) = 0;

        void Finalize() {
            this->ComputeFitStatistics();
        }

        virtual std::string ToString() {

        }

    protected:

        void ComputeFitStatistics() {
            degrees_of_freedom = y_m.size() - this->parameters_m.size();
            REAL_T ymean = ad::Mean<REAL_T > (y_m);
            REAL_T ssregression = 0;
            REAL_T sstotal = 0;
            REAL_T sse = 0;
            this->residuals_m = std::vector<REAL_T > (x_m.size());
            for (int i = 0; i< this->residuals_m.size(); i++) {
                REAL_T eval = this->Evaluate(x_m[i]);
                this->residuals_m[i] = (eval - y_m[i])*(eval - y_m[i]);
                sse += this->residuals_m[i];
                ssregression += (eval - ymean)*(eval - ymean);
                sstotal += (y_m[i] - ymean)*(y_m[i] - ymean);
            }

            this->residualStandardError = std::sqrt(sse / (static_cast<REAL_T> (this->degrees_of_freedom)));
            this->RSquared = 1.0 - sse / sstotal;
            this->AdjustedRSquared = 1.0 - (sse /
                    (static_cast<REAL_T> (y_m.size()) -
                    static_cast<REAL_T> (this->parameters_m.size())))
                    / (sstotal / (static_cast<REAL_T> (y_m.size()) - static_cast<REAL_T> (1)));

        }



    };

    /**
     *Regression object for 
     * y = mx+b
     */
    template<class REAL_T>
    class LinearRegression : public RegressionObject<REAL_T> {
        ad::Variable<REAL_T> mp;
        ad::Variable<REAL_T> bp;

        REAL_T m;
        REAL_T b;
    public:

        LinearRegression(const std::vector<REAL_T> &x
                , const std::vector<REAL_T> &y) : RegressionObject<REAL_T>(x, y) {
        }

        void Initialize() {
            //            this->predicted_m = std::vector<Variable<REAL_T> >(this->x_m.size());
            REAL_T sumX = REAL_T(0); //this->sum_m.X();
            REAL_T sumY = REAL_T(0); //this->sum_m.Y();
            REAL_T sumXX = REAL_T(0); //this->sum_m.X() * this->sum_m.X();
            REAL_T sumXY = REAL_T(0); //this->sum_m.X() * this->sum_m.Y();
            for (int i = 0; i <this->x_m.size(); i++) {
                REAL_T x = this->x_m[i];
                REAL_T y = this->y_m[i];
                sumX += x;
                sumY += y;
                REAL_T xx = x * x;
                sumXX += xx;
                REAL_T xy = x * y;
                sumXY += xy;
            }
            REAL_T sxx = sumXX - (sumX * sumX) / static_cast<REAL_T> (this->x_m.size());
            REAL_T sxy = sumXY - (sumX * sumY) / static_cast<REAL_T> (this->x_m.size());
            REAL_T xbar = sumX / static_cast<REAL_T> (this->x_m.size());
            REAL_T ybar = sumY / static_cast<REAL_T> (this->x_m.size());


            mp.SetName("m");
            bp.SetName("b");
            bp = sxy / sxx;
            mp = ybar - bp * xbar;
            this->Register(mp);
            this->Register(bp);
            this->SetVerbose(false);
        }

        void ObjectiveFunction(ad::Variable<REAL_T> &f) {
            f = 0.0;
            ad::Variable<REAL_T> temp;
            ad::Variable<REAL_T> sum = static_cast<REAL_T> (0);
            for (int i = 0; i< this->y_m.size(); i++) {
                temp = this->mp * this->x_m[i] + this->bp;
                f += (temp - this->y_m[i])*(temp - this->y_m[i]);

            }
        }

        void Finalize() {
            this->m = mp.GetValue();
            this->b = bp.GetValue();
            this->ComputeFitStatistics();


        }

        const REAL_T Evaluate(const REAL_T &x) {
            return m * x + b;
        }

        std::string ToString() {
            std::stringstream ss;

            ss << "Model: y = mx+b\n";
            ss << "Paramters:\n";
            ss << "m: " << this->m << "\n";
            ss << "b: " << this->b << "\n";
            ss << "Fit Statistics:\n";
            ss << "Residual Std. Error: " << this->residualStandardError << "\n";
            ss << "Degrees of Freedom: " << this->degrees_of_freedom << "\n";
            ss << "R-squared: " << this->RSquared << "\n";
            ss << "Adjusted R-squared: " << this->AdjustedRSquared;
            return ss.str();

        }



    };

    template<class REAL_T>
    class PolynomialRegression : public RegressionObject<REAL_T> {
        uint32_t order_m;
        std::vector<ad::Variable<REAL_T> > coefficients_m;
    public:

        PolynomialRegression(uint32_t order, const std::vector<REAL_T> &x
                , const std::vector<REAL_T> &y) : order_m(order), RegressionObject<REAL_T>(x, y) {

        }

        void Initialize() {
            this->SetVerbose(false);
            coefficients_m = std::vector<ad::Variable<REAL_T> >(this->order_m);


            size_t itemCount = this->y_m.size();


            std::vector<REAL_T> data(2 * itemCount);
            for (int item = 0; item < itemCount; item++) {
                data[0 + item] = this->x_m[item];
                data[1 * itemCount + item] = this->y_m[item];

            }

            size_t equations = order_m + 1;
            size_t coefficients = order_m + 2;

            std::vector<REAL_T> result(equations + 1);
            std::vector<REAL_T> matrix(equations * coefficients);
            REAL_T sumX = 0.0;
            REAL_T sumY = 0.0;

            for (size_t item = 0; item < itemCount; item++) {
                sumX += this->x_m[item];
                sumY += this->y_m[item];
                for (size_t eq = 0; eq < equations; eq++) {
                    for (size_t coe = 0; coe < coefficients - 1; coe++) {
                        matrix[eq * coefficients + coe] += std::pow(this->x_m[item], REAL_T(eq + coe));
                    }

                    matrix[eq * coefficients + coefficients - 1] += this->y_m[item]
                            * std::pow(this->x_m[item], REAL_T(eq));
                }
            }



            std::vector<REAL_T> subMatrix;
            size_t r, c;
            Sub(matrix, equations, coefficients, subMatrix, r, c);

            for (size_t eq = 1; eq < equations; eq++) {
                matrix[eq * coefficients] = REAL_T(0);
                for (size_t coe = 1; coe < coefficients; coe++) {
                    matrix[eq * coefficients + coe] = subMatrix[(eq - 1) * c + (coe - 1)];
                }
            }


            for (int eq = (equations - 1); eq > -1; eq--) {
                REAL_T value = matrix[eq * coefficients + (coefficients - 1)];
                for (size_t coe = eq; coe < coefficients - 1; coe++) {
                    value -= matrix[eq * coefficients + coe] * coefficients_m[coe].GetValue();
                }
                coefficients_m[eq] = value / matrix[eq * coefficients + eq];
            }

            for (int i = 0; i < order_m; i++) {
                this->Register(coefficients_m[i]);
            }


        }

        void ObjectiveFunction(ad::Variable<REAL_T> &f) {
            f = 0.0;
            ad::Variable<REAL_T> temp;
            ad::Variable<REAL_T> sum = static_cast<REAL_T> (0);
            for (int i = 0; i< this->y_m.size(); i++) {
                temp = 0;
                for (int j = 0; j < order_m; j++) {
                    temp += coefficients_m[j] * std::pow(this->x_m[i], j);
                }
                f += (temp - this->y_m[i])*(temp - this->y_m[i]);
            }
            //            f =static_cast<REAL_T>(this->y_m.size())*std::log(f);
        }

        const REAL_T Evaluate(const REAL_T &x) {
            REAL_T temp = 0;
            for (int j = 0; j < order_m; j++) {
                temp += coefficients_m[j].GetValue() * std::pow(x, j);
            }
            return temp;
        }

        std::string ToString() {
            std::stringstream ss;

            ss << "Model: y = a_n*x^n a_n-1*x^n-1 +...a_1*x+a0  \n";
            ss << "Paramters:\n";
            for (int i = 0; i < order_m; i++) {
                ss << coefficients_m[i] << " ";
            }
            ss << "\n";
            ss << "Fit Statistics:\n";
            ss << "Residual Std. Error: " << this->residualStandardError << "\n";
            ss << "Degrees of Freedom: " << this->degrees_of_freedom << "\n";
            ss << "R-squared: " << this->RSquared << "\n";
            ss << "Adjusted R-squared: " << this->AdjustedRSquared;
            return ss.str();

        }



    private:

        void Sub(const std::vector<REAL_T> &m, size_t rows, size_t cols, std::vector<REAL_T> &result, size_t &nrows, size_t &ncols) {
            size_t equations = rows;
            size_t coefficients = cols;
            nrows = (equations - 1);
            ncols = (coefficients - 1);
            result = std::vector<REAL_T > ((equations - 1)* (coefficients - 1));

            for (size_t eq = 1; eq < equations; eq++) {
                REAL_T div = m[eq * cols];

                if (div == REAL_T(0)) {
                    div = std::numeric_limits<REAL_T>::epsilon();
                }
                REAL_T factor = m[0] / div;
                for (size_t coe = 1; coe < coefficients; coe++) {
                    result[(eq - 1) * ncols + (coe - 1)] = m[ coe] - m[eq * cols + coe] * factor;
                }
            }


            if (equations == 1) {
                return;
            }

            if (result[0] == 0) {
                bool found = false;
                for (int i = 0; i < result.size(); i++) {
                    if (result[i * cols] != 0) {
                        found = true;
                        std::vector<REAL_T> temp(cols); // = Row(m, 0);

                        for (int c = 0; c < cols; c++) {
                            temp[c] = m[c];
                        }

                        for (size_t j = 0; j < (coefficients - 1); j++) {
                            result[j] = result[i * (coefficients - 1) + j];
                        }
                        for (size_t j = 0; j < (equations - 1); j++) {
                            result[i * (coefficients - 1) + j] = temp[j];
                        }
                        break;
                    }
                }
                if (!found) {
                    std::cout << "Equation has no solution!";
                    return;
                }
            }

            size_t r, c;
            std::vector<REAL_T> subMatrix;
            Sub(result, nrows, ncols, subMatrix, r, c);
            for (size_t eq = 1; eq < equations - 1; eq++) {
                result[eq * (coefficients - 1)] = 0;
                for (size_t coe = 1; coe < coefficients - 1; coe++) {
                    result[eq * ncols + coe] = subMatrix[(eq - 1) * c + (coe - 1)];
                }
            }



        }


    };

    template<class REAL_T>
    class LogrithmicRegression : public RegressionObject<REAL_T> {
        ad::Variable<REAL_T> a;
        ad::Variable<REAL_T> b;
    public:

        LogrithmicRegression(const std::vector<REAL_T> &x
                , const std::vector<REAL_T> &y) : RegressionObject<REAL_T>(x, y) {

        }
                
                
        void Initialize() {
            a =1.0;
            b =1.0;
            this->Register(a);
            this->Register(b);
        }

        void ObjectiveFunction(ad::Variable<REAL_T> &f) {
            f = 0.0;
            ad::Variable<REAL_T> temp;
            
            for (int i = 0; i< this->y_m.size(); i++) {
                temp = this->a +this->b*std::log(this->x_m[i]);
                f += (temp - this->y_m[i])*(temp - this->y_m[i]);
            }
            //            f =static_cast<REAL_T>(this->y_m.size())*std::log(f);
        }

        const REAL_T Evaluate(const REAL_T &x) {

        }


    };

    template<class REAL_T>
    class ExponetialRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class ExponetialRecoveryRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class PowerRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class GammaRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class GaussianRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class GaussianOffsetRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class RodbardRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    class InverseRegression : public RegressionObject<REAL_T> {
    };

    template<class REAL_T>
    void Exponential(const std::vector<REAL_T> &x, const std::vector<REAL_T> &y) {

    }





}




#endif	/* REGRESSION_HPP */

