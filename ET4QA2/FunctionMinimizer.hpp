/* 
 * File:   FunctionMinimizer.hpp
 * Author: matthewsupernaw
 *
 * Created on April 17, 2014, 1:10 PM
 */

#ifndef FUNCTIONMINIMIZER_HPP
#define	FUNCTIONMINIMIZER_HPP


//#include "support/admb/build/dist/include/fvar.hpp"
#include <valarray>
#include <vector>
#include <iomanip>
#include <sys/timeb.h>
#include <sstream>

#if defined(WIN32) || defined(WIN64)

#ifndef BOLD
#define BOLD ""
#endif

#ifndef DEFAULT_IO
#define DEFAULT_IO ""
#endif

#ifndef BLACK
#define BLACK ""
#endif

#ifndef RED
#define RED ""
#endif

#ifndef GREEN
#define GREEN ""
#endif

#ifdef BROWN
#define BROWN ""
#endif

#ifndef BLUE
#define BLUE ""
#endif

#ifndef MAGENTA
#define MAGENTA ""
#endif

#ifndef CYAN
#define CYAN ""
#endif

#ifndef GRAY
#define GRAY ""
#endif

#else

#ifndef BOLD
#define BOLD "\033[1m"
#endif

#ifndef DEFAULT_IO
#define DEFAULT_IO "\033[0m"
#endif

#ifndef BLACK
#define BLACK "\033[0;30m"
#endif

#ifndef RED
#define RED "\033[0;31m"
#endif

#ifndef GREEN
#define GREEN "\033[0;32m"
#endif

#ifdef BROWN
#define BROWN "\033[0;33m"
#endif

#ifndef BLUE
#define BLUE "\033[0;34m"
#endif

#ifndef MAGENTA
#define MAGENTA "\033[0;35m"
#endif

#ifndef CYAN
#define CYAN "\033[0;36m"
#endif

#ifndef GRAY
#define GRAY "\033[0;37m"
#endif

#endif

#define HAVE_ADMB

#ifdef HAVE_ADMB
#include <fvar.hpp>
#endif
namespace ad {

    //#define HAVE_GSL


#ifdef HAVE_GSL
#include <gsl/gsl_multimin.h>
    extern "C"
    double function_value_callback(const gsl_vector* x, void* params);

    void function_gradient_callback(const gsl_vector* x, void* params, gsl_vector* gradJ);

    void function_value_and_gradient_callback(const gsl_vector* x, void* params,
            double* J, gsl_vector* gradJ);
#endif

    /**
     *A derivative based function minimizer.
     */
    template<class T>
    class FunctionMinimizer {
    public:

        enum MinimizerType {
            DUBOUT_LBFGS = 0,
            NEWTON,
#ifdef HAVE_ADMB
            ADMB_AUTODIFF_MINIMIZER,
#endif
                    
#ifdef HAVE_GSL
            GSL_CONJUGATE_FR,
            GSL_CONJUGATE_PR,
            GSL_BFGS,
            GSL_BFGS2,
            GSL_STEEPEST_DESCENT
#endif
        };
    protected:



#ifdef HAVE_GSL
        friend double function_value_callback(const gsl_vector* x, void* params);

        friend void function_gradient_callback(const gsl_vector* x, void* params, gsl_vector* grad);

        friend void function_value_and_gradient_callback(const gsl_vector* x, void* params,
                double* fx, gsl_vector* grad);

#endif

        MinimizerType minimizer_type_m;
        std::vector<ad::Variable<T>* > active_parameters_m;
        std::vector<ad::Variable<T>* > parameters_m;
        std::vector<unsigned int> phases_m;
        std::vector<T> lower_bounds_m;
        std::vector<T> upper_bounds_m;

        T tolerance_m;
        unsigned int max_iterations_m;
        unsigned int iteration_m;
        unsigned int phase_m;
        unsigned int max_phase_m;
        std::vector<bool> is_constrained_m;

        size_t sum_time_in_user_function_m;
        size_t average_time_in_user_function_m;
        size_t sum_time_in_grad_calc_m;
        size_t average_time_in_grad_calc_m;
        size_t function_calls_m;
        size_t gradient_calls_m;

        bool has_constraints_m;
        T function_value_m;
        ad::Variable<T> function_result_m;
        std::valarray<T> gradient_m;

        bool verbose_m;
        unsigned int iprint_m;
        size_t max_history_m;

        T max_c;
        int unrecorded_calls_m;

    public:

        /**
         * Default constructor.
         */
        FunctionMinimizer()
        : tolerance_m(T(1e-4)),
        max_iterations_m(1000),
        is_constrained_m(false),
        verbose_m(true),
        iprint_m(10),
        max_history_m(1000),
        unrecorded_calls_m(0),
        minimizer_type_m(DUBOUT_LBFGS),
        max_c(std::numeric_limits<T>::min()) {

        }

        virtual ~FunctionMinimizer() {
            //            std::ofstream hess;
            //            std::ofstream grad;
            //            std::setprecision(50);
            //            grad.open("gradient.txt");
            //            hess.open("hessian.txt");
            //            for (int i = 0; i < this->active_parameters_m.size(); i++) {
            //                grad<<this->gradient_m[i]<<"\n";
            //                for (int j = 0; j < this->active_parameters_m.size(); j++) {
            //                    hess<<this->function_result_m.WRT(*this->active_parameters_m[i],*this->active_parameters_m[j]).GetValue()<<",";
            //                }
            //                hess<<"\n";
            //            }
            //            hess.close();
        }

        /**
         * Returns the iteration interval in which runtime information is sent
         * to stdout when verbose is set to true. 
         * 
         * @return 
         */
        unsigned int GetPrintInterval() const {
            return iprint_m;
        }

        /**
         * Sets the iteration interval in which runtime information is sent
         * to stdout when verbose is set to true. 
         * 
         * @param iprint
         */
        void SetPrintInterval(unsigned int iprint) {
            this->iprint_m = iprint;
        }

        /**
         * Returns the maximum history used by the l-bfgs algorithm.
         * @return 
         */
        size_t GetMaxHistory() const {
            return max_history_m;
        }

        /**
         * Sets the maximum history used by the l-bfgs algorithm.
         * 
         * @param max_history
         */
        void SetMaxHistory(size_t max_history) {
            this->max_history_m = max_history;
        }

        /**
         *
         * Returns the maximum number of iterations used in each call to 
         * the l-bfgs algorithm.
         * 
         * @return 
         */
        unsigned int GetMaxIterations() const {
            return max_iterations_m;
        }

        /**
         * Sets the maximum number of iterations used in each call to 
         * the l-bfgs algorithm.
         * 
         * @param max_iterations
         */
        void SetMaxIterations(unsigned int max_iterations) {
            this->max_iterations_m = max_iterations;
        }

        /**
         * Returns the tolerance for this minimizer.
         * @return 
         */
        T GetTolerance() const {
            return tolerance_m;
        }

        /**
         * Set the derivative tolerance for this minimizer. 
         * Default is 1e-4.
         * @param tolerance
         */
        void SetTolerance(T tolerance) {
            this->tolerance_m = tolerance;
        }

        /**
         * Is verbose on or off;
         * 
         * @return 
         */
        bool IsVerbose() const {
            return verbose_m;
        }

        /**
         * Set whether or not runtime information should be sent to stdout.
         * @param verbose
         */
        void SetVerbose(bool verbose) {
            this->verbose_m = verbose;
        }

        /**
         * Current phase.
         * 
         * @return 
         */
        unsigned int Phase() {
            return this->phase_m;
        }

        /**
         * Registers a unconstrained Variable<T> to the minimizer.
         * 
         * @param var
         * @param phase  -phase in which this Variable<T> becomes active.
         */
        void Register(ad::Variable<T> &var, unsigned int phase = 1) {
            this->parameters_m.push_back(&var);
            this->phases_m.push_back(phase);
            this->is_constrained_m.push_back(false);

        }

        /**
         * Registers a bounded Variable<T> to the minimizer. 
         * @param var 
         * @param lower_-bound -lower boundary
         * @param upper -bound -upper boundary
         * @param phase -phase in which this Variable<T> becomes active.
         */
        void Register(ad::Variable<T> &var, T lower_bound, T upper_bound, unsigned int phase = 1) {
            this->parameters_m.push_back(&var);
            this->phases_m.push_back(phase);
            this->lower_bounds_m.push_back(lower_bound);
            this->upper_bounds_m.push_back(upper_bound);
            this->is_constrained_m.push_back(true);

        }

        /**
         * Returns the last computed gradient vector.
         * @return 
         */
        std::vector<T> GetGradient() const {
            return gradient_m;
        }

        /**
         * Returns the last computed function value.
         * @return 
         */
        ad::Variable<T> GetFunctionValue() const {
            return function_result_m;
        }

        /**
         * Abstract function. Called before minimization begins.
         */
        virtual void Initialize() {

        }

        /**
         * Abstract function. The function to be minimized.
         * @param f -the value that is minimized.
         */
        virtual void ObjectiveFunction(ad::Variable<T> &f) {

        }

        /**
         * Abstract function. Called after the first phase and between phases 
         * in the minimization.
         */
        virtual void TransitionPhase() {

        }

        /**
         * Starts the minimizer. Initializes runtime members and handles 
         * phasing. Also, makes sure that bounded parameters are properly 
         * initialized. If a bounded parameter has an initial value outside of 
         * the specified bounds, the initial value is set to the center of the
         * bounded values.  
         * 
         * 
         * @return 
         */
        bool Run(MinimizerType type = DUBOUT_LBFGS) {
            this->minimizer_type_m = type;
            this->Initialize();
            this->max_phase_m = 1;
            this->max_c = 0.0;
            this->function_calls_m = 0;
            this->gradient_calls_m = 0;
            this->sum_time_in_user_function_m = 0;
            this->average_time_in_user_function_m = 0;
            this->sum_time_in_grad_calc_m = 0;
            this->average_time_in_grad_calc_m = 0;
            this->has_constraints_m = false;

            bool ret = false;

            for (int i = 0; i < parameters_m.size(); i++) {
                if (this->is_constrained_m[i]) {
                    this->has_constraints_m = true;
                    if (lower_bounds_m[i] > upper_bounds_m[i]) {
                        T temp = lower_bounds_m[i];
                        lower_bounds_m[i] = upper_bounds_m[i];
                        upper_bounds_m[i] = temp;
                    }

                    if (parameters_m[i]->GetValue() < lower_bounds_m[i]) {
                        parameters_m[i]->SetValue((lower_bounds_m[i] +(upper_bounds_m[i] - lower_bounds_m[i]) / 2.0));
                    } else if (parameters_m[i]->GetValue() > upper_bounds_m[i]) {
                        parameters_m[i]->SetValue((lower_bounds_m[i] +(upper_bounds_m[i] - lower_bounds_m[i]) / 2.0));
                    }
                }
            }

            //size_t max_phase = 1;
            for (int i = 0; i < this->phases_m.size(); i++) {
                if (this->phases_m[i] > max_phase_m) {
                    max_phase_m = phases_m[i];
                }
            }

            for (int p = 0; p < max_phase_m; p++) {
                this->phase_m = (p + 1);
                this->active_parameters_m.erase(active_parameters_m.begin(), active_parameters_m.end());

                for (int i = 0; i < this->parameters_m.size(); i++) {
                    if (this->phases_m[i] <= (p + 1)) {
                        this->parameters_m[i]->SetAsIndependent(true);
                        this->active_parameters_m.push_back(this->parameters_m[i]);
                    }
                }
                this->gradient_m.resize(this->active_parameters_m.size(), 0.0);
                //                std::cout << this->gradient_m.size() << "<<---" << std::flush;

                switch (this->minimizer_type_m) {
                    case DUBOUT_LBFGS:
                        ret = this->LBFGS(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
                    case NEWTON:
                        //                        ret = this->Newton(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
#ifdef HAVE_ADMB
                    case ADMB_AUTODIFF_MINIMIZER:
                        ret = this->ADMB_Minimizer(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
#endif
#ifdef HAVE_GSL
                    case GSL_CONJUGATE_FR:
                        ret = this->GSL_Multimin(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
                    case GSL_CONJUGATE_PR:
                        ret = this->GSL_Multimin(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
                    case GSL_BFGS:
                        ret = this->GSL_Multimin(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
                    case GSL_BFGS2:
                        ret = this->GSL_Multimin(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
                    case GSL_STEEPEST_DESCENT:
                        ret = this->GSL_Multimin(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
#endif
                    default:
                        ret = this->LBFGS(this->active_parameters_m, this->GetMaxIterations(), this->GetTolerance());
                        break;
                }

                //                this->Print(this->function_result_m, this->gradient_m, active_parameters_m, "Verbose:\nTransition");

                this->TransitionPhase();
            }


            //            this->LBFGS(this->parameters_m, this->GetMaxIterations(), this->GetTolerance());
            if (this->verbose_m) {
                this->ObjectiveFunction(function_result_m);
                this->Print(this->function_result_m, this->gradient_m, active_parameters_m, "Verbose:\nFinal Statistics");
            }
            this->Finalize();

            return ret;
        }

        /**
         * Abstract function. Called after minimizer is complete.
         */
        virtual void Finalize() {

        }

    protected:

        std::vector<ad::Variable<T>*> GetActiveParameters() const {
            return active_parameters_m;
        }

        /**
         * Computes the gradient with respect to active parameters. Also tracks 
         * the number of gradient function calls and the average time spent computing 
         * gradients. 
         * @param fx
         * @param parameters
         * @param gradient
         */
        virtual void Gradient(const ad::Variable<T> &fx, const std::vector<ad::Variable<T>* > &parameters, std::valarray<T> &gradient) {

            //            fx.Gradient(parameters,gradient);
            for (int i = 0; i < parameters.size(); i++) {
                // std::cout<<"Gradient i = "<<i<<std::endl;
                gradient[i] = fx.WRT(*parameters[i]);

                this->gradient_m[i] = gradient[i];
                // this->gradient_m[i] = gradient[i];
                if (std::fabs(gradient[i]) > max_c) {
                    max_c = std::fabs(gradient[i]);
                }
            }

        }

        const ad::Variable<T> GetCurrentFunctionValue() {
            ad::Variable<T> ret;
            this->ObjectiveFunction(ret);
            return ret;
        }

        const std::valarray<T> CalculateGradient() {
            std::valarray<T> gradient(this->active_parameters_m.size());
            ad::Variable<T> f;
            this->ObjectiveFunction(f);

            for (int i = 0; i < this->active_parameters_m.size(); i++) {
                gradient[i] = f.WRT(*active_parameters_m[i]);

            }

            return gradient;
        }

        const std::valarray<std::valarray<T> > EstimatedHessian() {
            Variable<T>::SetRecording(true);
            std::valarray<std::valarray<T> > hessian(
                    std::valarray<T > (this->active_parameters_m.size()), this->active_parameters_m.size());
            ad::Variable<T> f;
            for (size_t i = 0; i < active_parameters_m.size(); i++) {
                //   gradient(i) = f.NthPartial(parameters(i), 1);
                if (this->IsVerbose())
                    std::cout << "Estimating hessian row " << i << "\n";
                for (size_t j = 0; j < active_parameters_m.size(); j++) {
                    hessian[i] [j] = EstimateHessianElement(f, i, j);
                }
            }
            return hessian;
        }

    private:

        /**
         * Finite differences for second derivative.
         */
        const T EstimateHessianElement(ad::Variable<T> &f, size_t row, size_t column) {
            T h = T(0.0001);
            T hh = T(1.0) / h;
            size_t n = active_parameters_m.size();
            T ff;



            std::vector<T> xp2h(n);
            std::vector<T> xph(n);
            std::vector<T> xm2h(n);
            std::vector<T> xmh(n);

            ad::Variable<T> fp2h;
            ad::Variable<T> fph;
            ad::Variable<T> fm2h;
            ad::Variable<T> fmh;

            //updaste parameters
            for (size_t j = 0; j < n; j++) {
                if (j != column) {


                    xp2h[j] = active_parameters_m[j]->GetValue();
                    xph[j] = active_parameters_m[j]->GetValue();
                    xm2h[j] = active_parameters_m[j]->GetValue();
                    xmh[j] = active_parameters_m[j]->GetValue();
                } else {

                    xp2h[j] = active_parameters_m[j]->GetValue() + T(2.0) * h;
                    xph[j] = active_parameters_m[j]->GetValue() + h;
                    xm2h[j] = active_parameters_m[j]->GetValue() - T(2.0) * h;
                    xmh[j] = active_parameters_m[j]->GetValue() - h;
                }
            }

            for (size_t i = 0; i < n; i++) {
                active_parameters_m[i]->SetValue(xp2h[i]);
            }
            this->ObjectiveFunction(fp2h); // dfdx(f, xp2h, xiref);

            for (size_t i = 0; i < n; i++) {
                active_parameters_m[i]->SetValue(xph[i]);
            }
            this->ObjectiveFunction(fph); //dfdx(f, xph, xiref);

            for (size_t i = 0; i < n; i++) {
                active_parameters_m[i]->SetValue(xm2h[i]);
            }
            this->ObjectiveFunction(fm2h); //dfdx(f, xm2h, xiref);

            for (size_t i = 0; i < n; i++) {
                active_parameters_m[i]->SetValue(xmh[i]);
            }
            this->ObjectiveFunction(fmh); //dfdx(f, xmh, xiref);
            //            std::cout<<"here...\n";

            ff = (-1 * fp2h.WRT(*active_parameters_m[row]) +
                    8.0 * fph.WRT(*active_parameters_m[row]) -
                    8.0 * fmh.WRT(*active_parameters_m[row]) +
                    fm2h.WRT(*active_parameters_m[row])) / 12.0 * hh;


            return ff;

        }

        void CallGradient(ad::Variable<T> &fx, std::vector<ad::Variable<T>* > &parameters, std::valarray<T> &gradient) {
            this->gradient_calls_m++;
            this->max_c = 0;
            size_t start = this->GetMilliCount();
            Gradient(fx, parameters, gradient);
            size_t end = this->GetMilliCount();
            this->sum_time_in_grad_calc_m += (end - start);
            this->average_time_in_grad_calc_m = sum_time_in_grad_calc_m / this->gradient_calls_m;
        }

        /**
         * Intermediate function to call the ObjectiveFunction.
         * Tracks the amount of objective function calls and the average time 
         * spent in the objective function.
         * 
         * @param f
         */
        void CallObjectiveFunction(ad::Variable<T> &f) {
            //std::cout<<"called "<<__func__<<":"<<__LINE__<<std::endl;
            this->function_calls_m++;
            if (ad::Variable<T>::IsRecording()) {
                this->unrecorded_calls_m++;
            }
            clock_t start = GetMilliCount();
            this->ObjectiveFunction(f);
            clock_t end = GetMilliCount();
            //            this->function_result_m = ad::ADNumber<T > (f);

            sum_time_in_user_function_m += (end - start);
            average_time_in_user_function_m = sum_time_in_user_function_m / function_calls_m;
        }

        /**
         * \ingroup Matrix
         * Returns the determinant of Matrix m.
         * @param m
         * @return 
         */

        T Det(const std::valarray<std::valarray<T> > &m) {//this version is less expensive for ad types.


            T d = T(1);
            T ratio;


            size_t rows = m.size();
            size_t cols = m[0].size();

            if (rows == 1 && cols == 1) {
                return m[0][0];
            }

            if (rows == 2 && cols == 2) {
                d = (m[0][0] * m[1][1] -
                        m[0][1] * m[1][0]);
                return d;
            }


            //BELOW Here IS EXTREMELY SLOW!!
            size_t n = m.size();

            std::valarray<std::valarray<T> > matrix(std::valarray<T > (m[0].size(), 0), m.size());
            size_t i, j, k;

            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    if (j > i) {
                        ratio = m[j][i] / m[i][i];
                        for (k = 0; k < n; k++) {
                            T temp = matrix[j][k] - ratio * matrix[i][k];
                            //matrix[j][k] -= ratio * matrix[i][k];


                            //  matrix.SetValue(j, k, temp);
                            matrix[i][k] = temp;
                        }
                    }
                }
            }

            for (i = 0; i < n; i++) {
                d *= matrix[i][i]; //[i][i];
            }
            return d;

        }

        int Pivot(std::valarray<std::valarray<T> > &m, int row) {
            unsigned int k = row;
            T max = T(-1);
            T temp;

            for (size_t i = row; i < m.size(); i++) {
                if ((temp = std::fabs(m[i][ row])) > max && temp != 0.0) {
                    max = temp;
                    k = i;
                }

                if (m[k][row] == T(0)) {
                    return -1;
                }

                if (k != row) {
                    T tempVal;
                    for (size_t j = 0; j < m[0].size(); j++) {

                        tempVal = m[k][j];
                        m[k][j] = m[row][j];
                        m[row][ j] = tempVal;
                    }

                    return k;
                }
            }
            return 0;
        }

        void Inverse(std::valarray<std::valarray<T> >& m, std::valarray<std::valarray<T> >& ret) {
            //            std::valarray<std::valarray<T> > ret(std::valarray<T > (m[0].size(), 0),m.size());
            //ret =std::valarray<std::valarray<T> >(std::valarray<T > (m[0].size(), 1.0),m.size());
            std::valarray<std::valarray<T> > temp(std::valarray<T > (m[0].size(), 0), m.size());
            T d = Det(m);

            for (int i = 0; i < m.size(); i++)
                for (int j = 0; j < m.size(); j++)
                    ret[i][i] = 0.0;

            if (m.size() == 2 && m[0].size() == 2) {

                ret[0][0] = m[1][1] / d;
                ret[0][1] = m[0][1] / d * T(-1);
                ret[1][0] = m[1][0] / d * T(-1);
                ret[1][1] = m[0][0] / d;
                //                std::cout << ret[1][0] << " " << ret[1][1]<<std::endl;
                //                exit(0);
                //                return ret;
                return;
            }


            if (m.size() == m[0].size() && d != T(0)) {

                //                ret = std::valarray<std::valarray<T> >(std::valarray<T > (m[0].size(), 1.0), m.size());

                size_t i, j, k;
                T a1, a2, tempVal;


                for (k = 0; k < m.size(); k++) {

                    int index = Pivot(temp, k);
                    if (index == -1) {
                        //                        throw ArrayException("MATRIX Exception: call to Inverse() is not allowed, matrix is singular.");
                    }


                    if (index != 0) {
                        for (j = 0; j < m[0].size(); j++) {
                            tempVal = ret[k][j];
                            ret[k][j] = ret[index][j];
                            ret[index][j] = tempVal;
                        }
                    }

                    a1 = temp[k][k];

                    for (j = 0; j < m.size(); j++) {
                        tempVal = temp[k][j] / a1;
                        temp[k][j] = tempVal;
                        tempVal = ret[k][j] / a1;
                        ret[k][j] = tempVal;
                    }

                    for (i = 0; i < m.size(); i++) {
                        if (i != k) {
                            a2 = temp[i][k];
                            for (j = 0; j < m[0].size(); j++) {

                                tempVal = temp[i][j] - (a2 * temp[k][j]);
                                temp[i][j] = tempVal;
                                tempVal = ret[i][j] - (a2 * ret[k][j]);
                                ret[i][j] = tempVal;
                            }
                        }
                    }
                }

                //                return ret;
            } else {

                for (int i = 0; i < m.size(); i++)
                    ret[i][i] = 1.0;

                //                return ret;

            }
        }

        //        /**
        //         * A modified version of Charles Dubout's LBFGS algorithm released under the 
        //         * terms of GPL. 
        //         * http://www.idiap.ch/~cdubout/code/lbfgs.cpp
        //         * 
        //         * 
        //         * @param parameters
        //         * @param results
        //         * @param iterations
        //         * @param tolerance
        //         * @return 
        //         */
        //        bool LBFGS(std::vector<ad::Variable<T>* > &parameters, size_t iterations = 10000, T tolerance = (T(1e-5))) {
        //            //            std::cout<<"starting minimizer....\n";
        //            size_t max_history = max_history_m;
        //
        //            int maxLineSearches_ = 1000;
        //            std::valarray<T> x(parameters.size());
        //            //current gradient
        //            std::valarray<T> g(parameters.size());
        //
        //
        //            std::valarray<T> ng(x.size());
        //
        //            //initial evaluation
        //            ad::Variable<T> fx(0.0);
        //            T fx_value = 0.0;
        //
        //
        //            //Call the objective function and collect stats..
        //            this->CallObjectiveFunction(fx);
        //            this->function_value_m = fx.GetValue();
        //            //            ad::Variable<T> nfx(fx);
        //            //Historical evaluations
        //            std::valarray<T> px(parameters.size());
        //            std::valarray<T> pg(parameters.size());
        //            std::valarray<std::valarray<T> > dxs(std::valarray<T > (max_history), parameters.size());
        //            std::valarray<std::valarray<T> > dgs(std::valarray<T > (max_history), parameters.size());
        //
        //            //set parameters
        //            for (size_t i = 0; i < g.size(); i++) {
        //                x[i] = parameters[i]->GetValue();
        //            }
        //
        //            std::valarray<T> z(parameters.size());
        //
        //
        //            this->CallGradient(fx, parameters, g);
        //
        //
        //            T step;
        //            T relative_tolerance;
        //            T norm_g;
        //            const size_t nop = this->active_parameters_m.size();
        //
        //            for (int i = 0; i < iterations; ++i) {
        //
        //                iteration_m = i + 1;
        //
        //                norm_g = this->Norm(g);
        //
        //                // Backtracking using Wolfe's first condition (Armijo condition)
        //                step = i ? 1.0 : (1.0 / norm_g); //was just norm2
        //
        //                relative_tolerance = tolerance * std::max<T > (T(1.0), norm_g);
        //
        //                if (this->verbose_m && ((i % this->iprint_m) == 0)) {
        //                    this->Print(fx, g, parameters, "Verbose:\nMethod: L-BFGS");
        //                }
        //
        //                if (norm_g < relative_tolerance) {
        //
        //                    if (this->verbose_m) {
        //                        this->Print(fx, g, parameters, "Successful Convergence!");
        //                    }
        //                    return true;
        //                }
        //
        //                z = g;
        //
        //
        //                if (i > 0) {
        //
        //                    size_t h = std::min<size_t > (i, max_history);
        //                    size_t end = (i - 1) % h;
        //
        //                    //update histories
        //                    for (size_t r = 0; r < nop; r++) {
        //                        // std::cout<<r<<" "<<g[r] <<" - "<< pg[r]<<"\n";
        //                        dxs[r][end] = parameters[r]->GetValue() - px[r];
        //                        dgs[r][end] = g[r] - pg[r];
        //                    }
        //
        //                    std::valarray<T> p(h);
        //                    std::valarray<T>a(h);
        //
        //                    for (size_t j = 0; j < h; ++j) {
        //                        const size_t k = (end - j + h) % h;
        //                        p[k] = 1.0 / Dot(Column(dxs, k), Column(dgs, k));
        //
        //                        a[k] = p[k] * Dot(Column(dxs, k), z);
        //                        z -= a[k] * Column(dgs, k);
        //                    }
        //                    // Scaling of initial Hessian (identity matrix)
        //                    z *= Dot(Column(dxs, end), Column(dgs, end)) / Dot(Column(dgs, end), Column(dgs, end));
        //
        //                    for (size_t j = 0; j < h; ++j) {
        //                        const size_t k = (end + j + 1) % h;
        //                        const T b = p[k] * Dot(Column(dgs, k), z);
        //                        z += Column(dxs, k) * (a[k] - b);
        //
        //
        //                    }
        //
        //                }//end if(i>0)
        //
        //
        //
        //
        //                T descent = 0;
        //                for (size_t j = 0; j < nop; j++) {
        //                    px[j] = parameters[j]->GetValue();
        //                    x[j] = px[j];
        //                    pg[j] = g[j];
        //                    descent += z[j] * g[j];
        //                }//end for
        //
        //
        //                descent *= T(-1.0); // * Dot(z, g);
        //                if (descent > T(-0.000000001) * relative_tolerance /* tolerance relative_tolerance*/) {
        //
        //                    z = g;
        //                    iterations -= i;
        //                    i = 0;
        //                    step = 1.0;
        //                    descent = -1.0 * Dot(z, g);
        //                }//end if
        //
        //
        //
        //                bool down = false;
        //
        //                int ls;
        //
        //                fx_value = fx.GetValue();
        //
        //                ad::Variable<T>::SetRecordExpression(false);
        ////                std::cout << "expression recording is off...\n" << std::flush;
        //                for (ls = 0; ls < maxLineSearches_; ++ls) {
        //                    // Tentative solution, gradient and loss
        //                    std::valarray<T> nx = x - step * z;
        ////                                        std::cout<<"line search "<<ls<<"    \n";
        //                    bool bounded = false;
        //                    for (size_t j = 0; j < nop; j++) {
        //
        //                        if (this->has_constraints_m) {
        //
        //                            if (this->is_constrained_m[j] && (parameters[j]->GetValue() == lower_bounds_m[j]
        //                                    || parameters[j]->GetValue() == upper_bounds_m[j] || std::fabs(g[j]) <= tolerance)) {
        //                                bounded = true;
        //                                continue;
        //                            } else {
        //                                bounded = false;
        //                            }
        //
        //                            if (!this->is_constrained_m[j]) {
        //                                if (std::fabs(g[j]) <= tolerance) {
        //                                    bounded = true;
        //                                    continue;
        //                                } else {
        //                                    bounded = false;
        //                                }
        //
        //                            }
        //
        //                            if (this->is_constrained_m[j] && nx[j]< this->lower_bounds_m[j]) {
        //                                parameters[j]->SetValue(lower_bounds_m[j]);
        //                                continue;
        //                            } else if (this->is_constrained_m[j] && nx[j]> this->upper_bounds_m[j]) {
        //                                parameters[j]->SetValue(upper_bounds_m[j]);
        //                                continue;
        //                            } else {
        //                                parameters[j]->SetValue(nx[j]);
        //                                continue;
        //                            }
        //                        } else {
        //                            parameters[j]->SetValue(nx[j]);
        //                            continue;
        //                        }
        //                    }
        //
        //                    if (bounded) {
        //                        this->Print(fx, g, this->parameters_m, "Constrained Solution Found!");
        //                        return true;
        //                    }
        //
        //
        //                    this->CallObjectiveFunction(fx);
        //                    //this->CallGradient(fx, parameters, ng);
        //
        //                    //                    if (nfx.GetValue() != nfx.GetValue()) {
        //                    //                        return false;
        //                    //                    }
        //
        //                    //                    if (this->verbose_m) {
        //                    //                        std::cout << BOLD << "Line Search Value{" << ls << "}" << DEFAULT_IO << " = " << nfx << "\r";
        //                    //                    }
        //
        //
        //
        //
        //                    if (fx.GetValue() <= fx_value + tolerance * T(0.0001) * step * descent) { // First Wolfe condition
        //
        //                        ad::Variable<T>::SetRecordExpression(true);
        ////                        std::cout << "expression recording is on...\n" << std::flush;
        //                        this->CallObjectiveFunction(fx);
        //                        this->CallGradient(fx, parameters, ng);
        //
        //                        if ((-1.0 * Dot(z, ng) >= 0.9 * descent) || down) { // Second Wolfe condition
        //
        //
        //                            x = nx;
        //                            g = ng;
        //                            //                            fx = nfx;
        //                            this->function_value_m = fx.GetValue();
        //                            break;
        //                        } else {
        //                            ad::Variable<T>::SetRecordExpression(false);
        ////                            std::cout << "expression recording is off...\n" << std::flush;
        //                            step *= 100.0;
        //                        }
        //                    } else {
        //                        step /= 100.0;
        //                        down = true;
        //                    }
        //                }
        //
        //
        //
        //                if (ls == maxLineSearches_) {
        //                    std::cout << "Max line searches!\n";
        //                    return false;
        //                }
        //
        //            }
        //            return false;
        //        }

        bool NelderMead(std::vector<ad::Variable<T>* > &parameters, size_t iterations = 10000, T tolerance = (T(1e-5))) {

        }

        /**
         * A modified version of Charles Dubout's LBFGS algorithm released under the 
         * terms of GPL. 
         * http://www.idiap.ch/~cdubout/code/lbfgs.cpp
         * 
         * 
         * @param parameters
         * @param results
         * @param iterations
         * @param tolerance
         * @return 
         */
        bool LBFGS(std::vector<ad::Variable<T>* > &parameters, size_t iterations = 10000, T tolerance = (T(1e-4))) {

            size_t max_history = max_history_m;

            int maxLineSearches_ = 1000;



            std::valarray<T> x(parameters.size());
            //current gradient
            std::valarray<T> g(parameters.size());


            std::valarray<T> ng(x.size());

            //initial evaluation
            ad::Variable<T> fx(0.0);



            //Call the objective function and collect stats..
            this->CallObjectiveFunction(fx);
            this->function_value_m = fx.GetValue();
            //            ad::Variable<T> nfx(fx);
            //Historical evaluations
            std::valarray<T> px(parameters.size());
            std::valarray<T> pg(parameters.size());
            std::valarray<std::valarray<T> > dxs(std::valarray<T > (max_history), parameters.size());
            std::valarray<std::valarray<T> > dgs(std::valarray<T > (max_history), parameters.size());

            //set parameters
            for (size_t i = 0; i < g.size(); i++) {
                x[i] = parameters[i]->GetValue();
            }

            std::valarray<T> z(parameters.size());


            this->CallGradient(fx, parameters, g);


            T step = 0.1;
            T relative_tolerance;
            T norm_g;
            const size_t nop = this->active_parameters_m.size();
            std::valarray<T> p;
            std::valarray<T>a;
            for (int i = 0; i < iterations; ++i) {






                iteration_m = i + 1;

                norm_g = this->Norm(g);
                //
                //                // Backtracking using Wolfe's first condition (Armijo condition)
                //                step = i ? 1.0 : (1.0 / norm_g); //was just norm2

                relative_tolerance = tolerance * std::max<T > (T(1.0), norm_g);

                if (this->verbose_m && ((i % this->iprint_m) == 0)) {
                    this->Print(fx, g, parameters, "Verbose:\nMethod: L-BFGS");
                }

                if (norm_g < relative_tolerance) {

                    if (this->verbose_m) {
                        this->Print(fx, g, parameters, "Successful Convergence!");
                    }
                    return true;
                }

                z = g;


                if (i > 0) {

                    size_t h = std::min<size_t > (i, max_history);
                    size_t end = (i - 1) % h;

                    //update histories
                    for (size_t r = 0; r < nop; r++) {
                        // std::cout<<r<<" "<<g[r] <<" - "<< pg[r]<<"\n";
                        dxs[r][end] = parameters[r]->GetValue() - px[r];
                        dgs[r][end] = g[r] - pg[r];
                    }

                    p.resize(h);
                    a.resize(h);

                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end - j + h) % h;
                        p[k] = 1.0 / Dot(Column(dxs, k), Column(dgs, k));

                        a[k] = p[k] * Dot(Column(dxs, k), z);
                        z -= a[k] * Column(dgs, k);
                    }
                    // Scaling of initial Hessian (identity matrix)
                    z *= Dot(Column(dxs, end), Column(dgs, end)) / Dot(Column(dgs, end), Column(dgs, end));

                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end + j + 1) % h;
                        const T b = p[k] * Dot(Column(dgs, k), z);
                        z += Column(dxs, k) * (a[k] - b);


                    }

                }//end if(i>0)




                T descent = 0;
                for (size_t j = 0; j < nop; j++) {
                    px[j] = parameters[j]->GetValue();
                    x[j] = px[j];
                    pg[j] = g[j];
                    descent += z[j] * g[j];
                }//end for


                descent *= T(-1.0); // * Dot(z, g);
                if (descent > T(-0.0000000001) * relative_tolerance /* tolerance relative_tolerance*/) {

                    z = g;
                    iterations -= i;
                    i = 0;
                    step = 1.0;
                    descent = -1.0 * Dot(z, g);
                }//end if



                bool down = false;

                int ls;



                ad::Variable<T>::SetRecording(false);
                for (ls = 0; ls < maxLineSearches_; ++ls) {
                    // Tentative solution, gradient and loss
                    std::valarray<T> nx = x - step * z;

                    for (size_t j = 0; j < nop; j++) {

                        parameters[j]->SetValue(nx[j]);
                    }



                    this->CallObjectiveFunction(fx);

                    if (fx.GetValue() != fx.GetValue()) {
                        return false;
                    }

                    //                                        if (this->verbose_m) {
                    //                                            std::cout << BOLD << "Line Search Value{" << ls << "}" << DEFAULT_IO << " = " << fx << "                                \r";
                    //                                        }
                    //


                    if (fx.GetValue() <= this->function_value_m + tolerance * T(0.0001) * step * descent) { // First Wolfe condition

                        ad::Variable<T>::SetRecording(true);
                        this->CallObjectiveFunction(fx);
                        this->CallGradient(fx, parameters, ng);

                        if (down || (-1.0 * Dot(z, ng) >= 0.9 * descent)) { // Second Wolfe condition
                            x = nx;
                            g = ng;
                            //                            fx = fx;
                            this->function_value_m = fx.GetValue();
                            break;
                        } else {

                            step *= 10.0;
                        }
                    } else {
                        step /= 10.0;
                        down = true;
                    }
                }



                if (ls == maxLineSearches_) {
                    std::cout << "Max line searches!\n";
                    return false;
                }

            }
            return false;
        }


#ifdef HAVE_ADMB

        bool ADMB_Minimizer(std::vector<ad::Variable<T>* > &parameters, size_t iterations = 10000, T tolerance = (T(1e-4))) {
            const int nvar = parameters.size();
            double f;
            independent_variables x(1, nvar);
            dvector g(1, nvar);
            std::valarray<T> gradient(parameters.size());
            ad::Variable<T> fx;
            fmm fmc(nvar);
            fmc.crit = tolerance;
            fmc.iprint = 0;
            std::cout << fmc << "\n";
            //            fmc.dfn =0.01;
            //            exit(0);
            //            fmc.dfn    = 0.01;
            for (int i = 0; i < parameters.size(); i++) {
                x[i + 1] = parameters[i]->GetValue();
            }
            //            int iteration = this->iteration_m;
            this->CallObjectiveFunction(fx);
            f = fx.GetValue();
            this->CallGradient(fx, parameters, gradient);

            for (int i = 0; i < gradient.size(); i++) {
                g[i + 1] = gradient[i];
            }
            this->iteration_m = 0;
            //            this->verbose_m =false;
            while (fmc.ireturn >= 0) // Begin loop for minimization
            {
                this->iteration_m = fmc.itn;


                fmc.fmin(f, x, g); // Calls the minimization routine


                if (fmc.ireturn > 0) // Loop for evaluating the function and
                { // derivatives


                    for (int i = 0; i < parameters.size(); i++) {
                        *parameters[i] = x[i + 1];
                    }

                    if (fmc.ireturn == 1) {
                        Variable<T>::SetRecording(false);
                        this->CallObjectiveFunction(fx);
                        f = fx.GetValue();
                    } else {
                        Variable<T>::SetRecording(true);
                        this->CallObjectiveFunction(fx);
                        f = fx.GetValue();
                        this->Gradient(fx, parameters, gradient);

                        for (int i = 0; i < gradient.size(); i++) {
                            g[i + 1] = gradient[i];
                        }
                    }
                }

                                if (this->verbose_m && ((this->iteration_m % this->iprint_m) == 0)) {
                                    this->Print(fx, gradient, parameters, "Verbose:\nMethod: ADMB");
                                }

                //                iteration++;
            }

            return fmc.ireturn;

        }

#endif

#ifdef HAVE_GSL

        bool GSL_Multimin(std::vector<ad::Variable<T>* > &parameters, size_t iterations = 10000, T tolerance = (T(1e-3))) {

            std::string method = "unknown";
            const double initial_step_size = 0.01;
            const double line_search_tolerance = 1.0e-2;
            const double converged_gradient_norm = 1.0e-3;

            const gsl_multimin_fdfminimizer_type* minimizer_type
                    = gsl_multimin_fdfminimizer_vector_bfgs2;

            switch (this->minimizer_type_m) {

                case GSL_CONJUGATE_FR:
                    minimizer_type
                            = gsl_multimin_fdfminimizer_conjugate_fr;
                    method = "GSL_CONJUGATE_FR";
                    break;
                case GSL_CONJUGATE_PR:
                    minimizer_type
                            = gsl_multimin_fdfminimizer_conjugate_pr;
                    method = "GSL_CONJUGATE_PR";
                    break;
                case GSL_BFGS:
                    minimizer_type
                            = gsl_multimin_fdfminimizer_vector_bfgs;
                    method = "GSL_BFGS";
                    break;
                case GSL_BFGS2:
                    minimizer_type
                            = gsl_multimin_fdfminimizer_vector_bfgs2;
                    method = "GSL_BFGS2";
                    break;
                case GSL_STEEPEST_DESCENT:
                    minimizer_type
                            = gsl_multimin_fdfminimizer_steepest_descent;
                    method = "GSL_STEEPEST_DESCENT";
                    break;
                case NEWTON:
                    break;
                case DUBOUT_LBFGS:
                    break;


            }

            gsl_multimin_function_fdf my_function;
            my_function.n = parameters.size();
            my_function.f = ad::function_value_callback;
            my_function.df = ad::function_gradient_callback;
            my_function.fdf = ad::function_value_and_gradient_callback;
            my_function.params = reinterpret_cast<void*> (this);


            gsl_vector *x;
            x = gsl_vector_alloc(my_function.n);
            for (unsigned int i = 0; i < parameters.size(); i++) {
                gsl_vector_set(x, i, parameters[i]->GetValue());
            }


            gsl_multimin_fdfminimizer* minimizer
                    = gsl_multimin_fdfminimizer_alloc(minimizer_type, my_function.n);

            gsl_multimin_fdfminimizer_set(minimizer, &my_function, x,
                    initial_step_size, line_search_tolerance);


            size_t iter = 0;
            int status;
            this->CallObjectiveFunction(this->function_result_m);
            this->CallGradient(this->function_result_m, this->active_parameters_m, this->gradient_m);
            this->Print(this->function_result_m, this->gradient_m, parameters, "Verbose:\nMethod: " + method);

            do {

                ++iter;

                status = gsl_multimin_fdfminimizer_iterate(minimizer);


                if (status != GSL_SUCCESS) break;

                status = gsl_multimin_test_gradient(minimizer->gradient, converged_gradient_norm);


                if (this->verbose_m && ((iter % this->iprint_m) == 0)) {
                    this->Print(this->function_result_m, this->gradient_m, parameters, "Verbose:\nMethod: " + method);
                }

            } while (status == GSL_CONTINUE && iter < iterations);


            gsl_multimin_fdfminimizer_free(minimizer);
            gsl_vector_free(x);


            if (status == GSL_SUCCESS) {
                this->Print(this->function_result_m, this->gradient_m, parameters, "Verbose:\nMethod: GSL");
                return true;
            } else {
                this->Print(this->function_result_m, this->gradient_m, parameters, "Verbose:\nMethod: GSL");
                return false;
            }

        }
#endif

        /**
         * Compute the dot product of two vectors.
         * @param a
         * @param b
         * @return 
         */
        const T Dot(const std::valarray<T> &a, const std::valarray<T> &b) {
            T ret = 0;
            for (size_t i = 0; i < a.size(); i++) {
                ret += a[i] * b[i];
            }
            return ret;
        }

        /**
         * Compute the Norm of the vector v.
         *  
         * @param v
         * @return 
         */
        const T Norm(std::valarray<T> &v) {

            T ret = (T) 0.0;
            unsigned int i;
            for (i = 0; i < v.size(); i++) {
                ret += v[i] * v[i];

            }
            return std::sqrt(ret);
        }

        const std::valarray<T > GetGradient() {

        }

        /**
         * Get current time in milliseconds. Used for runtime statistics.
         */
        int GetMilliCount() {
#if defined(WIN32) || defined(WIN64) 
            retrun GetTickCount();
#else
            timeb tb;
            ftime(&tb);
            int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
            return nCount;
#endif
        }

        /**
         * Print current minimizer state to stdout.
         * 
         * @param ret
         * @param gradient
         * @param parameters
         * @param message
         */
        void Print(const ad::Variable<T> &ret, const std::valarray<T> &gradient, const std::vector<ad::Variable<T>* > &parameters, std::string message = "") {


            std::cout << BOLD << message << DEFAULT_IO << std::endl;
            std::cout << "Phase: " << BOLD << this->phase_m
                    << DEFAULT_IO << " of " << BOLD << this->max_phase_m
                    << DEFAULT_IO << "\n";
            std::cout << "Iteration: " << BOLD
                    << this->iteration_m << DEFAULT_IO << std::endl;
            std::cout << "Function Calls: "
                    << BOLD << this->function_calls_m << DEFAULT_IO << " (" << this->unrecorded_calls_m << " unrecorded line searches)" << std::endl;
            std::cout << "Average Time in Objective Function: "
                    << BOLD
                    << static_cast<long> (this->average_time_in_user_function_m)
                    << " ms\n" << DEFAULT_IO;
            std::cout << "Average Time Calculating Gradients: " << BOLD
                    << static_cast<long> (this->average_time_in_grad_calc_m)
                    << " ms\n" << DEFAULT_IO;
            int prec = std::cout.precision();
            std::cout.precision(50);
            std::cout << "Function Value = " << BOLD << ret
                    << std::endl << DEFAULT_IO;
            std::cout.precision(prec);
            std::cout << "Active Parameters: " << BOLD << parameters.size()
                    << DEFAULT_IO << std::endl;
            std::cout << std::scientific;
            std::cout << "Tolerance: " << BOLD << this->GetTolerance() << DEFAULT_IO
                    << std::endl;
            //            std::cout << "Expression Size: " << BOLD << ret.ExpressionSize() << DEFAULT_IO
            //                    << std::endl;
            std::cout << "Maximum Gradient Component Magnitude: " << BOLD
                    << max_c << DEFAULT_IO << std::endl;

            double number_good = 0.0;
            for (size_t i = 0; i < this->active_parameters_m.size(); i++) {
                bool is_max_c = false;


                if (std::fabs(gradient[i]) == max_c) {

                    is_max_c = true;
                }

                //make the output nice and clean by converting to string and 
                //back to double so we can use scientific format for
                //multi precision types.
                std::stringstream ss;
                ss << gradient[i];
                std::string grad_string = ss.str();

                ss.str("");
                ss << parameters[i]->GetValue();
                std::string param_string = ss.str();

                ss.str("");
                ss << max_c;
                std::string max_c_string = ss.str();

                if ((i % 10) == 0) {
                    std::cout << BOLD << std::left << std::setw(20)
                            << "Parameter"
                            << std::left << std::setw(20) << "Value"
                            << "Gradient["
                            << BOLD << RED << "" << toValue<double>(max_c_string) << "" << DEFAULT_IO
                            << BOLD << std::left << std::setw(10) << "]"
                            << "Bounded" << std::endl;
                }
                std::cout << BLUE << BOLD << std::left << std::setw(20)
                        << this->active_parameters_m[i]->GetName() << DEFAULT_IO << std::left
                        << std::setw(20) << toValue<double>(param_string) << std::left;


                if (is_max_c) {
                    std::cout << RED << BOLD << std::setw(31) << toValue<double>(grad_string) << DEFAULT_IO;
                } else {
                    if (std::fabs(toValue<double>(grad_string)) <= this->tolerance_m) {
                        number_good++;
                        std::cout << BOLD << GREEN << std::setw(31) << toValue<double>(grad_string) << DEFAULT_IO;
                    } else {
                        std::cout << std::setw(31) << toValue<double>(grad_string);
                    }
                }

                if (this->active_parameters_m[i]->IsBounded()) {
                    std::cout << "True" << "[" << this->active_parameters_m[i]->GetMinBoundary() << ","
                            << this->active_parameters_m[i]->GetMaxBoundary() << "]" << std::endl;
                } else {
                    std::cout << "False" << std::endl;
                }

            }

            //            double percent = 100 * (number_good / (double) this->active_parameters_m.size());
            //            std::cout << BOLD << "|";
            //            for (double i = 0; i < percent - 1; i++) {
            //                std::cout << GREEN << "=";
            //            }
            //            if (percent > 1)
            //                std::cout << GREEN << ">";
            //            for (int i = 0; i < (100 - percent); i++) {
            //                std::cout << RED << "*";
            //            }
            //            std::cout << DEFAULT_IO;
            //            std::cout << BOLD << "|" << DEFAULT_IO;
            //            std::cout.precision(prec);
            //
            //            std::cout << "[" << (int) percent << "%]\n";
            std::cout << std::endl;
            std::cout << std::fixed;
        }

        /**
         * returns the a column of a matrix as a std::valarray.
         * @param matrix
         * @param column
         * @return 
         */
        const std::valarray<T> Column(std::valarray<std::valarray<T> > &matrix, size_t column) {

            std::valarray<T> ret(this->active_parameters_m.size());

            for (int i = 0; i < ret.size(); i++) {
                ret[i] = matrix[i][column];
            }
            return ret;
        }

        template<class TT >
        TT toValue(const std::string & val) const {
            std::istringstream ss(val);
            //   ss<<val;
            TT result;
            return ss >> result ? result : 0;
        }

    };



    //gsl
#ifdef HAVE_GSL

    extern "C"
    double function_value_callback(const gsl_vector* x, void* params) {
        FunctionMinimizer<double>* fm = reinterpret_cast<FunctionMinimizer<double>*> (params);
        for (int i = 0; i < fm->active_parameters_m.size(); i++) {
            fm->active_parameters_m[i]->SetValue(x->data[i]);
        }
        Variable<double>::SetRecording(false);
        Variable<double> fx;
        fm->CallObjectiveFunction(fx);
        return fx.GetValue();

    }

    extern "C"
    void function_gradient_callback(const gsl_vector* x, void* params, gsl_vector* gradJ) {
        FunctionMinimizer<double>* fm = reinterpret_cast<FunctionMinimizer<double>*> (params);
        for (int i = 0; i < fm->active_parameters_m.size(); i++) {
            fm->active_parameters_m[i]->SetValue(x->data[i]);
        }
        Variable<double>::SetRecording(true);
        Variable<double> fx;
        std::valarray<double> g(fm->active_parameters_m.size());
        fm->CallObjectiveFunction(fx);

        fm->CallGradient(fx, fm->active_parameters_m, g);


        for (int i = 0; i < fm->active_parameters_m.size(); i++) {
            gradJ->data[i] = g[i];
        }

    }

    extern "C"
    void function_value_and_gradient_callback(const gsl_vector* x, void* params,
            double* J, gsl_vector* gradJ) {
        FunctionMinimizer<double>* fm = reinterpret_cast<FunctionMinimizer<double>*> (params);
        for (int i = 0; i < fm->active_parameters_m.size(); i++) {
            fm->active_parameters_m[i]->SetValue(x->data[i]);
        }
        Variable<double>::SetRecording(true);
        ad::Variable<double> fx;
        std::valarray<double> g(fm->active_parameters_m.size());
        fm->CallObjectiveFunction(fx);

        fm->CallGradient(fx, fm->active_parameters_m, g);

        *J = fx.GetValue();

        for (int i = 0; i < fm->active_parameters_m.size(); i++) {
            gradJ->data[i] = g[i];
        }
    }

#endif

}


#endif	/* FUNCTIONMINIMIZER_HPP */

