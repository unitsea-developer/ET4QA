/* 
 * File:   CatchAtAge2.hpp
 * Author: matthewsupernaw
 *
 * Created on April 2, 2014, 1:25 PM
 */

#ifndef CATCHATAGE2_HPP
#define	CATCHATAGE2_HPP
#include <sstream>

#include "IOStream.hpp"
#include "FunctionMinimizer.hpp"
#include "ET4AD.hpp"

template<class T>
class CatchAtAge : public ad::FunctionMinimizer<T> {
    typedef ad::Variable<T> variable;

   
    std::string input_file_path;

    /*
    
     init_int nyrs
   init_int nages
   init_matrix obs_catch_at_age(1,nyrs,1,nages)
   init_vector effort(1,nyrs)
   init_number M
   vector relwt(2,nages);
     */

    int nyrs;
    int nages;
    std::vector<T> obs_catch_at_age;
    std::vector<T> effort;
    T M;
    std::vector<T> relwt;

public:

    //Estimated
    variable log_q;
    variable log_popscale;
    std::vector<variable> log_sel_coff;
    std::vector<variable> log_relpop;
    std::vector<variable> effort_devs;
    //Runtime
    std::vector<variable> log_sel;
    std::vector<variable> log_initpop;
    std::vector<variable> F;
    std::vector<variable> Z;
    std::vector<variable> S;
    std::vector<variable> N;
    std::vector<variable> C;
    std::vector<variable> predicted_N;
    std::vector<variable> ratio_N;


    void Initialize() {
         StreamedDataFile<T> input_file;
        input_file_path = "/Users/matthewsupernaw/NetBeansProjects/ET4AD/catage.dat__";
        input_file.Parse(input_file_path);

        this->nyrs = static_cast<int> (input_file.Next());
        this->nages = static_cast<int> (input_file.Next());

        this->obs_catch_at_age = std::vector<T > (nyrs * nages);
        this->effort = std::vector<T > (nyrs);
        this->log_q.SetName("log_q");
        this->log_q = -1.0;
        this->Register(log_q);

        this->log_popscale.SetName("log_popscale");
        this->log_popscale = 5.0;
        this->Register(log_popscale);

        this->log_initpop = std::vector<variable > (nyrs + nages - 1);
        for (int i = 0; i< this->log_initpop.size(); i++) {
            this->log_initpop[i] = variable();
        }

        std::stringstream ss;

        this->log_sel_coff = std::vector<variable > (nages - 1);
        for (int i = 0; i< this->log_sel_coff.size(); i++) {
            ss.str("");
            ss << "log_sel_coff[" << i << "]";
            log_sel_coff[i] = variable();
            log_sel_coff[i].SetBounds(-15.0, 15.0);
            log_sel_coff[i].SetName(ss.str());
            this->Register(log_sel_coff[i],2);
        }

        this->log_relpop = std::vector<variable > (nyrs + nages - 1);
        for (int i = 0; i< this->log_relpop.size(); i++) {
            ss.str("");
            ss << "log_relpop[" << i << "]";
            log_relpop[i] = variable();
            log_relpop[i].SetBounds(-15.0, 15.0);
            log_relpop[i].SetName(ss.str());
            this->Register(log_relpop[i],2);
        }


        this->effort_devs = std::vector<variable > (nyrs);
        for (int i = 0; i< this->effort_devs.size(); i++) {
            ss.str("");
            ss << "effort_devs[" << i << "]";
            effort_devs[i] = variable(0.0);
            effort_devs[i].SetBounds(-5.0, 5.0);
            effort_devs[i].SetName(ss.str());
            this->Register(effort_devs[i],3);

        }
        this->log_sel = std::vector<variable > (nages);
        this->predicted_N = std::vector<variable > (nages);
        this->ratio_N = std::vector<variable > (nages);

        this->F = std::vector<variable > (nyrs * nages);
        this->Z = std::vector<variable > (nyrs * nages);
        this->S = std::vector<variable > (nyrs * nages);
        this->N = std::vector<variable > (nyrs * nages);
        this->C = std::vector<variable > (nyrs * nages);


        for (int i = 0; i < nyrs; i++) {
            for (int j = 0; j < nages; j++) {
                this->obs_catch_at_age[i * nages + j] = input_file.Next();
                this->F[i * nages + j] = variable();
                this->Z[i * nages + j] = variable();
                this->S[i * nages + j] = variable();
                this->N[i * nages + j] = variable();
                this->C[i * nages + j] = variable();
            }
        }

        for (int i = 0; i < nyrs; i++) {
            this->effort[i] = input_file.Next();
        }

        this->M = input_file.Next();

        this->relwt = std::vector<T > (nages);
        T w = 1.0;
        T maxw = 0;
        for (int i = 1; i < nages; i++) {
            relwt[i] = w;
            w++;
            relwt[i] = std::pow(relwt[i], (T) .5);
            if (relwt[i] > maxw) {
                maxw = relwt[i];
            }
        }

        for (int i = 1; i < nages; i++) {
            relwt[i] = relwt[i] / maxw;
        }

    }

    void GetMortalityAndSurvivalRates() {

        int i, j;
        // calculate the selectivity from the sel_coffs
        for (j = 0; j < nages - 1; j++) {
            log_sel[j] = log_sel_coff[j];

        }

        // the selectivity is the same for the last two age classes
        log_sel[nages - 1] = log_sel_coff[nages - 2];


        // This is the same as F(i,j)=exp(q)*effert(i)*exp(log_sel(j));
        //F = outer_prod(mfexp(log_q) * effort, mfexp(log_sel));
        for (i = 0; i < nyrs; i++) {
            for (j = 0; j < nages; j++) {
                F[i * nages + j] = (std::mfexp(log_q) * effort[i]) * std::mfexp(log_sel[j]);
                //                std::cout<<F[i * nages + j]<<" "<<F[i * nages + j].wrt(this->log_q)<<"\n";
                //                exit(0);
            }
        }

        //        exit(0);


        if (this->Phase() == 3) {//active(effort_devs)) {
            for (i = 0; i < nyrs; i++) {
                for (j = 0; j < nages; j++) {
                    F[i * nages + j] = F[i * nages + j] * std::mfexp(effort_devs[i]);
                }
            }
        }


        // get the total mortality
        for (i = 0; i < nyrs; i++) {
            for (j = 0; j < nages; j++) {
                Z[i * nages + j] = F[i * nages + j] + M;
                S[i * nages + j] = std::mfexp(static_cast<T>(-1.0) * Z[i * nages + j]);
            }

        }
    }

    void GetNumberAtAge() {

        int i, j;
        for (int i = 0; i < log_initpop.size(); i++) {
            log_initpop[i] = log_relpop[i] + log_popscale;
        }

        for (i = 0; i < nyrs; i++) {
            N[i * nages] = std::mfexp(log_initpop[i]);
        }

        for (j = 1; j < nages; j++) {
            N[j] = std::mfexp(log_initpop[(nyrs) + j - 1]);
        }

        for (i = 0; i < nyrs - 1; i++) {
            for (j = 0; j < nages - 1; j++) {
                N[(i + 1) * nages + (j + 1)] = N[i * nages + j] * S[i * nages + j];
            }
        }

        // calculated predicted numbers at age for next year
        //        for (j = 0; j < nages - 1; j++) {
        //            predicted_N[j + 1] = N[((nyrs - 1) * nages) + j] * S[((nyrs - 1) * nages) + j];
        //            //            std::cout << predicted_N[j+1] << " "<<N[((nyrs-1)*nages ) + j] <<" "<<S[((nyrs-1)*nages ) + j]<<"\n";
        //            ratio_N[j + 1] = predicted_N[j + 1] / N[((nyrs - 1) * nages) + (j + 1)];
        //            //            std::cout<< ratio_N[j + 1]<<" ";
        //        }
        //        exit(0);
        // calculated predicted Biomass for next year for
        // adjusted profile likelihood
        //        pred_B = (predicted_N * relwt);

    }

    void GetCatchAtAge() {

        for (int i = 0; i < C.size(); i++) {
            C[i] = (F[i] / Z[i])*(((T) 1.0 - S[i]) * N[i]);
        }
    }

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

    void ObjectiveFunction(variable & f) {
        //
//                if(this->function_calls_m == 251){
//                    exit(0);
//                }
        f = 0.0;

        GetMortalityAndSurvivalRates();
        GetNumberAtAge();
        GetCatchAtAge();

        f += (T) .01 * norm2(log_relpop);
        //        std::cout << this->log_q << ":" << this->log_popscale << " " << f << "---" << f.wrt(this->log_q) << "\n";
        //        exit(0);

        variable avg_F = (T) 0.0;
        for (int i = 0; i < F.size(); i++) {
            avg_F += F[i];
        }
        avg_F /= (double) F.size();


        if (this->Phase() == this->max_phase_m) {

            // a very small penalty on the average fishing mortality
            f += (T) .001 * (std::log(avg_F / (T) .2) * std::log(avg_F / (T) .2));
        } else {
            f += (T) 1000. * (std::log(avg_F / (T) .2) * std::log(avg_F / (T) .2));
        }

        variable sum;
        for (int i = 0; i < C.size(); i++) {
            sum += ((C[i] - obs_catch_at_age[i])*(C[i] - obs_catch_at_age[i])) / ((T) 0.01 + C[i]);
            //            std::cout << this->log_q << ":" << this->log_popscale << " " << f << "---" << sum.wrt(this->log_q) << "\n";
        }

        f += (T) 0.5 * T(C.size() + effort_devs.size()) * std::log(sum + (T) 0.1 * norm2<variable > (effort_devs));
    
        //std::cout<<f.ExpressionSize()<<"\n";
        //        std::cout << this->log_q << ":" << this->log_popscale << " " << f << "---" << f.wrt(this->log_q) << "\n";
        //        exit(0);
//        f.Diff(this->log_popscale);
//        std::cout<<ad::Variable<T>::misses_g<<"\n";
//       std::cout<<"current ="<< f.gv.size()<<"\n";
//                std::cout<<f.Diff(this->log_popscale)<<" == "<<f.WRT(this->log_popscale)<<"\n";
    }

    template<class TT>
    const TT norm2(std::vector<TT> &vect) {
        TT ret = 0.0;
        for (int i = 0; i < vect.size(); i++) {
            ret = ret+ vect[i] * vect[i];
        }
        return ret;
    }

    void Finalize() {
       
        std::cout << BOLD << "Estimated number of fish:\n" << DEFAULT_IO;
        for (int i = 0; i < nyrs; i++) {
            for (int j = 0; j < nages; j++) {
                std::cout << N[(i) * nages + (j)] << "\t";
            }
            std::cout << std::endl;
        }

        std::cout << BOLD << "\nEstimated catch:\n" << DEFAULT_IO;
        for (int i = 0; i < nyrs; i++) {
            for (int j = 0; j < nages; j++) {
                std::cout << C[(i) * nages + (j)] << "\t";
            }
            std::cout << std::endl;
        }

        std::cout << BOLD << "\nObserved Catch:\n" << DEFAULT_IO;
        for (int i = 0; i < nyrs; i++) {
            for (int j = 0; j < nages; j++) {
                std::cout << this->obs_catch_at_age[(i) * nages + (j)] << "\t";
            }
            std::cout << std::endl;
        }


        std::cout << BOLD << "\nEstimated Mortality:\n" << DEFAULT_IO;
        for (int i = 0; i < nyrs; i++) {
            for (int j = 0; j < nages; j++) {
                std::cout << F[(i) * nages + (j)] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::flush;
        //        sleep(1);
        std::valarray<T> gradient = this->CalculateGradient();

        std::ofstream out;
        out.open("gradient.data");

        out << "ParName\tValue\tGradient\n";
        for (int i = 0; i < this->GetActiveParameters().size(); i++) {
            out << this->GetActiveParameters().at(i)->GetName()
                    << "\t" << this->GetActiveParameters().at(i)->GetValue()
                    << "\t" << gradient[i] << "\n";
        }

        //        std::valarray<std::valarray<double> > hessian = this->EstimatedHessian();



    }



private:




};


#endif	/* CATCHATAGE2_HPP */

