/* 
 * File:   Portfolio.hpp
 * Author: matthewsupernaw
 *
 * Created on June 2, 2014, 9:01 AM
 */

#ifndef PORTFOLIO_HPP
#define	PORTFOLIO_HPP

#include <vector>
#include <string>
#include <curl/curl.h>
#include "IOStream.hpp"

template<class REAL_T>
class Quote {
    std::string date_m;
    REAL_T open_m;
    REAL_T low_m;
    REAL_T high_m;
    REAL_T close_m;
    REAL_T volume_m;
    REAL_T adj_close_m;


public:

    static uint32_t HISTORY_SIZE;

    Quote(std::string date,
            REAL_T open,
            REAL_T high,
            REAL_T low,
            REAL_T close,
            REAL_T volume,
            REAL_T adj) : date_m(date),
    open_m(open),
    high_m(high),
    low_m(low),
    close_m(close),
    volume_m(volume),
    adj_close_m(adj) {

    }

    REAL_T GetAdj_close() const {
        return adj_close_m;
    }

    void SetAdj_close(REAL_T adj_close) {
        this->adj_close_m = adj_close;
    }

    REAL_T GetClose() const {
        return close_m;
    }

    void SetClose(REAL_T close) {
        this->close_m = close;
    }

    std::string GetDate() const {
        return date_m;
    }

    void SetDate(std::string date) {
        this->date_m = date;
    }

    REAL_T GetHi() const {
        return high_m;
    }

    void SetHi(REAL_T hi) {
        this->high_m = hi;
    }

    REAL_T GetLow() const {
        return low_m;
    }

    void SetLow(REAL_T low) {
        this->low_m = low;
    }

    REAL_T GetOpen() const {
        return open_m;
    }

    void SetOpen(REAL_T open) {
        this->open_m = open;
    }

    REAL_T GetVolume() const {
        return volume_m;
    }

    void SetVolume(REAL_T volume) {
        this->volume_m = volume;
    }

    std::string ToString() {
        std::stringstream ss;
        ss << this->GetDate() << ",";
        ss << this->GetOpen() << ",";
        ss << this->GetHi() << ",";
        ss << this->GetLow() << ",";
        ss << this->GetClose() << ",";
        ss << this->GetVolume() << ",";
        ss << this->GetAdj_close();
        return ss.str();
    }



};

template<class REAL_T>
uint32_t Quote<REAL_T>::HISTORY_SIZE = 120;

template<class REAL_T>
class Asset {
public:
    std::string name;
    std::vector<REAL_T> price_history;
    REAL_T expected_return; //mean return rate
    REAL_T mean; //mean return rate
    REAL_T return_variance;
    REAL_T return_volatility;

    std::list<Quote<REAL_T> > history;
public:

    Asset(const std::string &name, std::string data) : name(name) {

        std::vector<std::string> lines;
        util::Tokenize(data, lines, "\n");
        int count =0;
        
        for (int i = 1; i < lines.size()-1; i++) {

            if (count < Quote<REAL_T>::HISTORY_SIZE) {
                std::vector<std::string> quote;
                util::Tokenize(lines[i], quote, ",");
//                std::cout<</lines[i]<<"\n";
                history.push_front(
                        Quote<REAL_T > (quote[0],
                        util::StringToNumber<REAL_T > (quote[1]),
                        util::StringToNumber<REAL_T > (quote[2]),
                        util::StringToNumber<REAL_T > (quote[3]),
                        util::StringToNumber<REAL_T > (quote[4]),
                        util::StringToNumber<REAL_T > (quote[5]),
                        util::StringToNumber<REAL_T > (quote[6])));
                count++;

            }else{
                break;
            }
        }
        
    
        std::cout << history.size() << " data points.\n";
//                typedef typename std::list<Quote<REAL_T> >::iterator iter;
//        
//                for (iter i = history.begin(); i != history.end(); i++) {
//                    std::cout << i->ToString() << "\n";
//                }
//    exit(0);


        //        mean = 0;
        //        expected_return = 0;
        //        for (int i = 0; i < price_history.size(); i++) {
        //            mean += price_history[i];
        //            if (i > 0) {
        //                expected_return += (price_history[0] - price_history[i]) / price_history[0];
        //            }
        //        }
        //
        //        expected_return /= price_history.size();
        //        mean /= price_history.size();
        //        REAL_T expected;
        //        for (int i = 1; i < price_history.size(); i++) {
        //
        //            expected += (((price_history[0] - price_history[i]) / price_history[0]) - expected_return)*
        //                    (((price_history[0] - price_history[i]) / price_history[0]) - expected_return);
        //
        //        }
        //
        //        return_variance = expected / (price_history.size() - 1);
        //        return_volatility = std::sqrt(return_variance);

    }

    REAL_T GetExpectedReturn() const {
        return expected_return;
    }

    void SetExpectedReturn(REAL_T expected_return) {
        this->expected_return = expected_return;
    }

    std::vector<REAL_T> GetPriceHistory() const {
        return price_history;
    }

    void SetPriceHistory(std::vector<REAL_T> price_history) {
        this->price_history = price_history;
    }

    REAL_T GetReturnVariance() const {
        return return_variance;
    }

    void SetReturnVariance(REAL_T return_variance) {
        this->return_variance = return_variance;
    }

    REAL_T GetReturnVolatility() const {
        return return_volatility;
    }

    void SetReturnVolatility(REAL_T return_volatility) {
        this->return_volatility = return_volatility;
    }



private:

};

template<class REAL_T, class EVAL_T = REAL_T>
class Portfolio {
    //initial
    std::vector<Asset<REAL_T> > assets;
    std::vector<EVAL_T> correlation_coefficients; //n x n matrix
    EVAL_T risk_tolerance; //q

    //runtime
    EVAL_T expected_return;
    EVAL_T return_variance;
    EVAL_T return_volatility; //standard deviation
    std::vector<EVAL_T> weights;

    EVAL_T objective_function_value;

    void PuildCorrelationMatrix() {

        for (size_t i = 0; i < assets.size(); i++) {
            for (size_t j = 0; j < assets.size(); j++) {
                correlation_coefficients[i * assets.size() + j] =
                        covariance(assets[i], assets[j], assets[i].price_history.size()) / (assets[i].return_variance * assets[j].return_variance);
            }
        }
    }

    const REAL_T covariance(const Asset<REAL_T> &a, const Asset<REAL_T> &b, const size_t &n) {
        REAL_T sum = 0;
        for (size_t i = 0; i < n; i++) {
            sum += ((a.price_history[i] - a.expected_return)*(b.price_history[i] - b.expected_return)) / n;
        }
    }

public:

    void ObjectiveFunction(EVAL_T &f) {
        //        expected_return = 0;
        //        return_variance = 0;
        //
        //        for (size_t i = 0; i < assets.size(); i++) {
        //            expected_return += weights[i] * assets[i].expected_return;
        //            return_variance += (weights[i] * weights[i])*(assets[i].return_variance * assets[i].return_variance);
        //            for (size_t j = 0; j < assets.size(); j++) {
        //                if (j ! = i) {
        //                    return_variance +=
        //                            (weights[i] *
        //                            weights[j] *
        //                            assets[i].return_volatility *
        //                            assets[j].return_volatility *
        //                            correlation_coefficients[i * this->assets.size() + j]);
        //                }
        //
        //            }
        //
        //            f = return_variance - this->risk_tolerance*expected_return;

        //        }


    }


private:


};


#endif	/* PORTFOLIO_HPP */

