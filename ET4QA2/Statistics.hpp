/* 
 * File:   Statistics.hpp
 * Author: matthewsupernaw
 *
 * Created on June 10, 2014, 8:10 AM
 */

#ifndef STATISTICS_HPP
#define	STATISTICS_HPP

#include <vector>

#include "BigFloat.hpp"

namespace ad {

    template<class T>
    const T Sum(const std::vector<T> &y) {
        T sum = 0;
        for (uint32_t i = 0; i < y.size(); i++) {
            sum += y[i];
        }
        return sum;
    }

    template<class T>
    const T Mean(const std::vector<T> &y) {
        return Sum<T > (y) / static_cast<T>(y.size());
    }

    /**
     * Compute the geometric mean. 
     * 
     * The n-th root of y_0*y_1...y_n
     * 
     * @param y
     * @param substitution - value to replace zeros, default is 1.0.
     * @return 
     */
    template<class T>
    const T GeometricMean(const std::vector<T> &y, T substitution = 1.0) {
        T product = 1.0;
        for (uint32_t i = 0; i < y.size(); i++) {
            if (y[i] != 0.0) {
                product *= y[i];
            } else {
                product *= substitution;
            }
        }

        return std::pow(product, 1.0 / y.size());

    }

    template<class T>
    const T HarmonicMean(const std::vector<T> &y) {

        T sum = 0;
        for (uint32_t i = 0; i < y.size(); i++) {
            sum += (1.0 / y[i]);
        }
        return y.size() / sum;

    }

    template<class T>
    const T Median(const std::vector<T> &y) {
        std::vector<T> z = y;
        std::sort(z.begin(), z.end());
        size_t center = y.size() / 2;
        if (y.size() % 2 == 0) {
            return (z[center] + z[center - 1]) / 2;
        } else {
            return z[center];
        }
    }

    /**
     * Return the first quartile in the set.<br>
     *  
     *<br>Method:<br>
     *1. Use the median to divide the ordered data set into two halves. <br>
     *   Do not include the median in either half.<br>
     * 
     *2. The first quartile value is the median of the lower half of the data. <br>
     * 
     * @param y
     * @return 
     */
    template<class T>
    const T Q1(const std::vector<T> &y) {
        std::vector<T> z = y;
        std::sort(z.begin(), z.end());
        size_t center = y.size() / 2;
        if (y.size() % 2 == 0) {
            return (z[(center) / 2] + z[((center) / 2) - 1]) / 2;
        } else {
            return z[center / 2];
        }
    }

    /**
     * Return the first quartile in the set.<br>
     *  
     *<br>Method:<br>
     *1. Use the median to divide the ordered data set into two halves. <br>
     *   Do not include the median in either half.<br>
     * 
     *2. The third quartile value is the median of the lower half of the data. <br>
     * 
     * @param y
     * @return 
     */
    template<class T>
    const T Q3(const std::vector<T> &y) {
        std::vector<T> z = y;
        std::sort(z.begin(), z.end());
        size_t center = y.size() / 2;
        if (y.size() % 2 == 0) {
            return (z[(3 * center) / 2] + z[(3 * center) / 2 - 1]) / 2;
        } else {
            return z[center + 1 + center / 2];
        }
    }

    template<class T>
    const T Variance(const std::vector<T> &y) {
        T mean = Mean<T > (y);

        T sum = 0;
        for (uint32_t i = 0; i < y.size(); i++) {
            sum += (mean - y[i])*(mean - y[i]);
        }
        return sum / static_cast<T>(y.size());
    }

    template<class T>
    const T StandardDeviation(const std::vector<T> &y) {
        return std::sqrt(Variance<T > (y));
    }

    template<class T>
    const T Min(const std::vector<T> &y) {
        T min = std::numeric_limits<T>::max();

        for (uint32_t i = 0; i < y.size(); i++) {
            if (y[i] < min) {
                min = y[i];
            }
        }
        return min;
    }

    template<class T>
    const T Max(const std::vector<T> &y) {
        T max = std::numeric_limits<T>::min();

        for (uint32_t i = 0; i < y.size(); i++) {
            if (y[i] > max) {
                max = y[i];
            }
        }
        return max;
    }

    template<class T>
    const T Covariance(const std::vector<T> &x, const std::vector<T> &y) {
        T xmean = Mean<T > (x);
        T ymean = Mean<T > (y);

        T sum = 0;
        for (uint32_t i = 0; i < y.size(); i++) {
            sum += ((x[i] - xmean)*(y[i] - ymean));
        }
        return sum / static_cast<T>(x.size());

    }

    /**
     * The correlation coefficient will vary from -1 to +1.
     * A -1 indicates perfect negative correlation, and +1 indicates 
     * perfect positive correlation. 
     * @param x
     * @param y
     * @return 
     */
    template<class T>
    const T CorrelationCoefficient(const std::vector<T> &x, const std::vector<T> &y) {
        T covariance = Covariance<T > (x, y);
        std::cout << covariance << " / " << (StandardDeviation<T > (x) * StandardDeviation<T > (y)) << "\n";
        return covariance / (StandardDeviation<T > (x) * StandardDeviation<T > (y));
    }

    /**
     * Range of the data set y.
     * 
     * @param y
     * @param min
     * @param max
     */
    template<class T>
    void Range(const std::vector<T> &y, const T &min, const T &max) {
         max = std::numeric_limits<T>::min();
         min = std::numeric_limits<T>::max();
        for (uint32_t i = 0; i < y.size(); i++) {
            if (y[i] > max) {
                max = y[i];
            }
            if (y[i] < min) {
                min = y[i];
            }
        }
    }


}

#endif	/* STATISTICS_HPP */

