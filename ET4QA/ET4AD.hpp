
/*!
 *  Software to compute derivatives. 
 */

/*!
 *   This library is dual-licensed: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 3 as 
 *   published by the Free Software Foundation. For the terms of this 
 *   license, see licenses/gpl_v3.txt or <http://www.gnu.org/licenses/>.
 *
 *   You are free to use this library under the terms of the GNU General
 *   Public License, but WITHOUT ANY WARRANTY; without even the implied 
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *  Alternatively, you can license this library under a commercial
 *  license. 
 *
 *               ET4AD Commercial License (ET4ADCL)
 *               ==================================
 * ---------------------------------------------------------------------
 * 
 *  The license is non-exclusively granted to a single person or company,
 *  per payment of the license fee, for the lifetime of that person or
 *  company. The license is non-transferable.
 * 
 *  The ET4ADCL grants the licensee the right to use ET4AD in commercial
 *  closed-source projects. Modifications may be made to ET4AD with no
 *  obligation to release the modified source code. ET4AD (or pieces
 *  thereof) may be included in any number of projects authored (in whole
 *  or in part) by the licensee.
 * 
 *  The licensee may use any version of ET4AD, past, present or future,
 *  as is most convenient. This license does not entitle the licensee to
 *  receive any technical support, updates or bug fixes, except as such are
 *  made publicly available to all ET4AD users.
 * 
 *  The licensee may not sub-license ET4AD itself, meaning that any
 *  commercially released product containing all or parts of ET4AD must
 *  have added functionality beyond what is available in ET4AD;
 *  ET4AD itself may not be re-licensed by the licensee.
 * 
 *  To obtain a commercial license agreement, contact:
 * 
 *  Matthew Supernaw
 *  msupernaw@gmail.com
 * 
 */


/*! 
 * File:   ET4AD.hpp
 * Author: matthewsupernaw
 *
 * Created on May 11, 2014, 8:24 AM
 */


#ifndef AD_ET4AD_HPP
#define	AD_ET4AD_HPP

#include <stdint.h>
#include <cmath>
#include <vector>
#include <set>
#include <iostream>

#define USE_MM_CACHE_MAP
//#define USE_HASH_TABLE
//this is the most crucial part of the software....we need a fast map!
#ifdef USE_MM_CACHE_MAP
#include "support/cache-table-0.2/mm/cache_map.hpp"
#elif defined(USE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(USE_HASH_TABLE)
#include "support/hash_table/hashtable.h"
#else
#include <map>
#endif

#define REAL_T double

namespace ad {

    enum Operation {
        MINUS = 0,
        PLUS,
        MULTIPLY,
        DIVIDE,
        SIN,
        COS,
        TAN,
        ASIN,
        ACOS,
        ATAN,
        ATAN2, //atan(adnumber,adnumber)
        ATAN3, //atan(T,adnumber)
        ATAN4, //atan(adnumber,T)
        SQRT,
        POW, //pow(adnumber,adnumber)
        POW1, //pow(T,adnumber)
        POW2, //pow(adnumber,T)
        LOG,
        LOG10,
        EXP,
        SINH,
        COSH,
        TANH,
        ABS,
        FABS,
        FLOOR,
        CEIL,
        CONSTANT,
        VARIABLE,
        NONE
    };

    /**
     * Statement class is used store the expression in a post-order vector. 
     * Ultimately used for higher order derivatives and expression string 
     * building.
     * 
     * @param op
     */
    class Statement {
    public:

        Statement(const Operation &op) : op_m(op) {

        }

        Statement(const Operation &op, const REAL_T &value) : op_m(op), value_m(value) {

        }

        Statement(const Operation &op, const REAL_T &value, const uint32_t &id) : op_m(op), value_m(value), id_m(id) {

        }

        Operation op_m;
        REAL_T value_m;
        uint32_t id_m;

    };

    /*!
     * Creates a unique identifier.
     * @return 
     */
    class IDGenerator {
    public:
        static IDGenerator * instance();

        const uint32_t next() {
            return ++_id;
        }
    private:

        IDGenerator() : _id(0) {
        }

        uint32_t _id;
    };

    static IDGenerator* only_copy;

    inline IDGenerator *
    IDGenerator::instance() {

        if (!only_copy) {
            only_copy = new IDGenerator();
        }

        return only_copy;
    }

    class Variable;

    template<class A>
    struct ExpressionBase {
    public:

        ExpressionBase() : id_m(0) {

        }

        ExpressionBase(const uint32_t & id) : id_m(id) {

        }

        const A & Cast() const {
            return static_cast<const A&> (*this);
        }

        const REAL_T GetValue() const {
            return Cast().GetValue();

        }

        uint32_t GetId() const {
            return this->id_m;
        }

        /**
         * Compute or get the stored derivative with respect to the unique 
         * identifier id.
         * 
         * @param id
         * @param found
         * @return 
         */
        virtual const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0;
        };






        uint32_t id_m;

    private:

        ExpressionBase& operator=(const ExpressionBase & exp) const {

            return *this;
        };
    };

    struct Literal : public ExpressionBase<Literal> {

        Literal(const REAL_T & value) : value_m(value) {

        }

        const REAL_T GetValue() const {
            return this->value_m;

        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0;
        };

        const REAL_T &value_m;


    };

    template <class LHS, class RHS>
    struct Add : public ExpressionBase<Add<LHS, RHS> > {

        Add(const ExpressionBase<LHS>& lhs, const ExpressionBase<RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() + rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found) + rhs_m.Derivative(id, found);
        }



    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class LHS, class RHS>
    inline
    Add<LHS, RHS> operator+(const ExpressionBase<LHS>& a,
    const ExpressionBase<RHS>& b) {
        return Add<LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    struct Minus : public ExpressionBase<Minus<LHS, RHS> > {

        Minus(const ExpressionBase<LHS>& lhs, const ExpressionBase<RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() - rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found) - rhs_m.Derivative(id, found);
        }

    protected:


    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class LHS, class RHS>
    inline
    Minus<LHS, RHS> operator-(const ExpressionBase<LHS>& lhs,
    const ExpressionBase<RHS>& rhs) {
        return Minus<LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    template <class LHS, class RHS>
    struct Multiply : public ExpressionBase<Multiply<LHS, RHS> > {

        Multiply(const ExpressionBase<LHS>& lhs, const ExpressionBase<RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            //            std::cout << "* " << a_.Derivative(id, found) << "*" << b_ << " +" <<a_ << " * " <<  b_.Derivative(id, found) << " \n";
            return lhs_m.Derivative(id, found) * rhs_m.GetValue() + lhs_m.GetValue() * rhs_m.Derivative(id, found);
        }

    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class LHS, class RHS>
    inline
    Multiply<LHS, RHS> operator*(const ExpressionBase<LHS>& lhs,
    const ExpressionBase<RHS>& rhs) {
        return Multiply<LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    template <class LHS, class RHS>
    struct Divide : public ExpressionBase<Divide<LHS, RHS> > {

        Divide(const ExpressionBase<LHS>& lhs, const ExpressionBase<RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() / rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return (lhs_m.Derivative(id, found) * rhs_m.GetValue() - lhs_m.GetValue() * rhs_m.Derivative(id, found)) / (rhs_m.GetValue() * rhs_m.GetValue());
        }

    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class LHS, class RHS>
    inline
    Divide<LHS, RHS> operator/(const ExpressionBase<LHS>& lhs,
    const ExpressionBase<RHS>& rhs) {
        return Divide<LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    template <class LHS>
    struct LiteralAdd : public ExpressionBase<LiteralAdd<LHS> > {

        LiteralAdd(const ExpressionBase<LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), value_m(lhs_m.GetValue() + rhs) {
        }

        const REAL_T GetValue() const {
            return value_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found);
        }


    private:
        const LHS& lhs_m;
        REAL_T value_m;
    };

    template <class LHS>
    inline
    LiteralAdd<LHS> operator+(const ExpressionBase<LHS>& lhs,
    const REAL_T& rhs) {
        return LiteralAdd<LHS > (lhs.Cast(), rhs);
    }

    template <class RHS>
    inline
    LiteralAdd<RHS> operator+(const REAL_T& lhs,
    const ExpressionBase<RHS>& rhs) {
        return LiteralAdd<RHS > (rhs.Cast(), lhs);
    }

    template <class LHS>
    struct LiteralMinus : public ExpressionBase<LiteralMinus<LHS> > {

        LiteralMinus(const ExpressionBase<LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), value_m(lhs_m.GetValue() - rhs) {
        }

        const REAL_T GetValue() const {
            return value_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found);
        }


    private:
        const LHS& lhs_m;
        REAL_T value_m;
    };

    template <class RHS>
    struct MinusLiteral : public ExpressionBase<MinusLiteral<RHS> > {

        MinusLiteral(const REAL_T & lhs, const ExpressionBase<RHS>& rhs)
        : rhs_m(rhs.Cast()), value_m(lhs - rhs_m.GetValue()) {
        }

        const REAL_T GetValue() const {
            return value_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return -1.0 * rhs_m.Derivative(id, found);
        }


    private:
        const RHS& rhs_m;
        REAL_T value_m;
    };

    template <class LHS>
    inline
    LiteralMinus<LHS> operator-(const ExpressionBase<LHS>& lhs,
    const REAL_T& rhs) {
        return LiteralMinus<LHS > (lhs.Cast(), rhs);
    }

    template <class RHS>
    inline
    MinusLiteral<RHS> operator-(const REAL_T& lhs,
    const ExpressionBase<RHS>& rhs) {
        return MinusLiteral<RHS > (lhs, rhs.Cast());
    }

    template <class LHS>
    struct LiteralTimes : public ExpressionBase<LiteralTimes<LHS> > {

        LiteralTimes(const ExpressionBase<LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return rhs_m * lhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found) * rhs_m;
        }


    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class LHS>
    struct TimesLiteral : public ExpressionBase<TimesLiteral<LHS> > {

        TimesLiteral(const ExpressionBase<LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return rhs_m * lhs_m.Derivative(id, found);
        }


    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class LHS>
    inline
    TimesLiteral<LHS> operator*(const ExpressionBase<LHS>& lhs,
    const REAL_T& rhs) {
        return TimesLiteral<LHS > (lhs.Cast(), rhs);
    }

    template <class RHS>
    inline
    LiteralTimes<RHS> operator*(const REAL_T& lhs,
    const ExpressionBase<RHS>& rhs) {
        return LiteralTimes<RHS > (rhs.Cast(), lhs);
    }

    template <class RHS>
    struct LiteralDivide : public ExpressionBase<LiteralDivide<RHS> > {

        LiteralDivide(const REAL_T & lhs, const ExpressionBase<RHS>& rhs)
        : rhs_m(rhs.Cast()), lhs_m(lhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m / rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return ( -rhs_m.Derivative(id, found) * lhs_m) / (rhs_m.GetValue() * rhs_m.GetValue());
        }


    private:
        const RHS& rhs_m;
        REAL_T lhs_m;
    };

    template <class LHS>
    struct DivideLiteral : public ExpressionBase<DivideLiteral<LHS> > {

        DivideLiteral(const ExpressionBase<LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() / rhs_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return (lhs_m.Derivative(id, found) * rhs_m) / (rhs_m * rhs_m);
        }


    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class LHS>
    inline
    DivideLiteral<LHS> operator/(const ExpressionBase<LHS>& lhs,
    const REAL_T& rhs) {
        return DivideLiteral<LHS > (lhs.Cast(), rhs);
    }

    template <class RHS>
    inline
    LiteralDivide<RHS> operator/(const REAL_T& lhs,
    const ExpressionBase<RHS>& rhs) {
        return LiteralDivide<RHS > (lhs, rhs.Cast());
    }

    template <class EXPR>
    struct Sin : public ExpressionBase<Sin<EXPR> > {

        Sin(const ExpressionBase<EXPR>& a)
        : expr_m(a.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::sin(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * std::cos(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Cos : public ExpressionBase<Cos<EXPR> > {

        Cos(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::cos(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx*-1.0 * std::sin(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Tan : public ExpressionBase<Tan<EXPR> > {

        Tan(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::tan(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return dx * ((1.0 / std::cos(v))*(1.0 / std::cos(v)));
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct ASin : public ExpressionBase<ASin<EXPR> > {

        ASin(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::asin(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * 1.0 / std::pow((1.0 - std::pow(v, 2.0)), 0.5));
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct ACos : public ExpressionBase<ACos<EXPR> > {

        ACos(const ExpressionBase<EXPR>& a)
        : expr_m(a.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::acos(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx* -1.0 / std::pow((1.0 - std::pow(v, 2.0)), 0.5));
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct ATan : public ExpressionBase<ATan<EXPR> > {

        ATan(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::atan(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * 1.0 / (v * v + 1.0));
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR1, class EXPR2>
    struct ATan2 : public ExpressionBase<ATan2<EXPR1, EXPR2> > {

        ATan2(const ExpressionBase<EXPR1>& expr1, const ExpressionBase<EXPR2>& expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::atan2(expr1_m.GetValue(), expr2_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            REAL_T dx = expr1_m.Derivative(id, found);

            if (found) {
                REAL_T v = expr1_m.GetValue();
                REAL_T v2 = expr2_m.GetValue();

                return ((v2 * dx) / (v * v + (v2 * v2)));
            } else {
                return 0.0;
            }
        }

    private:

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    template <class EXPR1>
    struct ATan2Literal : public ExpressionBase<ATan2Literal<EXPR1> > {

        ATan2Literal(const ExpressionBase<EXPR1>& expr1, const REAL_T & expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2) {
        }

        const REAL_T GetValue() const {
            return std::atan2(expr1_m.GetValue(), expr2_m);
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            REAL_T dx = expr1_m.Derivative(id, found);

            if (found) {
                REAL_T v = expr1_m.GetValue();
                REAL_T v2 = expr2_m;

                return ((v2 * dx) / (v * v + (v2 * v2)));
            } else {
                return 0.0;
            }
        }

    private:

        const EXPR1& expr1_m;
        const REAL_T& expr2_m;
    };

    template <class EXPR2>
    struct LiteralATan2 : public ExpressionBase<LiteralATan2<EXPR2> > {

        LiteralATan2(const REAL_T& expr1, const ExpressionBase<EXPR2>& expr2)
        : expr1_m(expr1), expr2_m(expr2.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::atan2(expr1_m, expr2_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            REAL_T dx = 0;

            if (found) {
                REAL_T v = expr1_m;
                REAL_T v2 = expr2_m.GetValue();

                return ((v2 * dx) / (v * v + (v2 * v2)));
            } else {
                return 0.0;
            }
        }

    private:

        const REAL_T& expr1_m;
        const ExpressionBase<EXPR2>& expr2_m;
    };

    template <class EXPR>
    struct Sqrt : public ExpressionBase<Sqrt<EXPR> > {

        Sqrt(const ExpressionBase<EXPR>& a)
        : expr_m(a.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::sqrt(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * .5 / std::sqrt(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR1, class EXPR2>
    struct Pow : public ExpressionBase<Pow<EXPR1, EXPR2> > {

        Pow(const ExpressionBase<EXPR1>& expr1, const ExpressionBase<EXPR2>& expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::pow(expr1_m.GetValue(), expr2_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            REAL_T dx = expr1_m.Derivative(id, found);
            REAL_T dx2 = expr2_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr1_m.GetValue();
                REAL_T v2 = expr2_m.GetValue();

                return (dx * v2) *std::pow(v, (v2 - 1.0));
            } else {
                return 0.0;
            }
        }

    private:

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    template <class EXPR>
    struct PowLiteral : public ExpressionBase<PowLiteral<EXPR> > {

        PowLiteral(const ExpressionBase<EXPR>& expr, const REAL_T & literal)
        : expr_m(expr.Cast()), literal_m(literal) {
        }

        const REAL_T GetValue() const {
            return std::pow(expr_m.GetValue(), literal_m);
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                REAL_T v2 = literal_m;

                return (dx * v2) *std::pow(v, (v2 - 1.0));
            } else {
                return 0.0;
            }
        }

    private:

        const EXPR& expr_m;
        const REAL_T& literal_m;
    };

    template <class EXPR>
    struct LiteralPow : public ExpressionBase<LiteralPow<EXPR> > {

        LiteralPow(const REAL_T& literal, const ExpressionBase<EXPR>& expr)
        : literal_m(literal), expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::pow(literal_m, expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            //          REAL_t dx = 0.0;
            //            if (found) {
            //                REAL_t v = a_;
            //                REAL_t v2 = b_.GetValue();
            //
            //                return (dx * v2) *std::pow(v, (v2 - 1.0));
            //            } else {
            //                return 0.0;
            //            }
            return 0.0;
        }

    private:

        const REAL_T& literal_m;
        const ExpressionBase<EXPR>& expr_m;
    };

    template <class EXPR>
    struct Log : public ExpressionBase<Log<EXPR> > {

        Log(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::log(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);

            if (found) {
                return (dx * 1.0) / expr_m.GetValue();
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Log10 : public ExpressionBase<Log10<EXPR> > {

        Log10(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::log10(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return (dx * 1.0) / (expr_m.GetValue() * std::log(10.0));
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Exp : public ExpressionBase<Exp<EXPR> > {

        Exp(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::exp(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            bool f = false;
            REAL_T dx = expr_m.Derivative(id, f);
            if (f) {
                found = true;
                return dx * std::exp(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Sinh : public ExpressionBase<Sinh<EXPR> > {

        Sinh(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::sinh(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * std::cosh(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Cosh : public ExpressionBase<Cosh<EXPR> > {

        Cosh(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::cosh(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * std::sinh(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Tanh : public ExpressionBase<Tanh<EXPR> > {

        Tanh(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::tanh(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return dx * (1.0 / std::cosh(v))*(1.0 / std::cosh(v));
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Fabs : public ExpressionBase<Fabs<EXPR> > {

        Fabs(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::fabs(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * v) / std::fabs(v);
            } else {
                return 0.0;
            }
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Floor : public ExpressionBase<Floor<EXPR> > {

        Floor(const ExpressionBase<EXPR>& a)
        : expr_m(a.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::floor(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0.0;
        }


    private:
        const EXPR& expr_m;
    };

    template <class EXPR>
    struct Ceil : public ExpressionBase<Ceil<EXPR> > {

        Ceil(const ExpressionBase<EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::ceil(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0.0;
        }


    private:
        const EXPR& expr_m;
    };

    class Variable : public ExpressionBase<Variable> {

        struct eqstr {

            bool operator()(uint32_t s1, uint32_t s2) const {
                return s1>s2;
            }
        };

        REAL_T value_m;
        std::string name_m;
        bool bounded_m;
        REAL_T min_boundary_m;
        REAL_T max_boundary_m;
        bool is_independent_m;
        static std::set<uint32_t> independent_variables_g;
        static bool is_recording_g;


        //Gradient related containers.
        typedef std::set<uint32_t>::iterator ind_iterator;

#ifdef USE_MM_CACHE_MAP
        typedef mm::cache_map<uint32_t, double> GradientMap;
        typedef mm::cache_map<uint32_t, mm::cache_map<uint32_t, double> > MixedPartialsMap;
        typedef GradientMap::const_iterator const_grads_iterator;
        typedef GradientMap::iterator grads_iterator;
        typedef MixedPartialsMap::const_iterator const_patrials_iterator;
        typedef MixedPartialsMap::iterator partials_iterator;
#elif defined(USE_TR1_UNORDERED_MAP)
        typedef std::tr1::unordered_map<uint32_t, double> GradientMap;
        typedef std::tr1::unordered_map<uint32_t, std::tr1::unordered_map<uint32_t, double> > MixedPartialsMap;
        typedef GradientMap::const_iterator const_grads_iterator;
        typedef GradientMap::iterator grads_iterator;
        typedef MixedPartialsMap::const_iterator const_patrials_iterator;
        typedef MixedPartialsMap::iterator partials_iterator;
#elif defined(USE_HASH_TABLE)
        typedef HashTable GradientMap;
#else
        typedef std::map<uint32_t, double> GradientMap;
        typedef std::map<uint32_t, std::map<uint32_t, double> > MixedPartialsMap;
        typedef GradientMap::const_iterator const_grads_iterator;
        typedef GradientMap::iterator grads_iterator;
        typedef MixedPartialsMap::const_iterator const_patrials_iterator;
        typedef MixedPartialsMap::iterator partials_iterator;

#endif
        GradientMap gradients;




    public:

        /**
         * Default constructor.
         */
        Variable() : ExpressionBase(0), value_m(0.0), bounded_m(false) {
        }

        /**
         * Construct a variable with a initial value.optional parameter to make
         * it a independent variable.
         * 
         * @param value
         * @param is_independent
         */
        Variable(const REAL_T& value, bool is_independent = false) : ExpressionBase(0), value_m(value) {
            if (is_independent) {
                this->id_m = IDGenerator::instance()->next();
                this->SetAsIndependent(is_independent);
                this->bounded_m = false;
                this->min_boundary_m = std::numeric_limits<REAL_T>::min();
                this->max_boundary_m = std::numeric_limits<REAL_T>::max();
                Variable::independent_variables_g.insert(this->GetId());
            }

        }

        /**
         * Copy constructor. 
         * 
         * @param rhs
         */
        Variable(const Variable& orig) {
            value_m = orig.GetValue();
            this->id_m = orig.GetId();
            this->gradients.insert(orig.gradients.begin(), orig.gradients.end());
            this->is_independent_m = orig.is_independent_m;
            this->bounded_m = orig.bounded_m;
            this->min_boundary_m = orig.min_boundary_m;
            this->max_boundary_m = orig.max_boundary_m;

        }

        /**
         * Constructs a variable from expression expr.
         * @param rhs
         */
        template<class T>
        Variable(const ExpressionBase<T>& expr) {

            this->bounded_m = false;
            this->min_boundary_m = std::numeric_limits<REAL_T>::min();
            this->max_boundary_m = std::numeric_limits<REAL_T>::max();
            if (Variable::is_recording_g) {
                this->id_m = expr.GetId();
                ind_iterator it;
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = expr.Derivative(*it, found);
                }
            }
            value_m = expr.GetValue();

        }

        /**
         * Control function. If true, derivatives are computed. If false
         * expressions are only evaluated.
         * 
         * @param is_recording
         */
        static void SetRecording(const bool &is_recording) {
            Variable::is_recording_g = is_recording;
        }

        /**
         * If true, derivatives are calculated and stored. If false,
         * expressions are evaluated for value only.
         * 
         * @return 
         */
        static bool IsRecording() {
            return Variable::is_recording_g;
        }

        /**
         * Returns the value of this variable.
         * 
         * @return 
         */
        const REAL_T GetValue() const {
            return this->value_m;
        }

        /**
         * Sets the value of this variable. If the variable is bounded, 
         * the value will be set between the min and max boundary. If the
         * value is less than the minimum boundary, the variables value is 
         * set to minimum boundary. If the  value is greater than the maximum
         * boundary, the variables value is set to maximum boundary. If the 
         * value is signaling nan, the value is set to the mid point between
         * the min and max boundary. 
         * 
         * @param value
         */
        void SetValue(const REAL_T &value) {
            if (this->bounded_m) {
                if (value != value) {//nan
                    this->value_m = this->min_boundary_m + (this->max_boundary_m - this->min_boundary_m) / 2.0;

                    return;
                }

                if (value<this->min_boundary_m) {
                    this->value_m = this->min_boundary_m;
                } else if (value>this->max_boundary_m) {
                    this->value_m = this->max_boundary_m;


                } else {
                    this->value_m = value;
                }
            } else {
                this->value_m = value;
            }
        }

        /**
         * Returns the name of this variable. Names are not initialized and 
         * if it is not set this function will return a empty string.
         * @return 
         */
        std::string GetName() const {
            return name_m;
        }

        /**
         * Set the name of this variable.
         * 
         * @param name
         */
        void SetName(std::string name) {
            this->name_m = name;
        }

        /**
         * Returns true if this variable is bounded. Otherwise,
         * false.
         * 
         * @return 
         */
        bool IsBounded() {
            return this->bounded_m;
        }

        /**
         * Set the min and max boundary for this variable. 
         * @param min
         * @param max
         */
        void SetBounds(const REAL_T& min, const REAL_T& max) {
            this->bounded_m = true;
            this->max_boundary_m = max;
            this->min_boundary_m = min;
            this->SetValue(this->GetValue());
        }

        /**
         * Return the minimum boundary for this variable. The default value 
         * is the results of std::numeric_limits<REAL_t>::min().
         * 
         * @return 
         */
        const REAL_T GetMinBoundary() {
            return this->min_boundary_m;
        }

        /**
         * Return the maximum boundary for this variable. The default value 
         * is the results of std::numeric_limits<REAL_t>::max().
         * 
         * @return 
         */
        const REAL_T GetMaxBoundary() {
            return this->max_boundary_m;
        }

        /**
         * Make this Variable an independent variable. If set to true,
         * the Variables unique identifier is registered in the static set
         * of independent variables. During function evaluations, gradient 
         * information is accumulated in post-order wrt to variables in the 
         * static set of independent variables. If set false, this Variables 
         * unique identifier will be removed from the set, if it existed.
         * 
         * 
         * To access the derivative wrt to a Variable, call function:
         * 
         * const REAL_t wrt(const Variable & ind)
         * 
         * 
         * @param is_independent
         */
        void SetAsIndependent(const bool &is_independent) {
            if (this->GetId() == 0) {
                this->id_m = IDGenerator::instance()->next();
            }
            if (this->is_independent_m && !is_independent) {
                Variable::independent_variables_g.erase(this->GetId());
                this->is_independent_m = false;
            }

            if (!this->is_independent_m && is_independent) {
                Variable::independent_variables_g.insert(this->GetId());
                this->is_independent_m = true;
            }
        }

        /**
         * Returns the derivative with respect to a variables who's
         * unique identifier is equal to the parameter id.
         * 
         * @param id
         * @param found
         * @return 
         */
        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {


            if (this->GetId() == id) {
                found = true;
                return 1.0;

            } else {
                const_grads_iterator git = this->gradients.find(id);

                if (git != this->gradients.end()) {
                    found = true;
                    return git->second;
                } else {
                    return 0.0;
                }
            }

        }

        /**
         * Finds the derivative in the encapsulated gradient map 
         * w.r.t a variable. 
         * 
         * @param ind
         * @return 
         */
        const REAL_T WRT(const Variable & ind) {
            const_grads_iterator git = this->gradients.find(ind.GetId());
            if (git != this->gradients.end()) {
                return git->second;
            } else {
                return 0;
            }
        }

        /**
         * Finds the derivative in the encapsulated gradient map 
         * w.r.t a variable. 
         * @param ind
         * @return 
         */
        const REAL_T WRT(const Variable & ind) const {
            const_grads_iterator git = this->gradients.find(ind.GetId());
            if (git != this->gradients.end()) {
                return git->second;
            } else {
                return 0;
            }
        }

        /**
         * Set the variable equal to the real type rhs.
         * 
         * The gradient map will be cleared.
         * 
         * @param rhs
         * @return 
         */
        Variable& operator=(const REAL_T& rhs) {
            this->gradients.clear();
            this->SetValue(rhs);
            return *this;
        }

        /**
         * Sets the value of this Variable to that of rhs. Derivatives 
         * are stored in the encapsulated gradient map.
         * 
         * @param rhs
         * @return 
         */
        Variable& operator=(const Variable& rhs) {
            //            this->id = rhs.getId();
            value_m = rhs.GetValue();
            //            this->gradients.clear();
            if (ad::Variable::is_recording_g) {
                //                this->gradients.insert(rhs.gradients.begin(), rhs.gradients.end());
                ind_iterator it;
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
            }
            this->is_independent_m = rhs.is_independent_m;
            return *this;
        }

        /**
         * Set the Variables value to the result of the expression rhs. 
         * Derivatives are calculated and stored in the encapsulated 
         * gradient map.
         * @param rhs
         * @return 
         */
        template<class T>
        Variable& operator=(const ExpressionBase<T>& rhs) {
            this->id_m = rhs.GetId();
            if (Variable::is_recording_g) {
                ind_iterator it; // = this->gradients.lower_bound();

                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
            }
            value_m = rhs.GetValue();
            return *this;
        }

        /**
         * 
         * @param rhs
         * @return 
         */
        template<class T>
        Variable& operator+=(const ExpressionBase<T>& rhs) {
            return *this = (*this +rhs);
        }

        Variable& operator+=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = this->gradients[(*it)] + rhs.Derivative(*it, found);
                }
            }
            this->value_m += rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable& operator-=(const ExpressionBase<T>& rhs) {
            return *this = (*this -rhs);
        }

        Variable& operator-=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;

                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = (this->gradients[(*it)] - rhs.Derivative(*it, found));
                }
            }
            this->value_m -= rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable& operator*=(const ExpressionBase<T>& rhs) {
            return *this = (*this * rhs);
        }

        Variable& operator*=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;

                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = this->gradients[(*it)] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
            }
            this->value_m *= rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable& operator/=(const ExpressionBase<T>& rhs) {
            return *this = (*this / rhs);
        }

        Variable& operator/=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;

                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients[(*it)] = (this->gradients[(*it)] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
            }
            this->value_m /= rhs.GetValue();
            return *this;
        }

        // And likewise for a Literal on the rhs

        Variable& operator+=(const REAL_T& rhs) {
            value_m += rhs;
            return *this;
        }

        Variable& operator-=(const REAL_T& rhs) {
            value_m -= rhs;
            return *this;
        }

        Variable& operator*=(const REAL_T& rhs) {
            return *this = (*this * rhs);
        }

        Variable& operator/=(const REAL_T& rhs) {
            return *this = (*this / rhs);
        }


    };

    //list of independent variables by unique identifier.
    std::set<uint32_t> Variable::independent_variables_g;
    bool Variable::is_recording_g = true;

    template<class T, class TT>
    inline const int operator==(const ad::ExpressionBase<T>& lhs, const ad::ExpressionBase<TT>& rhs) {
        return lhs.GetValue() == rhs.GetValue();
    }

    template<class T, class TT>
    inline const int operator!=(const ad::ExpressionBase<T>& lhs, const ad::ExpressionBase<TT>& rhs) {
        return lhs.GetValue() != rhs.GetValue();
    }

    template<class T, class TT>
    inline const int operator<(const ad::ExpressionBase<T>& lhs, const ad::ExpressionBase<TT>& rhs) {
        return lhs.GetValue() < rhs.GetValue();
    }

    template<class T, class TT>
    inline const int operator>(const ad::ExpressionBase<T>& lhs, const ad::ExpressionBase<TT>& rhs) {
        return lhs.GetValue() > rhs.GetValue();
    }

    template<class T, class TT>
    inline const int operator<=(const ad::ExpressionBase<T>& lhs, const ad::ExpressionBase<TT>& rhs) {
        return lhs.GetValue() <= rhs.GetValue();
    }

    template<class T, class TT>
    inline const int operator>=(const ad::ExpressionBase<T>& lhs, const ad::ExpressionBase<TT>& rhs) {
        return lhs.GetValue() >= rhs.GetValue();
    }

    template<class T>
    inline const int operator==(const T &lhs, const ad::ExpressionBase<T>& rhs) {
        return lhs == rhs.GetValue();
    }

    template<class T>
    inline const int operator!=(const T &lhs, const ad::ExpressionBase<T>& rhs) {
        return lhs != rhs.GetValue();
    }

    template<class T>
    inline const int operator<(const T &lhs, const ad::ExpressionBase<T>& rhs) {
        return lhs < rhs.GetValue();
    }

    template<class T>
    inline const int operator>(const T &lhs, const ad::ExpressionBase<T>& rhs) {
        return lhs > rhs.GetValue();
    }

    template<class T>
    inline const int operator<=(const T &lhs, const ad::ExpressionBase<T>& rhs) {
        return lhs <= rhs.GetValue();
    }

    template<class T>
    inline const int operator>=(const T &lhs, const ad::ExpressionBase<T>& rhs) {
        return lhs >= rhs.GetValue();
    }

    template<class T>
    inline const int operator==(const ad::ExpressionBase<T>& lhs, const T &rhs) {
        return lhs.GetValue() == rhs;
    }

    template<class T>
    inline const int operator!=(const ad::ExpressionBase<T>& lhs, const T &rhs) {
        return lhs.GetValue() != rhs;
    }

    template<class T>
    inline const int operator<(const ad::ExpressionBase<T>& lhs, const T &rhs) {
        return lhs.GetValue() <= rhs;
    }

    template<class T>
    inline const int operator>(const ad::ExpressionBase<T>& lhs, const T &rhs) {
        return lhs.GetValue() > rhs;
    }

    template<class T>
    inline const int operator<=(const ad::ExpressionBase<T>& lhs, const T &rhs) {
        return lhs.GetValue() <= rhs;
    }

    template<class T>
    inline const int operator>=(const ad::ExpressionBase<T>& lhs, const T &rhs) {
        return lhs.GetValue() >= rhs;
    }


}




namespace std {

    /**
     * Write the expression value to std::ostream out.
     * @param out
     * @param exp
     * @return 
     */
    template<class A>
    inline std::ostream& operator<<(std::ostream &out, const ad::ExpressionBase<A> &exp) {
        out << exp.GetValue();
        return out;
    }

    /**
     * Override for the sin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Sin<EXPR> sin(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Sin<EXPR > (expr.Cast());
    }

    /**
     * Override for the cos function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Cos<EXPR> cos(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Cos<EXPR > (expr.Cast());
    }

    /**
     * Override for the tan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Tan<EXPR> tan(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Tan<EXPR > (expr.Cast());
    }

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::ASin<EXPR> asin(const ad::ExpressionBase<EXPR>& expr) {
        return ad::ASin<EXPR > (expr.Cast());
    }

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::ACos<EXPR> acos(const ad::ExpressionBase<EXPR>& expr) {
        return ad::ACos<EXPR > (expr.Cast());
    }

    /**
     * Override for the atan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::ATan<EXPR> atan(const ad::ExpressionBase<EXPR>& a) {
        return ad::ATan<EXPR > (a.Cast());
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class EXPR1, class EXPR2>
    inline
    ad::ATan2<EXPR1, EXPR2> atan2(const ad::ExpressionBase<EXPR1>& expr1,
    const ad::ExpressionBase<EXPR2>& expr2) {
        return ad::ATan2<EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class EXPR>
    inline
    ad::ATan2Literal<EXPR> atan2(const ad::ExpressionBase<EXPR>& expr,
    const REAL_T& val) {
        return ad::ATan2Literal<EXPR > (expr.Cast(), val);
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class EXPR>
    inline
    ad::LiteralATan2<EXPR> atan2(const REAL_T& val,
    const ad::ExpressionBase<EXPR>& expr) {
        return ad::LiteralATan2<EXPR > (val, expr.Cast());
    }

    /**
     * Override for the sqrt function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Sqrt<EXPR> sqrt(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Sqrt<EXPR > (expr.Cast());
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class EXPR1, class EXPR2>
    inline
    ad::Pow<EXPR1, EXPR2> pow(const ad::ExpressionBase<EXPR1>& expr1,
    const ad::ExpressionBase<EXPR2>& expr2) {
        return ad::Pow<EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class EXPR>
    inline
    ad::PowLiteral<EXPR> pow(const ad::ExpressionBase<EXPR>& expr,
    const REAL_T& val) {
        return ad::PowLiteral<EXPR > (expr.Cast(), val);
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param val
     * @param expr2
     * @return 
     */
    template <class EXPR>
    inline
    ad::LiteralPow<EXPR> pow(const REAL_T& val,
    const ad::ExpressionBase<EXPR>& expr) {
        return ad::LiteralPow<EXPR > (val, expr.Cast());
    }

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Log<EXPR> log(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Log<EXPR > (expr.Cast());
    }

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Log10<EXPR> log10(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Log10<EXPR > (expr.Cast());
    }

    /**
     * Override for the exp function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Exp<EXPR> exp(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Exp<EXPR > (expr.Cast());
    }

    /**
     * Override for the sinh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Sinh<EXPR> sinh(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Sinh<EXPR > (expr.Cast());
    }

    /**
     * Override for the cosh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Cosh<EXPR> cosh(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Cosh<EXPR > (expr.Cast());
    }

    /**
     * Override for the tanh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Tanh<EXPR> tanh(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Tanh<EXPR > (expr.Cast());
    }

    /**
     * Override for the fabs function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Fabs<EXPR> fabs(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Fabs<EXPR > (expr.Cast());
    }

    /**
     * Override for the floor function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Floor<EXPR> floor(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Floor<EXPR > (expr.Cast());
    }

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class EXPR>
    inline ad::Ceil<EXPR> ceil(const ad::ExpressionBase<EXPR>& expr) {
        return ad::Ceil<EXPR > (expr.Cast());
    }



}

#endif	/* ETAD_HPP */

