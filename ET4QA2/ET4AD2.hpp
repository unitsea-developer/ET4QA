
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
#include <stack>

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

//#define REAL_T double

namespace ad2 {

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
    template<class REAL_T>
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

    template<class REAL_T>
    class Variable;

    template<class REAL_T, class A>
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

        virtual void Push(std::vector<Statement<REAL_T> > &statements) const {

        }

        virtual void MakeSet(std::set<uint32_t> &independent_variables)const {

        }





        uint32_t id_m;

    private:

        ExpressionBase& operator=(const ExpressionBase & exp) const {

            return *this;
        };
    };

    template<class REAL_T>
    struct Literal : public ExpressionBase<REAL_T, Literal<REAL_T> > {

        Literal(const REAL_T & value) : value_m(value) {

        }

        const REAL_T GetValue() const {
            return this->value_m;

        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0;
        };

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, value_m));
        }

        const REAL_T &value_m;


    };

    template <class REAL_T, class LHS, class RHS>
    struct Add : public ExpressionBase<REAL_T, Add<REAL_T, LHS, RHS> > {

        Add(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() + rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found) + rhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (PLUS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
            this->rhs_m.MakeSet(independent_variables);
        }


    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class REAL_T, class LHS, class RHS>
    inline
    Add<REAL_T, LHS, RHS> operator+(const ExpressionBase<REAL_T, LHS>& a,
    const ExpressionBase<REAL_T, RHS>& b) {
        return Add<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class REAL_T, class LHS, class RHS>
    struct Minus : public ExpressionBase<REAL_T, Minus<REAL_T, LHS, RHS> > {

        Minus(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() - rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found) - rhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MINUS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
            this->rhs_m.MakeSet(independent_variables);
        }

    protected:


    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class REAL_T, class LHS, class RHS>
    inline
    Minus<REAL_T, LHS, RHS> operator-(const ExpressionBase<REAL_T, LHS>& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return Minus<REAL_T, LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    template <class REAL_T, class LHS, class RHS>
    struct Multiply : public ExpressionBase<REAL_T, Multiply<REAL_T, LHS, RHS> > {

        Multiply(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            //            std::cout << "* " << a_.Derivative(id, found) << "*" << b_ << " +" <<a_ << " * " <<  b_.Derivative(id, found) << " \n";
            return lhs_m.Derivative(id, found) * rhs_m.GetValue() + lhs_m.GetValue() * rhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MULTIPLY));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
            this->rhs_m.MakeSet(independent_variables);
        }

    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class REAL_T, class LHS, class RHS>
    inline
    Multiply<REAL_T, LHS, RHS> operator*(const ExpressionBase<REAL_T, LHS>& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return Multiply<REAL_T, LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    template <class REAL_T, class LHS, class RHS>
    struct Divide : public ExpressionBase<REAL_T, Divide<REAL_T, LHS, RHS> > {

        Divide(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() / rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return (lhs_m.Derivative(id, found) * rhs_m.GetValue() - lhs_m.GetValue() * rhs_m.Derivative(id, found)) / (rhs_m.GetValue() * rhs_m.GetValue());
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (DIVIDE));
        }

    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
    };

    template <class REAL_T, class LHS, class RHS>
    inline
    Divide<REAL_T, LHS, RHS> operator/(const ExpressionBase<REAL_T, LHS>& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return Divide<REAL_T, LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    template <class REAL_T, class LHS>
    struct LiteralAdd : public ExpressionBase<REAL_T, LiteralAdd<REAL_T, LHS> > {

        LiteralAdd(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() + rhs_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (PLUS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
        }



    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class REAL_T, class LHS>
    inline
    LiteralAdd<REAL_T, LHS> operator+(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return LiteralAdd<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    inline
    LiteralAdd<REAL_T, RHS> operator+(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return LiteralAdd<REAL_T, RHS > (rhs.Cast(), lhs);
    }

    template <class REAL_T, class LHS>
    struct LiteralMinus : public ExpressionBase<REAL_T, LiteralMinus<REAL_T, LHS> > {

        LiteralMinus(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() - rhs_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return lhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (MINUS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
        }


    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class REAL_T, class RHS>
    struct MinusLiteral : public ExpressionBase<REAL_T, MinusLiteral<REAL_T, RHS> > {

        MinusLiteral(const REAL_T & lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : rhs_m(rhs.Cast()), lhs_m(lhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m - rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return -1.0 * rhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, lhs_m));
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MINUS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->rhs_m.MakeSet(independent_variables);
        }

    private:
        const RHS& rhs_m;
        REAL_T lhs_m;
    };

    template <class REAL_T, class LHS>
    inline
    LiteralMinus<REAL_T, LHS> operator-(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return LiteralMinus<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    inline
    MinusLiteral<REAL_T, RHS> operator-(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return MinusLiteral<REAL_T, RHS > (lhs, rhs.Cast());
    }

    template <class REAL_T, class LHS>
    struct LiteralTimes : public ExpressionBase<REAL_T, LiteralTimes<REAL_T, LHS> > {

        LiteralTimes(const REAL_T & lhs, const ExpressionBase<REAL_T, LHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()) {
        }

        const REAL_T GetValue() const {
            return lhs_m * rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return rhs_m.Derivative(id, found) * lhs_m;
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, lhs_m));
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MULTIPLY));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->rhs_m.MakeSet(independent_variables);
        }


    private:
        const LHS& rhs_m;
        REAL_T lhs_m;
    };

    template <class REAL_T, class LHS>
    struct TimesLiteral : public ExpressionBase<REAL_T, TimesLiteral<REAL_T, LHS> > {

        TimesLiteral(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return rhs_m * lhs_m.Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (MULTIPLY));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
        }

    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class REAL_T, class LHS>
    inline
    TimesLiteral<REAL_T, LHS> operator*(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return TimesLiteral<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    inline
    LiteralTimes<REAL_T, RHS> operator*(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return LiteralTimes<REAL_T, RHS > (lhs, rhs.Cast());
    }

    template <class REAL_T, class RHS>
    struct LiteralDivide : public ExpressionBase<REAL_T, LiteralDivide<REAL_T, RHS> > {

        LiteralDivide(const REAL_T & lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : rhs_m(rhs.Cast()), lhs_m(lhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m / rhs_m.GetValue();
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return ( -rhs_m.Derivative(id, found) * lhs_m) / (rhs_m.GetValue() * rhs_m.GetValue());
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, lhs_m));
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (DIVIDE));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->rhs_m.MakeSet(independent_variables);
        }


    private:
        const RHS& rhs_m;
        REAL_T lhs_m;
    };

    template <class REAL_T, class LHS>
    struct DivideLiteral : public ExpressionBase<REAL_T, DivideLiteral<REAL_T, LHS> > {

        DivideLiteral(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {
        }

        const REAL_T GetValue() const {
            return lhs_m.GetValue() / rhs_m;
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return (lhs_m.Derivative(id, found) * rhs_m) / (rhs_m * rhs_m);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (DIVIDE));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->lhs_m.MakeSet(independent_variables);
        }


    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
    };

    template <class REAL_T, class LHS>
    inline
    DivideLiteral<REAL_T, LHS> operator/(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return DivideLiteral<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    inline
    LiteralDivide<REAL_T, RHS> operator/(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return LiteralDivide<REAL_T, RHS > (lhs, rhs.Cast());
    }

    template <class REAL_T, class EXPR>
    struct Sin : public ExpressionBase<REAL_T, Sin<REAL_T, EXPR> > {

        Sin(const ExpressionBase<REAL_T, EXPR>& a)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (SIN));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Cos : public ExpressionBase<REAL_T, Cos<REAL_T, EXPR> > {

        Cos(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (COS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Tan : public ExpressionBase<REAL_T, Tan<REAL_T, EXPR> > {

        Tan(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (TAN));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct ASin : public ExpressionBase<REAL_T, ASin<REAL_T, EXPR> > {

        ASin(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::asin(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * static_cast<REAL_T> (1.0) / std::pow((static_cast<REAL_T> (1.0) - std::pow(v, static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
            } else {
                return 0.0;
            }
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ASIN));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct ACos : public ExpressionBase<REAL_T, ACos<REAL_T, EXPR> > {

        ACos(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::acos(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * static_cast<REAL_T> (-1.0) / std::pow((static_cast<REAL_T> (1.0) - std::pow(v, static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
            } else {
                return 0.0;
            }
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ACOS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct ATan : public ExpressionBase<REAL_T, ATan<REAL_T, EXPR> > {

        ATan(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ATAN));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR1, class EXPR2>
    struct ATan2 : public ExpressionBase<REAL_T, ATan2<REAL_T, EXPR1, EXPR2> > {

        ATan2(const ExpressionBase<REAL_T, EXPR1>& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr1_m.Push(statements);
            this->expr2_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ATAN));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr1_m.MakeSet(independent_variables);
            this->expr2_m.MakeSet(independent_variables);
        }


    private:

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    template <class REAL_T, class EXPR1>
    struct ATan2Literal : public ExpressionBase<REAL_T, ATan2Literal<REAL_T, EXPR1> > {

        ATan2Literal(const ExpressionBase<REAL_T, EXPR1>& expr1, const REAL_T & expr2)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr1_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, expr2_m));
            statements.push_back(Statement<REAL_T > (ATAN2));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr1_m.MakeSet(independent_variables);
        }

    private:

        const EXPR1& expr1_m;
        const REAL_T& expr2_m;
    };

    template <class REAL_T, class EXPR2>
    struct LiteralATan2 : public ExpressionBase<REAL_T, LiteralATan2<REAL_T, EXPR2> > {

        LiteralATan2(const REAL_T& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, expr1_m));
            expr2_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ATAN2));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr2_m.MakeSet(independent_variables);
        }
    private:

        const REAL_T& expr1_m;
        const ExpressionBase<REAL_T, EXPR2>& expr2_m;
    };

    template <class REAL_T, class EXPR>
    struct Sqrt : public ExpressionBase<REAL_T, Sqrt<REAL_T, EXPR> > {

        Sqrt(const ExpressionBase<REAL_T, EXPR>& a)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (SQRT));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }
    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR1, class EXPR2>
    struct Pow : public ExpressionBase<REAL_T, Pow<REAL_T, EXPR1, EXPR2> > {

        Pow(const ExpressionBase<REAL_T, EXPR1>& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
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

                return (dx * v2) *std::pow(v, (v2 - static_cast<REAL_T> (1.0)));
            } else {
                return 0.0;
            }
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr1_m.Push(statements);
            this->expr2_m.Push(statements);
            statements.push_back(Statement<REAL_T > (POW));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr1_m.MakeSet(independent_variables);
            this->expr2_m.MakeSet(independent_variables);
        }

    private:

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    template <class REAL_T, class EXPR>
    struct PowLiteral : public ExpressionBase<REAL_T, PowLiteral<REAL_T, EXPR> > {

        PowLiteral(const ExpressionBase<REAL_T, EXPR>& expr, const REAL_T & literal)
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

                return (dx * v2) *std::pow(v, (v2 - static_cast<REAL_T> (1.0)));
            } else {
                return 0.0;
            }
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, literal_m));
            statements.push_back(Statement<REAL_T > (POW));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:

        const EXPR& expr_m;
        const REAL_T& literal_m;
    };

    template <class REAL_T, class EXPR>
    struct LiteralPow : public ExpressionBase<REAL_T, LiteralPow<REAL_T, EXPR> > {

        LiteralPow(const REAL_T& literal, const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, literal_m));
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (POW));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:

        const REAL_T& literal_m;
        const ExpressionBase<REAL_T, EXPR>& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Log : public ExpressionBase<REAL_T, Log<REAL_T, EXPR> > {

        Log(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (LOG));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Log10 : public ExpressionBase<REAL_T, Log10<REAL_T, EXPR> > {

        Log10(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (LOG10));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Exp : public ExpressionBase<REAL_T, Exp<REAL_T, EXPR> > {

        Exp(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (EXP));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Sinh : public ExpressionBase<REAL_T, Sinh<REAL_T, EXPR> > {

        Sinh(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (SINH));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Cosh : public ExpressionBase<REAL_T, Cosh<REAL_T, EXPR> > {

        Cosh(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (COSH));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Tanh : public ExpressionBase<REAL_T, Tanh<REAL_T, EXPR> > {

        Tanh(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (TANH));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Fabs : public ExpressionBase<REAL_T, Fabs<REAL_T, EXPR> > {

        Fabs(const ExpressionBase<REAL_T, EXPR>& expr)
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

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (FABS));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Floor : public ExpressionBase<REAL_T, Floor<REAL_T, EXPR> > {

        Floor(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::floor(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0.0;
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (FLOOR));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }


    private:
        const EXPR& expr_m;
    };

    template <class REAL_T, class EXPR>
    struct Ceil : public ExpressionBase<REAL_T, Ceil<REAL_T, EXPR> > {

        Ceil(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        const REAL_T GetValue() const {
            return std::ceil(expr_m.GetValue());
        }

        const REAL_T Derivative(const uint32_t &id, bool &found, uint32_t index = 0) const {
            return 0.0;
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CEIL));
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            this->expr_m.MakeSet(independent_variables);
        }

    private:
        const EXPR& expr_m;
    };


}




namespace std2 {

    /**
     * Write the expression value to std::ostream out.
     * @param out
     * @param exp
     * @return 
     */
    template<class REAL_T, class A>
    inline std::ostream& operator<<(std::ostream &out, const ad::ExpressionBase<REAL_T, A> &exp) {
        out << exp.GetValue();
        return out;
    }

    /**
     * Override for the sin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Sin<REAL_T, EXPR> sin(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Sin<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the cos function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Cos<REAL_T, EXPR> cos(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Cos<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the tan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Tan<REAL_T, EXPR> tan(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Tan<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::ASin<REAL_T, EXPR> asin(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::ASin<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::ACos<REAL_T, EXPR> acos(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::ACos<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the atan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::ATan<REAL_T, EXPR> atan(const ad::ExpressionBase<REAL_T, EXPR>& a) {
        return ad::ATan<REAL_T, EXPR > (a.Cast());
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    ad::ATan2<REAL_T, EXPR1, EXPR2> atan2(const ad::ExpressionBase<REAL_T, EXPR1>& expr1,
    const ad::ExpressionBase<REAL_T, EXPR2>& expr2) {
        return ad::ATan2<REAL_T, EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    ad::ATan2Literal<REAL_T, EXPR> atan2(const ad::ExpressionBase<REAL_T, EXPR>& expr,
    const REAL_T& val) {
        return ad::ATan2Literal<REAL_T, EXPR > (expr.Cast(), val);
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    ad::LiteralATan2<REAL_T, EXPR> atan2(const REAL_T& val,
    const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::LiteralATan2<REAL_T, EXPR > (val, expr.Cast());
    }

    /**
     * Override for the sqrt function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Sqrt<REAL_T, EXPR> sqrt(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Sqrt<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    ad::Pow<REAL_T, EXPR1, EXPR2> pow(const ad::ExpressionBase<REAL_T, EXPR1>& expr1,
    const ad::ExpressionBase<REAL_T, EXPR2>& expr2) {
        return ad::Pow<REAL_T, EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    ad::PowLiteral<REAL_T, EXPR> pow(const ad::ExpressionBase<REAL_T, EXPR>& expr,
    const REAL_T& val) {
        return ad::PowLiteral<REAL_T, EXPR > (expr.Cast(), val);
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param val
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    ad::LiteralPow<REAL_T, EXPR> pow(const REAL_T& val,
    const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::LiteralPow<REAL_T, EXPR > (val, expr.Cast());
    }

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Log<REAL_T, EXPR> log(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Log<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Log10<REAL_T, EXPR> log10(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Log10<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the exp function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Exp<REAL_T, EXPR> exp(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Exp<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the sinh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Sinh<REAL_T, EXPR> sinh(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Sinh<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the cosh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Cosh<REAL_T, EXPR> cosh(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Cosh<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the tanh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Tanh<REAL_T, EXPR> tanh(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Tanh<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the fabs function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Fabs<REAL_T, EXPR> fabs(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Fabs<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the floor function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Floor<REAL_T, EXPR> floor(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Floor<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline ad::Ceil<REAL_T, EXPR> ceil(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Ceil<REAL_T, EXPR > (expr.Cast());
    }



}

namespace ad2 {

    /**
     * Variable class used for calculations and storage of computed derivatives.
     * 
     * Derivatives are only computed for variables marked as independent 
     * variables. To make a variable independent use:
     * 
     * Conclassor:
     * 
     * Variable(const REAL_T& value, bool is_independent = false)
     * 
     * or
     * 
     * Member Function:
     * 
     * void SetAsIndependent(const bool &is_independent)
     * 
     * 
     * 
     */
    template<class REAL_T>
    class Variable : public ExpressionBase<REAL_T, Variable<REAL_T> > {

        class eqstr {

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

        static bool is_supporting_arbitrary_order;

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
        GradientMap gradients_m;

        typedef std::vector<Statement<REAL_T> > ExpressionStatements;
        ExpressionStatements statements_m;

    public:

        /**
         * Default conclassor.
         */
        Variable() : ExpressionBase<REAL_T, Variable<REAL_T> >(0), value_m(0.0), bounded_m(false), is_independent_m(false) {
            if (Variable::is_recording_g) {
                if (Variable::IsSupportingArbitraryOrder()) {
                    this->statements_m.push_back(Statement<REAL_T > (VARIABLE, this->GetValue(), this->GetId()));
                }
            }
        }

        /**
         * Conclass a variable with a initial value. Optional parameter to make
         * it a independent variable.
         * 
         * @param value
         * @param is_independent
         */
        Variable(const REAL_T& value, bool is_independent = false) : is_independent_m(is_independent), ExpressionBase<REAL_T, Variable<REAL_T> >(0), value_m(value) {
            if (is_independent) {
                this->id_m = IDGenerator::instance()->next();
                this->SetAsIndependent(is_independent);
                this->bounded_m = false;
                this->min_boundary_m = std::numeric_limits<REAL_T>::min();
                this->max_boundary_m = std::numeric_limits<REAL_T>::max();
                Variable::independent_variables_g.insert(this->GetId());
                if (Variable::IsSupportingArbitraryOrder()) {
                    this->statements_m.push_back(Statement<REAL_T > (VARIABLE, this->GetValue(), this->GetId()));
                }
            }

        }

        /**
         * Copy conclassor. 
         * 
         * @param rhs
         */
        Variable(const Variable& orig) {
            value_m = orig.GetValue();
            this->id_m = orig.GetId();
#if defined(USE_HASH_TABLE)
            this->gradients_m = HashTable(orig.gradients_m);
#else
            this->gradients_m.insert(orig.gradients_m.begin(), orig.gradients_m.end());
#endif
            this->is_independent_m = orig.is_independent_m;
            this->bounded_m = orig.bounded_m;
            this->min_boundary_m = orig.min_boundary_m;
            this->max_boundary_m = orig.max_boundary_m;
            if (Variable::IsSupportingArbitraryOrder()) {
                orig.Push(this->statements_m);
            }

        }

        /**
         * Conclasss a variable from expression expr.
         * @param rhs
         */
        template<class T>
        Variable(const ExpressionBase<REAL_T, T>& expr) {

            this->bounded_m = false;
            this->min_boundary_m = std::numeric_limits<REAL_T>::min();
            this->max_boundary_m = std::numeric_limits<REAL_T>::max();
            if (Variable::is_recording_g) {
                this->id_m = expr.GetId();
                ind_iterator it;
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
#if defined(USE_HASH_TABLE)
                    this->gradients_m.Insert(*it)->value = expr.Derivative(*it, found);
#else
                    this->gradients_m[(*it)] = expr.Derivative(*it, found);
#endif
                }

                if (Variable::IsSupportingArbitraryOrder()) {
                    expr.Push(this->statements_m);
                }
            }
            value_m = expr.GetValue();

        }

        ~Variable() {
#if !defined(USE_HASH_TABLE)

            this->gradients_m.clear();

#endif
        }

        size_t Size() {
            return this->statements_m.size();
        }

        const REAL_T Diff(const Variable &wrt) {

            if (this->statements_m.size() == 0) {
                return 0.0;
            }
            std::vector<std::pair<double, double> > v;
            v.reserve(this->statements_m.size());
            std::stack<std::pair<double, double>,
                    std::vector<std::pair<double, double> > > stack;
            //            ad::Stack<std::pair<T, T> > stack;
            bool found = false;
            int size = this->statements_m.size();
            std::pair<REAL_T, REAL_T> lhs = std::pair<REAL_T, REAL_T > (0, 0);
            std::pair<REAL_T, REAL_T> rhs = std::pair<REAL_T, REAL_T > (0, 0);

            Statement<REAL_T>* edata = (Statement<REAL_T>*)this->statements_m.data();
            std::cout << this->statements_m.size() << "\n\n";
            for (int i = 0; i < size; i++) {


                REAL_T temp = 0;



                switch (edata[i].op_m) {

                    case CONSTANT:
                        stack.push(std::pair<REAL_T, REAL_T > (edata[i].value_m, 0.0));
                        break;
                    case VARIABLE:
                        if (edata[i].id_m == wrt.GetId() && wrt.GetId() > 0) {
                            found = true;
                            //f(x) = x
                            //f'(x) = 1
                            stack.push(std::pair<REAL_T, REAL_T > (edata[i].value_m, 1.0));
                        } else {//constant
                            //f(x) = C
                            //f'(x) = 0
                            stack.push(std::pair<REAL_T, REAL_T > (edata[i].value_m, 0.0));
                        }
                        break;
                    case PLUS:
                        rhs = stack.top();
                        stack.pop();
                        lhs = stack.top();
                        stack.pop();
                        stack.push(std::pair<REAL_T, REAL_T > (lhs.first + rhs.first, lhs.second + rhs.second));
                        // ret = lhs.second + rhs.second;
                        break;
                    case MINUS:
                        rhs = stack.top();
                        stack.pop();
                        lhs = stack.top();
                        stack.pop();
                        stack.push(std::pair<REAL_T, REAL_T > (lhs.first - rhs.first, lhs.second - rhs.second));
                        // ret = lhs.second-rhs.second;
                        break;
                    case MULTIPLY:
                        rhs = stack.top();
                        stack.pop();
                        lhs = stack.top();
                        stack.pop();
                        temp = lhs.second * rhs.first + lhs.first * rhs.second;
                        //                          std::cout << "* " <<lhs.second <<"*"<< rhs.first<<" +"<< lhs.first<<" * "<<rhs.second << " \n";

                        stack.push(std::pair<REAL_T, REAL_T > (lhs.first * rhs.first, temp));
                        // ret =temp;

                        break;
                    case DIVIDE:
                        rhs = stack.top();
                        stack.pop();
                        lhs = stack.top();
                        stack.pop();
                        temp = (lhs.second * rhs.first - lhs.first * rhs.second) / (rhs.first * rhs.first);
                        stack.push(std::pair<REAL_T, REAL_T > (lhs.first / rhs.first, temp));
                        // ret = temp;
                        break;

                    case SIN:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            stack.push(std::pair<REAL_T, REAL_T > (std::sin(lhs.first), lhs.second * std::cos(lhs.first)));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::sin(lhs.first), 0));
                        }

                        break;
                    case COS:

                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            stack.push(std::pair<REAL_T, REAL_T > (std::cos(lhs.first), lhs.second * (-1.0) * std::sin(lhs.first)));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::cos(lhs.first), 0));
                        }

                        break;
                    case TAN:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = lhs.second * ((1.0 / std::cos(lhs.first))*(1.0 / std::cos(lhs.first)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::tan(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::tan(lhs.first), 0));
                        }

                        break;
                    case ASIN:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * 1.0 / std::pow((1.0 - std::pow(lhs.first, 2.0)), 0.5));
                            stack.push(std::pair<REAL_T, REAL_T > (std::asin(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::asin(lhs.first), 0));
                        }

                        break;
                    case ACOS:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * (-1.0) / std::pow(((1.0) - std::pow(lhs.first, (2.0))), (0.5)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::acos(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::acos(lhs.first), (0)));
                        }

                        break;
                    case ATAN:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * (1.0) / (lhs.first * lhs.first + (1.0)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::atan(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::atan(lhs.first), (0)));
                        }

                        break;
                    case ATAN2:
                        rhs = stack.top();
                        stack.pop();
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (rhs.first * lhs.second / (lhs.first * lhs.first + (rhs.first * rhs.first)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::atan2(lhs.first, rhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::atan2(lhs.first, rhs.first), (0)));
                        }

                        break;
                    case SQRT:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = lhs.second * (.5) / std::sqrt(lhs.first);
                            stack.push(std::pair<REAL_T, REAL_T > (std::sqrt(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::sqrt(lhs.first), (0)));
                        }

                        break;
                    case POW:
                        rhs = stack.top();
                        stack.pop();
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * rhs.first) *
                                    std::pow(lhs.first, (rhs.first - (1.0)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::pow(lhs.first, rhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::pow(lhs.first, rhs.first), (0)));
                        }

                        break;
                    case LOG:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * (1.0)) / lhs.first;
                            stack.push(std::pair<REAL_T, REAL_T > (std::log(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::log(lhs.first), (0)));
                        }

                        break;
                    case LOG10:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * (1.0)) / (lhs.first * std::log((10.0)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::log10(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::log10(lhs.first), (0)));
                        }

                        break;
                    case EXP:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = lhs.second * std::exp(lhs.first);
                            stack.push(std::pair<REAL_T, REAL_T > (std::exp(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::exp(lhs.first), (0)));
                        }

                        break;
                    case SINH:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = lhs.second * std::cosh(lhs.first);
                            stack.push(std::pair<REAL_T, REAL_T > (std::sinh(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::sinh(lhs.first), (0)));
                        }

                        break;
                    case COSH:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = lhs.second * std::sinh(lhs.first);
                            stack.push(std::pair<REAL_T, REAL_T > (std::cosh(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::cosh(lhs.first), (0)));
                        }

                        break;
                    case TANH:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = lhs.second * ((1.0) / std::cosh(lhs.first))*((1.0) / std::cosh(lhs.first));
                            stack.push(std::pair<REAL_T, REAL_T > (std::tanh(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::tanh(lhs.first), (0)));
                        }

                        break;
                    case FABS:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * lhs.first) /
                                    std::fabs(lhs.first);
                            stack.push(std::pair<REAL_T, REAL_T > (std::fabs(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::fabs(lhs.first), (0)));
                        }

                        break;
                    case ABS:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * lhs.first) /
                                    std::fabs(lhs.first);
                            stack.push(std::pair<REAL_T, REAL_T > (std::fabs(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::fabs(lhs.first), (0)));
                        }

                        break;
                    case FLOOR:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {

                            temp = (0); //lhs.second * T(std::floor(lhs.first));
                            stack.push(std::pair<REAL_T, REAL_T > (std::floor(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::floor(lhs.first), (0)));
                        }

                        break;
                    case NONE:
                        std::cout << "nothing to do here.\n";
                        break;
                    default:
                        break;

                }


            }


            return stack.top().second;
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
         * If set to true, the expression is recorded into a vector of 
         * statements that can be manipulated to compute derivatives of 
         * arbitrary order. This functionality is costly and is rarely required,
         * thus default setting is false.
         * 
         * @param support_arbitrary_order
         */
        static void SetSupportArbitraryOrder(const bool &support_arbitrary_order) {
            Variable::is_supporting_arbitrary_order = support_arbitrary_order;
        }

        /*
         * If true, expressions are recorded into a vector of statements that 
         * can be used later to compute derivatives of arbitrary order, or
         * be used to build expression strings.
         * 
         */
        static bool IsSupportingArbitraryOrder() {
            return Variable::is_supporting_arbitrary_order;
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

        bool IsIndependent() {
            return this->is_independent_m;
        }

        /**
         * Returns the derivative with respect to a variables who's
         * unique identifier is equal to the parameter id.
         * 
         * @param id
         * @param found
         * @return 
         */
        const REAL_T Derivative(const uint32_t &id, bool &found) const {


            if (this->GetId() == id) {
                found = true;
                return 1.0;

            } else {

#if defined(USE_HASH_TABLE)
                HashTable::Cell* entry = this->gradients_m.Lookup(id);

                if (entry != NULL) {
                    return entry->value;
                } else {
                    return 0.0;
                }

#else


                const_grads_iterator git = this->gradients_m.find(id);


                if (git != this->gradients_m.end()) {
                    found = true;
                    return git->second;
                } else {
                    return 0.0;
                }

#endif
            }

        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            if (this->statements_m.size() > 0) {
                statements.insert(statements.end(), this->statements_m.begin(), statements_m.end());
            } else {
                statements.push_back(Statement<REAL_T > (VARIABLE, this->GetValue(), this->GetId()));
            }
        }

        void MakeSet(std::set<uint32_t> &independent_variables)const {
            if (this->GetId() != 0) {
                independent_variables.insert(this->GetId());
            }
            for (const_grads_iterator git = this->gradients_m.begin();
                    git != this->gradients_m.end(); ++git) {
                independent_variables.insert(git->first);
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

#if defined(USE_HASH_TABLE)
            HashTable::Cell* entry = this->gradients_m.Lookup(ind.GetId());

            if (entry != NULL) {
                return entry->value;
            } else {
                return 0.0;
            }
#else
            const_grads_iterator git = this->gradients_m.find(ind.GetId());
            if (git != this->gradients_m.end()) {
                return git->second;
            } else {
                return 0;
            }
#endif
        }

        /**
         * Finds the derivative in the encapsulated gradient map 
         * w.r.t a variable. 
         * @param ind
         * @return 
         */
        const REAL_T WRT(const Variable & ind) const {

#if defined(USE_HASH_TABLE)
            HashTable::Cell* entry = this->gradients_m.Lookup(ind.GetId());

            if (entry != NULL) {
                return entry->value;
            } else {
                return 0.0;
            }
#else
            const_grads_iterator git = this->gradients_m.find(ind.GetId());
            if (git != this->gradients_m.end()) {
                return git->second;
            } else {
                return 0;
            }
#endif
        }

        const std::string Serialize() {

        }

        void Deserialize(const std::string &str) {

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
            this->SetValue(rhs);
#if defined(USE_HASH_TABLE)
            this->gradients_m.Clear();
#else
            this->gradients_m.clear();
            if (Variable::IsRecording()) {
                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.clear();
                    this->statements_m.push_back(Statement<REAL_T > (VARIABLE, this->GetValue(), this->GetId()));
                }
            }
#endif

            return *this;
        }

        /**
         * Sets the value of this Variable to that of rhs. Derivatives 
         * are stored in the encapsulated gradient map.
         * 
         * @param rhs
         * @return 
         */
        Variable& operator=(const Variable& other) {
            //                        this->id = other.getId();
            value_m = other.GetValue();
            //            this->gradients.clear();
            if (ad::Variable<REAL_T>::is_recording_g) {

#if defined(USE_HASH_TABLE)
                ind_iterator it;
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = other.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                this->gradients.insert(rhs.gradients.begin(), rhs.gradients.end());
                ind_iterator it;
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m[(*it)] = other.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }

                std::set<uint32_t> iv;
                other.MakeSet(iv);

#endif

                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.clear();
                    other.Push(this->statements_m);
                }

            }
            this->is_independent_m = other.is_independent_m;
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
        Variable& operator=(const ExpressionBase<REAL_T, T>& expr) {
            this->id_m = expr.GetId();
            if (Variable::is_recording_g) {
                ind_iterator it; // = this->gradients.lower_bound();
#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = expr.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m[(*it)] = expr.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#endif
                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.clear();
                    expr.Push(this->statements_m);
                }

                std::set<uint32_t> iv;
                expr.MakeSet(iv);

            }
            value_m = expr.GetValue();
            return *this;
        }

        /**
         * 
         * @param rhs
         * @return 
         */
        template<class T>
        Variable& operator+=(const ExpressionBase<REAL_T, T>& rhs) {
//            return *this = (*this +rhs);
                        if (Variable::is_recording_g) {
                            ind_iterator it;
            #if defined(USE_HASH_TABLE)
                            for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                                bool found = false;
                                this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value + rhs.Derivative(*it, found);
                            }
            #else
                            for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                                bool found = false;
                                this->gradients_m[(*it)] = this->gradients_m[(*it)] + rhs.Derivative(*it, found);
                            }
            #endif
            
                            if (Variable::is_supporting_arbitrary_order) {
                                rhs.Push(this->statements_m);
                                this->statements_m.push_back(Statement<REAL_T>(PLUS));
                            }
            
                        }
                        this->value_m += rhs.GetValue();
            return *this;
        }

        Variable& operator+=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;
#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value + rhs.Derivative(*it, found);
                }
#else
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m[(*it)] = this->gradients_m[(*it)] + rhs.Derivative(*it, found);
                }
#endif

                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }

            }
            this->value_m += rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable& operator-=(const ExpressionBase<REAL_T, T>& rhs) {
            return *this = (*this -rhs);
            //            if (Variable::is_recording_g) {
            //                ind_iterator it;
            //#if defined(USE_HASH_TABLE)
            //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m.Insert((*it))->value = (this->gradients_m.Lookup(*it)->value - rhs.Derivative(*it, found));
            //                }
            //#else
            //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] - rhs.Derivative(*it, found));
            //                }
            //#endif
            //                if (Variable::is_supporting_arbitrary_order) {
            //                    rhs.Push(this->statements_m);
            //                    this->statements_m.push_back(Statement(MINUS));
            //                }
            //            }
            //            this->value_m -= rhs.GetValue();
            //            return *this;
        }

        Variable& operator-=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;
#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert((*it))->value = (this->gradients_m.Lookup(*it)->value - rhs.Derivative(*it, found));
                }
#else
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] - rhs.Derivative(*it, found));
                }
#endif
                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (MINUS));
                }
            }
            this->value_m -= rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable& operator*=(const ExpressionBase<REAL_T, T>& rhs) {
            return *this = (*this * rhs);
            //            if (Variable::is_recording_g) {
            //                ind_iterator it;
            //#if defined(USE_HASH_TABLE)
            //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#else
            //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#endif
            //                if (Variable::is_supporting_arbitrary_order) {
            //                    rhs.Push(this->statements_m);
            //                    this->statements_m.push_back(Statement(MULTIPLY));
            //                }
            //            }
            //            this->value_m *= rhs.GetValue();
            return *this;
        }

        Variable& operator*=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;
#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m[(*it)] = this->gradients_m[(*it)] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#endif
                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (MULTIPLY));
                }
            }
            this->value_m *= rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable& operator/=(const ExpressionBase<REAL_T, T>& rhs) {
            return *this = (*this / rhs);
            //            if (Variable::is_recording_g) {
            //                ind_iterator it;
            //#if defined(USE_HASH_TABLE)
            //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m.Insert(*it)->value = (this->gradients_m.Lookup(*it)->value * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#else
            //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#endif
            //                if (Variable::is_supporting_arbitrary_order) {
            //                    rhs.Push(this->statements_m);
            //                    this->statements_m.push_back(Statement(DIVIDE));
            //                }
            //            }
            //            this->value_m /= rhs.GetValue();
            return *this;
        }

        Variable& operator/=(Variable& rhs) {
            if (Variable::is_recording_g) {
                ind_iterator it;
#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = (this->gradients_m.Lookup(*it)->value * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#endif
                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (DIVIDE));
                }
            }
            this->value_m /= rhs.GetValue();
            return *this;
        }

        // And likewise for a Literal on the rhs

        Variable& operator+=(const REAL_T& rhs) {
            value_m += rhs;
            if (Variable::IsRecording()) {
                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }
            }
            return *this;
        }

        Variable& operator-=(const REAL_T& rhs) {
            if (Variable::IsRecording()) {
                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (MINUS));
                }
            }
            value_m -= rhs;
            return *this;
        }

        Variable& operator*=(const REAL_T& rhs) {
            if (Variable::IsRecording()) {
                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (MULTIPLY));
                }
            }
            return *this = (*this * rhs);
        }

        Variable& operator/=(const REAL_T& rhs) {
            if (Variable::IsRecording()) {
                if (Variable::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (DIVIDE));
                }
            }
            return *this = (*this / rhs);
        }


    };

    //list of independent variables by unique identifier.
    template<class REAL_T>
    std::set<uint32_t> Variable<REAL_T>::independent_variables_g;

    template<class REAL_T>
    bool Variable<REAL_T>::is_recording_g = true;

    template<class REAL_T>
    bool Variable<REAL_T>::is_supporting_arbitrary_order = false;

    template<class REAL_T, class T, class TT>
    inline const int operator==(const ad::ExpressionBase<REAL_T, T>& lhs, const ad::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() == rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator!=(const ad::ExpressionBase<REAL_T, T>& lhs, const ad::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() != rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator<(const ad::ExpressionBase<REAL_T, T>& lhs, const ad::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() < rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator>(const ad::ExpressionBase<REAL_T, T>& lhs, const ad::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() > rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator<=(const ad::ExpressionBase<REAL_T, T>& lhs, const ad::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() <= rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator>=(const ad::ExpressionBase<REAL_T, T>& lhs, const ad::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() >= rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator==(const REAL_T &lhs, const ad::ExpressionBase<REAL_T, T>& rhs) {
        return lhs == rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator!=(const REAL_T &lhs, const ad::ExpressionBase<REAL_T, T>& rhs) {
        return lhs != rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator<(const REAL_T &lhs, const ad::ExpressionBase<REAL_T, T>& rhs) {
        return lhs < rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator>(const REAL_T &lhs, const ad::ExpressionBase<REAL_T, T>& rhs) {
        return lhs > rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator<=(const REAL_T &lhs, const ad::ExpressionBase<REAL_T, T>& rhs) {
        return lhs <= rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator>=(const REAL_T &lhs, const ad::ExpressionBase<REAL_T, T>& rhs) {
        return lhs >= rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator==(const ad::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() == rhs;
    }

    template<class REAL_T, class T>
    inline const int operator!=(const ad::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() != rhs;
    }

    template<class REAL_T, class T>
    inline const int operator<(const ad::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() <= rhs;
    }

    template<class REAL_T, class T>
    inline const int operator>(const ad::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() > rhs;
    }

    template<class REAL_T, class T>
    inline const int operator<=(const ad::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() <= rhs;
    }

    template<class REAL_T, class T>
    inline const int operator>=(const ad::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() >= rhs;
    }


}



/**
 * Array components.
 */



#endif	/* ETAD_HPP */

