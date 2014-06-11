
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
#include <boost/container/set.hpp>
#include <tr1/unordered_set>
#define USE_MM_CACHE_MAP
//#define USE_HASH_TABLE
//this is the most crucial part of the software....we need a fast map!
#ifdef USE_MM_CACHE_MAP
#include "support/cache-table-0.2/mm/cache_map.hpp"
#include "support/cache-table-0.2/mm/cache_set.hpp"
#elif defined(USE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(USE_HASH_TABLE)
#include "support/hash_table/hashtable.h"
#elif defined(USE_JUDY_ARRAY)
#include "support/judy-code/src/Judy.h"
#else
#include <map>
#endif
#include <google/dense_hash_set>
#include <boost/unordered_set.hpp>
//typedef google::dense_hash_set<uint32_t> IdsSet;
typedef boost::unordered_set<uint32_t> IdsSet;



#include "BigFloat.hpp"

namespace ad {

    /**
     * Operation values used for recording expressions into a post-order 
     * expression tree. These operations are used primarily for supporting 
     * operations such as arbitray order derivatives,uncertainty calculations
     * and expression string building.
     */
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
        ATAN2,
        ATAN3,
        ATAN4,
        SQRT,
        POW,
        POW1,
        POW2,
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

        Statement() : op_m(CONSTANT), value_m(0), id_m(0) {

        }

        Statement(const Operation &op) : op_m(op), value_m(0), id_m(0) {

        }

        Statement(const Operation &op, const REAL_T &value) : op_m(op), value_m(value), id_m(0) {

        }

        Statement(const Operation &op, const REAL_T &value, const uint32_t &id) : op_m(op), value_m(value), id_m(id) {

        }

        Operation op_m;
        REAL_T value_m;
        uint32_t id_m;

    };

    /**
     * Interface class for storing variable information. The point of this class 
     * is to provide flexibility of the storage for the variables information. 
     * For instance, it may be desired to store the data in a database or on 
     * the disk rather than in random access memory. 
     */
    template<class REAL_T>
    class VariableStorage {
        template<class REAL_TT, int> friend class Variable;
        uint32_t id_m;
    public:

        VariableStorage() : id_m(0) {

        }

        virtual ~VariableStorage() {

        }

        uint32_t GetId() const {
            return id_m;
        }

        void SetId(uint32_t id) {
            this->id_m = id;
        }




        /**
         * Abstract function to set the value for this storage object.
         * 
         * @param value
         */
        virtual void SetValue(const REAL_T &value) = 0;

        /**
         * Abstract function to get the value for this storage object.
         * 
         * @return 
         */
        virtual const REAL_T GetValue() = 0;

        /**
         * Abstract function to return a derivative value for the independent 
         * variable with id.
         * 
         * @param id
         * @return 
         */
        virtual const std::pair<bool, REAL_T> GetDerivative(const uint32_t &id) = 0;

        /**
         * Abstract function to set the derivative value with respect to 
         * independent variable with id.
         * 
         * @param id
         * @param value
         * @return 
         */
        virtual void SetDerivative(const std::pair<bool, REAL_T> &entry, const uint32_t &id) = 0;

        /**
         * Abstract function to add a independent variable id to 
         * this storage set.
         * 
         * @param id
         * @return 
         */
        virtual void AddId(const uint32_t &id) = 0;

        /**
         * Abstract function to merge add id's from another VariableStorage 
         * to this one.
         * 
         * @param id
         * @return 
         */
        virtual void Merge(VariableStorage<REAL_T>* other) = 0;

        /**
         * Abstract function to return the size of the derivative info list. 
         * These are the computed intermediately stored derivatives. 
         * 
         */
        virtual const uint32_t DerivativeInfoSize() = 0;

        /**
         * Abstract function to clear the derivative info for this storage.
         * 
         * @return 
         */
        virtual const void ClearDerivativeInfo() = 0;

        /**
         * Abstract function to add a expression statement to the statement list. 
         * Statements are the recorded nodes in the expressions computational graph
         * and are stored in post-order(Post-order expression tree).
         * @param statement
         */
        virtual void AddStatement(const Statement<REAL_T> &statement) = 0;

        /**
         * Abstract function to return the size of the statement list size.
         * @return 
         */
        virtual const uint32_t ExpressionSize() = 0;

        /**
         * Abstract function to retrieve a Statement object at a given index.
         * 
         * @param index
         * @return 
         */
        virtual Statement<REAL_T> StatementAt(const uint32_t &index) = 0;

        /**
         * Abstract function to clear the statement list.
         */
        virtual void ClearExpression() = 0;

        /**
         * Rests this storage.
         */
        virtual void Reset() = 0;

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

        const uint32_t current() {
            return _id;
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

    template<class REAL_T, int group = 0 >
    class Variable;

    /**
     * Base class for expression types.
     */
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

        inline const REAL_T GetValue() const {
            return Cast().GetValue();
        }

        inline const uint32_t GetId() const {
            return id_m;
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            Cast().GetIdRange(min, max);
        }

        /**
         * Compute or get the stored derivative with respect to the unique 
         * identifier id.
         * 
         * @param id
         * @param found
         * @return 
         */
        const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return Cast().Derivative(id, found);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            Cast().Push(statements);
        }

        void PushIds(IdsSet & ids) const {
            Cast().PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            Cast().PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > &ids) const {
            Cast().PushIds(ids);
        }

        uint32_t id_m;

    private:

        ExpressionBase& operator=(const ExpressionBase & exp) const {

            return *this;
        };
    };

    /**
     * Expression template for constant values.
     * 
     * @param value
     */
    template<class REAL_T>
    struct Constant : public ExpressionBase<REAL_T, Constant<REAL_T> > {

        Constant(const REAL_T & value) : value_m(value) {

        }

        inline const REAL_T GetValue() const {
            return this->value_m;

        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return 0;
        };

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, value_m));
        }

        inline void PushIds(IdsSet & ids) const {

        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {

        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {

        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {

        }

        const REAL_T &value_m;


    };

    /**
     * 
     * Expression template for addition.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS, class RHS>
    struct Add : public ExpressionBase<REAL_T, Add<REAL_T, LHS, RHS> > {

        Add(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), value_m(lhs_m.GetValue() + rhs_m.GetValue()) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return lhs_m.Derivative(id, found) + rhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (PLUS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }




    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for addition of two expression templates.
     * @param a
     * @param b
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline
    Add<REAL_T, LHS, RHS> operator+(const ExpressionBase<REAL_T, LHS>& a,
    const ExpressionBase<REAL_T, RHS>& b) {
        return Add<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    /**
     * Expression template for subtraction.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS, class RHS>
    struct Minus : public ExpressionBase<REAL_T, Minus<REAL_T, LHS, RHS> > {

        Minus(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), value_m(lhs_m.GetValue() - rhs_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return lhs_m.Derivative(id, found) - rhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MINUS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }


    protected:


    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for expression template subtraction.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline const
    Minus<REAL_T, LHS, RHS> operator-(const ExpressionBase<REAL_T, LHS>& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return Minus<REAL_T, LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    /**
     * Expression template for multiplying two expression templates.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS, class RHS>
    struct Multiply : public ExpressionBase<REAL_T, Multiply<REAL_T, LHS, RHS> > {

        Multiply(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), value_m(lhs_m.GetValue() * rhs_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            //            std::cout << "* " << a_.Derivative(id, found) << "*" << b_ << " +" <<a_ << " * " <<  b_.Derivative(id, found) << " \n";
            return lhs_m.Derivative(id, found) * rhs_m.GetValue() + lhs_m.GetValue() * rhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MULTIPLY));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }



    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for expression template multiplication.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline const
    Multiply<REAL_T, LHS, RHS> operator*(const ExpressionBase<REAL_T, LHS>& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return Multiply<REAL_T, LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    /**
     * Expression template for division.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS, class RHS>
    struct Divide : public ExpressionBase<REAL_T, Divide<REAL_T, LHS, RHS> > {

        Divide(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), value_m(lhs_m.GetValue() / rhs_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return (lhs_m.Derivative(id, found) * rhs_m.GetValue() - lhs_m.GetValue() * rhs_m.Derivative(id, found)) / (rhs_m.GetValue() * rhs_m.GetValue());
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (DIVIDE));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
            this->rhs_m.PushIds(ids);
        }


    private:

        const LHS& lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for expression template division.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline const
    Divide<REAL_T, LHS, RHS> operator/(const ExpressionBase<REAL_T, LHS>& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return Divide<REAL_T, LHS, RHS > (lhs.Cast(), rhs.Cast());
    }

    /**
     * Expression template for adding a constant to a expression template.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS>
    struct ConstantAdd : public ExpressionBase<REAL_T, ConstantAdd<REAL_T, LHS> > {

        ConstantAdd(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs_m.GetValue() + rhs) {
        }

        inline const REAL_T GetValue() const {
            //            this->value_m = lhs_m.GetValue() + rhs_m;
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return lhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (PLUS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
        }

        inline inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
        }



    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
        REAL_T value_m;
    };

    /**
     * Operator for adding a constant value to a expression templates.
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS>
    inline const
    ConstantAdd<REAL_T, LHS> operator+(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return ConstantAdd<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    inline const
    ConstantAdd<REAL_T, RHS> operator+(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return ConstantAdd<REAL_T, RHS > (rhs.Cast(), lhs);
    }

    /**
     * Expression template for subtracting a constant from a expression 
     * template.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS>
    struct ConstantMinus : public ExpressionBase<REAL_T, ConstantMinus<REAL_T, LHS> > {

        ConstantMinus(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs_m.GetValue() - rhs_m) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return lhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (MINUS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
        }

        inline inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
        }




    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
        REAL_T value_m;
    };

    /**
     * Expression template for subtracting a expression template from a 
     * constant.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class RHS>
    struct MinusConstant : public ExpressionBase<REAL_T, MinusConstant<REAL_T, RHS> > {

        MinusConstant(const REAL_T & lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : rhs_m(rhs.Cast()), lhs_m(lhs), value_m(lhs_m - rhs_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return -1.0 * rhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, lhs_m));
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MINUS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->rhs_m.PushIds(ids);
        }


    private:
        const RHS& rhs_m;
        const REAL_T lhs_m;
        const REAL_T value_m;

    };

    /**
     * Operator for subtracting a constant from a expression template.
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS>
    inline const
    ConstantMinus<REAL_T, LHS> operator-(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return ConstantMinus<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    /**
     * Operator for subtracting a expression template from a constant.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS>
    inline const
    MinusConstant<REAL_T, RHS> operator-(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return MinusConstant<REAL_T, RHS > (lhs, rhs.Cast());
    }

    /**
     * Expression template for multiplying a constant by a expression template.
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS>
    struct ConstantTimes : public ExpressionBase<REAL_T, ConstantTimes<REAL_T, LHS> > {

        ConstantTimes(const REAL_T & lhs, const ExpressionBase<REAL_T, LHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m * rhs_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return rhs_m.Derivative(id, found) * lhs_m;
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, lhs_m));
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (MULTIPLY));
        }

        inline void PushIds(IdsSet & ids) const {
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->rhs_m.PushIds(ids);
        }



    private:
        const REAL_T lhs_m;
        const LHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Expression template for multiplying a expression template by a constant.
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS>
    struct TimesConstant : public ExpressionBase<REAL_T, TimesConstant<REAL_T, LHS> > {

        TimesConstant(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs_m.GetValue() * rhs_m) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return rhs_m * lhs_m.Derivative(id, found);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (MULTIPLY));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
        }


    private:
        const LHS& lhs_m;
        const REAL_T rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for multiplying a expression template by a constant.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS>
    inline const
    TimesConstant<REAL_T, LHS> operator*(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return TimesConstant<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    /**
     * Operator for multiplying a constant by a expression template.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS>
    inline const
    ConstantTimes<REAL_T, RHS> operator*(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return ConstantTimes<REAL_T, RHS > (lhs, rhs.Cast());
    }

    /**
     * Expression template for dividing a constant by a expression template.
     * 
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class RHS>
    struct ConstantDivide : public ExpressionBase<REAL_T, ConstantDivide<REAL_T, RHS> > {

        ConstantDivide(const REAL_T & lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : rhs_m(rhs.Cast()), lhs_m(lhs), value_m(lhs_m / rhs_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return ( -rhs_m.Derivative(id, found) * lhs_m) / (rhs_m.GetValue() * rhs_m.GetValue());
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            rhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, lhs_m));
            this->rhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (DIVIDE));
        }

        inline void PushIds(IdsSet & ids) const {
            this->rhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->rhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->rhs_m.PushIds(ids);
        }



    private:
        const RHS& rhs_m;
        REAL_T lhs_m;
        REAL_T value_m;
    };

    /**
     * Expression template for dividing a expression template by a constant.
     * @param lhs
     * @param rhs
     */
    template <class REAL_T, class LHS>
    struct DivideConstant : public ExpressionBase<REAL_T, DivideConstant<REAL_T, LHS> > {

        DivideConstant(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs_m.GetValue() / rhs_m) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return (lhs_m.Derivative(id, found) * rhs_m) / (rhs_m * rhs_m);
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            lhs_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->lhs_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, rhs_m));
            statements.push_back(Statement<REAL_T > (DIVIDE));
        }

        inline void PushIds(IdsSet & ids) const {
            this->lhs_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->lhs_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->lhs_m.PushIds(ids);
        }


    private:
        const LHS& lhs_m;
        REAL_T rhs_m;
        REAL_T value_m;
    };

    /**
     * Operator for dividing a expression template by a constant.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class LHS>
    inline const
    DivideConstant<REAL_T, LHS> operator/(const ExpressionBase<REAL_T, LHS>& lhs,
    const REAL_T& rhs) {
        return DivideConstant<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    /**
     * Operator for dividing a constant by a expression template.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS>
    inline const
    ConstantDivide<REAL_T, RHS> operator/(const REAL_T& lhs,
    const ExpressionBase<REAL_T, RHS>& rhs) {
        return ConstantDivide<REAL_T, RHS > (lhs, rhs.Cast());
    }

    /**
     * Expression template for taking the sine of an expression.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    struct Sin : public ExpressionBase<REAL_T, Sin<REAL_T, EXPR> > {

        Sin(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::sin(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * std::cos(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (SIN));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }



    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for taking the cosine of an expression.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Cos : public ExpressionBase<REAL_T, Cos<REAL_T, EXPR> > {

        Cos(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::cos(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx*-1.0 * std::sin(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (COS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }



    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the tangent of a expression template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Tan : public ExpressionBase<REAL_T, Tan<REAL_T, EXPR> > {

        Tan(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::tan(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return dx * ((1.0 / std::cos(v))*(1.0 / std::cos(v)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (TAN));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the inverse Sine of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct ASin : public ExpressionBase<REAL_T, ASin<REAL_T, EXPR> > {

        ASin(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::asin(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * static_cast<REAL_T> (1.0) / std::pow((static_cast<REAL_T> (1.0) - std::pow(v, static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ASIN));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }



    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the inverse cosine of an expression
     * template.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    struct ACos : public ExpressionBase<REAL_T, ACos<REAL_T, EXPR> > {

        ACos(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::acos(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * static_cast<REAL_T> (-1.0) / std::pow((static_cast<REAL_T> (1.0) - std::pow(v, static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ACOS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }



    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the inverse tangent of an expression 
     * template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct ATan : public ExpressionBase<REAL_T, ATan<REAL_T, EXPR> > {

        ATan(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * 1.0 / (v * v + 1.0));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ATAN));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both arguments are expression templates. 
     * 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    struct ATan2 : public ExpressionBase<REAL_T, ATan2<REAL_T, EXPR1, EXPR2> > {

        ATan2(const ExpressionBase<REAL_T, EXPR1>& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan2(expr1_m.GetValue(), expr2_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            REAL_T dx = expr1_m.Derivative(id, found);

            if (found) {
                REAL_T v = expr1_m.GetValue();
                REAL_T v2 = expr2_m.GetValue();

                return ((v2 * dx) / (v * v + (v2 * v2)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr1_m.GetIdRange(min, max);
            expr2_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr1_m.Push(statements);
            this->expr2_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ATAN));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr1_m.PushIds(ids);
            this->expr2_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr1_m.PushStorage(ids);
            this->expr2_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr1_m.PushIds(ids);
            this->expr2_m.PushIds(ids);
        }


    private:

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    /**
     *Expression template for computing the two argument inverse tangent, where
     * the second argument is a constant. 
     * 
     */
    template <class REAL_T, class EXPR1>
    struct ATan2Constant : public ExpressionBase<REAL_T, ATan2Constant<REAL_T, EXPR1> > {

        ATan2Constant(const ExpressionBase<REAL_T, EXPR1>& expr1, const REAL_T & expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2) {
        }

        inline const REAL_T GetValue() const {
            return std::atan2(expr1_m.GetValue(), expr2_m);
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            REAL_T dx = expr1_m.Derivative(id, found);

            if (found) {
                REAL_T v = expr1_m.GetValue();
                REAL_T v2 = expr2_m;

                return ((v2 * dx) / (v * v + (v2 * v2)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr1_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr1_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, expr2_m));
            statements.push_back(Statement<REAL_T > (ATAN2));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr1_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr1_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr1_m.PushIds(ids);
        }


    private:

        const EXPR1& expr1_m;
        const REAL_T& expr2_m;
    };

    /**
     *Expression template for computing the two argument inverse tangent,
     * where the first argument is a constant. 
     * 
     */
    template <class REAL_T, class EXPR2>
    struct ConstantATan2 : public ExpressionBase<REAL_T, ConstantATan2<REAL_T, EXPR2> > {

        ConstantATan2(const REAL_T& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan2(expr1_m, expr2_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            REAL_T dx = 0;

            if (found) {
                REAL_T v = expr1_m;
                REAL_T v2 = expr2_m.GetValue();

                return ((v2 * dx) / (v * v + (v2 * v2)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr2_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, expr1_m));
            expr2_m.Push(statements);
            statements.push_back(Statement<REAL_T > (ATAN2));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr2_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr2_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr2_m.PushIds(ids);
        }


    private:

        const REAL_T& expr1_m;
        const ExpressionBase<REAL_T, EXPR2>& expr2_m;
    };

    /**
     * Expression template for computing the square root of an expression
     * template.
     * @param a
     */
    template <class REAL_T, class EXPR>
    struct Sqrt : public ExpressionBase<REAL_T, Sqrt<REAL_T, EXPR> > {

        Sqrt(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::sqrt(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * .5 / std::sqrt(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (SQRT));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:
        const EXPR& expr_m;
    };

    /**
     *Expression template for computing the power of a expression template,
     * where both arguments are expression templates. 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    struct Pow : public ExpressionBase<REAL_T, Pow<REAL_T, EXPR1, EXPR2> > {

        Pow(const ExpressionBase<REAL_T, EXPR1>& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr1_m.GetValue(), expr2_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
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

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr1_m.GetIdRange(min, max);
            expr2_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr1_m.Push(statements);
            this->expr2_m.Push(statements);
            statements.push_back(Statement<REAL_T > (POW));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr1_m.PushIds(ids);
            this->expr2_m.Push(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr1_m.PushStorage(ids);
            this->expr2_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr1_m.PushIds(ids);
            this->expr2_m.Push(ids);
        }


    private:

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    /**
     * Expression template for computing the power of a expression template, 
     * where the second argument(the power) is a constant.
     * 
     * @param expr
     * @param constant
     */
    template <class REAL_T, class EXPR>
    struct PowConstant : public ExpressionBase<REAL_T, PowConstant<REAL_T, EXPR> > {

        PowConstant(const ExpressionBase<REAL_T, EXPR>& expr, const REAL_T & constant)
        : expr_m(expr.Cast()), constant_m(constant) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr_m.GetValue(), constant_m);
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                REAL_T v2 = constant_m;

                return (dx * v2) *std::pow(v, (v2 - static_cast<REAL_T> (1.0)));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CONSTANT, constant_m));
            statements.push_back(Statement<REAL_T > (POW));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:

        const EXPR& expr_m;
        const REAL_T& constant_m;
    };

    /**
     * Expression template to compute the the value of a constant raised to 
     * an expression template.
     * 
     * @param constant
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct ConstantPow : public ExpressionBase<REAL_T, ConstantPow<REAL_T, EXPR> > {

        ConstantPow(const REAL_T& constant, const ExpressionBase<REAL_T, EXPR>& expr)
        : constant_m(constant), expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(constant_m, expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
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

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.push_back(Statement<REAL_T > (CONSTANT, constant_m));
            this->expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (POW));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:

        const REAL_T& constant_m;
        const ExpressionBase<REAL_T, EXPR>& expr_m;
    };

    /**
     * Expression template to compute the log of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Log : public ExpressionBase<REAL_T, Log<REAL_T, EXPR> > {

        Log(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::log(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);

            if (found) {
                return (dx) / expr_m.GetValue();
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (LOG));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template to compute the log base 10 of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Log10 : public ExpressionBase<REAL_T, Log10<REAL_T, EXPR> > {

        Log10(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::log10(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return (dx * 1.0) / (expr_m.GetValue() * std::log(10.0));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (LOG10));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template to compute e raised to a expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Exp : public ExpressionBase<REAL_T, Exp<REAL_T, EXPR> > {

        Exp(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::exp(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            bool f = false;
            REAL_T dx = expr_m.Derivative(id, f);
            if (f) {
                found = true;
                return dx * std::exp(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (EXP));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template used to protect overflow in exp calculations. 
     * 
     * Author: Dave Fournier.
     * Original implementation in ADMB.
     * 
     * Source: http://admb-project.org/documentation/api/mfexp_8cpp_source.html
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct MFExp : public ExpressionBase<REAL_T, MFExp<REAL_T, EXPR> > {

        MFExp(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(Compute(expr.GetValue())) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
            //
            //            return std::exp(expr_m.GetValue());
        }

        const REAL_T Compute(const REAL_T & value) {
            REAL_T x = value;
            REAL_T b = REAL_T(60);
            if (x <= b && x >= REAL_T(-1) * b) {
                return std::exp(x);
            } else if (x > b) {
                return std::exp(b)*(REAL_T(1.) + REAL_T(2.) * (x - b)) / (REAL_T(1.) + x - b);
            } else {
                return std::exp(REAL_T(-1) * b)*(REAL_T(1.) - x - b) / (REAL_T(1.) + REAL_T(2.) * (REAL_T(-1) * x - b));
            }
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * this->GetValue();
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (EXP));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:
        const EXPR& expr_m;
        REAL_T value_m;
    };

    /**
     * Expression template for computing the hyperbolic sine of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Sinh : public ExpressionBase<REAL_T, Sinh<REAL_T, EXPR> > {

        Sinh(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::sinh(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * std::cosh(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (SINH));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the hyperbolic cosine of an 
     * expression template.
     */
    template <class REAL_T, class EXPR>
    struct Cosh : public ExpressionBase<REAL_T, Cosh<REAL_T, EXPR> > {

        Cosh(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::cosh(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                return dx * std::sinh(expr_m.GetValue());
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (COSH));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }


    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for computing the hyperbolic tangent of an expresison
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Tanh : public ExpressionBase<REAL_T, Tanh<REAL_T, EXPR> > {

        Tanh(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::tanh(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return dx * (1.0 / std::cosh(v))*(1.0 / std::cosh(v));
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (TANH));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for handling the absolute value of a expression 
     * template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Fabs : public ExpressionBase<REAL_T, Fabs<REAL_T, EXPR> > {

        Fabs(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::fabs(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {

            REAL_T dx = expr_m.Derivative(id, found);
            if (found) {
                REAL_T v = expr_m.GetValue();
                return (dx * v) / std::fabs(v);
            } else {
                return 0.0;
            }
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (FABS));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for handling the floor of an expression template.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    struct Floor : public ExpressionBase<REAL_T, Floor<REAL_T, EXPR> > {

        Floor(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::floor(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return 0.0;
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (FLOOR));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:
        const EXPR& expr_m;
    };

    /**
     * Expression template for handling the ceiling of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    struct Ceil : public ExpressionBase<REAL_T, Ceil<REAL_T, EXPR> > {

        Ceil(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::ceil(expr_m.GetValue());
        }

        inline const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return 0.0;
        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {
            expr_m.GetIdRange(min, max);
        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            expr_m.Push(statements);
            statements.push_back(Statement<REAL_T > (CEIL));
        }

        inline void PushIds(IdsSet & ids) const {
            this->expr_m.PushIds(ids);
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {
            this->expr_m.PushStorage(ids);
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            this->expr_m.PushIds(ids);
        }

    private:
        const EXPR& expr_m;
    };


}



/**
 * Utility functions added to the standard name space. Mostly cmath overloads 
 * for expressions.
 */
namespace std {

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
    inline const ad::Sin<REAL_T, EXPR> sin(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Sin<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the cos function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Cos<REAL_T, EXPR> cos(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Cos<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the tan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Tan<REAL_T, EXPR> tan(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Tan<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::ASin<REAL_T, EXPR> asin(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::ASin<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::ACos<REAL_T, EXPR> acos(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::ACos<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the atan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::ATan<REAL_T, EXPR> atan(const ad::ExpressionBase<REAL_T, EXPR>& a) {
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
    ad::ATan2Constant<REAL_T, EXPR> atan2(const ad::ExpressionBase<REAL_T, EXPR>& expr,
    const REAL_T& val) {
        return ad::ATan2Constant<REAL_T, EXPR > (expr.Cast(), val);
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
    ad::ConstantATan2<REAL_T, EXPR> atan2(const REAL_T& val,
    const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::ConstantATan2<REAL_T, EXPR > (val, expr.Cast());
    }

    /**
     * Override for the sqrt function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Sqrt<REAL_T, EXPR> sqrt(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
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
    ad::PowConstant<REAL_T, EXPR> pow(const ad::ExpressionBase<REAL_T, EXPR>& expr,
    const REAL_T& val) {
        return ad::PowConstant<REAL_T, EXPR > (expr.Cast(), val);
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
    ad::ConstantPow<REAL_T, EXPR> pow(const REAL_T& val,
    const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::ConstantPow<REAL_T, EXPR > (val, expr.Cast());
    }

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Log<REAL_T, EXPR> log(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Log<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Log10<REAL_T, EXPR> log10(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Log10<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the exp function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Exp<REAL_T, EXPR> exp(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Exp<REAL_T, EXPR > (expr.Cast());
    }

    template<class REAL_T, class EXPR>
    inline const const ad::MFExp<REAL_T, EXPR> mfexp(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::MFExp<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the sinh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Sinh<REAL_T, EXPR> sinh(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Sinh<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the cosh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Cosh<REAL_T, EXPR> cosh(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Cosh<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the tanh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Tanh<REAL_T, EXPR> tanh(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Tanh<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the fabs function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Fabs<REAL_T, EXPR> fabs(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Fabs<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the floor function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Floor<REAL_T, EXPR> floor(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Floor<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const ad::Ceil<REAL_T, EXPR> ceil(const ad::ExpressionBase<REAL_T, EXPR>& expr) {
        return ad::Ceil<REAL_T, EXPR > (expr.Cast());
    }



}

namespace ad {

    /**
     * Default implementation of the abstract VariableStorage class.
     */
    template<class REAL_T>
    class DefaultStorage : public VariableStorage<REAL_T> {
        typedef std::vector<std::pair<bool, REAL_T> > GradientVector;
        GradientVector g;

        typedef IdsSet indepedndent_variables;
        typedef IdsSet::iterator indepedndent_variables_iterator;
        typedef IdsSet::const_iterator const_indepedndent_variables_iterator;
        indepedndent_variables ids_m;

        typedef std::vector<Statement<REAL_T> > ExpressionStatements;
        ExpressionStatements statements_m;

        REAL_T value_m;
    public:

        DefaultStorage() : value_m(0) {

        }

        void SetValue(const REAL_T &value) {
            this->value_m = value;
        }

        const REAL_T GetValue() {
            return this->value_m;
        }

        const std::pair<bool, REAL_T> GetDerivative(const uint32_t &id) {

            if (id> this->g.size()) {
                return g[id];
            } else {
                return std::pair<bool, REAL_T > (false, 0);
            }

        }

        void SetDerivative(const std::pair<bool, REAL_T> &entry, const uint32_t &id) {

            if (g.size() < id) {
                g.resize(id + 1);
            }
            g[id] = entry;
        }

        void AddId(const uint32_t &id) {
            this->ids_m.insert(id);
        }

        void Merge(VariableStorage<REAL_T>* other) {
            DefaultStorage<REAL_T>* ds = (DefaultStorage<REAL_T>*)other;
            this->ids_m.insert(ds->ids_m.begin(), ds->ids_m.end());
        }

        const uint32_t DerivativeInfoSize() {
            return this->g.size();
        }

        const void ClearDerivativeInfo() {
            this->g.clear();
        }

        void AddStatement(const Statement<REAL_T> &statement) {
            this->statements_m.push_back(statement);
        }

        const uint32_t ExpressionSize() {
            return this->statements_m.size();
        }

        Statement<REAL_T> StatementAt(const uint32_t &index) {
            //            if(index >= this->statements_m.size()){
            //                std::cout<<"Error Statement<REAL_T> StatementAt(const uint32_t &index), index out of bounds."
            //            }

            return this->statements_m[index];
        }

        void ClearExpression() {
            this->statements_m.clear();
        }

        void Reset() {
            this->g.clear();
            this->ids_m.clear();
            this->statements_m.clear();
        }

    };

    /**
     * Variable class used for calculations and storage of computed derivatives.
     * 
     * Derivatives are only computed for variables marked as independent 
     * variables. To make a variable independent use:
     * 
     * constructor:
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
    template<class REAL_T, int group>
    class Variable : public ExpressionBase<REAL_T, Variable<REAL_T, group> > {
        VariableStorage<REAL_T>* storage;
        REAL_T value_m;

        std::string name_m;
        bool bounded_m;
        REAL_T min_boundary_m;
        REAL_T max_boundary_m;
        bool is_independent_m;
        uint32_t iv_id_m; //id when is a independent variable.

        //        uint32_t iv_min; //if this variable is caching derivatives, this is the min independent variable.
        //        uint32_t iv_max; //if this variable is caching derivatives, this is the max independent variable.


        //        static std::set<uint32_t> independent_variables_g;
        static bool is_recording_g;

        static bool is_supporting_arbitrary_order;


        typedef std::vector<std::pair<bool, REAL_T> > GradientVector;
        GradientVector g;
        //        GradientMap gradients_m;

        typedef IdsSet indepedndent_variables;
        typedef IdsSet::iterator indepedndent_variables_iterator;
        typedef IdsSet::const_iterator const_indepedndent_variables_iterator;
        indepedndent_variables ids_m;
        //        std::vector<bool> has_m;
        typedef std::vector<Statement<REAL_T> > ExpressionStatements;
        ExpressionStatements statements_m;

        template <typename TT >
        TT SwapBytes(const TT &u) {

            union {
                TT u;
                unsigned char u8[sizeof (TT)];
            } source, dest;

            source.u = u;

            for (size_t k = 0; k < sizeof (TT); k++)
                dest.u8[k] = source.u8[sizeof (TT) - k - 1];

            return dest.u;
        }


    public:
        static uint32_t misses_g;

        /**
         * Default conclassor.
         */
        Variable() : ExpressionBase<REAL_T, Variable<REAL_T> >(0),
        value_m(0.0),
        bounded_m(false),
        is_independent_m(false),
        storage(new DefaultStorage<REAL_T>()),
        iv_id_m(0) {

            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();
            //this->ids_m.set_empty_key(NULL);
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
        Variable(const REAL_T& value, bool is_independent = false) : storage(new DefaultStorage<REAL_T>()), value_m(value), is_independent_m(is_independent), iv_id_m(0) {
            //this->ids_m.set_empty_key(NULL);
            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();
            if (is_independent) {
                //                this->id_m = IDGenerator::instance()->next();
                this->SetAsIndependent(is_independent);
                this->bounded_m = false;
                this->min_boundary_m = std::numeric_limits<REAL_T>::min();
                this->max_boundary_m = std::numeric_limits<REAL_T>::max();
                //                Variable::independent_variables_g.insert(this->GetId());
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
        Variable(const Variable& orig) : storage(new DefaultStorage<REAL_T>()) {

            value_m = orig.GetValue();
            this->id_m = orig.GetId();
            this->iv_id_m = orig.iv_id_m;
            //            orig.PushIds(ids_m);
            //this->ids_m.set_empty_key(NULL);
            ids_m.insert(orig.ids_m.begin(), orig.ids_m.end()); // = orig.ids_m;
            //            g = orig.g;
            g.reserve(orig.g.size());
            g.insert(g.begin(), orig.g.begin(), orig.g.end());
            //            has_m = orig.has_m;
#if defined(USE_HASH_TABLE)
            this->gradients_m = HashTable(orig.gradients_m);
#else
            //            this->gradients_m.insert(orig.gradients_m.begin(), orig.gradients_m.end());
#endif
            this->is_independent_m = orig.is_independent_m;
            this->bounded_m = orig.bounded_m;
            this->min_boundary_m = orig.min_boundary_m;
            this->max_boundary_m = orig.max_boundary_m;
            //            iv_min = orig.iv_min;
            //            iv_max = orig.iv_max;
            if (Variable::IsSupportingArbitraryOrder()) {
                orig.Push(this->statements_m);
            }

        }

        /**
         * Conclasss a variable from expression expr.
         * @param rhs
         */
        template<class T>
        Variable(const ExpressionBase<REAL_T, T>& expr) : storage(new DefaultStorage<REAL_T>()) {

            //has_m.resize(IDGenerator::instance()->current() + 1);

            this->bounded_m = false;
            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();
            if (Variable::is_recording_g) {
                //                this->id_m = expr.GetId();
                //                ind_iterator it;
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //#if defined(USE_HASH_TABLE)
                //                    this->gradients_m.Insert(*it)->value = expr.Derivative(*it, found);
                //#else
                //                    this->gradients_m[(*it)] = expr.Derivative(*it, found);
                //#endif
                //                
                //            }

                //this->ids_m.set_empty_key(NULL);
                expr.PushIds(ids_m);
                //expr.GetIdRange(this->iv_min, this->iv_max);
                indepedndent_variables_iterator it;
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = false;
                    this->g[*it].first = true;
                    this->g[*it].second = expr.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }
                if (Variable::IsSupportingArbitraryOrder()) {
                    expr.Push(this->statements_m);
                }
            }
            value_m = expr.GetValue();

        }

        ~Variable() {
            delete storage;
#if !defined(USE_HASH_TABLE)
            //            if (this->is_independent_m) {
            //                this->independent_variables_g.erase(this->GetId());
            //            }
            //            this->statements_m.clear();
            //            this->gradients_m.clear();

#endif
        }

        size_t Size() {
            return this->statements_m.size();
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
        inline const REAL_T GetValue() const {
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
            if (this->iv_id_m == 0) {
                this->iv_id_m = IDGenerator::instance()->next();
                if (Variable<REAL_T>::IsSupportingArbitraryOrder()) {
                    this->statements_m.clear();
                    this->statements_m.push_back(Statement<REAL_T > (VARIABLE, this->GetValue(), this->GetId()));
                }
                //                this->id_m = iv_id_m;
                //                if (this->iv_min == std::numeric_limits<uint32_t>::max() && this->iv_max == std::numeric_limits<uint32_t>::min()) {
                //                    this->iv_max = this->GetId();
                //                    this->iv_min = this->GetId();
                //                }
            }
            if (this->is_independent_m && !is_independent) {
                //                Variable::independent_variables_g.erase(this->GetId());
                this->id_m = 0;
                this->is_independent_m = false;
            }

            if (!this->is_independent_m && is_independent) {
                //                Variable::independent_variables_g.insert(this->GetId());
                this->id_m = this->iv_id_m;
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



            if (this->Cast().id_m == id) {
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

                //                if(id < this->g.size()){
                //                    if(g[id] != 0){
                //                        found = true;
                //                    }
                //                    
                //                    return g[id];
                //                }

                //                std::find(this->ids_m.begin(),this->ids_m.end(),id );

                //                if (ad::binary_search(ids_m.begin(), ids_m.end(), id)) {
                //                    found = true;
                //                    return g[id];
                //                } else {
                //                    return 0.0;
                //                }
                //                
                if (id < g.size()) {
                    if (g[id].first) {
                        found = true;
                        return g[id].second;
                    } else {
                        return 0.0;
                    }
                } else {
                    return 0.0;
                }
                //                if (this->ids_m.find(id) != this->ids_m.end()) {
                //                    found = true;
                //                    return g[id].second;
                //                } else {
                //                    return 0.0;
                //                }

                //                 const_grads_iterator git = this->gradients_m.find(id);
                //                if (git != this->gradients_m.end()) {
                //                    found = true;
                //                    return git->second;
                //                } else {
                //                    Variable<REAL_T>::misses_g++;
                //                    return 0.0;
                //                }
            }

#endif


        }

        void GetIdRange(uint32_t &min, uint32_t & max) const {

            if (this->iv_min < min) {
                min = this->iv_min;
            }

            if (this->iv_max > max) {
                max = this->iv_max;
            }

            if (this->GetId() > 0) {
                if (this->GetId() < min) {
                    min = this->GetId();
                }

                //                if (this->GetId()) {
                //                    max = this->GetId();
                //                }
            }
            //
            //                        if (min == std::numeric_limits<uint32_t>::max()) {
            //                            min = max;
            //                        }

            //            std::cout << min << " " << max << "\n";



        }

        void Push(std::vector<Statement<REAL_T> > &statements) const {
            statements.insert(statements.end(), this->statements_m.begin(), statements_m.end());
            //            if (this->statements_m.size() > 0) {
            //                statements.insert(statements.end(), this->statements_m.begin(), statements_m.end());
            //            } else {
            //                statements.push_back(Statement<REAL_T > (VARIABLE, this->GetValue(), this->GetId()));
            //            }
        }

        inline void PushIds(IdsSet &ids) const {
            if (this->GetId() != 0) {
                //                uint32_t id = uint32_t
                ids.insert(this->GetId());
            } else {
                ids.insert(this->ids_m.begin(), this->ids_m.end());
            }
        }

        inline void PushStorage(VariableStorage<REAL_T> * ids) const {

            if (Variable<REAL_T, group>::IsRecording()) {

                if (this->GetId() != 0) {
                    ids->AddId(this->GetId());
                } else {
                    ids->Merge(this->storage);
                }


                for (int i = 0; i < this->storage->DerivativeInfoSize(); i++) {
                    ids->SetDerivative(this->storage->GetDerivative(i), i);
                }

                if (Variable<REAL_T, group>::IsSupportingArbitraryOrder()) {
                    for (int i = 0; i < this->storage->ExpressionSize(); i++) {
                        ids->AddStatement(this->storage->StatementAt(i));
                    }
                }

            }
        }

        inline void PushIds(std::vector < std::pair<bool, REAL_T> > & ids) const {
            for (int i = 1; i< this->g.size(); i++) {
                if (g[i].first) {
                    ids[i].first = true;
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

#if defined(USE_HASH_TABLE)
            HashTable::Cell* entry = this->gradients_m.Lookup(ind.GetId());

            if (entry != NULL) {
                return entry->value;
            } else {
                return 0.0;
            }
#else

            if (ind.GetId() < g.size()) {
                if (g[ind.GetId()].first) {

                    return g[ind.GetId()].second;
                } else {
                    return 0.0;
                }
            } else {
                return 0.0;
            }
            //            if (this->ids_m.find(ind.GetId()) != this->ids_m.end()) {
            //                //                found = true;
            //                return g[ind.GetId()].second;
            //            } else {
            //                return 0.0;
            //            }

            //            const_grads_iterator git = this->gradients_m.find(ind.GetId());
            //            if (git != this->gradients_m.end()) {
            //                return git->second;
            //            } else {
            //                return 0;
            //            }
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

            if (ind.GetId() < g.size()) {
                if (g[ind.GetId()].first) {

                    return g[ind.GetId()].second;
                } else {
                    return 0.0;
                }
            } else {
                return 0.0;
            }
            //            if (this->ids_m.find(ind.GetId()) != this->ids_m.end()) {
            //                //                found = true;
            //                return g[ind.GetId()].second;
            //            } else {
            //                return 0.0;
            //            }
            //            const_grads_iterator git = this->gradients_m.find(ind.GetId());
            //            if (git != this->gradients_m.end()) {
            //                return git->second;
            //            } else {
            //                return 0;
            //            }
#endif
        }

        void Serialize(std::ostream &out) {

            bool little_endian = true;

            int num = 1;
            if (*(char *) &num == 1) {
                little_endian = true;
            } else {
                little_endian = false;
            }

            //id stuff
            uint32_t id = this->GetId();


            if (!little_endian) {
                id = SwapBytes<uint32_t > (id);
            }

            out.write(reinterpret_cast<const char*> (&id), sizeof (uint32_t));

            //name
            uint32_t namesize = (uint32_t)this->GetName().size();

            if (!little_endian) {
                namesize = SwapBytes<uint32_t > (namesize);
            }
            out.write(reinterpret_cast<const char*> (&namesize), sizeof (uint32_t));


            out.write(this->GetName().data(), namesize);
            //is independent
            if (this->IsIndependent()) {
                out << '1';
            } else {
                out << '0';
            }
            //value stuff

            //bounded
            if (this->IsBounded()) {
                out << '1';
            } else {
                out << '0';
            }
            //min
            REAL_T value = this->GetMinBoundary();
            if (!little_endian) {
                value = SwapBytes<REAL_T > (value);
            }

            out.write(reinterpret_cast<const char*> (&value), sizeof ( REAL_T));

            //max
            value = this->GetMaxBoundary();

            if (!little_endian) {
                value = SwapBytes<REAL_T > (value);
            }
            out.write(reinterpret_cast<const char*> (&value), sizeof ( REAL_T));

            //actual value
            value = this->GetValue();

            if (!little_endian) {
                value = SwapBytes<REAL_T > (value);
            }
            out.write(reinterpret_cast<const char*> (&value), sizeof ( REAL_T));

            //derivative info
            uint32_t dsize = this->g.size();


            if (!little_endian) {
                dsize = SwapBytes<uint32_t > (dsize);
            }
            out.write(reinterpret_cast<const char*> (&dsize), sizeof (dsize));


            for (int i = 0; i < this->g.size(); i++) {

                bool b = this->g[i].first;
                value = this->g[i].second;
                if (b) {
                    out << '1';
                } else {
                    out << '0';
                }
                //
                //                if (!little_endian) {
                //                    id = SwapBytes<uint32_t > (id);
                //                }
                //
                //                out.write(reinterpret_cast<const char*> (&id), sizeof (id));
                //                
                if (!little_endian) {
                    value = SwapBytes<REAL_T > (value);
                }

                out.write(reinterpret_cast<const char*> (&value), sizeof ( REAL_T));

            }

            //statements
            uint32_t ssize = this->statements_m.size();

            if (!little_endian) {
                ssize = SwapBytes<uint32_t > (ssize);
            }

            out.write(reinterpret_cast<const char*> (&ssize), sizeof (ssize));

            uint32_t sid;
            uint32_t sop;
            REAL_T sval;

            for (int i = 0; i < this->statements_m.size(); i++) {
                sid = statements_m[i].id_m;
                if (!little_endian) {
                    sid = SwapBytes<uint32_t > (sid);
                }
                out.write(reinterpret_cast<const char*> (&sid), sizeof ( uint32_t));

                sop = statements_m[i].op_m;
                if (!little_endian) {
                    sop = SwapBytes<uint32_t > (sop);
                }
                out.write(reinterpret_cast<const char*> (&sop), sizeof ( uint32_t));

                sval = statements_m[i].value_m;
                if (!little_endian) {
                    sval = SwapBytes<REAL_T > (sval);
                }
                out.write(reinterpret_cast<const char*> (&sval), sizeof ( REAL_T));
            }





        }

        const Variable<REAL_T, group> Deserialize(std::istream &in) {


            Variable v;

            bool little_endian = true;

            int num = 1;
            if (*(char *) &num == 1) {
                little_endian = true;
            } else {
                little_endian = false;
            }

            if (!in.good()) {
                std::cout << "Error in Deserialize(std::istream &in), cannot read stream!";
                return v;
            }



            //read id
            uint32_t* id;
            char idc[sizeof (uint32_t) ];
            in.read(idc, sizeof (uint32_t));
            id = reinterpret_cast<uint32_t*> (idc);

            if (!little_endian) {
                v.id_m = SwapBytes<uint32_t > (*id);
            } else {
                v.id_m = *id;
            }

            //read name
            uint32_t* ns;
            uint32_t name_size;
            char nsc[sizeof (uint32_t) ];
            in.read(nsc, sizeof (uint32_t));
            ns = reinterpret_cast<uint32_t*> (nsc);

            if (!little_endian) {
                name_size = SwapBytes<uint32_t > (*ns);
            } else {
                name_size = *ns;
            }

            char namec[name_size + 1];
            namec[name_size] = '\0';
            in.read(namec, name_size);
            v.SetName(std::string(namec));




            //is independent
            char inde;
            in >> inde;

            if (inde == '1') {
                v.is_independent_m = true;
            } else {
                v.is_independent_m = false;
            }

            //is bounded
            char bounded;
            in >> bounded;

            if (bounded == '1') {
                v.bounded_m = true;
            } else {
                v.bounded_m = false;
            }


            //values min,max,actual
            REAL_T* value_ptr;
            REAL_T value;
            char valc[sizeof (REAL_T)];
            in.read(valc, sizeof (REAL_T));
            value_ptr = reinterpret_cast<REAL_T*> (valc);
            if (!little_endian) {
                value = SwapBytes<uint32_t > (*value_ptr);
            } else {
                value = *value_ptr;
            }
            v.min_boundary_m = value;

            in.read(valc, sizeof (REAL_T));
            value_ptr = reinterpret_cast<REAL_T*> (valc);

            if (!little_endian) {
                value = SwapBytes<uint32_t > (*value_ptr);
            } else {
                value = *value_ptr;
            }
            v.max_boundary_m = value;


            in.read(valc, sizeof (REAL_T));
            value_ptr = reinterpret_cast<REAL_T*> (valc);

            if (!little_endian) {
                value = SwapBytes<uint32_t > (*value_ptr);
            } else {
                value = *value_ptr;
            }
            v.value_m = value;


            //gradients
            uint32_t* gsize;
            uint32_t gs;
            char gsizec[sizeof (uint32_t) ];
            in.read(gsizec, sizeof (uint32_t));
            gsize = reinterpret_cast<uint32_t*> (gsizec);
            if (!little_endian) {
                gs = SwapBytes<uint32_t > (*gsize);
            } else {
                gs = *gsize;
            }



            v.g.resize(gs);


            for (uint32_t i = 0; i < gs; i++) {

                char c = '0';
                in >> c;

                bool b = false;
                if (c == '1') {
                    b = true;
                }


                in.read(valc, sizeof (REAL_T));
                value_ptr = reinterpret_cast<REAL_T*> (valc);
                if (!little_endian) {
                    value = SwapBytes<uint32_t > (*value_ptr);
                } else {
                    value = *value_ptr;
                }

                //                std::cout << v[] << "\n";
                v.g[i] = std::pair<bool, REAL_T > (b, value);

            }

            //statements

            uint32_t* ssize;
            uint32_t ss;
            char ssizec[sizeof (uint32_t) ];
            in.read(ssizec, sizeof (uint32_t));
            ssize = reinterpret_cast<uint32_t*> (ssizec);
            if (!little_endian) {
                ss = SwapBytes<uint32_t > (*ssize);
            } else {
                ss = *ssize;
            }

            v.statements_m.resize(ss);

            for (uint32_t i = 0; i < ss; i++) {


                uint32_t* sid_ptr;
                uint32_t sid;
                char sidc[sizeof (uint32_t) ];
                in.read(sidc, sizeof (uint32_t));
                sid_ptr = reinterpret_cast<uint32_t*> (sidc);
                if (!little_endian) {
                    sid = SwapBytes<uint32_t > (*sid_ptr);
                } else {
                    sid = *sid_ptr;
                }


                uint32_t* sop_ptr;
                uint32_t sop;
                char sopc[sizeof (uint32_t) ];
                in.read(sopc, sizeof (uint32_t));
                sop_ptr = reinterpret_cast<uint32_t*> (sopc);
                if (!little_endian) {
                    sop = SwapBytes<uint32_t > (*sop_ptr);
                } else {
                    sop = *sop_ptr;
                }


                REAL_T* value_ptr;
                REAL_T value;
                char valc[sizeof (REAL_T)];
                in.read(valc, sizeof (REAL_T));
                value_ptr = reinterpret_cast<REAL_T*> (valc);
                if (!little_endian) {
                    value = SwapBytes<uint32_t > (*value_ptr);
                } else {
                    value = *value_ptr;
                }
                //                Operation op = 
                v.statements_m[i] = Statement<REAL_T > (static_cast<Operation> (sop), value, sid);
            }



            return v;

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
            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();

#if defined(USE_HASH_TABLE)
            this->gradients_m.Clear();
#else
            //            this->gradients_m.clear();
            this->g.clear();
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
            //has_m.resize(IDGenerator::instance()->current() + 1);
            this->SetValue(other.GetValue());

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
                //                ind_iterator it;
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = other.Derivative(*it, found);
                //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                }
                // other.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                other.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it].first = true;
                    this->g[*it].second = other.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }


#endif

                if (Variable::is_supporting_arbitrary_order) {

                    std::vector<Statement<REAL_T> > temp_stmnt;
                    other.Push(temp_stmnt);
                    this->statements_m = temp_stmnt;
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

            //            this->id_m = expr.GetId();
            //has_m.resize(IDGenerator::instance()->current() + 1);
            if (Variable::is_recording_g) {
                //                ind_iterator it; // = this->gradients.lower_bound();
#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = expr.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                                    bool found = false;
                //                                    this->gradients_m[(*it)] = expr.Derivative(*it, found);
                //                                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                                }
#endif

//                                ids_m.clear();
                //                ids_m.set_empty_key(NULL);
                expr.PushIds(ids_m);
                // expr.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = false;
                    this->g[*it].first = true;
                    this->g[*it].second = expr.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }
                //                expr.GetIdRange(this->iv_min, this->iv_max);
                if (Variable::is_supporting_arbitrary_order) {
                    std::vector<Statement<REAL_T> > temp_stmnt;
                    expr.Push(temp_stmnt);
                    this->statements_m = temp_stmnt;
                }


            }
            this->SetValue(expr.GetValue());
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
            //has_m.resize(IDGenerator::instance()->current() + 1);
            if (Variable::is_recording_g) {

#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value + rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] + rhs.Derivative(*it, found);
                //                }
#endif

                // rhs.GetIdRange(this->iv_min, this->iv_max);

                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it].first = true;
                    this->g[*it].second = this->g[*it].second + rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }

                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }

            }
            this->value_m += rhs.GetValue();
            return *this;
        }

        Variable& operator+=(Variable& rhs) {
            if (Variable::is_recording_g) {
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value + rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] + rhs.Derivative(*it, found);
                //                }
#endif

                // rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it].first = true;
                    this->g[*it].second = this->g[*it].second + rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }

                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }

            }
            //            rhs.GetIdRange(this->iv_min, this->iv_max);
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
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert((*it))->value = (this->gradients_m.Lookup(*it)->value - rhs.Derivative(*it, found));
                }
#else
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] - rhs.Derivative(*it, found));
                //                }
#endif
                // rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it] = this->g[*it] - rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }
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
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                }
#endif
                //
                //                rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it] = this->g[*it] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }

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
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = (this->gradients_m.Lookup(*it)->value * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable::independent_variables_g.begin(); it != Variable::independent_variables_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                }
#endif
                //                rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_variables_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it] = (this->g[*it] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                    this->gradients_m[*it] = g[*it];
                }
                if (Variable::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (DIVIDE));
                }
            }
            this->value_m /= rhs.GetValue();
            return *this;
        }

        // And likewise for a Constant on the rhs

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

        const REAL_T Diff(const Variable &wrt) {

            if (this->statements_m.size() == 0) {
                return 0.0;
            }
            std::vector<std::pair<REAL_T, REAL_T> > v;
            v.reserve(this->statements_m.size());
            std::stack<std::pair<REAL_T, REAL_T>,
                    std::vector<std::pair<REAL_T, REAL_T> > > stack;
            //            ad::Stack<std::pair<T, T> > stack;
            bool found = false;
            int size = this->statements_m.size();
            std::pair<REAL_T, REAL_T> lhs = std::pair<REAL_T, REAL_T > (0, 0);
            std::pair<REAL_T, REAL_T> rhs = std::pair<REAL_T, REAL_T > (0, 0);

            //            Statement<REAL_T>* edata = (Statement<REAL_T>*)this->statements_m.data();
            //            std::cout << this->statements_m.size() << "\n\n";
            for (int i = 0; i < size; i++) {


                REAL_T temp = 0;



                switch (statements_m[i].op_m) {

                    case CONSTANT:
                        stack.push(std::pair<REAL_T, REAL_T > (statements_m[i].value_m, 0.0));
                        break;
                    case VARIABLE:
                        if (statements_m[i].id_m == wrt.GetId() && wrt.GetId() > 0) {
                            found = true;
                            //f(x) = x
                            //f'(x) = 1
                            stack.push(std::pair<REAL_T, REAL_T > (statements_m[i].value_m, 1.0));
                        } else {//constant
                            //f(x) = C
                            //f'(x) = 0
                            stack.push(std::pair<REAL_T, REAL_T > (statements_m[i].value_m, 0.0));
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
                        //                        std::cout<<lhs.first<<" "<<rhs.first<<"\n";
                        //                        temp = (lhs.second * rhs.first - lhs.first * rhs.second) / (rhs.first * rhs.first);
                        stack.push(std::pair<REAL_T, REAL_T > (lhs.first / rhs.first, 1));
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
                            temp = (lhs.second * 1.0 / std::pow((static_cast<REAL_T> (1.0) - std::pow(lhs.first, static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
                            stack.push(std::pair<REAL_T, REAL_T > (std::asin(lhs.first), temp));
                        } else {
                            stack.push(std::pair<REAL_T, REAL_T > (std::asin(lhs.first), 0));
                        }

                        break;
                    case ACOS:
                        lhs = stack.top();
                        stack.pop();
                        if (found) {
                            temp = (lhs.second * (-1.0) / std::pow((static_cast<REAL_T> (1.0) - std::pow(lhs.first, static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
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
                                    std::pow(lhs.first, (rhs.first - static_cast<REAL_T> (1.0)));
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




    };

    //list of independent variables by unique identifier.
    //    template<class REAL_T>
    //    std::set<uint32_t> Variable<REAL_T>::independent_variables_g;

    template<class REAL_T, int group>
    uint32_t Variable<REAL_T, group>::misses_g = 0;

    template<class REAL_T, int group>
    bool Variable<REAL_T, group>::is_recording_g = true;

    template<class REAL_T, int group>
    bool Variable<REAL_T, group>::is_supporting_arbitrary_order = false;

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

namespace array {







}

#endif	/* ETAD_HPP */

