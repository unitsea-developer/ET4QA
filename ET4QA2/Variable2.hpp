/* 
 * File:   Variable2.hpp
 * Author: matthewsupernaw
 *
 * Created on June 5, 2014, 2:25 PM
 */

#ifndef VARIABLE2_HPP
#define	VARIABLE2_HPP

#include "ET4AD.hpp"


namespace ad {

    template<class REAL_T, int group>
    class Variable2 : public ExpressionBase<REAL_T, Variable2<REAL_T, group> > {
        VariableStorage<REAL_T>* storage;
        REAL_T value_m;

        std::string name_m;
        bool bounded_m;
        REAL_T min_boundary_m;
        REAL_T max_boundary_m;
        bool is_independent_m;
        uint32_t iv_id_m; //id when is a independent Variable2.

        //        uint32_t iv_min; //if this Variable2 is caching derivatives, this is the min independent Variable2.
        //        uint32_t iv_max; //if this Variable2 is caching derivatives, this is the max independent Variable2.


        //        static std::set<uint32_t> independent_Variable2s_g;
        static bool is_recording_g;

        static bool is_supporting_arbitrary_order;


        typedef std::vector<std::pair<bool, REAL_T> > GradientVector;
        GradientVector g;
        //        GradientMap gradients_m;

        typedef IdsSet indepedndent_Variable2s;
        typedef IdsSet::iterator indepedndent_Variable2s_iterator;
        typedef IdsSet::const_iterator const_indepedndent_Variable2s_iterator;
        indepedndent_Variable2s ids_m;
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
        Variable2() : ExpressionBase<REAL_T, Variable2<REAL_T> >(0),
        value_m(0.0),
        bounded_m(false),
        is_independent_m(false),
        storage(new DefaultStorage<REAL_T>()){
//            this->storage->SetValue(0);
        }

        /**
         * Conclass a Variable2 with a initial value. Optional parameter to make
         * it a independent Variable2.
         * 
         * @param value
         * @param is_independent
         */
        Variable2(const REAL_T& value, bool is_independent = false) : storage(new DefaultStorage<REAL_T>()), value_m(value), is_independent_m(is_independent), iv_id_m(0) {
            //this->ids_m.set_empty_key(NULL);
            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();
            if (is_independent) {
                //                this->id_m = IDGenerator::instance()->next();
                this->SetAsIndependent(is_independent);
                this->bounded_m = false;
                this->min_boundary_m = std::numeric_limits<REAL_T>::min();
                this->max_boundary_m = std::numeric_limits<REAL_T>::max();
                //                Variable2::independent_Variable2s_g.insert(this->GetId());
                if (Variable2::IsSupportingArbitraryOrder()) {
                    this->statements_m.push_back(Statement<REAL_T > (Variable2, this->GetValue(), this->GetId()));
                }
            }

        }

        /**
         * Copy conclassor. 
         * 
         * @param rhs
         */
        Variable2(const Variable2& orig) : storage(new DefaultStorage<REAL_T>()) {
            orig.PushStorage(this->storage);
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
            if (Variable2::IsSupportingArbitraryOrder()) {
                orig.Push(this->statements_m);
            }

        }

        /**
         * Conclasss a Variable2 from expression expr.
         * @param rhs
         */
        template<class T>
        Variable2(const ExpressionBase<REAL_T, T>& expr) : storage(new DefaultStorage<REAL_T>()) {
            expr.PushStorage(this->storage);
            //has_m.resize(IDGenerator::instance()->current() + 1);

            this->bounded_m = false;
            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();
            if (Variable2::is_recording_g) {
                //                this->id_m = expr.GetId();
                //                ind_iterator it;
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
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
                indepedndent_Variable2s_iterator it;
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = false;
                    this->g[*it].first = true;
                    this->g[*it].second = expr.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }
                if (Variable2::IsSupportingArbitraryOrder()) {
                    expr.Push(this->statements_m);
                }
            }
            value_m = expr.GetValue();

        }

        ~Variable2() {
            delete storage;
#if !defined(USE_HASH_TABLE)
            //            if (this->is_independent_m) {
            //                this->independent_Variable2s_g.erase(this->GetId());
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
            Variable2::is_recording_g = is_recording;
        }

        /**
         * If true, derivatives are calculated and stored. If false,
         * expressions are evaluated for value only.
         * 
         * @return 
         */
        static bool IsRecording() {
            return Variable2::is_recording_g;
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
            Variable2::is_supporting_arbitrary_order = support_arbitrary_order;
        }

        /*
         * If true, expressions are recorded into a vector of statements that 
         * can be used later to compute derivatives of arbitrary order, or
         * be used to build expression strings.
         * 
         */
        static bool IsSupportingArbitraryOrder() {
            return Variable2::is_supporting_arbitrary_order;
        }

        /**
         * Returns the value of this Variable2.
         * 
         * @return 
         */
        inline const REAL_T GetValue() const {
            return this->value_m;
        }

        /**
         * Sets the value of this Variable2. If the Variable2 is bounded, 
         * the value will be set between the min and max boundary. If the
         * value is less than the minimum boundary, the Variable2s value is 
         * set to minimum boundary. If the  value is greater than the maximum
         * boundary, the Variable2s value is set to maximum boundary. If the 
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
         * Returns the name of this Variable2. Names are not initialized and 
         * if it is not set this function will return a empty string.
         * @return 
         */
        std::string GetName() const {
            return name_m;
        }

        /**
         * Set the name of this Variable2.
         * 
         * @param name
         */
        void SetName(std::string name) {
            this->name_m = name;
        }

        /**
         * Returns true if this Variable2 is bounded. Otherwise,
         * false.
         * 
         * @return 
         */
        bool IsBounded() {
            return this->bounded_m;
        }

        /**
         * Set the min and max boundary for this Variable2. 
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
         * Return the minimum boundary for this Variable2. The default value 
         * is the results of std::numeric_limits<REAL_t>::min().
         * 
         * @return 
         */
        const REAL_T GetMinBoundary() {
            return this->min_boundary_m;
        }

        /**
         * Return the maximum boundary for this Variable2. The default value 
         * is the results of std::numeric_limits<REAL_t>::max().
         * 
         * @return 
         */
        const REAL_T GetMaxBoundary() {
            return this->max_boundary_m;
        }

        /**
         * Make this Variable2 an independent Variable2. If set to true,
         * the Variable2s unique identifier is registered in the static set
         * of independent Variable2s. During function evaluations, gradient 
         * information is accumulated in post-order wrt to Variable2s in the 
         * static set of independent Variable2s. If set false, this Variable2s 
         * unique identifier will be removed from the set, if it existed.
         * 
         * 
         * To access the derivative wrt to a Variable2, call function:
         * 
         * const REAL_t wrt(const Variable2 & ind)
         * 
         * 
         * @param is_independent
         */
        void SetAsIndependent(const bool &is_independent) {
            if (this->iv_id_m == 0) {
                this->iv_id_m = IDGenerator::instance()->next();
                if (Variable2<REAL_T>::IsSupportingArbitraryOrder()) {
                    this->statements_m.clear();
                    this->statements_m.push_back(Statement<REAL_T > (Variable2, this->GetValue(), this->GetId()));
                }
                //                this->id_m = iv_id_m;
                //                if (this->iv_min == std::numeric_limits<uint32_t>::max() && this->iv_max == std::numeric_limits<uint32_t>::min()) {
                //                    this->iv_max = this->GetId();
                //                    this->iv_min = this->GetId();
                //                }
            }
            if (this->is_independent_m && !is_independent) {
                //                Variable2::independent_Variable2s_g.erase(this->GetId());
                this->id_m = 0;
                this->is_independent_m = false;
            }

            if (!this->is_independent_m && is_independent) {
                //                Variable2::independent_Variable2s_g.insert(this->GetId());
                this->id_m = this->iv_id_m;
                this->is_independent_m = true;
            }
        }

        bool IsIndependent() {
            return this->is_independent_m;
        }

        /**
         * Returns the derivative with respect to a Variable2s who's
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
                //                    Variable2<REAL_T>::misses_g++;
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
            //                statements.push_back(Statement<REAL_T > (Variable2, this->GetValue(), this->GetId()));
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

        inline void PushStorage(Variable2Storage<REAL_T> * ids) const {

            if (Variable2<REAL_T, group>::IsRecording()) {

                if (this->GetId() != 0) {
                    ids->AddId(this->GetId());
                } else {
                    ids->Merge(this->storage);
                }


                for (int i = 0; i < this->storage->DerivativeInfoSize(); i++) {
                    ids->SetDerivative(this->storage->GetDerivative(i), i);
                }

                if (Variable2<REAL_T, group>::IsSupportingArbitraryOrder()) {
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
         * w.r.t a Variable2. 
         * 
         * @param ind
         * @return 
         */
        const REAL_T WRT(const Variable2 & ind) {

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
         * w.r.t a Variable2. 
         * @param ind
         * @return 
         */
        const REAL_T WRT(const Variable2 & ind) const {

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

        const Variable2<REAL_T, group> Deserialize(std::istream &in) {


            Variable2 v;

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
         * Set the Variable2 equal to the real type rhs.
         * 
         * The gradient map will be cleared.
         * 
         * @param rhs
         * @return 
         */
        Variable2& operator=(const REAL_T& rhs) {
            this->SetValue(rhs);
            //            iv_min = std::numeric_limits<uint32_t>::max();
            //            iv_max = std::numeric_limits<uint32_t>::min();

#if defined(USE_HASH_TABLE)
            this->gradients_m.Clear();
#else
            //            this->gradients_m.clear();
            this->g.clear();
            if (Variable2::IsRecording()) {
                if (Variable2::is_supporting_arbitrary_order) {
                    this->statements_m.clear();
                    this->statements_m.push_back(Statement<REAL_T > (Variable2, this->GetValue(), this->GetId()));
                }
            }
#endif

            return *this;
        }

        /**
         * Sets the value of this Variable2 to that of rhs. Derivatives 
         * are stored in the encapsulated gradient map.
         * 
         * @param rhs
         * @return 
         */
        Variable2& operator=(const Variable2& other) {
            //                        this->id = other.getId();
            //has_m.resize(IDGenerator::instance()->current() + 1);
            this->SetValue(other.GetValue());
            this->storage->Reset();
            other.PushStorage(this->storage);
            //            this->gradients.clear();
            if (ad::Variable2<REAL_T>::is_recording_g) {

#if defined(USE_HASH_TABLE)
                ind_iterator it;
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = other.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                this->gradients.insert(rhs.gradients.begin(), rhs.gradients.end());
                //                ind_iterator it;
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = other.Derivative(*it, found);
                //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                }
                // other.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
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

                if (Variable2::is_supporting_arbitrary_order) {

                    std::vector<Statement<REAL_T> > temp_stmnt;
                    other.Push(temp_stmnt);
                    this->statements_m = temp_stmnt;
                }

            }
            this->is_independent_m = other.is_independent_m;
            return *this;
        }

        /**
         * Set the Variable2s value to the result of the expression rhs. 
         * Derivatives are calculated and stored in the encapsulated 
         * gradient map.
         * @param rhs
         * @return 
         */
        template<class T>
        Variable2& operator=(const ExpressionBase<REAL_T, T>& expr) {
            expr.PushStorage(this->storage);
            //            this->id_m = expr.GetId();
            //has_m.resize(IDGenerator::instance()->current() + 1);
            if (Variable2::is_recording_g) {
                //                ind_iterator it; // = this->gradients.lower_bound();
#if defined(USE_HASH_TABLE)
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = expr.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                                    bool found = false;
                //                                    this->gradients_m[(*it)] = expr.Derivative(*it, found);
                //                                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                                }
#endif

                //                ids_m.clear();
                //                ids_m.set_empty_key(NULL);
                expr.PushIds(ids_m);
                // expr.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = false;
                    this->g[*it].first = true;
                    this->g[*it].second = expr.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }
                //                expr.GetIdRange(this->iv_min, this->iv_max);
                if (Variable2::is_supporting_arbitrary_order) {
                    std::vector<Statement<REAL_T> > temp_stmnt;
                    expr.Push(temp_stmnt);
                    this->statements_m = temp_stmnt;
                }


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
        Variable2& operator+=(const ExpressionBase<REAL_T, T>& rhs) {
            //            return *this = (*this +rhs);
            //has_m.resize(IDGenerator::instance()->current() + 1);
            if (Variable2::is_recording_g) {

#if defined(USE_HASH_TABLE)
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value + rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] + rhs.Derivative(*it, found);
                //                }
#endif

                // rhs.GetIdRange(this->iv_min, this->iv_max);

                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it].first = true;
                    this->g[*it].second = this->g[*it].second + rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }

                if (Variable2::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }

            }
            this->value_m += rhs.GetValue();
            return *this;
        }

        Variable2& operator+=(Variable2& rhs) {
            if (Variable2::is_recording_g) {
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value + rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] + rhs.Derivative(*it, found);
                //                }
#endif

                // rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it].first = true;
                    this->g[*it].second = this->g[*it].second + rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }

                if (Variable2::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }

            }
            //            rhs.GetIdRange(this->iv_min, this->iv_max);
            this->value_m += rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable2& operator-=(const ExpressionBase<REAL_T, T>& rhs) {
            return *this = (*this -rhs);
            //            if (Variable2::is_recording_g) {
            //                ind_iterator it;
            //#if defined(USE_HASH_TABLE)
            //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m.Insert((*it))->value = (this->gradients_m.Lookup(*it)->value - rhs.Derivative(*it, found));
            //                }
            //#else
            //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] - rhs.Derivative(*it, found));
            //                }
            //#endif
            //                if (Variable2::is_supporting_arbitrary_order) {
            //                    rhs.Push(this->statements_m);
            //                    this->statements_m.push_back(Statement(MINUS));
            //                }
            //            }
            //            this->value_m -= rhs.GetValue();
            //            return *this;
        }

        Variable2& operator-=(Variable2& rhs) {
            if (Variable2::is_recording_g) {
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert((*it))->value = (this->gradients_m.Lookup(*it)->value - rhs.Derivative(*it, found));
                }
#else
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] - rhs.Derivative(*it, found));
                //                }
#endif
                // rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it] = this->g[*it] - rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }
                if (Variable2::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (MINUS));
                }
            }
            this->value_m -= rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable2& operator*=(const ExpressionBase<REAL_T, T>& rhs) {
            return *this = (*this * rhs);
            //            if (Variable2::is_recording_g) {
            //                ind_iterator it;
            //#if defined(USE_HASH_TABLE)
            //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#else
            //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#endif
            //                if (Variable2::is_supporting_arbitrary_order) {
            //                    rhs.Push(this->statements_m);
            //                    this->statements_m.push_back(Statement(MULTIPLY));
            //                }
            //            }
            //            this->value_m *= rhs.GetValue();
            return *this;
        }

        Variable2& operator*=(Variable2& rhs) {
            if (Variable2::is_recording_g) {
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = this->gradients_m.Lookup(*it)->value * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = this->gradients_m[(*it)] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                }
#endif
                //
                //                rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it] = this->g[*it] * rhs.GetValue() + this->GetValue() * rhs.Derivative(*it, found);
                    //                    this->gradients_m[*it] = g[*it];
                }

                if (Variable2::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (MULTIPLY));
                }
            }
            this->value_m *= rhs.GetValue();
            return *this;
        }

        template<class T>
        Variable2& operator/=(const ExpressionBase<REAL_T, T>& rhs) {
            return *this = (*this / rhs);
            //            if (Variable2::is_recording_g) {
            //                ind_iterator it;
            //#if defined(USE_HASH_TABLE)
            //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m.Insert(*it)->value = (this->gradients_m.Lookup(*it)->value * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#else
            //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
            //                    bool found = false;
            //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
            //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
            //                }
            //#endif
            //                if (Variable2::is_supporting_arbitrary_order) {
            //                    rhs.Push(this->statements_m);
            //                    this->statements_m.push_back(Statement(DIVIDE));
            //                }
            //            }
            //            this->value_m /= rhs.GetValue();
            return *this;
        }

        Variable2& operator/=(Variable2& rhs) {
            if (Variable2::is_recording_g) {
                //has_m.resize(IDGenerator::instance()->current() + 1);

#if defined(USE_HASH_TABLE)
                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                    bool found = false;
                    this->gradients_m.Insert(*it)->value = (this->gradients_m.Lookup(*it)->value * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                }
#else
                //                for (it = Variable2::independent_Variable2s_g.begin(); it != Variable2::independent_Variable2s_g.end(); ++it) {
                //                    bool found = false;
                //                    this->gradients_m[(*it)] = (this->gradients_m[(*it)] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                //                    //                std::cout<<"diff wrt "<<(*it)<<" = " <<rhs.Derivative(*it, found);
                //                }
#endif
                //                rhs.GetIdRange(this->iv_min, this->iv_max);
                //                for (uint32_t i = this->iv_min; i < this->iv_max + 1; i++) {
                indepedndent_Variable2s_iterator it;
                rhs.PushIds(ids_m);
                g.resize(IDGenerator::instance()->current() + 1);
                //                for (uint32_t i = this->iv_min; i< this->iv_max + 1; i++) {
                for (it = this->ids_m.begin(); it != ids_m.end(); ++it) {
                    bool found = true;
                    this->g[*it] = (this->g[*it] * rhs.GetValue() - this->GetValue() * rhs.Derivative(*it, found)) / (rhs.GetValue() * rhs.GetValue());
                    //                    this->gradients_m[*it] = g[*it];
                }
                if (Variable2::is_supporting_arbitrary_order) {
                    rhs.Push(this->statements_m);
                    this->statements_m.push_back(Statement<REAL_T > (DIVIDE));
                }
            }
            this->value_m /= rhs.GetValue();
            return *this;
        }

        // And likewise for a Constant on the rhs

        Variable2& operator+=(const REAL_T& rhs) {
            value_m += rhs;
            if (Variable2::IsRecording()) {
                if (Variable2::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (PLUS));
                }
            }
            return *this;
        }

        Variable2& operator-=(const REAL_T& rhs) {

            if (Variable2::IsRecording()) {
                if (Variable2::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (MINUS));
                }
            }
            value_m -= rhs;
            return *this;
        }

        Variable2& operator*=(const REAL_T& rhs) {
            if (Variable2::IsRecording()) {
                if (Variable2::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (MULTIPLY));
                }
            }
            return *this = (*this * rhs);
        }

        Variable2& operator/=(const REAL_T& rhs) {
            if (Variable2::IsRecording()) {
                if (Variable2::is_supporting_arbitrary_order) {
                    this->statements_m.push_back(Statement<REAL_T > (CONSTANT, rhs));
                    this->statements_m.push_back(Statement<REAL_T > (DIVIDE));
                }
            }
            return *this = (*this / rhs);
        }

        const REAL_T Diff(const Variable2 &wrt) {

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
                    case Variable2:
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

    //list of independent Variable2s by unique identifier.
    //    template<class REAL_T>
    //    std::set<uint32_t> Variable2<REAL_T>::independent_Variable2s_g;

    template<class REAL_T, int group>
    uint32_t Variable2<REAL_T, group>::misses_g = 0;

    template<class REAL_T, int group>
    bool Variable2<REAL_T, group>::is_recording_g = true;

    template<class REAL_T, int group>
    bool Variable2<REAL_T, group>::is_supporting_arbitrary_order = false;
}




#endif	/* VARIABLE2_HPP */

