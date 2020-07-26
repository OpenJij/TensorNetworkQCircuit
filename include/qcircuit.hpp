//Copyright (c) 2020 Jij Inc.


#pragma once

#include <itensor/all.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <iostream>
#include <array>
#include "circuit_topology.hpp"
#include "quantum_gate.hpp"
#include "qcircuit_exception.hpp"

namespace qcircuit {
    using namespace itensor;

    class QCircuit;
    Cplx overlap(QCircuit circuit1, const std::vector<ITensor>& op, QCircuit circuit2,
                 const Args& args = Args::global());

    /**
     * @brief Class to store and modify wave function.
     */
    class QCircuit {
    private:
        std::vector<Index> a;   //!< @brief Link (bond) indices.
        std::vector<Index> s;   //!< @brief Physical (on-site) indices.
        std::vector<ITensor> M; //!< @brief Tensor entities.
        ITensor Psi;            //!< @brief TPS wave function.

        const CircuitTopology topology; //!< @brief Circuit topology.

        std::pair<std::size_t, std::size_t> cursor; //!< @brief Position of cursor, which must be laid across two neighboring sites.

        std::mt19937 random_engine;

        Args default_args = Args(); //!< @brief Default arguments for ITensor functions.

    public:
        /** @brief Constructor to initialize TPS wave function.
         *
         * Cursor position will be set at (0, 1).
         *
         * @param topology Circuit Topology.
         * @param init_qubits Initial qubit states for each site.
         * @param physical_indices Physical (on-site) indices.
         * If not specified, the indices will be initialized with new IDs.
         * This argument is mainly used to share physical indices among
         * some "replica" wave functions in the same circuit.
         **/
        QCircuit(const CircuitTopology& topology,
                 const std::vector<std::pair<std::complex<double>, std::complex<double>>>& init_qubits,
                 const std::vector<Index>& physical_indices = std::vector<Index>()) :
            a(), s(physical_indices), topology(topology), random_engine(std::random_device()()) {

            if(!topology.isConnectedGraph()) {
                throw QCircuitException("Invalid circuit topology : Some nodes are unreachable");
            }

            /* Initialize link indices */
            a.reserve(topology.numberOfLinks());
            for(size_t i = 0; i < topology.numberOfLinks();i++) {
                a.emplace_back(1, "LinkInd");
            }

            /* Initialize physical indices */
            if(s.empty()) {
                s.reserve(this->size());
                for(size_t i = 0;i < this->size();i++) {
                    s.emplace_back(2, "SiteInd");
                }
            }

            /* Initialize tensor entities */
            M.reserve(this->size());
            for(size_t i = 0;i < this->size();i++) {
                auto neighbors = topology.neighborsOf(i);

                std::vector<Index> ind_list;
                ind_list.reserve(neighbors.size()+1);
                std::vector<IndexVal> ind_val_list;
                ind_val_list.reserve(neighbors.size()+1);

                //fill ind_list
                ind_list.push_back(s[i]);
                for(const auto& neighbor : neighbors){
                    ind_list.push_back(a[neighbor.link]);
                }

                //fill ind_val_list
                ind_val_list.push_back(s[i](1));
                for(const auto& neighbor : neighbors){
                    ind_val_list.push_back(a[neighbor.link](1));
                }

                //insert into M
                M.emplace_back(ind_list);
                ind_val_list[0] = s[i](1);
                M[i].set(ind_val_list, init_qubits[i].first);
                ind_val_list[0] = s[i](2);
                M[i].set(ind_val_list, init_qubits[i].second);
            }

            /* set cursor position */
            {
                cursor.first = 0;

                /* find the minimum numbered bit of neighbors of the bit 0. */
                size_t index = this->size();
                for(const auto& neighbor : topology.neighborsOf(cursor.first)) {
                    if(neighbor.site < index) {
                        index = neighbor.site;
                    }
                }
                cursor.second = index;
            }

            Psi = M[cursor.first]*M[cursor.second];
        }

        /** @brief Constructor to initialize TPS wave function with |000 ... 000>.
         *
         * Cursor position will be set at (0, 1).
         *
         * @param topology Circuit Topology.
         * @param physical_indices Physical (on-site) indices.
         * If not specified, the indices will be initialized with new IDs.
         * This argument is mainly used to share physical indices among
         * some "replica" wave functions in the same circuit.
         **/
        QCircuit(const CircuitTopology& topology,
                 const std::vector<Index>& physical_indices = std::vector<Index>()) :
            QCircuit(topology,
                     std::vector<std::pair<std::complex<double>, std::complex<double>>>(topology.numberOfBits(), std::make_pair(1.0, 0.0)),
                     physical_indices) {}

        /** @brief returns number of qubits */
        size_t size() const {
            return this->topology.numberOfBits();
        }

        /** @brief decompose and truncate the system wave-function at cursor position.  */
        void decomposePsi(const Args& args) {
            ITensor U,S,V;

            // fetch nodelist in index "first"
            auto list = topology.neighborsOf(cursor.first);

            /* get link index corresponding to cursor position
             * (this can be combined with the std::remove_if below)
             */
            auto second_itr = std::find_if(list.begin(), list.end(),
                                              [&](auto x) {return x.site == cursor.second;});
            auto link_ind = (*second_itr).link;

            // remove the element which connects to node second
            auto pend = std::remove_if(list.begin(), list.end(), [&](const auto& t){return t.site == cursor.second;});
            std::vector<CircuitTopology::Neighbor> newlist;
            for(auto it = list.begin(); it != pend; ++it){
                newlist.push_back(*it);
            }

            // set U
            std::vector<Index> ind_list;
            ind_list.reserve(newlist.size()+1);
            ind_list.push_back(s[cursor.first]);
            for(const auto& elem : newlist){
                ind_list.push_back(a[elem.link]);
            }

            U = ITensor(ind_list);

            Spectrum spec = svd(Psi, U, S, V, args);
            a[link_ind] = commonIndex(U,S);
            S /= norm(S); //normalize
            M[cursor.first] = U;
            M[cursor.second] = S*V;
        }

        void decomposePsi() {
            decomposePsi(default_args);
        }

        /** @brief shifts cursor position to specified neighboring site `dest`.
         *
         * `direction` specifies which site of the cursor position to be used as the "head" of move.
         * `direction == 0` means not specifying the direction.
         */
        Spectrum shiftCursorTo(size_t dest, int direction, const Args& args) {
            assert(dest != cursor.first);
            assert(dest != cursor.second);
            assert(direction == 0 || direction == 1 || direction == 2);

            if(direction == 0) {
                //first
                for(auto&& i : topology.neighborsOf(cursor.first)) {
                    if(i.site == dest) {
                        direction = 1;
                        break;
                    }
                }
                //second
                for(auto&& i : topology.neighborsOf(cursor.second)) {
                    if(i.site == dest) {
                        direction = 2;
                        break;
                    }
                }
            }

            Spectrum spec;

            if(direction == 1) {
                ITensor U,S,V;

                //fetch nodelist in index "second"
                auto list = this->topology.neighborsOf(cursor.second);
                //get link index
                size_t link_ind = 0;
                for(auto&& i : list){
                    if(i.site == cursor.first){
                        link_ind = i.link;
                    }
                }
                //remove the element which connects to node first
                auto pend = std::remove_if(list.begin(), list.end(), [&](const auto& t){return t.site == cursor.first;});
                std::vector<CircuitTopology::Neighbor> newlist;
                for(auto it = list.begin(); it != pend; ++it){
                    newlist.push_back(*it);
                }

                //set U
                std::vector<Index> ind_list;
                ind_list.reserve(newlist.size()+1);
                ind_list.push_back(s[cursor.first]);
                for(const auto& elem : newlist){
                    ind_list.push_back(a[elem.link]);
                }

                U = ITensor(ind_list);

                spec = svd(Psi, U, S, V, args);

                a[link_ind] = commonIndex(S,V);
                S /= norm(S); //normalize
                M[cursor.second] = V;
                Psi = M[dest]*U*S;

                cursor.second = cursor.first;
                cursor.first = dest;
            } else if(direction == 2) {
                ITensor U,S,V;

                //fetch nodelist in index "first"
                auto list = topology.neighborsOf(cursor.first);
                //get link index
                size_t link_ind = 0;
                for(auto&& i : list) {
                    if(i.site == cursor.second){
                        link_ind = i.link;
                    }
                }
                //remove the element which connects to node second
                auto pend = std::remove_if(list.begin(), list.end(), [&](const auto& t){return t.site == cursor.second;});
                std::vector<CircuitTopology::Neighbor> newlist;
                for(auto it = list.begin(); it != pend; ++it){
                    newlist.push_back(*it);
                }

                //set U
                std::vector<Index> ind_list;
                ind_list.reserve(newlist.size()+1);
                ind_list.push_back(s[cursor.first]);
                for(const auto& elem : newlist){
                    ind_list.push_back(a[elem.link]);
                }

                U = ITensor(ind_list);

                spec = svd(Psi, U, S, V, args);

                a[link_ind] = commonIndex(U,S);
                S /= norm(S); //normalize
                M[cursor.first] = U;
                Psi = S*V*M[dest];

                cursor.first = cursor.second;
                cursor.second = dest;
            } else {
                assert(false && "cannot move to this direction");
            }

            return spec;
        }

        /** @brief moves cursor along give path.
         *
         * If the cursor is positioned at (0, 1) and calling
         * `moveCursorAlong({2, 3, 4, 5})`,
         * the cursor will be moved to (4, 5) along the path.
         */
        void moveCursorAlong(const std::vector<size_t>& path, const Args& args) {
            for(auto site : path) {
                shiftCursorTo(site, 0, args);
            }
        }

        void moveCursorAlong(const std::vector<size_t>& path) {
            moveCursorAlong(path, default_args);
        }

        /** @brief moves cursor to given destination along the shortest path. */
        void moveCursorTo(size_t destination1, size_t destination2, const Args& args) {
            if(!topology.hasLinkBetween(destination1, destination2)) {
                std::stringstream ss;
                ss << "There is no link between " << destination1 << " and " << destination2;
                throw QCircuitException(ss.str());
            }

            if( (destination1 == cursor.first || destination1 == cursor.second) &&
                (destination2 == cursor.first || destination2 == cursor.second) ) {
                return; // No need to move
            }

            auto path = topology.getRoute(cursor, std::make_pair(destination1, destination2));
            moveCursorAlong(path, args);


            /* Here, only one of the destination sites is covered by the cursor.
             * So move the cursor to the remaining destination.
             */
            int direction = (cursor.first == destination1 || cursor.first == destination2) ?
                1 : 2;
            int dest = (cursor.first == destination1 || cursor.second == destination1) ?
                destination2 : destination1;
            shiftCursorTo(dest ,direction, args);
        }

        void moveCursorTo(size_t destination1, size_t destination2) {
            moveCursorTo(destination1, destination2, default_args);
        }

        /** @brief applies `op` at cursor position. */
        void applyAtCursor(const ITensor& op) {
            assert(op.inds().size() == 4);
            for(auto&& elem : op.inds()){
                //check if
                assert(elem == s[cursor.first] || elem == s[cursor.second] || elem == prime(s[cursor.first]) || elem == prime(s[cursor.second]));
            }

            this->Psi = op * prime(Psi, s[cursor.first], s[cursor.second]);
        }

        /**
         * @brief applies one-site gates `gate1` and `gate2` onto the gate position.
         *
         * Cursor position will be automatically moved.
         */
        void apply(const OneSiteGate& gate1,
                   const OneSiteGate& gate2,
                   const Args& args) {
            moveCursorTo(gate1.site, gate2.site, args);
            auto op = generateTensorOp(gate1) * generateTensorOp(gate2);
            applyAtCursor(op);
        }

        void apply(const OneSiteGate& gate1,
                   const OneSiteGate& gate2) {
            apply(gate1, gate2, default_args);
        }

        /**
         * @brief applies one-site gates `gate1` onto the gate position and
         * identity gate onto a dummy site.
         *
         * Cursor position will be automatically moved.
         */
        void apply(const OneSiteGate& gate1, const Args& args) {
            Id gate2(topology.neighborsOf(gate1.site)[0].site);

            this->apply(gate1, gate2, args);
        }

        void apply(const OneSiteGate& gate1) {
            apply(gate1, default_args);
        }

        /**
         * @brief applies two-site gate `gate` onto the gate position.
         *
         * Cursor position will be automatically moved.
         */
        void apply(const TwoSiteGate& gate, const Args& args) {
            moveCursorTo(gate.site1, gate.site2, args);
            auto op = generateTensorOp(gate);
            applyAtCursor(op);
        }

        void apply(const TwoSiteGate& gate) {
            apply(gate, default_args);
        }

        /** @brief returns the tensor operator corresponding to `gate`. */
        ITensor generateTensorOp(const Gate& gate) const {
            return gate.op(s);
        }

        /**
          * @brief returns probability the qubit at `site` to be observed as value "0",
          * i.e. `<psi|Proj_0(i)|psi>` as "Born rule".
          */
        double probabilityOfZero(size_t site) const {
            std::vector<ITensor> op;
            op.reserve(this->size());
            for(size_t i = 0;i < this->size();i++) {
                if(i == site) {
                    op.push_back(generateTensorOp(Proj_0(i)));
                } else {
                    op.push_back(generateTensorOp(Id(i)));
                }
            }

            /*
             * <psi|Proj_0(i)|psi> = |w_{i,0}|^2 is real by definition,
             * where w_{i,0} is the coefficient corresponding to the base |....0....>.
             *                                                     (ith qubit) ^
             */
            return std::real(overlap(*this, op, *this));
        }

        /** @brief observes the qubit state at `site` and returns the projected qubit value (0 or 1). */
        int observeQubit(size_t site, const Args& args) {
            auto prob0 = probabilityOfZero(site);

            std::uniform_real_distribution<> dist(0.0, 1.0);
            int state = (dist(random_engine) < prob0) ? 0 : 1; // measurement

            size_t neighbor = topology.neighborsOf(site)[0].site; // Dummy site to which Id operator is applied.
            if(state == 0) {
                apply(Proj_0(site), Id(neighbor), args);
            } else {
                apply(Proj_1(site), Id(neighbor), args);
            }
            this->normalize();

            return state;
        }

        int observeQubit(size_t site) {
            return observeQubit(site, default_args);
        }

        /** @brief reset the qubit state at `site` to |0>. */
        void resetQubit(size_t site, const Args& args) {
            auto prob0 = probabilityOfZero(site);

            size_t neighbor = topology.neighborsOf(site)[0].site; // Dummy site to which Id operator is applied.
            if(prob0 > 0.0) {
                apply(Proj_0(site), Id(neighbor), args);
            } else {
                apply(Proj_1(site), Id(neighbor), args);
            }
            this->normalize();
        }

        void resetQubit(size_t site) {
            resetQubit(site, default_args);
        }

        /** @brief sets cutoff. */
        QCircuit& setCutoff(double cutoff) {
            default_args.add("Cutoff", cutoff);
            /*
             * According to ITensor implementation,
             * existing value will be just replaced.
             */

            return *this; // for method chaining
        }

        /** @brief returns cutoff. If not being set, returns 0.  */
        double getCutoff() const {
            return default_args.getReal("Cutoff", 0);
        }

        /** @brief sets cutoff. */
        QCircuit& setMaxDim(int max_dim) {
            default_args.add("MaxDim", max_dim);
            /*
             * According to ITensor implementation,
             * existing value will be just replaced.
             */

            return *this; // for method chaining
        }

        /** @brief returns maximum bond dimension to be kept.
            If not being set, returns 0.  */
        int getMaxDim() const {
            return default_args.getInt("MaxDim", 0);
        }

        void normalize() {
            Psi /= norm(Psi);
        }

        void primeAll() {
            for(auto& elem : this->s){
                elem = prime(elem);
            }

            for(auto& elem : this->a){
                elem = prime(elem);
            }

            for(auto& elem : this->M){
                elem = prime(elem);
            }

            Psi = prime(Psi);
        }

        const ITensor& Mref(size_t i) const {
            assert(0 <= i && i < this->size());
            return this->M[i];
        }

        const std::vector<ITensor>& Mref() const {
            return this->M;
        }

        const ITensor& Psiref() const {
            return this->Psi;
        }

        const Index& site(size_t i) const {
            assert(0 <= i && i < this->size());
            return this->s[i];
        }

        const std::vector<Index>& site() const {
            return this->s;
        }

        const std::pair<size_t, size_t>& getCursor() const {
            return this->cursor;
        }

        void PrintMat() const {
            for(auto&& elem : M){
                PrintData(elem);
            }
            std::cout << "-----------" << std::endl;
            PrintData(Psi);
        }

        void PrintCursor() const {
            std::cout << "(" << cursor.first << "," << cursor.second << ")" << std::endl;
        }
    };

    inline Cplx overlap(QCircuit circuit1,
                        const std::vector<ITensor>& op,
                        QCircuit circuit2,
                        const Args& args) {
        /* Prototype declaration of this function is placed at the top. */
        assert(op.size() == circuit1.size() && op.size() == circuit2.size());

        circuit1.decomposePsi(args);
        circuit2.decomposePsi(args);

        circuit2.primeAll();

        ITensor ret_t = dag(circuit1.Mref(0))*op[0]*circuit2.Mref(0);
        for(size_t i = 1; i < circuit1.size(); i++) {
            //reduction
            ret_t = dag(circuit1.Mref(i))*op[i]*ret_t*circuit2.Mref(i);
        }

        return ret_t.cplx();
    }

} // namespace qcircuit
