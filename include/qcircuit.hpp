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
        std::vector<Index> s;     //!< @brief Physical (on-site) indices.
        std::vector<ITensor> M;   //!< @brief Tensor entities.
        std::vector<ITensor> SV;  //!< @brief Singular value entities.
        ITensor Psi;  //!< @brief TPS wave function

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
            s(physical_indices), topology(topology), random_engine(std::random_device()()) {

            if(!topology.isConnectedGraph()) {
                throw QCircuitException("Invalid circuit topology : Some nodes are unreachable");
            }

            /* Initialize link indices */
            std::vector<std::pair<Index, Index>> a;
            a.reserve(topology.numberOfLinks());
            SV.reserve(topology.numberOfLinks());
            for(size_t i = 0; i < topology.numberOfLinks();i++) {
                auto index1 = Index(1, "LinkInd");
                auto index2 = Index(1, "LinkInd");
                a.emplace_back(index1, index2);
                SV.emplace_back(index1, index2);
                SV.back().set(index1=1, index2=1, 1.0);
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
                    if(i < neighbor.site) {
                        ind_list.push_back(a[neighbor.link].first);
                    } else {
                        ind_list.push_back(a[neighbor.link].second);
                    }
                }

                //fill ind_val_list
                ind_val_list.push_back(s[i](1));
                for(const auto& neighbor : neighbors){
                    if(i < neighbor.site) {
                        ind_val_list.push_back((a[neighbor.link].first)(1));
                    } else {
                        ind_val_list.push_back((a[neighbor.link].second)(1));
                    }
                }

                //insert into M
                M.emplace_back(ind_list);
                ind_val_list[0] = s[i](1);
                M[i].set(ind_val_list, init_qubits[i].first);
                ind_val_list[0] = s[i](2);
                M[i].set(ind_val_list, init_qubits[i].second);
            }

            /* set cursor position */
            cursor.first = 0;

            /* find the minimum numbered bit of neighbors of the bit 0. */
            size_t cursor_second_index = this->size();
            for(const auto& neighbor : topology.neighborsOf(cursor.first)) {
                if(neighbor.site < cursor_second_index) {
                    cursor_second_index = neighbor.site;
                }
            }
            cursor.second = cursor_second_index;

            updatePsi();
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

        /** @brief update `Psi` (canonical center) with current cursor position */
        void updatePsi() {
            size_t link_index = topology.getLinkIdBetween(cursor.first, cursor.second);

            Psi = M[cursor.first]*SV[link_index]*M[cursor.second];

            // add singular-value matrices at egdes
            for(auto&& neighbor : topology.neighborsOf(cursor.first)) {
                if(neighbor.link != link_index) {
                    Psi = Psi*SV[neighbor.link];
                }
            }
            for(auto&& neighbor : topology.neighborsOf(cursor.second)) {
                if(neighbor.link != link_index) {
                    Psi = Psi*SV[neighbor.link];
                }
            }
        }

        /** @brief decomposes wave-function at cursor position and truncates bonds.
         *
         * Decomposition is achieved with the following way.
         * Below, symbols "o" means site tensors
         * and symbols "*" means singular-value tensors.
         *
         * 1. Decompose Psi:
         @verbatim
                          Psi
                +---------------------+                     U   S   V
            ----|--*---o---*---o---*--|----  =>  -----------o---*---o-----------
                +---------------------+
         @endverbatim
         * 2. Since Psi includes singular-value tensors at each edge,
         *    U and V must be reconstructed NOT to include them.
         *    To factor out them, first, insert singular-value tensors and their inverses (marked as #):
         @verbatim
                       U   S   V                            U   S   V
            -----------o---*---o----------   =>  ---*---#---o---*---o---#---*---
         @endverbatim
         * 3. Then, regard the part "--#--o--" as U (and also V):
         @verbatim
                       U   S   V                                S
            ---*---#---o---*---o---#---*---  =>  ---*-----o-----*-----o-----*---
                  ^^^^^^^     ^^^^^^^                     U           V
                   new U       new V
         @endverbatim
         *
         */
        Spectrum decomposePsi(const Args& args) {
            const double SINGULAR_VALUE_THRESHOLD = 1e-16;
            // Very small singular values could cause numerical instability
            // when calculating inverse of them.
            // This value is used as the threshold to discard such small singular values.

            /* Prepare indices to be free ones of U */
            auto link_index = topology.getLinkIdBetween(cursor.first, cursor.second);

            std::vector<Index> outer_indices_U;
            outer_indices_U.push_back(s[cursor.first]);
            for(auto&& neighbor : topology.neighborsOf(cursor.first)) {
                if(neighbor.link != link_index) {
                    outer_indices_U.push_back(commonIndex(Psi, SV[neighbor.link]));
                }
            }

            ITensor U, S, V;
            U = ITensor(outer_indices_U);

            Spectrum spec = svd(Psi, U, S, V, args);

            S /= norm(S); // normalization

            /* Reconstruct U */
            for(auto&& neighbor : topology.neighborsOf(cursor.first)) {
                if(neighbor.link != link_index) {
                    Index index_i = uniqueIndex(SV[neighbor.link], U);
                    Index index_j = commonIndex(U, SV[neighbor.link]); // to be contracted

                    ITensor inv(prime(index_j), index_i);
                    for(auto i : range1(std::min(dim(index_i), dim(index_j)))) {
                        // Singular-value matrix has only diagonal elements, so just invert them

                        auto x = SV[neighbor.link].real(index_i=i, index_j=i);
                        if(x < SINGULAR_VALUE_THRESHOLD) {
                            break; // SVs are in descending order, so discard remaining ones
                        }
                        inv.set(prime(index_j)=i, index_i=i, 1.0/x);
                    }
                    U.prime(index_j);
                    U = U*inv;
                }
            }

            /* Reconstruct V */
            for(auto&& neighbor : topology.neighborsOf(cursor.second)) {
                if(neighbor.link != link_index) {
                    Index index_i = uniqueIndex(SV[neighbor.link], V);
                    Index index_j = commonIndex(V, SV[neighbor.link]); // to be contracted

                    ITensor inv(prime(index_j), index_i);
                    for(auto i : range1(std::min(dim(index_i), dim(index_j)))) {
                        // Singular-value matrix has only diagonal elements, so just invert them

                        auto x = SV[neighbor.link].real(index_i=i, index_j=i);
                        if(x < SINGULAR_VALUE_THRESHOLD) {
                            break; // SVs are in descending order, so discard remaining ones
                        }
                        inv.set(prime(index_j)=i, index_i=i, 1.0/x);

                    }
                    V.prime(index_j);
                    V = V*inv;
                }
            }

            SV[link_index] = S;
            M[cursor.first] = U;
            M[cursor.second] = V;

            return spec;
        }

        Spectrum decomposePsi() {
            return decomposePsi(default_args);
        }

        /** @brief shifts cursor position to specified neighboring site `dest`.
         *
         * `direction` specifies which site of the cursor position to be used as the "head" of move.
         * `direction == 0` means not specifying the direction.
         */
        Spectrum shiftCursorTo(size_t dest, int direction, const Args& args) {
            const static int AUTO_HEAD = 0;
            const static int FIRST_AS_HEAD = 1;
            const static int SECOND_AS_HEAD = 2;

            assert(dest != cursor.first);
            assert(dest != cursor.second);
            assert(direction == AUTO_HEAD || direction == FIRST_AS_HEAD || direction == SECOND_AS_HEAD);

            if(direction == AUTO_HEAD) {
                //first
                for(auto&& i : topology.neighborsOf(cursor.first)) {
                    if(i.site == dest) {
                        direction = FIRST_AS_HEAD;
                        break;
                    }
                }
                //second
                for(auto&& i : topology.neighborsOf(cursor.second)) {
                    if(i.site == dest) {
                        direction = SECOND_AS_HEAD;
                        break;
                    }
                }
            }


            Spectrum spec = decomposePsi();


            switch(direction) {
            case FIRST_AS_HEAD:
                cursor.second = cursor.first;
                cursor.first = dest;
                break;
            case SECOND_AS_HEAD:
                cursor.first = cursor.second;
                cursor.second = dest;
                break;
            default:
                assert(false && "cannot move to this direction");
                break;
            }

            updatePsi();

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
          * @brief returns probability the qubit at `site` to be observed as value `expected` (0 or 1),
          * i.e. `<psi|Proj_X(i)|psi>` as "Born rule".
          */
        double probabilityOf(size_t site, int expected) const {
            assert(expected == 0 | expected == 1);

            std::vector<ITensor> op;
            op.reserve(this->size());
            for(size_t i = 0;i < this->size();i++) {
                if(i == site) {
                    if(expected == 0) {
                        op.push_back(generateTensorOp(Proj_0(i)));
                    } else {
                        op.push_back(generateTensorOp(Proj_1(i)));
                    }
                } else {
                    op.push_back(generateTensorOp(Id(i)));
                }
            }

            /*
             * <psi|Proj_X(i)|psi> = |w_{i,X}|^2 is real by definition,
             * where w_{i,X} is the coefficient corresponding to the base |....X....>.
             *                                                     (ith qubit) ^
             */
            return std::real(overlap(*this, op, *this));
        }

        /**
          * @brief returns probability the qubit at `site` to be observed as value 0.
          */
        double probabilityOfZero(size_t site) const {
            return probabilityOf(site, 0);
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

            for(auto& elem : this->M){
                elem = prime(elem);
            }

            for(auto& elem : this->SV){
                elem = prime(elem);
            }

            Psi = prime(Psi);
        }

        const CircuitTopology& getTopology() const {
            return this->topology;
        }

        const ITensor& Mref(size_t i) const {
            assert(0 <= i && i < this->size());
            return this->M[i];
        }

        const std::vector<ITensor>& Mref() const {
            return this->M;
        }

        const ITensor& SVref(size_t i) const {
            assert(0 <= i && i < this->topology.numberOfLinks());
            return this->SV[i];
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

    /* Prototype declaration of this function is placed at the top. */
    inline Cplx overlap(QCircuit circuit1,
                        const std::vector<ITensor>& op,
                        QCircuit circuit2,
                        const Args& args) {
        assert(op.size() == circuit1.size() && op.size() == circuit2.size());

        auto topology = circuit1.getTopology();

        circuit1.decomposePsi(args);
        circuit2.decomposePsi(args);

        circuit2.primeAll();

        ITensor ret_t(1.0); // Rank zero tensor with amplitude 1.0
        for(size_t i = 0; i < circuit1.size(); i++) {
            //reduction
            ret_t = dag(circuit1.Mref(i))*op[i]*ret_t*circuit2.Mref(i);
            for(auto&& neighbor : topology.neighborsOf(i)) {
                if(i < neighbor.site) {
                    ret_t = dag(circuit1.SVref(neighbor.link))*ret_t*circuit2.SVref(neighbor.link);
                }
            }
        }

        return ret_t.cplx();
    }

} // namespace qcircuit
