//Copyright (c) 2020 Jij Inc.

#pragma once
#include "itensor/all.h"

namespace qcircuit {
    using namespace itensor;

    /**
     * @brief Abstract class to represent a quantum gate.
     */
    class Gate {
    public:
        /**
         * @brief returns corresponding tensor operator.
         */
        virtual ITensor op(const std::vector<Index>& slist) const = 0;
    };

    /**
     * @brief Abstract class to represent a one-site quantum gate.
     *
     * This class holds only the site ID on which the gate acts.
     */
    class OneSiteGate : public Gate {
    public:
        const size_t site; //!< @brief Site (physical index) ID

        OneSiteGate(size_t site) : site(site) {}
    };

    /**
     * @brief Identity gate.
     */
    class Id : public OneSiteGate {
    public:
        Id(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(1), prime(s)(1), 1);
            ret.set(s(2), prime(s)(2), 1);
            return ret;
        }
    };

    /**
     * @brief Pauli X gate.
     */
    class X : public OneSiteGate {
    public:
        X(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(1), prime(s)(2), 1);
            ret.set(s(2), prime(s)(1), 1);
            return ret;
        }
    };

    /**
     * @brief Pauli Y gate.
     */
    class Y : public OneSiteGate {
    public:
        Y(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(1), prime(s)(2), -1_i);
            ret.set(s(2), prime(s)(1), 1_i);
            return ret;
        }
    };

    /**
     * @brief Pauli Z gate.
     */
    class Z : public OneSiteGate {
    public:
        Z(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(1), prime(s)(1), 1);
            ret.set(s(2), prime(s)(2), -1);
            return ret;
        }
    };

    /**
     * @brief Projection onto |0>, i.e. |0><0|.
     */
    class Proj_0 : public OneSiteGate {
    public:
        Proj_0(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(1), prime(s)(1), 1);
            return ret;
        }
    };

    /**
     * @brief Projection onto |1>, i.e. |1><1|.
     */
    class Proj_1 : public OneSiteGate {
    public:
        Proj_1(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(2), prime(s)(2), 1);
            return ret;
        }
    };

    /**
     * @brief Map from |0> to |1> (Pauli S-), i.e. |1><0|.
     */
    class Proj_0_to_1 : public OneSiteGate {
    public:
        Proj_0_to_1(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(2), prime(s)(1), 1);
            return ret;
        }
    };

    /**
     * @brief Map from |1> to |0> (Pauli S+), i.e. |0><1|.
     */
    class Proj_1_to_0 : public OneSiteGate {
    public:
        Proj_1_to_0(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s(1), prime(s)(2), 1);
            return ret;
        }
    };

    /**
     * @brief Hadamard gate.
     */
    class H : public OneSiteGate {
    public:
        H(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return (1.0/sqrt(2.0))*(Proj_0(site).op(slist) + Proj_0_to_1(site).op(slist))
                 + (1.0/sqrt(2.0))*(Proj_1(site).op(slist) - Proj_1_to_0(site).op(slist));
        }
    };






    /**
     * @brief Abstract class to represent a two-site quantum gate.
     *
     * This class holds only the site IDs on which the gate acts.
     */
    class TwoSiteGate : public Gate {
    public:
        const size_t site1; //!< @brief Site (physical index) ID
        const size_t site2; //!< @brief Site (physical index) ID

        TwoSiteGate(size_t site1, size_t site2) : site1(site1), site2(site2) {}
    };

    /**
     * @brief Controlled NOT gate (cX gate).
     */
    class CNOT : public TwoSiteGate {
    public:
        CNOT(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*X(site2).op(slist);
        }
    };

    /**
     * @brief cY gate.
     */
    class CY : public TwoSiteGate {
    public:
        CY(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*Y(site2).op(slist);
        }
    };

    /**
     * @brief cZ gate.
     */
    class CZ : public TwoSiteGate {
    public:
        CZ(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*Z(site2).op(slist);
        }
    };
}
