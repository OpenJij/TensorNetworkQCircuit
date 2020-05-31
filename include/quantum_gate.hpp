//Copyright (c) 2020 Jij Inc.

#pragma once
#include "itensor/all.h"

namespace qcircuit {
    using namespace itensor;

    class Gate {
    public:
        virtual ITensor op(const std::vector<Index>& slist) const = 0;
    };

    /**
     * @brief Class to represent a quantum gate acting on a site
     */
    class OneSiteGate : public Gate {
    public:
        const size_t site; //!< @brief Site (physical index) ID

        OneSiteGate(size_t site) : site(site) {}
    };

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

    class H : public OneSiteGate {
    public:
        H(size_t site) : OneSiteGate(site) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return (1.0/sqrt(2.0))*(Proj_0(site).op(slist) + Proj_0_to_1(site).op(slist))
                 + (1.0/sqrt(2.0))*(Proj_1(site).op(slist) - Proj_1_to_0(site).op(slist));
        }
    };


    class TwoSiteGate : public Gate {
    public:
        const size_t site1; //!< @brief Site (physical index) ID
        const size_t site2; //!< @brief Site (physical index) ID

        TwoSiteGate(size_t site1, size_t site2) : site1(site1), site2(site2) {}
    };

    class CNOT : public TwoSiteGate {
    public:
        CNOT(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*X(site2).op(slist);
        }
    };

    class CY : public TwoSiteGate {
    public:
        CY(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*Y(site2).op(slist);
        }
    };

    class CZ : public TwoSiteGate {
    public:
        CZ(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*Z(site2).op(slist);
        }
    };
}
