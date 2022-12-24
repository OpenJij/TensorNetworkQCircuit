//Copyright (c) 2020 Jij Inc.

#pragma once
#include <complex>
#include <itensor/all.h>

namespace qcircuit {
    using namespace itensor;

    /**
     * @brief Abstract class to represent a quantum gate.
     */
    class Gate {
    public:
        virtual ~Gate() {}

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
        virtual ~OneSiteGate() {}
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
            ret.set(s=1, prime(s)=1, 1);
            ret.set(s=2, prime(s)=2, 1);
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
            ret.set(s=1, prime(s)=2, 1);
            ret.set(s=2, prime(s)=1, 1);
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
            ret.set(s=1, prime(s)=2, -1_i);
            ret.set(s=2, prime(s)=1, 1_i);
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
            ret.set(s=1, prime(s)=1, 1);
            ret.set(s=2, prime(s)=2, -1);
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
            ret.set(s=1, prime(s)=1, 1);
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
            ret.set(s=2, prime(s)=2, 1);
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
            ret.set(s=2, prime(s)=1, 1);
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
            ret.set(s=1, prime(s)=2, 1);
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
                 + (1.0/sqrt(2.0))*(Proj_1_to_0(site).op(slist) - Proj_1(site).op(slist));
        }
    };

    /**
     * @brief Phase gate.
     */
    class P : public OneSiteGate {
    public:
        const double theta;

        P(size_t site, double theta) : OneSiteGate(site), theta(theta) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s = slist[site];
            ITensor ret(s, prime(s));
            ret.set(s=1, prime(s)=1, 1);
            ret.set(s=2, prime(s)=2, std::exp(1_i * theta));
            return ret;
        }
    };

    /**
     * @brief Universal unitary gate.
     *
     * U(theta, phi, lambda) = R_z(phi)R_y(theta)R_z(lambda),
     * where R_z(phi) = exp(-i theta Z / 2),  R_y(phi) = exp(-i theta Y / 2).
     * Z and Y are the Pauli matrices for each axis.
     * This can express any element of SU(2).
     */
    class UniversalUnitary : public OneSiteGate {
    public:
      const double theta;
      const double phi;
      const double lambda;

      UniversalUnitary(size_t site,
                       double theta, double phi, double lambda) :
       OneSiteGate(site), theta(theta), phi(phi), lambda(lambda) {}

      ITensor op(const std::vector<Index>& slist) const override {
          std::complex<double> alpha = std::exp(-1_i*(phi+lambda)/2)*cos(theta/2);
          std::complex<double> beta = -std::exp(-1_i*(phi-lambda)/2)*sin(theta/2);

          auto s = slist[site];
          ITensor ret(s, prime(s));
          ret.set(s=1, prime(s)=1, alpha);
          ret.set(s=1, prime(s)=2, beta);
          ret.set(s=2, prime(s)=1, -std::conj(beta));
          ret.set(s=2, prime(s)=2, std::conj(alpha));

          return ret;
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
        virtual ~TwoSiteGate() {}
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

    /**
     * @brief Controlled Phase gate.
     */
    class CP : public TwoSiteGate {
    public:
        const double theta;

        CP(size_t site1, size_t site2, double theta) : TwoSiteGate(site1, site2), theta(theta) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*P(site2, theta).op(slist);
        }
    };

    /**
     * @brief Controlled Universal Unitary gate.
     */
    class CUniversalUnitary : public TwoSiteGate {
    public:
        const double theta;
        const double phi;
        const double lambda;

        CUniversalUnitary(size_t site1, size_t site2, double theta, double phi, double lambda) : TwoSiteGate(site1, site2), theta(theta), phi(phi), lambda(lambda) {}

        ITensor op(const std::vector<Index>& slist) const override {
            return Proj_0(site1).op(slist)*Id(site2).op(slist)
                + Proj_1(site1).op(slist)*UniversalUnitary(site2, theta, phi, lambda).op(slist);
        }
    };

    /**
     * @brief Swap gate.
     */
    class Swap : public TwoSiteGate {
    public:
        Swap(size_t site1, size_t site2) : TwoSiteGate(site1, site2) {}

        ITensor op(const std::vector<Index>& slist) const override {
            auto s1 = slist[site1];
            auto s2 = slist[site2];
            ITensor ret(s1, prime(s1), s2, prime(s2));
            ret.set(s1=1, prime(s1)=1, s2=1, prime(s2)=1, 1);
            ret.set(s1=2, prime(s1)=2, s2=2, prime(s2)=2, 1);
            ret.set(s1=1, prime(s1)=2, s2=2, prime(s2)=1, 1);
            ret.set(s1=2, prime(s1)=1, s2=1, prime(s2)=2, 1);

            return ret;
        }
    };
}
