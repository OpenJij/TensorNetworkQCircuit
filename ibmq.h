//Copyright (c) 2019 Kohji Nishimura
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#pragma once

#include "itensor/all.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <iostream>
#include <array>

namespace ibmq {

    using namespace itensor;

    /******
     * 1-indexed vector
     ******/

    template<typename _Tp, typename _Alloc = std::allocator<_Tp>>
        class myvector : public std::vector<_Tp, _Alloc>
    {
        using _Base = std::vector<_Tp, _Alloc>;

        public:
        myvector() : _Base(){}

        explicit
            myvector(size_t __n, const _Alloc& __a = _Alloc()) : _Base(__n, __a){}

        myvector(const myvector& __x) : _Base(__x){}
        myvector(myvector&& __x) noexcept : _Base(__x){}

        myvector(const _Base& __x) : _Base(__x){}
        myvector(_Base&& __x) noexcept : _Base(__x){}

        myvector& operator=(const myvector& __x){
            this->clear();
            for(auto& elem : __x){
                this->push_back(elem);
            }
            return *this;
        }

        //index from one

        _Tp& operator[](size_t __n){
            return _Base::operator[](__n - 1);
        }

        const _Tp& operator[](size_t __n) const{
            return _Base::operator[](__n - 1);
        }
    };

    /**
     * @brief IBM Q circuit
     */
    class IBMQPeps{
        private:
            myvector<Index> a;   //Links
            myvector<Index> s;   // Sites
            myvector<ITensor> M; // ITensors
            ITensor Psi; //\Psi

            std::pair<std::size_t, std::size_t> cursor;

            //[(link_ind, node_ind), ...]
            myvector<myvector<std::pair<size_t, size_t>>> link_node;

            void generate_link(size_t ind1, size_t ind2){
                //generate link between ind1 and ind2
                a.push_back(Index(1, "LinkInd"));
                size_t link_ind = a.size(); //1-indexed

                link_node[ind1].push_back(std::make_pair(link_ind, ind2));
                link_node[ind2].push_back(std::make_pair(link_ind, ind1));
            }

        public:
            constexpr static size_t NUM_BITS = 53;

            IBMQPeps(const std::array<std::pair<double, double>, NUM_BITS>& init_state, const myvector<Index> site_ind = myvector<Index>()) : s(site_ind), link_node(NUM_BITS){
                //generate tensors
                generate_link(1,2);
                generate_link(2,3);
                generate_link(3,4);
                generate_link(4,5);

                generate_link(1,6);
                generate_link(5,7);
                generate_link(6,8);
                generate_link(7,12);

                generate_link(8,9);
                generate_link(9,10);
                generate_link(10,11);
                generate_link(11,12);

                generate_link(8,13);
                generate_link(12,14);
                generate_link(13,15);
                generate_link(14,16);
                generate_link(15,17);
                generate_link(16,19);

                generate_link(10,18);

                generate_link(17,20);
                generate_link(19,21);
                generate_link(20,22);
                generate_link(21,23);
                generate_link(22,24);
                generate_link(23,28);

                generate_link(18,26);

                generate_link(24,25);
                generate_link(25,26);
                generate_link(26,27);
                generate_link(27,28);

                generate_link(24,29);
                generate_link(28,30);
                generate_link(29,31);
                generate_link(30,35);

                generate_link(31,32);
                generate_link(32,33);
                generate_link(33,34);
                generate_link(34,35);

                generate_link(31,36);
                generate_link(35,37);
                generate_link(36,38);
                generate_link(37,39);
                generate_link(38,40);
                generate_link(39,42);

                generate_link(33,41);

                generate_link(40,43);
                generate_link(42,44);
                generate_link(43,45);
                generate_link(44,46);
                generate_link(45,47);
                generate_link(46,51);

                generate_link(41,49);

                generate_link(47,48);
                generate_link(48,49);
                generate_link(49,50);
                generate_link(50,51);

                generate_link(47,52);
                generate_link(51,53);

                //sites
                if(s.size() == 0){
                    for(auto i : range1(NUM_BITS)){
                        s.push_back(Index(2, "SiteInd"));
                    }
                }

                //tensors
                for(auto i : range1(NUM_BITS)){
                    switch(link_node[i].size()){
                        case 1:
                            M.push_back(ITensor(s[i], a[link_node[i][1].first]));
                            M[i].set(a[link_node[i][1].first](1), s[i](1), init_state[i-1].first);
                            M[i].set(a[link_node[i][1].first](1), s[i](2), init_state[i-1].second);
                            break;
                        case 2:
                            M.push_back(ITensor(s[i], a[link_node[i][1].first], a[link_node[i][2].first]));
                            M[i].set(a[link_node[i][1].first](1), a[link_node[i][2].first](1), s[i](1), init_state[i-1].first);
                            M[i].set(a[link_node[i][1].first](1), a[link_node[i][2].first](1), s[i](2), init_state[i-1].second);
                            break;
                        case 3:
                            M.push_back(ITensor(s[i], a[link_node[i][1].first], a[link_node[i][2].first], a[link_node[i][3].first]));
                            M[i].set(a[link_node[i][1].first](1), a[link_node[i][2].first](1), a[link_node[i][3].first](1), s[i](1), init_state[i-1].first);
                            M[i].set(a[link_node[i][1].first](1), a[link_node[i][2].first](1), a[link_node[i][3].first](1), s[i](2), init_state[i-1].second);
                            break;
                        default:
                            assert(false);
                    }
                }

                cursor.first = 1;
                cursor.second = 2;

                Psi = M[cursor.first]*M[cursor.second];
            }

            void decomposePsi(Args args = Args::global()){
                ITensor U,S,V;

                //fetch nodelist in index "first"
                auto list = link_node[cursor.first];
                //get link index
                size_t link_ind = 0;
                for(auto&& i : list){
                    if(i.second == cursor.second){
                        link_ind = i.first;
                    }
                }
                //remove the element which connects to node second
                auto pend = std::remove_if(list.begin(), list.end(), [&](std::pair<size_t, size_t> t){return t.second == cursor.second;});
                myvector<std::pair<size_t, size_t>> newlist;
                for(auto it = list.begin(); it != pend; ++it){
                    newlist.push_back(*it);
                }
                switch(newlist.size()){
                    case 0:
                        U = ITensor(s[cursor.first]);
                        break;
                    case 1:
                        U = ITensor(s[cursor.first], a[newlist[1].first]);
                        break;
                    case 2:
                        U = ITensor(s[cursor.first], a[newlist[1].first], a[newlist[2].first]);
                        break;
                    default:
                        assert(false);
                }

                Spectrum spec = svd(Psi, U, S, V, args);
				a[link_ind] = commonIndex(U,S);
				S /= norm(S); //normalize
				M[cursor.first] = U;
				M[cursor.second] = S*V;
            }

            Spectrum shift_to(size_t ind, const Args& args = Args::global()){
                assert(ind != cursor.first);
                assert(ind != cursor.second);

                //search ind from list
                int flag = 0; //0 not found 1 connected with cursor.first  2 connected with cursor.second

                //first
                for(auto&& i : link_node[cursor.first]){
                    if(i.second == ind){
                        flag = 1;
                        break;
                    }
                }
                //second
                for(auto&& i : link_node[cursor.second]){
                    if(i.second == ind){
                        flag = 2;
                        break;
                    }
                }

                Spectrum spec;

                if(flag == 1){
                    ITensor U,S,V;

                    //fetch nodelist in index "second"
                    auto list = link_node[cursor.second];
                    //get link index
                    size_t link_ind = 0;
                    for(auto&& i : list){
                        if(i.second == cursor.first){
                            link_ind = i.first;
                        }
                    }
                    //remove the element which connects to node first
                    auto pend = std::remove_if(list.begin(), list.end(), [&](std::pair<size_t, size_t> t){return t.second == cursor.first;});
                    myvector<std::pair<size_t, size_t>> newlist;
                    for(auto it = list.begin(); it != pend; ++it){
                        newlist.push_back(*it);
                    }
                    switch(newlist.size()){
                        case 0:
                            V = ITensor(s[cursor.second]);
                            break;
                        case 1:
                            V = ITensor(s[cursor.second], a[newlist[1].first]);
                            break;
                        case 2:
                            V = ITensor(s[cursor.second], a[newlist[1].first], a[newlist[2].first]);
                            break;
                        default:
                            assert(false);
                    }

                    spec = svd(Psi, U, S, V, args);

                    a[link_ind] = commonIndex(S,V);
                    S /= norm(S); //normalize
                    M[cursor.second] = V;
                    Psi = M[ind]*U*S;

                    cursor.second = cursor.first;
                    cursor.first = ind;

                }
                else if(flag == 2){
                    ITensor U,S,V;

                    //fetch nodelist in index "first"
                    auto list = link_node[cursor.first];
                    //get link index
                    size_t link_ind = 0;
                    for(auto&& i : list){
                        if(i.second == cursor.second){
                            link_ind = i.first;
                        }
                    }
                    //remove the element which connects to node second
                    auto pend = std::remove_if(list.begin(), list.end(), [&](std::pair<size_t, size_t> t){return t.second == cursor.second;});
                    myvector<std::pair<size_t, size_t>> newlist;
                    for(auto it = list.begin(); it != pend; ++it){
                        newlist.push_back(*it);
                    }
                    switch(newlist.size()){
                        case 0:
                            U = ITensor(s[cursor.first]);
                            break;
                        case 1:
                            U = ITensor(s[cursor.first], a[newlist[1].first]);
                            break;
                        case 2:
                            U = ITensor(s[cursor.first], a[newlist[1].first], a[newlist[2].first]);
                            break;
                        default:
                            assert(false);
                    }

                    spec = svd(Psi, U, S, V, args);

                    a[link_ind] = commonIndex(U,S);
                    S /= norm(S); //normalize
                    M[cursor.first] = U;
                    Psi = S*V*M[ind];

                    cursor.first = cursor.second;
                    cursor.second = ind;
                }
                else{
                    assert(false && "cannot move to this direction");
                }

                return spec;

            }

            void apply(const ITensor& Op){
				assert(Op.inds().size() == 4);
				for(auto&& elem : Op.inds()){
					//check if 
					assert(elem == s[cursor.first] || elem == s[cursor.second] || elem == prime(s[cursor.first]) || elem == prime(s[cursor.second]));
				}
				this->Psi = Op * prime(Psi,s[cursor.first],s[cursor.second]);
            }

            inline void normalize(){
                Psi /= norm(Psi);
            }

            void primeAll(){
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

            inline const ITensor& Mref(size_t i) const{
				assert(1<=i && i<=NUM_BITS);
                return this->M[i];
            }

            inline const myvector<ITensor>& Mref() const{
                return this->M;
            }

            inline const ITensor& Psiref() const{
                return this->Psi;
            }

            inline const Index& site(size_t i) const{
				assert(1<=i && i<=NUM_BITS);
                return this->s[i];
            }

            inline const myvector<Index>& site() const{
                return this->s;
            }

            inline const std::pair<size_t, size_t>& pos() const{
                return this->cursor;
            }

            void PrintMat(){
                for(auto&& elem : M){
                    PrintData(elem);
                }
                std::cout << "-----------" << std::endl;
                PrintData(Psi);
            }

            void PrintCursor(){
                std::cout << "(" << cursor.first << "," << cursor.second << ")" << std::endl;
            }




    };

	//ITensor Hamiltonian
	inline const ITensor Id(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(1), 1);
		ret.set(s(2), prime(s)(2), 1);
		return ret;
	}

    //pauli x
	inline const ITensor X(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(2), 1);
		ret.set(s(2), prime(s)(1), 1);
		return ret;
	}

    //pauli y
	inline const ITensor Y(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(2), -1_i);
		ret.set(s(2), prime(s)(1), 1_i);
		return ret;
	}

    //pauli z
	inline const ITensor Z(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(1), 1);
		ret.set(s(2), prime(s)(2), -1);
		return ret;
	}

    //|0><0|
	inline const ITensor proj_0(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(1), 1);
		return ret;
	}

    //|1><1|
	inline const ITensor proj_1(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(2), prime(s)(2), 1);
		return ret;
	}

    //|1><0|
	inline const ITensor proj_0_to_1(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(2), prime(s)(1), 1);
		return ret;
	}

    //|0><1|
	inline const ITensor proj_1_to_0(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(2), 1);
		return ret;
	}

    //Hadamard
	inline const ITensor H(const Index& s){
		return (1.0/sqrt(2.0))*(proj_0(s) + proj_0_to_1(s)) + (1.0/sqrt(2.0))*(proj_1(s) - proj_1_to_0(s));
	}

    //CNOT
	inline const ITensor CNOT(const Index& s1, const Index& s2){
		return proj_0(s1)*Id(s2) + proj_1(s1)*X(s2);
	}

    //Control-Y
	inline const ITensor CY(const Index& s1, const Index& s2){
		return proj_0(s1)*Id(s2) + proj_1(s1)*Y(s2);
	}

    //Control-Z
	inline const ITensor CZ(const Index& s1, const Index& s2){
		return proj_0(s1)*Id(s2) + proj_1(s1)*Z(s2);
	}




	Cplx overlap(const IBMQPeps& t_peps1, const myvector<ITensor>& op, const IBMQPeps& t_peps2, const Args& args = Args::global()){
		assert(op.size() == IBMQPeps::NUM_BITS);

		//create copy
		IBMQPeps peps1 = t_peps1;
		IBMQPeps peps2 = t_peps2;

		peps1.decomposePsi(args);
		peps2.decomposePsi(args);

		peps2.primeAll();

		ITensor ret_t = dag(peps1.Mref(1))*op[1]*peps2.Mref(1);
		for(size_t i=2; i<= IBMQPeps::NUM_BITS; i++){
			//reduction
			ret_t = dag(peps1.Mref(i))*op[i]*ret_t*peps2.Mref(i);
		}
		
		return ret_t.cplx();
	}



















} // namespace ibmq
