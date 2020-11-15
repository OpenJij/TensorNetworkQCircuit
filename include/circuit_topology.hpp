//Copyright (c) 2020 Jij Inc.


#pragma once

#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <sstream>
#include <cassert>
#include "qcircuit_exception.hpp"

namespace qcircuit {

    /**
     * @brief class to represent circuit topology
     */
    class CircuitTopology {
    public:
        struct Neighbor {
            size_t site; //!< @brief Site (physical index) ID
            size_t link; //!< @brief Link (bond index) ID

            Neighbor(size_t site, size_t link) : site(site), link(link) {}
        };

    private:
        const size_t num_bits; //!< @brief Number of qubits
        size_t num_links; //!< @brief Number of links
        std::vector<std::vector<Neighbor>> neighbors_list; //!< @brief List of neighboring sites of each site.

    public:
        CircuitTopology(size_t num_bits) : num_bits(num_bits),
                                           num_links(0),
                                           neighbors_list(num_bits) {}

        /**
         * @brief generates a link between `site1` and `site2`.
         */
        void generateLink(size_t site1, size_t site2) {
            /* Check invalid site index */
            if(site1 >= num_bits || site2 >= num_bits) {
                std::stringstream ss;
                ss << "Link can't be generated between (" << site1 << ", " << site2
                   << ") : Index exceeds the number of qubits or ";
                throw QCircuitException(ss.str());
            }
            if(site1 < 0 || site2 < 0) {
                std::stringstream ss;
                ss << "Link can't be generated between (" << site1 << ", " << site2
                   << ") : Negative index specified";
                throw QCircuitException(ss.str());
            }
            if(site1 == site2) {
                std::stringstream ss;
                ss << "Link can't be generated between (" << site1 << ", " << site2
                   << ") : Same indices specified";
                throw QCircuitException(ss.str());
            }

            /* Check whether the specified link already exists. */
            {
                auto itr = std::find_if(neighbors_list[site1].begin(), neighbors_list[site1].end(),
                                        [&](const Neighbor& nei) { return nei.site == site2; });
                if(itr != neighbors_list[site1].end()) {
                    std::stringstream ss;
                    ss << "Link can't be generated between (" << site1 << ", " << site2
                       << ") : Link already exists";
                    throw QCircuitException(ss.str());
                }
            }

            neighbors_list[site1].emplace_back(site2, num_links);
            neighbors_list[site2].emplace_back(site1, num_links);
            num_links++;
        }

        /**
         * @brief returns whether a link exists between `site1` and `site2`.
         */
        bool hasLinkBetween(size_t site1, size_t site2) const {
            auto v = this->neighbors_list[site1];
            auto find_result = std::find_if(v.begin(), v.end(), [&site2] (const Neighbor& x) {
                                                                    return x.site == site2;
                                                                });
            return find_result != v.end();
        }

        /**
         * @brief returns a link ID between `site1` and `site2`.
         */
        size_t getLinkIdBetween(size_t site1, size_t site2) const {
            auto v = this->neighbors_list[site1];
            auto find_result = std::find_if(v.begin(), v.end(), [&site2] (const Neighbor& x) {
                                                                    return x.site == site2;
                                                                });
            if(find_result == v.end()) {
                std::stringstream ss;
                ss << "There is no link between (" << site1 << ", " << site2 << ")";
                throw QCircuitException(ss.str());
            }

            return (*find_result).link;
        }

        size_t numberOfLinks() const {
            return this->num_links;
        }

        size_t numberOfBits() const {
            return this->num_bits;
        }

        const std::vector<Neighbor>& neighborsOf(size_t index) const {
            return this->neighbors_list[index];
        }


        /**
         * @brief Searches the shortest route from `origin` to `destination`.
         *
         * Breadth First Search algorithm is used.
         */
        std::vector<size_t> getRoute(std::pair<size_t, size_t> origin, std::pair<size_t, size_t> destination) const {
            static const int NOT_YET_REACHED = -1;
            std::vector<int> back_to(num_bits, NOT_YET_REACHED);
            // back_to[i] points which neighboring site we go back at ith site
            // to reach `origin` in the shortest path.

            // If not necessary to move (at least one of destination is the same as origins)
            if((origin.first == destination.first ||  origin.first == destination.second) ||
               (origin.second == destination.first || origin.second == destination.second)) {
                return std::vector<size_t>();
            }

            std::queue<size_t> queue;
            queue.push(origin.first);
            back_to[origin.first] = origin.first;
            queue.push(origin.second);
            back_to[origin.second] = origin.second;

            while(!queue.empty()) {
                auto site = queue.front();
                queue.pop(); // pop() returns nothing.

                if(site == destination.first || site == destination.second) {
                    break;
                }

                for(auto neighbor : neighbors_list[site]) {
                    if(back_to[neighbor.site] == NOT_YET_REACHED) {
                        back_to[neighbor.site] = site;
                        queue.push(neighbor.site);
                    }
                }
            }

            size_t site;
            std::vector<size_t> result;
            if(back_to[destination.first] != NOT_YET_REACHED) {
                site = destination.first;
            } else if(back_to[destination.second] != NOT_YET_REACHED) {
                site = destination.second;
            } else {
                std::stringstream ss;
                ss << "Path to (" << destination.first << ", "
                   << destination.second << ") not found";
                throw QCircuitException(ss.str());
            }

            while(site != origin.first && site != origin.second) {
                result.push_back(site);
                site = back_to[site];
            }
            std::reverse(result.begin(), result.end());
            return result;
        }


        /**
         * @brief Find the shortest swapping path to move `target` to
         * one neighboring site of `origin`
         *
         * Breadth First Search algorithm is used. Returned vector contains `target`
         * and does not contain `origin`.
         * If `target` is already a neighbor of `origin`, returned vector contains only `target`.
         *
         * @note The implementation of this function is quite similar to `getRoute()` so
         * this function could be combined with it.
         */
        std::vector<size_t> getSwapPath(size_t origin, size_t target) const {
            static const int NOT_YET_REACHED = -1;
            std::vector<int> back_to(num_bits, NOT_YET_REACHED);
            // back_to[i] points which neighboring site we go back at ith site
            // to reach one neighbor of `destination` in the shortest path.

            if(origin == target) {
                assert(false && "Unintended call (unnecessary to swap)");
            }

            std::queue<size_t> queue;
            queue.push(origin);
            back_to[origin] = origin;

            while(!queue.empty()) {
                auto site = queue.front();
                queue.pop(); // pop() returns nothing.

                if(site == target) {
                    break;
                }

                for(auto neighbor : neighbors_list[site]) {
                    if(back_to[neighbor.site] == NOT_YET_REACHED) {
                        back_to[neighbor.site] = site;
                        queue.push(neighbor.site);
                    }
                }
            }

            size_t site = target;
            std::vector<size_t> result;

            while(site != origin) {
                result.push_back(site);
                site = back_to[site];
            }

            return result;
        }


        /**
         * @brief returns true if the circuit is a connected graph,
         * i.e. there is a path between every pair of qubits (vertices).
         *
         * Breadth First Search algorithm is used.
         */
        bool isConnectedGraph() const {
            std::vector<bool> reached(num_bits, false);

            std::queue<size_t> queue;
            reached[0] = true;
            queue.push(0);

            while(!queue.empty()) {
                auto site = queue.front();
                queue.pop(); // pop() returns nothing.

                for(auto neighbor : neighbors_list[site]) {
                    if(!reached[neighbor.site]) {
                        reached[neighbor.site] = true;
                        queue.push(neighbor.site);
                    }
                }
            }

            return std::all_of(reached.begin(), reached.end(),
                               [](size_t x) { return x; });
        }

        /**
         * @brief Returns  DOT language (Graphviz style) string representing the graph topology.
         *
         * To generate PDF file from exported DOT file, run the following:  
         * `dot -Tpdf <name>.dot -o <name>.pdf`
         */
        std::string convertToDotString(const std::string& layout = "neato",
                             const std::string& shape = "circle") const {
            std::stringstream stream;

            stream << "// Convert to pdf:" << std::endl;
            stream << "// dot -Tpdf <name>.dot -o <name>.pdf" << std::endl;
            stream << std::endl;
            stream << "graph {" << std::endl;
            stream << "    graph[layout=" << layout << "]" << std::endl;
            stream << "    node[shape=" << shape << "]" << std::endl << std::endl;
            for(size_t i = 0;i < num_bits;i++) {
                for(const auto& neighbor : neighbors_list[i]) {
                    if(i > neighbor.site) {
                        stream << "    " << i << " -- " << neighbor.site << ";" << std::endl;
                    }
                }
            }
            stream << "}" << std::endl;

            return stream.str();
        }
    };
} // namespace qcircuit
