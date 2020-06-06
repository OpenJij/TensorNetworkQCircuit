//Copyright (c) 2020 Jij Inc.


#pragma once

#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <fstream>

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
        std::vector<std::vector<Neighbor>> neighbors_list; //!< @brief List of neighboring sites of each site.
        size_t num_links; //!< @brief Number of links
        const size_t num_bits; //!< @brief Number of qubits

    public:
        CircuitTopology(size_t num_bits) : num_bits(num_bits),
                                           num_links(0),
                                           neighbors_list(num_bits) {}

        void generateLink(size_t site1, size_t site2) {
            // TODO: assert or throw Exception when index exceeds number of bits.
            // TODO: assert or throw Exception when multiple link is specified.

            neighbors_list[site1].emplace_back(site2, num_links);
            neighbors_list[site2].emplace_back(site1, num_links);
            num_links++;
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

            // if not necessary to move
            if((origin.first == destination.first &&  origin.second == destination.second) ||
               (origin.first == destination.second && origin.second == destination.first)) {
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
                result.push_back(destination.second);
                site = destination.first;
            } else if(back_to[destination.second] != NOT_YET_REACHED) {
                result.push_back(destination.first);
                site = destination.second;
            } else {
                // TODO: assert or throw Exception.
                // If reach this else branch, given sites are unreachable.
            }

            while(site != origin.first && site != origin.second) {
                result.push_back(site);
                site = back_to[site];
            }
            std::reverse(result.begin(), result.end());
            return result;
        }

        /**
         * @brief Exports graph topology as DOT language (Graphviz style) to given file.
         *
         * To generate PDF file from exported DOT file, run the following:  
         * `dot -Tpdf <name>.dot -o <name>.pdf`
         */
        void exportDotString(const std::string& filename,
                             const std::string& layout = "neato",
                             const std::string& shape = "circle") const {
            std::ofstream stream(filename);

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
        }
    };
} // namespace qcircuit
