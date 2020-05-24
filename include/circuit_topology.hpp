//Copyright (c) 2020 Jij Inc.


#pragma once

#include <vector>

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

        void generate_link(size_t site1, size_t site2) {
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
    };
} // namespace qcircuit
