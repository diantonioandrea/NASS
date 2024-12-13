/**
 * @file NASS_Matrix.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Matrix.hpp implementation.
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <print>

#include "../include/Reals.hpp"
#include "../include/Matrix.hpp"

namespace nass {
    namespace internal {

        void Pr_RmNN_0(const real_t* Rm0, const natural_t& N0, const natural_t& N1) {
            std::println("--- Matrix.");

            for(natural_t N2 = 0; N2 < N0; ++N2) {
                for(natural_t N3 = 0; N3 < N1 - 1; ++N3)
                    Pr_R_0(Rm0[N3 * N0 + N2]);
                
                Pn_R_0(Rm0[(N1 - 1) * N0 + N2]);
            }

            std::println("---");
        }

    }
}