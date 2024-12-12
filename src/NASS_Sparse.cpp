/**
 * @file NASS_Sparse.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Sparse.hpp implementations.
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <fstream>
#include <sstream>
#include <cassert>

#include "../include/Sparse.hpp"

namespace nass {
    namespace internal {

        /**
         * @brief Loads a (CSC) sparse matrix from a file.
         * 
         * @param St String [St]
         * @return std::tuple<natural_t, natural_t*, natural_t*, real_t*> (CSC) sparse matrix [Spc].
         */
        std::tuple<natural_t, natural_t*, natural_t*, real_t*> Spc_St_NNvNvRv(const std::string& St0) {
            natural_t N0, N1, N2, N3 = 0, N4 = 0;
            real_t R0;
            
            // Line and file.
            std::string St1;
            std::ifstream F0(St0);

            // First line.
            std::getline(F0, St1);

            // First line stream.
            std::istringstream Ss0(St1);

            // Parameters.
            Ss0 >> N0; // Rows.
            Ss0 >> N1; // Columns.
            Ss0 >> N2; // Entries.

            // Square check.
            assert(N0 == N1);

            // Initialization.
            natural_t* Nv0 = new natural_t[N1 + 1];
            natural_t* Nv1 = new natural_t[N2];
            real_t* Rv0 = new real_t[N2];

            // Reading.
            while(std::getline(F0, St1)) {

                // Line stream.
                std::istringstream Ss1(St1);

                // Coordinates and entry.
                Ss1 >> N1; --N1; // Row.
                Ss1 >> N2; --N2; // Column.
                Ss1 >> R0; // Entry.

                // Sparse check.
                if(std::abs(R0) < real_tol)
                    continue;

                for(; N3 < N2; ++N3)
                    Nv0[N3 + 1] = N4;

                Nv1[N4] = N1;
                Rv0[N4] = R0;

                ++N4;
            }

            // Last entries.
            for(; N3 < N0; ++N3)
                Nv0[N3 + 1] = N4;

            return {N0, Nv0, Nv1, Rv0};
        }

    }
}