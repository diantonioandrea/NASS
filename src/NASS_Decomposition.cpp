/**
 * @file NASS_Decomposition.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../include/Vectors.hpp"
#include "../include/Decomposition.hpp"

namespace nass {
    namespace internal {
        
        /**
         * @brief (Thin) QR decomposition.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param Rmt1 Real matrix [Rm], target [t].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
        void TQR_RmtRmtNN_0(real_t* Rmt0, real_t* Rmt1, const natural_t& N0, const natural_t& N1) {
            
            // Column.
            real_t* Rv0 = new real_t[N0];

            // First column.

            // Copy.
            Cp_RvtRvN_0(Rv0, Rmt1, N0);

            // Direction.
            Rv0[0] -= Nr_RvN_R(Rv0, N0);

            // Normalization.
            Nrz_RvtN_0(Rv0, N0);

            // Q.
            for(natural_t N2 = 0; N2 < N0; ++N2) {
                for(natural_t N3 = 0; N3 < N2; ++N3)
                    Rmt0[N2 * N0 + N3] = -2.0 * Rv0[N2] * Rv0[N3];

                for(natural_t N3 = N2 + 1; N3 < N0; ++N3)
                    Rmt0[N2 * N0 + N3] = -2.0 * Rv0[N2] * Rv0[N3];

                Rmt0[N2 * (N0 + 1)] = 1.0 - 2.0 * Rv0[N2] * Rv0[N3];
            }

            // R.
            // [!]

            // Other columns.
            for(natural_t N2 = 1; N2 < N1; ++N2) {

                // Copy.
                Cp_RvtRvN_0(Rv0 + N2, Rmt1 + N2 * (N0 + 1), N0 - N2);

                // Direction.
                // [!]

                // Normalization.
                Nrz_RvtN_0(Rv0 + N2, N0 - N2);

                // Q.
                // [!]

                // R.
                // [!]
            }

            delete[] Rv0;
        }

    }
}