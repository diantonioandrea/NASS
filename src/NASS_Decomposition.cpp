/**
 * @file NASS_Decomposition.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <cmath>

#ifndef NDEBUG // Assertions.
#include <cassert>
#endif

#include "../include/Vectors.hpp"
#include "../include/Matrix.hpp"
#include "../include/Decomposition.hpp"

namespace nass {
    namespace internal {

        /**
         * @brief QR pivoting.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param Nvt0 Natural vector [Nv], target [t].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @param N2 Natural number [N].
         */
        void Pv_RmtNvtNNN_0(real_t* Rmt0, natural_t* Nvt0, const natural_t& N0, const natural_t& N1, const natural_t& N2) {
            real_t R0 = Nr_RvN_R(Rmt0 + N2 * (N0 + 1), N0 - N2);
            natural_t N3 = N2;

            for(natural_t N4 = N2 + 1; N4 < N1; ++N4) {
                const real_t R1 = Nr_RvN_R(Rmt0 + N4 * N0 + N2, N0 - N2);

                if(R1 > R0) {
                    R0 = R1;
                    N3 = N4;
                }
            }

            if(N3 != N2) { // [!] Pivoting can be made implicit.
                for(natural_t N4 = 0; N4 < N0; ++N4) {
                    real_t R1 = Rmt0[N2 * N0 + N4];
                    Rmt0[N2 * N0 + N4] = Rmt0[N3 * N0 + N4];
                    Rmt0[N3 * N0 + N4] = R1;
                }

                natural_t N4 = Nvt0[N2];
                Nvt0[N2] = Nvt0[N3];
                Nvt0[N3] = N4;
            }
        }


        /**
         * @brief Left Householder product.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @param N2 Natural number [N].
         */
        void Lh_RmtRvNNN_0(real_t* Rmt0, const real_t* Rv0, const natural_t& N0, const natural_t& N1, const natural_t& N2) {
            for(natural_t N3 = N0 - N2; N3 < N1; ++N3) {
                const natural_t N4 = N3 * N0;
                const natural_t N5 = N0 - N2;

                const real_t R0 = 2.0 * Dt_RvRvN_R(Rmt0 + N4 + N5, Rv0 + N5, N2);

                for(natural_t N6 = N0 - N2; N6 < N0; ++N6)
                    Rmt0[N4 + N6] -= Rv0[N6] * R0;
            }
        }


        /**
         * @brief Right Householder product.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
        void Rh_RmtRvNN_0(real_t* Rmt0, const real_t* Rv0, const natural_t& N0, const natural_t& N1) {
            for(natural_t N3 = 0; N3 < N0; ++N3) {
                real_t R0 = 0.0;

                for(natural_t N5 = N0 - N1; N5 < N0; ++N5)
                    R0 += Rmt0[N5 * N0 + N3] * Rv0[N5];

                R0 *= 2.0;

                for(natural_t N5 = N0 - N1; N5 < N0; ++N5)
                    Rmt0[N5 * N0 + N3] -= R0 * Rv0[N5];
            }
        }


        /**
         * @brief Left Householder product.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
        void Lh_RvtRvNN_0(real_t* Rvt0, const real_t* Rv0, const natural_t& N0, const natural_t& N1) {
            const real_t R0 = 2.0 * Dt_RvRvN_R(Rvt0 + N0 - N1, Rv0 + N0 - N1, N1);

            for(natural_t N2 = N0 - N1; N2 < N0; ++N2)
                Rvt0[N2] -= R0 * Rv0[N2];
        }

        
        /**
         * @brief (Thin) QR decomposition.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param Rmt1 Real matrix [Rm], target [t].
         * @param Rmt1 Real matrix [Rm], target [t].
         * @param Nvt0 Natural vector [Nv], target [t].
         * @param Rvt0 Real vector [Rm], target [t].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @param N2 Natural number [N].
         */
        void TQR_RmtRmtNvtRvtNN_0(real_t* Rmt0, real_t* Rmt1, natural_t* Nvt0, real_t* Rvt0, const natural_t& N0, const natural_t& N1) {
            #ifndef NDEBUG // Integrity check.
            assert(N1 < N0 - 1);
            #endif

            // Pivoting.
            for(natural_t N2 = 0; N2 < N1; ++N2)
                Nvt0[N2] = N2;

            // Column.
            real_t* Rv0 = new real_t[N0];

            // First column.

            // // Pivoting.
            // Pv_RmtNvtNNN_0(Rmt1, Nvt0, N0, N1, 0);

            // Copy.
            Cp_RvtRvN_0(Rv0, Rmt1, N0);

            // Direction.
            Rv0[0] -= std::copysign(1.0, Rv0[0]) * Nr_RvN_R(Rv0, N0);

            // Normalization.
            Nrz_RvtN_0(Rv0, N0);

            // Q.
            for(natural_t N2 = 0; N2 < N0; ++N2) {
                for(natural_t N3 = 0; N3 < N2; ++N3)
                    Rmt0[N2 * N0 + N3] = -2.0 * Rv0[N2] * Rv0[N3];

                for(natural_t N3 = N2 + 1; N3 < N0; ++N3)
                    Rmt0[N2 * N0 + N3] = -2.0 * Rv0[N2] * Rv0[N3];

                Rmt0[N2 * (N0 + 1)] = 1.0 - 2.0 * Rv0[N2] * Rv0[N2];
            }

            // R.
            Lh_RmtRvNNN_0(Rmt1, Rv0, N0, N1, N0);

            // Vector.
            Lh_RvtRvNN_0(Rvt0, Rv0, N0, N0);

            // Other columns.
            for(natural_t N2 = 1; N2 < N1; ++N2) {

                // // Pivoting.
                // Pv_RmtNvtNNN_0(Rmt1, Nvt0, N0, N1, N2);

                // Copy.
                Cp_RvtRvN_0(Rv0 + N2, Rmt1 + N2 * (N0 + 1), N0 - N2);

                // Direction.
                Rv0[N2] -= std::copysign(1.0, Rv0[N2]) * Nr_RvN_R(Rv0 + N2, N0 - N2);

                // Normalization.
                Nrz_RvtN_0(Rv0 + N2, N0 - N2);

                // Q.
                Rh_RmtRvNN_0(Rmt0, Rv0, N0, N0 - N2);

                // R.
                Lh_RmtRvNNN_0(Rmt1, Rv0, N0, N1, N0 - N2);

                // Vector.
                Lh_RvtRvNN_0(Rvt0, Rv0, N0, N0 - N2);
            }

            delete[] Rv0;
        }

    }
}