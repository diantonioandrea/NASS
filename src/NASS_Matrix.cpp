/**
 * @file NASS_Matrix.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Matrix.hpp implementation.
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <cstdlib>
#include <ctime>
#include <print>

#include <cmath>
#define _USE_MATH_DEFINES

#include "../include/Reals.hpp"
#include "../include/Matrix.hpp"

namespace nass {
    namespace internal {

        /**
         * @brief 
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rm0 Real matrix [Rm].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
        void Ml_RvtRmRvNN_0(real_t* Rvt0, const real_t* Rm0, const real_t* Rv0, const natural_t& N0, const natural_t& N1) {
            #pragma omp parallel
            {
                for(natural_t N2 = 0; N2 < N1; ++N2) {
                    const real_t R0 = Rv0[N2];

                    #pragma omp for
                    for(natural_t N3 = 0; N3 < N0; ++N3)
                        Rvt0[N3] += Rm0[N2 * N0 + N3] * R0;
                }
            }
        }


        /**
         * @brief Products between two real_t*.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param TRm0 Real matrix [Rm], transposed.
         * @param Rm0 Real matrix [Rm].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @param N2 Natural number [N].
         */
        void Ml_RmtTRmRmNNN_0(real_t* Rmt0, const real_t* TRm0, const real_t* Rm0, const natural_t& N0, const natural_t& N1, const natural_t& N2) {
            for(natural_t N3 = 0; N3 < N2; ++N3)
                for(natural_t N4 = 0; N4 < N1; ++N4)
                    for(natural_t N5 = 0; N5 < N0; ++N5)
                        Rmt0[N3 * N1 + N4] += TRm0[N4 * N0 + N5] * Rm0[N3 * N0 + N5];
        }


        /**
         * @brief Products between two real_t*.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param Rm0 Real matrix [Rm].
         * @param Rm1 Real matrix [Rm].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @param N2 Natural number [N].
         */
        void Ml_RmtRmRmNNN_0(real_t* Rmt0, const real_t* Rm0, const real_t* Rm1, const natural_t& N0, const natural_t& N1, const natural_t& N2) {
            for(natural_t N3 = 0; N3 < N2; ++N3)
                for(natural_t N4 = 0; N4 < N1; ++N4)
                    for(natural_t N5 = 0; N5 < N0; ++N5)
                        Rmt0[N3 * N0 + N5] += Rm0[N4 * N0 + N5] * Rm1[N3 * N1 + N4];
        }


        /**
         * @brief Gauss embedding.
         * 
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
        [[nodiscard]] real_t* Gs_NN_Rm(const natural_t& N0, const natural_t& N1) {
            std::srand(std::time(nullptr));

            real_t* Rm0 = new real_t[N0 * N1];

            #pragma omp parallel for
            for(natural_t N2 = 0; N2 < N0 * N1; ++N2) {
                const real_t R0 = static_cast<real_t>(std::rand()) / RAND_MAX;

                Rm0[N2] = std::sqrt(-2.0 * std::log(R0)) * std::cos(2.0 * M_PI * R0);
            }

            return Rm0;
        }


        /**
         * @brief Prints a (matrix) real_t*.
         * 
         * @param Rm0 Real matrix [Rm].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
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