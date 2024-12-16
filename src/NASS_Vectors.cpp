/**
 * @file NASS_Vectors.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Vectors.hpp implementations.
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <print>
#include <cmath>
#include <cstring>

#include "../include/Vectors.hpp"

namespace nass {
    namespace internal {
        
        /**
         * @brief Copy a real_t* into a real_t*, both of size N0.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         */
        void Cp_RvtRvN_0(real_t* Rvt0, const real_t* Rv0, const natural_t& N0) {
            #ifdef _NEON

            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_0, Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0));
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_1, Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1));
            }

            for(natural_t N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1)
                Rvt0[N1] = Rv0[N1];

            #else
            #ifdef _OPENMP

            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0 - 1; N1 += 2) {
                Rvt0[N1] = Rv0[N1];
                Rvt0[N1 + 1] = Rv0[N1 + 1];
            }

            for(natural_t N1 = N0 - (N0 % 2); N1 < N0; ++N1)
                Rvt0[N1] = Rv0[N1];

            #else
            memcpy(Rvt0, Rv0, N0 * sizeof(real_t));
            #endif
            #endif
        }


        /**
         * @brief Dot product between two real_t*, both of size N0.
         * 
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         * @param N0 Natural number [N].
         * @return real_t Real number [R].
         */
        real_t Dt_RvRvN_R(const real_t* Rv0, const real_t* Rv1, const natural_t& N0) {
            #ifdef _NEON
            
            reals_t Rs0 = Ex_R_Rs(0.0), Rs1 = Ex_R_Rs(0.0);
            
            #pragma omp parallel for reduction(Rd_Rs: Rs0, Rs1)
            for(natural_t N1 = 0; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                const reals_t Rs00 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs01 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1);

                const reals_t Rs10 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs11 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_1);

                Rs0 = Ad_RsRs_Rs(Rs0, Ml_RsRs_Rs(Rs00, Rs10));
                Rs1 = Ad_RsRs_Rs(Rs1, Ml_RsRs_Rs(Rs01, Rs11));
            }

            real_t R0 = Rd_Rs_R(Rs0) + Rd_Rs_R(Rs1);

            for(natural_t N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1)
                R0 += Rv0[N1] * Rv1[N1];

            return R0;

            #else

            real_t R0 = 0.0, R1 = 0.0;

            #pragma omp parallel for reduction(+: R0, R1)
            for(natural_t N1 = 0; N1 < N0 - 1; N1 += 2) {
                R0 += Rv0[N1] * Rv1[N1];
                R1 += Rv0[N1 + 1] * Rv1[N1 + 1];
            }

            for(natural_t N1 = N0 - (N0 % 2); N1 < N0; ++N1)
                R0 += Rv0[N1] * Rv1[N1];

            return R0 + R1;

            #endif
        }


        /**
         * @brief Non-parallel dot product between two real_t*, both of size N0.
         * 
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         * @param N0 Natural number [N].
         * @return real_t Real number [R].
         */
        real_t NPDt_RvRvN_R(const real_t* Rv0, const real_t* Rv1, const natural_t& N0) {
            natural_t N1 = 0;

            #ifdef _NEON
            
            reals_t Rs0 = Ex_R_Rs(0.0), Rs1 = Ex_R_Rs(0.0);

            for(; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                const reals_t Rs00 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs01 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1);

                const reals_t Rs10 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs11 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_1);

                Rs0 = Ad_RsRs_Rs(Rs0, Ml_RsRs_Rs(Rs00, Rs10));
                Rs1 = Ad_RsRs_Rs(Rs1, Ml_RsRs_Rs(Rs01, Rs11));
            }

            real_t R0 = Rd_Rs_R(Rs0) + Rd_Rs_R(Rs1);

            for(; N1 < N0; ++N1)
                R0 += Rv0[N1] * Rv1[N1];

            return R0;

            #else

            real_t R0 = 0.0, R1 = 0.0;

            for(; N1 < N0 - 1; N1 += 2) {
                R0 += Rv0[N1] * Rv1[N1];
                R1 += Rv0[N1 + 1] * Rv1[N1 + 1];
            }

            for(; N1 < N0; ++N1)
                R0 += Rv0[N1] * Rv1[N1];

            return R0 + R1;

            #endif
        }


        /**
         * @brief Norm of a real_t*.
         * 
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         * @return real_t 
         */
        real_t Nr_RvN_R(const real_t* Rv0, const natural_t& N0) {
            #ifdef _NEON
            
            reals_t Rs0 = Ex_R_Rs(0.0), Rs1 = Ex_R_Rs(0.0);
            
            #pragma omp parallel for reduction(Rd_Rs: Rs0, Rs1)
            for(natural_t N1 = 0; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                const reals_t Rs2 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs3 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1);

                Rs0 = Ad_RsRs_Rs(Rs0, Ml_RsRs_Rs(Rs2, Rs2));
                Rs1 = Ad_RsRs_Rs(Rs1, Ml_RsRs_Rs(Rs3, Rs3));
            }

            real_t R0 = Rd_Rs_R(Rs0) + Rd_Rs_R(Rs1);

            for(natural_t N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1) {
                const real_t R3 = Rv0[N1];

                R0 += R3 * R3;
            }

            return std::sqrt(R0);

            #else

            real_t R0 = 0.0, R1 = 0.0;

            #pragma omp parallel for reduction(+: R0, R1)
            for(natural_t N1 = 0; N1 < N0 - 1; N1 += 2) {
                const real_t R3 = Rv0[N1], R4 = Rv0[N1 + 1];

                R0 += R3 * R3;
                R1 += R4 * R4;
            }

            for(natural_t N1 = N0 - (N0 % 2); N1 < N0; ++N1) {
                const real_t R3 = Rv0[N1];

                R0 += R3 * R3;
            }

            return std::sqrt(R0 + R1);

            #endif
        }


        /**
         * @brief Normalization of a real_t*.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param N0 Natural number [N].
         */
        void Nrz_RvtN_0(real_t* Rvt0, const natural_t& N0) {
            const real_t R0 = Nr_RvN_R(Rvt0, N0);

            #ifdef _NEON
            
            const reals_t Rs0 = Ex_R_Rs(R0);
            
            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                const reals_t Rs1 = Ld_Rv_Rs(Rvt0 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs2 = Ld_Rv_Rs(Rvt0 + N1 + MEMORY_OFFSET_1);

                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_0, Dv_RsRs_Rs(Rs1, Rs0));
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_1, Dv_RsRs_Rs(Rs2, Rs0));
            }

            for(natural_t N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1)
                Rvt0[N1] /= R0;

            #else

            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0 - 1; N1 += 2) {
                Rvt0[N1] /= R0;
                Rvt0[N1 + 1] /= R0;
            }

            for(natural_t N1 = N0 - (N0 % 2); N1 < N0; ++N1)
                Rvt0[N1] /= R0;

            #endif
        }


        /**
         * @brief Projection onto a real_t*.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         * @param N0 Natural number [N].
         */
        void Prj_RvtRvRvN_0(real_t* Rvt0, const real_t* Rv0, const real_t* Rv1, const natural_t& N0) {
            const real_t R0 = Dt_RvRvN_R(Rv0, Rv1, N0);

            #ifdef _NEON
            
            const reals_t Rs0 = Ex_R_Rs(R0);
            
            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                const reals_t Rs1 = Ld_Rv_Rs(Rvt0 + N1 + MEMORY_OFFSET_0);
                const reals_t Rs2 = Ld_Rv_Rs(Rvt0 + N1 + MEMORY_OFFSET_1);

                const reals_t Rs3 = Ml_RsRs_Rs(Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0), Rs0);
                const reals_t Rs4 = Ml_RsRs_Rs(Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1), Rs0);

                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_0, Sb_RsRs_Rs(Rs1, Rs3));
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_1, Sb_RsRs_Rs(Rs2, Rs4));
            }

            for(natural_t N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1) {
                Rvt0[N1] -= R0 * Rv0[N1];
            }

            #else

            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0 - 1; N1 += 2) {
                Rvt0[N1] -= R0 * Rv0[N1];
                Rvt0[N1 + 1] -= R0 * Rv0[N1 + 1];
            }

            for(natural_t N1 = N0 - (N0 % 2); N1 < N0; ++N1) {
                Rvt0[N1] -= R0 * Rv0[N1];
            }

            #endif
        }


        /**
         * @brief Prints a (row) real_t*.
         * 
         * @param Rrv0 Real row vector [Rrv].
         * @param N0 Natural number [N].
         */
        void Pr_RrvN_0(const real_t* Rrv0, const natural_t& N0) {
            std::println("--- Row vector.");

            for(natural_t N1 = 0; N1 < N0; ++N1) {
                Pr_R_0(Rrv0[N1]);
            }

            std::println("\n---");
        }


        /**
         * @brief Prints a (column) real_t*.
         * 
         * @param Rcv0 Real column vector [Rcv].
         * @param N0 Natural number [N].
         */
        void Pr_RcvN_0(const real_t* Rcv0, const natural_t& N0) {
            std::println("--- Column vector.");

            for(natural_t N1 = 0; N1 < N0; ++N1) {
                Pn_R_0(Rcv0[N1]);
            }

            std::println("---");
        }

    }
}