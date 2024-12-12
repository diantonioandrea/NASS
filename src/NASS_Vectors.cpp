/**
 * @file NASS_Vectors.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Vectors.hpp implementations.
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <cassert>
#include <iostream>

#include "../include/Vectors.hpp"

#ifdef _OPENMP
namespace nla {
    namespace internal {
        void Rd_RstRs_0(reals_t &Rst0, const reals_t &Rs0) { Rst0 = Ad_RsRs_Rs(Rst0, Rs0); }
        void In_Rst_0(reals_t &Rst0) { Rst0 = Ex_R_Rs(0.0); }

        #pragma omp declare reduction(Rd_Rs: reals_t: Rd_RstRs_0(omp_out, omp_in)) initializer(In_Rst_0(omp_out))
    }
}
#endif

namespace nla {
    namespace internal {
        
        /**
         * @brief Copy a real_t* into a real_t*, both of size N0.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         */
        void Cp_RvtRvN_0(real_t* Rvt0, const real_t* Rv0, const natural_t& N0) {
            #ifdef __ARM_NEON

            natural_t N1;

            #pragma omp parallel for
            for(N1 = 0; N1 < N0 - LOOP_OFFSET + 1; N1 += LOOP_OFFSET) {
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_0, Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0));
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_1, Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1));
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_2, Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_2));
                St_RvtRs_0(Rvt0 + N1 + MEMORY_OFFSET_3, Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_3));
            }

            for(N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1)
                Rvt0[N1] = Rv0[N1];

            #else
            memcpy(Rvt0, Rv0, N0 * sizeof(real_t));
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
            #ifndef NDEBUG // Integrity check.
            assert(N0 > LOOP_OFFSET);
            #endif

            natural_t N1;

            #ifdef __ARM_NEON
            if(N0 > LOOP_OFFSET) {
                reals_t Rs0 = Ex_R_Rs(0.0);
                reals_t Rs1 = Ex_R_Rs(0.0);
                reals_t Rs2 = Ex_R_Rs(0.0);
                reals_t Rs3 = Ex_R_Rs(0.0);
                
                #pragma omp parallel for reduction(Rd_Rs: Rs0, Rs1, Rs2, Rs3)
                for(N1 = 0; N1 < N0 - MEMORY_OFFSET_3 + 1; N1 += LOOP_OFFSET) {
                    const reals_t Rs00 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_0);
                    const reals_t Rs01 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_1);
                    const reals_t Rs02 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_2);
                    const reals_t Rs03 = Ld_Rv_Rs(Rv0 + N1 + MEMORY_OFFSET_3);

                    const reals_t Rs10 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_0);
                    const reals_t Rs11 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_1);
                    const reals_t Rs12 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_2);
                    const reals_t Rs13 = Ld_Rv_Rs(Rv1 + N1 + MEMORY_OFFSET_3);

                    const reals_t Rs20 = Ml_RsRs_Rs(Rs00, Rs10);
                    const reals_t Rs21 = Ml_RsRs_Rs(Rs01, Rs11);
                    const reals_t Rs22 = Ml_RsRs_Rs(Rs02, Rs12);
                    const reals_t Rs23 = Ml_RsRs_Rs(Rs03, Rs13);

                    Rs0 = Ad_RsRs_Rs(Rs0, Rs20);
                    Rs1 = Ad_RsRs_Rs(Rs1, Rs21);
                    Rs2 = Ad_RsRs_Rs(Rs2, Rs22);
                    Rs3 = Ad_RsRs_Rs(Rs3, Rs23);
                }

                real_t R0 = Rd_Rs_R(Rs0) + Rd_Rs_R(Rs1) + Rd_Rs_R(Rs2) + Rd_Rs_R(Rs3);

                for(N1 = N0 - (N0 % LOOP_OFFSET); N1 < N0; ++N1)
                    R0 += Rv0[N1] * Rv1[N1];

                return R0;
            }
            #endif

            real_t R0 = 0.0;
            real_t R1 = 0.0;
            real_t R2 = 0.0;
            real_t R3 = 0.0;

            #pragma omp parallel for reduction(+: R0, R1, R2, R3)
            for(N1 = 0; N1 < N0 - 3; N1 += 4) {
                R0 += Rv0[N1 + 0] * Rv1[N1 + 0]; // Superfluous.
                R1 += Rv0[N1 + 1] * Rv1[N1 + 1];
                R2 += Rv0[N1 + 2] * Rv1[N1 + 2];
                R3 += Rv0[N1 + 3] * Rv1[N1 + 3];
            }

            for(N1 = N0 - (N0 % 4); N1 < N0; ++N1)
                R0 += Rv0[N1 + 0] * Rv1[N1 + 0];

            return R0 + R1 + R2 + R3;
        }

    }
}