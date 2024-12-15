/**
 * @file NASS_Solver.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief
 * @date 2024-12-13
 *
 * @copyright Copyright (c) 2024
 *
 */

#if !defined(SPARSE_SKETCH) && !defined(GAUSS_SKETCH)
#define SPARSE_SKETCH
#endif

#include "../include/Vectors.hpp"
#include "../include/Matrix.hpp"
#include "../include/Sparse.hpp"
#include "../include/Decomposition.hpp"
#include "../include/Solver.hpp"

namespace nass {
    namespace internal {

        /**
         * @brief Residual estimate.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rm0 Real matrix [Rm].
         * @param Rv0 Real vector [Rv].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         */
        void REst_RvtRmRvNN_0(real_t* Rvt0, const real_t* Rm0, const real_t* Rv0, const natural_t& N0, const natural_t& N1) {
            for(natural_t N2 = 0; N2 < N1; ++N2)
                for(natural_t N3 = 0; N3 < N0; ++N3)
                    Rvt0[N3] -= Rm0[N2 * N0 + N3] * Rv0[N2];
        }


        /**
         * @brief Sketched GMRES.
         *
         * @param Rvt0 Real vector [Rv], target [t].
         * @param N0 Natural number [N].
         * @param Nv0 Natural vector [Nv].
         * @param Nv1 Natural vector [Nv].
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         * @param N1 Natural number [N].
         * @param N2 Natural number [N].
         */
        real_t sGMRES_RvNNvNvRvRvNN_0(real_t* Rvt0, const natural_t& N0, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rv1, const natural_t& N1, const natural_t& N2) {

            #ifndef NVERBOSE
            std::println("--- sGMRES.");
            std::println("Parameters: {}, {}, {}", N0, N1, N2);
            std::println("---");
            #endif


            // PREPARATION.


            // Embedding.
            const natural_t N3 = 2 * (N1 + 1);

            #if defined(SPARSE_SKETCH)
            const auto [Nv2, Nv3, Rv2] = Sec_NN_NvNvRv(N1, N0);
            #elif defined(GAUSS_SKETCH)
            real_t* Rm0 = Gs_NN_Rm(N3, N0);
            #endif

            // Basis.
            real_t* Rm1 = new real_t[N0 * N1];

            // LS matrix.
            real_t* Rm2 = new real_t[N0 * N1];

            // Sketch matrix.
            real_t* Rm3 = new real_t[N3 * N1];

            // QR.
            real_t* Rm4 = new real_t[N3 * N3];
            natural_t* Nv4 = new natural_t[N1];

            // Minimizer.
            real_t* Rv3 = new real_t[N1];

            // Residual and residual sketch.
            real_t* Rv4 = new real_t[N0];
            real_t* Rv5 = new real_t[N3];

            // LS solution.
            real_t* Rv6 = new real_t[N1];

            // Residual estimates.
            real_t* Rv7 = new real_t[N3];


            // sGMRES.

            // Residual.
            RMlc_RvtNNvNvRvRvRv_0(Rv4, N0, Nv0, Nv1, Rv0, Rvt0, Rv1);

            // Residual sketch.
            #if defined(SPARSE_SKETCH)
            Mlc_RvtNNvNvRvRv_0(Rv5, N0, Nv2, Nv3, Rv2, Rv4);
            #elif defined(GAUSS_SKETCH)
            Ml_RvtRmRvNN_0(Rv5, Rm0, Rv4, N3, N0);
            #endif

            Cp_RvtRvN_0(Rv7, Rv5, N3);

            // First basis column.
            Cp_RvtRvN_0(Rm1, Rv4, N0);
            Nrz_RvtN_0(Rm1, N0);

            // First LS column.
            Mlc_RvtNNvNvRvRv_0(Rm2, N0, Nv0, Nv1, Rv0, Rm1);

            // Truncated Arnoldi.
            const natural_t N4 = N2 > N1 ? N1 : N2;

            // Truncated Arnoldi, first columns.
            for(natural_t N5 = 1; N5 < N4; ++N5) {

                // Copy.
                Cp_RvtRvN_0(Rm1 + N5 * N0, Rm2 + (N5 - 1) * N0, N0);

                // Orthogonalization.
                for(natural_t N6 = 1; N6 <= N5; ++N6)
                    Prj_RvtRvRvN_0(Rm1 + N5 * N0, Rm1 + (N5 - N6) * N0, Rm1 + N5 * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm1 + N5 * N0, N0);

                // LS matrix.
                Mlc_RvtNNvNvRvRv_0(Rm2 + N5 * N0, N0, Nv0, Nv1, Rv0, Rm1 + N5 * N0);
            }

            // Truncated Arnoldi, last columns.
            for(natural_t N5 = N4; N5 < N1; ++N5) {

                // Copy.
                Cp_RvtRvN_0(Rm1 + N5 * N0, Rm2 + (N5 - 1) * N0, N0);

                // Orthogonalization.
                for(natural_t N6 = 1; N6 <= N4; ++N6)
                    Prj_RvtRvRvN_0(Rm1 + N5 * N0, Rm1 + (N5 - N6) * N0, Rm1 + N5 * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm1 + N5 * N0, N0);

                // LS matrix.
                Mlc_RvtNNvNvRvRv_0(Rm2 + N5 * N0, N0, Nv0, Nv1, Rv0, Rm1 + N5 * N0);
            }

            // Sketch matrix.
            #if defined(SPARSE_SKETCH)
            Mlc_RmtNNNvNvRvRmN_0(Rm3, N3, N0, Nv2, Nv3, Rv2, Rm2, N1);
            #elif defined(GAUSS_SKETCH)
            Ml_RmtRmRmNNN_0(Rm3, Rm0, Rm2, N3, N0, N1);
            #endif

            // QR.
            TQR_RmtRmtNvtRvtNN_0(Rm4, Rm3, Nv4, Rv5, N3, N1);

            // (Reduced) LS problem, backward substitution.
            for(natural_t N5 = N1; N5 > 0; --N5) {
                real_t R0 = 0.0;

                for(natural_t N6 = N5; N6 < N1; ++N6)
                    R0 += Rm3[N6 * N3 + N5 - 1] * Rv6[N6];

                Rv6[N5 - 1] = (Rv5[N5 - 1] - R0) / Rm3[(N5 - 1) * (N3 + 1)];
            }

            // Pivoting.
            for(natural_t N5 = 0; N5 < N1; ++N5) {
                const real_t R0 = Rv6[N5];
                Rv6[N5] = Rv6[Nv4[N5]];
                Rv6[Nv4[N5]] = R0;
            }

            // Residual estimation.
            REst_RvtRmRvNN_0(Rv7, Rm4, Rv5, N3, N1);
            const real_t R0 = Nr_RvN_R(Rv7, N3);

            // Solution update.
            Ml_RvtRmRvNN_0(Rvt0, Rm1, Rv6, N0, N1);


            // CLEAN-UP AND RETURN.


            #if defined(SPARSE_SKETCH)
            delete[] Nv2; delete[] Nv3; delete[] Rv2;
            #elif defined(GAUSS_SKETCH)
            delete[] Rm0;
            #endif

            delete[] Rm1;
            delete[] Rm2;
            delete[] Rm3;
            delete[] Rm4; delete[] Nv4;
            delete[] Rv3;
            delete[] Rv4; delete[] Rv5;
            delete[] Rv6;
            delete[] Rv7;

            return R0;
        }

    }
}