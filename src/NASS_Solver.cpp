/**
 * @file NASS_Solver.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../include/Matrix.hpp"
#include "../include/Sparse.hpp"
#include "../include/Decomposition.hpp"
#include "../include/Solver.hpp"

namespace nass {
    namespace internal {
        
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
        void sGMRES_RvNNvNvRvRvNN_0(real_t* Rvt0, const natural_t& N0, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rv1, const natural_t& N1, const natural_t& N2) {

            // PREPARATION.


            // Embedding.
            const natural_t N3 = 2 * (N1 + 1);
            const auto [Nv2, Nv3, Rv2] = Sec_NN_NvNvRv(N1, N0);

            // Basis.
            real_t* Rm0 = new real_t[N0 * N1];

            // LS matrix.
            real_t* Rm1 = new real_t[N0 * N1];

            // Sketch matrix.
            real_t* Rm2 = new real_t[N3 * N1];

            // (Thin) QR matrices.
            real_t* Rm3 = new real_t[N3 * N3];

            // Minimizer.
            real_t* Rv3 = new real_t[N1];

            // Residual and residual sketch.
            real_t* Rv4 = new real_t[N0];
            real_t* Rv5 = new real_t[N3];

            // LS solution.
            real_t* Rv6 = new real_t[N1];

            
            // sGMRES.


            // Residual.
            RMlc_RvtNNvNvRvRvRv_0(Rv4, N0, Nv0, Nv1, Rv0, Rvt0, Rv1);

            // Residual sketch.
            Mlc_RvtNNvNvRvRv_0(Rv5, N0, Nv2, Nv3, Rv2, Rv4);

            // First basis column.
            Cp_RvtRvN_0(Rm0, Rv4, N0);
            Nrz_RvtN_0(Rm0, N0);

            // First LS column.
            Mlc_RvtNNvNvRvRv_0(Rm1, N0, Nv0, Nv1, Rv0, Rm0);

            // Truncated Arnoldi, first columns.
            for(natural_t N4 = 1; N4 <= N2; ++N4) {
                
                // Copy.
                Cp_RvtRvN_0(Rm0 + N4 * N0, Rm1 + (N4 - 1) * N0, N0);

                // Orthogonalization.
                for(natural_t N5 = 1; N5 <= N4; ++N5)
                    Prj_RvtRvRvN_0(Rm0 + N4 * N0, Rm0 + (N4 - N5) * N0, Rm1 + (N4 - 1) * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm0 + N4 * N0, N0);

                // LS matrix.
                Mlc_RvtNNvNvRvRv_0(Rm1 + N4 * N0, N0, Nv0, Nv1, Rv0, Rm0 + N4 * N0);
            }
            
            // Truncated Arnoldi, last columns.
            for(natural_t N4 = N2 + 1; N4 < N1; ++N4) {
                
                // Copy.
                Cp_RvtRvN_0(Rm0 + N4 * N0, Rm0 + (N4 - 1) * N0, N0);

                // Orthogonalization.
                for(natural_t N5 = 1; N5 <= N2; ++N5)
                    Prj_RvtRvRvN_0(Rm0 + N4 * N0, Rm0 + (N4 - N5) * N0, Rm1 + (N4 - 1) * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm0 + N4 * N0, N0);

                // LS matrix.
                Mlc_RvtNNvNvRvRv_0(Rm1 + N4 * N0, N0, Nv0, Nv1, Rv0, Rm0 + N4 * N0);
            }

            // Sketch matrix.
            Mlc_RmtNNvNvRvRmN_0(Rm2, N3, Nv2, Nv3, Rv2, Rm1, N1);

            // (Thin) QR.
            TQR_RmtRmtRvtNN_0(Rm3, Rm2, Rv5, N3, N1);

            // (Reduced) LS problem, backward substitution.
            for(natural_t N4 = N1; N4 > 0; --N4) {
                real_t R0 = 0.0;

                for(natural_t N5 = N4; N5 < N1; ++N5)
                    R0 += Rm2[N5 * N3 + N4 - 1] * Rv6[N5];
                
                Rv6[N4 - 1] = (Rv5[N4 - 1] - R0) / Rm2[(N4 - 1) * (N3 + 1)];
            }

            // Residual estimation.
            // [!]

            // Solution update.
            Ml_RvtRmRvNN_0(Rvt0, Rm1, Rv6, N0, N1);


            // CLEAN-UP AND RETURN.


            delete[] Nv2; delete[] Nv3; delete[] Rv2;
            delete[] Rm0;
            delete[] Rm1;
            delete[] Rm2;
            delete[] Rm3;
            delete[] Rv3;
            delete[] Rv4; delete[] Rv5;
        }

    }
}