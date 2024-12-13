/**
 * @file NASS_Solver.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief 
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../include/Sparse.hpp"
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

            // Embedding.
            const natural_t N3 = 2 * (N1 + 1);
            const auto [Nv2, Nv3, Rv2] = Sec_NN_NvNvRv(N1, N0);

            // Basis.
            real_t* Mv0 = new real_t[N0 * N1];

            // LS matrix.
            real_t* Mv1 = new real_t[N0 * N1];

            // Sketch matrix.
            real_t* Mv2 = new real_t[N3 * N1];

            // (Thin) QR matrices.
            real_t* Mv3 = new real_t[N3 * N1];
            real_t* Mv4 = new real_t[N1 * N1];

            // Minimizer.
            real_t* Rv3 = new real_t[N1];

            // Residual and residual sketch.
            real_t* Rv4 = new real_t[N0];
            real_t* Rv5 = new real_t[N3];


            // Residual.
            // [!]

            // Residual sketch.
            Mlc_RvtNNvNvRvRv_0(Rv5, N0, Nv2, Nv3, Rv2, Rv4);

            // First basis column.
            Cp_RvtRvN_0(Mv0, Rv4, N0);
            Nrz_RvN_0(Mv0, N0);

            // First LS column.
            Mlc_RvtNNvNvRvRv_0(Mv1, N0, Nv0, Nv1, Rv0, Mv0);


            // ...


            // Clean-up.
            delete[] Nv2; delete[] Nv3; delete[] Rv2;
            delete[] Mv0;
            delete[] Mv1;
            delete[] Mv2;
            delete[] Mv3; delete[] Mv4;
            delete[] Rv3;
            delete[] Rv4; delete[] Rv5;
        }

    }
}