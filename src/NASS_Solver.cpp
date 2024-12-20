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

#if defined(SPARSE_SKETCH) && defined(GAUSS_SKETCH)
#error "Unsafe constant definition."
#endif

#if defined(GAUSS_SKETCH) && defined(MEMORY_PRIORITY)
#error "Unsafe constant definition."
#endif

#ifndef NVERBOSE
#include <chrono>
#endif

#ifndef NDEBUG
#include <cassert>
#endif

#include "../include/Vectors.hpp"
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
         * @return std::array<real_t, 2> Real numbers [R].
         */
        std::array<real_t, 2> sGMRES_RvNNvNvRvRvNN_RR(real_t* Rvt0, const natural_t& N0, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rv1, const natural_t& N1, const natural_t& N2) {
            #ifndef NDEBUG // Integrity check.
            assert(N1 > 0);
            assert(N2 > 0);
            assert(N2 <= N1);
            #endif

            #ifndef NVERBOSE
            #ifndef MEMORY_PRIORITY
            std::println("--- sGMRES.");
            #else
            std::println("--- sGMRES, memory priority.");
            #endif
            std::println("Parameters: {}, {}, {}", N0, N1, N2);
            std::println("Timings:");
            #endif


            // PREPARATION.


            // Embedding.
            const natural_t N3 = 2 * (N1 + 1);


            #ifndef NVERBOSE
            auto T0 = std::chrono::high_resolution_clock::now();
            #endif


            #if defined(SPARSE_SKETCH)
            const auto [Nv2, Nv3, Rv2] = Sec_NN_NvNvRv(N1, N0);
            #elif defined(GAUSS_SKETCH)
            real_t* Rm0 = Gs_NN_Rm(N3, N0);
            #endif


            #ifndef NVERBOSE
            auto T1 = std::chrono::high_resolution_clock::now();

            std::println("\tSketch generation: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif


            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            #ifndef MEMORY_PRIORITY
            real_t* Rm1 = new real_t[N0 * N1]; // Basis.
            real_t* Rm2 = new real_t[N0 * N1]; // LS matrix.
            #else 
            real_t* Rm1 = new real_t[N0 * (N2 + 1)]; // Partial basis.
            #endif

            // Sketched LS matrix.
            real_t* Rm3 = new real_t[N3 * N1];

            // QR.
            real_t* Rm4 = new real_t[N3 * N1];
            natural_t* Nv4 = new natural_t[N1];

            // Minimizer.
            real_t* Rv3 = new real_t[N1];

            // Residual and residual sketch.
            real_t* Rv4 = new real_t[N0];
            real_t* Rv5 = new real_t[N3];

            // LS solution and permuted version.
            real_t* Rv6 = new real_t[N1];
            real_t* Rv7 = new real_t[N1];

            // Residual estimates.
            real_t* Rv8 = new real_t[N3];


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            std::println("\tAllocation: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif


            // sGMRES.

            // Residual.
            RMlc_RvtNNvNvRvRvRv_0(Rv4, N0, Nv0, Nv1, Rv0, Rvt0, Rv1);

            // Residual sketch.
            #if defined(SPARSE_SKETCH)
            Mlc_RvtNNvNvRvRv_0(Rv5, N0, Nv2, Nv3, Rv2, Rv4);
            #elif defined(GAUSS_SKETCH)
            Ml_RvtRmRvNN_0(Rv5, Rm0, Rv4, N3, N0);
            #endif

            Cp_RvtRvN_0(Rv8, Rv5, N3);

            // Arnoldi.


            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            // First basis column.
            Cp_RvtRvN_0(Rm1, Rv4, N0);
            Nrz_RvtN_0(Rm1, N0);

            #ifndef MEMORY_PRIORITY

            // First LS column.
            Mlc_RvtNNvNvRvRv_0(Rm2, N0, Nv0, Nv1, Rv0, Rm1);

            // Truncated Arnoldi, first part.
            for(natural_t N4 = 1; N4 <= N2; ++N4) {

                // Copy.
                Cp_RvtRvN_0(Rm1 + N4 * N0, Rm2 + (N4 - 1) * N0, N0);

                // Orthogonalization.
                for(natural_t N5 = 0; N5 < N4; ++N5)
                    Prj_RvtRvRvN_0(Rm1 + N4 * N0, Rm1 + N5 * N0, Rm1 + N4 * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm1 + N4 * N0, N0);

                // LS matrix.
                Mlc_RvtNNvNvRvRv_0(Rm2 + N4 * N0, N0, Nv0, Nv1, Rv0, Rm1 + N4 * N0);
            }

            // Truncated Arnoldi, second part.
            for(natural_t N4 = N2 + 1; N4 < N1; ++N4) {

                // Copy.
                Cp_RvtRvN_0(Rm1 + N4 * N0, Rm2 + (N4 - 1) * N0, N0);

                // Orthogonalization.
                for(natural_t N5 = N4 - N2; N5 < N4; ++N5)
                    Prj_RvtRvRvN_0(Rm1 + N4 * N0, Rm1 + N5 * N0, Rm1 + N4 * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm1 + N4 * N0, N0);

                // LS matrix.
                Mlc_RvtNNvNvRvRv_0(Rm2 + N4 * N0, N0, Nv0, Nv1, Rv0, Rm1 + N4 * N0);
            }

            #else

            // Truncated Arnoldi and sketching, first part.
            for(natural_t N4 = 1; N4 <= N2; ++N4) {

                // LS column.
                Mlc_RvtNNvNvRvRv_0(Rm1 + N4 * N0, N0, Nv0, Nv1, Rv0, Rm1 + (N4 - 1) * N0);

                // Sketch application.
                Mlc_RvtNNvNvRvRv_0(Rm3 + (N4 - 1) * N3, N0, Nv2, Nv3, Rv2, Rm1 + N4 * N0);

                // Orthogonalization.
                for(natural_t N5 = 0; N5 < N4; ++N5)
                    Prj_RvtRvRvN_0(Rm1 + N4 * N0, Rm1 + N5 * N0, Rm1 + N4 * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm1 + N4 * N0, N0);
            }

            // Truncated Arnoldi and sketching, second part.
            for(natural_t N4 = N2 + 1; N4 < N1; ++N4) {

                // Indices.
                const natural_t N5 = N4 % (N2 + 1);
                const natural_t N6 = (N4 - 1) % (N2 + 1);

                // Zeroing.
                #pragma omp parallel for
                for(natural_t N7 = 0; N7 < N0; ++N7)
                    Rm1[N5 * N0 + N7] = 0.0;

                // LS column.
                Mlc_RvtNNvNvRvRv_0(Rm1 + N5 * N0, N0, Nv0, Nv1, Rv0, Rm1 + N6 * N0);

                // Sketch application.
                Mlc_RvtNNvNvRvRv_0(Rm3 + (N4 - 1) * N3, N0, Nv2, Nv3, Rv2, Rm1 + N5 * N0);

                // Orthogonalization.
                for(natural_t N7 = 0; N7 <= N2; ++N7) {
                    if(N5 == N7)
                        continue;

                    Prj_RvtRvRvN_0(Rm1 + N5 * N0, Rm1 + N7 * N0, Rm1 + N5 * N0, N0);
                }

                // Normalization.
                Nrz_RvtN_0(Rm1 + N5 * N0, N0);
            }

            // Zeroing.
            #pragma omp parallel for
            for(natural_t N7 = 0; N7 < N0; ++N7)
                Rm1[(N1 % (N2 + 1)) * N0 + N7] = 0.0;

            // LS column.
            Mlc_RvtNNvNvRvRv_0(Rm1 + (N1 % (N2 + 1)) * N0, N0, Nv0, Nv1, Rv0, Rm1 + ((N1 - 1) % (N2 + 1)) * N0);

            // Sketch application.
            Mlc_RvtNNvNvRvRv_0(Rm3 + (N1 - 1) * N3, N0, Nv2, Nv3, Rv2, Rm1 + (N1 % (N2 + 1)) * N0);

            #endif


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            #ifndef MEMORY_PRIORITY
            std::println("\tArnoldi: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #else
            std::println("\tArnoldi (1) and sketch application: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif
            #endif


            #ifndef MEMORY_PRIORITY
            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            // Sketch matrix.
            #if defined(SPARSE_SKETCH)
            Mlc_RmtNNNvNvRvRmN_0(Rm3, N3, N0, Nv2, Nv3, Rv2, Rm2, N1);
            #elif defined(GAUSS_SKETCH)
            Ml_RmtRmRmNNN_0(Rm3, Rm0, Rm2, N3, N0, N1);
            #endif


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            std::println("\tSketch application: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif
            #endif


            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            // QR.
            TQR_RmtRmtNvtNN_0(Rm4, Rm3, Nv4, N3, N1);


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            std::println("\tQR decomposition: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif


            // Condition number estimate.
            real_t R0 = std::abs(Rm3[0]), R1 = std::abs(Rm3[0]);

            for(natural_t N4 = 1; N4 < N1; ++N4) {
                const real_t R2 = std::abs(Rm3[N4 * (N3 + 1)]);

                if(R2 < R0)
                    R0 = R2;

                if(R2 > R1)
                    R1 = R2;
            }


            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            // (Reduced) LS problem, Qt.
            Mqt_RvtRmNN_0(Rv5, Rm4, N3, N1, N1);

            // (Reduced) LS problem, backward substitution.
            for(natural_t N4 = N1; N4 > 0; --N4) {
                real_t R2 = 0.0;

                for(natural_t N5 = N4; N5 < N1; ++N5)
                    R2 += Rm3[N5 * N3 + N4 - 1] * Rv6[N5];

                Rv6[N4 - 1] = (Rv5[N4 - 1] - R2) / Rm3[(N4 - 1) * (N3 + 1)];
            }

            // Pivoting.
            for(natural_t N4 = 0; N4 < N1; ++N4)
                Rv7[Nv4[N4]] = Rv6[N4];

            // Residual estimation, Q.
            Mq_RvtRmNN_0(Rv5, Rm4, N3, N1, N1);

            // Residual estimation, subtraction.
            for(natural_t N4 = 0; N4 < N3; ++N4)
                Rv8[N4] -= Rv5[N4];

            // Residual estimation. [!] Q needs to be truncated.
            const real_t R2 = Nr_RvN_R(Rv8, N3);


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            std::println("\tLS problem: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif


            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            #ifndef MEMORY_PRIORITY

            // Solution update.
            Ml_RvtRmRvNN_0(Rvt0, Rm1, Rv7, N0, N1);

            #else 

            // First basis column.
            Cp_RvtRvN_0(Rm1, Rv4, N0);
            Nrz_RvtN_0(Rm1, N0);

            // Solution update.
            #pragma omp parallel for
            for(natural_t N5 = 0; N5 < N0; ++N5)
                Rvt0[N5] += Rm1[N5] * Rv7[0];

            // Truncated Arnoldi and solution update, first part.
            for(natural_t N4 = 1; N4 <= N2; ++N4) {

                // Zeroing.
                #pragma omp parallel for
                for(natural_t N5 = 0; N5 < N0; ++N5)
                    Rm1[N4 * N0 + N5] = 0.0;

                // LS column.
                Mlc_RvtNNvNvRvRv_0(Rm1 + N4 * N0, N0, Nv0, Nv1, Rv0, Rm1 + (N4 - 1) * N0);

                // Orthogonalization.
                for(natural_t N5 = 0; N5 < N4; ++N5)
                    Prj_RvtRvRvN_0(Rm1 + N4 * N0, Rm1 + N5 * N0, Rm1 + N4 * N0, N0);

                // Normalization.
                Nrz_RvtN_0(Rm1 + N4 * N0, N0);

                // Solution update.
                const real_t R3 = Rv7[N4];

                #pragma omp parallel for
                for(natural_t N5 = 0; N5 < N0; ++N5)
                    Rvt0[N5] += Rm1[N4 * N0 + N5] * R3;
            }

            // Truncated Arnoldi, second part.
            for(natural_t N4 = N2 + 1; N4 < N1; ++N4) {

                // Indices.
                const natural_t N5 = N4 % (N2 + 1);
                const natural_t N6 = (N4 - 1) % (N2 + 1);

                // Zeroing.
                #pragma omp parallel for
                for(natural_t N7 = 0; N7 < N0; ++N7)
                    Rm1[N5 * N0 + N7] = 0.0;

                // LS column.
                Mlc_RvtNNvNvRvRv_0(Rm1 + N5 * N0, N0, Nv0, Nv1, Rv0, Rm1 + N6 * N0);

                // Orthogonalization.
                for(natural_t N7 = 0; N7 <= N2; ++N7) {
                    if(N5 == N7)
                        continue;

                    Prj_RvtRvRvN_0(Rm1 + N5 * N0, Rm1 + N7 * N0, Rm1 + N5 * N0, N0);
                }

                // Normalization.
                Nrz_RvtN_0(Rm1 + N5 * N0, N0);

                // Solution update.
                const real_t R3 = Rv7[N4];

                #pragma omp parallel for
                for(natural_t N7 = 0; N7 < N0; ++N7)
                    Rvt0[N7] += Rm1[N5 * N0 + N7] * R3;
            }

            #endif


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            #ifndef MEMORY_PRIORITY
            std::println("\tSolution update: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #else
            std::println("\tArnoldi (2) and solution update: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif
            #endif


            // CLEAN-UP AND RETURN.


            #ifndef NVERBOSE
            T0 = std::chrono::high_resolution_clock::now();
            #endif


            #if defined(SPARSE_SKETCH)
            delete[] Nv2; delete[] Nv3; delete[] Rv2;
            #elif defined(GAUSS_SKETCH)
            delete[] Rm0;
            #endif

            delete[] Rm1;

            #ifndef MEMORY_PRIORITY
            delete[] Rm2;
            #endif

            delete[] Rm3;
            delete[] Rm4; delete[] Nv4;
            delete[] Rv3;
            delete[] Rv4; delete[] Rv5;
            delete[] Rv6;
            delete[] Rv7;
            delete[] Rv8;


            #ifndef NVERBOSE
            T1 = std::chrono::high_resolution_clock::now();

            std::println("\tDeallocation: {}", std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0));
            #endif


            #ifndef NVERBOSE
            std::println("---");
            #endif

            return {R2, R1 / R0};
        }

    }
}