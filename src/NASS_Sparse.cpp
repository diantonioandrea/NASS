/**
 * @file NASS_Sparse.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Sparse.hpp implementations.
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <cstdlib>
#include <ctime>
#include <cmath>
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
        [[nodiscard]] std::tuple<natural_t, natural_t*, natural_t*, real_t*> Spc_St_NNvNvRv(const std::string& St0) {
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

            // File closing.
            F0.close();

            return {N0, Nv0, Nv1, Rv0};
        }


        /**
         * @brief Loads a (CSR) sparse matrix from a file.
         * 
         * @param St String [St]
         * @return std::tuple<natural_t, natural_t*, natural_t*, real_t*> (CSC) sparse matrix [Spc].
         */
        [[nodiscard]] std::tuple<natural_t, natural_t*, natural_t*, real_t*> Spr_St_NNvNvRv(const std::string& St0) {
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
            natural_t* Nv0 = new natural_t[N0 + 1];
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

                for(; N3 < N1; ++N3)
                    Nv0[N3 + 1] = N4;

                Nv1[N4] = N2;
                Rv0[N4] = R0;

                ++N4;
            }

            // Last entries.
            for(; N3 < N0; ++N3)
                Nv0[N3 + 1] = N4;

            // File closing.
            F0.close();

            return {N0, Nv0, Nv1, Rv0};
        }


        /**
         * @brief Multiplies a (CSC) sparse matrix by a real_t*.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param N0 Natural number [N].
         * @param Nv0 Natural vector [Nv].
         * @param Nv1 Natural vector [Nv].
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         */
        void Mlc_RvtNNvNvRvRv_0(real_t* Rvt0, const natural_t& N0, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rv1) {
            for(natural_t N1 = 0; N1 < N0; ++N1)
                for(natural_t N2 = Nv0[N1]; N2 < Nv0[N1 + 1]; ++N2)
                    Rvt0[Nv1[N2]] += Rv0[N2] * Rv1[N1];
        }


        /**
         * @brief Multiplies a (CSR) sparse matrix by a real_t*.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param N0 Natural number [N].
         * @param Nv0 Natural vector [Nv].
         * @param Nv1 Natural vector [Nv].
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         */
        void Mlr_RvtNNvNvRvRv_0(real_t* Rvt0, const natural_t& N0, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rv1) {

            #pragma omp parallel for
            for(natural_t N1 = 0; N1 < N0; ++N1)
                for(natural_t N2 = Nv0[N1]; N2 < Nv0[N1 + 1]; ++N2)
                    Rvt0[N1] += Rv0[N2] * Rv1[Nv1[N2]];
        }


        /**
         * @brief Evaluates the residual of a (CSC) sparse linear system.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param N0 Natural number [N].
         * @param Nv0 Natural vector [Nv].
         * @param Nv1 Natural vector [Nv].
         * @param Rv0 Real vector [Rv].
         * @param Rv1 Real vector [Rv].
         * @param Rv2 Real vector [Rv].
         */
        void RMlc_RvtNNvNvRvRvRv_0(real_t* Rvt0, const natural_t& N0, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rv1, const real_t* Rv2) {
            Cp_RvtRvN_0(Rvt0, Rv2, N0);

            for(natural_t N1 = 0; N1 < N0; ++N1)
                for(natural_t N2 = Nv0[N1]; N2 < Nv0[N1 + 1]; ++N2)
                    Rvt0[Nv1[N2]] -= Rv0[N2] * Rv1[N1];
        }


        /**
         * @brief Multiplies a (CSC) sparse matrix by a real_t*.
         * 
         * @param Rmt0 Real matrix [Rm], target [t].
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @param Nv0 Natural vector [Nv].
         * @param Nv1 Natural vector [Nv].
         * @param Rv0 Real vector [Rv].
         * @param Rm0 Real matrix [Rm].
         * @param N2 Natural number [N].
         */
        void Mlc_RmtNNNvNvRvRmN_0(real_t* Rmt0, const natural_t& N0, const natural_t& N1, const natural_t* Nv0, const natural_t* Nv1, const real_t* Rv0, const real_t* Rm0, const natural_t& N2) {

            #pragma omp parallel for
            for(natural_t N3 = 0; N3 < N2; ++N3)
                Mlc_RvtNNvNvRvRv_0(Rmt0 + N3 * N0, N1, Nv0, Nv1, Rv0, Rm0 + N3 * N1);
        }


        /**
         * @brief (CSC) sparse embedding.
         * 
         * @param N0 Natural number [N].
         * @param N1 Natural number [N].
         * @return std::tuple<natural_t*, natural_t*, real_t*> 
         */
        [[nodiscard]] std::tuple<natural_t*, natural_t*, real_t*> Sec_NN_NvNvRv(const natural_t& N0, const natural_t& N1) {
            std::srand(std::time(nullptr));
            
            const natural_t N2 = 2 * (N0 + 1);
            const natural_t N3 = static_cast<natural_t>(std::ceil(2.0 * std::log(N0 + 1.0)));
            const real_t R0 = 2.0 / std::sqrt(static_cast<real_t>(N2));

            natural_t* Nv0 = new natural_t[N1 + 1];
            natural_t* Nv1 = new natural_t[N3 * N1];
            real_t* Rv0 = new real_t[N3 * N1];
            
            #pragma omp parallel for
            for(natural_t N4 = 0; N4 < N1 + 1; ++N4)
                Nv0[N4] = N4 * N3;

            #pragma omp parallel for
            for(natural_t N4 = 0; N4 < N3 * N1; ++N4) {
                Nv1[N4] = std::rand() % N2;
                Rv0[N4] = (static_cast<real_t>(std::rand()) / RAND_MAX - 0.5) * R0;
            }

            return {Nv0, Nv1, Rv0};
        }

    }
}