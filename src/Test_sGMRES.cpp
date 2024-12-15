/**
 * @file Test_sGMRES.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief sGMRES testing.
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "./Test.hpp"

#include "../include/Sparse.hpp"
#include "../include/Solver.hpp"
using namespace nass;

int main(int argc, char** argv) {

    // Seeding.
    std::srand(std::time(nullptr));

    // Arguments.
    if(argc < 3) {
        std::println("Usage: {} St [String, path] N [Natural, Subspace dimension] N? [Natural, Arnoldi]", argv[0]);
        return -1;
    }

    // Parameters.
    const natural_t N1 = std::atoi(argv[2]);
    const natural_t N2 = argc > 3 ? std::atoi(argv[3]) : 4;

    // Sparse matrix.
    auto [N0, Nv0, Nv1, Rv0] = internal::Spc_St_NNvNvRv(argv[1]);

    // Solution.
    real_t* Rv1 = new real_t[N0];

    // Expected solution.
    real_t* Rv2 = new real_t[N0];

    // RHS.
    real_t* Rv3 = new real_t[N0];

    // Residual.
    real_t* Rv4 = new real_t[N0];

    // Expected solution.
    for(natural_t N3 = 0; N3 < N0; ++N3)
        Rv2[N3] = 1.0; // static_cast<real_t>(std::rand()) / RAND_MAX;

    // RHS.
    internal::Mlc_RvtNNvNvRvRv_0(Rv3, N0, Nv0, Nv1, Rv0, Rv2);

    // sGMRES.
    const real_t R0 = internal::sGMRES_RvNNvNvRvRvNN_0(Rv1, N0, Nv0, Nv1, Rv0, Rv3, N1, N2);

    // Residual.
    internal::RMlc_RvtNNvNvRvRvRv_0(Rv4, N0, Nv0, Nv1, Rv0, Rv1, Rv3);

    // Relative residual.
    const real_t R1 = internal::Nr_RvN_R(Rv4, N0) / internal::Nr_RvN_R(Rv3, N0);

    // Output.
    std::println("Residual, estimate: {:.3e}", R0);
    std::println("Residual, relative: {:.3e}", R1);

    // Clean-up.
    delete[] Nv0; delete[] Nv1; delete[] Rv0;
    delete[] Rv1;
    delete[] Rv3;
    delete[] Rv4;

    return 0;
}