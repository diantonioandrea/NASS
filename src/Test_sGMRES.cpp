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

    // Parameters.
    const natural_t N1 = 20; // Subspace size.
    const natural_t N2 = 4; // Arnoldi truncation.

    // Arguments.
    if(argc < 2) {
        std::println("Usage: {} St [String, path]", argv[0]);
        return -1;
    }

    // Sparse matrix.
    auto [N0, Nv0, Nv1, Rv0] = internal::Spc_St_NNvNvRv(argv[1]);

    // Solution.
    real_t* Rv1 = new real_t[N0];

    // RHS.
    real_t* Rv2 = new real_t[N0];

    // Residual.
    real_t* Rv3 = new real_t[N0];

    // Random RHS filling.
    for(natural_t N1 = 0; N1 < N0; ++N1)
        Rv2[N1] = static_cast<real_t>(std::rand()) / RAND_MAX;

    // sGMRES.
    internal::sGMRES_RvNNvNvRvRvNN_0(Rv1, N0, Nv0, Nv1, Rv0, Rv2, N1, N2);

    // Residual.
    internal::RMlc_RvtNNvNvRvRvRv_0(Rv3, N0, Nv0, Nv1, Rv0, Rv1, Rv2);

    // Relative residual.
    const real_t R0 = internal::Nr_RvN_R(Rv3, N0) / internal::Nr_RvN_R(Rv2, N0);

    // Output.
    std::println("Relative residual: {:.3e}", R0);

    // Clean-up.
    delete[] Nv0; delete[] Nv1; delete[] Rv0;
    delete[] Rv1;
    delete[] Rv2;
    delete[] Rv3;

    return 0;
}