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


    // TIMED.


    // Start.
    auto T0 = high_resolution_clock::now();

    // Sparse matrix.
    auto [N0, Nv0, Nv1, Rv0] = internal::Spc_St_NNvNvRv(argv[1]);

    // End.
    auto T1 = high_resolution_clock::now();


    // TIMED.


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


    // TIMED.


    // Start.
    auto T2 = high_resolution_clock::now();

    // sGMRES.
    const auto [R0, R1] = internal::sGMRES_RvNNvNvRvRvNN_RR(Rv1, N0, Nv0, Nv1, Rv0, Rv3, N1, N2);

    // End.
    auto T3 = high_resolution_clock::now();


    // TIMED.


    // Durations.
    auto D0 = duration_cast<seconds>(T1 - T0);
    auto D1 = duration_cast<seconds>(T3 - T2);

    // Residual.
    internal::RMlc_RvtNNvNvRvRvRv_0(Rv4, N0, Nv0, Nv1, Rv0, Rv1, Rv3);

    // Relative residual.
    const real_t R2 = internal::Nr_RvN_R(Rv4, N0) / internal::Nr_RvN_R(Rv3, N0);

    // Output.
    std::println("Residual, estimate:  {:.3e}", R0);
    std::println("Condition, estimate: {:.3e}\n", R1);
    std::println("Residual, relative:  {:.3e}\n", R2);
    std::println("Elapsed:\n\tLoading: {}\n\tsGMRES: {}", D0, D1);

    // Clean-up.
    delete[] Nv0; delete[] Nv1; delete[] Rv0;
    delete[] Rv1;
    delete[] Rv3;
    delete[] Rv4;

    return 0;
}