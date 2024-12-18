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
    for(natural_t N3 = 0; N3 < N0; ++N3) {
        const real_t R0 = static_cast<real_t>(std::rand()) / RAND_MAX;

        Rv2[N3] = std::sqrt(-2.0 * std::log(R0)) * std::cos(2.0 * M_PI * R0);
    }

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
    auto D0 = duration_cast<milliseconds>(T1 - T0);
    auto D1 = duration_cast<milliseconds>(T3 - T2);

    // Residual.
    internal::RMlc_RvtNNvNvRvRvRv_0(Rv4, N0, Nv0, Nv1, Rv0, Rv1, Rv3);

    // Relative residual.
    const real_t R2 = internal::Nr_RvN_R(Rv4, N0);
    const real_t R3 = R2 / internal::Nr_RvN_R(Rv3, N0);

    // Output.
    std::println("--- sGMRES testing.");
    std::println("Results:\n\tResidual: {:.3e}\n\tResidual, relative: {:.3e}", R2, R3);
    std::println("Estimates:\n\tResidual: {:.3e}\n\tCondition: {:.3e}", R0, R1);
    std::println("Timings:\n\tLoading: {}\n\tsGMRES: {}", D0, D1);
    std::println("---");

    // Clean-up.
    delete[] Nv0; delete[] Nv1; delete[] Rv0;
    delete[] Rv1;
    delete[] Rv3;
    delete[] Rv4;

    return 0;
}