/**
 * @file Solver.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Sketched GMRES.
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_SOLVER_HPP
#define NASS_SOLVER_HPP

#include "./Core.hpp"

namespace nass {
    namespace internal {

        // sGMRES.

        void sGMRES_RvNNvNvRvRvNN_0(real_t*, const natural_t&, const natural_t*, const natural_t*, const real_t*, const real_t*, const natural_t&, const natural_t&);

    }
}

#endif