/**
 * @file Decomposition.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief (Thin) QR decomposition.
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

        // Products.

        // [!]

        // (Thin) QR.

        void TQR_RmtRmtNN_0(real_t*, real_t*, const natural_t&, const natural_t&);

    }
}

#endif