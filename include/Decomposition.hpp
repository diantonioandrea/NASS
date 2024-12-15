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

        // Pivoting.

        void Pv_RmtNvtNNN_0(real_t*, natural_t*, const natural_t&, const natural_t&, const natural_t&);

        // Products.

        void Lh_RmtRvNNN_0(real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);
        void Rh_RmtRvNNN_0(real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);

        void Lh_RvtRvNN_0(real_t*, const real_t*, const natural_t&, const natural_t&);

        // (Thin) QR.

        void TQR_RmtRmtNvtRvtNN_0(real_t*, real_t*, natural_t*, real_t*, const natural_t&, const natural_t&);

    }
}

#endif