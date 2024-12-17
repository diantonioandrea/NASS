/**
 * @file Decomposition.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief (Thin) QR decomposition.
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_DECOMPOSITION_HPP
#define NASS_DECOMPOSITION_HPP

#include "./Core.hpp"

namespace nass {
    namespace internal {

        // Pivoting.

        void Pv_RmtNvtNNN_0(real_t*, natural_t*, const natural_t&, const natural_t&, const natural_t&);

        // Products.

        void Lh_RmtRvNNN_0(real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);

        void Mq_RvtRmNN_0(real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);
        void Mqt_RvtRmNN_0(real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);

        // (Thin) QR.

        void TQR_RmtRmtNvtNN_0(real_t*, real_t*, natural_t*, const natural_t&, const natural_t&);

    }
}

#endif