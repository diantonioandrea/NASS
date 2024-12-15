/**
 * @file Matrix.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Real matrices.
 * @date 2024-12-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_MATRICES_HPP
#define NASS_MATRICES_HPP

#include "./Core.hpp"

namespace nass {
    namespace internal {

        // Products.
        
        void Ml_RvtRmRvNN_0(real_t*, const real_t*, const real_t*, const natural_t&, const natural_t&);

        void Ml_RmtTRmRmNNN_0(real_t*, const real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);
        void Ml_RmtRmRmNNN_0(real_t*, const real_t*, const real_t*, const natural_t&, const natural_t&, const natural_t&);

        // Output.

        void Pr_RmNN_0(const real_t*, const natural_t&, const natural_t&);

    }
}

#endif