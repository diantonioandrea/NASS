/**
 * @file Vectors.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Real vectors.
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_VECTORS_HPP
#define NASS_VECTORS_HPP

#include <cstring>

#include "./Reals.hpp"

namespace nass {
    namespace internal {

        // Copy.
        
        void Cp_RvtRvN_0(real_t*, const real_t*, const natural_t&);

        // Operations.

        // Dot product.

        real_t Dt_RvRvN_R(const real_t*, const real_t*, const natural_t&);

    }
}

#endif