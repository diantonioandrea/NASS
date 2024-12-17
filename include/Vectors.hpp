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

#include "./Reals.hpp"

namespace nass {
    namespace internal {

        // Copy.
        
        void Cp_RvtRvN_0(real_t*, const real_t*, const natural_t&);

        // Operations.

        // Dot product.

        real_t Dt_RvRvN_R(const real_t*, const real_t*, const natural_t&);
        real_t NPDt_RvRvN_R(const real_t*, const real_t*, const natural_t&);

        // Norm.

        real_t Nr_RvN_R(const real_t*, const natural_t&);
        void Nrz_RvtN_0(real_t*, const natural_t&);

        real_t NPNr_RvN_R(const real_t*, const natural_t&);
        void NPNrz_RvtN_0(real_t*, const natural_t&);

        // Projections.
        
        void Prj_RvtRvRvN_0(real_t*, const real_t*, const real_t*, const natural_t&);

        // Output.

        void Pr_RrvN_0(const real_t*, const natural_t&);
        void Pr_RcvN_0(const real_t*, const natural_t&);

    }
}

#endif