/**
 * @file Sparse.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Real sparse matrices.
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_SPARSE_HPP
#define NASS_SPARSE_HPP

#include <string>
#include <tuple>

#include "./Vectors.hpp"

namespace nass {
    namespace internal {

        // Load. [!] These two methods can be unified.

        [[nodiscard]] std::tuple<natural_t, natural_t*, natural_t*, real_t*> Spc_St_NNvNvRv(const std::string&);
        [[nodiscard]] std::tuple<natural_t, natural_t*, natural_t*, real_t*> Spr_St_NNvNvRv(const std::string&);

        // Products.

        void Mlc_RvtNNvNvRvRv_0(real_t*, const natural_t&, const natural_t*, const natural_t*, const real_t*, const real_t*);
        void Mlr_RvtNNvNvRvRv_0(real_t*, const natural_t&, const natural_t*, const natural_t*, const real_t*, const real_t*);

    }
}

#endif