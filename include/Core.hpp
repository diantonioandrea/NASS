/**
 * @file Core.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Core declarations.
 * @date 2024-12-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_CORE_HPP
#define NASS_CORE_HPP

// Includes.

#include <cstddef>


// Neon support and floating-point types.

#if defined(LOOP_OFFSET) || defined(MEMORY_OFFSET) || defined(MEMORY_OFFSET_0) || defined(MEMORY_OFFSET_1) || defined(MEMORY_OFFSET_2) || defined(MEMORY_OFFSET_3)
#error "Unsafe constant definition."
#endif

#ifdef __ARM_NEON
#include <arm_neon.h>

#ifdef NEON16 // 16-bit floating-point.

#if defined(NEON32) || defined(NEON64)
#error "Unsafe constant definition."
#endif

#define LOOP_OFFSET 32
#define MEMORY_OFFSET 8

namespace nla {
    using real_t = float16_t;
    
    namespace internal {
        using reals_t = float16x8_t;

        constexpr real_t real_tol = 5.0E-3;
    }
}

#else
#ifdef NEON32 // 32-bit floating-point.

#ifdef NEON64
#error "Unsafe constant definition."
#endif

#define LOOP_OFFSET 16
#define MEMORY_OFFSET 4

namespace nla {
    using real_t = float32_t;
    
    namespace internal {
        using reals_t = float32x4_t;

        constexpr real_t real_tol = 1.0E-7;
    }
}

#else // 64-bit floating-point.

#ifndef NEON64
#define NEON64
#endif

#define LOOP_OFFSET 8
#define MEMORY_OFFSET 2

namespace nla {
    using real_t = float64_t;
    
    namespace internal {
        using reals_t = float64x2_t;

        constexpr real_t real_tol = 1.0E-14;
    }
}

#endif
#endif

// Multiple offsets for loop-unrolling.

#define MEMORY_OFFSET_0 MEMORY_OFFSET * 0 // Superfluous.
#define MEMORY_OFFSET_1 MEMORY_OFFSET * 1 // Superfluous.
#define MEMORY_OFFSET_2 MEMORY_OFFSET * 2
#define MEMORY_OFFSET_3 MEMORY_OFFSET * 3

#else // Missing Neon support.

namespace nla {
    using real_t = double;

    namespace internal {
        constexpr real_t real_tol = 1.0E-14;
    }
}

#endif


// Integral types.

namespace nla {
    using natural_t = std::size_t;
    using integer_t = std::ptrdiff_t;
}

#endif