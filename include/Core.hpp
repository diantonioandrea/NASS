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

#include <cstddef>

// OpenMP support.

#ifdef _OPENMP
#include <omp.h>
#endif

// Neon support and floating-point types.

#if defined(__ARM_NEON) && !defined(NNEON)
#ifndef _NEON
#define _NEON // Custom Neon flag.
#endif
#endif

#if defined(LOOP_OFFSET) || defined(MEMORY_OFFSET) || defined(MEMORY_OFFSET_0) || defined(MEMORY_OFFSET_1)
#error "Unsafe constant definition."
#endif

#ifdef _NEON
#include <arm_neon.h>

#ifdef NEON32 // 32-bit floating-point.

#ifdef NEON64
#error "Unsafe constant definition."
#endif

#define LOOP_OFFSET 8
#define MEMORY_OFFSET 4

namespace nass {
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

#define LOOP_OFFSET 4
#define MEMORY_OFFSET 2

namespace nass {
    using real_t = float64_t;
    
    namespace internal {
        using reals_t = float64x2_t;

        constexpr real_t real_tol = 1.0E-14;
    }
}

#endif

// Multiple offsets for loop-unrolling.

#define MEMORY_OFFSET_0 MEMORY_OFFSET * 0 // Superfluous.
#define MEMORY_OFFSET_1 MEMORY_OFFSET * 1 // Superfluous.

#else // Missing Neon support.

namespace nass {
    using real_t = double;

    namespace internal {
        constexpr real_t real_tol = 1.0E-14;
    }
}

#endif

// Integral types.

namespace nass {
    using natural_t = std::size_t;
    using integer_t = std::ptrdiff_t;
}

// Checks.

#ifndef NVERBOSE
#ifdef _OPENMP
#pragma message "OpenMP enabled."
#else
#pragma message "OpenMP NOT supported/enabled."
#endif

#ifdef _NEON
#pragma message "Neon enabled."
#else
#pragma message "Neon NOT supported/enabled."
#endif
#endif

#endif