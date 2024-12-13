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

// Checks.

#ifndef NDEBUG

#ifdef _OPENMP
#pragma message "OpenMP supported."
#else
#pragma message "OpenMP NOT supported."
#endif

#ifdef __ARM_NEON
#pragma message "Neon supported."
#else
#pragma message "Neon NOT supported."
#endif

#endif

// OpenMP support.

#ifdef _OPENMP
#include <omp.h>
#endif

// Neon support and floating-point types.

#ifdef __ARM_NEON
#ifndef _NEON
#define _NEON // Custom Neon flag.
#endif
#endif

#if defined(LOOP_OFFSET) || defined(MEMORY_OFFSET) || defined(MEMORY_OFFSET_0) || defined(MEMORY_OFFSET_1)
#error "Unsafe constant definition."
#endif

#ifdef _NEON
#include <arm_neon.h>

#ifdef NEON16 // 16-bit floating-point.

#if defined(NEON32) || defined(NEON64)
#error "Unsafe constant definition."
#endif

#define LOOP_OFFSET 16
#define MEMORY_OFFSET 8

namespace nass {
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

#endif