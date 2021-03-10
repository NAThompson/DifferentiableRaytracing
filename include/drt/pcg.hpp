#ifndef DRT_PCG_HPP
#define DRT_PCG_HPP
#include <cmath>
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

namespace drt {

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// TODO: Make sure this works with floats!!!
// Also, I suspect that this is not an optimal implementation.
// Something is wrong here.
template<typename Real>
Real pcg_real_01(pcg32_random_t* rng) {
    return ldexp(pcg32_random_r(rng), -32);
}

}
#endif