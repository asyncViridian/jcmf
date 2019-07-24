/* prng.h
**
*/

#ifndef __PRNG_H__
#define __PRNG_H__

#define MT19937AR_N 624

#define DEV_URANDOM "/dev/urandom"

// Seed the pseudo-random number generator
void init_rand();

// Get a uniform random double from [0,1) 
double get_rand();

#endif /* __PRNG_H__ */
