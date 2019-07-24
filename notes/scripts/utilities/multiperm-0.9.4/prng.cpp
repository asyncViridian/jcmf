/* prng.cpp
**
*/

#include <stdlib.h>
#include <stdio.h>
#include "mt19937ar.h"
#include "prng.h"

#ifdef USE_R
#include "R.h"
#endif

void init_rand()
{
#ifdef USE_R
    // Get R's Random Number Generator ready
    GetRNGstate();
#else
    FILE *dev_rand = fopen(DEV_URANDOM, "r");
    if(!dev_rand)
    {
	fprintf(stderr, "FATAL ERROR. Couldn't open '%s' for seeding the "
		"Pseudo-Random Number Generator.\n", DEV_URANDOM);
	exit(-1);
    }

    unsigned long seed = 0;

    if (fread(&seed, sizeof(seed), 1, dev_rand) != 1)
    {
	fprintf(stderr, "FATAL ERROR. Couldn't read enough bytes from '%s' "
		"for seeding the Pseudo-Random Number Generator.\n",
		DEV_URANDOM);
	exit(-1);
	
    }

    init_genrand(seed);
#endif
}


double get_rand()
{
#ifdef USE_R
    return unif_rand();
#else
    return genrand_real2();    
#endif
}
