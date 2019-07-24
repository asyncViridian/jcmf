#include "chi_square.h"
#include <math.h>

#ifdef USE_R
#include <R.h>
#include <Rembedded.h>
#include <Rinternals.h>
#include <Rmath.h>
#endif


void chi_square::compute(const sequence &seq_a, const sequence &seq_b,
			 double &chisq, double &degrees_of_freedom,
			 double &pvalue)
{
    const int *seq_a_stats = seq_a.getDinuclCount();
    const int *seq_b_stats = seq_b.getDinuclCount();

    chisq = 0.0;
    for(int i=0; i<NUM_DINUCLEOTIDES; ++i)
    {
        chisq += (double) (seq_b_stats[i] - seq_a_stats[i]) * (seq_b_stats[i] - seq_a_stats[i])
                 / (double) seq_a_stats[i];
    }

#ifdef USE_R
    // We actually have 12 degrees of freedom, not 15
    degrees_of_freedom = 12;
    pvalue = pchisq(chisq, degrees_of_freedom, false, false);
#else
    // No pvalue calculation
    degrees_of_freedom = -1.0;
    pvalue = -1.0;
#endif

    return;
}
