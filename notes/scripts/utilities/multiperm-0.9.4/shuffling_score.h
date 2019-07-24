#ifndef __SHUFFLING_SCORE_H__
#define __SHUFFLING_SCORE_H__

#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include "alignment_block.h"

class shuffling_score
{
 public:

    typedef enum {
	ALGO_OUR_SHUFFLE,
	ALGO_RNAZ_SHUFFLE
    } algo;
	    
    
    static void score(const alignment_block &original_block,
		      double &shufflability,
		      double &shuffle_score,
		      double &mean_pvalue,
		      double &sd_pvalue,
		      algo alg = ALGO_OUR_SHUFFLE);
    
    static void rnaz_score(const alignment_block &original_block,
			   double &shuffle_score,
			   double &mean_pvalue,
			   double &sd_pvalue);        
};

#endif /* __SHUFFLING_SCORE_H__ */
