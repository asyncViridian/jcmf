#ifndef __RNAZ_RANDOM_SHUFFLE_H__
#define __RNAZ_RANDOM_SHUFFLE_H__

#include<vector>
#include "alignment_block.h"

#define TMP_IN_FILENAME "tmp.maf"
#define TMP_OUT_FILENAME "tmp_rnaz_rand.maf"

class rnaz_random_shuffle
{
 public:
    static void shuffle(const alignment_block &in_block,
			std::vector<alignment_block> &out_blocks);
};


#endif /* __RNAZ_RANDOM_SHUFFLE_H__ */
