#include <vector>
#include "sequence.h"
#include "shuffling_score.h"
#include "multiperm.h"
#include "chi_square.h"
#include <ostream>
#include <math.h>
#include "rnaz_random_shuffle.h"
#include "assert.h"


using namespace std;

void shuffling_score::score(const alignment_block &original_block,
			    double &shufflability,
			    double &shuffle_score,
			    double &mean_pvalue,
			    double &sd_pvalue,
			    algo alg)
{
    const int num_tries = (alg == ALGO_RNAZ_SHUFFLE) ? 25 : 1000;
    int num_duplicates = 0;
#ifdef USE_R
    double sum_pvalues = 0.0;
    double sum_sq_pvalues = 0.0;
    int num_pvalues = 0;
#endif

    shufflability = 0.0;
    vector<alignment_block> blocks;

    for(int i=0; i<num_tries; ++i)
    {
	alignment_block perm_block;

	if (alg == ALGO_RNAZ_SHUFFLE)
	{
	    vector<alignment_block> rnaz_blocks;
	    rnaz_random_shuffle::shuffle(original_block, rnaz_blocks);
	    perm_block = rnaz_blocks[0];	    
	}
	else // (alg == ALGO_OUR_SHUFFLE)
	{
	    unsigned int offset;
	    int num_conservation_patterns;
	    string col_ranks;
	    string violations;
	    string orig_column_labels;
	    string perm_column_labels;
	    string edge_output;
	    multiperm::our_conserve_permute(original_block, perm_block,
					    offset,
					    num_conservation_patterns,
					    col_ranks, violations,
					    orig_column_labels,
					    perm_column_labels,
					    edge_output, false);	
	}
	
		
	if (original_block.is_duplicate(perm_block))
	    num_duplicates++;

	shufflability += original_block.score(perm_block);

#ifdef USE_R
	vector<sequence>::const_iterator seq_itr =
	    original_block.sequences.begin();
	vector<sequence>::const_iterator perm_seq_itr =
	    perm_block.sequences.begin();
	
	for(; (seq_itr < original_block.sequences.end() &&
	       perm_seq_itr < perm_block.sequences.end());
	    ++seq_itr, ++perm_seq_itr)
	{
	    double chisq = 0.0;
	    double df = 0;
	    double pvalue = 0.0;
	    
	    chi_square::compute(*seq_itr, *perm_seq_itr, chisq, df, pvalue);
	    
	    sum_pvalues += pvalue;
	    sum_sq_pvalues += pvalue * pvalue;
	    num_pvalues++;
	}
#endif

	blocks.push_back(perm_block);
    }

    sort(blocks.begin(), blocks.end());

    assert(blocks.size() > 0);
    int num_unique=1;

    vector<alignment_block>::iterator prev = blocks.begin();
    vector<alignment_block>::iterator cur = prev;
    for(++cur; cur < blocks.end(); ++cur, ++prev)
    {
	if (prev->is_duplicate(*cur))
	    continue;

	num_unique++;
    }

#ifdef USE_R
    mean_pvalue = sum_pvalues / (double) num_pvalues;
    sd_pvalue = sqrt((sum_sq_pvalues - num_pvalues * mean_pvalue * mean_pvalue)
		     /(num_pvalues - 1));
#else
    mean_pvalue = -1.0;
    sd_pvalue = -1.0;
#endif
    
    shuffle_score = (double) num_unique / (double) num_tries;
    shufflability /= (double) num_tries;
}


void shuffling_score::rnaz_score(const alignment_block &original_block,
				 double &shuffle_score,
				 double &mean_pvalue,
				 double &sd_pvalue)
{
    const int num_tries = 25;
    int num_duplicates = 0;
    double sum_pvalues = 0.0;
    double sum_sq_pvalues = 0.0;
    int num_pvalues = 0;

    for(int i=0; i<num_tries; ++i)
    {
	vector<alignment_block> rnaz_blocks;
	rnaz_random_shuffle::shuffle(original_block, rnaz_blocks);
	alignment_block &perm_block = rnaz_blocks[0];
	
	if (original_block.is_duplicate(perm_block))
	    num_duplicates++;
	
	vector<sequence>::const_iterator seq_itr =
	    original_block.sequences.begin();
	vector<sequence>::const_iterator perm_seq_itr =
	    perm_block.sequences.begin();
	
	for(; (seq_itr < original_block.sequences.end() &&
	       perm_seq_itr < perm_block.sequences.end());
	    ++seq_itr, ++perm_seq_itr)
	{
	    double chisq = 0.0;
	    double df = 0;
	    double pvalue = 0.0;
	    
	    chi_square::compute(*seq_itr, *perm_seq_itr, chisq, df, pvalue);
	    
	    sum_pvalues += pvalue;
	    sum_sq_pvalues += pvalue * pvalue;
	    num_pvalues++;
	}
    }

    //printf("sum_pvalues = %f, sum_sq_pvalues = %f, num_pvalues = %d\n", sum_pvalues, sum_sq_pvalues, num_pvalues);
    
    shuffle_score = 1.0 - (double) num_duplicates / (double) num_tries;
    mean_pvalue = sum_pvalues / (double) num_pvalues;
    sd_pvalue = sqrt((sum_sq_pvalues - num_pvalues * mean_pvalue * mean_pvalue)
		     /(num_pvalues - 1));
}

