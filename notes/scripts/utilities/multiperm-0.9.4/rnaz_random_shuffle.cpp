#include "rnaz_random_shuffle.h"
#include "maf_file.h"
#include <signal.h>


using namespace std;

void rnaz_random_shuffle::shuffle(const alignment_block &in_block,
				  vector<alignment_block> &out_blocks)
{
#if USE_RNAZ
    
    // Write file that is the input to the RNAz random shuffle
    FILE *fout = fopen(TMP_IN_FILENAME, "w");

    if (!fout)
    {
	fprintf(stderr, "FATAL ERROR. Could not open file '%s' for writing. Exiting.\n", TMP_IN_FILENAME);
    }
    
    fprintf(fout, "a score=0\n");
	
    vector<sequence>::const_iterator seq_itr = in_block.sequences.begin();
    for(; seq_itr < in_block.sequences.end(); ++seq_itr)
    {
	char buf[MAX_LINE_SIZE];
	seq_itr->getMAFline(buf);
	fprintf(fout, "%s", buf);
    }
    fprintf(fout, "\n");
    fclose(fout);

    // Invoke the RNAz Perl script for random shuffling
    // The location of the script, RNAZ_RANDOMIZE_ALN, is defined in the Makefile
    int ret = system(RNAZ_RANDOMIZE_ALN " -l 1 tmp.maf > " TMP_OUT_FILENAME);

    if (WIFSIGNALED(ret) &&
	(WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT))
    {
	raise(WTERMSIG(ret));
	return;
    }
    
    // Read in RNAz permuted random sequence
    maf_file::read(TMP_OUT_FILENAME, out_blocks);

#endif
}
