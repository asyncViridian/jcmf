#include <string>
#include <map>
#include <vector>
#include "assert.h"
#include "sequence.h"
#include "alignment_block.h"
#include "chi_square.h"
#include "base_file.h"
#include "maf_file.h"
#include "clustalw_file.h"
#include "multiperm.h"
#include "rnaz_random_shuffle.h"
#include "options.h"
#include "shuffling_score.h"
#include "prng.h"
#include <signal.h>
#include <stdarg.h>

#ifdef USE_R
#include <R.h>
#include <Rembedded.h>
#include <Rinternals.h>
#endif



using namespace std;


// Global options
options g_options;

inline void verbose(FILE *stream, const char *format, ...)
{
    if (g_options.verbose())
    {
	va_list ap;
	va_start(ap, format);
	vfprintf(stream, format, ap);
	va_end(ap);
    }
}


static void sigint_handler(int sig)
{
    fprintf(stderr, "\nFATAL ERROR. Caught termination signal. Exiting!\n");

#ifdef USE_R
    // R related cleanup
    // PutRNGstate();
    // Rf_endEmbeddedR(1);
#endif    

    exit(-1);
}


#ifdef USE_R
void my_R_WriteConsole(char *buf, int len)
{
    printf("R message: %s", buf);
}
#endif
     

int main(int argc, char *argv[])
{
    // Parse command line options
    g_options.parse(argc, argv);

#ifdef USE_R
    // Initialize R (RHOME is defined in the Makefile)
    if (setenv("R_HOME", RHOME, 1) !=0 )
    {
	fprintf(stderr, "FATAL ERROR. Couldn't set R_HOME in environment\n");
	exit(-1);
    }
    char *R_argv[] = {"Rembedded", "--gui=none", "--silent",
		      "--slave", "--vanilla"};
    int R_argc = sizeof(R_argv) / sizeof(R_argv[0]);
    Rf_initEmbeddedR(R_argc, R_argv);
#endif
    
    // Seed the PRNG
    init_rand();
    
    // Set up signal handler, for cleanup
    signal(SIGINT, sigint_handler);
    
    // Open the multiple sequence alignment file
    base_file *fin;
    if (g_options.clustalw())
    {
	// Open the Clustal W file
	fin = new clustalw_file(g_options.maf_filename().c_str());
    }
    else
    {
	// Open the MAF file
	fin = new maf_file(g_options.maf_filename().c_str());
    }

    // Open the output files
    vector<FILE *> fout;
    for(int i=0; i<g_options.num(); ++i)
    {
	char filename[256];

	char *maf_file = (char *) g_options.maf_filename().c_str();
	char *last_slash = strrchr(maf_file, '/');
	if (last_slash) maf_file = last_slash + 1;
	
	snprintf(filename, sizeof(filename), "perm_%03d_%s", (i+1), maf_file);

	FILE *out = fopen(filename, "w");

	if (!out)
	{
	    fprintf(stderr, "FATAL ERROR. Unable to open output file '%s'. "
		    "Exiting.\n", filename);
	    exit(-1);
	}

	fout.push_back(out);
    }
    
    // Run through all the multiple alignments
    alignment_block block;
    string comments;
    bool verb = g_options.verbose();
    while(fin->get_next_block(block, comments))
    {
	if (!g_options.discardgapseqs() ||
	    (g_options.discardgapseqs() &&
	     block.get_num_gaps() == 0))
	{
	    double our_shufflability = 0.0;
	    double our_score = 0.0;
	    double our_mean_pvalue = 0.0;
	    double our_sd_pvalue = 0.0;
	    double rnaz_shufflability = 0.0;
	    double rnaz_score = 0.0;
	    double rnaz_mean_pvalue = 0.0;
	    double rnaz_sd_pvalue = 0.0;
	    if (g_options.experimentalscore())
	    {
		// Compute the shuffling score for this block
		// The score ranges from 0.0 to 1.0, with 1.0
		// being the best (lots of random permutations possible)
		shuffling_score::score(block, our_shufflability, our_score,
				       our_mean_pvalue, our_sd_pvalue,
				       shuffling_score::ALGO_OUR_SHUFFLE);
	    }
	    if (g_options.runrnaz())
	    {		
		shuffling_score::score(block, rnaz_shufflability, rnaz_score,
				       rnaz_mean_pvalue, rnaz_sd_pvalue,
				       shuffling_score::ALGO_RNAZ_SHUFFLE);
	    }

	    for(int i=0; i<g_options.num(); ++i)
	    {
		fprintf(fout[i], "%s", comments.c_str());
		
		// Write out original alignment
		verbose(fout[i], "# Original alignment block:\n");
		if (verb) block.print_sequences(fout[i]);

		// Compute permutations
		alignment_block perm_block;
		unsigned int offset;
		int num_conservation_patterns;
		string col_ranks;
		string violations;
		string orig_column_labels;
		string perm_column_labels;
		string edge_output;
		multiperm::our_conserve_permute(block, perm_block, offset,
						num_conservation_patterns,
						col_ranks, violations,
						orig_column_labels,
						perm_column_labels,
						edge_output, true);	
		
		if (verb) block.print_col_ranks(fout[i], col_ranks);
		verbose(fout[i], "%s", orig_column_labels.c_str()); 
		verbose(fout[i], "# \n");	
		verbose(fout[i], "%s# \n", edge_output.c_str());

		if (g_options.experimentalscore())
		{
		    if (our_mean_pvalue >= 0.0)
		    {
			verbose(fout[i],
				"# With 1000 multiperm permutations, "
				"Mean pvalue = %.3f "
				"(standard deviation = %.3f), "
				"Shuffling Score = %.3f, "
				"Shufflability = %.3f\n",
				our_mean_pvalue, our_sd_pvalue, our_score, 
				our_shufflability);
		    }
		    else
		    {
			verbose(fout[i],
				"# With 1000 multiperm permutations, "
				"Shuffling Score = %.3f, "
				"Shufflability = %.3f\n",
				our_score, 
				our_shufflability);			
		    }
		}

		if (g_options.runrnaz())
		{
		    if (rnaz_mean_pvalue >= 0.0)
		    {
			verbose(fout[i],
				"# With 25 RNAz permutations, "
				"Mean pvalue = %.3f "
				"(standard deviation = %.3f), "
				"Shuffling Score = %.3f, "
				"Shufflability = %.3f\n",
				rnaz_mean_pvalue, rnaz_sd_pvalue, rnaz_score,
				rnaz_shufflability);
		    }
		    else
		    {
			verbose(fout[i],
				"# With 25 RNAz permutations, "
				"Shuffling Score = %.3f, "
				"Shufflability = %.3f\n",
				rnaz_score,
				rnaz_shufflability);			
		    }
		}

		verbose(fout[i], "# Random column offset: %d\n", offset);
		verbose(fout[i], "# Reference sequence length: %d\n",
			block.get_reference().getText().size());
		verbose(fout[i], "# Number of conservation patterns: %d\n# \n",
			num_conservation_patterns);

		verbose(fout[i], "# Random alignment block (using our shuffle):\n");
		perm_block.print_sequences(fout[i], 's', verb);
		if (verb) perm_block.print_col_ranks(fout[i], col_ranks);
		verbose(fout[i], "%s", perm_column_labels.c_str());

		if (verb) block.print_violations(fout[i], violations);

		verbose(fout[i], "# \n");

		if (verb)
		    multiperm::print_statistics(fout[i], block, perm_block,
						"Original / Our Permuted");
		
		// Run the RNAz random shuffling algorithm, if necessary
		if (g_options.runrnaz())
		{
		    vector<alignment_block> rnaz_blocks;
		    // Shuffle using RNAz algorithm
		    rnaz_random_shuffle::shuffle(block, rnaz_blocks);
		    verbose(fout[i], "# Random alignment block (using the RNAz algorithm):\n");
		    if (verb) rnaz_blocks.front().print_sequences(fout[i]);
		    verbose(fout[i], "# \n");
		    if (verb)
			multiperm::print_statistics(fout[i], block,rnaz_blocks.front(),
						    "Original / RNAz Permuted");
		}
		
		verbose(fout[i], "# ---\n");
		fprintf(fout[i], "\n");
	    }
	}
    }

    // Close the output files
    for(int i=0; i<g_options.num(); ++i)
    {
	fclose(fout[i]);
    }
    delete fin;

#ifdef USE_R
    // R related cleanup
    PutRNGstate();
    Rf_endEmbeddedR(1);
#endif

    return 0;
}
