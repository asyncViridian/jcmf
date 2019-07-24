#include "assert.h"
#include<string>
#include<map>
#include<vector>
#include "sequence.h"
#include "alignment_block.h"
#include "chi_square.h"
#include "our_shuffle.h"
#include "multiperm.h"
#include "options.h"
#include "math.h"

using namespace std;

extern options g_options;


void multiperm::our_permute(const alignment_block &original_block,
			    const vector<int> &col_ranks,
			    alignment_block &permuted_block,
			    unsigned int &offset,
			    string &violations,
			    string &orig_column_labels,
			    string &perm_column_labels,
			    string &edge_output,
			    bool print_edges)
{
    // Set the type of the multiple alignment
    permuted_block.set_type(original_block.get_type());
    
    // Get the reference sequence
    const string &s = original_block.get_reference().getText();
    
    // Compute our permutation of the reference sequence
    string L;
    vector<int> path;
    our_shuffle::dinuclShuffle(s, original_block, col_ranks, L, path,
			       offset, violations, orig_column_labels,
			       perm_column_labels, edge_output, print_edges);

    // Go through the sequences and rearrange columns according to the
    // permutation just constructed
    vector<sequence>::const_iterator seq_itr =
	original_block.sequences.begin();
    for(; seq_itr < original_block.sequences.end(); ++seq_itr)
    {
	const char *str = seq_itr->getText().c_str();
	int len = seq_itr->getText().length();
	
	// Permute
	string perm;
	for(int i=0; i<len; ++i)
	{
	    perm.append(1,str[path[i]]);
	}
		
	const sequence seq_perm(seq_itr->getSrc(),
				seq_itr->getStart(),
				seq_itr->getSize(), seq_itr->getStrand(),
				seq_itr->getSrcSize(), perm);

	permuted_block.sequences.push_back(seq_perm);
    }
}


void multiperm::print_statistics(FILE *fout,
				 const alignment_block &original_block,
				 const alignment_block &permuted_block,
				 const char *src_label)
{
    vector<sequence>::const_iterator seq_itr =
	original_block.sequences.begin();
    vector<sequence>::const_iterator perm_seq_itr =
	permuted_block.sequences.begin();
    
    seq_itr->print_dinucleotide_header(fout, src_label);
    
    for(; (seq_itr < original_block.sequences.end() &&
	   perm_seq_itr < permuted_block.sequences.end());
	++seq_itr, ++perm_seq_itr)
    {
	seq_itr->print_stats(fout, *perm_seq_itr);
    }
    fprintf(fout, "# \n");
}



void multiperm::our_conserve_permute(const alignment_block &original_block,
				     alignment_block &permuted_block,
				     unsigned int &offset,
				     int &num_conservation_patterns,
				     string &col_ranks_str,
				     string &violations,
				     string &orig_column_labels,
				     string &perm_column_labels,
				     string &edge_output,
				     bool print_edges)
{
    map<string, vector<int> > groups;
    vector<int> col_ranks;
    group_columns(original_block, groups);
    rank_columns(original_block, groups, col_ranks_str, col_ranks);

    num_conservation_patterns = groups.size();
    
    map<string, vector<int> >::const_iterator group_itr = groups.begin();

    our_permute(original_block, col_ranks, permuted_block, offset,
		violations, orig_column_labels, perm_column_labels,
		edge_output, print_edges);
    
    // Compute dinucleotide statistics
    vector<sequence>::iterator seq_itr = permuted_block.sequences.begin();
    for(; seq_itr < permuted_block.sequences.end(); ++seq_itr)
    {
	seq_itr->init_dinucleotide_stats();
    }
    permuted_block.finalize();
}


void multiperm::get_mask(const alignment_block &block,
			 const int col,
			 string &mask)
{
    map<char, int> seen;
    int counter = 0;
    vector<sequence>::const_iterator seq_itr = block.sequences.begin();
    for(; seq_itr < block.sequences.end(); ++seq_itr)
    {
	if (g_options.conservation() == CONSERVATION_GAPS ||
	    g_options.conservation() == CONSERVATION_LEVEL0 ||
	    g_options.conservation() == CONSERVATION_LEVEL1)
	{
	    if (seq_itr->getText()[col] == '-')
		mask.append(1,'-');
	    else
		mask.append(1,'N');
	}
	else if (g_options.conservation() == CONSERVATION_FULL)
	{
	    if (seq_itr->getText()[col] == '-')
		mask.append(1,'-');
	    else if (seen.find(seq_itr->getText()[col]) != seen.end())
	    {
		char buf[2];
		snprintf(buf, sizeof(buf), "%d", seen[seq_itr->getText()[col]]);
		mask.append(buf);
	    }
	    else
	    {
		counter++;
		seen[seq_itr->getText()[col]] = counter;
		char buf[2];
		snprintf(buf, sizeof(buf), "%d", counter);
		mask.append(buf);		
	    }
	}
    }

    if (g_options.conservation() == CONSERVATION_LEVEL0 ||
	g_options.conservation() == CONSERVATION_LEVEL1)
    {
	// Calculates mean pairwise identity for each column
	// (just as in RNAz.pm) and add it to the mask
	int pairs = 0;
	int matches = 0;
	vector<sequence>::const_iterator i_seq_itr = block.sequences.begin();
	for(; i_seq_itr < block.sequences.end(); ++i_seq_itr)
	{
	    vector<sequence>::const_iterator j_seq_itr = i_seq_itr + 1;
	    for(; j_seq_itr < block.sequences.end(); ++j_seq_itr)
	    {
		char ch_i = i_seq_itr->getText()[col];
		char ch_j = j_seq_itr->getText()[col];
		
		if ((ch_i != '-') && (ch_j != '-'))
		{
		    if (ch_i == ch_j)
			matches++;
		    pairs++;
		}
	    }
	}
	
	char buf[16];
	if (pairs>0)
	{
	    
	    if (g_options.conservation() == CONSERVATION_LEVEL0)
	    {
		snprintf(buf, sizeof(buf), "%d",
			 (int) floor( (double) matches / (double) pairs + 0.5));
	    }
	    else if (g_options.conservation() == CONSERVATION_LEVEL1)
	    {
		snprintf(buf, sizeof(buf), "%.1f",
			 (double) matches / (double) pairs);
	    }
	}
	else
	{
	    snprintf(buf, sizeof(buf), "%.0f", 1.0);
	}
	mask.append(buf);
    }
}


void multiperm::group_columns(const alignment_block &block,
			      map<string, vector<int> > &groups)
{
    // Get the first sequence
    vector<sequence>::const_iterator seq_itr = block.sequences.begin();

    for(unsigned int col=0; col<seq_itr->getText().length(); ++col)
    {
	string mask;
	get_mask(block, col, mask);

	groups[mask].push_back(col);
    }
}


void multiperm::rank_columns(const alignment_block &block,
			     const map<string, vector<int> > &groups,
			     string &col_ranks_str,
			     vector<int> &col_ranks)
{
    multimap<int, string> ranking;

    map<string, vector<int> >::const_iterator itr = groups.begin();

    for(; itr != groups.end(); ++itr)
    {
	ranking.insert(pair<int, string>(itr->second.size(), itr->first));
    }

    map<string, char> rankchar;
    map<string, int> ranknum;
    
    const char rankstr[] =
	"0123456789"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"abcdefghijklmnopqrstuvwxyz";

    const int rankstr_len = strlen(rankstr)-1;

    int rank = 1;
    
    multimap<int, string>::reverse_iterator rank_itr = ranking.rbegin();
    for(; rank_itr != ranking.rend(); ++rank_itr, ++rank)
    {
	if (rank > rankstr_len)
	    rankchar[rank_itr->second] = '*';
	else
	    rankchar[rank_itr->second] = rankstr[rank];

	ranknum[rank_itr->second] = rank;
    }

    // Get the first sequence
    vector<sequence>::const_iterator seq_itr = block.sequences.begin();

    for(unsigned int col=0; col<seq_itr->getText().length(); ++col)
    {
	string mask;
	get_mask(block, col, mask);

	col_ranks_str.append(1, rankchar[mask]);
	col_ranks.push_back(ranknum[mask]);
    }    
}


void multiperm::select_columns(const alignment_block &block,
			       const vector<int> &columns,
			       alignment_block &selected_block)
{
    // Prepare selected_block
    vector<sequence>::const_iterator seq_itr = block.sequences.begin();
    int num_cols = seq_itr->getText().length();
    for(; seq_itr < block.sequences.end(); ++seq_itr)
    {
	sequence seq(seq_itr->getSrc(), seq_itr->getStart(),
		     seq_itr->getSize(), seq_itr->getStrand(),
		     seq_itr->getSrcSize(), "");
	selected_block.sequences.push_back(seq);
    }

    // Select columns
    unsigned int col_index = 0;
    for(int i=0; i<num_cols && col_index < columns.size() ; ++i)
    {
	if (i == columns[col_index])
	{
	    seq_itr = block.sequences.begin();
	    vector<sequence>::iterator selected_seq_itr =
		selected_block.sequences.begin();
	    for(; (seq_itr < block.sequences.end() &&
		   selected_seq_itr < selected_block.sequences.end());
		++seq_itr, ++selected_seq_itr)
	    {
		selected_seq_itr->appendToText(seq_itr->getText()[i]);
	    }
	    col_index++;
	}
    }
    selected_block.finalize();
}


void multiperm::replace_columns(alignment_block &block,
				const vector<int> &columns,
				const alignment_block &replacement_block)
{
    if (block.sequences.size() != replacement_block.sequences.size())
    {
	fprintf(stderr, "FATAL ERROR. Attempting to perform in-place replacement of block with %d sequences with block with %d sequences. Exiting.\n", block.sequences.size(), replacement_block.sequences.size());
	exit(-1);
    }

    // Replace columns
    vector<sequence>::iterator seq_itr = block.sequences.begin();
    int num_cols = seq_itr->getText().length();
    unsigned int col_index = 0;
    for(int i=0; i<num_cols && col_index < columns.size() ; ++i)
    {
	if (i == columns[col_index])
	{
	    seq_itr = block.sequences.begin();
	    vector<sequence>::const_iterator replacement_seq_itr =
		replacement_block.sequences.begin();
	    for(; (seq_itr < block.sequences.end() &&
		   replacement_seq_itr < replacement_block.sequences.end());
		++seq_itr, ++replacement_seq_itr)
	    {
		char ch = replacement_seq_itr->getText()[col_index]; 
		seq_itr->replaceInText(i, ch);
	    }
	    col_index++;
	}
    }
}

