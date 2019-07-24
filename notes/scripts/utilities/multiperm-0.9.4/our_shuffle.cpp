#include "our_shuffle.h"
#include "assert.h"
#include "prng.h"
#include <math.h>
#include <string.h>
#ifdef DEBUG
#include <iostream>
#endif

using namespace std;

extern options g_options;

const char * const our_shuffle::base62_str =
"0123456789"
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz";


void our_shuffle::dinuclShuffle(const std::string &s,
				const alignment_block &original_block,
				const std::vector<int> &col_ranks,
				std::string &L,
				std::vector<int> &path,
				unsigned int &offset,
				string &violations,
				string &orig_column_labels,
				string &perm_column_labels,
				string &edge_output,
				bool print_edges)
{
    edge_map edges;

#ifdef DEBUG
    cout << s << ' ' << s.size() << ' ' << col_ranks.size() << "\n";
    for(int ii =0; ii< 100; ii++){cout << col_ranks[ii] << ',';};
#endif

    assert(s.size() == col_ranks.size());

    unsigned int n = s.size();
    if (g_options.norandomoffset())
	offset = 0;
    else
	offset = (unsigned int) floor(get_rand() * n);
    assert(offset < n);

    int col = (0 + offset) % n;
    char last_nucl = s[col];
    
    char last_non_maj_gaps_nucl = '\0';
    bool last_maj_gaps = original_block.is_majority_gaps_column(col);
    bool majhdl = !(g_options.nomajgaphandling());

    // Construct the graph
    for(unsigned int i=1; i<n; ++i)
    {
	int col = (i + offset) % n;
	
        int cur_rank = col_ranks[col];
	bool cur_maj_gaps = original_block.is_majority_gaps_column(col);

	// If there was a transition from majority gaps to not...
	if (majhdl && last_non_maj_gaps_nucl && last_maj_gaps && !cur_maj_gaps)
	    edges[last_non_maj_gaps_nucl][cur_rank].push_back(edge(s[col],col));
	else // no transition
	    edges[last_nucl][cur_rank].push_back(edge(s[col],col));

	last_nucl = s[col];
	last_maj_gaps = cur_maj_gaps;
	if (!cur_maj_gaps)
	    last_non_maj_gaps_nucl = s[col];
    }

    if (print_edges)
    {
	edge_output.append("# Original edge list:\n");
	print_edge_list(edges, edge_output);
    }

    // Shuffle the edges
    shuffleEdges(edges);

    if (print_edges)
    {
	edge_output.append("# \n# Shuffled edge list:\n");
	print_edge_list(edges, edge_output);
    }

    L.resize(n, ' ');
    path.resize(n, 0);
    violations.resize(n, ' ');
    
    col = (0 + offset) % n;
    last_nucl = s[col];
    L[col] = last_nucl;
    path[col] = col;
    violations[col] = '^';

    last_non_maj_gaps_nucl = '\0';
    last_maj_gaps = original_block.is_majority_gaps_column(0);

    for(unsigned int i=1; i<n; ++i)
    {
	int col = (i + offset) % n;
	
	int cur_rank = col_ranks[col];
	bool cur_maj_gaps = original_block.is_majority_gaps_column(col);

	char nucl;
	int index;
	bool violation = false;

	// If there was a transition from majority gaps to not...
	if (majhdl && last_non_maj_gaps_nucl && last_maj_gaps && !cur_maj_gaps)
	    pick_edge(last_non_maj_gaps_nucl, cur_rank, edges, nucl, index, violation);
	else // no transition
	    pick_edge(last_nucl, cur_rank, edges, nucl, index, violation);
	
	L[col] = nucl;
	path[col] = index;
	violations[col] = violation ? '*' : ' '; 

	last_nucl = nucl;
	last_maj_gaps = cur_maj_gaps;
	if (!cur_maj_gaps)
	    last_non_maj_gaps_nucl = nucl;
    }

    if (print_edges)
	label_columns(original_block, path, orig_column_labels,
		      perm_column_labels);
}



void our_shuffle::shuffleEdges(edge_map &edges)
{
    edge_map::iterator edge_itr = edges.begin();

    for(; edge_itr != edges.end(); ++edge_itr)
    {
	map<int, vector<edge> >::iterator edge_rank_itr = 
	    edge_itr->second.begin();

	for(; edge_rank_itr != edge_itr->second.end(); ++edge_rank_itr)
	{
	    vector<edge> &edge_vec = edge_rank_itr->second;

	    int n = edge_vec.size();

	    if (n <= 1)
		continue; // nothing to do

	    for(int i=n; --i; )
	    {
		int j = (int) (get_rand() * (i + 1));
		assert(j <= i);

		char tmp_nucl = edge_vec[j].nucl;
		int tmp_col = edge_vec[j].col;

		edge_vec[j].nucl = edge_vec[i].nucl;
		edge_vec[j].col = edge_vec[i].col;

		edge_vec[i].nucl = tmp_nucl;
		edge_vec[i].col = tmp_col;
	    }   
	}
    }
}


void our_shuffle::pick_edge(const char last_nucl,
			    const int cur_rank,
			    edge_map &edges,
			    char &nucl,
			    int &index,
			    bool &violation)
{
    vector<edge> &edge_vec = edges[last_nucl][cur_rank];

    if (edge_vec.empty())
    {
	pick_random_edge(last_nucl, cur_rank, edges, nucl, index);
	violation = true;
	return;
    }

    vector<edge>::iterator edge_itr = edge_vec.begin();

    nucl = edge_itr->nucl;
    index = edge_itr->col;

    edge_vec.erase(edge_itr);

    violation = false;
}


void our_shuffle::pick_random_edge(
    const char last_nucl,
    const int cur_rank,
    edge_map &edges,
    char &nucl,
    int &index)
{
    edge_map::iterator edge_itr = edges.begin();
    int total = 0;
    
    for(; edge_itr != edges.end(); ++edge_itr)
	total += edge_itr->second[cur_rank].size();

    double rnd_num = get_rand();
    
    int running_total = 0;
    edge_itr = edges.begin();
    for(; edge_itr != edges.end(); ++edge_itr)
    {
	if (edge_itr->second[cur_rank].empty())
	    continue;

	running_total += edge_itr->second[cur_rank].size();

	if (rnd_num < (double) running_total / (double) total)
	{
	    vector<edge> &edge_vec = edge_itr->second[cur_rank];

	    vector<edge>::iterator itr = edge_vec.begin();
	    
	    nucl = itr->nucl;
	    index = itr->col;

	    edge_vec.erase(itr);

	    return;
	}   
    }

    // This should never happen!
    assert(0);
}


void our_shuffle::print_edge_list(edge_map &edges,
				  string &edge_output)
{
    edge_map::iterator edge_itr = edges.begin();
    for(; edge_itr != edges.end(); ++edge_itr)
    {
	string top;
	string bottom;
	top.append("# ");
	top.append(1, edge_itr->first);
	top.append(": ");
	bottom.append("#    ");

	map<int, vector<edge> >::iterator edge_rank_itr = 
	    edge_itr->second.begin();
	for(; edge_rank_itr != edge_itr->second.end(); ++edge_rank_itr)
	{
	    unsigned int rank = edge_rank_itr->first;
	    char rank_ch = rank >= strlen(base62_str) ? '*':base62_str[rank];
	    
	    vector<edge> &edge_vec = edge_rank_itr->second;

	    for(unsigned int i=0; i<edge_vec.size(); ++i)
	    {
		top.append(1, edge_vec[i].nucl);
		bottom.append(1, rank_ch);
	    }
	    top.append(1, ' ');
	    bottom.append(1, ' ');
	}
	edge_output.append(top);
	edge_output.append(1, '\n');
	edge_output.append(bottom);
	edge_output.append(1, '\n');
    }
}


void our_shuffle::label_columns(const alignment_block &block,
				const vector<int> &path,
				string &orig_column_labels,
				string &perm_column_labels)
{
    unsigned int n = block.sequences.front().getText().length();
    vector<int> column_labels(n,0);

    int label = 1;
    for(unsigned int i=0; i<n; ++i)
    {
	if (column_labels[i] != 0)
	    continue;

	column_labels[i] = label;

	for(unsigned int j=i+1; j<n; ++j)
	{
	    if (block.is_duplicate_column(i,j))
		column_labels[j] = label;
	}

	label++;
    }

    string orig_labels_top;
    string orig_labels_bottom;

    int base = strlen(base62_str);
    
    for(unsigned int i=0; i<column_labels.size(); ++i)
    {
	int label = column_labels[i];

	if (label >= base*base)
	{
	    orig_labels_top.append(1,'*');
	    orig_labels_bottom.append(1, '*');
	}
	
	orig_labels_top.append(1, base62_str[label / base]);
	orig_labels_bottom.append(1, base62_str[label % base]);
    }
    
    string perm_labels_top;
    string perm_labels_bottom;

    for (unsigned int i=0; i<path.size(); ++i)
    {
	int label = column_labels[path[i]];

	if (label >= base*base)
	{
	    perm_labels_top.append(1,'*');
	    perm_labels_bottom.append(1, '*');
	}
	
	perm_labels_top.append(1, base62_str[label / base]);
	perm_labels_bottom.append(1, base62_str[label % base]);
    }
    
    orig_column_labels =
	"#                 COLUMN...                                " +
	orig_labels_top + "\n" +
	"#                 ...NUMBER                                " +
	orig_labels_bottom + "\n";
    
    perm_column_labels =
	"#                 COLUMN...                                " +
	perm_labels_top + "\n"
	"#                 ...NUMBER                                " +
	perm_labels_bottom + "\n";
}


