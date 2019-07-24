#ifndef __OUR_SHUFFLE_H__
#define __OUR_SHUFFLE_H__

#include <string>
#include <vector>
#include <map>
#include "alignment_block.h"
#include "options.h"

class our_shuffle
{
 public:
    struct edge {
	edge(char n, int c) :
	    nucl(n), col(c) {}
	
	char nucl;
	int col;
    };

    typedef std::map<char, std::map<int, std::vector<edge> > > edge_map;
    
    static void dinuclShuffle(const std::string &s,
			      const alignment_block &original_block,
			      const std::vector<int> &col_ranks,
			      std::string &L,
			      std::vector<int> &path,
			      unsigned int &offset,
			      std::string &violations,
			      std::string &orig_column_labels,
			      std::string &perm_column_labels,
			      std::string &edge_output,
			      bool print_edges = false);

 private:
    static void shuffleEdges(edge_map &edges);

    static void pick_edge(const char last_nucl,
			  const int cur_rank,
			  edge_map &edges,
			  char &nucl,
			  int &index,
			  bool &violation);
    
    static void pick_random_edge(const char last_nucl,
				 const int cur_rank,
				 edge_map &edges,
				 char &nucl,
				 int &index);

    static void print_edge_list(edge_map &edges,
				std::string &edge_output);
    
    static void label_columns(const alignment_block &block,
			      const std::vector<int> &path,
			      std::string &orig_column_labels,
			      std::string &perm_column_labels);

    static const char * const base62_str;
};

#endif /* __OUR_SHUFFLE_H__ */
