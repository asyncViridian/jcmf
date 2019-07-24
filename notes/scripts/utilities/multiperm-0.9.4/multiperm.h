#ifndef __MULTIPERM_H__
#define __MULTIPERM_H__

#include<string>
#include <cstring>
#include<map>
#include<vector>
#include "alignment_block.h"

class multiperm
{
 public:
    static void our_permute(const alignment_block &original_block,
			    const std::vector<int> &col_ranks,
			    alignment_block &permuted_block,
			    unsigned int &offset,
			    std::string &violations,
			    std::string &orig_column_labels,
			    std::string &perm_column_labels,
			    std::string &edge_output,
			    bool print_edges = false);

    static void print_statistics(FILE *fout,
				 const alignment_block &original_block,
				 const alignment_block &permuted_block,
				 const char *src_label = "");

    static void our_conserve_permute(const alignment_block &original_block,
				     alignment_block &permuted_block,
				     unsigned int &offset,
				     int &num_conservation_patterns,
				     std::string &col_ranks_str,
				     std::string &violations,
				     std::string &edge_output,
				     std::string &orig_column_labels,
				     std::string &perm_column_labels,
				     bool print_edges = false);

 private:
    static void get_mask(const alignment_block &block,
			 const int col,
			 std::string &mask);
    
    static void group_columns(const alignment_block &block,
			      std::map<std::string,std::vector<int> > &groups);

    static void rank_columns(const alignment_block &block,
		       const std::map<std::string, std::vector<int> > &groups,
		       std::string &col_ranks_str,
		       std::vector<int> &col_ranks);

    static void select_columns(const alignment_block &block,
			       const std::vector<int> &columns,
			       alignment_block &selected_block);

    static void replace_columns(alignment_block &block,
				const std::vector<int> &columns,
				const alignment_block &replacement_block);
};


#endif /* __MULTIPERM_H__ */
