#ifndef __ALIGNMENT_BLOCK_H__
#define __ALIGNMENT_BLOCK_H__

#include <map>
#include <vector>
#include "sequence.h"
#include "options.h"

class alignment_block
{
 public:

    typedef enum {
	MULT_ALIGN_MAF,
	MULT_ALIGN_CLUSTALW
    } align_type;

    alignment_block() :
	m_type(MULT_ALIGN_MAF),
	m_total_num_gaps(0),
	m_total_num_nucleotides(0) {}

    alignment_block(align_type type) :
	m_type(type),
	m_total_num_gaps(0),
	m_total_num_nucleotides(0) {}    
    
    std::vector<sequence> sequences;

    void finalize();

    void clear() {
	sequences.clear();
	m_total_num_gaps = 0;
	m_total_num_nucleotides = 0;
	m_reference.clear();
    }

    void print_sequences(FILE *fout, const char first_char = '#',
			 bool verbose = true) const;

    void print_col_ranks(FILE *fout, std::string &col_ranks) const;

    void print_violations(FILE *fout, std::string &violations) const;

    align_type get_type() const {return m_type;}
    
    void set_type(align_type type) {m_type = type;}
    
    int get_num_gaps() const {return m_total_num_gaps;}

    int get_num_nucleotides() const {return m_total_num_nucleotides;}

    const sequence &get_reference() const {return m_reference;}

    bool is_duplicate(alignment_block &block) const;

    bool is_duplicate_column(const unsigned int i, const unsigned int j) const;

    bool is_majority_gaps_column(const unsigned int i) const;

    bool is_all_gaps_column(const unsigned int i) const;

    double score(const alignment_block &block) const;

    bool operator <(const alignment_block &rhs) const;

 private:
    void compute_first_sequence();
    
    void compute_consensus_sequence();

    align_type m_type;
    
    int m_total_num_gaps;
    int m_total_num_nucleotides;

    sequence m_reference;
};


#endif /* __ALIGNMENT_BLOCK_H__ */
