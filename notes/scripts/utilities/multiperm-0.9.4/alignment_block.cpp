#include "alignment_block.h"
#include "config.h"

using namespace std;

extern options g_options;

bool alignment_block::operator <(const alignment_block &rhs) const
{
    vector<sequence>::const_iterator seq_itr = sequences.begin();
    vector<sequence>::const_iterator rhs_seq_itr = rhs.sequences.begin();
    
    for(; seq_itr < sequences.end() && rhs_seq_itr < rhs.sequences.end();
	++seq_itr, ++rhs_seq_itr)
    {
	const string &lhs_text = seq_itr->getText();
	const string &rhs_text = rhs_seq_itr->getText();

	if (lhs_text < rhs_text)
	    return true;
	if (lhs_text > rhs_text)
	    return false;
    }

    // The two blocks are equal
    return false;
}


void alignment_block::finalize()
{
    // Called from clustalw_file:: or maf_file::get_next_block
    // after a complete block has been read.
    // Compute block statistics & reference sequence.

    ////////////
    //
    // 14 Sep 2010: quick fix - columns containing *only* dashes
    // cause a crash, since the reference seq excludes them, but col_ranks doesn't;
    // so we just delete any such columns (maybe there's a better choice). 
    // no idea why hg18 17-way alignments should have them.  --- wlr
    
    int len = sequences[0].getText().length();
    for(int i=0; i < len; ++i)
    {
      while(is_all_gaps_column(i))
	{
	  vector<sequence>::iterator seq_itr = sequences.begin();	
	  for(; seq_itr < sequences.end(); ++seq_itr)
	    {
	      seq_itr->eraseInText(i,1);
	    }
	}
    }

    //
    ////////////

    vector<sequence>::const_iterator seq_itr = sequences.begin();	
    for(; seq_itr < sequences.end(); ++seq_itr)
    {
	const char *str = seq_itr->getText().c_str();
	int len = seq_itr->getText().length();
	
	for(int i=0; i < len; ++i)
	{
	    if (str[i] == '-')
		m_total_num_gaps++;
	    
	    m_total_num_nucleotides++;
	}
    }

    // Determine the reference sequence
    if (g_options.refseq() == REFSEQ_FIRST)
	compute_first_sequence();
    else // (g_options.refseq() == REFSEQ_CONSENSUS)
	compute_consensus_sequence();
}

bool alignment_block::is_all_gaps_column(const unsigned int i) const
{
    vector<sequence>::const_iterator itr = sequences.begin();
    for(; itr < sequences.end(); ++itr)
    {
	if (itr->getText()[i] != '-')
	  return false;
    }
    return true;
}

void alignment_block::print_sequences(FILE *fout, const char first_char, bool verbose) const
{
    if (m_type == MULT_ALIGN_MAF)
    {
	vector<sequence>::const_iterator seq_itr = sequences.begin();
	for(; seq_itr < sequences.end(); ++seq_itr)
	{
	    fprintf(fout, "%c %-25s %10d %6d %c %10d %s\n", first_char,
		    seq_itr->getSrc().c_str(),
		    seq_itr->getStart(),
		    seq_itr->getSize(),
		    seq_itr->getStrand(),
		    seq_itr->getSrcSize(),
		    seq_itr->getText().c_str());
	}
	if (verbose)
	    fprintf(fout, "# %25s %30s %s\n", m_reference.getSrc().c_str(), "",
		    m_reference.getText().c_str());
    }
    else if (m_type == MULT_ALIGN_CLUSTALW)
    {
	vector<sequence>::const_iterator seq_itr = sequences.begin();
	unsigned int len = seq_itr->getText().length();
	unsigned int num_chunks = len / 60;
	if (len % 60)
	    num_chunks++;

	for(unsigned int i=0; i<num_chunks; ++i)
	{
	    unsigned int start = i*60;
	    unsigned int end = (i+1)*60;
	    if (end > len)
		end = len;

	    for(seq_itr=sequences.begin(); seq_itr<sequences.end(); ++seq_itr)
	    {
		fprintf(fout, "%-35s %s\n",
			seq_itr->getSrc().c_str(),
			seq_itr->getText().substr(start, end-start).c_str());
	    }
	    fprintf(fout, "\n");
	}
    }
}


void alignment_block::print_col_ranks(FILE *fout, string &col_ranks) const
{
    fprintf(fout, "# %25s %30s %s\n", "CONSERVATION PATTERNS", "", col_ranks.c_str());
}


void alignment_block::print_violations(FILE *fout, string &violations) const
{
    fprintf(fout, "# %25s %30s %s\n", "DINUCL CONSERV VIOLATIONS", "", violations.c_str());
}


void alignment_block::compute_first_sequence()
{
    vector<sequence>::const_iterator seq_itr = sequences.begin();
    unsigned int len = seq_itr->getText().length();
    
    m_reference.setSrc("REFERENCE (FIRST)");

    for(unsigned int col=0; col<len; ++col)
    {
	seq_itr = sequences.begin();
	for(; seq_itr < sequences.end(); ++seq_itr)
	{
	    char ch = seq_itr->getText()[col];
	    
	    if (ch != '-')
	    {
		m_reference.appendToText(ch);
		break;
	    }
	}
    }
}


void alignment_block::compute_consensus_sequence()
{
    vector<sequence>::const_iterator seq_itr = sequences.begin();
    unsigned int len = seq_itr->getText().length();

    m_reference.setSrc("REFERENCE (CONSENSUS)");
    		       
    for(unsigned int col=0; col<len; ++col)
    {
	int max = 0;
	map<char,int> freq;
	seq_itr = sequences.begin();
	for(; seq_itr < sequences.end(); ++seq_itr)
	{
	    char ch = seq_itr->getText()[col];

	    if (ch == '-')
		continue;
	    
	    if (freq.find(ch) == freq.end())
		freq[ch] = 1;
	    else
		freq[ch]++;

	    if (freq[ch] > max)
		max = freq[ch];
	}

	seq_itr = sequences.begin();
	for(; seq_itr < sequences.end(); ++seq_itr)
	{
	    char ch = seq_itr->getText()[col];

	    if (ch != '-' && freq[ch] == max)
	    {
		m_reference.appendToText(ch);
		break;
	    }
	}
    }
}


bool alignment_block::is_duplicate(alignment_block &block) const
{
    if (sequences.size() != block.sequences.size())
	return false;
    
    vector<sequence>::const_iterator this_itr = sequences.begin();
    vector<sequence>::const_iterator that_itr = block.sequences.begin();
    
    for(; this_itr < sequences.end() && that_itr < block.sequences.end();
	++this_itr, ++that_itr)
    {
	if (this_itr->getText() != that_itr->getText())
	    return false;
    }
    return true;
}


bool alignment_block::is_duplicate_column(const unsigned int i,
					  const unsigned int j) const
{
    vector<sequence>::const_iterator itr = sequences.begin();
    
    for(; itr < sequences.end(); ++itr)
    {
	if (itr->getText()[i] != itr->getText()[j])
	    return false;
    }
    return true;   
}

bool alignment_block::is_majority_gaps_column(const unsigned int i) const
{
    int num_gaps = 0;
    int num_non_gaps = 0;
    vector<sequence>::const_iterator itr = sequences.begin();
    for(; itr < sequences.end(); ++itr)
    {
	if (itr->getText()[i] == '-')
	    num_gaps++;
	else
	    num_non_gaps++;
    }
    if (num_gaps > num_non_gaps)
	return true;

    return false;
}

double alignment_block::score(const alignment_block &block) const
{
    if (sequences.size() != block.sequences.size())
    {
	fprintf(stderr, "FATAL ERROR. Attempting to score blocks "
		"of unequal sizes. Exiting.\n");
	exit(-1);      
    }
    
    vector<sequence>::const_iterator this_itr = sequences.begin();
    vector<sequence>::const_iterator that_itr = block.sequences.begin();

    int nucl_change = 0;
    int num_nucl = 0;

    for(; this_itr < sequences.end() && that_itr < block.sequences.end();
	++this_itr, ++that_itr)
    {
	const string &this_seq = this_itr->getText();
	const string &that_seq = that_itr->getText();

	if (this_seq.length() != that_seq.length())
	{
	    fprintf(stderr, "FATAL ERROR. Attempting to score blocks "
		    "containing sequences of unequal lengths. Exiting.\n");
	    exit(-1);
	}

	unsigned int len = this_seq.length();

	for(unsigned int i = 0; i<len; ++i)
	{
	    if (this_seq[i] == '-')
		continue;

	    if (this_seq[i] != that_seq[i])
		nucl_change++;

	    num_nucl++;
	}
    }

    double score = (double) nucl_change / (double) num_nucl;

    return score;
}

