#include<map>
#include<string>
#include "sequence.h"
#include "chi_square.h"
#include "prng.h"


sequence::sequence() :
    m_src(""),
    m_start(0),
    m_size(0),
    m_strand('+'),
    m_srcSize(0),
    m_text("")
{
}

void sequence::clear()
{
    m_src = "";
    m_start = 0;
    m_size = 0;
    m_strand = '+';
    m_srcSize = 0;
    m_text = "";
}


sequence::sequence(std::string src, int start, int size,
		   char strand, int srcSize, std::string text) :
    m_src(src),
    m_start(start),
    m_size(size),
    m_strand(strand),
    m_srcSize(srcSize),
    m_text(text)
{
    int len = getText().length();
    
    for(int i=0; i<len; ++i)
    {
	m_text[i] = toupper(m_text[i]);

	// If the nucleotide is something like 'N', 'R', etc, just
	// randomly pick a nucleotide
	if (m_text[i] != 'A' && m_text[i] != 'C'
	    && m_text[i] != 'G' && m_text[i] != 'T'
	    && m_text[i] != '-')
	{
	    double rand_num = get_rand();
	    if (rand_num < 0.25)
		m_text[i] = 'A';
	    else if (rand_num < 0.5)
		m_text[i] = 'C';
	    else if (rand_num < 0.75)
		m_text[i] = 'G';
	    else
		m_text[i] = 'T';
	}
    }

    init_nucleotide_stats();
    init_dinucleotide_stats();    
}


inline const int sequence::id(const char di[2]) const
{
    switch (di[0])
    {
    case 'A':
	switch (di[1])
	{
	case 'A':
	    return 0;
	case 'C':
	    return 1;
	case 'G':
	    return 2;
	case 'T':
	    return 3;
	}
	break;
    case 'C':
	switch (di[1])
	{
	case 'A':
	    return 4;
	case 'C':
	    return 5;
	case 'G':
	    return 6;
	case 'T':
	    return 7;
	}
	break;	
    case 'G':
	switch (di[1])
	{
	case 'A':
	    return 8;
	case 'C':
	    return 9;
	case 'G':
	    return 10;
	case 'T':
	    return 11;
	}
    case 'T':
	switch (di[1])
	{
	case 'A':
	    return 12;
	case 'C':
	    return 13;
	case 'G':
	    return 14;
	case 'T':
	    return 15;
	}
    }

    // This should never happen, so abort!
    fprintf(stderr,
	    "ERROR: Tried to look up dinucleotide '%s', which doesn't exist!\n",
	    di);
    exit(-1);
}

inline const char *sequence::di(const int id) const
{
    switch (id)
    {
    case  0:
	return "AA";
    case  1:
	return "AC";
    case  2:
	return "AG";
    case  3:
	return "AT";
    case  4:
	return "CA";
    case  5:
	return "CC";
    case  6:
	return "CG";
    case  7:
	return "CT";
    case  8:
	return "GA";
    case  9:
	return "GC";
    case 10:
	return "GG";
    case 11:
	return "GT";
    case 12:
	return "TA";
    case 13:
	return "TC";
    case 14:
	return "TG";
    case 15:
	return "TT";
    }

    // This should never happen, so abort!
    fprintf(stderr,
	    "ERROR: Tried to look up dinucleotide id '%d', which doesn't exist!\n",
	    id);
    exit(-1);

}


void sequence::init_nucleotide_stats()
{
    // Initialize
    for(int i=0; i<NUM_NUCLEOTIDES; ++i)
	m_nucl_count[i] = 0;
    
    const char *seq_str = getText().c_str();
    int len = getText().length();
    
    for(int i=0; i<len; ++i)
    {
	switch (seq_str[i])
	{
	case 'A':
	    m_nucl_count[0]++;
	    break;
	case 'C':
	    m_nucl_count[1]++;
	    break;
	case 'G':
	    m_nucl_count[2]++;
	    break;
	case 'T':
	    m_nucl_count[3]++;
	    break;
	case '-':
	    // ignore
	    break;
	default:
	    fprintf(stderr, "FATAL ERROR: unknown nucleotide '%c'\n",
		    seq_str[i]);
	    exit(-1);
	    break;
	}
    }
}


void sequence::init_dinucleotide_stats()
{
    // We add a pseudocount to the dinucleotide count to avoid
    // having zeroes that mess up the chisq statistics
    const int pseudocount = 1;

    // Initialize
    for(int i=0; i<NUM_DINUCLEOTIDES; ++i)
	m_dinucl_count[i] = pseudocount;

    char dinucleotide[3];

    const char *seq_str = getText().c_str();

    int len = getText().length();
    
    for(int i=0; i<len-1; ++i)
    {
	while (seq_str[i] == '-' && i<len-1)
	    i++;
	if (i >= len-1)
	    break;
	dinucleotide[0] = seq_str[i];

	while (seq_str[i+1] == '-' && i<len-1)
	    i++;
	if (i >= len-1)
	    break;
	dinucleotide[1] = seq_str[i+1];
	
	dinucleotide[2] = '\0';

	m_dinucl_count[id(dinucleotide)]++;
    }
}


void sequence::print_stats(FILE *fout) const
{
    fprintf(fout, "# %25s %12s ", getSrc().c_str(), ""); 

    for(int i=0; i<NUM_DINUCLEOTIDES; ++i)
    {
	fprintf(fout, "%s:%3d ", di(i), m_dinucl_count[i]);
    }

    return;
}


void sequence::print_stats(FILE *fout, const sequence &seq) const
{
    fprintf(fout, "# %25s %10d %c ", getSrc().c_str(), getStart(), getStrand());
    
    double chisq = 0.0;
    double df = 0;
    double pvalue = 0.0;
    
    chi_square::compute(*this, seq, chisq, df, pvalue);

    if (df >= 0.0)
	fprintf(fout, "%6.2f  %2d  %1.3f   ", chisq, (int) df, pvalue);
    else
	fprintf(fout, "%6.2f   -    -     ", chisq);

    for(int i=0; i<NUM_DINUCLEOTIDES; ++i)
    {
	fprintf(fout, "%3d/%3d ", getDinuclCount()[i], seq.getDinuclCount()[i]);
    }
    fprintf(fout, " |  ");
    for(int i=0; i<NUM_NUCLEOTIDES; ++i)
    {
	fprintf(fout, "%3d/%3d ", getNuclCount()[i], seq.getNuclCount()[i]);
    }
    fprintf(fout, "\n");
    
    return;
}


void sequence::print_dinucleotide_header(FILE *fout, const char *src_label) const
{
    fprintf(fout, "# %25s %12s ", src_label, "");

    fprintf(fout, " chisq  df  pvalue  ");
    
    for(int i=0; i<NUM_DINUCLEOTIDES; ++i)
    {
	fprintf(fout, "   %s   ", di(i));
    }
    fprintf(fout, " |     A       C       G       T    \n");
    
    return;    
}

void sequence::getMAFline(char *line) const
{
    sprintf(line, "s %s %d %d %c %d %s\n", m_src.c_str(), m_start,
	    m_size, m_strand, m_srcSize, m_text.c_str());

}
