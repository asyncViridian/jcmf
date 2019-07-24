#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include<string>

#define NUM_NUCLEOTIDES 4
#define NUM_DINUCLEOTIDES 16

class sequence
{
 public:

    sequence();
    
    sequence(std::string src,   // "source", typically species & chromosome from which the sequence derives
	     int start,         // start location within source (1-based, perhaps, but doesn't matter here)
	     int size,          // number of nucleotides; excluding alignment gaps
	     char strand,       // from which strand? (+/-)
	     int srcSize,       // may be total size of src seq, but is often 0 (doesn't matter here)
	     std::string text   // the actual nuc sequence, *including* gaps
	     );

    void clear();
	    
    const int id(const char di[2]) const;
    const char *di(const int id) const;
    
    void print_stats(FILE *fout) const;

    void print_stats(FILE *fout, const sequence &seq) const;

    void print_dinucleotide_header(FILE *fout, const char *src_label = "") const;

    const std::string &getSrc() const {return m_src;}
    const int getStart() const {return m_start;}
    const int getSize() const {return m_size;}
    const char getStrand() const {return m_strand;}
    const int getSrcSize() const {return m_srcSize;}
    const std::string &getText() const {return m_text;}

    void setSrc(std::string src) {m_src = src;}
    void setStart(int start) {m_start = start;}
    void setSize(int size) {m_size = size;}
    void setStrand(char strand) {m_strand = strand;}
    void setSrcSize(int srcSize) {m_srcSize = srcSize;}
    void setText(std::string text) {m_text = text;} 
    
    void appendToText(const char ch) {m_text.append(1,ch);}
    void appendToText(const char *str) {m_text.append(str);}

    void eraseInText(const int col, const int n) {m_text.erase(col, n);}
    void replaceInText(const int col, const char ch) {m_text[col] = ch;}

    const int *getNuclCount() const {return m_nucl_count;}
    const int *getDinuclCount() const {return m_dinucl_count;}

    void getMAFline(char *line) const;

    void init_nucleotide_stats();
    void init_dinucleotide_stats();
    
 private:

    std::string m_src;
    int m_start;
    int m_size;
    char m_strand;
    int m_srcSize;
    std::string m_text;

    // dinucleotide statistics
    int m_nucl_count[NUM_NUCLEOTIDES];
    int m_dinucl_count[NUM_DINUCLEOTIDES];
    
};


#endif /* __SEQUENCE_H__ */
