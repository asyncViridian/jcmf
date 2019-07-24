#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include <string>

typedef enum
{
    REFSEQ_FIRST,
    REFSEQ_CONSENSUS
}
refseq_t;

typedef enum
{
    CONSERVATION_NONE,
    CONSERVATION_GAPS,
    CONSERVATION_FULL,
    CONSERVATION_LEVEL0,
    CONSERVATION_LEVEL1
}
conservation_t;

class options
{
 public:
    options();

    bool clustalw() {return m_clustalw;}
    int num() {return m_num;}
    bool discardgapseqs() {return m_discardgapseqs;}
    bool experimentalscore() {return m_experimental_score;}
    bool runrnaz() {return m_runrnaz;}
    refseq_t refseq() {return m_refseq;}
    conservation_t conservation() {return m_conservation;}
    bool nomajgaphandling() {return m_nomajgaphandling;}
    bool norandomoffset() {return m_norandomoffset;}
    bool verbose() {return m_verbose;}
    std::string &maf_filename() {return m_maf_filename;}
    
    void parse(int argc, char * const argv[]);
    
 private:
    void print_version() const;
    void print_usage_and_exit(int exit_code) const;

    bool m_clustalw;
    int m_num;
    bool m_discardgapseqs;
    bool m_experimental_score; 
    bool m_runrnaz;
    refseq_t m_refseq;
    conservation_t m_conservation;
    bool m_nomajgaphandling;
    bool m_norandomoffset;
    bool m_verbose;
    std::string m_maf_filename;
    
    static struct option g_long_options[];
    static char g_short_options[];
};


#endif /* __OPTIONS_H__ */
