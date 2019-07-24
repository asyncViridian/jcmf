#include "clustalw_file.h"
#include <strings.h>

using namespace std;

clustalw_file::clustalw_file(const char *filename) :
    base_file(filename)
{
    m_file = fopen(filename, "r");

    if (!m_file)
    {
	fprintf(stderr, "FATAL ERROR. Unable to open Clustal W file '%s'. "
		"Exiting.\n", filename);
	exit(-1);
    }
}

clustalw_file::~clustalw_file()
{
    fclose(m_file);
}

bool clustalw_file::get_next_block(alignment_block &block,
				   string &comments)
{
    char buf[MAX_LINE_SIZE];

    block.clear();
    block.set_type(alignment_block::MULT_ALIGN_CLUSTALW);
    comments.clear();

    if (!fgets(buf, sizeof(buf), m_file))
	return false; // Second invocation of get_next_block(), exiting.

    if (strncasecmp(buf, "CLUSTAL", 7))
    {
	fprintf(stderr, "FATAL ERROR. Invalid format for Clustal W file. "
		"The file does not begin with 'CLUSTAL'. Exiting.\n");
	exit(-1);
    }
    
    comments.append(buf);

    bool start = true;
    while(fgets(buf, sizeof(buf), m_file) != NULL)
    {
	if (start && isspace(buf[0]))
	{
	    comments.append(buf);
	    continue;
	}
	else if (isspace(buf[0]))
	{
	    break;
	}

	start = false;
	
	char src[MAX_LINE_SIZE];
	int start = 0;
	int size = 0;
	char strand = ' ';
	int srcSize = 0;
	char text[MAX_LINE_SIZE];
	
	int num_elements = sscanf(buf, "%s %s", src, text);

	if (num_elements != 2)
	{
	    fprintf(stderr,
		    "FATAL ERROR. Invalid Clustal W format, line:\n%s\n", buf);
	    exit(-1);
	}
	
	sequence seq(src, start, size, strand, srcSize, text);
	
	block.sequences.push_back(seq);
    }
    

    unsigned int seq_num = 0;
    while(fgets(buf, sizeof(buf), m_file) != NULL)
    {
	if (isspace(buf[0]))
	{
	    seq_num = 0;
	    continue;
	}

	char src[MAX_LINE_SIZE];
	char text[MAX_LINE_SIZE];
	
	sscanf(buf, "%s %s", src, text);

	if (block.sequences[seq_num].getSrc() != src)
	{
	    fprintf(stderr, "FATAL ERROR. Bad format for Clustal W file. "
		    "Expecting '%s', read '%s'. Exiting.\n",
		    block.sequences[seq_num].getSrc().c_str(), src);
	    exit(-1);
	}

	block.sequences[seq_num].appendToText(text);
	seq_num++;
    }

    block.finalize();

    return true;
}

