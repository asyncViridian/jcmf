#include "maf_file.h"

using namespace std;

maf_file::maf_file(const char *filename) :
    base_file(filename)
{
    m_file = fopen(filename, "r");

    if (!m_file)
    {
	fprintf(stderr, "FATAL ERROR. Unable to open MAF file '%s'. Exiting.\n", filename);
	exit(-1);
    }
}

maf_file::~maf_file()
{
    fclose(m_file);
}

bool maf_file::get_next_block(alignment_block &block,
			      string &comments)
{
    char buf[MAX_LINE_SIZE];

    block.clear();
    comments.clear();

    while(fgets(buf, sizeof(buf), m_file) != NULL)
    {
	if (buf[0] == '\n')
	{
	    if (block.sequences.size() > 0)
	    {
		// We're done with this block
		block.finalize();
		return true;
	    }
	}
	else if (buf[0] == 's')
	{
	    char type;
	    char src[MAX_LINE_SIZE];
	    int start;
	    int size;
	    char strand;
	    int srcSize;
	    char text[MAX_LINE_SIZE];

	    int num_elements = sscanf(buf, "%c %s %d %d %c %d %s",
				&type, src, &start, &size,
				&strand, &srcSize, text);

	    if (num_elements != 7)
	    {
		fprintf(stderr, "FATAL ERROR. Invalid MAF format, line:\n%s\n",
			buf);
		exit(-1);
	    }

	    sequence seq(src, start, size, strand, srcSize, text);

	    block.sequences.push_back(seq);
	}
	else
	{
	    comments.append(buf);
	}
    }

    return false;
}

void maf_file::read(const char *filename,
		    vector<alignment_block> &blocks,
		    const char *src_prefix)
{
    FILE *fin = fopen(filename, "r");

    if (!fin)
    {
	fprintf(stderr, "FATAL ERROR. Unable to open MAF file '%s'. Exiting.\n", filename);
	exit(-1);
    }
    
    char buf[MAX_LINE_SIZE];
    
    alignment_block block;

    while(fgets(buf, sizeof(buf), fin) != NULL)
    {
	if (buf[0] == '\n')
	{
	    // We're done with this block
	    block.finalize();
	    blocks.push_back(block);
	    block.clear();
	    continue;
	}
	else if (buf[0] == 's')
	{
	    char type;
	    char src[MAX_LINE_SIZE];
	    int start;
	    int size;
	    char strand;
	    int srcSize;
	    char text[MAX_LINE_SIZE];

	    int prefix_len = 0;
	    if (src_prefix)
	    {
		prefix_len = strlen(src_prefix);
		strncpy(src, src_prefix, prefix_len);
	    }
	    
	    sscanf(buf, "%c %s %d %d %c %d %s",
		   &type, (src + prefix_len), &start, &size,
		   &strand, &srcSize, text);

	    sequence seq(src, start, size, strand, srcSize, text);

	    block.sequences.push_back(seq);
	}
	
    }
    
    fclose(fin);
}
