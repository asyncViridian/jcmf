#ifndef __MAF_FILE_H__
#define __MAF_FILE_H__

#include <string>
#include <cstring>
#include <vector>
#include "alignment_block.h"
#include "base_file.h"

#define MAX_LINE_SIZE (10*1024)


class maf_file : public base_file
{
 public:
    maf_file(const char *filename);

    ~maf_file();

    bool get_next_block(alignment_block &block,
			std::string &comments);
    
    static void read(const char *filename,
		     std::vector<alignment_block> &blocks,
		     const char *src_prefix = NULL);

 private:
    FILE *m_file;
};

#endif /* __MAF_FILE_H__ */
