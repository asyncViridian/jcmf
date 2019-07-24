#ifndef __CLUSTALW_FILE_H__
#define __CLUSTALW_FILE_H__

#include <string>
#include <vector>
#include "alignment_block.h"
#include "base_file.h"

#define MAX_LINE_SIZE (10*1024)


class clustalw_file : public base_file
{
 public:
    clustalw_file(const char *filename);

    ~clustalw_file();

    bool get_next_block(alignment_block &block,
			std::string &comments);
    
  private:
    FILE *m_file;
};

#endif /* __CLUSTALW_FILE_H__ */
