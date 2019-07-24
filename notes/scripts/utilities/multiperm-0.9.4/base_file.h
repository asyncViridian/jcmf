/* base_file.h
**
**
*/

#ifndef __BASE_FILE_H__
#define __BASE_FILE_H__


class base_file
{
 public:
    base_file(const char *filename) {};

    virtual ~base_file() {};

    virtual bool get_next_block(alignment_block &block,
				std::string &comments) = 0;
};



#endif /* __BASE_FILE_H__ */
