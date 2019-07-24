#ifndef __CHI_SQUARE_H__
#define __CHI_SQUARE_H__

#include "sequence.h"

class chi_square
{
 public:
    static void compute(const sequence &seq_a, const sequence &seq_b,
			double &chisq, double &degrees_of_freedom,
			double &pvalue);
};

#endif /* __CHI_SQUARE_H__ */
