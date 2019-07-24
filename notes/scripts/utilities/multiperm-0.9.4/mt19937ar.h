/* mt19937ar.h
**
**
*/

#ifndef __MT19937AR_H__
#define __MT19937AR_H__

void init_genrand(unsigned long s);

void init_by_array(unsigned long init_key[], int key_length);

unsigned long genrand_int32(void);

double genrand_real2(void);

#endif /* __MT19937AR_H__ */
