bin_PROGRAMS = multiperm

multiperm_SOURCES = mt19937ar.h mt19937ar.cpp prng.h prng.cpp options.h options.cpp sequence.h sequence.cpp alignment_block.h alignment_block.cpp rnaz_random_shuffle.h rnaz_random_shuffle.cpp chi_square.h chi_square.cpp base_file.h maf_file.h maf_file.cpp clustalw_file.h clustalw_file.cpp our_shuffle.h our_shuffle.cpp multiperm.h multiperm.cpp shuffling_score.h shuffling_score.cpp main.cpp 

multiperm_CPPFLAGS = -O3 -g -Wall

if USE_R
multiperm_CPPFLAGS += -I$(RHOME)/include -DRHOME=\"$(RHOME)\" -DUSE_R
multiperm_LDFLAGS = -L$(RHOME)/lib -lRblas -lR
endif

if USE_RNAZ
multiperm_CPPFLAGS += -DRNAZ_RANDOMIZE_ALN=\"$(RNAZ_RANDOMIZE_ALN)\" -DUSE_RNAZ
endif
