#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <strings.h>
#include "options.h"
#include "config.h"

using namespace std;

struct option options::g_long_options[] = {
    {"clustalw", 1, 0, 'w'},
    {"num", 1, 0, 'n'},
    {"refseq", 1, 0, 'r'},
    {"conservation", 1, 0, 'c'},
    {"nomajgaphandling", 0, 0, 'g'},
    {"norandomoffset", 0, 0, 'o'},
    {"discardgapseqs", 0, 0, 'd'},
    {"expscore", 0, 0, 'e'},
#if USE_RNAZ
    {"runrnaz", 0, 0, 'z'},
#endif
    {"verbose", 0, 0, 'v'},
    {"version", 0, 0, 's'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
};


#if USE_RNAZ
char options::g_short_options[] = "wn:r:c:godezvsh";
#else
char options::g_short_options[] = "wn:r:c:godevsh";
#endif


options::options() :
    m_clustalw(false),
    m_num(1),
    m_discardgapseqs(false),
    m_experimental_score(false),
    m_runrnaz(false),
    m_refseq(REFSEQ_CONSENSUS),
    m_conservation(CONSERVATION_LEVEL1),
    m_nomajgaphandling(false),
    m_norandomoffset(false),
    m_verbose(false),
    m_maf_filename("")
{
}
    

void options::parse(int argc, char * const argv[])
{
    // We'll handle the error messages
    opterr = 0;
    
    while(1)
    {
	int options_index;
	int c = getopt_long(argc, argv, g_short_options,
			    g_long_options, &options_index);

	if (c == -1)
	    break;
	
	switch(c)
	{
	case 'w':
	    m_clustalw = true;
	    m_verbose = false;
	    break;
	    
	case 'n':
	    if (optarg)
		m_num = atoi(optarg);
	    break;
		
	case 'r':
	    if (optarg && !strcasecmp(optarg, "first"))
		m_refseq = REFSEQ_FIRST;
	    else if (optarg && !strcasecmp(optarg, "consensus"))
		m_refseq = REFSEQ_CONSENSUS;
	    else
		print_usage_and_exit(-1);
	    break;

	case 'c':
	    if (optarg && !strcasecmp(optarg, "none"))
		m_conservation = CONSERVATION_NONE;
	    else if (optarg && !strcasecmp(optarg, "gaps"))
		m_conservation = CONSERVATION_GAPS;
	    else if (optarg && !strcasecmp(optarg, "full"))
		m_conservation = CONSERVATION_FULL;	    
	    else if (optarg && !strcasecmp(optarg, "level0"))
		m_conservation = CONSERVATION_LEVEL0;
	    else if (optarg && !strcasecmp(optarg, "level1"))
		m_conservation = CONSERVATION_LEVEL1;
	    else
		print_usage_and_exit(-1);
	    break;

	case 'g':
	    m_nomajgaphandling = true;
	    break;

	case 'o':
	    m_norandomoffset = true;
	    break;
	    
	case 'd':
	    m_discardgapseqs = true;
	    break;

	case 'e':
	    if (!m_clustalw)
	    {
		m_experimental_score = true;
		m_verbose = true;
	    }
	    break;

#if USE_RNAZ
	case 'z':
	    if (!m_clustalw)
	    {
		m_runrnaz = true;
		m_verbose = true;
	    }
	    break;
#endif
	    
	case 'v':
	    if (!m_clustalw)
	    {
		m_verbose = true;
	    }
	    break;

        case 's':
            print_version();
            exit(0);

	case 'h':
	    print_usage_and_exit(0);
	    break;
	    
	case '?':
	default:
	    print_usage_and_exit(-1);
	    break;
	}
    }

    if (optind != argc - 1)
	print_usage_and_exit(-1);

    m_maf_filename = argv[optind];
}

void options::print_version() const
{
    const char *version_string =
    "\n"
    "Multiple Sequence Alignment Random Shuffle\n"
    PACKAGE_STRING "\n"
    "Report bugs to: " PACKAGE_BUGREPORT "\n"
    "\n";
    fprintf(stderr, "%s", version_string);
}

void options::print_usage_and_exit(int exit_code) const
{
    const char *usage =
    "\n"
    "Usage:  multiperm [OPTIONS] ALIGNMENT_FILE\n"
    "\n"
    "Input:  A file containing (one or more) multiple sequence alignments\n"
    "        in MAF format. If the '-w' option is specified, the file should\n"
    "        be in Clustal W format instead.\n"
    "\n"
    "Output: N files named perm_XXX_[ALIGNMENT_FILE] where [ALIGNMENT_FILE]\n"
    "        is the input filename and XXX is a number. The output files are\n"
    "        in MAF format unless the '-w' option is specified, in which case\n"
    "        they are in Clustal W format.\n"
    "\n"	
    "Options:\n"
    "-w, --clustalw            The input file contains a multiple sequence\n"
    "                          alignment in Clustal W format. Turns off\n"
    "                          'verbose'.\n"
    "\n"
    "-n, --num=N               Number of times the random shuffle is\n"
    "                          performed for each multiple alignment\n"
    "                          This determines the number of output files.\n"
    "                          Default n: 1\n"
    "\n"
    "-r, --refseq=REFSEQ       Use  --refseq=first  to perform the\n"
    "                          shuffle with the first sequence in the\n"
    "                          multiple alignment.\n"
    "                          Use  --refseq=consensus  to perform the\n"
    "                          shuffle with the consensus sequence.\n"
    "                          Default refseq: consensus\n"
    "\n"
    "-c, --conservation=CONS   Use  --conservation=none  to not impose any\n"
    "                          conservation restrictions on the shuffle.\n"
    "                          Use  --conservation=gaps  to only conserve\n"
    "                          gap structure.\n"
    "                          Use  --conservation=full  to preserve local\n"
    "                          conservation patterns in addition to gap\n"
    "                          structure.\n"
    "                          Use  --conservation=level0  to coarse-grain\n"
    "                          preserve local conservation patterns and\n"
    "                          strictly preserve gap structure (modeled\n"
    "                          after rnazRandomizeAln.pl's level 0)\n"    
    "                          Use  --conservation=level1  to coarse-grain\n"
    "                          preserve local conservation patterns and\n"
    "                          strictly preserve gap structure (modeled\n"
    "                          after rnazRandomizeAln.pl's level 1)\n"
    "                          Default: level1\n"
    "\n"
    "-g, --nomajgaphandling    Turn off the special handling of majority\n"
    "                          gap regions of the multiple alignment.\n"
    "\n"
    "-o, --norandomoffset      Turn off random column offsets. Makes all\n"
    "                          shuffles start at the first column.\n"
    "\n"
    "-d, --discardgapseqs      Discard all multiple alignments in the\n"
    "                          input file that have a gap in any sequence.\n"
    "                          Default: no discarding.\n"
    "\n"
    "-e, --expscore            Perform 1000 random permutations to\n"
    "                          obtain the experimental shuffling score.\n"
    "                          Turns on 'verbose'.\n"
    "\n"
#if USE_RNAZ
    "-z, --runrnaz             Run 'rnazRandomizeAln.pl -l1' on the same\n"
    "                          alignments, for comparison. Turns on 'verbose'.\n"
    "\n"
#endif
    "-v, --verbose             Output additional information as MAF comments.\n"
    "\n"
    "-s, --version             Print package version number and exit.\n"
    "\n"
    "-h, --help                This help screen.\n"	
    "\n";

    print_version();

    fprintf(stderr, "%s", usage);

    exit(exit_code);
}
