package cmf;

/**
 *
 * cmfinder04.pl
 */
public class cmfinder {

    // default parameters
    int CAND = 40;
    int MAXSPAN1 = 100;
    int MINSPAN1 = 30;
    int MAXSPAN2 = 100;
    int MINSPAN2 = 40;
    int CLUSTER = 3;
    double FRACTION = 0.8;
    int SINGLE = 5;
    int DOUBLE = 5;
    int verbose = 0;
    int help = 0;
    int COMBINE = 0;
    String cand_weight_option = "";
    String cmfinderBaseExe = "cmfinder04";
    // ($likeold,$skipClustalw,$useOldCmfinder,$simpleMotifsAlreadyDone,$justGetCmfinderCommand,
    // $copyCmfinderRunsFromLog,$amaa,$version,$filterNonFrag,$fragmentary,$commaSepEmFlags,
    // $commaSepSummarizeFlags,$commaSepCandfFlags,$saveTimer,$allCpus,$cpu,$candsParallel,$outFileSuffix,
    // $columnOnlyBasePairProbs);
    int emulate_apparent_bug_in_resolve_overlap = 1;
    
    //parameters from comb_motif.pl

}
