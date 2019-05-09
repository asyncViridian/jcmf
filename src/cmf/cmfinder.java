package cmf;

/**
 *
 * cmfinder04.pl
 */
import static cmf.Io.*;
import static cmf.utilities.*;
import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.json.*;

public class cmfinder {

    public static String bin_path = "cmfinder/bin";
    // default parameters
    static String CMFINDER_PACKAGE_VERSION = "0.4.1.15";
    static int CAND = 40;
    static int MAXSPAN1 = 100;
    static int MINSPAN1 = 30;
    static int MAXSPAN2 = 100;
    static int MINSPAN2 = 40;
    static int CLUSTER = 3;
    static double FRACTION = 0.8;
    static int SINGLE = 5;
    static int DOUBLE = 5;
    static boolean verbose = true;
    static int help = 0;
    static boolean COMBINE = false;
    static String cand_weight_option = "";
    static String cmfinderBaseExe = "cmfinder04";
    static boolean emulate_apparent_bug_in_resolve_overlap = false;

    //parameters from comb_motif.pl
    static int comb_max_gap = 100; // motif instances within this distance of each other can be considered close enough for merging
    static double comb_len_energy_threshold = 0.1;
    static int comb_min_overlap = 2;
    static double comb_min_num = 2.5;
    static int comb_max_len = 200;
    static boolean output_more_like_old_cmfinder_pl = false; //perl using 0, debug flag
    static String motifList;
    static int minCandScoreInFinal = 0; // be more like old cmfinder for now
    static String emSeq;
    static boolean skipClustalw;
    static boolean useOldCmfinder; // set in the static block or main method
    static String outFileSuffix;

    //getOptions parameters from json file
    static JSONObject jo;

    static ArrayList<String> cmfinder_inf11FlagsList = new ArrayList<String>();
    static ArrayList<String> summarizeFlagsList = new ArrayList<String>();

    static String cmfinder_inf11Flags = "";
    static String summarizeFlagsStr = "";
    static String candfExtraFlags = "";
    static String saveTimerFlag = "";
    static String saveTimer03Flag = "";

    static HashMap<String, Seq> unaligned_seqs;

    //usage for aligments: try_merge
    static HashMap<String, Alignment> alignments = new HashMap<>();

    //SEQ setting in main method
    static String SEQ;
    static String seqForExpectationMaximization = SEQ;
    static String dummyCmfileParamForCmfinder = "";
    static String tempFileListFileName = SEQ + ".file-list";

    static ArrayList<String> cands = new ArrayList();

    static {

        unaligned_seqs = read_fasta(seqForExpectationMaximization);
        try {
            jo = read_json_file("./src/cmf/cmfinder_param.json");

            //get some value
            useOldCmfinder = jo.optBoolean("useOldCmfinder");
            motifList = jo.optString("motifList");
            COMBINE = jo.getBoolean("combine");
            skipClustalw = jo.getBoolean("skipClustalw");

            if ((!jo.optString("emSeq").isEmpty())
                      && (!jo.optString("emSeq").equals(null))) {
                seqForExpectationMaximization = jo.optString("emSeq").toString();
            }

            //set value outFileSuffix
            if ((outFileSuffix.isEmpty()) || (outFileSuffix.equals(null))) {
                outFileSuffix = "";
            } else {
                outFileSuffix = jo.optString("outFileSuffix=s");
            }

            if (jo.optBoolean("useOldCmfinder")) {
                cmfinderBaseExe = "cmfinder";
                dummyCmfileParamForCmfinder = SEQ + ".temp.cm";
            }

            int cpu = jo.getInt("cpu");
            cpu = jo.optBoolean("allCpus") ? -1 : cpu;
            if (cpu != 0) {
                cmfinder_inf11FlagsList.add("--cpu " + cpu);
            }
            /*
            skip following perl
            if ($likeold && !$useOldCmfinder) {
                push @cmfinder_inf11FlagsList,"--enone";
                push @cmfinder_inf11FlagsList,"--p56";
                push @cmfinder_inf11FlagsList,"--degen-rand";
                push @cmfinder_inf11FlagsList,"--ints-like-03";
                push @cmfinder_inf11FlagsList,"--min-seq-weight 0";
                push @cmfinder_inf11FlagsList,"--no-elim-iden-seq";
                push @cmfinder_inf11FlagsList,"--no-elim-iden-subseq";
                push @cmfinder_inf11FlagsList,"--min-cand-score-in-final 0";
            }
             */
            cmfinder_inf11FlagsList.add("--min-cand-score-in-final " + minCandScoreInFinal);
            if (jo.optBoolean("amma")) {
                cmfinder_inf11FlagsList.add("--amaa");
            }
            if (jo.optBoolean("filterNonFrag")) {
                cmfinder_inf11FlagsList.add("--filter-non-frag");
            }
            cmfinder_inf11Flags = String.join(" ", cmfinder_inf11FlagsList);
            //
            Iterator<?> key1 = jo.getJSONObject("commaSepEmFlags").keys();
            while (key1.hasNext()) {
                String key = (String) key1.next();
                String value = jo.getJSONObject("commaSepEmFlags").get(key).toString();
                if (!value.equals("null") && !value.equals("false")) {
                    //check if key has <x> part ?
                    cmfinder_inf11Flags
                              = cmfinder_inf11Flags
                              + ((cmfinder_inf11Flags != null && !cmfinder_inf11Flags.isEmpty()) ? " " : "")
                              + ((key.indexOf("<") > 1 ? "--" + key.substring(0, key.indexOf("<") - 1) : "--" + key)
                              + " "
                              + value);
                }
            }
            if (jo.optBoolean("columnOnlyBasePairProbs")) {
                cmfinder_inf11Flags = cmfinder_inf11Flags + " --column-only-base-pair-probs";
            }
            System.out.println(cmfinder_inf11Flags);

            if (jo.optBoolean("fragmentary")) {
                cmfinder_inf11FlagsList.add("--fragmentary");
                summarizeFlagsList.add("--fragmentary");
            }
            summarizeFlagsStr = String.join(" ", summarizeFlagsList);
            Iterator<?> key2 = jo.getJSONObject("commaSepSummarizeFlags").keys();
            while (key2.hasNext()) {
                String key = (String) key2.next();
                String value = jo.getJSONObject("commaSepSummarizeFlags").get(key).toString();
                if (!value.equals("null") && !value.equals("false")) {
                    //check if key has <x> part ?
                    summarizeFlagsStr
                              = summarizeFlagsStr
                              + ((summarizeFlagsStr != null && !summarizeFlagsStr.isEmpty()) ? " " : "")
                              + ((key.indexOf("<") > 1 ? "--" + key.substring(0, key.indexOf("<") - 1) : "--" + key)
                              + " "
                              + value);
                }
            }
            System.out.println(summarizeFlagsStr);

            Iterator<?> key3 = jo.getJSONObject("commaSepCandfFlags").keys();
            while (key3.hasNext()) {
                String key = (String) key3.next();
                String value = jo.getJSONObject("commaSepCandfFlags").get(key).toString();
                if (!value.equals("null") && !value.equals("false")) {
                    //check if key has <x> part ?
                    candfExtraFlags
                              = candfExtraFlags
                              + ((candfExtraFlags != null && !candfExtraFlags.isEmpty()) ? " " : "")
                              + ((key.indexOf("<") > 1 ? "--" + key.substring(0, key.indexOf("<") - 1) : "--" + key)
                              + " "
                              + value);
                }
            }
            System.out.println(candfExtraFlags);

            cand_weight_option = (jo.optString("w=s").equals(null))
                      ? cand_weight_option : jo.optString("w=s");
            System.out.println(cand_weight_option);
            if (!cand_weight_option.equals("")) {
                cand_weight_option = "-w " + cand_weight_option;
            }

            //saveTimer is a file, such as "./src/cmf/test"
            if (!jo.optString("saveTimer").equals(null)
                      && !jo.optString("saveTimer").equals("")) {
                System.out.println(jo.optString("saveTimer"));
                //delete file first
                deleteFile(jo.optString("saveTimer"));
                saveTimerFlag = "--timer-append " + jo.optString("saveTimer");
                saveTimer03Flag = "-t " + jo.optString("saveTimer");
            }

            if (jo.optBoolean("justGetCmfinderCommand")) {
                //just output in screen
                String cmfinder_cmd = bin_path + "/"
                          + cmfinderBaseExe + " "
                          + cmfinder_inf11Flags + " "
                          + cand_weight_option;
                System.out.println("-justGetCmfinderCommand:" + cmfinder_cmd);
                //exit
                try {
                    throw new MyException("Exit due to justGetCmfinderCommand setting.");
                } catch (MyException ex) {
                    Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            //check motifList file, slight different from perl
            if ((!motifList.equals(null)) && (!motifList.isEmpty())) {
                //check file existing otherwise IOException 
                findFile(bin_path, motifList);
            }

            if (jo.optBoolean("copyCmfinderRunsFromLog")) {
                /* below find the SEQ's folder                   
                    my $dir=$SEQ;
                    $dir =~ s/\/[^\/]+$//g;
                    print "dir=$dir\n";                    
                 */
                String dir = Paths.get(SEQ).toAbsolutePath().getParent().toString();
                System.out.println("dir=" + dir);
                //int gotCmfinderCmd = 0;  this indicate is stupid, remove
                findFile(bin_path, jo.optString("copyCmfinderRunsFromLog"));

                Pattern p_cmfinder = Pattern.compile("^\\S*cmfinder(|_inf11)\\s");
                Pattern p_a = Pattern.compile("-a\\s(\\S+)\\s"); // get -a arhument \s is space
                Pattern p_o = Pattern.compile("-o\\s(\\S+)\\s"); //get output filename
                Pattern p4 = Pattern.compile("[.]cm[.]");

                try (Stream<String> stream = Files.lines(Paths.get(SEQ))) {
                    stream.forEach(line -> {
                        //System.out.println(line);
                        Matcher m1 = p_cmfinder.matcher(line);
                        if (m1.find()) {
                            String version = m1.group(0);
                            String inputMsa = "";
                            String inputFasta = "";
                            String outputFileBase = "";
                            Matcher m2 = p_a.matcher(line);
                            if (m2.find()) {
                                inputMsa = m2.group(0);
                            } else {
                                try {
                                    throw new MyException("debug: p2 no match, exit.");
                                } catch (MyException ex) {
                                    System.err.println(ex);
                                }
                            }
                            Matcher m3 = p_o.matcher(line);
                            if (m3.find()) {
                                //get filename part only
                                outputFileBase = new File(m3.group(0)).getName();
                            } else {
                                try {
                                    throw new MyException("debug: p2 no match, exit.");
                                } catch (MyException ex) {
                                    System.err.println(ex);
                                }
                            }
                            String[] paceSepList = line.split("\\s+");
                            String cm = paceSepList[paceSepList.length - 1]; //get last element of array
                            if (Pattern.compile("[.]cm[.]").matcher(cm).find()) {
                                if (version.equals("_inf11")) {
                                    try {
                                        throw new MyException("cmfinder_inf11 should not have cmfile");
                                    } catch (MyException ex) {
                                        System.err.println(ex);
                                    }
                                }
                                inputFasta = cm;
                            } else {
                                if (version.length() == 0) {
                                    try {
                                        throw new MyException("cmfinder (0.3) should have cmfile");
                                    } catch (MyException ex) {
                                        Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
                                    }
                                }
                                inputFasta = cm;
                            }
                            System.out.println("ORIGINAL CMD:" + line);
                            if (!fileExists(inputMsa)) {
                                System.out.println("\t" + "SKIPPING: input msa " + inputMsa
                                          + " doesn't exist, so I assume it isn't good, and anyway it'd be a hassle to run.");
                            } else {
                                String output_motif_file = bin_path + "/" + cmfinderBaseExe;
                                String[] cmd = {output_motif_file,
                                    cmfinder_inf11Flags, cand_weight_option,
                                    "-o " + dir + "/" + outputFileBase,
                                    "-a " + inputMsa, inputFasta};
                                try {
                                    RunCmfinder(cmd, output_motif_file);
                                } catch (MyException ex) {
                                    Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
                                }
                                //gotCmfinderCmd = 1;  //lamdba can't change variable
                            }
                        } else {
                            try {
                                throw new MyException("didn't find any cmfinder commands");
                            } catch (MyException ex) {
                                Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    });
                }
                try {
                    //need end whole program here
                    throw new MyException("Exit due to copyCmfinderRunsFromLog.");
                } catch (MyException ex) {
                    Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
                }
            } // copyCmfinderRunsFromLog end

            if (!jo.getBoolean("simpleMotifsAlreadyDone")) {
                int runCandsParallel = 0;
                if (jo.getInt("cpu") >= 2
                          || jo.getBoolean("candsParallel")
                          || jo.getBoolean("allCpus")
                          || checkCPU() >= 2) {
                    runCandsParallel = 1;
                }

                HashMap<String, String[]> candsJobs = new HashMap<>();
                SINGLE = jo.getInt("s1=i");
                DOUBLE = jo.getInt("s2=i");
                if (SINGLE > 0) {
                    candsJobs.put("candf",
                              new String[]{bin_path + "/candf", candfExtraFlags, saveTimer03Flag, " -c " + CAND,
                                  " -o " + SEQ + ".cand.h1" + outFileSuffix, " -M " + MAXSPAN1, " -m " + MINSPAN1,
                                  " -s 1 ", " -S 1 ", SEQ});
                    candsJobs.put("cands",
                              new String[]{bin_path + "/cands", saveTimer03Flag,
                                  " -l " + tempFileListFileName + ".single" + outFileSuffix,
                                  " -n " + SINGLE, " -f " + FRACTION, SEQ, SEQ + ".cand.h1" + outFileSuffix});
                }
                if (DOUBLE > 0) {
                    candsJobs.put("candf",
                              new String[]{bin_path + "/candf", candfExtraFlags, saveTimer03Flag, " -c " + CAND,
                                  " -o " + SEQ + ".cand.h2" + outFileSuffix, " -M " + MAXSPAN2, " -m " + MINSPAN2,
                                  " -s 2 ", " -S 2 ", SEQ});
                    candsJobs.put("cands",
                              new String[]{bin_path + "/cands", saveTimer03Flag,
                                  " -l " + tempFileListFileName + ".single" + outFileSuffix,
                                  " -n " + DOUBLE, " -f " + FRACTION, SEQ, SEQ + ".cand.h2" + outFileSuffix});
                }

                //java 8: ExecutorService's invokeAll()
                if (runCandsParallel > 0) {
                    System.out.println(
                              "Parallel running "
                              + candsJobs.keySet().stream().collect(Collectors.joining(", "))
                    //join() will remove last delimiter, so no worry here
                    );
                    ExecutorService executor = Executors.newFixedThreadPool(checkCPU());
                    List<Callable<cmdOut>> callables
                              = //   remove below hard coded from map, using map to stream                       
                              //                                    = Arrays.asList(
                              //                                    () -> runCmd(candsJobs.get("candf")),
                              //                                    () -> runCmd(candsJobs.get("cands"))
                              //                          );
                              candsJobs.entrySet().stream()
                                        .map(e -> (Callable<cmdOut>) (() -> runCmd(e.getValue())))
                                        .collect(Collectors.toList());

                    executor.invokeAll(callables)
                              .stream()
                              .map(future -> {
                                  try {
                                      return future.get();
                                  } catch (Exception e) {
                                      throw new IllegalStateException(e);
                                  }
                              })
                              .forEach(System.out::println);

                    executor.shutdown();
                } else {
                    System.out.println(
                              "Serial running "
                              + candsJobs.keySet().stream().collect(Collectors.joining(", "))
                    //join() will remove last delimiter, so no worry here
                    );
                    //just run serially in stream
                    candsJobs.entrySet()
                              .stream().forEach(e -> runCmd(e.getValue()));
                }
                // need verify input arg is full path
                cands.add(tempFileListFileName + ".single" + outFileSuffix);
                cands.add(tempFileListFileName + ".double" + outFileSuffix);
            }

            ArrayList<String> motifFiles = new ArrayList();
            for (String cand : cands) {
                String pattern = SEQ + "(\\.\\d)?\\.cand.*\\_\\d+$";
                Pattern p = Pattern.compile(pattern);
                Matcher m = p.matcher(cand);
                if (m.find()) {
                    String pref = m.group(0);
                    String align = cand;
                    align = align.replace("cand", "align");
                    String motif = cand;
                    motif = motif.replace("cand", "motif");
                    if (jo.getBoolean("simpleMotifsAlreadyDone")) {
                        // perl logic error?
                        /*
                        if (-e $motif) {push @motifFiles,$motif;  }
                         */
                        if (motif.isEmpty()) {
                            motifFiles.add(motif);
                        }
                    } else {
                        String[] cmd1 = {bin_path + "/canda", saveTimer03Flag, cand, SEQ, align};
                        runCmd(cmd1);
                        String[] cmd2 = {bin_path + "/" + cmfinderBaseExe,
                            saveTimerFlag, cmfinder_inf11Flags, cand_weight_option,
                            "-o " + motif, "-a " + align,
                            seqForExpectationMaximization, dummyCmfileParamForCmfinder};
                        if (RunCmfinder(cmd2, motif)) {
                            //produced acceptable output
                            motifFiles.add(motif);
                            if ((!motifList.equals(null)) && (!motifList.isEmpty())) {
                                System.out.println("motifList=" + motif);
                            }
                        } else {
                            // oh well, nothing acceptable
                        }
                    }
                } else {
                    throw new MyException("weird cand file name " + cand);
                }
            }

            //combine
            if (COMBINE) {
                CombMotif(cand_weight_option, seqForExpectationMaximization, motifFiles);
            }

            //to do
        } catch (InterruptedException | MyException | JSONException | IOException ex) {
            System.out.println(ex);
        }
    }

    //read json file
    public static JSONObject read_json_file(String file_name) throws JSONException, IOException {
        Path p = Paths.get(file_name);
        //System.out.println(p.toAbsolutePath().toString());
        byte[] jsonByte = Files.readAllBytes(p.toAbsolutePath());
        JSONObject jsonObject = new JSONObject(new String(jsonByte));
        return jsonObject;
    }

    // arg: tempFileListFileName must with absolute path and file name
    public static ArrayList<String> GetCandFiles(String tempFileListFileName) throws MyException, IOException {
        ArrayList<String> list = new ArrayList();
        Stream<String> stream;
        if (!fileExists(tempFileListFileName)) {
            throw new MyException("cannot open " + tempFileListFileName);
        } else {
            stream = Files.lines(Paths.get(tempFileListFileName));
        }
        stream.forEach(e -> list.add(e));
        return list;
    }

    public static int resolve_overlap(Alignment alignment1, Alignment alignment2) {

        int cost1 = 0;
        int cost2 = 0;

        HashMap<String, AlignSeq> align1 = alignment1.getSeqs();
        HashMap<String, AlignSeq> align2 = alignment2.getSeqs();

        for (HashMap.Entry<String, AlignSeq> entry : align1.entrySet()) {
            if (align2.containsKey(entry.getKey())) {
                AlignSeq motif1 = entry.getValue();
                AlignSeq motif2 = align2.get(entry.getKey());
                String seq1 = motif1.getAlignSeq();
                String seq2 = motif2.getAlignSeq();
                String ss1 = motif1.getSs();
                String ss2 = motif2.getSs();
                HashMap<Integer, Integer> map1 = new HashMap<>();
                HashMap<Integer, Integer> map2 = new HashMap<>();
                if (emulate_apparent_bug_in_resolve_overlap) {
                    //do nothing in perl
                } else {
                    map1 = motif1.getAlignMap();
                    map2 = motif2.getAlignMap();
                }
                if (motif2.getStart() <= motif1.getEnd()) {
                    int olap_start = motif2.getStart();
                    int olap_end = motif1.getEnd();
                    HashMap<Integer, Integer> pt2 = pair_table(ss2);
                    int conflict2 = 0;
                    for (int i = 0; i < seq2.length(); i++) {
                        if (!map2.containsKey(i)) {
                            continue;
                        }
                        if (map2.get(i) > olap_end) {
                            break;
                        }
                        if (!pt2.containsKey(i)) {
                            continue;
                        }
                        if (pt2.get(i) >= 0) {
                            conflict2++;
                        }
                    }
                    HashMap<Integer, Integer> pt1 = pair_table(ss1);
                    int conflict1 = 0;
                    for (int i = seq1.length() - 1; i >= 0; i--) {
                        if (!map1.containsKey(i)) {
                            continue;
                        }
                        if (map1.get(i) < olap_start) {
                            break;
                        }
                        if (pt1.get(i) >= 0) {
                            conflict1++;
                        }
                    }
                    cost1 += conflict1 * motif1.getWeight();
                    cost2 += conflict2 * motif2.getWeight();
                }
            }

        }

        return (cost1 < cost2) ? 1 : 2;
    }

    public static Merge_Motif_Output merge_motif(AlignSeq motif1, AlignSeq motif2, String whole_seq, int olap_own) {
        //assume that $motif1 should be before $motif2;
        Merge_Motif_Output mm = new Merge_Motif_Output(0, 0, "", "", "", "", "", "");

        if (motif1 == null) {
            return new Merge_Motif_Output(
                      motif2.getStart(),
                      motif2.getEnd(),
                      "", "", "", "",
                      motif2.getAlignSeq(),
                      motif2.getAlignSs());
        }
        if (motif2 == null) {
            return new Merge_Motif_Output(
                      motif1.getStart(),
                      motif1.getEnd(),
                      "", "", "", "",
                      motif1.getAlignSeq(),
                      motif1.getAlignSs());
        }

        // Two motifs can't be merged.    
        if (motif1.getStart() > motif2.getStart()) {
            return new Merge_Motif_Output(-1, -1, "", "", "", "", "", "");
        }

        String gap_seq = "";
        String gap_ss = "";
        String seq1 = motif1.getAlignSeq();
        String seq2 = motif2.getAlignSeq();
        String ss1 = motif1.getAlignSs();
        String ss2 = motif2.getAlignSs();
        int start = motif1.getStart();
        int end = motif2.getEnd();

        // motif2 is contained in motif1
        if (motif2.getEnd() < motif1.getEnd()) {
            //return ($motif1->start, $motif1->end, $seq1, $ss1, "", "", "","");
            return new Merge_Motif_Output(-1, -1, "", "", "", "", "", "");
        }

        HashMap<Integer, Integer> map1 = motif1.getAlignMap();
        HashMap<Integer, Integer> map2 = motif2.getAlignMap();

        if (motif2.getStart() <= motif1.getEnd()) {
            int olap_start = motif2.getStart();
            int olap_end = motif1.getEnd();
            //assign the overlap region to motif1, remove it from motif2
            if (olap_own == 1) {
                HashMap<Integer, Integer> pt = pair_table(ss2);
                for (int i = 0; i < seq2.length(); i++) {
                    if (!map2.containsKey(i)) {
                        continue;
                    }
                    if (map2.get(i) > olap_end) {
                        break;
                    }
                    //replace with string index i to .
                    //perl using 4th version substr, str, offset, length, replacement
                    //index same as java starting from 0
                    seq2 = seq2.substring(0, i - 1) + "." + seq2.substring(i + 1);
                    ss2 = ss2.substring(0, i - 1) + "." + ss2.substring(i + 1);
                    if (!pt.containsKey(i)) {
                        continue;
                    }
                    if (pt.get(i) >= 0) {
                        ss2 = ss2.substring(0, pt.get(i) - 1) + "." + ss2.substring(pt.get(i) + 1);
                    }
                }
            } else {
                HashMap<Integer, Integer> pt = pair_table(ss1);
                for (int i = seq1.length() - 1; i >= 0; i--) {
                    if (!map1.containsKey(i)) {
                        continue;
                    }
                    if (map1.get(i) < olap_start) {
                        break;
                    }
                    seq1 = seq1.substring(0, i - 1) + "." + seq1.substring(i + 1);
                    ss1 = ss1.substring(0, i - 1) + "." + ss1.substring(i + 1);
                    if (pt.get(i) >= 0) {
                        ss1 = ss1.substring(0, pt.get(i) - 1) + "." + ss1.substring(pt.get(i) + 1);
                    }
                }
            }
            gap_seq = "";
            gap_ss = "";
        } //no overlap
        else {
            if (motif2.getStart() == motif1.getEnd() + 1) {
                gap_seq = "";
                gap_ss = "";
            } else {
                //substring (begin_index,end_index)
                //length end_index-begin_index, so end_index byte is excluded
                gap_seq = whole_seq.substring(motif1.getEnd(), motif2.getStart());

                if (gap_seq.length() > comb_max_gap) {
                    if (motif1.getScore() > motif2.getScore()) {
                        start = motif1.getStart();
                        end = motif1.getEnd();
                        gap_seq = "";
                        gap_ss = "";
                        seq2 = make_string(seq2.length(), ".");
                        ss2 = make_string(seq2.length(), "."); // check this is a bug? using ss2 instead?
                    } else {
                        start = motif2.getStart();
                        end = motif2.getEnd();
                        gap_seq = "";
                        gap_ss = "";
                        seq1 = make_string(seq1.length(), ".");
                        ss1 = make_string(seq1.length(), "."); // check this is a bug? using ss1 instead?
                    }
                }
                gap_ss = make_string(gap_seq.length(), ".");

            }

        }
        gap_seq = gap_seq.replaceAll("(?i)t", "U");  // sed ~s /[tT]/U/g;  case-insentive replace

        return new Merge_Motif_Output(start, end, seq1, ss1, gap_seq, gap_ss, seq2, ss2);
    }

    public static Alignment merge_alignment(
              Alignment alignment1,
              Alignment alignment2,
              HashMap<String, Seq> seqs,
              String out) throws MyException {

        HashMap<String, AlignSeq> align1 = alignment1.getSeqs();
        HashMap<String, AlignSeq> align2 = alignment2.getSeqs();
        String ss_cons1 = alignment1.getSsCons();
        String ss_cons2 = alignment2.getSsCons();
        String rf1 = alignment1.getRf();
        String rf2 = alignment2.getRf();

        HashMap<String, Merge_Motif_Output> motif_overlap = new HashMap<>();
        int max_gap_len = 0;
        int max_m1_len = 0;
        int max_m2_len = 0;
        double avg_gap_len = 0;
        int start;
        int end;
        String seq1;
        String ss1;
        String gap_seq;
        String gap_ss;
        String seq2;
        String ss2;
        HashMap<String, AlignSeq> merged_motif = new HashMap<>();
        String merged_ss_cons;
        String merged_rf;
        int olap_own = resolve_overlap(alignment1, alignment2);

        HashMap<String, Integer> ids = new HashMap<>();
        int max_id = 0;

        for (Map.Entry<String, AlignSeq> entry : align1.entrySet()) {
            align_seq_map(entry.getValue()); //entry AlignSeq's two maps updated
            ids.put(entry.getKey(), entry.getValue().getId());
            if (ids.get(entry.getKey()) > max_id) {
                max_id = ids.get(entry.getKey());
            }
        }

        for (Map.Entry<String, AlignSeq> entry : align2.entrySet()) {
            align_seq_map(entry.getValue()); //entry AlignSeq's two maps updated
            if (ids.containsKey(entry.getKey())) {
                continue;
            }
            max_id++;
            ids.put(entry.getKey(), max_id);
        }

        for (Map.Entry<String, Integer> entry : ids.entrySet()) {
            String acc;
            double weight;
            if (align2.containsKey(entry.getKey())) {
                acc = align2.get(entry.getKey()).getAcc();
            } else {
                acc = align1.get(entry.getKey()).getAcc();
            }
            if (!seqs.containsKey(entry.getKey())) {
                throw new MyException(entry.getKey() + " does not exist in input sequences ");
            }
            if (!align2.containsKey(entry.getKey())) {
                weight = align1.get(entry.getKey()).getWeight();
                Merge_Motif_Output mm = merge_motif(
                          align1.get(entry.getKey()),
                          null,
                          seqs.get(entry.getKey()).getSeq(),
                          olap_own
                );
                start = mm.getStart();
                end = mm.getEnd();
                seq1 = mm.getSeq1();
                ss1 = mm.getSs1();
                gap_seq = mm.getGapSeq();
                gap_ss = mm.getGapSs();
                seq2 = mm.getSeq2();
                ss2 = mm.getSs2();
            } else if (!align1.containsKey(entry.getKey())) {
                weight = align2.get(entry.getKey()).getWeight();
                Merge_Motif_Output mm = merge_motif(
                          null,
                          align2.get(entry.getKey()),
                          seqs.get(entry.getKey()).getSeq(),
                          olap_own
                );
                start = mm.getStart();
                end = mm.getEnd();
                seq1 = mm.getSeq1();
                ss1 = mm.getSs1();
                gap_seq = mm.getGapSeq();
                gap_ss = mm.getGapSs();
                seq2 = mm.getSeq2();
                ss2 = mm.getSs2();
            } else {
                weight = align2.get(entry.getKey()).getWeight()
                          + align1.get(entry.getKey()).getWeight();
                Merge_Motif_Output mm = merge_motif(
                          align1.get(entry.getKey()),
                          align2.get(entry.getKey()),
                          seqs.get(entry.getKey()).getSeq(),
                          olap_own
                );
                start = mm.getStart();
                end = mm.getEnd();
                seq1 = mm.getSeq1();
                ss1 = mm.getSs1();
                gap_seq = mm.getGapSeq();
                gap_ss = mm.getGapSs();
                seq2 = mm.getSeq2();
                ss2 = mm.getSs2();
            }

            if (start >= 0 && end >= 0) {
                Merge_Motif_Output mm = new Merge_Motif_Output(
                          start, end, seq1, ss1, gap_seq, gap_ss, seq2, ss2
                );
                mm.setWeight(weight / 2.0d);
                motif_overlap.put(entry.getKey(), mm);
                if (max_gap_len < gap_seq.length()) {
                    max_gap_len = gap_seq.length();
                }
                if (max_gap_len < seq1.length()) {
                    max_gap_len = seq1.length();
                }
                if (max_gap_len < seq2.length()) {
                    max_gap_len = seq2.length();
                }
            }
        } // loop ids is done

        if (motif_overlap.size() < 2) {
            System.out.println("Merge_alignment returns null");
            return null;
        }

        //now out file should be using gap file extension
        String out_file = out;
        if (out_file == null || out_file.isEmpty() || out_file.equals("")) {
            out_file = "out.gap";
        }
        Path p_out = Paths.get(out_file);
        int gap_count = 0;
        try {
            Files.createFile(p_out);
            //now loop motif_overlap
            for (Map.Entry<String, Merge_Motif_Output> entry : motif_overlap.entrySet()) {
                String gap = motif_overlap.get(entry.getKey()).getGapSeq();
                if (gap.length() > 0) {
                    Files.write(p_out,
                              (">" + entry.getKey() + System.lineSeparator()).getBytes(),
                              StandardOpenOption.APPEND);
                    Files.write(p_out,
                              (gap + System.lineSeparator()).getBytes(),
                              StandardOpenOption.APPEND);
                    gap_count++;
                }
            }
        } catch (FileAlreadyExistsException ignored) {
            System.err.println("File already exists, use existing.");
        } catch (IOException e) {
            System.err.println("Can't open file " + out_file + ": " + e);
        }

        if (gap_count > 1) {
            if (jo.optBoolean("skipClustalw")) {
                if (output_more_like_old_cmfinder_pl) {
                    System.out.println("Can't exec \"clustalw\": No such file or directory at " + bin_path + "/merge_motif.pl line 290.");
                    System.out.println("FATAL: Alignment file " + out + ".gap.aln could not be opened for reading");
                    System.out.println("Illegal division by zero at " + bin_path + "/io.pl line 353.");
                } else {
                    System.out.println("Aborting merge_motif since -skipClustalw");
                }
                return null;
            } else {
                String[] cmd1 = {findFile(bin_path, "clustalw"),
                    "-infile=" + out_file,
                    "-outfile=" + out_file + ".aln"};
                runCmd(cmd1, 1800); //giving runCmd result

                String[] cmd2 = {findFile(bin_path, "sreformat"),
                    "stockholm",
                    out_file + ".aln",
                    ">",
                    out_file + ".align"};
                runCmd(cmd2, 1800);//giving runCmd result
                Alignment gap_sto = read_stockholm(out_file + ".align");
                HashMap<String, AlignSeq> gap_align = gap_sto.getSeqs();

                for (Map.Entry<String, Merge_Motif_Output> entry : motif_overlap.entrySet()) {
                    if (gap_align.containsKey(entry.getKey())) {
                        gap_seq = gap_align.get(entry.getKey()).getAlignSeq();
                        entry.getValue().setGapSeq(gap_seq);
                        entry.getValue().setGapSs("");
                        if (max_gap_len < gap_seq.length()) {
                            max_gap_len = gap_seq.length();
                        }
                    }
                }
                //delete out.gap out_file
                deleteFile(out_file);
            }
        }

        merged_ss_cons = ss_cons1 + pad_string("", max_gap_len, ".", 1) + ss_cons2;
        merged_rf = rf1 + pad_string("", max_gap_len, ".", 1) + rf2;

        //  my ($align_seq,$align_ss,$score1,$score2,$score,$weight,$desc);
        String align_seq;
        String align_ss;
        float score1 = 0;
        float score2 = 0;
        float score;
        double weight;
        String desc;

        for (Map.Entry<String, Merge_Motif_Output> entry : motif_overlap.entrySet()) {
            gap_seq = pad_string(entry.getValue().getGapSeq(), max_gap_len, ".", 1);
            gap_ss = pad_string(entry.getValue().getGapSs(), max_gap_len, ".", 1);
            align_seq = pad_string(entry.getValue().getSeq1(), max_m1_len, ".", 1)
                      + gap_seq + pad_string(entry.getValue().getSeq2(), max_m2_len, ".", 1);
            align_ss = pad_string(entry.getValue().getSs1(), max_m1_len, ".", 1)
                      + gap_ss + pad_string(entry.getValue().getSs2(), max_m2_len, ".", 1);
            if (align1.containsKey(entry.getKey())) {
                score1 = align1.get(entry.getKey()).getScore();
            }
            if (align2.containsKey(entry.getKey())) {
                score2 = align2.get(entry.getKey()).getScore();
            }
            start = entry.getValue().getStart();
            end = entry.getValue().getEnd();
            score = score1 + score2 - remove_gap(gap_seq).length();
            weight = entry.getValue().getWeight();
            desc = String.format("%3d", start)
                      + ".." + String.format("%3d", end) + "\t"
                      + String.format("%.3f", score);

            merged_motif.put(entry.getKey(),
                      new AlignSeq(entry.getKey(), ids.get(entry.getKey()),
                                start, end, desc, score, weight, align_seq, align_ss));
        }

        //final return value
        return new Alignment(merged_motif, alignment1.getFlags(), merged_ss_cons, merged_rf);
    }

    public static void CombMotif(String cand_weight_option, String seq_file, ArrayList<String> motifFilesRef)
              throws MyException {

        ArrayList<String> align_files = motifFilesRef;
        ArrayList<String> all_files = align_files;
        HashMap<String, HashMap<String, String>> all_stats = new HashMap<>();

        for (String ff : align_files) {
            Alignment align = read_stockholm(ff);
            if (align.getWeight() >= comb_min_num) {
                HashMap<String, String> stat = RunSummarize(ff);
                all_stats.put(ff, stat);
                alignments.put(ff, align);
            }
        }
        //try merging all pairs of motifs, and see how they fit together
        HashMap<String, MergeMotif> merge_motif = new HashMap<>();
        for (String f1 : align_files) {
            for (String f2 : align_files) {
                // java This method returns 0 if two Strings are equal or if both are null, 
                // a negative number if the first String comes before the argument, and 
                // a number greater than zero if the first String comes after the argument String
                if (f1.compareTo(saveTimerFlag) < 0) {
                    merge_motif = try_merge(f1, f2, merge_motif);
                }
            }
        }

        HashMap<String, MergeMotif> processed = new HashMap<>();
        HashMap<String, String> merged_files = new HashMap<>();
        //add crete a reverse sort merge_motif_reverse

        HashMap<String, MergeMotif> merge_motif_reverse_sorted
                  = merge_motif.entrySet().stream()
                            .sorted(Map.Entry.comparingByValue(
                                      Comparator.comparingDouble(MergeMotif::getWeight)
                                                .reversed()))
                            .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                                      (oldValue, newValue) -> oldValue, LinkedHashMap::new));

        while (merge_motif.size() > 0) {
            if (verbose) {
                System.out.println("entering for my \\$id, with list:");
                //print key value based on object's weight reverse sorted
                merge_motif_reverse_sorted.entrySet().stream()
                          .forEachOrdered(e -> System.out.println(e.getKey()));
            }
            for (Map.Entry<String, MergeMotif> entry : merge_motif_reverse_sorted.entrySet()) {
                //ake motifs whose combination has the biggest weight first
                if (verbose) {
                    System.out.println("while keys merge_motif : id=" + entry.getKey());
                }
                //MergeMotif m = merge_motif.get(entry.getKey());
                MergeMotif m = entry.getValue();
                merge_motif.remove(entry.getKey());

                if (processed.containsKey(entry.getKey())) {
                    if (verbose) {
                        System.out.println("while keys merge_motif NEXT : exists \\$processed{"
                                  + entry.getKey() + "}");
                    }
                    continue;
                }
                processed.put(entry.getKey(), m);
                String f1 = m.getMotif1();
                String f2 = m.getMotif2();

                if (m.getWeight() <= 0) {
                    if (verbose) {
                        System.out.println("while keys merge_motif NEXT : weight <= 0");
                    }
                }

                if (m.getWeight() > 0) { //has to be favorable
                    if (merged_files.containsKey(f1) || merged_files.containsKey(f2)) {
                        //# ??  apparently each motif can only appear in one merging?
                        if (verbose) {
                            System.out.println("while keys merge_motif NEXT : exists \\$merged_files{"
                                      + f1 + "}=" + (merged_files.containsKey(f1) ? "true" : "false")
                                      + " || exists \\$merged_files{" + f2 + "}="
                                      + (merged_files.containsKey(f2) ? "true" : "false"));
                        }
                        continue;
                    }
                    if (m.getGap() > comb_max_gap) {
                        double m_gap = m.getGap();
                        if (verbose) {
                            System.out.println("while keys merge_motif NEXT : \\$m->gap > \\$comb_max_gap :"
                                      + m_gap + ">" + comb_max_gap);
                        }
                        continue;
                    }

                    String f = entry.getKey();
                    int found = 0;
                    for (String tmp : all_files) {
                        if (tmp.equals(f)) {
                            found = 1;
                            break;
                        }
                    }

                    if (found != 0) {
                        if (verbose) {
                            System.out.println("while keys merge_motif NEXT : not found.  all_files = "
                                      + String.join(" ", all_files));
                        }
                        continue;  //# file has already been made
                    }

                    System.out.println("( near merge_motif: "
                              + entry.getKey() + ", " + m.getNum_seq() + ", \t" + f1 + "\t" + f2 + "\t,"
                              + m.getWeight() + ",t" + m.getGap() + ",\t" + m.getOverlap());

                    // this used to be a call to the merge_motif.pl script
                    if (output_more_like_old_cmfinder_pl) {
                        if (verbose) {
                            System.out.println(bin_path + "/merge_motif.pl " + SEQ
                                      + " " + f1 + " " + f2 + " " + f + ".temp");
                        }
                    } else {
                        if (verbose) {
                            System.out.println("call to merge_motif "
                                      + f1 + " " + f2 + " " + f + ".temp");
                        }

                        if (!alignments.containsKey(f1) || !alignments.containsKey(f2)) {
                            throw new MyException("internal error");
                        }

                        String f_temp = f + ".temp";
                        Alignment new_alignment
                                  = merge_alignment(alignments.get(f1), alignments.get(f2), unaligned_seqs, f_temp);
                        if (new_alignment != null) {
                            if (!skipClustalw) {
                                throw new MyException("unexpected");
                            }
                            if (output_more_like_old_cmfinder_pl) {
                                if (verbose) {
                                    String cmfile = seq_file;
                                    cmfile = cmfile.replace(".motif.", ".cm.");  //$cmfile =~ s/[.]motif[.]/.cm./g;
                                    System.out.println(bin_path + "/cmfinder      -o "
                                              + f + " -a " + f_temp + " " + seq_file + " " + cmfile);
                                    System.out.println("FATAL: Alignment file "
                                              + f_temp + " could not be opened for reading");
                                    System.out.println("while keys merge_motif NEXT : !-s " + f);
                                }
                            } else {
                                if (verbose) {
                                    System.out.println("while keys merge_motif NEXT : -skipClustalw");
                                }
                            }
                            continue;
                        }
                        write_stockholm(new_alignment, f_temp);

                        //2019-04-24 continue
                        String[] a = {bin_path + "/" + cmfinderBaseExe,
                            saveTimerFlag
                    // todo @cmfinder_inf11FlagsList
                    };
                        if (!RunCmfinder(a, f)) {
                            //couldn't produce acceptable output
                            //perl next = java continue?
                            continue;
                        }
                    }
                }
            }
        }
    }

    public static HashMap<String, MergeMotif>
              try_merge(String f1, String f2, HashMap<String, MergeMotif> merge_motif_ref)
              throws MyException {
        if (merge_motif_ref == null) { //perl !defined not isEmpty()
            throw new MyException("merge_motif_ref Hashmap is empty");
        }

        String index = String.join(".", my_strcmp(f1, f2));
        System.out.println("try_merge " + f1 + " " + f2);
        if (merge_motif_ref.containsKey(index)) {
            if (verbose) {
                System.out.println();
            }
            return null;
        }

        double num_overlap = 0;
        double start1 = 0;
        double start2 = 0;
        double start1_score = 0;
        double start2_score = 0;
        int gap1 = 0;
        int gap2 = 0;
        int overlap1 = 0;
        int overlap2 = 0;

        if (!alignments.containsKey(f1) || !alignments.containsKey(f2)) {
            if (verbose) {
                System.out.println("alignments don't containt key for " + f1 + " or " + f2);
            }
            return null;
        }

        HashMap<String, AlignSeq> align1 = alignments.get(f1).getSeqs();
        HashMap<String, AlignSeq> align2 = alignments.get(f2).getSeqs();

        for (HashMap.Entry<String, AlignSeq> entry : align1.entrySet()) {
            if (!align2.containsKey(entry.getKey())) {
                continue;
            }
            //below, we know that hit id $id is common to both alignments $f1 and $f2
            //figure out if the hits for id=$id in the two motifs overlap
            AlignSeq motif1 = align1.get(entry.getKey());
            AlignSeq motif2 = align2.get(entry.getKey());

            //we only look for motifs on the forward strand
            if ((motif1.getStart() > motif1.getEnd())
                      || (motif2.getStart() > motif2.getEnd())) {
                throw new MyException("motif start is bigger than end unexpected");
            }

            int len1 = motif1.getEnd() - motif1.getStart();
            int len2 = motif2.getEnd() - motif2.getStart();

            if (motif1.getStart() > motif2.getStart()) {
                if (motif1.getEnd() < motif2.getEnd()) {
                    //overlap
                    //motif1 contained within motif2
                    num_overlap += motif1.getWeight() * motif2.getWeight();
                } else {
                    int overlap = motif2.getEnd() - motif1.getStart();
                    if (overlap > -comb_max_gap) {
                        if (((len1 - overlap < 25) || (len2 - overlap < 25))
                                  && ((overlap > 0.9 * len1) || (overlap > 0.9 * len2))) {
                            //the motif instances almost entirely overlap (have up leave at most 25 nucs 
                            //and 10% of each instance un-overlapped
                            num_overlap += motif1.getWeight() * motif2.getWeight();
                        } else {
                            //he motif instances are near to one another or somewhat overlapping, let's 
                            //record the weight for this orientation (motif2 instance coords > motif1 instance coords)
                            start2 += (motif1.getWeight() + motif2.getWeight()) / 2;
                            start2_score += motif1.getScore() * motif1.getWeight() + motif2.getScore() * motif2.getWeight();
                            if (overlap > 0) {
                                //the motif instances actually overlap
                                overlap2 += overlap * motif1.getWeight() * motif2.getWeight();
                            } else if (-overlap < comb_max_gap) {//I think this is always true because of the test some lines up
                                gap2 += -overlap * motif1.getWeight() * motif2.getWeight();
                            }
                        }
                    }
                }
            } else { //this 'else' case is the mirror image of the 'if' case
                if (motif1.getEnd() < motif2.getEnd()) {
                    int overlap = motif1.getEnd() - motif2.getStart();
                    if (overlap > -comb_max_gap) {
                        if (((len1 - overlap < 25) || (len2 - overlap < 25))
                                  && ((overlap > 0.9 * len1) || overlap > 0.9 * len2)) {
                            num_overlap += motif1.getWeight() * motif2.getWeight();
                        } else {
                            start1 += (motif1.getWeight() + motif2.getWeight()) / 2;
                            start1_score += motif1.getScore() * motif1.getWeight() + motif2.getScore() * motif2.getWeight();
                            if (overlap > 0) {
                                overlap1 += overlap * motif1.getWeight() * motif2.getWeight();
                            } else if (-overlap1 < comb_max_gap) {
                                gap1 += -overlap * motif1.getWeight() * motif2.getWeight();
                            }
                        }
                    }
                } else {
                    //overlap
                    // motif2 contained within motif1
                    num_overlap += motif1.getWeight() * motif2.getWeight();
                }
            }
        }

        if (verbose) {
            System.out.println("start1=" + start1 + " start2=" + start2 + " num_overlap=" + num_overlap);
        }
        if (num_overlap > start1 && num_overlap > start2) {
            if (verbose) {
                System.out.println("num_overlap>start1 && num_overlap>start2 :"
                          + num_overlap + ">" + start1 + " && " + num_overlap + ">" + start2);
            }
            return null;
        }
        if (start1_score > start2_score) {
            //# is motif1 more often to the left of motif2 or vice versa?  Where's the bias?
            if (start1 < comb_min_overlap) {
                //# not enough sequences favor this orientation in absolute terms (compared to the constant 
                //$comb_min_overlap) for it to be worthwhile  
                if (verbose) {
                    System.out.println("start1<comb_min_overlap :" + start1 + "<" + comb_min_overlap);
                }
                return null;
            }

            //	#motif1 is before motif2
            merge_motif_ref.put(index,
                      new MergeMotif(
                                f1,
                                f2,
                                start1,
                                start1_score,
                                gap1 / start1,
                                overlap1 / start1,
                                start1_score - gap1 / 2 - overlap1));
            System.out.println("accept id=" + index + ":"
                      + f1 + "-" + f2 + ":" + start1_score + " " + gap1 / 2 + " " + overlap1 + " " + start1
            );
        } else { //# this 'else' case mirrors the 'if' case
            if (start2 < comb_min_overlap) {
                if (verbose) {
                    System.out.println("start2 < comb_min_overlap :" + start2 + "<" + comb_min_overlap);
                }
                return null;
            }
            //#motif2 is before motif1
            merge_motif_ref.put(index,
                      new MergeMotif(
                                f2,
                                f1,
                                start2,
                                start2_score,
                                gap2 / start2,
                                overlap2 / start2,
                                start2_score - gap2 / 2 - overlap2)
            );
            System.out.println("taccept id=" + index + ":"
                      + f2 + "-" + f1 + ":" + start2_score + " " + gap2 / 2 + " " + overlap2 + " " + start2);
        }

        return merge_motif_ref;
    }

    public static String[] my_strcmp(String s1, String s2) {
        String[] t1 = s1.split("\\.");
        String[] t2 = s2.split("\\.");
        String prefix = "";
        for (int i = 0; i < t1.length && i < t2.length; i++) {
            if (!t1[i].equals(t2[i])) {
                break;
            }
            if (prefix.equals("")) {
                prefix = t1[i];
            } else {
                prefix = prefix + "." + t1[i];
            }
        }
        String suffix1 = String.join(".", t1);
        String suffix2 = String.join(".", t2);
        return new String[]{prefix, suffix1, suffix2};
    }

    public static boolean RunCmfinder(String[] cmd, String output_motif_file) throws MyException {
        if (output_motif_file == null || output_motif_file.isEmpty()) {
            throw new MyException("must pass output_motif_file");
        }
        if (output_more_like_old_cmfinder_pl) {
            System.out.println(String.join(" ", cmd));
        } else {
            System.out.println("Running:" + String.join(" ", cmd));
        }

        cmdOut co = runCmd(cmd);  //30 minutes default
        boolean status = co.getextVal();
        int exitCode = co.getExitCode();

        if (status) {
            //made an alignment
            return true;
        } else {
            if (exitCode >= 2 && exitCode <= 8) {
                //couldn't produce acceptable output
                return false;
            } else {
                if (useOldCmfinder) {
                    System.out.println("cmfinder returned error, but I'm ignoring it because I think it's benign");
                    deleteFile(output_motif_file);
                    return true;
                } else {
                    throw new MyException("problem running cmd:" + String.join(" ", cmd));
                }
            }
        }
    }

    // note we assume the file is in bin_path folder
    public static HashMap<String, String> RunSummarize(String file_name) throws MyException {
        String f = findFile(bin_path, file_name);
        if (false) {  // prel if(0)
            String[] cmd = {findFile(bin_path, "summarize"), f};
            System.out.println("Running " + String.join(" ", cmd));
            cmdOut exec1 = runCmd(cmd, 300);
            ArrayList<String> summary = exec1.getOut(); //ArrayList to String
            if (!exec1.getextVal()) {
                throw new MyException(String.join(" ", cmd) + " failed.");
            }
        }
        String cmfinder = findFile(bin_path, "cmfinder04");
        String[] cmd2 = {cmfinder, saveTimerFlag, "--summarize", summarizeFlagsStr, f};
        System.out.println("Running " + String.join(" ", cmd2));
        cmdOut exec2 = runCmd(cmd2, 300);
        ArrayList<String> summary2 = exec2.getOut();
        if (!exec2.getextVal()) {
            throw new MyException(String.join(" ", cmd2) + " failed.");
        }
        //he matcher.group() function expects to take a single integer argument: The capturing group index, starting from 1. 
        //The index 0 is special, which means "the entire match". 
        Pattern p = Pattern.compile("\\s*(\\S+)=(\\S+)\\s*");
        HashMap<String, String> result = new HashMap<>();
        summary2.forEach(e -> {
            Matcher m = p.matcher(e);
            if (m.matches()) {
                result.put(m.group(1), m.group(2));
            }
        }
        );
        return result;
    }

    public void print_version() {
        System.out.println("CMFINDER_PACKAGE_VERSION=" + CMFINDER_PACKAGE_VERSION);
    }

    public static void print_help() {
        String help_content = "CMFINDER [options] SEQ\n"
                  + "Options:\n"
                  + "    -c <number>      \n"
                  + "     The maximum number of candidates in each sequence. Default 40. No bigger than 100.\n"
                  + "    -m <number>      \n"
                  + "     The minimum length of candidates. Default 30\n"
                  + "    -M <number>      \n"
                  + "     The maximum length of candidates. Default 100\n"
                  + "    -f <number>      \n"
                  + "     The fraction of the sequences expected to contain the motif. Default 0.80\n"
                  + "    -s1 <number>     \n"
                  + "     The max number of output single stem-loop motifs\n"
                  + "    -s2 <number>    \n"
                  + "     The max number of output double stem-loop motifs    \n"
                  + "    -minspan1 <number>\n"
                  + "     minimum span of a candidate sub-sequence in the heuristics to come up with an initial alignment for single-hairpin (h1) motifs\n"
                  + "    -maxspan1 <number>\n"
                  + "     like -minspan1, but maximum\n"
                  + "    -minspan2 <number>\n"
                  + "     like -minspan1, but for double-hairpin (h2) motifs\n"
                  + "    -maxspan2 <number>\n"
                  + "     like -minspan2, but maximal\n"
                  + "    -combine         \n"
                  + "     Combine the output motifs. Default False\n"
                  + "    -motifList <file> \n"
                  + "     Produce a list of motifs generated, one motif per line.\n"
                  + "    -o <number>\n"
                  + "     Minimum overlap for combining motifs\n"
                  + "    -n <number>      \n"
                  + "     Minimum number of sequences (weighted) for combining motifs\n"
                  + "    -emSeq <file>\n"
                  + "     Use the sequences in this fasta file for the expectation maximization step (i.e., the C executable cmfinder), but not for the earlier steps related to finding candidate motifs.  The reason for this distinction is that it is somewhat easier to add weighting to the cmfinder program, than the various canda, candf, cands and align programs.\n"
                  + "    -likeold         \n"
                  + "     Behave as much as possible like the old CMfinder, e.g., passing --enone, --p56 and --degen-rand to cmfinder_inf11.  It's not possible to produce identical results to CMfinder 0.3, but these flags make it more similar.\n"
                  + "    -fragmentary\n"
                  + "     Pass --fragmentary for cmfinder\n"
                  + "    -amaa            \n"
                  + "     Pass --amaa to cmfinder (align max align accuracy)\n"
                  + "    -useOldCmfinder  \n"
                  + "     Run the old cmfinder executable, mainly to test whether we get different results because of this perl script, or the cmfinder_inf11 executable.\n"
                  + "    -skipClustalw    \n"
                  + "     Do not run clustalw, like older installations lacking this program.\n"
                  + "    -justGetCmfinderCommand    \n"
                  + "     Print the command to run for the cmfinder executable, with appropriate partial flags.  This can be used to realign an existing .sto file, for example.\n"
                  + "    -copyCmfinderRunsFromLog <log-file> \n"
                  + "     For debugging.  Reads a log file that contains cmfinder commands, and re-runs them with new CMfinder.\n"
                  + "    -commaSepEmFlags x<flags>\n"
                  + "     List of flags and arguments to pass to the EM-step cmfinder exe.  There's an 'x' at the beginning of the flags, so that perl doesn't interpret the flags as flags for it.  It's comma-separated where on the command line it would be space separated.  I think commas are safe, and mean that I don't have to worry about quoting stuff.  e.g., -commaSepEmFlags x--fragmentary,--filter-non-frag,--filter-non-frag-pad,10 would pass this to the cmfinder program: \"--fragmentary --filter-non-frag --filter-non-frag-pad 10\", i.e., just replace commas with spaces.\n"
                  + "    -commaSepSummarizeFlags x<flags>\n"
                  + "     Flags to pass to the --summarize command.  Same syntax as for --commaSepEmFlags\n"
                  + "    -commaSepCandfFlags x<flags>\n"
                  + "     Flags to pass to the candf command.  Same syntax as for --commaSepEmFlags\n"
                  + "    -minCandScoreInFinal <number>    \n"
                  + "     Pass --min-cand-score-in-final <number> to cmfinder.  WARNING: there's a difference between using this flag (where even intermediate motifs will avoid these hits) and taking the low-scoring instances out of the final alignments (which might be combinations of motifs in which the sequence would have been lower-scoring).\n"
                  + "    -filterNonFrag\n"
                  + "     Pass --filter-non-frag to cmfinder\n"
                  + "    -columnOnlyBasePairProbs\n"
                  + "     Pass --column-only-base-pair-probs to cmfinder\n"
                  + "    -saveTimer <file>\n"
                  + "     create tab-delimited <file> containing timing stats on various sub-processes of this script.  the first tab-delimited field is the description of the sub-process, the second field is the total CPU time (user+sys) and the third field is the wall-clock time.  Sub-processes can occur in multiple lines if they are run multiple timers, so the caller should add them.  Due to my laziness, the time of the clustalw program (if used) is not counted.\n"
                  + "    -cpu <num>\n"
                  + "     use <num> CPUs for functionality that can use multi-CPUs (currently only the internal cmsearch commands in cmfinder04)\n"
                  + "    -allCpus\n"
                  + "     equivalent to -cpu X , where X is the number of available processors.\n"
                  + "    -candsParallel\n"
                  + "     run the two cands jobs in parallel, even if -cpu 1\n"
                  + "    -outFileSuffix <string>\n"
                  + "     add <string> to the output file names.  this is useful if you want to run the script in multiple ways in the same directory.\n"
                  + "    -h               \n"
                  + "     Show this list\n"
                  + "    -version\n"
                  + "     Show package version\n";
        System.err.println(help_content);
    }

    public static void main(String[] args) {

        //test case insensitive replaceAll
        //        String target = "FOOfBar";
        //        target = target.replaceAll("(?i)f", "U");
        //        System.out.println(target);
        //test findfiles
        //        Path p = Paths.get(bin_path + "/clustalw");
        //        System.out.println(p.toAbsolutePath().toString());
        //        findFiles(bin_path, "clustalw").forEach(System.out::println);
        //        System.out.println("find it:" + findFile(bin_path, "clustalw"));
        //test my_strcmp
        //System.out.println(Arrays.toString(my_strcmp("abc.def.ghi", "abc.def.ghi")));
        SEQ = "";
    }
}
