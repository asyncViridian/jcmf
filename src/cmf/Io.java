package cmf;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.*;

/**
 * replace bin/io.pl
 */
public class Io {

    /**
     * Remove gap characters '.' and '-' from sequence string
     *
     * @param s sequence
     * @return de-gapped sequence
     */
    public static String remove_gap(String s) {
        String strNew = s.replace(".", "").replace("-", "");
        return strNew;
    }

    /**
     * Create n copies of ch
     *
     * @param n number of times to repeat ch
     * @param ch string to repeat
     * @return ch repeated n times
     */
    public static String make_string(int n, String ch) {
        //String strNew = new String(new char[n]).replace("\0", s); // before java 8
        String strNew = String.join("", Collections.nCopies(n, ch));  // Java 8, Java 11 repeat()
        return strNew;
    }

    /**
     * Adds up to (n - seq.length) copies of ch to either the left or the right
     * side of seq, such that the resulting sequence is n units long
     *
     * @param seq sequence to pad onto
     * @param n length to make the result sequence reach
     * @param ch sequence to pad with
     * @param dir direction to pad in: if 0, left; if 1, right
     * @return padded sequence
     */
    public static String pad_string(String seq, int n, String ch, int dir) {
        String strNew = seq;
        int l = n - seq.length();
        if (l > 0) { //pad right
            if (dir == 1) {
                strNew = seq + make_string(l, ch);
            } // pad left
            else if (dir == 0) {
                strNew = make_string(l, ch) + seq;
            }
        }
        return strNew;
    }

    /**
     * Reads the contents of the given file(name) into a mapping from sequence
     * name to sequence struct
     *
     * @param file_name filename to read
     * @return map from sequence nickname to struct representing sequence
     */
    public static HashMap<String, Seq> read_fasta(String file_name) {
        Path fasta = Paths.get(file_name);
        HashMap<String, Seq> seqs = new HashMap<>();
        int i = 0;
        //pattern for > line
        Pattern p = Pattern.compile("^>(\\S+)\\s+(\\S+.*)");
        Matcher m;
        String acc = null;
        String desc = null;
        StringBuilder seq = new StringBuilder();
        try {
            List<String> lines = Files.readAllLines(fasta);
            //acc = lines.get(0).substring(1);
            //desc = lines.get(0).substring(1);
            for (String line : lines) {
                if (line.length() > 0 && line.charAt(0) == '>') {
                    // begin next entry
                    if (i > 0) {
                        seqs.put(acc, new Seq(acc, i, desc, seq.toString()));
                    }
                    //try to find > line match the pattern with space, then break acc and desc
                    if ((m = p.matcher(line)).find()) {  //find() vs matches, later is whole string
                        acc = m.group(1);
                        desc = m.group(2);
                        //System.out.println("id= " + (i + 1) + " break acc=" + acc + " desc=" + desc);
                    } else {
                        acc = line.substring(1);
                        desc = line.substring(1);
                    }
                    i++;
                    seq = new StringBuilder();
                } else {
                    // continue an entry
                    seq.append(line);
                }
            }
            //update map value
            seqs.put(acc, new Seq(acc, i, desc, seq.toString()));
        } catch (IOException e) {
            System.out.println("Can't open file: " + e);
        }

        return seqs;
    }

    /**
     * Output fasta to a file or screen based on fasta data input order
     *
     * @param file_name filename to read
     * @param seqs is a Map<String, Seq>
     */
    public static void write_fasta(String file_name, HashMap<String, Seq> seqs) {
        //sort using treemap
        TreeMap<Integer, Seq> tree_seqs = new TreeMap<>();
        //move into TreeMap
        seqs.entrySet().stream()
                .forEach(e -> tree_seqs.put(e.getValue().getId(), e.getValue())
                );
        // now try to write file or print to screen
        if (file_name != null && !file_name.isEmpty()) {
            Path fasta_out = Paths.get(file_name);
            try {
                Files.write(fasta_out,
                        (Iterable<String>) tree_seqs.entrySet().stream()
                                .map(e -> ">" + e.getValue().getAcc() + "\t"
                                + ((e.getValue().getDesc().equals(e.getValue().getAcc())) ? "" : e.getValue().getDesc())
                                + System.lineSeparator()
                                + e.getValue().getSeq()
                                )::iterator);
            } catch (IOException e) {
                System.out.println("Can't open file for writing " + e);
            }
        } else {
            //print to screen treemap
            tree_seqs.entrySet().stream().forEach(e -> System.out.println(
                    ">" + e.getValue().getAcc() + "\t"
                    + ((e.getValue().getDesc().equals(e.getValue().getAcc())) ? "" : e.getValue().getDesc()) + "\n"
                    + e.getValue().getSeq()
            //e.getValue().getSeqString()
            )
            );
        }
    }

    //STOCKHOLM methods
    /**
     * Output AlignSeq object
     *
     * @param a obj
     * @return a
     *
     */
    public static AlignSeq align_seq_map(AlignSeq a) {
        String seq = a.getAlignSeq();
        Integer j = a.getStart();
        Boolean dir = a.getStart() < a.getEnd();
        String[] letters = seq.split("");
        //System.out.println(Arrays.toString(letters));
        Pattern p = Pattern.compile("\\.|-");
        for (int i = 0; i < seq.length(); i++) {
            if (!p.matcher(letters[i]).find()) {
                a.getAlignMap().put(i, j);
                a.getRevMap().put(j, i);
                if (dir) {
                    j++;
                } else {
                    j--;
                }
            }
        }
        return a;
    }

    /**
     * Output HashMap
     *
     * @param ss_str from AlignSeq align_ss value
     * @return pt
     *
     */
    public static HashMap<Integer, Integer> pair_table(String ss_str) {
        HashMap<Integer, Integer> pt = new HashMap<>();
        Integer i = 0;
        Integer j = 0;
        if (ss_str != null && !ss_str.isEmpty()) {
            int len = ss_str.length();
            String[] ss = ss_str.split("");
            ArrayList<Integer> stack = new ArrayList<>();
            while (i < len) { //int compare Integer is fine
                while (ss[i].equals("<")) {
                    stack.set(j, i);
                    j++;
                    i++;
                }
                while (i < len && ss[i].equals(">")) {
                    j--;
                    pt.put(i, stack.get(j));
                    pt.put(stack.get(j), i);
                    i++;
                }
                if (i < len && !ss[i].equals("<") && !ss[i].equals(">")) {
                    i++;
                }
            }
        } else {
            try { // better not using easy runtime exception
                throw new MyException("!defined(ss_str), which is odd");
            } catch (MyException ex) {
                System.out.println(ex.getMessage());
            }
        }
        return pt;
    }

    /**
     * read Stockholm file output Alignment obj
     *
     * @param file_name filename to read
     * @return alignment
     */
    public static Alignment read_stockholm(String file_name) {

        String id = "";
        String acc;
        int weight;
        float score;
        String desc;
        Integer start;
        Integer end;

        String ss_cons = "";
        String rf = "";
        HashMap<String, AlignSeq> seqs = new HashMap<>();
        HashMap<String, Integer> flags = new HashMap<>();
        AlignSeq new_seq;

        Path RFAM = Paths.get(file_name);
        Matcher m;

        try {
            List<String> lines = Files.readAllLines(RFAM);
            int new_id = 0;
            int last_count = 0;
            String last_acc = "";
            for (String line : lines) {
                //main process start

                //1st line check
                if (Pattern.compile("^\\#\\s+").matcher(line).matches()
                        || Pattern.compile("^\\s*$").matcher(line).matches()
                        || Pattern.compile("^\\#=GF").matcher(line).matches()) {
                    last_acc = "";
                    last_count = 0;
                    continue;
                }

                //2nd line check
                if ((m = Pattern.compile("\\#=GS\\s+(\\S+)\\s+WT\\s+(\\S+)")
                        .matcher(line))
                        .find()) {
                    acc = m.group(1);
                    weight = Integer.parseInt(m.group(2)); //parseInt takes String
                    if (!last_acc.equals(acc)) {
                        last_acc = acc;
                        last_count = 0;
                        id = acc;
                    } else {
                        last_count++;
                        id = acc + last_count;
                    }
                    if (!seqs.containsKey(id)) {
                        new_seq = new AlignSeq(acc, weight, new_id);
                        flags.put("WGT", 1);
                        new_id++;
                        seqs.put(id, new_seq);
                    } else {
                        seqs.get(id).setWeight(weight);
                    }
                } else if ((m = Pattern.compile("\\#=GS\\s+(\\S+)\\s+DE\\s+(.+)$")
                        .matcher(line))
                        .find()) {
                    acc = m.group(1);
                    desc = m.group(2);
                    flags.put("DE", 1);
                    if (!last_acc.equals(acc)) {
                        last_acc = acc;
                        last_count = 0;
                        id = acc;
                    } else {
                        last_count++;
                        id = acc + last_count;
                    }
                    if (!seqs.containsKey(id)) {
                        new_seq = new AlignSeq(acc, 1, new_id, desc, "desc");  //weight int?
                        new_id++;
                        seqs.put(id, new_seq);
                    } else {
                        seqs.get(id).setDesc(desc);
                    }
                    if ((m = Pattern.compile("(\\d+)\\.\\.\\s*(\\d+)\\s*(\\S*)$")
                            .matcher(desc))
                            .find()) {
                        start = Integer.parseInt(m.group(1));  //  \d+ to Integer
                        end = Integer.parseInt(m.group(2));
                        flags.put("START", 1);
                        flags.put("END", 1);
                        seqs.get(id).setStart(start);
                        seqs.get(id).setEnd(end);
                        if (!m.group(3).equals("")) {  // extract (\S*) string
                            score = Float.valueOf(m.group(3));
                            flags.put("SC", 1);
                            seqs.get(id).setScore(score);
                        }
                    }
                } else if ((m = Pattern.compile("\\#=GC\\s+SS_cons\\s+(\\S+)")
                        .matcher(line))
                        .find()) {
                    flags.put("SS_cons", 1);
                    ss_cons = ss_cons + m.group(1);
                } else if ((m = Pattern.compile("\\#=GC\\s+RF\\s+(\\S+)")
                        .matcher(line))
                        .find()) {
                    flags.put("RF", 1);
                    rf = rf + m.group(1);
                } else if ((m = Pattern.compile("^\\#=GR\\s+(\\S+)\\s+SS\\s+(\\S+)")
                        .matcher(line))
                        .find()) {
                    acc = m.group(1);
                    String ss = m.group(2);
                    flags.put("SS", 1);
                    if (!seqs.containsKey(id)) {
                        new_seq = new AlignSeq(acc, 1, new_id, ss, "align_ss");
                        new_id++;
                        seqs.put(id, new_seq);
                    } else {
                        if (seqs.get(id).getAlignSs() != null) {
                            String old_ss = seqs.get(id).getAlignSs();
                            seqs.get(id).setAlignSs(old_ss + ss);
                        } else {
                            seqs.get(id).setAlignSs(ss);
                        }
                    }
                } else if ((Pattern.compile("^\\#")
                        .matcher(line))
                        .matches()) {
                    System.out.println("Unrecognized : " + line + "\n");
                } else if ((m = Pattern.compile("^(\\S+)\\s+(\\S+)")
                        .matcher(line))
                        .find()) {
                    acc = m.group(1);
                    String seq = m.group(2);
                    if (!last_acc.equals(acc)) {
                        last_acc = acc;
                        last_count = 0;
                        id = acc;
                    } else {
                        last_count++;
                        id = acc + last_count;
                    }
                    if (!seqs.containsKey(id)) {
                        new_seq = new AlignSeq(acc, 1, new_id, seq, "align_seq");
                        new_id++;
                        seqs.put(id, new_seq);
                    } else {
                        if (seqs.get(id).getAlignSeq() != null) {
                            String old_seq = seqs.get(id).getAlignSeq();
                            seqs.get(id).setAlignSeq(old_seq + seq);
                        } else {
                            seqs.get(id).setAlignSeq(seq);
                        }
                    }
                } else {
                    last_acc = "";
                    last_count = 0;
                }
            }
        } catch (IOException e) {
            System.out.println("Can't open file: " + e);
        }

        float sum_score = 0;
        int sum_weight = 0;
        int sum_len = 0;

//       lamdba can't sum local variable for sum_ variables 
//       seqs.forEach((k, v) -> {
//            v.setAlignSeq(remove_gap(v.getAlignSeq()));
//            if (flags.containsKey("SS") && flags.get("SS") == 1) {
//                v.setAlignSs(remove_gap(v.getAlignSs()));
//            }
//        });
//        
        for (Map.Entry<String, AlignSeq> entry : seqs.entrySet()) {
            entry.getValue().setAlignSeq(remove_gap(entry.getValue().getAlignSeq()));
            if (flags.containsKey("SS") && flags.get("SS") == 1) {
                entry.getValue().setAlignSs(remove_gap(entry.getValue().getAlignSs()));
            }
            if (flags.containsKey("SC")) {
                sum_score += (entry.getValue().getScore())
                        * entry.getValue().getWeight();
            }
            sum_weight += entry.getValue().getWeight();  //weight is int?
            sum_len += (entry.getValue().getSeq().length())
                    * entry.getValue().getWeight();
        }

        if (sum_weight == 0) {
            sum_weight = seqs.size();
        }

        return new Alignment(
                seqs, flags, ss_cons, rf, sum_score / sum_weight, sum_weight, sum_len / sum_weight);
    }

    /*
    **
    * @Param file_name
    * @return Alignment
     */
    public static Alignment read_rfam(String rfam_file) {
        Alignment alignment = read_stockholm(rfam_file);
        HashMap<String, AlignSeq> seqs = alignment.getSeqs();
        Matcher m;
        for (Map.Entry<String, AlignSeq> entry : seqs.entrySet()) {
            if ((m = Pattern.compile("(\\S+)\\/(\\d+)-(\\d+)")
                    .matcher(entry.getKey())).find()) {
                entry.getValue().setStart(Integer.parseInt(m.group(2)));
                entry.getValue().setEnd(Integer.parseInt(m.group(3)));
                seqs.put(m.group(1), entry.getValue());
                seqs.remove(entry.getKey());
            }
        }
        alignment.setSeqs(seqs);
        return alignment;
    }

    /*
    **
    * @Param file_name
    * 
     */
    public static void write_stockholm(Alignment alignment, String file_name) {
        if (file_name != null && !file_name.isEmpty()) {
            Path out = Paths.get(file_name);
            try {
                Files.write(
                        out,
                        ("# STOCKHOLM 1.0" + System.lineSeparator() + System.lineSeparator()).getBytes()
                );
                //different variables
                int line_len = 80;
                HashMap<String, AlignSeq> seqs = alignment.getSeqs();
                if (seqs.isEmpty()) {
                    throw new MyException("Empty alignments seqs");
                }
                HashMap<String, Integer> flags = alignment.getFlags();
                String ss_cons = alignment.getSsCons();
                String rf = alignment.getRf();
                //get max_name_length
                int max_name_length = seqs.entrySet().stream()
                        .map(entry -> entry.getValue().getAcc().length())
                        .max(Comparator.comparing(Integer::valueOf))
                        .get();
                max_name_length++;
                // printf("'%-20s'\n", "Hello"); 'Hello   
                String name_format = "%-" + max_name_length + "s";  // indeed a right padding

                if (flags.containsKey("WGT")) {
                    Files.write(out,
                            //print base on AlignSeq.id value
                            (Iterable<String>) () -> seqs.entrySet().stream()
                                    //.sorted((e1, e2)->
                                    //{return Integer.compare(e1.getValue().getId(), e2.getValue().getId());})
                                    .sorted(Comparator.comparing(e -> e.getValue().getId()))
                                    .map(e -> "#=GS "
                                    + String.format(name_format, e.getValue().getAcc()) //right space padding
                                    + "WT\t"
                                    + e.getValue().getWeight()
                                    + System.lineSeparator()
                                    ).iterator(),
                            StandardOpenOption.APPEND
                    );
                    //"\n\n"
                    Files.write(out, System.lineSeparator().getBytes(), StandardOpenOption.APPEND);
                    Files.write(out, System.lineSeparator().getBytes(), StandardOpenOption.APPEND);
                }

                if (flags.containsKey("DE")) {
                    Files.write(out,
                            //print base on AlignSeq.id value
                            (Iterable<String>) () -> seqs.entrySet().stream()
                                    //.sorted((e1, e2)->
                                    //{return Integer.compare(e1.getValue().getId(), e2.getValue().getId());})
                                    .sorted(Comparator.comparing(e -> e.getValue().getId()))
                                    .map(e -> "#=GS "
                                    + String.format(name_format, e.getValue().getAcc()) //right space padding
                                    + "DE\t"
                                    + e.getValue().getWeight()
                                    + System.lineSeparator()
                                    ).iterator(),
                            StandardOpenOption.APPEND
                    );
                    //"\n\n"
                    Files.write(out, System.lineSeparator().getBytes(), StandardOpenOption.APPEND);
                    Files.write(out, System.lineSeparator().getBytes(), StandardOpenOption.APPEND);
                }

                int ss_len = seqs.entrySet().iterator().next().getValue().getAlignSs().length();
                int seq_len = ss_len;
                String gap1 = "            ";
                String gap2 = "     ";

                for (int len = 0; len < seq_len; len += line_len) {
                    int l = line_len;
                    l = (ss_len - len < l) ? ss_len - len : l;
                    final int bindex = len;    //begin index , note perl offset = begin index
                    final int eindex = len + l; //end index, note perl length+begin index = end index

                    Files.write(out,
                            //print base on AlignSeq.id value
                            (Iterable<String>) () -> seqs.entrySet().stream()
                                    .sorted(Comparator.comparing(e -> e.getValue().getId()))
                                    .map(e -> String.format(name_format, e.getValue().getAcc()) //right space padding
                                    + gap1 + e.getValue().getAlignSeq().substring(bindex, eindex) //using substring here
                                    + System.lineSeparator()
                                    + ((flags.containsKey("SS")) ? "#=GR "
                                    + String.format(name_format, e.getValue().getAcc())
                                    + "SS" + gap2
                                    + e.getValue().getAlignSs().substring(bindex, eindex)
                                    + System.lineSeparator()
                                    : "")
                                    ).iterator(),
                            StandardOpenOption.APPEND);

                    if (flags.containsKey("SS_cons")) {
                        Files.write(out,
                                (String.format(name_format, "#=GC SS_cons")
                                        + gap1
                                        + ss_cons.substring(bindex, eindex)
                                        + System.lineSeparator()).getBytes(),
                                StandardOpenOption.APPEND);
                    }

                    if (flags.containsKey("RF")) {
                        Files.write(out,
                                (String.format(name_format, "#=GC RF")
                                        + gap1
                                        + rf.substring(bindex, eindex)
                                        + System.lineSeparator()).getBytes(),
                                StandardOpenOption.APPEND);
                    }
                    Files.write(out, System.lineSeparator().getBytes(), StandardOpenOption.APPEND);
                }
                Files.write(out, ("//" + System.lineSeparator()).getBytes(), StandardOpenOption.APPEND);
            } catch (IOException e) {
                System.out.println("Can't open file for writing " + e);
            } catch (MyException ex) {
                System.out.println(ex.getMessage());
            }
        } else {
            System.out.println("Write Stockholm failed, the file name is null or empty.");
        }
    }

    // we add nice extra method, name can explain
    public static String rpad(String s, int n) {
        return String.format("%-" + n + "s", s);
    }

    public static String lpad(String s, int n) {
        return String.format("%" + n + "s", s);
    }

    public static void main(String[] args) {
        HashMap<String, Seq> result = Io.read_fasta("test/cmf/data/example.fasta");
        //using stream print out Map
        //        result.entrySet().stream().
        //                forEach(e -> System.out.println(e.getValue().getSeqString()));
        write_fasta("", result);
        write_fasta("test/cmf/data/test.fasta", result);

        //test Files.write()
        /*
        Path path = Paths.get("test/try.txt");
        try {
            Files.write(path, "some test content...\n".getBytes());
            Iterable<String> iterable = Arrays.asList("line1", "line2");
            Files.write(path, iterable, StandardOpenOption.APPEND);
            byte[] bytes = Files.readAllBytes(path);
            System.out.println(new String(bytes));
        } catch (IOException ex) {
            Logger.getLogger(Io.class.getName()).log(Level.SEVERE, null, ex);
        }
         */
    }

}
