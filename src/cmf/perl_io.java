package cmf;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.regex.*;

/**
 * replace bin/io.pl
 */
public class perl_io {

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
    public static Map<String, Seq> read_fasta(String file_name) {
        Path fasta = Paths.get(file_name);
        Map<String, Seq> seqs = new HashMap<String, Seq>();
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
                    if ((m = p.matcher(line)).find()) {
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
            System.out.println(e);
        }

        return seqs;
    }

    /**
     * Output fasta to a file or screen based on fasta data input order
     *
     * @param file_name filename to read
     * @param segs is a Map<String, Seq>
     */
    public static void write_fasta(String file_name, Map<String, Seq> seqs) {
        //sort using treemap
        TreeMap<Integer, Seq> tree_seqs = new TreeMap<Integer, Seq>();
        //move into TreeMap
        seqs.entrySet().stream()
                .forEach(e -> tree_seqs.put(new Integer(e.getValue().getId()), e.getValue())
                );
        // now try to write file or print to screen
        if (file_name != null && !file_name.isEmpty()) {
            Path fasta_out = Paths.get(file_name);
            try {
                Files.write(fasta_out,
                        (Iterable<String>) tree_seqs.entrySet().stream()
                                .map(e -> ">" + e.getValue().getAcc() + "\t"
                                + ((e.getValue().getDesc().equals(e.getValue().getAcc())) ? "" : e.getValue().getDesc()) + "\n"
                                + e.getValue().getSeq()
                                )::iterator);
            } catch (IOException e) {
                System.out.println(e);
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

    public static void main(String[] args) {
        Map<String, Seq> result = perl_io.read_fasta("test/cmf/data/example.fasta");
        //using stream print out Map
        //        result.entrySet().stream().
        //                forEach(e -> System.out.println(e.getValue().getSeqString()));
        write_fasta("", result);
        write_fasta("test/cmf/data/test.fasta", result);

    }

}
