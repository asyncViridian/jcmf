package cmf;

import java.io.IOException;
import java.nio.file.*;
import java.util.Collections;
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
     * @param n  number of times to repeat ch
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
     * @param n   length to make the result sequence reach
     * @param ch  sequence to pad with
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
     * Reads the contents of the given file(name) into a mapping from sequence name to sequence struct
     *
     * @param file_name filename to read
     * @return map from sequence nickname to struct representing sequence
     */
    public static Map<String, Seq> read_fasta(String file_name) {
        Path fasta = Paths.get(file_name);
        Map<String, Seq> seqs = new HashMap<String, Seq>();
        int i = 0;
        //pattern for > line
        Pattern p1 = Pattern.compile("^>(\\S+)\\s*$");
        Pattern p2 = Pattern.compile("^>(\\S+)\\s+(\\S+.*)");
        Matcher m1;
        Matcher m2;
        String acc = "";
        String desc = "";
        StringBuilder seq = new StringBuilder();
        try {
            List<String> lines = Files.readAllLines(fasta);
            acc = lines.get(0).split(">")[1];
            for (String line : lines) {
                if (line.length() > 0 && line.charAt(0) == '>') {
                    // begin next entry
                    seqs.put(acc, new Seq(acc, i, acc, seq.toString()));
                    acc = line.split(">")[1];
                    i++;
                    seq = new StringBuilder();
                } else {
                    // continue an entry
                    seq.append(line);
                }
//                //System.out.println(line);
//                m1 = p1.matcher(line);
//                m2 = p2.matcher(line);
//                if (m1.find() || m2.find()) {
//
//                    if (m1.find()) {
//                        acc = m1.group(0);
//                        desc = "";
//
//                    } else if (m2.find()) {
//                        acc = m2.group(0);
//                        desc = m2.group(1);
//
//                    } else {
//                        seq = seq.append(line);
//                    }
//
//                }
            }
         seqs.put(acc, new Seq(acc, i, acc, seq.toString()));
        } catch (IOException e) {
            System.out.println(e);
        }

        return seqs;
    }
    public static void main(String[] args) {
        Map<String, cmf.Seq> result = perl_io.read_fasta("test/cmf/data/example.fasta");
    }

}
