package cmf;

import java.util.Collections;

/**
 *
 * replace bin/io.pl
 */
public class perl_io {

    public static String remove_gap(String s) {
        String strNew = s.replace(".", "").replace("-", "");
        return strNew;
    }

    public static String make_string(int n, String ch) {
        //String strNew = new String(new char[n]).replace("\0", s); // before java 8
        String strNew = String.join("", Collections.nCopies(n, ch));  // Java 8, Java 11 repeat()
        return strNew;
    }

    public static String pad_string(String seq, int n, String ch, int dir) {
        String strNew = seq;
        int l = n - seq.length();
        if (l > 0) { //pad right
            if (dir == 1) {
                strNew = seq + make_string(n, ch);
            } // pad left
            else if (dir == 0) {
                strNew = make_string(n, ch) + seq;
            }
        }
        return strNew;
    }

    //FASTA FORMAT
    // some comment
    
    
    public static void main(String[] args) {
        System.out.println(pad_string("this is a good test", 35, "hello", 0));
    }

}
