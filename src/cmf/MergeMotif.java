package cmf;

/**
 *
 * @MergeMotif for merge_motif return multiple values
 */
public class MergeMotif {

//    String motif1;
//    String motif2;
//    int num_seq;
//    double score;
//    int gao;
//    String overlap;
//    int weight;
    private int start;
    private int end;
    private String seq1;
    private String ss1;
    private String gap_seq;
    private String gap_ss;
    private String seq2;
    private String ss2;

    MergeMotif(
            int start,
            int end,
            String seq1,
            String ss1,
            String gap_seq,
            String gap_ss,
            String seq2,
            String ss2) {
        this.start = start;
        this.end = end;
        this.seq1 = seq1;
        this.ss1 = ss1;
        this.gap_seq = gap_seq;
        this.gap_ss = gap_ss;
        this.seq2 = seq2;
        this.ss2 = ss2;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getSeq1() {
        return seq1;
    }

    public void setSeq1(String seq1) {
        this.seq1 = seq1;
    }

    public String getSs1() {
        return ss1;
    }

    public void setSs1(String ss1) {
        this.ss1 = ss1;
    }

    public String getGapSeq() {
        return gap_seq;
    }

    public void setGapSeq(String gap_seq) {
        this.gap_seq = gap_seq;
    }

    public String getGapSs() {
        return gap_ss;
    }

    public void setGapSs(String gap_ss) {
        this.gap_ss = gap_ss;
    }

    public String getSeq2() {
        return seq2;
    }

    public void setSeq2(String seq2) {
        this.seq2 = seq2;
    }

    public String getSs2() {
        return ss2;
    }

    public void setSs2(String ss2) {
        this.ss2 = ss2;
    }
}
