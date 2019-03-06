package cmf;

import java.util.HashMap;

/**
 *
 * struct 'Alignment'
 *
 * STOCKHOLM FORMAT
 */
public class Alignment {

    private HashMap<String, AlignSeq> seqs;
    private HashMap<String, Integer> flags;
    private String ss_cons;
    private String rf;
    private float score;
    private double weight;
    private float len;

    public Alignment(
            HashMap<String, AlignSeq> v_seqs,
            HashMap<String, Integer> v_flags,
            String v_ss_cons,
            String v_rf,
            float v_score,
            double v_weight,
            float v_len) {
        seqs = v_seqs;
        flags = v_flags;
        ss_cons = v_ss_cons;
        rf = v_rf;
        score = v_score;
        weight = v_weight;
        len = v_len;
    }

    public Alignment(
            HashMap<String, AlignSeq> v_seqs,
            HashMap<String, Integer> v_flags,
            String v_ss_cons,
            String v_rf) {
        seqs = v_seqs;
        flags = v_flags;
        ss_cons = v_ss_cons;
        rf = v_rf;
    }

    public HashMap<String, AlignSeq> getSeqs() {
        return seqs;
    }

    public void setSeqs(HashMap<String, AlignSeq> seqs) {
        this.seqs = seqs;
    }

    public HashMap<String, Integer> getFlags() {
        return flags;
    }

    public void setFlags(HashMap<String, Integer> flags) {
        this.flags = flags;
    }

    public String getSsCons() {
        return ss_cons;
    }

    public String getRf() {
        return rf;
    }

    public float getScore() {
        return score;
    }

    public double getWeight() {
        return weight;
    }

    public float getint() {
        return len;
    }
}
