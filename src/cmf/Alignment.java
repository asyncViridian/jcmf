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
    private int weight;
    private float len;

    public Alignment(
            HashMap<String, AlignSeq> v_seqs,
            HashMap<String, Integer> v_flags,
            String v_ss_cons,
            String v_rf,
            float v_score,
            int v_weight,
            float v_len) {
        seqs = v_seqs;
        flags = v_flags;
        ss_cons = v_ss_cons;
        rf = v_rf;
        score = v_score;
        weight = v_weight;
        len = v_len;
    }

    public HashMap<String, AlignSeq> getSeqs() {
        return seqs;
    }

    public void setSeqs(HashMap<String, AlignSeq> v_seqs) {
        this.seqs = new HashMap<>(v_seqs);
    }

    public HashMap<String, Integer> getFlags() {
        return flags;
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

    public int getWeight() {
        return weight;
    }

    public float getint() {
        return len;
    }
}
