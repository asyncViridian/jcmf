package cmf;

import java.util.Map;

/**
 *
 * struct 'Alignment'
 *
 * STOCKHOLM FORMAT
 */
public class Alignment {

    Map<String, AlignSeq> seqs;
    Map<String, Integer> flags;
    String ss_cons;
    String rf;
    float score;
    int weight;
    float len;

    public Alignment(
            Map<String,AlignSeq> v_seqs,
            Map<String, Integer> v_flags,
            String v_ss_cons,
            String v_rf,
            float v_score,
            int v_weight,
            int v_len) {
        seqs = v_seqs;
        flags = v_flags;
        ss_cons = v_ss_cons;
        rf = v_rf;
        score = v_score;
        weight = v_weight;
        len = v_len;
    }

    public Map<String, AlignSeq> getSeqs() {
        return seqs;
    }

    public Map<String, Integer> getFlags() {
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
