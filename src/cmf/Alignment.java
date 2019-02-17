package cmf;

import java.util.Map;

/**
 *
 * struct 'Alignment'
 *
 * STOCKHOLM FORMAT
 */
public class Alignment {

    Map<String, Seq> seqs;
    Map<String, Integer> flags;
    String ss_cons;
    String rf;
    float score;
    String weight;
    int len;

    public Alignment(
            Map<String, Seq> v_seqs,
            Map<String, Integer> v_flags,
            String v_ss_cons,
            String v_rf,
            float v_score,
            String v_weight,
            int v_len) {
        seqs = v_seqs;
        flags = v_flags;
        ss_cons = v_ss_cons;
        rf = v_rf;
        score = v_score;
        weight = v_weight;
        len = v_len;
    }

    public Map<String, Seq> getSeqs() {
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

    public String getWeight() {
        return weight;
    }

    public int getint() {
        return len;
    }
}
