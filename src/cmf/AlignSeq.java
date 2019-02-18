package cmf;

import java.util.HashMap;

/**
 *
 * struct 'AlignSeq'
 *
 * STOCKHOLM FORMAT
 */
public class AlignSeq {

    String acc;
    int id;
    String desc;
    String weight;
    String seq;
    String align_seq;
    String ss;
    String align_ss;
    HashMap<Integer, Integer> align_map;
    HashMap<Integer, Integer> rev_map;
    Integer start;
    Integer end;
    float score;

    public AlignSeq(
            String v_acc,
            int v_id,
            String v_desc,
            String v_weight,
            String v_seq,
            String v_align_seq,
            String v_ss,
            String v_align_ss,
            HashMap<Integer, Integer> v_align_map,
            HashMap<Integer, Integer> v_rev_map,
            Integer v_start,
            Integer v_end,
            float v_score) {
        acc = v_acc;
        id = v_id;
        desc = v_desc;
        weight = v_weight;
        seq = v_seq;
        align_seq = v_align_seq;
        ss = v_ss;
        align_ss = v_align_ss;
        align_map = v_align_map;
        rev_map = v_rev_map;
        start = v_start;
        end = v_end;
        score = v_score;
    }

    public String getAcc() {
        return acc;
    }

    public void setAcc(String v_acc) {
        acc = v_acc;
    }

    public int getId() {
        return id;
    }

    public void getId(int v_id) {
        id = v_id;
    }

    public String getDesc() {
        return desc;
    }

    public void setDesc(String v_desc) {
        desc = v_desc;
    }

    public String getWeight() {
        return weight;
    }

    public void setWeight(String v_weight) {
        weight = v_weight;
    }

    public String getSeq() {
        return seq;
    }

    public void setSeq(String v_seq) {
        seq = v_seq;
    }

    public String getAlignSeq() {
        return align_seq;
    }

    public void setAlignSeq(String v_align_seq) {
        align_seq = v_align_seq;
    }

    public String getSs() {
        return ss;
    }

    public void setSs(String v_ss) {
        ss = v_ss;
    }

    public String getAlignSs() {
        return align_ss;
    }

    public void setAlignSs(String v_align_ss) {
        align_ss = v_align_ss;
    }

    public HashMap<Integer, Integer> getAlignMap() {
        return align_map;
    }

    //    The simplest way to clone a simple structure. 
    //    This work fine as long as type are immutable
    public void setAlignMap(HashMap<Integer, Integer> v_align_map) {
        align_map = new HashMap<>(v_align_map);
    }

    public HashMap<Integer, Integer> getRevMap() {
        return rev_map;
    }

    public void setRevMap(HashMap<Integer, Integer> v_rev_map) {
        rev_map = new HashMap<>(v_rev_map);
    }

    public Integer getStart() {
        return start;
    }

    public void getStart(Integer v_start) {
        start = v_start;
    }

    public Integer getEnd() {
        return end;
    }

    public void setEnd(Integer v_end) {
        end = v_end;
    }

    public float getScore() {
        return score;
    }

    public void setScore(float v_score) {
        score = v_score;
    }

}
