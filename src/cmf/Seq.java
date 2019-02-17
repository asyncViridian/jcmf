package cmf;

/**
 * struct 'Seq'
 * 
 * FASTA format
 */
public class Seq {

    String acc;
    int id;
    String desc;
    String seq;

    public Seq(String v_acc, int v_id, String v_desc, String v_seq) {
        acc = v_acc;
        id = v_id;
        desc = v_desc;
        seq = v_seq;
    }

    public String getSeqString() {
        return "id=" + id + " acc=" + acc + " desc=" + desc + " \n  seq=" + seq;
    }
    
    public String getAcc() {
        return acc;
    }

    public int getId() {
        return id;
    }

    public String getDesc() {
        return desc;
    }

    public String getSeq() {
        return seq;
    }

}
