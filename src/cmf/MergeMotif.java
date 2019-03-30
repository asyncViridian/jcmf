/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cmf;

/**
 *
 * @author james
 */
public class MergeMotif {

    private String motif1;
    private String motif2;
    private double num_seq;
    private double score;
    private double gap;
    private double overlap;
    private double weight;

    public MergeMotif(
            String motif1,
            String motif2,
            double num_seq,
            double score,
            double gap,
            double overlap,
            double weight
    ) {
        this.motif1 = motif1;
        this.motif2 = motif2;
        this.num_seq = num_seq;
        this.score = score;
        this.gap = gap;
        this.overlap = overlap;
        this.weight = weight;
    }

    public String getMotif1() {
        return motif1;
    }

    public void setMotif1(String motif1) {
        this.motif1 = motif1;
    }

    public String getMotif2() {
        return motif2;
    }

    public void setMotif2(String motif2) {
        this.motif2 = motif2;
    }

    public double getNum_seq() {
        return num_seq;
    }

    public void setNum_seq(double num_seq) {
        this.num_seq = num_seq;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public double getGap() {
        return gap;
    }

    public void setGap(double gap) {
        this.gap = gap;
    }

    public double getOverlap() {
        return overlap;
    }

    public void setOverlap(double overlap) {
        this.overlap = overlap;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

}
