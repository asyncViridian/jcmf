package blockmerge.util;

import java.math.BigInteger;
import java.util.List;

public class AlignmentBlock {
    private List<Sequence> sequences;

    public AlignmentBlock(List<Sequence> sequences) {
        this.sequences = sequences;
    }

    public static class Sequence {
        String src;
        BigInteger start;
        BigInteger size;
        /**
         * True means + strand, false means - strand
         */
        boolean strand;
        BigInteger srcSize;
        String contents;
        AdjacentDetail left;
        AdjacentDetail right;
        /**
         * True means this is the reference sequence (the first one listed in
         * a block)
         */
        boolean isReference;

        public Sequence(String src, BigInteger start, BigInteger size,
                        boolean strand, BigInteger srcSize, String contents,
                        AdjacentDetail left, AdjacentDetail right,
                        boolean isReference) {
            this.src = src;
            this.start = start;
            this.size = size;
            this.strand = strand;
            this.srcSize = srcSize;
            this.contents = contents;
            this.left = left;
            this.right = right;
            this.isReference = isReference;
        }

        /**
         * Construct a reference MAF sequence
         * @param onlyLine line for the sequence
         */
        public Sequence(String onlyLine) {
            String[] split = onlyLine.split("\\s");
            this.buildBasic(split);
            this.left = null;
            this.right = null;
            this.isReference = true;
        }

        /**
         * Construct a non-reference MAF sequence
         * @param line1 first line for the sequence
         * @param line2 second line for the sequence
         */
        public Sequence(String line1, String line2) {
            String[] split1 = line1.split("\\s");
            String[] split2 = line2.split("\\s");
            this.buildBasic(split1);
            this.left = new AdjacentDetail(split2[2],
                    new BigInteger(split2[3]));
            this.right = new AdjacentDetail(split2[4],
                    new BigInteger(split2[5]));
            this.isReference = false;
        }

        /**
         * Populate the fields available in all sequence types given the
         * first line for the sequence
         *
         * @param firstLine first line of the sequence, whitespace-split
         */
        private void buildBasic(String[] firstLine) {
            this.src = firstLine[1];
            this.start = new BigInteger(firstLine[2]);
            this.size = new BigInteger(firstLine[3]);
            this.strand = firstLine[4].equals("+");
            this.srcSize = new BigInteger(firstLine[5]);
            this.contents = firstLine[6];
        }
    }

    public static class AdjacentDetail {
        public enum Type {
            /**
             * the sequence before or after is contiguous with this block
             */
            C,
            /**
             * there are bases between the bases in this block and the one
             * before or after it
             */
            I,
            /**
             * this is the first sequence from this src chrom or scaffold
             */
            N,
            /**
             * this is the first sequence from this src chrom or scaffold but
             * it is bridged by another alignment from a different chrom or
             * scaffold
             */
            n,
            /**
             * there is missing data before or after this block (Ns in the
             * sequence)
             */
            M,
            /**
             * the sequence in this block has been used before in a previous
             * block (likely a tandem duplication)
             */
            T
        }

        Type type;
        BigInteger length;

        public AdjacentDetail(String type, BigInteger length) {
            this.type = Type.valueOf(type);
            this.length = length;
        }
    }
}
