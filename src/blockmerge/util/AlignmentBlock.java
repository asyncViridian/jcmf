package blockmerge.util;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;

public class AlignmentBlock {
    /**
     * All sequences and gapped sections in the alignment block.
     */
    public Map<String, Sequence> sequences = new HashMap<>();
    /**
     * The name of the species that this alignment block is "aligned to"
     */
    public String reference;
    /**
     * The score of the alignment block
     */
    public BigDecimal score;

    /**
     * Construct an AlignmentBlock from a set of Strings that comprise it,
     * starting from the first "a" line to the last "e" line
     *
     * @param lines
     */
    public AlignmentBlock(List<String> lines) {
        for (int i = 0; i < lines.size(); i++) {
            String line = lines.get(i);
            if (line.indexOf("##") == 0) {
                // this is a comment line, idc about this
                // TODO more general MAF support?
            } else if (line.charAt(0) == 'a') {
                // if this is the initial line
                // TODO make this more resilient for general MAF-format?
                this.score = new BigDecimal(line.split("\\s")[1].split("=")[1]);
            } else {
                // this is a non-initial line, so it corresponds sto a
                // sequence from some species
                String source = line.split("\\s+")[1];
                source = source.substring(0, source.indexOf('.'));
                if (line.charAt(0) == 's') {
                    if (i != lines.size() - 1 && lines.get(i + 1).charAt(
                            0) == 'i') {
                        // if there's a context line
                        this.sequences.put(source, new Sequence(line, lines.get(
                                i + 1)));
                        i++;
                    } else {
                        // no context line
                        this.sequences.put(source, new Sequence(line));
                        // make this the reference
                        this.reference = source;
                    }
                } else if (line.charAt(0) == 'e') {
                    // create a gap sequence
                    this.sequences.put(source, new Sequence(line));
                }
            }
        }
    }

    /**
     * Returns the set of species names that are represented in this
     * alignment block (sequences or gaps both)
     *
     * @return set of species in this block
     */
    public Set<String> getSpecies() {
        Set<String> result = new HashSet<>();
        for (Sequence s : this.sequences.values()) {
            result.add(s.src);
        }
        return result;
    }

    public static class Sequence {
        /**
         * The species that this sequence comes from
         */
        public String src;
        /**
         * The name of the chromosome (or scaffold) that this sequence comes
         * from
         */
        public String section;
        /**
         * The start index of this sequence in the chromosome relative to the
         * start position of the strand
         */
        public BigInteger start;
        /**
         * The number of non-gap bases in this sequence
         */
        public BigInteger size;
        /**
         * True means + strand, false means - strand
         */
        public boolean strand;
        /**
         * The total length of the source chromosome/scaffold that this
         * sequence comes from
         */
        public BigInteger srcSize;
        /**
         * The contents of this sequence, including gap characters etc.
         */
        public String contents;
        /**
         * The left-side context of the sequence. Is "C 0" if this sequence
         * is reference or gap.
         */
        public AdjacentDetail left;
        /**
         * The right-side context of the sequence. Is "C 0" if this sequence
         * is reference or gap.
         */
        public AdjacentDetail right;
        /**
         * True means this is the reference sequence (the first one listed in
         * a block). Iff true, then the AdjacentDetails must be "C 0 C 0"
         */
        public boolean isReference;
        /**
         * True means that this sequence is entirely empty throughout the
         * containing alignment block. Iff true, then has a non-null gapType,
         * the AdjacentDetails must be "C 0 C 0", and the contents is "".
         */
        public boolean isGap;
        /**
         * The type of gap that this sequence is (see GapType documentation
         * for more details)
         */
        public GapType gapType;

        /**
         * Construct a Sequence given each argument manually. Part of this
         * javadoc is copied from the UCSC MAF format specification.
         *
         * @param src         The name of one of the source sequences for the
         *                    alignment. For sequences that are resident in a
         *                    browser assembly, the form 'database
         *                    .chromosome' allows automatic creation of links
         *                    to other assemblies. Non-browser sequences are
         *                    typically reference by the species name alone.
         * @param section     The chromosome that this sequence comes from.
         * @param start       The start of the aligning region in the source
         *                    sequence. This is a zero-based number. If the
         *                    strand field is "-" then this is the start
         *                    relative to the reverse-complemented source
         *                    sequence.
         * @param size        The size of the aligning region in the source
         *                    sequence. This number is equal to the number of
         *                    non-dash characters in the alignment text field.
         * @param strand      Either "+" (true) or "-" (false). If "-", then
         *                    the alignment is to the reverse-complemented
         *                    source.
         * @param srcSize     The size of the entire source sequence, not
         *                    just the parts involved in the alignment.
         * @param contents    The nucleotides (or amino acids) in the
         *                    alignment and any insertions (dashes) as well.
         * @param left        Information about the context of the sequence
         *                    line: relationship between the sequence in this
         *                    block and the sequence in the block immediately
         *                    preceding it.
         * @param right       Information about the context of the sequence
         *                    line: relationship between the sequence in this
         *                    block and the sequence in the block immediately
         *                    following it.
         * @param isReference True iff this is the reference sequence for the
         *                    containing alignment block.
         * @param isGap       True iff this sequence is gapped for the
         *                    relevant alignment block.
         * @param gapType     Null iff isGap is true, otherwise equals the
         *                    type of gap that this sequence is.
         * @see
         * <a href="https://genome.ucsc.edu/FAQ/FAQformat.html#format5">the UCSC MAF spec</a>
         */
        public Sequence(String src, String section, BigInteger start,
                        BigInteger size, boolean strand, BigInteger srcSize,
                        String contents, AdjacentDetail left,
                        AdjacentDetail right, boolean isReference,
                        boolean isGap, GapType gapType) {
            this.src = src;
            this.section = section;
            this.start = start;
            this.size = size;
            this.strand = strand;
            this.srcSize = srcSize;
            this.contents = contents;
            this.left = left;
            this.right = right;
            this.isReference = isReference;
            this.isGap = isGap;
            this.gapType = gapType;
            checkRep();
        }

        /**
         * Construct a single-line MAF sequence
         *
         * @param onlyLine line for the sequence
         */
        public Sequence(String onlyLine) {
            String[] split = onlyLine.split("\\s+");
            this.buildBasic(split);
            this.left = new AdjacentDetail("C", BigInteger.valueOf(0));
            this.right = new AdjacentDetail("C", BigInteger.valueOf(0));
            if (onlyLine.charAt(0) == 's') {
                // Create a reference sequence
                this.isReference = true;
                this.isGap = false;
                this.gapType = null;
            } else if (onlyLine.charAt(0) == 'e') {
                // Create a gap sequence
                this.contents = "";
                this.isGap = true;
                this.gapType = GapType.valueOf(split[6]);
            }
            checkRep();
        }

        /**
         * Construct a non-reference MAF sequence (with context information)
         *
         * @param line1 first (s-) line for the sequence
         * @param line2 second (i-) line for the sequence
         */
        public Sequence(String line1, String line2) {
            String[] split1 = line1.split("\\s+");
            String[] split2 = line2.split("\\s+");
            this.buildBasic(split1);
            this.left = new AdjacentDetail(split2[2],
                                           new BigInteger(split2[3]));
            this.right = new AdjacentDetail(split2[4],
                                            new BigInteger(split2[5]));
            this.isReference = false;
            this.isGap = false;
            this.gapType = null;
        }

        /**
         * Populate the fields available in all non-gap sequence types given the
         * first line for the sequence
         *
         * @param firstLine first line of the sequence, whitespace-split
         */
        private void buildBasic(String[] firstLine) {
            this.src = firstLine[1].substring(0, firstLine[1].indexOf('.'));
            this.section = firstLine[1].substring(
                    firstLine[1].indexOf('.') + 1);
            this.start = new BigInteger(firstLine[2]);
            this.size = new BigInteger(firstLine[3]);
            this.strand = firstLine[4].equals("+");
            this.srcSize = new BigInteger(firstLine[5]);
            this.contents = firstLine[6];
        }

        private void checkRep() {
            // cannot both be a reference sequence and a gap sequence
            assert !(this.isReference && this.isGap);
            // check that isReference enforces the correct requirements
            if (this.isReference) {
                assert this.left.equals(
                        new AdjacentDetail("C", BigInteger.ZERO));
                assert this.right.equals(
                        new AdjacentDetail("C", BigInteger.ZERO));
            }
            // check that isGap enforces the correct requirements
            if (this.isGap) {
                assert this.left.equals(
                        new AdjacentDetail("C", BigInteger.ZERO));
                assert this.right.equals(
                        new AdjacentDetail("C", BigInteger.ZERO));
                assert this.gapType != null;
                assert this.contents.equals("");
            } else {
                assert this.gapType == null;
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
                 * this is the first sequence from this src chrom or scaffold
                 * but it is bridged by another alignment from a different
                 * chrom or scaffold
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

            public Type type;
            public BigInteger length;

            public AdjacentDetail(String type, BigInteger length) {
                this.type = Type.valueOf(type);
                this.length = length;
            }

            @Override
            public boolean equals(Object o) {
                if (!(o instanceof AdjacentDetail)) {
                    return false;
                }
                AdjacentDetail ad = (AdjacentDetail) o;
                return ad.type.equals(this.type) && ad.length.equals(
                        this.length);
            }

            @Override
            public int hashCode() {
                return Objects.hash(this.type, this.length);
            }
        }

        public enum GapType {
            /**
             * the sequence before and after is contiguous implying that this
             * region was either deleted in the source or inserted in the
             * reference sequence
             */
            C,
            /**
             * there are non-aligning bases in the source species between
             * chained alignment blocks before and after this block
             */
            I,
            /**
             * there are non-aligning bases in the source and more than 90%
             * of them are Ns in the source
             */
            M,
            /**
             * there are non-aligning bases in the source and the next
             * aligning block starts in a new chromosome or scaffold that is
             * bridged by a chain between still other blocks
             */
            n,
            /**
             * TODO
             * This shows up in the data but there's no documentation
             * on what this means. For now, treat it as a generic gap...
             */
            T
        }
    }

}
