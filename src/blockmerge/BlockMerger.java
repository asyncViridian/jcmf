package blockmerge;

import blockmerge.util.AlignmentBlock;
import blockmerge.util.MAFReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public class BlockMerger {

    /**
     * The type of block-merging to do.
     * <p>
     * Bases = each output alignment block contains a certain number of bases.
     * <p>
     * Blocks = each output alignment block contains a certain number of
     * original alignment blocks.
     */
    public static final MergeType type = MergeType.BLOCKS;
    /**
     * For merge type block: the number of blocks to contain in each output
     * alignment block
     */
    public static final int NUM_BLOCKS_PER_OUTPUT = 2;
    /**
     * For merge type bases: the number of bases to contain in each output
     * alignment block
     */
    public static final int NUM_BASES_PER_OUTPUT = 300;
    /**
     * For merge type bases: the number of bases to shift forwards by in each
     * block
     */
    public static final int NUM_BASES_INCREMENT = 50;
    /**
     * The maximum gap length between two blocks that still allows them to be
     * considered mergeable.
     */
    public static final int GAP_THRESHOLD = 300;

    // TODO add a parameter for "minimum relevant species"

    public static void main(String[] args) throws IOException {
        // TODO take this from input?
        String mafSrc = "toy.maf";
        //"multiz100way_chr12_62602752-62622213.maf";
        // TODO take this from input?
        String outputPrefix = "toy";
        //"m100_chr12_62602752-62622213";
        // TODO take this from input?
        String outputDir = "output";
        MAFReader reader = new MAFReader(mafSrc);

        if (type == MergeType.BASES) {
            // TODO do base-based merging
            throw new RuntimeException("not implemented yet");
        } else if (type == MergeType.BLOCKS) {
            LinkedList<AlignmentBlock> current = new LinkedList<>();
            BigInteger i = BigInteger.ONE;
            while (reader.hasNext()) {
                // add the next block
                current.add(reader.next());
                if (current.size() == 1) {
                    // abort if we only have one block
                    // (don't merge single block)
                    continue;
                }
                // keep the set of working blocks small
                // (we only write n blocks at a time)
                if (current.size() > NUM_BLOCKS_PER_OUTPUT) {
                    current.remove();
                }

                // find the set of mergeable species
                List<String> speciesToMerge = new LinkedList<>();
                for (String species : current.getFirst().getSpecies()) {
                    // initialize iterators through the blocks
                    Iterator<AlignmentBlock> secondIt = current.iterator();
                    secondIt.next();
                    Iterator<AlignmentBlock> firstIt = current.iterator();
                    // retrieve the list of species that can be included in
                    // this set of alignment blocks
                    boolean mergeable = true;
                    while (secondIt.hasNext()) {
                        if (mergeable && !isMergeable(firstIt.next(),
                                                      secondIt.next(),
                                                      species)) {
                            mergeable = false;
                        }
                    }
                    if (mergeable) {
                        speciesToMerge.add(species);
                    }
                }

                // create output FASTA file
                Path file = Paths.get(outputDir,
                                      outputPrefix + "_" + i.toString() +
                                              ".fasta");
                Files.deleteIfExists(file);
                Files.createFile(file);
                BufferedWriter writer = Files.newBufferedWriter(file);

                // write FASTA lines to disk
                for (String species : speciesToMerge) {
                    // write each species sequence
                    // write the species and originating chromosome
                    AlignmentBlock.Sequence first =
                            current.getFirst().sequences.get(
                            species);
                    AlignmentBlock.Sequence last =
                            current.getLast().sequences.get(
                            species);
                    // build the sequence name (species, chrom, strand, coords)
                    StringBuilder speciesHeader = new StringBuilder();
                    speciesHeader.append(species);
                    speciesHeader.append(".");
                    speciesHeader.append(first.section);
                    speciesHeader.append(
                            "(" + (first.strand ? "+" : "-") + ")");
                    speciesHeader.append(":");
                    speciesHeader.append(first.start);
                    speciesHeader.append("-");
                    speciesHeader.append(last.start.add(last.size));
                    writeToFasta(writer, speciesHeader.toString(), null);
                    for (AlignmentBlock block : current) {
                        // write the contents of each block for each species
                        AlignmentBlock.Sequence seq = block.sequences.get(
                                species);
                        if (!seq.isGap) {
                            // if we have contents write the contents directly
                            writeToFasta(writer, null, seq.contents);
                            if (!seq.right.length.equals(BigInteger.ZERO)) {
                                // if there is nonzero amount of context
                                // bases (to the right)
                                writeToFasta(writer, null,
                                             repeat("N", seq.right.length));
                            }
                        } else {
                            // if this is gap (unknown reference?) fill with Ns
                            writeToFasta(writer, null, repeat("N", seq.size));
                        }
                    }
                }

                // increment the file number
                i = i.add(BigInteger.ONE);

                writer.close();
            }
        }
    }

    /**
     * Return a string consisting of src repeated num times
     *
     * @param src String to repeat
     * @param num number of times to repeat this string
     * @return src repeated num times
     */
    public static String repeat(String src, BigInteger num) {
        StringBuilder builder = new StringBuilder();
        for (BigInteger i = BigInteger.ZERO; i.compareTo(num) < 0; i = i.add(
                BigInteger.ONE)) {
            builder.append(src);
        }
        return builder.toString();
    }

    /**
     * Append a FASTA-format sequence with the given title and the given
     * content to the end of the given file.
     *
     * @param writer   writer for file to append to
     * @param header   title of the sequence. Iff null then just append to
     *                 the last sequence, do not write a new sequence
     * @param sequence content of the sequence. Iff null then do not write
     *                 any sequence content.
     */
    public static void writeToFasta(BufferedWriter writer, String header,
                                    String sequence) {
        try {
            if (header != null) {
                writer.write(">" + header + "\n");
            }
            if (sequence != null) {
                writer.write(sequence + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Returns true iff the two given AlignmentBlocks can be merged for the
     * given species.
     *
     * @param first   the first alignment block to merge (the one that comes
     *                earlier)
     * @param second  the second alignment block to merge (the one that comes
     *                later)
     * @param species the species to attempt to merge these two blocks for
     * @return true iff the two given blocks can merge for the given species
     */
    public static boolean isMergeable(AlignmentBlock first,
                                      AlignmentBlock second, String species) {
        AlignmentBlock.Sequence firstSeq = first.sequences.get(species);
        AlignmentBlock.Sequence secondSeq = second.sequences.get(species);

        // test if they are in different non-scaffold chromosomes:
        String firstSec = firstSeq.section;
        String secondSec = secondSeq.section;
        if (!firstSec.contains("scaffold") && !secondSec.contains("scaffold")) {
            if (!firstSec.equals(secondSec)) {
                // If neither chromosome is a scaffold and they are in
                // different chromosomes, then they are not mergeable
                return false;
            }
        }

        // test if they are on different strands:
        if (firstSeq.strand != secondSeq.strand) {
            // If they are on different strands, they are not mergeable
            return false;
        }

        // test if they are at incomparable locations / overlapping:
        if (firstSeq.start.add(firstSeq.size).compareTo(secondSeq.start) < 0) {
            // if the ending of the first block is after
            // the beginning of the second block
            return false;
            // NOTE I don't think this will happen with my input MAFS
            // but it is probably still good to check for
        }

        // test if the gap between them is too long:
        if (firstSeq.isGap && firstSeq.gapType.equals(
                AlignmentBlock.Sequence.GapType.C)) {
            // if first is a gap sequence (with no bases in the section)
            // and it is too long to be included
            if (firstSeq.size.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                return false;
            }
        }
        if (secondSeq.isGap) {
            // if second is a gap sequence and it is too long to be included
            if (secondSeq.size.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                return false;
            }
        }
        if (!firstSeq.isGap && !secondSeq.isGap) {
            // if neither is a gap, but the intervening "gap" is too long
            if (firstSeq.right.length.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                // the context "after" the first block
                return false;
            } else if (secondSeq.left.length.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                // the context "before" the second block
                // should be the same as thing above???
                return false;
            }
            // TODO what else should I check here?
        }

        // TODO are there more criteria?

        return true;
    }

    public enum MergeType {
        BASES, BLOCKS
    }
}
