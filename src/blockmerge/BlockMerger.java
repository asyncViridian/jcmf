package blockmerge;

import blockmerge.util.AlignmentBlock;
import blockmerge.util.MAFReader;
import com.sun.scenario.effect.Merge;

import java.io.IOException;
import java.math.BigInteger;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

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

    public static void main(String[] args) throws IOException {
        // TODO take this from input?
        String mafSrc = "multiz100way_chr12_62602752-62622213.maf";
        // TODO take this from input?
        String outputPrefix = "multiz100way_chr12_62602752-62622213";
        // TODO take this from input?
        String outputDir = "output";
        MAFReader reader = new MAFReader(mafSrc);

        if (type == MergeType.BASES) {
            // TODO do base-based merging
            throw new RuntimeException("not implemented yet");
        } else if (type == MergeType.BLOCKS) {
            LinkedList<AlignmentBlock> current = new LinkedList<>();
            // TODO add block join logic
            while (reader.hasNext()) {
                // add the next block
                current.add(reader.next());
                // keep the set of working blocks small
                // (we only write n blocks at a time)
                if (current.size() > NUM_BLOCKS_PER_OUTPUT) {
                    current.remove();
                }
                // TODO write to disk
            }
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
        } else {
            // if first isn't a gap sequence, but the intervening (non-block)
            // space is too long
            // todo
        }
        if (secondSeq.isGap) {
            // if second is a gap sequence and it is too long to be included
            if (secondSeq.size.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                return false;
            }
        }
        // TODO

        // TODO are there more criteria?

        return true;
    }

    public enum MergeType {
        BASES, BLOCKS
    }
}
