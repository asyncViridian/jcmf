package blockmerge;

import blockmerge.util.AlignmentBlock;
import blockmerge.util.MAFReader;
import org.apache.commons.cli.*;

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
     * Options for the block merging
     */
    private static Options options;
    // TODO switch to using the Options values instead of class-level variables.

    /**
     * The type of block-merging to do.
     * <p>
     * Bases = each output alignment block contains a certain number of bases.
     * <p>
     * Blocks = each output alignment block contains a certain number of
     * original alignment blocks.
     */
    private static MergeType MERGE_TYPE;
    private static final MergeType MERGE_TYPE_DEFAULT = MergeType.BLOCKS;
    /**
     * For merge type block: the number of blocks to contain in each output
     * alignment block
     */
    private static int NUM_BLOCKS_PER_OUTPUT;
    private static final int NUM_BLOCKS_PER_OUTPUT_DEFAULT = 5;
    /**
     * For merge type bases: the number of bases to contain in each output
     * alignment block
     */
    private static int NUM_BASES_PER_OUTPUT;
    private static final int NUM_BASES_PER_OUTPUT_DEFAULT = 300;
    /**
     * For merge type bases: the number of bases to shift forwards by in each
     * block
     */
    private static int NUM_BASES_INCREMENT;
    private static final int NUM_BASES_INCREMENT_DEFAULT = 50;
    /**
     * The maximum gap length between two blocks that still allows them to be
     * considered mergeable.
     */
    private static int GAP_THRESHOLD;
    private static final int GAP_THRESHOLD_DEFAULT = 300;
    /**
     * The minimum number of species in a merged block that still allows the
     * block to be output (if a block has fewer species, it is considered
     * irrelevant/useless for the purposes of alignment)
     */
    private static int MIN_NUM_SPECIES;
    private static final int MIN_NUM_SPECIES_DEFAULT = 5;
    /**
     * The source MAF file with all alignment blocks we want to merge.
     */
    private static String srcName;
    /**
     * The directory to read input files from.
     */
    private static String srcDir;
    private static final String srcDir_DEFAULT = "";
    /**
     * The custom name part of the FASTA-format output files that we will write.
     */
    private static String outName;
    /**
     * The directory to output files into.
     */
    private static String outDir;
    private static final String outDir_DEFAULT = "";

    public static void main(String[] args) throws IOException, ParseException {
        // handle command-line argument processing :)
        // add all the arguments we need
        BlockMerger.options = new Options();
        {
            BlockMerger.options.addOption(
                    Option.builder("base")
                            .desc("Use base-counting based merging (by "
                                          + "default BlockMerger uses "
                                          + "block-based merging).")
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("numBlocksPerOutput")
                            .hasArg()
                            .desc("Number of source blocks to include in each" +
                                          " " +
                                          "merged block (if using " +
                                          "block-counting " +
                                          "merging). Defaults to " + NUM_BLOCKS_PER_OUTPUT_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("numBasesPerOutput")
                            .hasArg()
                            .desc("Number of bases to include in each "
                                          + "output block (if using "
                                          + "base-counting merging). " +
                                          "Defaults to " + NUM_BASES_PER_OUTPUT_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("numBasesIncrement")
                            .hasArg()
                            .desc("Number of bases to increment between "
                                          + "each output merged block (if" +
                                          " using base-counting merging)." +
                                          " Defaults to " + NUM_BASES_INCREMENT_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("gapThreshold")
                            .hasArg()
                            .desc("Maximum gap length between two source " +
                                          "alignment blocks before they " +
                                          "are determined unmergeable. " +
                                          "Defaults to " + GAP_THRESHOLD_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("minNumSpecies")
                            .hasArg()
                            .desc("Minimum number of species that are " +
                                          "mergeable to include in each " +
                                          "result alignment block. " +
                                          "Defaults to " + MIN_NUM_SPECIES_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder("sd")
                            .longOpt("srcDir")
                            .hasArg()
                            .desc("Directory to read input files from. " +
                                          "Defaults to "
                                          + (srcDir_DEFAULT.equals(
                                    "") ? "current directory" :
                                    "./" + srcDir_DEFAULT))
                            .build());
            BlockMerger.options.addOption(
                    Option.builder("s")
                            .longOpt("srcName")
                            .hasArg()
                            .desc("Input MAF file")
                            .required()
                            .build());
            BlockMerger.options.addOption(
                    Option.builder("od")
                            .longOpt("outDir")
                            .hasArg()
                            .desc("Directory to read input files from. " +
                                          "Defaults to "
                                          + (outDir_DEFAULT.equals(
                                    "") ? "current directory" :
                                    "./" + outDir_DEFAULT))
                            .build());
            BlockMerger.options.addOption(
                    Option.builder("o")
                            .longOpt("outName")
                            .hasArg()
                            .desc("Output FASTA files prefix")
                            .required()
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the merge type
            if (line.hasOption("base")) {
                BlockMerger.MERGE_TYPE = MergeType.BASES;
            } else {
                BlockMerger.MERGE_TYPE = MERGE_TYPE_DEFAULT;
            }
            // set the numBlocksPerOutput
            if (line.hasOption("numBlocksPerOutput")) {
                BlockMerger.NUM_BLOCKS_PER_OUTPUT = Integer.valueOf(
                        line.getOptionValue("numBlocksPerOutput"));
            } else {
                BlockMerger.NUM_BLOCKS_PER_OUTPUT =
                        NUM_BLOCKS_PER_OUTPUT_DEFAULT;
            }
            // set the numBasesPerOutput
            if (line.hasOption("numBasesPerOutput")) {
                BlockMerger.NUM_BASES_PER_OUTPUT = Integer.valueOf(
                        line.getOptionValue("numBasesPerOutput"));
            } else {
                BlockMerger.NUM_BASES_PER_OUTPUT = NUM_BASES_PER_OUTPUT_DEFAULT;
            }
            // set the numBasesIncrement
            if (line.hasOption("numBasesIncrement")) {
                BlockMerger.NUM_BASES_INCREMENT = Integer.valueOf(
                        line.getOptionValue("numBasesIncrement"));
            } else {
                BlockMerger.NUM_BASES_INCREMENT = NUM_BASES_INCREMENT_DEFAULT;
            }
            // set the gapThreshold
            if (line.hasOption("gapThreshold")) {
                BlockMerger.GAP_THRESHOLD = Integer.valueOf(
                        line.getOptionValue("gapThreshold"));
            } else {
                BlockMerger.GAP_THRESHOLD = GAP_THRESHOLD_DEFAULT;
            }
            // set the minNumSpecies
            if (line.hasOption("minNumSpecies")) {
                BlockMerger.MIN_NUM_SPECIES = Integer.valueOf(
                        line.getOptionValue("minNumSpecies"));
            } else {
                BlockMerger.MIN_NUM_SPECIES = MIN_NUM_SPECIES_DEFAULT;
            }
            // set the src directory
            if (line.hasOption("sd")) {
                BlockMerger.srcDir = line.getOptionValue("sd");
            } else {
                BlockMerger.srcDir = srcDir_DEFAULT;
            }
            // set the src filename
            BlockMerger.srcName = line.getOptionValue("s");
            // set the out directory
            if (line.hasOption("od")) {
                BlockMerger.outDir = line.getOptionValue("od");
            } else {
                BlockMerger.outDir = outDir_DEFAULT;
            }
            // set the out fileprefix
            BlockMerger.outName = line.getOptionValue("o");
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "BlockMerger [options] -s <filename> -o <fileprefix>",
                    options);
            return;
        }

        // start reading from the given files etc.
        MAFReader reader = new MAFReader(BlockMerger.srcDir,
                                         BlockMerger.srcName);

        if (MERGE_TYPE == MergeType.BASES) {
            // TODO write base-based merging...
            throw new RuntimeException("not implemented yet");
        } else if (MERGE_TYPE == MergeType.BLOCKS) {
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
                    while (secondIt.hasNext() && mergeable) {
                        if (!isMergeable(firstIt.next(),
                                         secondIt.next(),
                                         species)) {
                            mergeable = false;
                        }
                    }
                    if (mergeable) {
                        speciesToMerge.add(species);
                    }
                }

                // checks the # of merged species threshold
                if (speciesToMerge.size() < MIN_NUM_SPECIES) {
                    continue;
                }

                // create output FASTA file
                Path file = Paths.get(BlockMerger.outDir,
                                      BlockMerger.outName + "_" + i.toString() + ".fasta");
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
                // note that we have created a file
                System.out.println("Wrote " + file.toString());

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
    private static void writeToFasta(BufferedWriter writer, String header,
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
    private static boolean isMergeable(AlignmentBlock first,
                                       AlignmentBlock second, String species) {
        AlignmentBlock.Sequence firstSeq = first.sequences.get(species);
        AlignmentBlock.Sequence secondSeq = second.sequences.get(species);

        // test if the species is actually in both blocks
        if (firstSeq == null || secondSeq == null) {
            return false;
        }

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
        if (firstSeq.start.add(firstSeq.size).compareTo(secondSeq.start) > 0) {
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
            // TODO may need to add more checks here if bugs appear
        }

        // TODO check if there are more criteria?

        return true;
    }

    public enum MergeType {
        BASES, BLOCKS
    }
}
