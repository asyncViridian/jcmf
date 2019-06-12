package blockmerge;

import util.MAFAlignmentBlock;
import util.MAFReader;
import org.apache.commons.cli.*;
import util.SimpleScatterPlot;

import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
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

    private static final boolean REMOVE_GAPS = true;
    private static final boolean REMOVE_NEWLINES = true;

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
     * The maximum length of the reference genome in a resulting block that
     * will still allow the resulting block to be output.
     */
    private static int MAX_OUTPUT_LENGTH;
    private static int MAX_OUTPUT_LENGTH_DEFAULT = 5000;
    /**
     * The minimum length of the reference genome in a resulting block that
     * will still allow the resulting block to be output.
     */
    private static int MIN_OUTPUT_LENGTH;
    private static int MIN_OUTPUT_LENGTH_DEFAULT = 10;
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

    public static void main(String[] args) throws IOException {
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
                    Option.builder()
                            .longOpt("maxOutputLength")
                            .hasArg()
                            .desc("Maximum length of the reference genome " +
                                          "section in a resulting block that " +
                                          "will still allow the resulting " +
                                          "block to be output." +
                                          "Defaults to " + MAX_OUTPUT_LENGTH_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("minOutputLength")
                            .hasArg()
                            .desc("Minimum length of the reference genome " +
                                          "section in a resulting block that " +
                                          "will still allow the resulting " +
                                          "block to be output." +
                                          "Defaults to " + MIN_OUTPUT_LENGTH_DEFAULT)
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
                            .desc("Directory to write output files to. " +
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
            // set the maxOutputLength
            if (line.hasOption("maxOutputLength")) {
                BlockMerger.MAX_OUTPUT_LENGTH = Integer.valueOf(
                        line.getOptionValue("maxOutputLength"));
            } else {
                BlockMerger.MAX_OUTPUT_LENGTH = MAX_OUTPUT_LENGTH_DEFAULT;
            }
            // set the minOutputLength
            if (line.hasOption("minOutputLength")) {
                BlockMerger.MIN_OUTPUT_LENGTH = Integer.valueOf(
                        line.getOptionValue("minOutputLength"));
            } else {
                BlockMerger.MIN_OUTPUT_LENGTH = MIN_OUTPUT_LENGTH_DEFAULT;
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
            // TODO write percentage/threshold size based merging...
        } else if (MERGE_TYPE == MergeType.BLOCKS) {
            // Create overall statistics trackers
            // Track gap content (N bases) in each merged block
            // TODO make this filename argument-able???
            Path gapStatsFile = Paths.get(BlockMerger.outDir,
                                          "graph_postfilter_gapStatistics" +
                                                  ".png");
            SimpleScatterPlot gapStats = new SimpleScatterPlot(
                    gapStatsFile,
                    null,
                    "Merged sequence length",
                    "Gaps percentage");

            LinkedList<MAFAlignmentBlock> current = new LinkedList<>();
            BigInteger i = BigInteger.ONE;
            while (reader.hasNext()) {
                // add the next block
                current.add(reader.next());
                if (current.size() < NUM_BLOCKS_PER_OUTPUT) {
                    // abort if we have too few blocks merged
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
                    Iterator<MAFAlignmentBlock> secondIt = current.iterator();
                    secondIt.next();
                    Iterator<MAFAlignmentBlock> firstIt = current.iterator();
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

                // Filter:
                // checks the # of merged species threshold
                if (speciesToMerge.size() < MIN_NUM_SPECIES) {
                    continue;
                }

                // Filter:
                // require that the merged species includes human
                if (!speciesToMerge.contains("hg38")) {
                    continue;
                }

                // Filter:
                // check output against sequence length bounds
                // TODO: perhaps generalize to allow other human assemblies?
                // TODO: possibly integrate too-short blocks into other sections
                MAFAlignmentBlock.Sequence firstCheckSeqLen =
                        current.getFirst().sequences.get(
                                "hg38");
                MAFAlignmentBlock.Sequence lastCheckSeqLen =
                        current.getLast().sequences.get(
                                "hg38");
                BigInteger humanSeqLen = lastCheckSeqLen.start.add(
                        lastCheckSeqLen.size).subtract(firstCheckSeqLen.start);
                // check lower bound
                if (humanSeqLen.compareTo(
                        BigInteger.valueOf(MIN_OUTPUT_LENGTH)) < 0) {
                    continue;
                }
                // check upper bound
                if (humanSeqLen.compareTo(
                        BigInteger.valueOf(MAX_OUTPUT_LENGTH)) > 0) {
                    continue;
                }

                // Write to disk!:
                // create output FASTA file
                Path file = Paths.get(BlockMerger.outDir,
                                      BlockMerger.outName + "_" + i.toString() + ".fasta");
                Files.deleteIfExists(file);
                Files.createFile(file);
                BufferedWriter writer = Files.newBufferedWriter(file);

                // write FASTA lines to disk
                for (String species : speciesToMerge) {
                    // Initialize to track info for gapStats
                    BigInteger seqLength = BigInteger.ZERO;
                    BigInteger gapLength = BigInteger.ZERO;

                    // write each species sequence
                    // write the species and originating chromosome
                    MAFAlignmentBlock.Sequence first =
                            current.getFirst().sequences.get(
                                    species);
                    MAFAlignmentBlock.Sequence last =
                            current.getLast().sequences.get(
                                    species);

                    seqLength = last.start.add(last.size).subtract(first.start);

                    // if the merged sequence has a size of 0
                    if (seqLength.equals(BigInteger.ZERO)) {
                        // don't even print it
                        // skip this line of output
                        continue;
                    }

                    // build the sequence name
                    StringBuilder speciesHeader = new StringBuilder();
                    // write species
                    speciesHeader.append(species);
                    speciesHeader.append(":");
                    // write chromosome
                    speciesHeader.append(first.section);
                    speciesHeader.append(":");
                    // write strand
                    speciesHeader.append(first.strand ? "+" : "-");
                    speciesHeader.append(":");
                    // write start and end coords
                    speciesHeader.append(first.start);
                    speciesHeader.append("-");
                    speciesHeader.append(last.start.add(last.size));
                    speciesHeader.append(":");
                    // write individual start and end coords for each block
                    Iterator<MAFAlignmentBlock> it = current.iterator();
                    MAFAlignmentBlock b = it.next();
                    MAFAlignmentBlock.Sequence s = b.sequences.get(species);
                    speciesHeader.append(s.start);
                    speciesHeader.append("-");
                    speciesHeader.append(s.start.add(s.size));
                    while (it.hasNext()) {
                        b = it.next();
                        s = b.sequences.get(species);
                        speciesHeader.append(";");
                        speciesHeader.append(s.start);
                        speciesHeader.append("-");
                        speciesHeader.append(s.start.add(s.size));
                    }
                    writeFastaHeader(writer, speciesHeader.toString());

                    // write the contents of each block for each species
                    int numBlock = 1;
                    boolean lastBlockWasGap = false;
                    for (MAFAlignmentBlock block : current) {
                        // get the raw sequence and fill it up with N as needed
                        MAFAlignmentBlock.Sequence seq = block.sequences.get(
                                species);
                        if (!seq.isGap) {
                            // if we have contents write the contents directly
                            writeFastaContent(writer, seq.contents);
                            if (!seq.left.length.equals(
                                    BigInteger.ZERO) && numBlock != 1 && !lastBlockWasGap) {
                                // if there is nonzero amount of context
                                // bases (to the right)
                                // But we do not want to include gaps that are
                                // not in the "internals" of the merged area
                                writeFastaContent(writer, repeat("N",
                                                                 seq.left.length));
                                gapLength = gapLength.add(seq.left.length);
                            }
                            lastBlockWasGap = false;
                        } else {
                            // if this is gap (unknown reference?) fill with Ns
                            writeFastaContent(writer, repeat("N", seq.size));
                            gapLength = gapLength.add(seq.size);
                            lastBlockWasGap = true;
                        }
                        numBlock++;
                    }
                    writer.write("\n");
                    // write the stats for this sequence to tracker
                    // write the sequence length vs gap percentage
                    gapStats.addValue(new BigDecimal(seqLength),
                                      (new BigDecimal(gapLength))
                                              .multiply(BigDecimal.valueOf(100))
                                              .divide(new BigDecimal(seqLength),
                                                      RoundingMode.HALF_EVEN));
                }
                // note that we have created a file
                System.out.println("Wrote " + file.toString());

                // increment the file number
                i = i.add(BigInteger.ONE);

                writer.close();
            }

            // Output overall statistics graphics
            gapStats.write();
            System.out.println("Wrote statistics information");

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
     * Writes a line to the given writer starting with ">" followed by the
     * given header, then a newline.
     *
     * @param writer Writer to use.
     * @param header Content of the header to write.
     * @throws IOException if something goes wrong with using writer
     */
    private static void writeFastaHeader(BufferedWriter writer, String header)
            throws IOException {
        if (writer == null) {
            throw new IllegalArgumentException("null writer passed");
        } else if (header == null) {
            throw new IllegalArgumentException("null header passed");
        }
        writer.write(">" + header + "\n");
    }

    /**
     * Writes the given content (sequence data) with the given writer. Note
     * that this method will append newlines as determined by the variable
     * REMOVE_NEWLINES.
     *
     * @param writer  Writer to use.
     * @param content Content of the sequence to write.
     * @throws IOException if something goes wrong with using writer
     */
    private static void writeFastaContent(BufferedWriter writer,
                                          String content) throws IOException {
        if (writer == null) {
            throw new IllegalArgumentException("null writer passed");
        } else if (content == null) {
            throw new IllegalArgumentException("null header passed");
        }
        if (REMOVE_GAPS) {
            content = content.replace("-", "");
        }
        writer.write(content);
        if (!REMOVE_NEWLINES) {
            writer.write("\n");
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
    private static boolean isMergeable(MAFAlignmentBlock first,
                                       MAFAlignmentBlock second,
                                       String species) {
        MAFAlignmentBlock.Sequence firstSeq = first.sequences.get(species);
        MAFAlignmentBlock.Sequence secondSeq = second.sequences.get(species);

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
                MAFAlignmentBlock.Sequence.GapType.C)) {
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
