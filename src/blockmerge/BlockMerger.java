package blockmerge;

import util.*;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class BlockMerger {

    /**
     * Options for the block merging
     */
    private static Options options;

    /**
     * Remove "-" characters in a FASTA-format file output.
     */
    private static final boolean REMOVE_GAPS = true;
    /**
     * The reference sequence to use.
     */
    private static final String REF_SPECIES = "hg38";

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
    private static int MAX_OUTPUT_LENGTH_DEFAULT = 3000;
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
     * The custom name part of the output files that we will write.
     */
    private static String outName;
    /**
     * The directory to output files into.
     */
    private static String outDir;
    private static final String outDir_DEFAULT = "";
    /**
     * The file output format to use.
     */
    private static FileType outType;
    private static final FileType outType_DEFAULT = FileType.FASTA;

    /**
     * Exists to track the output of isMergeable over the course of a program
     * run.
     */
    private static SimpleBarChart disjointReasonsStats;
    private static SimpleBarChart disjointSpeciesStats;

    public static void main(String[] args) throws IOException {
        // handle command-line argument processing :)
        // add all the arguments we need
        BlockMerger.options = new Options();
        {
            BlockMerger.options.addOption(
                    Option.builder("t")
                            .longOpt("mergeType")
                            .hasArg()
                            .desc("" +
//                                          "'bases' (fixed number of bases
//                                          per merge), " +
                                          "'blocks' (fixed number of blocks " +
                                          "per merge) or 'fillblocks' " +
                                          "(contiguous blocks merged to " +
                                          "within a range of merge lenghs) to" +
                                          " determine the type of merge used " +
                                          "(by default BlockMerger uses "
                                          + "block-based merging).")
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("numBlocksPerOutput")
                            .hasArg()
                            .desc("Number of source blocks to include in each" +
                                          "merged block (if using " +
                                          "block-counting " +
                                          "merging). Defaults to " + NUM_BLOCKS_PER_OUTPUT_DEFAULT)
                            .build());
//            BlockMerger.options.addOption(
//                    Option.builder()
//                            .longOpt("numBasesPerOutput")
//                            .hasArg()
//                            .desc("Number of bases to include " +
//                                          "in each "
//                                          + "output block (if using "
//                                          + "base-counting merging). " +
//                                          "Defaults to " +
//                                          NUM_BASES_PER_OUTPUT_DEFAULT)
//                            .build());
//            BlockMerger.options.addOption(
//                    Option.builder()
//                            .longOpt("numBasesIncrement")
//                            .hasArg()
//                            .desc("Number of bases to " +
//                                          "increment between "
//                                          + "each output merged block (if" +
//                                          " using base-counting merging)." +
//                                          " Defaults to " +
//                                          NUM_BASES_INCREMENT_DEFAULT)
//                            .build());
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
                                          "block to be output. " +
                                          "Defaults to " + MAX_OUTPUT_LENGTH_DEFAULT)
                            .build());
            BlockMerger.options.addOption(
                    Option.builder()
                            .longOpt("minOutputLength")
                            .hasArg()
                            .desc("Minimum length of the reference genome " +
                                          "section in a resulting block that " +
                                          "will still allow the resulting " +
                                          "block to be output. " +
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
                            .desc("Output file prefix")
                            .required()
                            .build());
            BlockMerger.options.addOption(
                    Option.builder("ot")
                            .longOpt("outType")
                            .hasArg()
                            .desc("'FASTA' to output FASTA file, " +
                                          "'MAF' to output MAF file. " +
                                          "Defaults to FASTA format.")
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the merge type
            if (line.hasOption("t")) {
                String value = line.getOptionValue("t");
                if (value.toLowerCase().equals("bases")) {
                    BlockMerger.MERGE_TYPE = MergeType.BASES;
                } else if (value.toLowerCase().equals("fillblocks")) {
                    BlockMerger.MERGE_TYPE = MergeType.FILLBLOCKS;
                } else { //if (value.toLowerCase().equals("blocks"))
                    BlockMerger.MERGE_TYPE = MergeType.BLOCKS;
                }
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
//            // set the numBasesPerOutput
//            if (line.hasOption("numBasesPerOutput")) {
//                BlockMerger.NUM_BASES_PER_OUTPUT = Integer.valueOf(
//                        line.getOptionValue("numBasesPerOutput"));
//            } else {
//                BlockMerger.NUM_BASES_PER_OUTPUT =
//                NUM_BASES_PER_OUTPUT_DEFAULT;
//            }
//            // set the numBasesIncrement
//            if (line.hasOption("numBasesIncrement")) {
//                BlockMerger.NUM_BASES_INCREMENT = Integer.valueOf(
//                        line.getOptionValue("numBasesIncrement"));
//            } else {
//                BlockMerger.NUM_BASES_INCREMENT = NUM_BASES_INCREMENT_DEFAULT;
//            }
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
            // set the out type
            if (line.hasOption("ot")) {
                String value = line.getOptionValue("ot");
                if (value.toLowerCase().equals("fasta")) {
                    BlockMerger.outType = FileType.FASTA;
                } else { // if (value.toLowerCase().equals("maf")) {
                    BlockMerger.outType = FileType.MAF;
                }
            } else {
                BlockMerger.outType = BlockMerger.outType_DEFAULT;
            }
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
            // TODO write base-based merging...???
            throw new RuntimeException("not implemented");
        } else if (MERGE_TYPE == MergeType.BLOCKS || MERGE_TYPE == MergeType.FILLBLOCKS) {
            // Case where merge style is any of the two below:
            // Block-based merging (block-shift based, fixed number of blocks)
            // Fillblocks merging (block-shift based, but size threshold)

            // Create overall statistics trackers
            // Track gap content (N bases) in each merged block
            Path gapStatsFile = Paths.get(BlockMerger.outDir,
                                          "graph_merge_gapStatistics" +
                                                  ".png");
            SimpleScatterPlot gapStats = new SimpleScatterPlot(
                    gapStatsFile,
                    null,
                    "Merged sequence length (all species)",
                    "% gap length");
            // Track # of blocks used in each merge
            Path numBlocksFile = Paths.get(BlockMerger.outDir,
                                           "graph_merge_numBlocks" +
                                                   ".png");
            SimpleNumberHistogram numBlocksStats = new SimpleNumberHistogram(
                    numBlocksFile,
                    "",
                    "Number of blocks in mergeblock",
                    "% of mergeblocks",
                    20);
            // Track length of ref sequence involved
            Path refLengthFile = Paths.get(BlockMerger.outDir,
                                           "graph_merge_refLength" +
                                                   ".png");
            SimpleNumberHistogram refLengthStats = new SimpleNumberHistogram(
                    refLengthFile,
                    "",
                    "Length of reference sequence",
                    "% of merged results",
                    50);
            // Track causes of disjoint
            Path disjointReasonsFile = Paths.get(BlockMerger.outDir,
                                                 "graph_merge_disjointReasons" +
                                                         ".png");
            disjointReasonsStats = new SimpleBarChart(
                    disjointReasonsFile,
                    "",
                    "# Incidences");
            // Track species that are disjoint
            Path disjointSpeciesFile = Paths.get(BlockMerger.outDir,
                                                 "graph_merge_disjointSpecies" +
                                                         ".png");
            disjointSpeciesStats = new SimpleBarChart(
                    disjointSpeciesFile,
                    "",
                    "# Incidences");
            // Track block-combination limits
            Path mergeFailsFile = Paths.get(BlockMerger.outDir,
                                            "graph_merge_mergefails" +
                                                    ".png");
            SimpleBarChart mergeFailsStats = new SimpleBarChart(
                    mergeFailsFile,
                    "",
                    "# Incidences");

            LinkedList<MAFAlignmentBlock> toMerge = new LinkedList<>();
            BigInteger i = BigInteger.ONE;
            while ((MERGE_TYPE == MergeType.BLOCKS) ?
                    reader.hasNext() : // the condition if doing n-block merge
                    i.equals(BigInteger.ONE) // condition if doing fillblock
                            || toMerge.size() != 0) {
                // BLOCK-based merging rules:
                if (MERGE_TYPE == MergeType.BLOCKS) {
                    toMerge.add(reader.next());
                }
                if ((MERGE_TYPE == MergeType.BLOCKS) &&
                        toMerge.size() < NUM_BLOCKS_PER_OUTPUT) {
                    // abort if we have too few blocks merged
                    continue;
                }
                // keep the set of working blocks small
                // (we only write n blocks at a time)
                if ((MERGE_TYPE == MergeType.BLOCKS) &&
                        toMerge.size() > NUM_BLOCKS_PER_OUTPUT) {
                    toMerge.remove();
                }
                // FILLBLOCKS-based merging rules:
                BigInteger refSeqLength =
                        (toMerge.size() == 0) ? BigInteger.ZERO :
                                toMerge.getLast().sequences.get(REF_SPECIES)
                                        .start
                                        .add(
                                                toMerge.getLast().sequences.get(
                                                        REF_SPECIES).size)
                                        .subtract(
                                                toMerge.getFirst().sequences.get(
                                                        REF_SPECIES).start);
                // check upper bound
                if ((MERGE_TYPE == MergeType.FILLBLOCKS) &&
                        refSeqLength.compareTo(
                                BigInteger.valueOf(MAX_OUTPUT_LENGTH)) > 0) {
                    // need fewer blocks to reduce to maximum
                    mergeFailsStats.addValue("fblocks_above_maxlength");
                    toMerge.removeFirst();
                    continue;
                }
                // check lower bound
                if ((MERGE_TYPE == MergeType.FILLBLOCKS) &&
                        refSeqLength.compareTo(
                                BigInteger.valueOf(MIN_OUTPUT_LENGTH)) < 0) {
                    // need more blocks to hit minimum
                    mergeFailsStats.addValue("fblocks_below_minlength");
                    if (reader.hasNext()) {
                        toMerge.add(reader.next());
                    } else {
                        // there are no more blocks to add and no point removing
                        break;
                    }
                    continue;
                }

                // find the set of mergeable species
                List<String> speciesToMerge = new LinkedList<>();
                for (String species : toMerge.getFirst().getSpecies()) {
                    // initialize iterators through the blocks
                    Iterator<MAFAlignmentBlock> secondIt = toMerge.iterator();
                    secondIt.next();
                    Iterator<MAFAlignmentBlock> firstIt = toMerge.iterator();
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

                // Standard Filter:
                // checks the # of merged species threshold
                if (speciesToMerge.size() < MIN_NUM_SPECIES) {
                    mergeFailsStats.addValue("too_few_species_mergeable");
                    continue;
                }

                // Standard Filter:
                // require that the merged species includes the stated reference
                if (!speciesToMerge.contains(REF_SPECIES)) {
                    mergeFailsStats.addValue("ref_species_not_mergeable");
                    continue;
                }

                // Standard Filter:
                // check output against sequence length bounds
                MAFAlignmentBlock.Sequence firstCheckSeqLen =
                        toMerge.getFirst().sequences.get(
                                REF_SPECIES);
                MAFAlignmentBlock.Sequence lastCheckSeqLen =
                        toMerge.getLast().sequences.get(
                                REF_SPECIES);
                BigInteger humanSeqLen = lastCheckSeqLen.start.add(
                        lastCheckSeqLen.size).subtract(firstCheckSeqLen.start);
                // check lower bound
                if (humanSeqLen.compareTo(
                        BigInteger.valueOf(MIN_OUTPUT_LENGTH)) < 0) {
                    mergeFailsStats.addValue("merged_below_minlength");
                    continue;
                }
                // check upper bound
                if (humanSeqLen.compareTo(
                        BigInteger.valueOf(MAX_OUTPUT_LENGTH)) > 0) {
                    mergeFailsStats.addValue("merged_above_maxlength");
                    continue;
                }

                refLengthStats.addValue(new BigDecimal(refSeqLength));

                // Write to disk!:
                // create output file
                Path file = Paths.get(BlockMerger.outDir,
                                      BlockMerger.outName
                                              + "_"
                                              + i.toString()
                                              + "."
                                              + BlockMerger.outType.toString().toLowerCase());
                Files.deleteIfExists(file);
                Files.createFile(file);
                BufferedWriter writer = Files.newBufferedWriter(file);

                // build the sequence lines for each individual species first
                Map<String, StringBuilder> headers = new HashMap<>();
                Map<String, StringBuilder> sequences = new HashMap<>();
                Map<String, BigInteger> seqLengths = new HashMap<>();
                Map<String, BigInteger> gapLengths = new HashMap<>();
                // Generate the headers
                for (String species : speciesToMerge) {
                    // track info for gapStats later
                    BigInteger seqLength = BigInteger.ZERO;

                    // write the species and originating chromosome
                    MAFAlignmentBlock.Sequence first =
                            toMerge.getFirst().sequences.get(
                                    species);
                    MAFAlignmentBlock.Sequence last =
                            toMerge.getLast().sequences.get(
                                    species);

                    seqLength = last.start.add(last.size).subtract(first.start);

                    // if the merged sequence has a size of 0
                    if (seqLength.equals(BigInteger.ZERO)) {
                        headers.put(species, null);
                        seqLengths.put(species, seqLength);
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
                    Iterator<MAFAlignmentBlock> it = toMerge.iterator();
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
                    // we have the species header at this point
                    headers.put(species, speciesHeader);
                    seqLengths.put(species, seqLength);
                }

                // build the individual (aligned) full sequences
                // temporary variable to keep track of how merging is going
                Map<String, Boolean> lastBlockWasGap = new HashMap<>();
                for (String species : speciesToMerge) {
                    // initialize the sequence-related maps
                    sequences.put(species, new StringBuilder());
                    lastBlockWasGap.put(species, Boolean.FALSE);
                    gapLengths.put(species, BigInteger.ZERO);
                }
                for (int blockNum = 0; blockNum < toMerge.size(); blockNum++) {
                    MAFAlignmentBlock block = toMerge.get(blockNum);

                    // write the interval section between blocks if appropriate
                    if (blockNum != 0) {
                        // collect the list of intervals need to be written
                        Set<String> areGaps = new HashSet<>();
                        for (String species : speciesToMerge) {
                            if (!lastBlockWasGap.get(species)) {
                                // we will need to write some Ns
                                areGaps.add(species);
                            }
                        }
                        // write the intervals
                        for (String species : areGaps) {
                            MAFAlignmentBlock.Sequence seq =
                                    block.sequences.get(species);
                            for (String currentlyWriting : speciesToMerge) {
                                if (species.equals(currentlyWriting)) {
                                    // write the Ns
                                    sequences.get(currentlyWriting)
                                            .append(StringManip.repeat("N",
                                                                       seq.left.length));
                                    gapLengths.put(
                                            currentlyWriting,
                                            gapLengths.get(currentlyWriting)
                                                    .add(seq.left.length));
                                } else {
                                    // write the -s
                                    sequences.get(currentlyWriting)
                                            .append(StringManip.repeat("-",
                                                                       seq.left.length));
                                }
                            }
                        }
                    }

                    // now write each individual block contents
                    // we do the block-scale inserts the same way as
                    //     we do mid-block inserts
                    // list all gap/insert parts first!
                    // also update the "last block gap status"
                    Set<String> areGaps = new HashSet<>();
                    BigInteger alignmentLength = BigInteger.ZERO;
                    for (String species : speciesToMerge) {
                        if (block.sequences.get(species).isGap) {
                            areGaps.add(species);
                            lastBlockWasGap.put(species, Boolean.TRUE);
                        } else {
                            lastBlockWasGap.put(species, Boolean.FALSE);
                            alignmentLength = BigInteger.valueOf(
                                    block.sequences.get(species).contents
                                            .length());
                        }
                    }
                    // write the gap sections first
                    for (String species : areGaps) {
                        MAFAlignmentBlock.Sequence seq =
                                block.sequences.get(species);
                        for (String currentlyWriting : speciesToMerge) {
                            if (species.equals(currentlyWriting)) {
                                // Write the Ns for the given sequence
                                sequences.get(currentlyWriting)
                                        .append(StringManip.repeat("N",
                                                                   seq.size));
                                // and increment the gap length for seq
                                gapLengths.put(
                                        currentlyWriting,
                                        gapLengths.get(currentlyWriting)
                                                .add(seq.size));
                            } else {
                                // Write a series of unaligned section (-)
                                sequences.get(currentlyWriting)
                                        .append(StringManip.repeat("-",
                                                                   seq.size));
                                // we don't need to increment the gap length lol
                            }
                        }
                    }
                    // write the actual alignment section now
                    for (String species : speciesToMerge) {
                        if (areGaps.contains(species)) {
                            // write a series of unaligned section (-)
                            sequences.get(species)
                                    .append(StringManip.repeat("-",
                                                               alignmentLength));
                        } else {
                            // write the actual alignment contents
                            sequences.get(species)
                                    .append(block.sequences
                                                    .get(species).contents);
                        }
                    }
                }

                // Now write the data!!!
                // write an overall file header if necessary
                if (BlockMerger.outType == FileType.MAF) {
                    // MAF header
                    writer.write("##maf version=1 program=blockmerger.jar\n");
                    writer.write("a \n");
                }
                // write file lines for each species to disk
                for (String species : speciesToMerge) {
                    // Initialize to track info for gapStats
                    BigInteger seqLength = seqLengths.get(species);

                    // if the merged sequence has a size of 0
                    if (seqLength.equals(BigInteger.ZERO)) {
                        // don't even print it
                        // skip this line of output
                        continue;
                    }

                    // build the sequence name
                    StringBuilder speciesHeader = headers.get(species);

                    // build the contents of each block for each species
                    StringBuilder sequence = sequences.get(species);

                    // Write the file proper
                    if (BlockMerger.outType == FileType.FASTA) {
                        // Write a FASTA entry
                        writer.write(">" + speciesHeader.toString() + "\n");
                        String content = sequence.toString();
                        if (REMOVE_GAPS) {
                            content = content.replace("-", "");
                        }
                        writer.write(content);
                        writer.write("\n");
                    } else if (BlockMerger.outType == FileType.MAF) {
                        // A source of canonical seq information...
                        MAFAlignmentBlock.Sequence first =
                                toMerge.getFirst().sequences.get(
                                        species);

                        // write a MAF entry
                        // standard line begin
                        writer.write("s ");
                        // For MAF we ignore potential scaffold merges...
                        // Species and chromosome
                        String print = species + "." + first.section;
                        writer.write(print + StringManip.repeat(" ",
                                                                BigInteger.valueOf(
                                                                        24 - print.length()))
                                             + " ");
                        // start position
                        print = "" + first.start;
                        writer.write(print + StringManip.repeat(" ",
                                                                BigInteger.valueOf(
                                                                        12 - print.length()))
                                             + " ");
                        // number of actual bases used (excluding dash)
                        String degapped = sequence.toString()
                                .replace("-", "");
                        print = "" + degapped.length();
                        writer.write(print + StringManip.repeat(" ",
                                                                BigInteger.valueOf(
                                                                        6 - print.length()))
                                             + " ");
                        // strand used
                        writer.write((first.strand ? "+" : "-") + " ");
                        // Still ignore scaffold merges
                        // Total length of source sequence
                        print = "" + first.srcSize;
                        writer.write(print + StringManip.repeat(" ",
                                                                BigInteger.valueOf(
                                                                        14 - print.length()))
                                             + " ");
                        // The alignment sequence itself
                        writer.write(sequence.toString());
                        writer.write("\n");
                    }
                }
                // Add a closing newline
                if (BlockMerger.outType == FileType.MAF) {
                    writer.write("\n");
                }
                // Plot to gapStats
                for (String species : speciesToMerge) {
                    // write the sequence length vs gap percentage
                    // (one point per species-sequence, excludes 0-length seq)
                    if (seqLengths.get(species).equals(BigInteger.ZERO)) {
                        continue;
                    }
                    BigDecimal seqLength = new BigDecimal(
                            seqLengths.get(species));
                    BigDecimal gapLength = new BigDecimal(
                            gapLengths.get(species));
                    gapStats.addValue(seqLength,
                                      (gapLength)
                                              .multiply(BigDecimal.valueOf(100))
                                              .divide(seqLength,
                                                      RoundingMode.HALF_EVEN));
                }
                // write the mergeblocksize statistic
                // (one point per output file)
                numBlocksStats.addValue(new BigDecimal(toMerge.size()));
                // note that we have created a file
                System.out.println("Wrote " + file.toString());

                // increment the file number
                i = i.add(BigInteger.ONE);

                writer.close();

                // for fillblocks: shift forwards
                if ((MERGE_TYPE == MergeType.FILLBLOCKS)) {
                    if (reader.hasNext()) {
                        // absorb the next block if possible
                        toMerge.add(reader.next());
                    } else {
                        // otherwise leave the first one behind
                        toMerge.removeFirst();
                    }
                }
            }

            // Output overall statistics graphics
            gapStats.write();
            numBlocksStats.write();
            disjointReasonsStats.write();
            disjointSpeciesStats.write();
            mergeFailsStats.write();
            refLengthStats.write();
            System.out.println("Wrote statistics information");

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
            disjointReasonsStats.addValue("missing_from_block");
            disjointSpeciesStats.addValue(species);
            return false;
        }

        // test if they are in different non-scaffold chromosomes:
        String firstSec = firstSeq.section;
        String secondSec = secondSeq.section;
        if (!firstSec.contains("scaffold") && !secondSec.contains("scaffold")) {
            if (!firstSec.equals(secondSec)) {
                // If neither chromosome is a scaffold and they are in
                // different chromosomes, then they are not mergeable
                disjointReasonsStats.addValue("diff_chromosome");
                disjointSpeciesStats.addValue(species);
                return false;
            }
        }

        // test if they are on different strands:
        if (firstSeq.strand != secondSeq.strand) {
            // If they are on different strands, they are not mergeable
            disjointReasonsStats.addValue("diff_strand");
            disjointSpeciesStats.addValue(species);
            return false;
        }

        // test if they are at incomparable locations / overlapping:
        if (firstSeq.start.add(firstSeq.size).compareTo(secondSeq.start) > 0) {
            // if the ending of the first block is after
            // the beginning of the second block
            // Note that this works for reverse strand too,
            // by convention of the alignment
            disjointReasonsStats.addValue("reversed_order");
            disjointSpeciesStats.addValue(species);
            return false;
            // NOTE I don't think this will happen with my input MAFS
            // but it is probably still good to check for
        }

        // test if the gap between them is too long:
        if (firstSeq.isGap && firstSeq.gapType.equals(
                MAFAlignmentBlock.Sequence.GapType.C)) {
            // if first is a gap sequence (with no bases in the section)
            // and it is too long to be included
            disjointReasonsStats.addValue("large_gap");
            disjointSpeciesStats.addValue(species);
            if (firstSeq.size.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                return false;
            }
        }
        if (secondSeq.isGap) {
            // if second is a gap sequence and it is too long to be included
            disjointReasonsStats.addValue("large_gap");
            disjointSpeciesStats.addValue(species);
            if (secondSeq.size.compareTo(
                    BigInteger.valueOf(GAP_THRESHOLD)) > 0) {
                return false;
            }
        }
        if (!firstSeq.isGap && !secondSeq.isGap) {
            // if neither is a gap, but the intervening "gap" is too long
            disjointReasonsStats.addValue("large_gap");
            disjointSpeciesStats.addValue(species);
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
        BASES, BLOCKS, FILLBLOCKS
    }

}
