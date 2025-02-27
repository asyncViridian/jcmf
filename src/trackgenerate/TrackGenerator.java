package trackgenerate;

import org.apache.commons.cli.*;
import org.apache.commons.lang3.tuple.Pair;
import util.SimpleNumberHistogram;
import util.SimpleBarChart;
import util.SimpleScatterPlot;
import util.ScoredStockholmAlignmentBlock;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class TrackGenerator {
    private static Options options;
    private static String srcDir;
    private static String outputDir;
    private static Boolean writeHeader;

    public static void main(String[] args) throws IOException {
        // handle command-line argument processing :)
        // add all the arguments we need
        TrackGenerator.options = new Options();
        {
            TrackGenerator.options.addOption(
                    Option.builder("s")
                            .longOpt("srcDir")
                            .hasArg()
                            .desc("Input score files directory")
                            .required()
                            .build());
            TrackGenerator.options.addOption(
                    Option.builder("o")
                            .longOpt("outputDir")
                            .hasArg()
                            .desc("Output directory")
                            .required()
                            .build());
            TrackGenerator.options.addOption(
                    Option.builder("h")
                            .longOpt("writeHeader")
                            .desc("Include the header line")
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            TrackGenerator.srcDir = line.getOptionValue("s");
            // set the output filename
            TrackGenerator.outputDir = line.getOptionValue("o");
            // set the header line boolean
            TrackGenerator.writeHeader = line.hasOption("h");
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "TrackGenerator [options] -s <dirname> -o <filename>",
                    options);
            return;
        }

        // generate the output files
        // the multiple-block motifs file
        Path output = Paths.get(TrackGenerator.outputDir, "multi.bed");
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);
        // the single-block motifs file
        Path outputSingle = Paths.get(TrackGenerator.outputDir,
                                      "single.bed");
        Files.deleteIfExists(outputSingle);
        Files.createFile(outputSingle);
        BufferedWriter writerSingle = Files.newBufferedWriter(outputSingle);
        // histograms
        Path rnaScorePreFile = Paths.get(TrackGenerator.outputDir,
                                         "graph_prefilter_rnaScoreStats" +
                                                 ".png");
        SimpleNumberHistogram rnaScorePreStats = new SimpleNumberHistogram(
                rnaScorePreFile,
                "",
                "RNA posterior score",
                "% of motifs",
                20);
        Path rnaScorePostFile = Paths.get(TrackGenerator.outputDir,
                                          "graph_postfilter_rnaScoreStats" +
                                                  ".png");
        SimpleNumberHistogram rnaScorePostStats = new SimpleNumberHistogram(
                rnaScorePostFile,
                "",
                "RNA posterior score",
                "% of motifs",
                20);
        Path pairScorePreFile = Paths.get(TrackGenerator.outputDir,
                                          "graph_prefilter_pairScoreStats" +
                                                  ".png");
        SimpleNumberHistogram pairScorePreStats = new SimpleNumberHistogram(
                pairScorePreFile,
                "",
                "Pair posterior score",
                "% of motifs",
                20);
        Path pairScorePostFile = Paths.get(TrackGenerator.outputDir,
                                           "graph_postfilter_pairScoreStats" +
                                                   ".png");
        SimpleNumberHistogram pairScorePostStats = new SimpleNumberHistogram(
                pairScorePostFile,
                "",
                "Pair posterior score",
                "% of motifs",
                20);
        Path blockSpanPreFile = Paths.get(TrackGenerator.outputDir,
                                          "graph_prefilter_blockSpanStats" +
                                                  ".png");
        SimpleNumberHistogram blockSpanPreStats = new SimpleNumberHistogram(
                blockSpanPreFile,
                "",
                "Number of blocks spanned",
                "% of motifs",
                10);
        Path blockSpanPostFile = Paths.get(TrackGenerator.outputDir,
                                           "graph_postfilter_blockSpanStats" +
                                                   ".png");
        SimpleNumberHistogram blockSpanPostStats = new SimpleNumberHistogram(
                blockSpanPostFile,
                "",
                "Number of blocks spanned",
                "% of motifs",
                10);
        Path sequenceLengthPreFile = Paths.get(TrackGenerator.outputDir,
                                               "graph_prefilter_sequenceLengthStats" +
                                                       ".png");
        SimpleNumberHistogram sequenceLengthPreStats =
                new SimpleNumberHistogram(
                        sequenceLengthPreFile,
                        "",
                        "Length of merged sequence",
                        "% of motifs",
                        20);
        Path sequenceLengthPostFile = Paths.get(TrackGenerator.outputDir,
                                                "graph_postfilter_sequenceLengthStats" +
                                                        ".png");
        SimpleNumberHistogram sequenceLengthPostStats =
                new SimpleNumberHistogram(
                        sequenceLengthPostFile,
                        "",
                        "Length of merged sequence",
                        "% of motifs",
                        20);
        Path motifSpeciesPreFile = Paths.get(TrackGenerator.outputDir,
                                             "graph_prefilter_motifSpeciesStats" +
                                                     ".png");
        SimpleBarChart motifSpeciesPreStats =
                new SimpleBarChart(
                        motifSpeciesPreFile,
                        "",
                        "Number of Motifs"
                );
        Path motifSpeciesPostFile = Paths.get(TrackGenerator.outputDir,
                                              "graph_postfilter_motifSpeciesStats" +
                                                      ".png");
        SimpleBarChart motifSpeciesPostStats =
                new SimpleBarChart(
                        motifSpeciesPostFile,
                        "",
                        "Number of Motifs"
                );
        Path pairVsRnaScorePreFile = Paths.get(TrackGenerator.outputDir,
                                               "graph_prefilter_pairVsRnaStats" +
                                                       ".png");
        SimpleScatterPlot pairVsRnaScorePreStats =
                new SimpleScatterPlot(
                        pairVsRnaScorePreFile,
                        null,
                        "RNA Posterior Score",
                        "Pair Posterior Score"
                );
        Path pairVsRnaScorePostFile = Paths.get(TrackGenerator.outputDir,
                                                "graph_postfilter_pairVsRnaStats" +
                                                        ".png");
        SimpleScatterPlot pairVsRnaScorePostStats =
                new SimpleScatterPlot(
                        pairVsRnaScorePostFile,
                        null,
                        "RNA Posterior Score",
                        "Pair Posterior Score"
                );

        // write the header lines
        if (writeHeader) {
            // TODO set a more useful name field?
            writer.write("track"
                                 + " name=" + "motifsFoundMultiBlock"
                                 + " type=" + "bedDetail"
                                 + "\n");
            // TODO set a more useful name field?
            writerSingle.write("track"
                                       + " name=" + "motifsFoundSingleBlock"
                                       + " type=" + "bedDetail"
                                       + "\n");
        }

        // read the given directory...
        Path source = Paths.get(TrackGenerator.srcDir);
        for (File f : source.toFile().listFiles()) {
            // for each file in the given directory,
            // check if it should be added to the BED and add it if necessary
            ScoredStockholmAlignmentBlock block =
                    ScoredStockholmAlignmentBlock.constructFromScore(f);

            // something was malformatted?
            if (block == null) {
                continue;
            }

            // check if block contains human DNA
            if (!block.containsSpecies("hg38")) {
                continue;
            }

            BigDecimal pairScore = block.pairScore;
            BigDecimal rnaScore = block.rnaScore;

            //  Get prefilter stats
            rnaScorePreStats.addValue(rnaScore);
            pairScorePreStats.addValue(pairScore);
            blockSpanPreStats.addValue(
                    BigDecimal.valueOf(block.motifInNumBlocks("hg38")));
            Pair<BigInteger, BigInteger> interval
                    = block.getInterval("hg38");
            sequenceLengthPreStats.addValue(
                    new BigDecimal(
                            interval.getRight().subtract(interval.getLeft()))
            );
            for (String species : block.getSpecies()) {
                motifSpeciesPreStats.addValue(species);
            }
            pairVsRnaScorePreStats.addValue(rnaScore, pairScore);

            // score filter
            if (rnaScore.compareTo(BigDecimal.valueOf(50L)) < 0) {
                continue;
            }

            // Begin postfilter details
            // Get postfilter stats
            rnaScorePostStats.addValue(rnaScore);
            pairScorePostStats.addValue(pairScore);
            blockSpanPostStats.addValue(
                    BigDecimal.valueOf(block.motifInNumBlocks("hg38")));
            sequenceLengthPostStats.addValue(
                    new BigDecimal(
                            interval.getRight().subtract(interval.getLeft()))
            );
            for (String species : block.getSpecies()) {
                motifSpeciesPostStats.addValue(species);
            }
            pairVsRnaScorePostStats.addValue(rnaScore, pairScore);

            // TODO set a better score filter and value guideline?
            BigDecimal max = BigDecimal.valueOf(125L);

            // construct the line that we will output for this motif
            String chr = block.getChromosome("hg38");
            String name = f.getName();
            // set cap of the computed scaled score to 1000
            BigDecimal score = rnaScore.divide(max, RoundingMode.HALF_EVEN)
                    .multiply(BigDecimal.valueOf(1000));
            if (score.compareTo(BigDecimal.valueOf(1000)) > 0) {
                score = BigDecimal.valueOf(1000);
            }
            String thingToWrite = chr + "\t"
                    + interval.getKey() + "\t"
                    + interval.getValue() + "\t"
                    + name + "\t"
                    + score.intValue() + "\n";
            // determine which BED track we want to output to
            if (block.motifInNumBlocks("hg38") == 1) {
                // add it to the single-block BED
                writerSingle.write(thingToWrite);
            } else {
                // add it to the multiblock BED
                writer.write(thingToWrite);
            }
        }
        writer.close();
        writerSingle.close();
        System.out.println("Wrote BED files");

        // Output histograms
        rnaScorePreStats.write();
        rnaScorePostStats.write();
        pairScorePreStats.write();
        pairScorePostStats.write();
        blockSpanPreStats.write();
        blockSpanPostStats.write();
        sequenceLengthPreStats.write();
        sequenceLengthPostStats.write();
        motifSpeciesPreStats.write();
        motifSpeciesPostStats.write();
        pairVsRnaScorePreStats.write();
        pairVsRnaScorePostStats.write();
        System.out.println("Wrote statistics graphs");
    }
}
