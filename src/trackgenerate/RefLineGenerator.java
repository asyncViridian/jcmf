package trackgenerate;

import org.apache.commons.cli.*;
import util.MAFReader;
import util.SimpleHistogram;

import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class RefLineGenerator {
    private static Options options;
    private static String srcMaf;
    private static String outDir;
    private static Boolean writeHeader;

    public static void main(String[] args) throws IOException {
        // handle command-line argument processing :)
        // add all the arguments we need
        RefLineGenerator.options = new Options();
        {
            RefLineGenerator.options.addOption(
                    Option.builder("s")
                            .longOpt("srcMaf")
                            .hasArg()
                            .desc("Input MAF to make a reference for")
                            .required()
                            .build());
            RefLineGenerator.options.addOption(
                    Option.builder("o")
                            .longOpt("outDir")
                            .hasArg()
                            .desc("Output directory")
                            .required()
                            .build());
            RefLineGenerator.options.addOption(
                    Option.builder("h")
                            .longOpt("writeHeader")
                            .desc("Include the header line")
                            .build());
            // TODO add options??? do I even need more options???
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            RefLineGenerator.srcMaf = line.getOptionValue("s");
            // set the output filename
            RefLineGenerator.outDir = line.getOptionValue("o");
            // set the header line boolean
            RefLineGenerator.writeHeader = line.hasOption("h");
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "RefLineGenerator [options] -s <filename> -o <filename>",
                    options);
            return;
        }

        // generate the output files
        // the multiple-block motifs file
        Path output = Paths.get(RefLineGenerator.outDir, "ref.bed");
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);
        // histograms
        Path blockLengthStatsFile = Paths.get(RefLineGenerator.outDir,
                                              "graph_refBlockLengthStats" +
                                                      ".png");
        SimpleHistogram blockLengthStats = new SimpleHistogram(
                blockLengthStatsFile,
                "",
                "Block length",
                "% of blocks",
                20);

        // write the header lines
        if (writeHeader) {
            // TODO set a more useful name field?
            writer.write("track"
                                 + " name=" + "referenceBlocksTrack"
                                 + " type=" + "bedDetail"
                                 + "\n");
        }

        // read the given MAF
        MAFReader reader = new MAFReader("", RefLineGenerator.srcMaf);
        while (reader.hasNext()) {
            // write block length statistics for the next block
            blockLengthStats.addValue(
                    new BigDecimal(reader.next().sequences.get("hg38").size));
        }

        // output BED tracks
        writer.close();
        System.out.println("Wrote BED files");

        // Output histograms
        blockLengthStats.write();
        System.out.println("Wrote statistics graphs");
    }
}
