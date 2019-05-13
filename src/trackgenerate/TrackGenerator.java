package trackgenerate;

import javafx.util.Pair;
import org.apache.commons.cli.*;
import util.StockholmAlignmentBlock;

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
    private static String outFile;
    private static String outFileSingle;

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
                            .longOpt("outFile")
                            .hasArg()
                            .desc("Output BED filename (including extension)")
                            .required()
                            .build());
            TrackGenerator.options.addOption(
                    Option.builder("os")
                            .longOpt("outFileSingle")
                            .hasArg()
                            .desc("Output BED filename for single-block " +
                                          "motifs(including extension)")
                            .required()
                            .build());
            // TODO add options??? do I even need more options???
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            TrackGenerator.srcDir = line.getOptionValue("s");
            // set the output filename
            TrackGenerator.outFile = line.getOptionValue("o");
            // set the single-block output filename
            TrackGenerator.outFileSingle = line.getOptionValue("os");
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "TrackGenerator [options] -s <dirname> -o <filename>",
                    options);
            return;
        }

        // generate the output file
        Path output = Paths.get(TrackGenerator.outFile);
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);
        Path outputSingle = Paths.get(TrackGenerator.outFileSingle);
        Files.deleteIfExists(outputSingle);
        Files.createFile(outputSingle);


        // TODO write the header lines.
        // TODO write header lines in bedDetail format

        // read the given directory...
        Path source = Paths.get(TrackGenerator.srcDir);
        for (File f : source.toFile().listFiles()) {
            // for each file in the given directory,
            // check if it should be added to the BED and add it if necessary
            StockholmAlignmentBlock block =
                    StockholmAlignmentBlock.constructFromScore(f);

            // something was malformatted?
            if (block == null) {
                continue;
            }

            // check if block contains human DNA
            if (!block.containsSpecies("hg38")) {
                continue;
            }

            // score filter
            if (block.pairScore.compareTo(BigDecimal.valueOf(40L)) < 0) {
                continue;
            }

            // TODO set a better score filter and value guideline
            BigDecimal max = BigDecimal.valueOf(125L);

            // TODO check for if it is in a single block
            // add it to the BED
            Pair<BigInteger, BigInteger> interval
                    = block.getInterval("hg38");
            String chr = block.getChromosome("hg38");
            String name = f.getName();
            BigDecimal pairScore = block.pairScore;
            writer.write(chr + "\t"
                                 + interval.getKey() + "\t"
                                 + interval.getValue() + "\t"
                                 + name + "\t"
                                 + pairScore.divide(max,
                                                    RoundingMode.HALF_EVEN).multiply(
                    BigDecimal.valueOf(1000)) + "\n");
        }
        writer.close();
    }
}
