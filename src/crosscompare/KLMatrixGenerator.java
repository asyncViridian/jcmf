package crosscompare;

import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

public class KLMatrixGenerator {

    private static Options options;
    private static String srcDir;
    private static String outputFile;
    private static int numSamples;

    public static void main(String[] args) throws IOException {
        long startTime = System.nanoTime();
        // handle command-line argument processing :)
        // add all the arguments we need
        KLMatrixGenerator.options = new Options();
        {
            KLMatrixGenerator.options.addOption(
                    Option.builder("s")
                            .longOpt("srcDir")
                            .hasArg()
                            .desc("Input CM files directory. Required.")
                            .required()
                            .build());
            KLMatrixGenerator.options.addOption(
                    Option.builder("o")
                            .longOpt("outputFile")
                            .hasArg()
                            .desc("Output filename. Required.")
                            .required()
                            .build());
            KLMatrixGenerator.options.addOption(
                    Option.builder("n")
                            .hasArg()
                            .desc("The number of samples to use for " +
                                          "divergence calculation. " +
                                          "Required.")
                            .required()
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            KLMatrixGenerator.srcDir = line.getOptionValue("s");
            // set the output filename
            KLMatrixGenerator.outputFile = line.getOptionValue("o");
            // set the number of samples to use
            KLMatrixGenerator.numSamples
                    = Integer.valueOf(line.getOptionValue("n"));
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "KLMatrixGenerator -s <dirname> -o " +
                            "<filename> -n <numsamples>",
                    options);
            return;
        }

        // Get the input & output files set up
        Path source = Paths.get(KLMatrixGenerator.srcDir);
        Path output = Paths.get(KLMatrixGenerator.outputFile);
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);

        // Read all the input files into list
        Object[] files = Arrays.stream(source.toFile().listFiles())
                .filter(f -> f.getName().endsWith(".cm"))
                .toArray();
        // Get the max filename length so we can make pretty whitespace ...
        int maxFilenameLength = 0;
        for (int i = 0; i < files.length; i++) {
            int currFilenameLength = ((File) files[i]).getName().length();
            if (currFilenameLength > maxFilenameLength) {
                maxFilenameLength = currFilenameLength;
            }
        }

        String[] cmd = {"/bin/sh", "-c", "ls > hello"};
        Runtime.getRuntime().exec(cmd);

        System.out.println(
                "KLMatrixGenerator finished running in "
                        + (System.nanoTime() - startTime) / 1000000 + "ms");
    }
}
