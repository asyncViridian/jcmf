package crosscompare;

import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.UUID;
import java.util.concurrent.*;

public class KLMatrixGenerator {

    private static Options options;
    /**
     * The directory that this program reads input CM files from.
     */
    private static String srcDir;
    private static String outputFile;
    private static int numSamples;
    /**
     * The command that this program will run to get a K-L divergence score
     * between two CMs.
     */
    private static String scoreCalcCmd = "./KLDivergenceEstimate.sh";
    /**
     * The number of output scores expected from the calc cmd (separated by
     * newlines)
     */
    private static int expectedOutputs = 1;

    public static void main(String[] args)
            throws IOException, InterruptedException {
        long startTime = System.nanoTime();
        // handle command-line argument processing :)
        // add all the arguments we need
        KLMatrixGenerator.options = new Options();
        {
            KLMatrixGenerator.options.addOption(
                    Option.builder("s")
                            .longOpt("srcDir")
                            .hasArg()
                            .desc("Input CM files directory. " +
                                          "Note that this program will do a " +
                                          "matrix comparison on ALL CM files " +
                                          "within this directory. It will " +
                                          "ONLY include files that have a " +
                                          "file extension of '.cm'. Required.")
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
            KLMatrixGenerator.options.addOption(
                    Option.builder("e")
                            .hasArg()
                            .desc("The command/executable to use for finding " +
                                          "divergence score. Must take in 4 " +
                                          "args in the given order: a 'temp' " +
                                          "working directory name, a CM " +
                                          "filepath, another CM filepath, and" +
                                          " a number of samples to use. " +
                                          "Defaults to '" + scoreCalcCmd + "'")
                            .build());
            KLMatrixGenerator.options.addOption(
                    Option.builder("es")
                            .hasArg()
                            .desc("The expected number of output floats from " +
                                          "the command/executable used for " +
                                          "scoring. One output file will be " +
                                          "produced per score. Defaults to " + expectedOutputs)
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
            // set the executable
            if (line.hasOption("e")) {
                KLMatrixGenerator.scoreCalcCmd =
                        line.getOptionValue("e");
                KLMatrixGenerator.expectedOutputs =
                        Integer.valueOf(line.getOptionValue("es"));
            }
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

        // Check that the script being used exists in the same directory
        if (Files.notExists(Paths.get(scoreCalcCmd))) {
            System.err.println("Please ensure that the command " +
                                       "'" + scoreCalcCmd + "' " +
                                       "can be executed from the working " +
                                       "directory.");
            return;
        }

        // Get the input & output files set up
        Path source = Paths.get(KLMatrixGenerator.srcDir);
        Path output = Paths.get(KLMatrixGenerator.outputFile);
        Path[] outputs = new Path[expectedOutputs];
        outputs[0] = output;
        if (expectedOutputs > 1) {
            for (int i = 1; i < expectedOutputs; i++) {
                outputs[i] = Paths
                        .get(outputFile
                                     .substring(0,
                                                outputFile.lastIndexOf('.'))
                                     + "_" + i
                                     + outputFile
                                .substring(
                                        outputFile.lastIndexOf('.')));
            }
        }
        for (Path o : outputs) {
            Files.deleteIfExists(o);
            Files.createFile(o);
        }

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
        maxFilenameLength += 2;

        // Start constructing our matrix
        Future<BigDecimal[]>[][] scoreFutures =
                new Future[files.length][files.length];
        ExecutorService executor = Executors.newFixedThreadPool(8);
        for (int i = 0; i < files.length; i++) {
            String pathP = ((File) files[i]).getPath();
            for (int j = 0; j < files.length; j++) {
                String pathQ = ((File) files[j]).getPath();
                scoreFutures[i][j] = executor.submit(() -> {
                    // Create temp working directory for the estimate script
                    Path tempDir = Paths
                            .get("klmatrixtemp-" + UUID.randomUUID().toString());
                    Files.deleteIfExists(tempDir);
                    Files.createDirectory(tempDir);

                    // Run the K-L divergence score estimator script
                    // TODO make it more possible to parameterize the script???
                    String[] cmd = {scoreCalcCmd,
                            tempDir.toString(),
                            pathP,
                            pathQ,
                            "" + numSamples};
                    try {
                        // execute script and collect the output from it
                        Process p = Runtime.getRuntime().exec(cmd);
                        BufferedReader reader = new BufferedReader(
                                new InputStreamReader(p.getInputStream()));
                        StringBuilder pOutput = new StringBuilder();
                        String line;
                        while ((line = reader.readLine()) != null) {
                            pOutput.append(line);
                            pOutput.append("\n");
                        }
                        reader.close();
                        String[] pOutLines = pOutput.toString().split("\n");
                        System.out.println("line "+pOutput);
                        System.out.println("pOut "+Arrays.toString(pOutLines));
                        BigDecimal[] result = new BigDecimal[expectedOutputs];
                        System.out.println("bDec "+Arrays.toString(result));
                        for (int res = 0; res < result.length; res++) {
                            result[res] = BigDecimal.valueOf(
                                    Double.valueOf(pOutLines[res])
                            );
                        }
                        return result;
                    } catch (IOException e) {
                        e.printStackTrace();
                        // something went wrong!!!
                        // note that K-L divergence is never negative so this
                        // is a clear sign of Bad Things Happened(tm)
                        return new BigDecimal[]{BigDecimal.valueOf(-1)};
                    } finally {
                        // Files.deleteIfExists(tempDir);
                    }
                });
            }
        }
        // Collect the results of our matrix
        BigDecimal[][][] scores = new BigDecimal
                [expectedOutputs]
                [files.length]
                [files.length];
        // set them all to null for now
        for (BigDecimal[][] a : scores) {
            for (BigDecimal[] b : a) {
                for (int i = 0; i < b.length; i++) {
                    b[i] = null;
                }
            }
        }
        int done = 0;
        while (done != (files.length * files.length)) {
            for (int i = 0; i < files.length; i++) {
                for (int j = 0; j < files.length; j++) {
                    if (scores[0][i][j] == null) {
                        if (scoreFutures[i][j].isDone()) {
                            try {
                                BigDecimal[] results = scoreFutures[i][j].get();
                                for (int o = 0; o < expectedOutputs; o++) {
                                    scores[o][i][j] = results[o];
                                }
                            } catch (InterruptedException | ExecutionException e) {
                                // this shouldn't happen if its done
                                // but we will mark results with a -1 anyway
                                for (int o = 0; o < expectedOutputs; o++) {
                                    scores[o][i][j] = BigDecimal.valueOf(-1);
                                }
                                e.printStackTrace();
                            } finally {
                                done++;
                                System.out.println(
                                        System.currentTimeMillis() / 1000 +
                                                " " + done);
                            }
                        }
                    }
                }
            }
        }
        executor.shutdown();

        // Write to the output file(s)!!!
        for (int fnum = 0; fnum < expectedOutputs; fnum++) {
            BufferedWriter writer = Files.newBufferedWriter(outputs[fnum]);
            // write the header line
            String corner = "#header" + fnum + " ";
            writer.write(corner + StringUtils.repeat(" ",
                                                     maxFilenameLength -
                                                             corner.length()));
            for (Object file : files) {
                File f = (File) file;
                writer.write(f.getName() + StringUtils.repeat(" ",
                                                              maxFilenameLength -
                                                                      f.getName().length()));
            }
            writer.write("\n");

            // write the contents
            for (int i = 0; i < files.length; i++) {
                File P = (File) files[i];
                writer.write(P.getName() + StringUtils.repeat(" ",
                                                              maxFilenameLength -
                                                                      P.getName().length()));
                for (int j = 0; j < files.length; j++) {
                    String score = scores[fnum][j][j].toString();
                    writer.write(score + StringUtils.repeat(" ",
                                                            maxFilenameLength -
                                                                    score.length()));
                }
                writer.write("\n");
            }
            writer.close();
        }

        System.out.println(
                "KLMatrixGenerator finished running in "
                        + (System.nanoTime() - startTime) / 1000000 + "ms");
    }
}
