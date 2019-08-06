package crosscompare;

import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;
import util.StringManip;

import java.io.*;
import java.math.BigDecimal;
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
    private static final String SCORE_CALCULATE_CMD =
            "./KLDivergenceEstimate.sh";

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

        // Check that the script being used exists in the same directory
        if (Files.notExists(Paths.get(SCORE_CALCULATE_CMD))) {
            System.err.println("Please ensure that the command " +
                                       "'" + SCORE_CALCULATE_CMD + "' " +
                                       "can be executed from the working " +
                                       "directory.");
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
        maxFilenameLength += 2;

        // Start constructing our matrix
        Future[][] scoreFutures = new Future[files.length][files.length];
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
                    String[] cmd = {SCORE_CALCULATE_CMD,
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
                        }
                        reader.close();
                        return BigDecimal.valueOf(
                                Double.valueOf(pOutput.toString()));
                    } catch (IOException e) {
                        e.printStackTrace();
                        // something went wrong!!!
                        // note that K-L divergence is never negative so this
                        // is a clear sign of Bad Things Happened(tm)
                        return BigDecimal.valueOf(-1);
                    }
//                    finally {
//                         TODO delete the tempDir???
//                         Files.deleteIfExists(tempDir);
//                    }
                });
            }
        }
        // Collect the results of our matrix
        BigDecimal[][] scores = new BigDecimal[files.length][files.length];
        for (BigDecimal[] a : scores) {
            for (int i = 0; i < a.length; i++) {
                a[i] = null;
            }
        }
        int done = 0;
        while (done != (files.length * files.length)) {
            TimeUnit.SECONDS.sleep(1);
            for (int i = 0; i < files.length; i++) {
                for (int j = 0; j < files.length; j++) {
                    if (scores[i][j] == null) {
                        if (scoreFutures[i][j].isDone()) {
                            try {
                                scores[i][j] =
                                        (BigDecimal) scoreFutures[i][j].get();
                            } catch (InterruptedException | ExecutionException e) {
                                // this shouldn't happen if its done
                                // but we will mark it with a neat -1 anyway
                                scores[i][j] = BigDecimal.valueOf(-1);
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

        // Write to the output file!!!
        // write the header line
        String corner = "#header ";
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
                String score = scores[i][j].toString();
                writer.write(score + StringUtils.repeat(" ",
                                                        maxFilenameLength -
                                                                score.length()));
            }
            writer.write("\n");
        }
        writer.close();

        System.out.println(
                "KLMatrixGenerator finished running in "
                        + (System.nanoTime() - startTime) / 1000000 + "ms");
    }
}
