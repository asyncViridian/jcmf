package crosscompare;

import org.apache.commons.cli.*;
import util.ScoredStockholmAlignmentBlock;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

public class ConsensusOverlapCompare {
    private static Options options;
    private static String srcDir;
    private static String outputFile;
    private static String chr;
    private static BigInteger startPoint;

    public static void main(String[] args) throws IOException {
        // handle command-line argument processing :)
        // add all the arguments we need
        ConsensusOverlapCompare.options = new Options();
        {
            ConsensusOverlapCompare.options.addOption(
                    Option.builder("s")
                            .longOpt("srcDir")
                            .hasArg()
                            .desc("Input score files directory. Required.")
                            .required()
                            .build());
            ConsensusOverlapCompare.options.addOption(
                    Option.builder("o")
                            .longOpt("outputFile")
                            .hasArg()
                            .desc("Output filename. Required.")
                            .required()
                            .build());
            ConsensusOverlapCompare.options.addOption(
                    Option.builder("c")
                            .longOpt("chr")
                            .hasArg()
                            .desc("The chromosome name (of hg38) to use. " +
                                          "Required.")
                            .required()
                            .build());
            ConsensusOverlapCompare.options.addOption(
                    Option.builder("p")
                            .longOpt("position")
                            .hasArg()
                            .desc("The start position in the chr to use. " +
                                          "Required.")
                            .required()
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            ConsensusOverlapCompare.srcDir = line.getOptionValue("s");
            // set the output filename
            ConsensusOverlapCompare.outputFile = line.getOptionValue("o");
            // set the chr
            ConsensusOverlapCompare.chr = line.getOptionValue("c");
            // set the chr
            ConsensusOverlapCompare.startPoint = BigInteger.valueOf(
                    Long.valueOf(line.getOptionValue("p")));
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "ConsensusOverlapCompare -s <dirname> -o " +
                            "<filename> -c <chromosome> -p <position>",
                    options);
            return;
        }
        // Get the input & output files set up
        Path source = Paths.get(ConsensusOverlapCompare.srcDir);
        Path output = Paths.get(ConsensusOverlapCompare.outputFile);
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);

        // Sort the files by blockgen number
        File[] sortedFiles = source.toFile().listFiles();
        Arrays.sort(sortedFiles,
                    (f1, f2) -> {
                        if (f1.getName().endsWith(
                                ".html") && f2.getName().endsWith(".html")) {
                            return 0;
                        }
                        if (f1.getName().endsWith(".html")) {
                            return 1;
                        }
                        if (f2.getName().endsWith(".html")) {
                            return -1;
                        }
                        ScoredStockholmAlignmentBlock block1 = null;
                        ScoredStockholmAlignmentBlock block2 = null;
                        try {
                            block1 =
                                    ScoredStockholmAlignmentBlock
                                            .constructFromScore(f1);
                            block2 =
                                    ScoredStockholmAlignmentBlock
                                            .constructFromScore(f2);
                        } catch (FileNotFoundException e) {
                            e.printStackTrace();
                        }
                        ScoredStockholmAlignmentBlock.Source hgsrc1 =
                                block1.sources.get("hg38");
                        ScoredStockholmAlignmentBlock.Source hgsrc2 =
                                block2.sources.get("hg38");

                        int startPos = hgsrc1.totalSpan.getLeft().
                                add(block1.intervals.get("hg38").getLeft())
                                .compareTo(hgsrc2.totalSpan.getLeft()
                                                   .add(block2.intervals.get(
                                                           "hg38").getLeft()));
                        int seqCompare = block1.SS_cons.compareTo(
                                block2.SS_cons);
                        int fileCompare=f1.getName().compareTo(f2.getName());

                        if (startPos != 0) {
                            return startPos;
                        } else if (seqCompare != 0) {
                            return seqCompare;
                        } else {
                            return fileCompare;
                        }
                    });

        // Start writing the file
        writer.write("<html><body style=\"font-family:monospace;\">");
        writer.write("<table><tr>");
        // Write the left column: filenames
        writer.write("<td><div style=\"overflow-x:auto;width:10vw;" +
                             "white-space:nowrap;\">");
        for (File f : sortedFiles) {
            // skip HTML files, we only want scored alignments
            if (f.getName().endsWith(".html")) {
                continue;
            }

            ScoredStockholmAlignmentBlock block =
                    ScoredStockholmAlignmentBlock.constructFromScore(f);
            ScoredStockholmAlignmentBlock.Source hgsrc = block.sources.get(
                    "hg38");

            // skip nonmatching chr
            if (!hgsrc.chr.equals(ConsensusOverlapCompare.chr)) {
                continue;
            }

            writer.write(f.getName());
            writer.write("<br/>");
        }
        writer.write("</div></td>");
        // Write the right column: consensus structure aligned
        writer.write("<td><div style=\"overflow-x:auto;width:85vw;" +
                             "white-space:nowrap;\">");
        for (File f : sortedFiles) {
            // skip HTML files, we only want scored alignments
            if (f.getName().endsWith(".html")) {
                continue;
            }

            ScoredStockholmAlignmentBlock block =
                    ScoredStockholmAlignmentBlock.constructFromScore(f);
            ScoredStockholmAlignmentBlock.Source hgsrc = block.sources.get(
                    "hg38");

            // skip nonmatching chr
            if (!hgsrc.chr.equals(ConsensusOverlapCompare.chr)) {
                continue;
            }

            // insert spaces corresponding to the actual position where it
            // starts
            for (int i = 0;
                 hgsrc.totalSpan.getLeft()
                         .add(block.intervals.get("hg38").getLeft())
                         .subtract(startPoint)
                         .compareTo(
                                 BigInteger.valueOf(i)) > 0;
                 i++) {
                writer.write("&nbsp;");
            }
            writer.write(block.SS_cons);
            writer.write("<br/>");
        }
        writer.write("</div></td>");
        writer.write("</tr></table>");
        writer.write("</body></html>");
        writer.close();
    }
}
