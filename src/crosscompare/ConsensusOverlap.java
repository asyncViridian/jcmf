package crosscompare;

import org.apache.commons.cli.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import util.ScoredStockholmAlignmentBlock;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class ConsensusOverlap {
    private static Options options;
    private static String srcDir;
    private static String outputFile;
    private static String chr;
    private static Boolean link;
    private static String relativeLinkPrefix = "";

    public static void main(String[] args) throws IOException {
        long startTime = System.nanoTime();
        // handle command-line argument processing :)
        // add all the arguments we need
        ConsensusOverlap.options = new Options();
        {
            ConsensusOverlap.options.addOption(
                    Option.builder("s")
                            .longOpt("srcDir")
                            .hasArg()
                            .desc("Input score files directory. Required.")
                            .required()
                            .build());
            ConsensusOverlap.options.addOption(
                    Option.builder("o")
                            .longOpt("outputFile")
                            .hasArg()
                            .desc("Output filename. Required.")
                            .required()
                            .build());
            ConsensusOverlap.options.addOption(
                    Option.builder("c")
                            .longOpt("chr")
                            .hasArg()
                            .desc("The chromosome name (of hg38) to use. " +
                                          "Required.")
                            .required()
                            .build());
            ConsensusOverlap.options.addOption(
                    Option.builder("l")
                            .longOpt("link")
                            .desc("Iff there should be a link to the details " +
                                          "page for a given motif entry.")
                            .build());
            ConsensusOverlap.options.addOption(
                    Option.builder("pre")
                            .longOpt("linkPrefix")
                            .desc("The prefix of the relative link to use. " +
                                          "Defaults to blank.")
                            .hasArg()
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            ConsensusOverlap.srcDir = line.getOptionValue("s");
            // set the output filename
            ConsensusOverlap.outputFile = line.getOptionValue("o");
            // set the chr
            ConsensusOverlap.chr = line.getOptionValue("c");
            // set whether to link
            ConsensusOverlap.link = line.hasOption("l");
            // set the link prefix
            if (line.hasOption("pre")) {
                ConsensusOverlap.relativeLinkPrefix
                        = line.getOptionValue("pre");
            }
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "ConsensusOverlap -s <dirname> -o " +
                            "<filename> -c <chromosome> -p <position>",
                    options);
            return;
        }
        // Get the input & output files set up
        Path source = Paths.get(ConsensusOverlap.srcDir);
        Path output = Paths.get(ConsensusOverlap.outputFile);
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);

        // Sort the files by motif starting position
        File[] allFiles = source.toFile().listFiles();
        List<File> sortedFiles = new ArrayList<>();
        // just pick out the files that are not html files
        for (int i = 0; i < allFiles.length; i++) {
            if (!allFiles[i].getName().endsWith(".html")
                    && allFiles[i].getName().contains(chr + ".")) {
                sortedFiles.add(allFiles[i]);
            }
        }
        // sort all the motif files
        // sorted by beginning coordinate of consensus motif structure
        sortedFiles.sort(
                (File f1, File f2) -> {
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
                    int fileCompare = f1.getName().compareTo(f2.getName());

                    if (startPos != 0) {
                        return startPos;
                    } else if (seqCompare != 0) {
                        return seqCompare;
                    } else {
                        return fileCompare;
                    }
                });

        // Collect all the species involved
        Set<String> species = new HashSet<>();

        // Group the files by overlapping motif structure
        BigInteger groupStartPos = BigInteger.valueOf(0);
        BigInteger groupEndPos = BigInteger.valueOf(0);
        List<Pair<List<File>, BigInteger>> motifGroupings = new LinkedList<>();
        List<File> currMotifGroup = new LinkedList<>();
        for (int i = 0; i < sortedFiles.size(); i++) {
            // read the actual file
            ScoredStockholmAlignmentBlock block = null;
            try {
                block = ScoredStockholmAlignmentBlock
                        .constructFromScore(sortedFiles.get(i));
            } catch (FileNotFoundException e) {
                throw new RuntimeException("file is now missing?");
            }
            // save the species that this block supports
            species.addAll(block.scores.keySet());
            // get the specific human block now
            ScoredStockholmAlignmentBlock.Source src =
                    block.sources.get("hg38");

            // get position indexes
            BigInteger startPos = src.totalSpan.getLeft()
                    .add(block.intervals.get("hg38").getLeft());
            BigInteger endPos = src.totalSpan.getLeft()
                    .add(block.intervals.get("hg38").getRight());

            // place it into a group or add a group to the list
            if (startPos.compareTo(groupEndPos) > 0) {
                // if the start of this motif is after the end of the last group
                if (currMotifGroup.size() != 0) {
                    // if there was a previous group, close it
                    Pair<List<File>, BigInteger> pair = new ImmutablePair<>(
                            currMotifGroup,
                            groupStartPos
                    );
                    motifGroupings.add(pair);
                }
                // start a new group
                groupStartPos = startPos;
                currMotifGroup = new LinkedList<>();
                currMotifGroup.add(sortedFiles.get(i));
                groupEndPos = endPos;
            } else {
                // add this motif to the group
                currMotifGroup.add(sortedFiles.get(i));
                groupEndPos = groupEndPos.max(endPos);
            }
        }
        // close the fencepost (finish off the last grouping)
        Pair<List<File>, BigInteger> pair = new ImmutablePair<>(
                currMotifGroup,
                groupStartPos
        );
        motifGroupings.add(pair);

        // Start writing the actual file
        writer.write("<html><body style=\"font-family:monospace;\">");
        writer.write("<p style=\"text-align:center;\">" +
                             ConsensusOverlap.chr + " motif clusters" +
                             "</p>");
        for (Pair<List<File>, BigInteger> currentlyPrinting : motifGroupings) {
            writer.write("<hr>");
            writer.write("<p>Motif cluster starting at "
                                 + currentlyPrinting.getRight() + "</p>");

            // one table for each motif grouping
            writer.write("<table><tr>");
            // Write the left column: filenames with links
            {
                writer.write("<td><div style=\"overflow-x:scroll;width:10vw;" +
                                     "white-space:nowrap;\">");
                writer.write("<table style=\"" +
                                     "border-collapse: collapse;" +
                                     "border-spacing: 0;\">");
                for (File f : currentlyPrinting.getLeft()) {
                    writer.write("<tr><td>&nbsp;</td></tr>");
                    writer.write("<tr><td style=\"white-space: nowrap;\">");
                    ScoredStockholmAlignmentBlock block =
                            ScoredStockholmAlignmentBlock.constructFromScore(f);
                    ScoredStockholmAlignmentBlock.Source hgsrc =
                            block.sources.get(
                                    "hg38");

                    // skip nonmatching chr
                    if (!hgsrc.chr.equals(ConsensusOverlap.chr)) {
                        continue;
                    }

                    if (ConsensusOverlap.link) {
                        writer.write("<a href=\""
                                             + ConsensusOverlap.relativeLinkPrefix + f.getName() + ".html\">");
                    }
                    writer.write(f.getName());
                    if (ConsensusOverlap.link) {
                        writer.write("</a>");
                    }
                    writer.write("</td></tr>");
                }
                writer.write("</table>");
                writer.write("</div></td>");
            }
            // Write the right column: consensus structure aligned
            {
                writer.write("<td><div style=\"overflow-x:scroll;width:85vw;" +
                                     "white-space:nowrap;\">");
                writer.write("<table style=\"" +
                                     "border-collapse: collapse;" +
                                     "border-spacing: 0;\">");
                for (File f : currentlyPrinting.getLeft()) {
                    ScoredStockholmAlignmentBlock block =
                            ScoredStockholmAlignmentBlock.constructFromScore(f);
                    ScoredStockholmAlignmentBlock.Source hgsrc =
                            block.sources.get(
                                    "hg38");

                    // skip nonmatching chr
                    if (!hgsrc.chr.equals(ConsensusOverlap.chr)) {
                        continue;
                    }

                    // write the sequence itself
                    writer.write("<tr><td style=\"white-space: nowrap;\">");
                    // insert spaces corresponding to the actual position
                    // where it starts
                    for (int i = 0;
                         hgsrc.totalSpan.getLeft()
                                 .add(block.intervals.get("hg38").getLeft())
                                 .subtract(currentlyPrinting.getRight())
                                 .compareTo(
                                         BigInteger.valueOf(i)) > 0;
                         i++) {
                        writer.write("&nbsp;");
                    }
                    writer.write(block.sources.get("hg38").sequence);
                    writer.write("</td></tr>");

                    // write the folding structure
                    writer.write("<tr><td style=\"white-space: nowrap;\">");
                    // insert spaces corresponding to the actual position
                    // where it starts
                    for (int i = 0;
                         hgsrc.totalSpan.getLeft()
                                 .add(block.intervals.get("hg38").getLeft())
                                 .subtract(currentlyPrinting.getRight())
                                 .compareTo(
                                         BigInteger.valueOf(i)) > 0;
                         i++) {
                        writer.write("&nbsp;");
                    }
                    writer.write(block.SS_cons);
                    writer.write("</td></tr>");
                }
                writer.write("</table>");
                writer.write("</div></td>");
            }
            // Write the scoretable: relative quality scores for all species
            {
                // TODO order the species display a more reasonable way
                List<String> lspecies = new ArrayList<>();
                lspecies.addAll(species);
                writer.write("<td><div style=\"overflow-x:scroll;width:85vw;" +
                                     "white-space:nowrap;\">");
                writer.write("<table style=\"" +
                                     "border-collapse: collapse;" +
                                     "border-spacing: 0;\">" +
                                     "<tr>");
                for (String s : lspecies) {
                    writer.write("<th style=\"white-space: nowrap;\">");
                    writer.write(s);
                    writer.write("</th>");
                }
                writer.write("</tr>");
                boolean first = true;
                for (File f : currentlyPrinting.getLeft()) {
                    ScoredStockholmAlignmentBlock block =
                            ScoredStockholmAlignmentBlock.constructFromScore(f);

                    if (!first) {
                        writer.write("<tr>" +
                                             "<td colspan=\""
                                             + lspecies.size() + "\">" +
                                             "&nbsp;</td>" +
                                             "</tr>");
                    }
                    writer.write("<tr>");
                    for (String s : lspecies) {
                        writer.write("<td>");
                        if (block.scores.containsKey(s)) {
                            writer.write("" + block.scores.get(s));
                        } else {
                            // there is no score
                            writer.write("&nbsp;");
                        }
                        writer.write("</td>");
                    }
                    writer.write("</tr>");
                    first = false;
                }
                writer.write("</table>");
                writer.write("</div></td>");
            }
            writer.write("</tr></table>");
        }
        writer.write("</body></html>");
        writer.close();

        System.out.println("ConsensusOverlap finished running in "
                                   + (System.nanoTime() - startTime) / 1000000 + "ms");
    }
}
