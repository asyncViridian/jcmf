package trackgenerate;

import org.apache.commons.cli.*;
import util.ScoredStockholmAlignmentBlock;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.MessageFormat;

public class PageGenerator {
    private static Options options;
    private static String srcFile;

    /**
     * Template string for the entire HTML.
     * <p>
     * {0} = filename that is being html-ified
     * <p>
     * {1} = Pair-posterior score
     * <p>
     * {2} = RNA-posterior score
     * <p>
     * {3} = Sequence data table HTML code (not including table-tag)
     * <p>
     * {4} = Alignment data table HTML code (not including table-tag)
     */
    private static final String HTML_TEMPLATE =
            "<!DOCTYPE html><html><head>" +
                    "<title>{0}</title>" +
                    "</head><body style=\"font-family:monospace;\"><div " +
                    "style=\"text-align:center;\">{0}<br/>" +
                    "<a href=\"{0}\">link to source score file</a>" +
                    "</div><hr/><div>" +
                    "Total pair posterior {1}" +
                    "<br/>" +
                    "Total RNA posterior {2}" +
                    "</div><hr/><div>" +
                    "Sequence segments:" +
                    "<br/>" +
                    "<table>{3}</table>" +
                    "</div><hr/><div>Alignment: (use shift-scroll to scroll " +
                    "sideways)<br/>" +
                    "<table><tr>{4}</tr></table>" +
                    "</div></body></html>";

    /**
     * Template string for a single table-row of sequence data table.
     * <p>
     * {0} = long species information string
     * <p>
     * {1} = coordinate span string
     * <p>
     * {2} = quality score
     */
    private static final String SEQ_DATA_TABLE_ROW_TEMPLATE =
            "<tr><td>{0}</td>" +
                    "<td>{1}</td>" +
                    "<td>{2}</td></tr>";

    /**
     * Template string for the alignment data table columns...
     * <p>
     * {0} = string for species name column. Consists of each species name
     * separated by two break tags.
     * <p>
     * {1} = string for alignment data column. Consists of alternating lines
     * of sequence data and motif-alignment data, separated by break tags.
     */
    private static final String ALIGN_DATA_TABLE_TEMPLATE =
            "<td>{0}</td>" +
                    "<td><div style=\"overflow-x:scroll;width:85vw;" +
                    "white-space:nowrap;\">{1}</div></td>";

    public static void main(String[] args) throws IOException {
        PageGenerator.options = new Options();
        {
            PageGenerator.options.addOption(
                    Option.builder("s")
                            .longOpt("srcFile")
                            .hasArg()
                            .desc("Input score file")
                            .required()
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            PageGenerator.srcFile = line.getOptionValue("s");
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "PageGenerator [options] -s <filename>",
                    options);
            return;
        }

        // generate the output HTML file
        Path output = Paths.get(PageGenerator.srcFile + ".html");
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);

        // read the input file
        File file = Paths.get(PageGenerator.srcFile).toFile();
        ScoredStockholmAlignmentBlock block =
                ScoredStockholmAlignmentBlock.constructFromScore(file);
        if (block == null) {
            // something was malformatted
            return;
        }

        // Generate the output HTML...

        // First generate both the component tables
        StringBuilder srcTableContents = new StringBuilder();
        StringBuilder alignTableSpecies = new StringBuilder();
        StringBuilder alignTableData = new StringBuilder();
        for (String species : block.getSpecies()) {
            ScoredStockholmAlignmentBlock.Source source =
                    block.sources.get(species);
            // Build the sequence summary data table
            srcTableContents.append(MessageFormat.format(
                    PageGenerator.SEQ_DATA_TABLE_ROW_TEMPLATE,
                    source.inputLine,
                    source.totalSpan.getLeft() + ".." + source.totalSpan.getRight(),
                    block.scores.get(species)
            ));
            // Build the alignment data table species column
            alignTableSpecies.append(species);
            alignTableSpecies.append("<br/><br/>");
            // Build the alignment data table data column
            String sequence = source.sequence;
            String motif = source.motif;
            StringBuilder sequenceLine = new StringBuilder();
            StringBuilder motifLine = new StringBuilder();
            boolean light = false;
            while (sequence.length() != 0) {
                if (sequence.length() != motif.length()) {
                    throw new IllegalStateException("mismatched seq-motif len");
                }

                sequenceLine.append("<span style=\"color:");
                motifLine.append("<span style=\"color:");
                if (light) {
                    sequenceLine.append("grey");
                    motifLine.append("grey");
                    // add light shaded
                } else {
                    sequenceLine.append("black");
                    motifLine.append("black");
                    // add dark shaded
                }
                sequenceLine.append("\">");
                motifLine.append("\">");
                sequenceLine.append(
                        sequence,
                        0,
                        Math.min(50, sequence.length()));
                motifLine.append(
                        motif,
                        0,
                        Math.min(50, motif.length()));
                sequenceLine.append("</span>");
                motifLine.append("</span>");

                // update cycle
                light = !light;
                sequence = sequence.substring(Math.min(50, sequence.length()));
                motif = motif.substring(Math.min(50, motif.length()));
            }
            alignTableData.append(sequenceLine.toString());
            alignTableData.append("<br/>");
            alignTableData.append(motifLine.toString());
            alignTableData.append("<br/>");
        }
        // add the consensus to the table
        alignTableSpecies.append(
                "<span style=\"font-weight:bold\">" +
                        "CONSENSUS</span><br/>" +
                        "<span style=\"font-weight:bold\">RF</span><br/>");
        // data column
        String SS_cons = block.SS_cons;
        String RF = block.RF;
        StringBuilder SSLine = new StringBuilder();
        StringBuilder RFLine = new StringBuilder();
        boolean light = false;
        while (SS_cons.length() != 0) {
            if (SS_cons.length() != RF.length()) {
                throw new IllegalStateException("mismatched SS_CON, RF len");
            }

            SSLine.append("<span style=\"font-weight:bold;color:");
            RFLine.append("<span style=\"font-weight:bold;color:");
            if (light) {
                SSLine.append("grey");
                RFLine.append("grey");
                // add light shaded
            } else {
                SSLine.append("black");
                RFLine.append("black");
                // add dark shaded
            }
            SSLine.append("\">");
            RFLine.append("\">");
            SSLine.append(
                    SS_cons,
                    0,
                    Math.min(50, SS_cons.length()));
            RFLine.append(
                    RF,
                    0,
                    Math.min(50, RF.length()));
            SSLine.append("</span>");
            RFLine.append("</span>");

            // update cycle
            light = !light;
            SS_cons = SS_cons.substring(Math.min(50, SS_cons.length()));
            RF = RF.substring(Math.min(50, RF.length()));
        }
        alignTableData.append(SSLine.toString());
        alignTableData.append("<br/>");
        alignTableData.append(RFLine.toString());
        alignTableData.append("<br/>");

        // build the overall alignment data table
        String alignTableContents = MessageFormat.format(
                PageGenerator.ALIGN_DATA_TABLE_TEMPLATE,
                alignTableSpecies.toString(),
                alignTableData.toString()
        );

        // Generate the overall HTML
        writer.write(
                MessageFormat.format(
                        PageGenerator.HTML_TEMPLATE,
                        file.getName(),
                        block.pairScore,
                        block.rnaScore,
                        srcTableContents.toString(),
                        alignTableContents
                ));
        writer.close();
        System.out.println("Wrote summary html for " + PageGenerator.srcFile);
    }
}
