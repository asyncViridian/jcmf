package trackgenerate;

import org.apache.commons.cli.*;

import java.nio.file.Path;
import java.nio.file.Paths;

public class PageGenerator {
    private static Options options;
    private static String srcFile;
    private static String outputDir;

    public static void main(String[] args) {
        PageGenerator.options = new Options();
        {
            PageGenerator.options.addOption(
                    Option.builder("s")
                            .longOpt("srcFile")
                            .hasArg()
                            .desc("Input score file")
                            .required()
                            .build());
            PageGenerator.options.addOption(
                    Option.builder("o")
                            .longOpt("outputDir")
                            .hasArg()
                            .desc("Output directory " +
                                          "to write the matching HTML file in")
                            .required()
                            .build());
        }
        // Parse the commandline arguments
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            // set the src dirname
            PageGenerator.srcFile = line.getOptionValue("s");
            // set the output filename
            PageGenerator.outputDir = line.getOptionValue("o");
        } catch (ParseException exp) {
            // something went wrong
            System.err.println("Parsing failed. Reason: " + exp.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(
                    "PageGenerator [options] -s <filename> -o <filename>",
                    options);
            return;
        }

        // generate the output HTML file
        Path output = Paths.get(PageGenerator.outputDir, "multi.bed");
    }
}
