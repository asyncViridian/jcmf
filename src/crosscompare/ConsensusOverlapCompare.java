package crosscompare;

import org.apache.commons.cli.Options;
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
    private static String srcDir = "notes/trackHub/hg38/bigscanf_src";
    private static String outputDir = "output";
    private static String chr = "chr12";
    private static BigInteger startPoint = BigInteger.valueOf(62602730);

    public static void main(String[] args) throws IOException {
        Path source = Paths.get(ConsensusOverlapCompare.srcDir);
        Path output = Paths.get(ConsensusOverlapCompare.outputDir,
                                "consensus_comparison.html");
        Files.deleteIfExists(output);
        Files.createFile(output);
        BufferedWriter writer = Files.newBufferedWriter(output);
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

                        int i = hgsrc1.totalSpan.getLeft().compareTo(
                                hgsrc2.totalSpan.getLeft());
                        if (i != 0) {
                            return i;
                        } else {
                            return block1.SS_cons.compareTo(block2.SS_cons);
                        }
                    });

        writer.write("<html><body style=\"font-family:monospace;\">");
        writer.write("<table><tr>");
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

            for (int i = 0;
                 hgsrc.totalSpan.getLeft().subtract(
                         startPoint).compareTo(BigInteger.valueOf(i)) > 0;
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
