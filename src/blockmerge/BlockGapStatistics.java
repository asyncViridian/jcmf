package blockmerge;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.nio.file.Path;

public class BlockGapStatistics {
    private XYSeries lengthToGaps = new XYSeries("Gaps");
    private Path output;

    public BlockGapStatistics(Path output) {
        this.output = output;
    }

    public void addGap(BigInteger seqLength, BigInteger gapLength) {
        this.lengthToGaps.add(seqLength,
                              (new BigDecimal(gapLength))
                                      .divide(new BigDecimal(seqLength),
                                              RoundingMode.CEILING));
    }

    public void write() throws IOException {
        JFreeChart chart = ChartFactory.createScatterPlot(
                null,
                "Merged sequence length",
                "Gaps percentage",
                new XYSeriesCollection(this.lengthToGaps)
        );
        ChartUtils.saveChartAsPNG(this.output.toFile(), chart, 800, 600);
    }
}
