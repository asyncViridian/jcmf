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
        // ignore cases where gap is greater than sequence length
        // TODO also there are other cases where gap is at the very end
        // TODO and therefore erroneously included. Figure this bug out.
        if (gapLength.compareTo(seqLength) > 0) {
            // TODO I might want to take a different approach here...
            return;
        }
        // Add to the dataset.
        this.lengthToGaps.add(seqLength,
                              (new BigDecimal(gapLength))
                                      .multiply(BigDecimal.valueOf(100))
                                      .divide(new BigDecimal(seqLength),
                                              RoundingMode.HALF_EVEN));
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
