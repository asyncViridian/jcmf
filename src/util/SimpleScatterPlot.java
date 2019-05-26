package util;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.io.IOException;
import java.math.BigDecimal;
import java.nio.file.Path;

public class SimpleScatterPlot implements StatsGraph {
    private XYSeries dataset = new XYSeries("Gaps");
    private Path output;
    private String title;
    private String xAxisLabel;
    private String yAxisLabel;

    public SimpleScatterPlot(Path output, String title, String xAxisLabel,
                             String yAxisLabel) {
        this.output = output;
        this.title = title;
        this.xAxisLabel = xAxisLabel;
        this.yAxisLabel = yAxisLabel;
    }

    public void addValue(BigDecimal seqLength, BigDecimal gapPercentage) {
        // Add to the dataset.
        this.dataset.add(seqLength, gapPercentage);
    }

    @Override
    public void write() throws IOException {
        JFreeChart chart = ChartFactory.createScatterPlot(
                this.title,
                this.xAxisLabel,
                this.yAxisLabel,
                new XYSeriesCollection(this.dataset)
        );
        ChartUtils.saveChartAsPNG(this.output.toFile(), chart, 800, 600);
    }
}
