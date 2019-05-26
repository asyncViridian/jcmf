package util;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.title.Title;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;

public class SimpleHistogram implements StatsGraph {
    private ArrayList<BigDecimal> values = new ArrayList<>();
    private String key;
    private String xLabel;
    private String yLabel;
    private Path output;
    private int numBins;

    public SimpleHistogram(Path output,
                           String key,
                           String xLabel,
                           String yLabel,
                           int numBins) {
        this.output = output;
        this.key = key;
        this.xLabel = xLabel;
        this.yLabel = yLabel;
        this.numBins = numBins;
    }

    public void addValue(BigDecimal score) {
        values.add(score);
    }

    @Override
    public void write() throws IOException {
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        double[] temp = new double[values.size()];
        Iterator<BigDecimal> it = values.iterator();
        int i = 0;
        while (it.hasNext()) {
            temp[i] = it.next().doubleValue();
            i++;
        }
        dataset.addSeries(key, temp, numBins);

        JFreeChart chart = ChartFactory.createHistogram(
                null,
                xLabel,
                yLabel,
                dataset);
        chart.addSubtitle(
                new TextTitle("average = " + this.getAverage().toString()));
        ChartUtils.saveChartAsPNG(this.output.toFile(), chart, 800, 600);
    }

    public BigDecimal getAverage() {
        return this.values.stream().reduce(BigDecimal.ZERO,
                                           BigDecimal::add)
                .divide(BigDecimal.valueOf(this.values.size()),
                        RoundingMode.HALF_EVEN);
    }
}
