package trackgenerate;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

import java.io.IOException;
import java.math.BigDecimal;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;

public class SimpleHistogram {
    private ArrayList<Double> values = new ArrayList<>();
    private String key;
    private String xLabel;
    private String yLabel;
    private Path output;
private int numBins;

    public SimpleHistogram(Path output, String key, String xLabel,
                           String yLabel, int numBins) {
        this.output = output;
        this.key = key;
        this.xLabel = xLabel;
        this.yLabel = yLabel;
        this.numBins=numBins;
    }

    public void addValue(BigDecimal score) {
        values.add(score.doubleValue());
    }

    public void write() throws IOException {
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        double[] temp = new double[values.size()];
        Iterator<Double> it = values.iterator();
        int i = 0;
        while (it.hasNext()) {
            temp[i] = it.next();
            i++;
        }
        dataset.addSeries(key, temp, numBins);

        JFreeChart chart = ChartFactory.createHistogram(
                null,
                xLabel,
                yLabel,
                dataset);
        ChartUtils.saveChartAsPNG(this.output.toFile(), chart, 800, 600);
    }
}
