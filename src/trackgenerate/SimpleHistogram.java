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
    private ArrayList<Double> scores = new ArrayList<>();
    private String key;
    private Path output;

    public SimpleHistogram(Path output, String key) {
        this.output = output;
        this.key = key;
    }

    public void addScore(BigDecimal score) {
        scores.add(score.doubleValue());
    }

    public void write() throws IOException {
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        double[] temp = new double[scores.size()];
        Iterator<Double> it = scores.iterator();
        int i = 0;
        while (it.hasNext()) {
            temp[i] = it.next();
            i++;
        }
        dataset.addSeries(key, temp, 20);

        JFreeChart chart = ChartFactory.createHistogram(
                null,
                "Number of motifs",
                "Score",
                dataset);
        ChartUtils.saveChartAsPNG(this.output.toFile(), chart, 800, 600);
    }
}
