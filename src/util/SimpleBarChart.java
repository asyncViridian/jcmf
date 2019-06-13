package util;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.DefaultCategoryDataset;

import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Path;
import java.util.HashMap;

public class SimpleBarChart implements StatsGraph {
    private HashMap<String, BigInteger> values = new HashMap<>();
    private String key;
    private String vertLabel;
    private Path output;

    public SimpleBarChart(Path output,
                          String key,
                          String vertLabel) {
        this.output = output;
        this.key = key;
        this.vertLabel = vertLabel;
    }

    public void addValue(String item) {
        if (!values.containsKey(item)) {
            values.put(item, BigInteger.ZERO);
        }
        values.put(item, values.get(item).add(BigInteger.ONE));
    }

    @Override
    public void write() throws IOException {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (String item : values.keySet()) {
            dataset.addValue(
                    values.get(item),
                    "rowkey_" + item,
                    item
            );
        }

        JFreeChart chart = ChartFactory.createBarChart(
                null,
                "categoryAxisLabel",
                "valueAxisLabel",
                dataset
        );
        ChartUtils.saveChartAsPNG(this.output.toFile(), chart, 800, 600);
    }
}
