package util;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

import java.awt.*;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.file.Path;
import java.util.*;
import java.util.List;

public class SimpleBarChart implements StatsGraph {
    private Map<String, BigInteger> values = new HashMap<>();
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
        List<Map.Entry<String, BigInteger>> list = new LinkedList<>(
                values.entrySet());
        list.sort(Comparator.comparing(Map.Entry::getValue));
        for (Map.Entry<String, BigInteger> entry : list) {
            String item = entry.getKey();
            dataset.addValue(
                    values.get(item),
                    "",
                    item
            );
        }

        JFreeChart chart = ChartFactory.createBarChart(
                null,
                "",
                vertLabel,
                dataset,
                PlotOrientation.HORIZONTAL,
                false, false, false
        );
        chart.getCategoryPlot().getDomainAxis().setTickLabelFont(
                new Font(null, Font.PLAIN, 10));
        ChartUtils.saveChartAsPNG(this.output.toFile(),
                                  chart,
                                  800,
                                  600);
    }
}
