package blockmerge.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class MAFReader {
    List<String> lines;

    public MAFReader(String filepath) throws IOException {
        this.lines = Files.readAllLines(Paths.get(filepath));
    }
}
