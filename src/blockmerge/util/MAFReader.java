package blockmerge.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class MAFReader {
    private String filepath;

    public MAFReader(String filepath) {
        this.filepath = filepath;
    }

    public AlignmentBlock getNextBlock() throws IOException {
        List<String> lines = new ArrayList<>();
        BufferedReader r = Files.newBufferedReader(Paths.get("data", filepath));
        Iterator<String> it = r.lines().iterator();
        while (it.hasNext()) {
            String s = it.next();
            if (!s.isEmpty()) {
                lines.add(s);
            } else {
                break;
            }
        }
        AlignmentBlock result = new AlignmentBlock(lines);
        System.out.print("done");
        return result;
    }

}
