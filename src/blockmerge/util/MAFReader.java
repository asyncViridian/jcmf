package blockmerge.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class MAFReader implements Iterator<AlignmentBlock> {
    private String filepath;
    private BufferedReader r;
    private Iterator<String> it;

    public MAFReader(String filepath) throws IOException {
        this.filepath = filepath;
        this.r = Files.newBufferedReader(Paths.get("data", filepath));
        this.it = this.r.lines().iterator();
    }


    public void reset() {
        this.it = r.lines().iterator();
    }

    @Override
    public boolean hasNext() {
        return this.it.hasNext();
    }

    @Override
    public AlignmentBlock next() {
        List<String> lines = new ArrayList<>();
        while (it.hasNext()) {
            String s = it.next();
            if (!s.isEmpty()) {
                lines.add(s);
            } else {
                break;
            }
        }
        AlignmentBlock result = new AlignmentBlock(lines);
        return result;
    }
}
