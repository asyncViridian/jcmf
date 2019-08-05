package blockmerge;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import util.StringManip;

import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Paths;

import static org.junit.jupiter.api.Assertions.*;

class BlockMergerTest {

    @BeforeEach
    void setUp() {
    }

    @AfterEach
    void tearDown() {
    }

    @Test
    void main() {
        try {
            Files.list(Paths.get("test/blockmerge/data"))
                    .map(p -> p.getFileName().toString())
                    .map(s -> s.substring(0, s.indexOf(".maf")))
                    .forEach(filename -> {
                        String[] args = {"-numBlocksPerOutput", "2",
                                "--minNumSpecies", "0",
                                "-sd", "test/blockmerge/data",
                                "-s", filename + ".maf",
                                "-od", "test/blockmerge/output",
                                "-o", filename,
                                "-ot", "maf"};
                        try {
                            BlockMerger.main(args);
                        } catch (Exception e) {
                            e.printStackTrace();
                            fail("check test integrity");
                        }

                        // TODO add check against expected behavior???
                    });
        } catch (Exception e) {
            e.printStackTrace();
            fail("check test integrity");
        }
    }

    @Test
    void testRepeat() {
        String[] s = {"", "a", "aB"};

        // Test repeating zero times
        for (int i = 0; i < s.length; i++) {
            assertEquals("", StringManip.repeat(s[i], BigInteger.ZERO));
        }

        // Test repeating one time
        for (int i = 0; i < s.length; i++) {
            assertEquals(s[i], StringManip.repeat(s[i], BigInteger.ONE));
        }

        // Test repeating two times
        for (int i = 0; i < s.length; i++) {
            assertEquals(s[i] + s[i],
                         StringManip.repeat(s[i], BigInteger.valueOf(2)));
        }
    }
}