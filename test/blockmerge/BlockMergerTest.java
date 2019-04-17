package blockmerge;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.math.BigInteger;

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

    }

    @Test
    void testRepeat() {
        String[] s = {"", "a", "aB"};

        // Test repeating zero times
        for (int i = 0; i < s.length; i++) {
            assertEquals("", BlockMerger.repeat(s[i], BigInteger.ZERO));
        }

        // Test repeating one time
        for (int i = 0; i < s.length; i++) {
            assertEquals(s[i], BlockMerger.repeat(s[i], BigInteger.ONE));
        }

        // Test repeating two times
        for (int i = 0; i < s.length; i++) {
            assertEquals(s[i] + s[i],
                         BlockMerger.repeat(s[i], BigInteger.valueOf(2)));
        }
    }
}