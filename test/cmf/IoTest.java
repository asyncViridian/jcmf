/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cmf;

import java.util.Map;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * @author james
 */
public class IoTest {

    public IoTest() {
    }

    @BeforeAll
    public static void setUpClass() {
    }

    @AfterAll
    public static void tearDownClass() {
    }

    /**
     * Test of remove_gap method, of class perl_io.
     */
    @Test
    public void testRemove_gap() {
        assertEquals("ACGT", Io.remove_gap("ACGT"));
        assertEquals("AGT", Io.remove_gap("A.GT"));
        assertEquals("CGT", Io.remove_gap("-CGT"));
        assertEquals("CG", Io.remove_gap("-CG."));
        assertEquals("", Io.remove_gap("----"));
        assertEquals("", Io.remove_gap("...."));
        assertEquals("", Io.remove_gap(".-.-"));
    }

    /**
     * Test of make_string method, of class perl_io.
     */
    @Test
    public void testMake_string() {
        assertEquals("", Io.make_string(0, "a"));
        assertEquals("", Io.make_string(1, ""));
        assertEquals("a", Io.make_string(1, "a"));
        assertEquals("ab", Io.make_string(1, "ab"));
        assertEquals("ab ab ab ", Io.make_string(3, "ab "));
    }

    /**
     * Test of pad_string method, of class perl_io.
     */
    @Test
    public void testPad_string() {
        assertEquals("", Io.pad_string("", 0, "", 0));
        assertEquals("A", Io.pad_string("A", 0, "", 0));
        assertEquals("A", Io.pad_string("A", 1, "", 0));
        assertEquals("A", Io.pad_string("A", 0, "a", 0));
        assertEquals("A", Io.pad_string("A", 1, "a", 0));
        assertEquals("aaaaA", Io.pad_string("A", 5, "a", 0));
        assertEquals("a", Io.pad_string("", 1, "a", 0));
        assertEquals("aaaaa", Io.pad_string("", 5, "a", 0));
        assertEquals("", Io.pad_string("", 0, "", 1));
        assertEquals("A", Io.pad_string("A", 0, "", 1));
        assertEquals("A", Io.pad_string("A", 1, "", 1));
        assertEquals("A", Io.pad_string("A", 0, "a", 1));
        assertEquals("A", Io.pad_string("A", 1, "a", 1));
        assertEquals("Aaaaa", Io.pad_string("A", 5, "a", 1));
        assertEquals("a", Io.pad_string("", 1, "a", 1));
        assertEquals("aaaaa", Io.pad_string("", 5, "a", 1));
        assertEquals("", Io.pad_string("", 0, "", 2));
        assertEquals("A", Io.pad_string("A", 0, "", 2));
        assertEquals("A", Io.pad_string("A", 1, "", 2));
        assertEquals("A", Io.pad_string("A", 0, "a", 2));
        assertEquals("A", Io.pad_string("A", 1, "a", 2));
        assertEquals("A", Io.pad_string("A", 5, "a", 2));
    }

    /**
     * Test of read_fasta method, of class perl_io.
     */
    @Test
    public void testRead_fasta() {
        Map<String, cmf.Seq> result = Io.read_fasta("test/cmf/data/example.fasta");
        assertEquals(result.size(), 55);
        // TODO implement the rest of the test
    }

}
