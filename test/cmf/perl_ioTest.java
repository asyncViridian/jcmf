/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cmf;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * @author james
 */
public class perl_ioTest {

    public perl_ioTest() {
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
        System.out.println("remove_gap");
        String s = "";
        String expResult = "";
        String result = perl_io.remove_gap(s);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of make_string method, of class perl_io.
     */
    @Test
    public void testMake_string() {
        System.out.println("make_string");
        int n = 0;
        String ch = "";
        String expResult = "";
        String result = perl_io.make_string(n, ch);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of pad_string method, of class perl_io.
     */
    @Test
    public void testPad_string() {
        System.out.println("pad_string");
        String seq = "";
        int n = 0;
        String ch = "";
        int dir = 0;
        String expResult = "";
        String result = perl_io.pad_string(seq, n, ch, dir);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of main method, of class perl_io.
     */
    @Test
    public void testMain() {
        System.out.println("main");
        String[] args = null;
        perl_io.main(args);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

}
