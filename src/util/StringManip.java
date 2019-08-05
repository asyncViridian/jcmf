package util;

import java.math.BigInteger;

public class StringManip {
    /**
     * Return a string consisting of src repeated num times
     *
     * @param src String to repeat
     * @param num number of times to repeat this string
     * @return src repeated num times
     */
    public static String repeat(String src, BigInteger num) {
        StringBuilder builder = new StringBuilder();
        for (BigInteger i = BigInteger.ZERO; i.compareTo(num) < 0; i = i.add(
                BigInteger.ONE)) {
            builder.append(src);
        }
        return builder.toString();
    }
}
