package cmf;

// RuntimeException are unchecked while Exception are checked 
public class MyException extends Exception {

    public MyException(String s) {
        super(s);
    }
}
