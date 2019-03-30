package cmf;

import java.util.ArrayList;

/**
 *
 * return output by execute shell command runCmd
 */
public class cmdOut {

    private ArrayList<String> out;
    private boolean extVal;
    private int exitCode;

    cmdOut(ArrayList<String> out, boolean extVal) {
        this.out = out;
        this.extVal = extVal;
    }

    cmdOut(ArrayList<String> out, boolean extVal, int exitCode) {
        this.out = out;
        this.extVal = extVal;
        this.exitCode = exitCode;
    }

    public boolean getextVal() {
        return extVal;
    }

    public boolean isExtVal() {
        return extVal;
    }

    public void setExtVal(boolean extVal) {
        this.extVal = extVal;
    }

    public int getExitCode() {
        return exitCode;
    }

    public void setExitCode(int exitCode) {
        this.exitCode = exitCode;
    }

    public ArrayList<String> getOut() {
        return out;
    }

}
