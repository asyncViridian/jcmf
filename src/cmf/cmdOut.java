package cmf;

import java.util.ArrayList;

/**
 *
 * return output by execute shell command runCmd
 */
public class cmdOut {

    private ArrayList<String> out;
    private boolean extVal;

    cmdOut(ArrayList<String> out, boolean extVal) {
        this.out = out;
        this.extVal = extVal;
    }
    
    public boolean getextVal(){
        return extVal;
    }
    
    public ArrayList<String> getOut(){
        return out;
    }

}
