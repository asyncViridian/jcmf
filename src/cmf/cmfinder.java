package cmf;

/**
 *
 * cmfinder04.pl
 */
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.*;

public class cmfinder {

    // default parameters
    int CAND = 40;
    int MAXSPAN1 = 100;
    int MINSPAN1 = 30;
    int MAXSPAN2 = 100;
    int MINSPAN2 = 40;
    int CLUSTER = 3;
    double FRACTION = 0.8;
    int SINGLE = 5;
    int DOUBLE = 5;
    int verbose = 0;
    int help = 0;
    int COMBINE = 0;
    String cand_weight_option = "";
    String cmfinderBaseExe = "cmfinder04";
    // ($likeold,$skipClustalw,$useOldCmfinder,$simpleMotifsAlreadyDone,$justGetCmfinderCommand,
    // $copyCmfinderRunsFromLog,$amaa,$version,$filterNonFrag,$fragmentary,$commaSepEmFlags,
    // $commaSepSummarizeFlags,$commaSepCandfFlags,$saveTimer,$allCpus,$cpu,$candsParallel,$outFileSuffix,
    // $columnOnlyBasePairProbs);
    int emulate_apparent_bug_in_resolve_overlap = 1;

    //read json file
    public static JSONObject read_json_file(String file_name) throws JSONException, IOException {
        Path p = Paths.get(file_name);
        byte[] jsonByte = Files.readAllBytes(p);
        JSONObject jsonObject = new JSONObject(new String(jsonByte));
        return jsonObject;
    }

    public static void main(String[] args) {
        try {
            JSONObject jo = read_json_file("src/cmf/cmfinder_param.json");
            System.out.println(jo.toString());
            System.out.println(jo.get("h"));
            System.out.println(jo.getJSONObject("commaSepEmFlags").get("prior <f>"));
        } catch (JSONException ex) {
            Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(cmfinder.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
