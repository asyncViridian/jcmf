package util;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Compares species ordering based on a known list of preferred ordering.
 * <p>
 * Any elements that are not in the list are ordered after all other
 * elements and are considered equal to each other. Null values are
 * considered to not be on the list.
 */
public class PhyloComparator implements Comparator<String> {
    /**
     * The preferred ordering of elements.
     * Does not include species outside of those used in MULTIZ 100.
     */
    private static List<String> pref = new ArrayList<>();

    static {
        pref.add("hg38");
        pref.add("panTro4");
        pref.add("gorGor3");
        pref.add("ponAbe2");
        pref.add("nomLeu3");
        pref.add("macFas5");
        pref.add("rheMac3");
        pref.add("papAnu2");
        pref.add("chlSab2");
        pref.add("saiBol1");
        pref.add("calJac3");
        pref.add("otoGar3");
        pref.add("tupChi1");
        pref.add("oryCun2");
        pref.add("ochPri3");
        pref.add("hetGla2");
        pref.add("cavPor3");
        pref.add("chiLan1");
        pref.add("octDeg1");
        pref.add("speTri2");
        pref.add("jacJac1");
        pref.add("mm10");
        pref.add("rn6");
        pref.add("micOch1");
        pref.add("criGri1");
        pref.add("mesAur1");
        pref.add("eriEur2");
        pref.add("sorAra2");
        pref.add("conCri1");
        pref.add("susScr3");
        pref.add("vicPac2");
        pref.add("camFer1");
        pref.add("turTru2");
        pref.add("orcOrc1");
        pref.add("panHod1");
        pref.add("bosTau8");
        pref.add("oviAri3");
        pref.add("capHir1");
        pref.add("pteAle1");
        pref.add("pteVam1");
        pref.add("eptFus1");
        pref.add("myoDav1");
        pref.add("myoLuc2");
        pref.add("equCab2");
        pref.add("cerSim1");
        pref.add("felCat8");
        pref.add("canFam3");
        pref.add("musFur1");
        pref.add("ailMel1");
        pref.add("odoRosDiv1");
        pref.add("lepWed1");
        pref.add("dasNov3");
        pref.add("oryAfe1");
        pref.add("chrAsi1");
        pref.add("echTel2");
        pref.add("triMan1");
        pref.add("loxAfr3");
        pref.add("monDom5");
        pref.add("sarHar1");
        pref.add("macEug2");
        pref.add("ornAna1");
        pref.add("anoCar2");
        pref.add("cheMyd1");
        pref.add("chrPic2");
        pref.add("pelSin1");
        pref.add("apaSpi1");
        pref.add("allMis1");
        pref.add("anaPla1");
        pref.add("galGal4");
        pref.add("melGal1");
        pref.add("colLiv1");
        pref.add("falChe1");
        pref.add("falPer1");
        pref.add("melUnd1");
        pref.add("amaVit1");
        pref.add("araMac1");
        pref.add("pseHum1");
        pref.add("ficAlb2");
        pref.add("taeGut1");
        pref.add("zonAlb1");
        pref.add("geoFor1");
        pref.add("xenTro7");
        pref.add("latCha1");
        pref.add("lepOcu1");
        pref.add("danRer10");
        pref.add("astMex1");
        pref.add("gadMor1");
        pref.add("gasAcu1");
        pref.add("oryLat2");
        pref.add("xipMac1");
        pref.add("tetNig2");
        pref.add("fr3");
        pref.add("takFla1");
        pref.add("oreNil2");
        pref.add("neoBri1");
        pref.add("hapBur1");
        pref.add("mayZeb1");
        pref.add("punNye1");
        pref.add("petMar2");
    }

    @Override
    public int compare(String a, String b) {
        if (pref.contains(a) && pref.contains(b)) {
            return Integer.compare(pref.indexOf(a), pref.indexOf(b));
        } else if (pref.contains(a) && !pref.contains(b)) {
            return -1;
        } else if (!pref.contains(a) && pref.contains(b)) {
            return 1;
        } else {
            return 0;
        }
    }
}
