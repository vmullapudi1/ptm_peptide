package utsw.joachimiak.vish;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class psBase {
	public static void main(String[] args) throws Exception {
		HashMap<Integer, ArrayList<Integer>> links = new HashMap<>();
		char[] tau2N4R =
				("MAEPR" + "QEFEV" + "MEDHA" + "GTYGL" + "GDRKD" + "QGGYT" + "MHQDQ" + "EGDTD" + "AGLKE" + "SPLQT" //0-49
						+ "PTEDG" + "SEEPG" + "SETSD" + "AKSTP" + "TAEDV" + "TAPLV" + "DEGAP" + "GKQAA" + "AQPHT" + "EIPEG" //50-99
						+ "TTAEE" + "AGIGD" + "TPSLE" + "DEAAG" + "HVTQA" + "RMVSK" + "SKDGT" + "GSDDK" + "KAKGA" + "DGKTK" //100-149
						+ "IATPR" + "GAAPP" + "GQKGQ" + "ANATR" + "IPAKT" + "PPAPK" + "TPPSS" + "GEPPK" + "SGDRS" + "GYSSP" //150-199
						+ "GSPGT" + "PGSRS" + "RTPSL" + "PTPPT" + "REPKK" + "VAVVR" + "TPPKS" + "PSSAK" + "SRLQT" + "APVPM" //200-249
						+ "PDLKN" + "VKSKI" + "GSTEN" + "LKHQP" + "GGGKV" + "QIINK" + "KLDLS" + "NVQSK" + "CGSKD" + "NIKHV" //250-299
						+ "PGGGS" + "VQIVY" + "KPVDL" + "SKVTS" + "KCGSL" + "GNIHH" + "KPGGG" + "QVEVK" + "SEKLD" + "FKDRV" //300-349
						+ "QSKIG" + "SLDNI" + "THVPG" + "GGNKK" + "IETHK" + "LTFRE" + "NAKAK" + "TDHGA" + "EIVYK" + "SPVVS" //350-399
						+ "GDTSP" + "RHLSN" + "VSSTG" + "SIDMV" + "DSPQL" + "ATLAD" + "EVSAS" + "LAKQG" + "L").toCharArray();

		for (int i = 0; i < tau2N4R.length; i++) {
			char j = tau2N4R[i];
			ArrayList<Integer> basicSites = new ArrayList<>();
			if (j == 'S' || j == 'T' || j == 'Y') {
				int beginning = i - 25 >= 0 ? (i - 25) : 0;
				int end = i + 25 < tau2N4R.length ? (i + 25) : tau2N4R.length - 1;
				for (int k = beginning; k <= end; k++) {
					char z = tau2N4R[k];
					if (z == 'R' || z == 'K' || z == 'H') {
						basicSites.add(k);
					}

				}
				links.put(i, basicSites);
			}
		}

		FileWriter f = new FileWriter("Tau Links25.csv", false);
		f.write("Protein 1,Protein 2,AbsPos1,AbsPos2\n");
		StringBuilder output = new StringBuilder(26736);
		for (int i : links.keySet()) {
			ArrayList<Integer> bases = links.get(i);
			for (int asdf : bases) {
				output.append("Tau2N4R,Tau2N4R," + (1 + i) + "," + (1 + asdf) + "\n");
			}
		}
		f.write(output.toString());
		f.close();
	}
}
