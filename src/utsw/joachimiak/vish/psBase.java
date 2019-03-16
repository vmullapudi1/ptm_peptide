package utsw.joachimiak.vish;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

public class psBase {
	public static void main(String[] args) throws Exception {
		HashMap<Integer, ArrayList<Integer>> links = new HashMap<>();
		final char[] tau2N4R =
				("MAEPR" + "QEFEV" + "MEDHA" + "GTYGL" + "GDRKD" + "QGGYT" + "MHQDQ" + "EGDTD" + "AGLKE" + "SPLQT" //0-49
						+ "PTEDG" + "SEEPG" + "SETSD" + "AKSTP" + "TAEDV" + "TAPLV" + "DEGAP" + "GKQAA" + "AQPHT" + "EIPEG" //50-99
						+ "TTAEE" + "AGIGD" + "TPSLE" + "DEAAG" + "HVTQA" + "RMVSK" + "SKDGT" + "GSDDK" + "KAKGA" + "DGKTK" //100-149
						+ "IATPR" + "GAAPP" + "GQKGQ" + "ANATR" + "IPAKT" + "PPAPK" + "TPPSS" + "GEPPK" + "SGDRS" + "GYSSP" //150-199
						+ "GSPGT" + "PGSRS" + "RTPSL" + "PTPPT" + "REPKK" + "VAVVR" + "TPPKS" + "PSSAK" + "SRLQT" + "APVPM" //200-249
						+ "PDLKN" + "VKSKI" + "GSTEN" + "LKHQP" + "GGGKV" + "QIINK" + "KLDLS" + "NVQSK" + "CGSKD" + "NIKHV" //250-299
						+ "PGGGS" + "VQIVY" + "KPVDL" + "SKVTS" + "KCGSL" + "GNIHH" + "KPGGG" + "QVEVK" + "SEKLD" + "FKDRV" //300-349
						+ "QSKIG" + "SLDNI" + "THVPG" + "GGNKK" + "IETHK" + "LTFRE" + "NAKAK" + "TDHGA" + "EIVYK" + "SPVVS" //350-399
						+ "GDTSP" + "RHLSN" + "VSSTG" + "SIDMV" + "DSPQL" + "ATLAD" + "EVSAS" + "LAKQG" + "L").toCharArray();
		final int beginning = 0;
		final int end = tau2N4R.length - 1;
		for (int windowEnd = 5; windowEnd < tau2N4R.length; windowEnd++) {
			int windowStart = 0;

			for (int i = beginning; i <= end; i++) {
				char j = tau2N4R[i];
				ArrayList<Integer> basicSites = new ArrayList<>();
				if (j == 'S' || j == 'T' || j == 'Y') {
					for (int k = i - windowEnd >= 0 ? i - windowEnd : 0; ((k <= end) && (k <= i + windowEnd)); k++) {
						char z = tau2N4R[k];
						if (z == 'R' || z == 'K' || z == 'H') {
							basicSites.add(k);
						}

					}
					links.put(i, basicSites);
				}
			}
			Files.createDirectories(Paths.get("contactOutput"));
			outputToFile(windowStart, windowEnd, links);
		}
	}

	private static void outputToFile(int windowStart, int windowEnd, HashMap<Integer, ArrayList<Integer>> links) throws IOException {
		FileWriter f = new FileWriter("contactOutput/Tau Links" + windowStart + "-" + windowEnd + ".csv", false);
		f.write("Protein1,Protein2,AbsPos1,AbsPos2\n");
		StringBuilder output = new StringBuilder(26736);
		for (int i : links.keySet()) {
			ArrayList<Integer> bases = links.get(i);
			for (int resIndex : bases) {
				output.append("sp|P10636-8|TAU_HUMAN,sp|P10636-8|TAU_HUMAN," + (1 + i) + "," + (1 + resIndex) + "\n");
			}
		}
		f.write(output.toString().trim());
		f.close();
		System.out.println("tau links output to file");
	}
}
