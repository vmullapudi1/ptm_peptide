package utsw.joachimiak.vish;

import java.io.*;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.stream.Stream;

public class CSVParse {
	private static final String sourceCSV = "Human_Tau_PSM.csv";
	private static final String outputCSV = "HumanPSMOutput.csv";
	private static final String tau2N4R = "MAEPR" + "QEFEV" + "MEDHA" + "GTYGL" + "GDRKD" + "QGGYT" + "MHQDQ" + "EGDTD" + "AGLKE" + "SPLQT" //0-49
			+ "PTEDG" + "SEEPG" + "SETSD" + "AKSTP" + "TAEDV" + "TAPLV" + "DEGAP" + "GKQAA" + "AQPHT" + "EIPEG" //50-99
			+ "TTAEE" + "AGIGD" + "TPSLE" + "DEAAG" + "HVTQA" + "RMVSK" + "SKDGT" + "GSDDK" + "KAKGA" + "DGKTK" //100-149
			+ "IATPR" + "GAAPP" + "GQKGQ" + "ANATR" + "IPAKT" + "PPAPK" + "TPPSS" + "GEPPK" + "SGDRS" + "GYSSP" //150-199
			+ "GSPGT" + "PGSRS" + "RTPSL" + "PTPPT" + "REPKK" + "VAVVR" + "TPPKS" + "PSSAK" + "SRLQT" + "APVPM" //200-249
			+ "PDLKN" + "VKSKI" + "GSTEN" + "LKHQP" + "GGGKV" + "QIINK" + "KLDLS" + "NVQSK" + "CGSKD" + "NIKHV" //250-299
			+ "PGGGS" + "VQIVY" + "KPVDL" + "SKVTS" + "KCGSL" + "GNIHH" + "KPGGG" + "QVEVK" + "SEKLD" + "FKDRV" //300-349
			+ "QSKIG" + "SLDNI" + "THVPG" + "GGNKK" + "IETHK" + "LTFRE" + "NAKAK" + "TDHGA" + "EIVYK" + "SPVVS" //350-399
			+ "GDTSP" + "RHLSN" + "VSSTG" + "SIDMV" + "DSPQL" + "ATLAD" + "EVSAS" + "LAKQG" + "L";
	private static String[] columnHeaders;

	public static void main(String[] args) {
		Stream<String[]> fileData = null;

		try {
			fileData = ingestFile();
		} catch (IOException ex) {
			System.err.println("Error: File error " + System.getenv() + "\n" + ex.getMessage());
			System.exit(1);
		}
		System.out.println("File ingestion was a success.");


		//DEBUG
		//fileData.forEach(lineArr->System.out.println(Arrays.toString(lineArr)));

		//Get all of the phosphorylated peptides' lines
		Stream<String[]> phosphoPeptides = getPhosphoPeptides(fileData);

		//DEBUG
		//phosphoPeptides.forEach((lineArr->System.out.println(Arrays.toString(lineArr))));
	}

	/**
	 * Ingest file ,using comma as the delimiting value between items and producing a Stream<String[]> containing the file's data
	 *
	 * @return A stream from the file containing each line of the CSV file as a string array
	 * @throws IOException           if there is another error in reading/opening/working with the file (From BufferedReader)
	 * @throws FileNotFoundException if the file is not found in the local directory
	 */
	private static Stream<String[]> ingestFile() throws IOException {
		if (Files.exists((FileSystems.getDefault().getPath(sourceCSV)))) {
			System.out.println("attempting to open input file...");
		} else {
			System.err.println("File not found");
			throw new FileNotFoundException("File was not found");
		}
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(new FileInputStream(sourceCSV)));
		columnHeaders = inputReader.readLine().split(",");
		return inputReader.lines().map(dataLine -> dataLine.split(","));
	}

	/**
	 * @param fileLines The Stream<String[]> containing the file data (each String[] is a line, each string is a element in that line
	 * @return A Stream<String[]> containing all of the peptides indicated by the "Modifications" column to contain one or more phosphorylated sites
	 */
	private static Stream<String[]> getPhosphoPeptides(Stream<String[]> fileLines) {
		int modificationIndex = Arrays.asList(columnHeaders).indexOf("Modifications");
		return fileLines.filter(s -> s[modificationIndex].contains("Phospho"));
	}
}
