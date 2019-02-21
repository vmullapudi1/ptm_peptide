package utsw.joachimiak.vish;

import java.io.*;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;

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
		ArrayList<Peptide> fileData = null;
		HashSet<String> fileID_list;
		Map<String, List<Peptide>> peptidesByFile;
		try {
			fileData = ingestFile();
		} catch (IOException ex) {
			System.err.println("Error: File error " + System.getenv() + "\n" + ex.getMessage());
			System.exit(1);
		}
		System.out.println("File ingestion was a success.");

		//@DEBUG Column header array generation
		//Arrays.stream(columnHeaders).forEach(System.out::println);

		//@DEBUG file ingestion
		fileData.forEach(System.out::println);

		//Split the peptides by which file they came from
		peptidesByFile = fileData.stream().collect(Collectors.groupingBy(Peptide::getFileID));

		//@DEBUG Splitting by file ID
		//System.out.println(peptidesByFile.toString().replaceAll(",","\n"));
	}


	/**
	 * Ingest file, using comma as the delimiting value between items and producing an ArrayList<String[]> containing the file's data	 *
	 * @return An ArrayList<String[]> from the file containing each line of the CSV file as a string array
	 * @throws IOException           if there is another error in reading/opening/working with the file (From BufferedReader)
	 * @throws FileNotFoundException if the file is not found in the local directory
	 */
	private static ArrayList<Peptide> ingestFile() throws IOException {
		ArrayList<Peptide> data;
		List<String> headerList;
		if (Files.exists((FileSystems.getDefault().getPath(sourceCSV)))) {
			System.out.println("attempting to open input file...");
		} else {
			System.err.println("File not found");
			throw new FileNotFoundException("File was not found");
		}
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(new FileInputStream(sourceCSV)));
		columnHeaders = inputReader.readLine().split(",");
		headerList = Arrays.asList(columnHeaders);
		int fileIDindex = headerList.indexOf("File ID");
		int sequenceIndex = headerList.indexOf("Annotated Sequence");
		int phosphoIndex = headerList.indexOf("Modifications");
		data = inputReader.lines().map(line -> line.split(","))
				.map(lineArr -> new Peptide(lineArr[fileIDindex], lineArr[sequenceIndex], lineArr[phosphoIndex].split(";")))
				.collect(Collectors.toCollection(ArrayList::new));
		return data;
	}

}

