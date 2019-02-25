package utsw.joachimiak.vish;

import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

class CSVParse {

	private static final String sourceCSV = "Human_Tau_PSM.csv";
	private static final String outputFileNameFormat = "HumanPSMOutput.csv";
	private static final String tau2N4R =
			"MAEPR" + "QEFEV" + "MEDHA" + "GTYGL" + "GDRKD" + "QGGYT" + "MHQDQ" + "EGDTD" + "AGLKE" + "SPLQT" //0-49
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
		//fileData.forEach(System.out::println);

		//Split the peptides by which file they came from
		peptidesByFile = fileData.stream().collect(Collectors.groupingBy(Peptide::getFileID));

		//@DEBUG Splitting by file ID
		//System.out.println(peptidesByFile.toString().replaceAll(",","\n"));

		//generate abundance and phosphorylation data for each file
		for(String s:peptidesByFile.keySet()) {
			List<Peptide> p = peptidesByFile.get(s);
			calcPeptideTauPhosLocalizations(p);
			double[][] abundance = generateAbundances(p);
		}
	}

	/**
	 * Ingest file, using comma as the delimiting value between items and producing an ArrayList<String[]> containing the file's data	 *
	 * @return An ArrayList<String[]> from the file containing each line of the CSV file as a string array
	 * @throws IOException           if there is another error in reading/opening/working with the file (From BufferedReader)
	 * @throws FileNotFoundException if the file is not found in the local directory
	 */
	private static ArrayList<Peptide> ingestFile() throws IOException {
		ArrayList<Peptide> data=new ArrayList<>();
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
		int fileIDIndex = headerList.indexOf("File ID");
		int sequenceIndex = headerList.indexOf("Annotated Sequence");
		int modIndex = headerList.indexOf("Modifications");
		int abundanceIndex=headerList.indexOf("Precursor Abundance");

		while(inputReader.ready()) {
			String[] lineArr = inputReader.readLine().replaceAll("[,]{2}",",-1,").split(",");
			data.add(new Peptide(lineArr[fileIDIndex], lineArr[sequenceIndex], lineArr[modIndex], Double.parseDouble(lineArr[abundanceIndex])));
		}
		return data;
	}

	/**
	 * @param peptideList List of peptides for which to convert the peptide phosphorylation index to the 2N4R Tau
	 *                    residue number (zero-indexed)
	 */
	private static void calcPeptideTauPhosLocalizations(@NotNull List<Peptide> peptideList) {
		for(Peptide p:peptideList){
			String[] peptideSites=p.getPhosphorylations();
			if (peptideSites.length == 0) {
				continue;
			}
			if (peptideSites[0].equals("-1")){
				continue;
			}

			String seq = p.getSequence().toUpperCase();
			int peptideTauIndex=tau2N4R.indexOf(seq);
			if (peptideTauIndex == -1) {
				//System.out.println(p.getSequence());
			}
			p.setTauIndex(peptideTauIndex);
			int[] tauLocal = new int[peptideSites.length];

			for (int i = 0; i < peptideSites.length; i++) {
				int site = Integer.parseInt(peptideSites[i].replaceAll("[\\D]", ""));
				tauLocal[i]=site-1+peptideTauIndex;
			}
			p.setTauPhosLocalization(tauLocal);
		}
	}

	/**
	 * Calculates the modified (phosphorylated) and unmodified abundances for each phosphorylation locus given in the peptideList
	 *
	 * @param peptideList the list of peptides from a particular MS run/sample
	 * @return a 2d int array, where for every residue of tau a modified and unmodified abundance is stored
	 */
	private static double[][] generateAbundances(List<Peptide> peptideList) {
		double[][] abundances = new double[tau2N4R.length()][2];
		ArrayList<Peptide> notInTauFL = new ArrayList<>();

		for (Peptide p : peptideList) {
			int index = p.getTauIndex();
			if (index == -1) {
				notInTauFL.add(p);
				continue;
			}
			double a = p.getAbundance();
			//add the abundance of the peptide to all the unmodified sites that the peptide covers
			for (int i = index; i < p.getSequence().length(); i++) {
				abundances[i][1] += a;
			}
			//add the peptide abundance to the modified abundance of all the phophorylated sites
			for (int i : p.getTauPhosLocalization()) {
				//if no phosphorylation, break
				if (i == -1) {
					break;
				}
				abundances[i][0] += a;
			}
		}

		return abundances;
	}
}

