package utsw.joachimiak.vish;

import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.nio.file.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.stream.Collectors;

class CSVParse {
	public static void main(String[] args) {
		//Will contain the parsed raw data from the CSV of the MS output
		ArrayList<Peptide> fileData = null;
		//Will contain the data as peptide fragment objects grouped by fileID
		Map<String, List<Peptide>> peptidesByFile;
		//Get the input parameters
		final String sourceCSV = getInputFile();
		final String outputFileNameFormat = getOutputFormat();
		final String PROTEIN_SEQ = getProtein();
		final String outputFolderName = getOutputFolderName();

		//Read the file, complain if it doesn't work
		try {
			fileData = ingestFile(sourceCSV);
		} catch (IOException e) {
			System.err.println("error-file could not be read or found\n" + e);
			System.exit(1);
		}

		//Split the peptides by which file they came from
		peptidesByFile = fileData.stream()
				.collect(Collectors.groupingBy(Peptide::getFileID));

		//generate abundance and phosphorylation data for each file
		long tID = System.currentTimeMillis();//This tID is used as a unique folder name so that each run goes into it's own folder

		for (String fileID : peptidesByFile.keySet()) {
			//Get the peptide list of that file ID
			List<Peptide> p = peptidesByFile.get(fileID);

			//Align the peptide to Tau to generate its index in the tau isoform
			calcProteinPhosLocalizations(p, PROTEIN_SEQ);

			//Calculate the phosphorylated and unphosphorylated abundances for each residue in Tau
			Abundance[] abundance = generateResidueAbundances(p, PROTEIN_SEQ.length());

			//Attempt to print the abundance data to a file, or complain and exit if it doesn't work
			try {
				outputCSV(fileID, abundance, tID, sourceCSV, outputFileNameFormat, outputFolderName);
				System.out.println("Output " + fileID + " to folder " + tID);
			} catch (IOException e) {
				System.err.println("Error writing to file" + e);
				System.exit(1);
			}
		}
	}

	/**
	 * Ingest file, using comma as the delimiting value between items and producing an ArrayList<String[]> containing the file's data	 *
	 * @return An ArrayList<Peptide> from the file containing each line of the CSV file as a string array
	 * @throws IOException           if there is another error in reading/opening/working with the file (From BufferedReader)
	 * @throws FileNotFoundException if the file is not found in the local directory
	 */
	private static ArrayList<Peptide> ingestFile(String sourceCSV) throws IOException {
		//The list that wil contain the data from the file about each peptide
		ArrayList<Peptide> data=new ArrayList<>();
		//The first line of the file containing the column headers of the CSV file
		List<String> headerList;

		//Check to see if the source file exists and can be read, or else throw an exception
		if (Files.exists((FileSystems.getDefault().getPath(sourceCSV))) && new File(sourceCSV).canRead()) {
			System.out.println("attempting to open input file...");
		} else {
			throw new FileNotFoundException();
		}

		//Attempt to read the file. Throw an exception if it fails.
		try (BufferedReader inputReader = Files.newBufferedReader(Paths.get(sourceCSV))) {
			//Get the column headers for use in locating the fields of interest
			headerList = Arrays.asList(inputReader.readLine().split(","));

			//Find the locations of the fields of interest in the CSV output
			int fileIDIndex = headerList.indexOf("File ID");
			int sequenceIndex = headerList.indexOf("Annotated Sequence");
			int modIndex = headerList.indexOf("Modifications");
			int abundanceIndex = headerList.indexOf("Abundance: F1: Sample") == -1 ? headerList.indexOf("Precursor Abundance") : headerList.indexOf("Abundance: F1: Sample");
			if (abundanceIndex == -1) {
				abundanceIndex = headerList.indexOf("Abundance");
			}
			//Read the whole file, creating new peptides from the CSV data
			while (inputReader.ready()) {
				//String[] lineArr = inputReader.readLine().split(",");
				String[] lineArr = inputReader.readLine().replaceAll(",", ", ").split(",");
				if(lineArr.length==8){
					lineArr=Arrays.copyOf(lineArr,9);
					lineArr[8]="0";
				}
				for (int k = 0; k < lineArr.length; k++) {
					if (lineArr[k].equals(" ")) {
						lineArr[k] = "0";
					}
				}
				if(fileIDIndex==-1){
					data.add(new Peptide("FileID Null", lineArr[sequenceIndex], lineArr[modIndex], Double.parseDouble(lineArr[abundanceIndex])));
				}
				else {
					data.add(new Peptide(lineArr[fileIDIndex], lineArr[sequenceIndex], lineArr[modIndex], Double.parseDouble(lineArr[abundanceIndex])));
				}
			}
			return data;
		} catch (IOException e) {
			System.err.println("Error occurred while parsing input file\n" + e.getMessage());
			throw e;
		}
	}

	/**
	 * @param peptideList List of peptides for which to convert the peptide phosphorylation index to the protein
	 *                    residue number (zero-indexed)
	 * @param protein The sequence of the protein against which to index the peptide
	 */
	private static void calcProteinPhosLocalizations(@NotNull List<Peptide> peptideList, String protein) {
		for(Peptide p:peptideList){
			String seq = p.getSequence().toUpperCase();

			int seqProteinIndex = protein.indexOf(seq);
			//If the peptide isn't in the protein, ignore it and go to the next peptide
			if (seqProteinIndex == -1) {
				continue;
			}

			p.setPeptideProteinIndex(seqProteinIndex);
			String[] peptideSites=p.getPhosphorylations();
			if (peptideSites == null || peptideSites.length == 0 || peptideSites[0].equals("-1")) {
				continue;
			}

			int[] proteinPhosphorylationLocalizations = new int[peptideSites.length];
			Arrays.fill(proteinPhosphorylationLocalizations, -1);

			for (int i = 0; i < peptideSites.length; i++) {
				String siteStr = peptideSites[i];
				if (siteStr.matches(".*\\d+.*")) {
					int site = Integer.parseInt(siteStr.replaceAll("[\\D]", ""));
					proteinPhosphorylationLocalizations[i] = site - 1 + seqProteinIndex;
				} else {
					p.containsUnlocalizedPhosphorylation = true;
				}
			}
			p.setProteinPhosLocalizations(proteinPhosphorylationLocalizations);
		}
	}

	/**
	 * Calculates the modified (phosphorylated) and unmodified abundances for each phosphorylation locus given in the peptideList
	 *
	 * @param peptideList the list of peptides from a particular MS run/sample
	 * @return a 2d int array, where for every residue of tau a modified and unmodified abundance is stored
	 */
	private static Abundance[] generateResidueAbundances(List<Peptide> peptideList, int proteinLength) {

		Abundance[] residueAbundances = new Abundance[proteinLength];
		for (int i = 0; i < residueAbundances.length; i++) {
			residueAbundances[i] = new Abundance(0, 0);
		}
		//ArrayList<Peptide> notInTauFL = new ArrayList<>();
		for (Peptide p : peptideList) {
			int index = p.getIndexInProtein();
			if (index == -1) {
				//notInTauFL.add(p);
				continue;
			}
			double a = p.getAbundance();
			if (a < 0) {
				a = 0;
			}
			//add the abundance of the peptide to all the unmodified sites that the peptide covers
			for (int i = index; i < index + p.getSequence().length(); i++) {
				residueAbundances[i].total += a;
			}

			if (p.getProteinPhosLocalizations().length != 0) {
				//add the peptide abundance to the modified abundance of all the phosphorylated sites
				for (int i : p.getProteinPhosLocalizations()) {
					//If the localization doesn exist e.g. unlocalized phosphorylation, ignore it
					if (i == -1) {
						continue;
					}
					residueAbundances[i].phosphorylated += a;
				}
			}
		}
		return residueAbundances;
	}

	/**
	 * @param fileID     THe fileID of the source of the peptides
	 * @param abundances the array containing phosphorylated and total abundances of a residue
	 *                   in the PTM data
	 */
	private static void outputCSV(String fileID, Abundance[] abundances, long tID, String sourceCSV, String outputFormat, String outputFolderName) throws IOException {
		//Create the output file name and folder, using the system time as folder
		//and the date and time as filename modifiers
		LocalDateTime date = LocalDateTime.now();
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("kk-mm-MM-dd-YY");
		String dateString = date.format(formatter);
		//todo allow user picking of output folder
		Path path = Paths.get("output/" + outputFolderName);
		//create the folder for the current run
		Files.createDirectories(path);
		String sourceFile = sourceCSV.replaceAll("[/.<>:\"\\\\|\\-?*[\\s]]", "").strip();
		//The file writer. Creates new file instead of appending old ones, failing if the file already exists
		BufferedWriter w = Files.newBufferedWriter(
				new File(("output/" + outputFolderName + "/" + fileID + "_" + dateString + sourceFile + outputFormat).strip())
						.toPath(), StandardOpenOption.CREATE_NEW);

		StringBuilder outputBuffer = new StringBuilder(7000);
		//for residue-based phosphorylation data
		outputBuffer.append("Phosphorylation site, Phosphorylation Abundance, Total Abundance, Modification Proportion\n");
		double phos;
		double tot;

		for (int i = 0; i < abundances.length; i++) {
			phos = abundances[i].phosphorylated;
			tot = abundances[i].total;
			outputBuffer.append(i + "," + phos + "," + tot + "," + (phos / tot) + "\n");
		}

		w.write(outputBuffer.toString());
		w.close();
	}

	/*
	 //TODO Implement
	private static void outputPeptidePhosphorylation(String fileID, Map<Peptide,Abundance> p,long tID, String sourceCSV, String outputFormat)throws IOException{
		//Create the output file name and folder, using the system time as folder
		//and the date and time as filename modifiers
		LocalDateTime date = LocalDateTime.now();
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("kk-mm-MM-dd-YY");
		String dateString = date.format(formatter);
		//todo allow user picking of output folder
		Path path = Paths.get("output/" + tID);
		//create the folder for the current run
		Files.createDirectories(path);

		//The file writer. Creates new file instead of appending old ones, failing if the file already exists
		BufferedWriter w = Files.newBufferedWriter(new File("output/" + tID + "/" + fileID + "_" + dateString +
				sourceCSV + outputFormat).toPath(), StandardOpenOption.CREATE_NEW);

		StringBuilder outputBuffer = new StringBuilder(7000);
		//for peptide modification
		outputBuffer.append("Peptide Fragment,Tau Index,Fragment Length,Phosphorylation Abundance,Total Abundance,Modification Proportion\n");
		Map<String,List<Peptide>> fragments=p.stream().collect(Collectors.groupingBy(Peptide::getSequence));
		for(String s:fragments.keySet()){
			double phos=0;
			double tot=0;
			for(Peptide pep:fragments.get(s)){
				tot+=pep.getAbundance();
				//TODO: This doesn't fix the non-assigned localizations.
				if(pep.getPhosphorylations().length!=0){
					phos+=pep.getAbundance();
				}
			}
			int index=tau2N4R.indexOf(s.toUpperCase());
			if(index!=-1) {
				outputBuffer.append(s + ","+index+","+s.length()+"," + phos + "," + tot +","+ (phos/tot)+"\n");
			}
		}
		w.write(outputBuffer.toString());
		w.close();
	}
	*/

	//todo implement gui file picker or command line input of file
	private static String getInputFile() {
		Scanner s = new Scanner(System.in);
		System.out.println("Enter the name/path of the file you want to parse");
		String f = s.next();
		s.close();
		return f;
	}

	private static String getOutputFolderName() {
		Scanner s = new Scanner(System.in);
		System.out.println("Enter the folder name the out;ut should be stored in");
		String f = s.next();
		s.close();
		return f;
	}

	//todo:allow custom output format from user
	private static String getOutputFormat() {
		return "_Output.csv";
	}

	/*
		Todo: allow specifying the input protein to align the peptide fragments against
			Using keyboard input/system.in or file?
	 */
	private static String getProtein() {
		//Full length 2n4r tau isoform
		return
				"MAEPR" + "QEFEV" + "MEDHA" + "GTYGL" + "GDRKD" + "QGGYT" + "MHQDQ" + "EGDTD" + "AGLKE" + "SPLQT" //0-49
						+ "PTEDG" + "SEEPG" + "SETSD" + "AKSTP" + "TAEDV" + "TAPLV" + "DEGAP" + "GKQAA" + "AQPHT" + "EIPEG" //50-99
						+ "TTAEE" + "AGIGD" + "TPSLE" + "DEAAG" + "HVTQA" + "RMVSK" + "SKDGT" + "GSDDK" + "KAKGA" + "DGKTK" //100-149
						+ "IATPR" + "GAAPP" + "GQKGQ" + "ANATR" + "IPAKT" + "PPAPK" + "TPPSS" + "GEPPK" + "SGDRS" + "GYSSP" //150-199
						+ "GSPGT" + "PGSRS" + "RTPSL" + "PTPPT" + "REPKK" + "VAVVR" + "TPPKS" + "PSSAK" + "SRLQT" + "APVPM" //200-249
						+ "PDLKN" + "VKSKI" + "GSTEN" + "LKHQP" + "GGGKV" + "QIINK" + "KLDLS" + "NVQSK" + "CGSKD" + "NIKHV" //250-299
						+ "PGGGS" + "VQIVY" + "KPVDL" + "SKVTS" + "KCGSL" + "GNIHH" + "KPGGG" + "QVEVK" + "SEKLD" + "FKDRV" //300-349
						+ "QSKIG" + "SLDNI" + "THVPG" + "GGNKK" + "IETHK" + "LTFRE" + "NAKAK" + "TDHGA" + "EIVYK" + "SPVVS" //350-399
						+ "GDTSP" + "RHLSN" + "VSSTG" + "SIDMV" + "DSPQL" + "ATLAD" + "EVSAS" + "LAKQG" + "L";
	}
}