package utsw.joachimiak.vish;

import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.nio.file.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

class CSVParse {
	private CSVParse() {
	}

	public static void main(String[] args) {
		//Will contain the parsed raw data from the CSV of the MS output
		ArrayList<Fragment> fileData = null;
		//Will contain the data as peptide fragment objects grouped by fileID
		Map<String, List<Fragment>> peptidesByFile;

		//Get the input parameters
		Scanner s = new Scanner(System.in);
		String sourceCSV = CSVParse.getInputFile(s);
		String outputFileNameFormat = CSVParse.getOutputFormat();
		String PROTEIN_SEQ = CSVParse.getProtein();
		String outputFolderName = CSVParse.getOutputFolderName(s);
		boolean byFragment = CSVParse.queryPeptideOrResidueAnalysis(s);
		boolean isPeptideFile = isPeptideFormat(s);
		outputFileNameFormat = byFragment ? "peptide_analysis" + outputFileNameFormat : "residue_analysis" +
				outputFileNameFormat;
		s.close();
		final long tID = System.currentTimeMillis();//This tID is used as a unique folder name so that each run goes into
		// it's own folder

		//Read the file, complain if it doesn't work
		try {
			fileData = CSVParse.ingestFile(sourceCSV, isPeptideFile);
		} catch (IOException e) {
			System.err.println("error-file could not be read or found\n" + e);
			System.exit(1);
		}
		//Split the peptides by which file they came from
		peptidesByFile = fileData.stream()
				.collect(Collectors.groupingBy(Fragment::getFileID));

		//generate abundance and phosphorylation data for each file
		for (Map.Entry<String, List<Fragment>> stringListEntry : peptidesByFile.entrySet()) {
			//Get the peptide list of that file ID
			List<Fragment> p = stringListEntry.getValue();

			//Calculate residue modification analysis
			//otherwise calculate the peptide modification analysis
			if (byFragment) {
				try {
					CSVParse.outputFragmentCSV(stringListEntry.getKey(), CSVParse.generateFragmentPhosAbundances(p), tID, sourceCSV,
							outputFileNameFormat, outputFolderName);
					System.out.println("Output " + stringListEntry.getKey() + " to folder " + "output/" + outputFolderName + tID);
				} catch (IOException e) {
					System.err.println("Error writing to file " + e);
					System.exit(1);
				}
			} else {
				//Align the peptide to Tau to generate its index in the tau isoform
				CSVParse.calcProteinPhosLocalizations(p, PROTEIN_SEQ);

				//Calculate the phosphorylated and unphosphorylated abundances for each residue in Tau
				Abundance[] residuePhosAbundances = CSVParse.generateResiduePhosAbundances(p, PROTEIN_SEQ.length());

				//Attempt to print the abundance data to a file, or complain and exit if it doesn't work
				try {
					CSVParse.outputResidueCSV(stringListEntry.getKey(), residuePhosAbundances, tID, sourceCSV, outputFileNameFormat,
							outputFolderName);
					System.out.println("Output " + stringListEntry.getKey() + " to folder " + "output/" +
							outputFolderName + tID);
				} catch (IOException e) {
					System.err.println("Error writing to file " + e);
					System.exit(1);
				}
			}
			}
		}

	/**
	 * Determines whether or not the uses wan't to perform by residue or by peptide analysis
	 * TODO convert to CLI
	 */
	private static boolean queryPeptideOrResidueAnalysis(Scanner s) {
		System.out.println("Should the output format be by peptide fragment instead of by residue?\n[Y/N]");
		String ans = s.next();
		if ("Y".equalsIgnoreCase(ans)) {
			return true;
		} else if ("N".equalsIgnoreCase(ans)) {
			return false;
		}
		System.err.println("Invalid input " + ans);
		System.exit(2);
		return false;
	}

	private static boolean isPeptideFormat(Scanner s) {
		System.out.println("Is this file a peptide file, without a separate fileID column, instead denoting different " +
				"samples using \"Abundance: [ID}\"?\n[Y/N]?");
		String ans = s.next();
		if ("Y".equalsIgnoreCase(ans)) {
			return true;
		} else if ("N".equalsIgnoreCase(ans)) {
			return false;
		}
		System.err.println("Invalid input " + ans);
		System.exit(2);
		return false;
	}

	/**
	 * Ingest file, using comma as the delimiting value between items and producing an ArrayList<String[]> containing
	 * the file's data
	 * @return An ArrayList of fragments from the file containing each line of the CSV file as a string array
	 * @throws IOException           if there is another error in reading/opening/working with the file
	 * (From BufferedReader)
	 * @throws FileNotFoundException if the file is not found in the local directory
	 */
	private static ArrayList<Fragment> ingestFile(String sourceCSV, boolean isPeptideFile) throws FileNotFoundException, IOException {
		//The list that wil contain the data from the file about each peptide
		ArrayList<Fragment> data = new ArrayList<>(3000);
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
			int abundanceIndex = headerList.indexOf("Abundance: F1: Sample") == -1 ?
					headerList.indexOf("Precursor Abundance") : headerList.indexOf("Abundance: F1: Sample");
			if (abundanceIndex == -1) {
				abundanceIndex = headerList.indexOf("Abundance");
			}
			HashMap<Integer, String> abundanceIndexWithFileID = new HashMap<>();
			if (isPeptideFile) {
				Matcher abundancePattern = Pattern.compile("Abundance: (\\w)*").matcher("");
				for (int i = 0; i < headerList.size(); i++) {
					abundancePattern.reset(headerList.get(i));
					if (abundancePattern.matches()) {
						String fileID = abundancePattern.group(1);
						abundanceIndexWithFileID.put(i, fileID);
					}
				}
			}
			//todo: for peptide file formats grab the fileID from the Abundance:F# column or Abundance: #W
			//Read the whole file, creating new peptides from the CSV data
			Pattern delimiter = Pattern.compile(",");
			Matcher m = delimiter.matcher("");
			while (inputReader.ready()) {
				String line = inputReader.readLine();
				m.reset(line);
				String[] lineArr = m.replaceAll(", ").split(",");
				if (lineArr.length < headerList.size()) {
					int lineArrlen = lineArr.length;
					lineArr = Arrays.copyOf(lineArr, headerList.size());
					Arrays.fill(lineArr, lineArrlen - 1, lineArr.length - 1, "0");
				}
				for (int k = 0; k < lineArr.length; k++) {
					if (" ".equals(lineArr[k])) {
						lineArr[k] = "0";
					}
				}
				if (isPeptideFile) {
					for (Map.Entry<Integer, String> fileEntry : abundanceIndexWithFileID.entrySet()) {
						if (lineArr[fileEntry.getKey()].equals("0")) {
							continue;
						}
						data.add(new Fragment(fileEntry.getValue(), lineArr[fileEntry.getKey()], lineArr[modIndex],
								Double.parseDouble(lineArr[fileEntry.getKey()])));
					}
				} else if (fileIDIndex == -1) {
					data.add(new Fragment("FileID F1", lineArr[sequenceIndex], lineArr[modIndex],
							Double.parseDouble(lineArr[abundanceIndex])));
				}
				else {
					if (lineArr[sequenceIndex].equals("0")) {
						continue;
					}
					data.add(new Fragment(lineArr[fileIDIndex], lineArr[sequenceIndex], lineArr[modIndex],
							Double.parseDouble(lineArr[abundanceIndex])));
				}
			}
			return data;
		} catch (IOException e) {
			System.err.println("Error occurred while parsing input file\n" + e.getMessage());
			throw e;
		}
	}

	/**
	 * @param fragmentList List of peptide fragments for which to convert the peptide phosphorylation index to the protein
	 *                    residue number (zero-indexed)
	 * @param protein The sequence of the protein against which to index the peptide
	 */
	private static void calcProteinPhosLocalizations(@NotNull List<Fragment> fragmentList, String protein) {
		Pattern phosLocalizer = Pattern.compile(".*\\d+.*");
		Pattern digitIsolator = Pattern.compile("[\\D]");
		Matcher finder = phosLocalizer.matcher("");
		Matcher replacer = digitIsolator.matcher("");
		for (Fragment p : fragmentList) {
			String seq = p.getSequence().toUpperCase();
			int seqProteinIndex = protein.indexOf(seq);

			//If the peptide isn't in the protein, ignore it and go to the next peptide
			if (seqProteinIndex == -1) {
				continue;
			}

			p.setPeptideProteinIndex(seqProteinIndex);
			String[] peptideSites = p.getPhosphorylations();
			if (peptideSites == null || peptideSites.length == 0 || "-1".equals(peptideSites[0])) {
				continue;
			}

			int[] proteinPhosphorylationLocalizations = new int[peptideSites.length];
			Arrays.fill(proteinPhosphorylationLocalizations, -1);

			for (int i = 0; i < peptideSites.length; i++) {
				String siteStr = peptideSites[i];
				finder.reset(siteStr);
				if (finder.find()) {
					replacer.reset(siteStr);
					int site = Integer.parseInt(replacer.replaceAll(""));
					proteinPhosphorylationLocalizations[i] = site - 1 + seqProteinIndex;
				} else {
					p.containsUnlocalizedPhosphorylation = true;
				}
			}
			p.setProteinPhosLocalizations(proteinPhosphorylationLocalizations);
		}
	}

	/**
	 * Calculates the modified (phosphorylated) and unmodified abundances for each phosphorylation locus given in the
	 * fragmentList
	 * @param fragmentList the list of peptide fragments from a particular MS run/sample
	 * @return a 2d int array, where for every residue of tau a modified and unmodified abundance is stored
	 */
	private static Abundance[] generateResiduePhosAbundances(List<Fragment> fragmentList, int proteinLength) {
		Abundance[] residueAbundances = new Abundance[proteinLength];
		for (int i = 0; i < residueAbundances.length; i++) {
			residueAbundances[i] = new Abundance(0, 0);
		}

		for (Fragment p : fragmentList) {
			int index = p.getIndexInProtein();
			//If the fragment isn't in the protein, ignore it
			if (index == -1) {
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
					//If the localization doesn't exist e.g. unlocalized phosphorylation, ignore it
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
	 * Combines all fragments in the input that have the same amino acid sequence together and sums their abundances,
	 * creating an accumulated fragment that contains to total phosporylated and unphosphorylated incidence of the fragment,
	 * and notes whether that fragment contains unlocalized phosphorylations
	 * @param fragmentList The fragments to be analyzed
	 * @return A HashMap<Fragment, Abundance> mapping each combined fragment to its abundance
	 */
	private static HashMap<Fragment, Abundance> generateFragmentPhosAbundances(List<Fragment> fragmentList) {
		Map<String, List<Fragment>> fragmentsBySeq = fragmentList.stream().collect(Collectors.groupingBy(f ->
				f.getSequence().toUpperCase()));
		HashMap<Fragment, Abundance> ans = new HashMap<>(fragmentsBySeq.keySet().size());
		String fileID;
		int protIndex;

		for (Map.Entry<String, List<Fragment>> stringListEntry : fragmentsBySeq.entrySet()) {
			List<Fragment> sameSeq = stringListEntry.getValue();
			fileID = sameSeq.get(0).getFileID();
			protIndex = CSVParse.getProtein().indexOf(stringListEntry.getKey());
			Abundance combinedAbundance = new Abundance(0, 0);
			Set<String> phosphorylations = new HashSet<>(10);
			Set<Integer> proteinPhosSites = new HashSet<>(10);
			boolean unlocalized = false;
			//if the fragment doesn't align to the protein, ignore it and continue
			if (protIndex == -1) {
				continue;
			}
			for (Fragment f : sameSeq) {
				if (f.getPhosphorylations().length > 0) {
					combinedAbundance.phosphorylated += f.getAbundance();
					phosphorylations.addAll(Arrays.asList(f.getPhosphorylations()));
					proteinPhosSites.addAll(IntStream.of(f.getProteinPhosLocalizations()).boxed().collect(Collectors.toList()));
				}
				combinedAbundance.total += f.getAbundance();
				if (f.containsUnlocalizedPhosphorylation) {
					unlocalized = true;
				}
			}

			ans.put(new Fragment(fileID, stringListEntry.getKey(), phosphorylations.toArray(new String[0]),
					proteinPhosSites.stream().filter(i -> i != null && i >= 0)
							.mapToInt(Integer::intValue)
							.toArray(), protIndex, unlocalized), combinedAbundance);
		}
		return ans;
	}

	/**
	 *	Outputs the data created by the residue-based phosphorylation analysis to a csv file
	 * @param fileID-The file id of the residues being output
	 * @param abundances An array containing the abundance data for each residue in the protein, indexed to the residue of
	 *                   the protein (zero-indexed)
	 * @param tID An id number used to temporally segregate runs
	 * @param sourceCSV the filename of the source CSV file
	 * @param outputFormat a string to modify the output filename with
	 * @param outputFolderName the folder name to output the run to. Is added to the tID to produce a more unique folder
	 *                         name.
	 * @throws IOException When an error occurs in opening the file to be written to or in writing the file
	 */
	private static void outputResidueCSV(String fileID, Abundance[] abundances, long tID,
										 String sourceCSV, String outputFormat,
										 String outputFolderName) throws IOException {
		//Create the output file name and folder, using the system time as folder
		//and the date and time as filename modifiers
		LocalDateTime date = LocalDateTime.now();
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("kk-mm-MM-dd-YY");
		String dateString = date.format(formatter);
		//todo allow user picking of output folder
		Path path = Paths.get("output/" + outputFolderName + tID);
		//create the folder for the current run
		Files.createDirectories(path);
		String sourceFile = sourceCSV.replaceAll("[/.<>:\"\\\\|\\-?*[\\s]]", "").strip();
		//The file writer. Creates new file instead of appending old ones, failing if the file already exists
		BufferedWriter w = Files.newBufferedWriter(
				new File(("output/" + outputFolderName + tID + "/" + fileID + "_" + dateString + sourceFile + outputFormat).strip())
						.toPath(), StandardOpenOption.CREATE_NEW);

		StringBuilder outputBuffer = new StringBuilder(7000);
		//for residue-based phosphorylation data
		outputBuffer.append("Phosphorylation site, Phosphorylation Abundance, Total Abundance, Modification Proportion\n");
		double phos;
		double tot;

		for (int i = 0; i < abundances.length; i++) {
			phos = abundances[i].phosphorylated;
			tot = abundances[i].total;
			String out = i + "," + phos + "," + tot + "," + (phos / tot) + "\n";
			outputBuffer.append(out);
		}

		w.write(outputBuffer.toString());
		w.close();
	}

	/**
	 * Outputs the data created by the protein-fragment based phosphorylation analysis to a CSV file
	 *
	 * @param fileID              he file id of the residues being output
	 * @param combinedFragments-A HashMap of the accumulated Fragment objects containing the aggregated data from all
	 *                            fragments
	 *                            of that fileId with the same protein sequence, along with the associated abundance
	 *                            values.
	 * @param tID                 An id number used to temporally segregate runs
	 * @param sourceCSV           the filename of the source CSV file
	 * @param outputFormat        a string to modify the output filename with
	 * @param outputFolderName    the folder name to output the run to. Is added to the tID to produce a more unique
	 *                            folder name.
	 * @throws IOException When an error occurs in opening the file to be written to or in writing the file
	 */
	private static void outputFragmentCSV(String fileID, HashMap<Fragment, Abundance> combinedFragments,
										  long tID, String sourceCSV, String outputFormat,
										  String outputFolderName) throws IOException {
		//Create the output file name and folder, using the system time as folder
		//and the date and time as filename modifiers
		LocalDateTime date = LocalDateTime.now();
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("kk-mm-MM-dd-YY");
		String dateString = date.format(formatter);
		//todo allow user picking of output folder
		Path path = Paths.get("output/" + outputFolderName + tID);
		//create the folder for the current run
		Files.createDirectories(path);
		String sourceFile = sourceCSV.replaceAll("[/.<>:\"\\\\|\\-?*[\\s]]", "").strip();
		//The file writer. Creates new file instead of appending old ones, failing if the file already exists
		BufferedWriter w = Files.newBufferedWriter(
				new File(("output/" + outputFolderName + tID + "/" + fileID + "_" + dateString + sourceFile + outputFormat).strip())
						.toPath(), StandardOpenOption.CREATE_NEW);

		StringBuilder outputBuffer = new StringBuilder(7000);
		outputBuffer.append("Fragment,Start Index,End Index,Modified Abundance,Unmodified Abundance,Modification Proportion\n");

		ArrayList<Fragment> sortedFragments = new ArrayList<>(combinedFragments.keySet());
		sortedFragments.sort(Comparator.comparing(Fragment::getIndexInProtein));

		for (Fragment f : sortedFragments) {
			double phos = combinedFragments.get(f).phosphorylated;
			double tot = +combinedFragments.get(f).total;
			String out = f.getSequence() + "," + f.getIndexInProtein() + "," +
					(f.getIndexInProtein() + f.getSequence().length()) + "," + phos + "," + tot + "," + (phos / tot) + "\n";
			outputBuffer.append(out);

		}
		w.write(outputBuffer.toString());
		w.close();
	}

	/**
	 * Gets the user's desired input's filename
	 *
	 * @param s The scanner object from which to retrieve the user's input
	 * @return a String containing the file path input by the user
	 */
	//todo implement command line input of file or a config file
	private static String getInputFile(@NotNull Scanner s) {
		System.out.println("Enter the name/path of the file you want to parse");
		return s.nextLine();
	}

	/**
	 * Get's the user's desired output folder name
	 * @param s the Scanner object from which to query the user's desired output folder name
	 * @return a String containing the desired output folder name
	 */
	//todo implement config file or commandline arg
	private static String getOutputFolderName(Scanner s) {
		System.out.println("Enter the folder name the output should be stored in");
		return s.next();
	}

	/**
	 * Gets the desired output string to append to the output filename
	 * @return A string containing the postfix for the output filename
	 */
	//todo:allow custom output format from user, either from CLI or config file
	private static String getOutputFormat() {
		return "_Output.csv";
	}


	//Todo: allow specifying the input protein to align the peptide fragments against
	// 	Using keyboard input/system.in or file/config file?

	/**
	 * @return The protein sequence against which to index peptides
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

		/*//1N4R P301S Tau
		return
						"MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPG" +
						"SETSDAKSTPTAEAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKT" +
						"KIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSR" +
						"SRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQ" +
						"PGGGKVQIINKKLDLSNVQSKCGSKDNIKHVSGGGSVQIVYKPVDLSKVTSKCGSLGNIH" +
						"HKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHG" +
						"AEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL";*/
	}
}