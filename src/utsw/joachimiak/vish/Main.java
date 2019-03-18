package utsw.joachimiak.vish;

import org.apache.poi.EncryptedDocumentException;
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

class Main {


	private static ArrayList<Double> abundances;
	private static ArrayList<ArrayList<Double>> abundancesAD;

	private static ArrayList<ArrayList<String>> extractFile(String data) throws EncryptedDocumentException, IOException {
		@SuppressWarnings("resource")
		Workbook workbook = WorkbookFactory.create(new File(data));
		ArrayList<ArrayList<String>> excelReadout = new ArrayList<ArrayList<String>>();
		// Retrieving the number of sheets in the Workbook
		System.out.println("Workbook has " + workbook.getNumberOfSheets() + " Sheets : ");
		Iterator<Sheet> sheetIterator = workbook.sheetIterator();
		System.out.println("Retrieving Sheets using Iterator");
		while (sheetIterator.hasNext()) {
			Sheet sheet = sheetIterator.next();
			System.out.println("=> " + sheet.getSheetName());
		}
		// Getting the Sheet at index zero
		Sheet sheet = workbook.getSheetAt(0);

		// Create a DataFormatter to format and get each cell's value as String
		DataFormatter dataFormatter = new DataFormatter();

		// 1. You can obtain a rowIterator and columnIterator and iterate over them
		//System.out.println("\n\nIterating over Rows and Columns using Iterator\n");
		Iterator<Row> rowIterator = sheet.rowIterator();

		Row rowZero = rowIterator.next();
		//System.out.println(row.getCell(0).toString());
		// Now let's iterate over the columns of the current row
		Iterator<Cell> titles = rowZero.cellIterator();
		//int numberOfCol = 0;
		int annSeq = 0;
		int mods = 0;
		int id = 0;
		int abund = 0;
		int l = 0;
		while (titles.hasNext()) {
			Cell cell = titles.next();
			//System.out.println(cell.getStringCellValue());
			//System.out.println(titles[l]);
			if (cell.getStringCellValue().contains("Annotated Sequence")) {
				annSeq = l;
			}
			if (cell.getStringCellValue().contains("Modifications") && !cell.getStringCellValue().contains("in Master Proteins")) {
				mods = l;
			}
			//for the mice data
			if (cell.getStringCellValue().contains("Abundance:557622")) {
				abund = l;
			}
				/*if (cell.getStringCellValue().contains("File ID"))	{
					id = l;
				}*/
			l++;
		}
		//System.out.println(l);
		abundances = new ArrayList<Double>();
		abundancesAD = new ArrayList<ArrayList<Double>>();
		abundancesAD.add(new ArrayList<Double>());
		abundancesAD.add(new ArrayList<Double>());
		abundancesAD.add(new ArrayList<Double>());
		abundancesAD.add(new ArrayList<Double>());
		abundancesAD.add(new ArrayList<Double>());
		abundancesAD.add(new ArrayList<Double>());
		while (rowIterator.hasNext()) {
			Row row = rowIterator.next();
			ArrayList<String> shortlist = new ArrayList<String>();
			String frag = row.getCell(annSeq).toString(); //reads out just the 3rd column
			String[] newfrag = frag.split("]");
			String fragseq = newfrag[1];
			int x = fragseq.indexOf("[");
			String peptideseq = fragseq.substring(1, x - 1);
			String pep = peptideseq.toUpperCase();
			shortlist.add(pep); //grabs only the aa sequence desired
			//String idTag = row.getCell(id).toString();
			//shortlist.add(idTag);
			String modification = row.getCell(mods).toString(); //reads of the 4th column
			String toAdd = null;
			if (modification.contains("Phospho")) {
				String[] mod = modification.replaceAll("[(Phospho)]", "").split("[);]");//modification.split("o ");
				String Mod = mod[1];//.split("]")[0];
				toAdd = Mod.substring(1);
			}
			shortlist.add(toAdd);
			String C1Abund = row.getCell(abund).toString();
			if (C1Abund.equals("")) {
				C1Abund = "0.0";
			}
			shortlist.add(C1Abund);

			excelReadout.add(shortlist);
		}

		return excelReadout;
	}

	//[TDHGAEIVYKSPVVSGDTSPR, S11; S15; S19, 662379.375]
	private static ArrayList<ArrayList<Integer>> findIndexes(String full, ArrayList<ArrayList<String>> week5) {
		ArrayList<ArrayList<Integer>> answer = new ArrayList<ArrayList<Integer>>();
		for (ArrayList<String> combo : week5) {
			String frag = combo.get(0);
			ArrayList<Integer> singlelist = new ArrayList<Integer>();
			int start = full.indexOf(frag);
			int end = start + frag.length() - 1;
			//String specificMod = combo.get(1);
			singlelist.add(start);
			singlelist.add(end);
			answer.add(singlelist);
			//System.out.println(start + " " + end);
		}
		return answer;
	}

	private static ArrayList<ArrayList<Integer>> findModSite(ArrayList<ArrayList<Integer>> fragmentLocation, ArrayList<ArrayList<String>> week5) {
		ArrayList<ArrayList<Integer>> answer = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < week5.size(); i++) {
			int start = fragmentLocation.get(i).get(0);
			String siteLoc = week5.get(i).get(1);
			ArrayList<Integer> toAdd = new ArrayList<Integer>();
			//System.out.println(siteLoc);
			//System.out.println(siteLoc);
			if (siteLoc == null) {

			} else if (!siteLoc.contains("/") && siteLoc.length() <= 3 && siteLoc.length() > 1) { //note: makes sure to filter out the ambiguous modifications
				int modLoc = start + Integer.parseInt(siteLoc.substring(1)) - 1;
				toAdd.add(modLoc);
			} else if (siteLoc.length() > 3) {
				String[] sites = siteLoc.split("; ");
				for (int k = 0; k < sites.length; k++) {
					if (sites[k].length() > 1 && !sites[k].contains("/")) {
						toAdd.add(start + Integer.parseInt(sites[k].substring(1)) - 1);
					}
				}
			}
			answer.add(toAdd);
			//System.out.println(toAdd);
		}
		return answer;
	}

	//this is going by aa, independent of fragments
	private static LinkedHashMap<String, ArrayList<Double>> calcAbundances(ArrayList<ArrayList<String>> week5, ArrayList<ArrayList<Integer>> modPlacementInTau, String full, ArrayList<ArrayList<Integer>> fragmentLocation, String type, Boolean method) {
		LinkedHashMap<String, ArrayList<Double>> answer = new LinkedHashMap<String, ArrayList<Double>>();
		for (int j = 0; j < modPlacementInTau.size(); j++) {        //looking at ever Tau aa
			for (int i : modPlacementInTau.get(j)) {
				answer.put(type + "  -" + i, new ArrayList<Double>());
				answer.put(type + " MOD -" + i, new ArrayList<Double>());
			}
		}

		for (int i = 0; i < full.length(); i++) {
			for (int j = 0; j < fragmentLocation.size(); j++) { //looking at every fragment
				ArrayList<Integer> aa = fragmentLocation.get(j);
				if (aa.get(0) <= i && i <= aa.get(1) && !aa.get(0).equals(-1) && !aa.get(1).equals(-1)) {            //looking to see if that Tau aa is in the fragment
					ArrayList<Integer> mod = modPlacementInTau.get(j);
					if (week5.get(j).size() == 3 && answer.containsKey(type + "  -" + i) || answer.containsKey(type + "~ -" + i)) {
						for (int site : mod) {
							if (site == i) { //looking to see if the mod site == current aa site
								ArrayList<Double> temp1 = answer.get(type + " MOD -" + i);
								temp1.add(Double.parseDouble(week5.get(j).get(2)));
								answer.put(type + " MOD -" + i, temp1);
								ArrayList<Double> temp2 = answer.get(type + "  -" + i);
								temp2.add(Double.parseDouble(week5.get(j).get(2)));
								answer.put(type + "  -" + i, temp2);
							}
						}
						ArrayList<Double> temp2 = answer.get(type + "  -" + i);
						temp2.add(Double.parseDouble(week5.get(j).get(2)));
						answer.put(type + "  -" + i, temp2);
					}
				}
			}
		}
		return answer;
		//[TDHGAEIVYKSPVVSGDTSPR, S11; S15; S19, 662379.375]
	}

	//TODO fix for peptide/PSM input file format
	//TODO change to separate phosphorylation localization based on fileID
	public static void main(String[] args) throws EncryptedDocumentException, IOException {
		//Tau 2N4R isoform
		String tauIsoformF =
				"MAEPR" + "QEFEV" + "MEDHA" + "GTYGL" + "GDRKD" + "QGGYT" + "MHQDQ" + "EGDTD" + "AGLKE" + "SPLQT" //0-49
						+ "PTEDG" + "SEEPG" + "SETSD" + "AKSTP" + "TAEDV" + "TAPLV" + "DEGAP" + "GKQAA" + "AQPHT" + "EIPEG" //50-99
						+ "TTAEE" + "AGIGD" + "TPSLE" + "DEAAG" + "HVTQA" + "RMVSK" + "SKDGT" + "GSDDK" + "KAKGA" + "DGKTK" //100-149
						+ "IATPR" + "GAAPP" + "GQKGQ" + "ANATR" + "IPAKT" + "PPAPK" + "TPPSS" + "GEPPK" + "SGDRS" + "GYSSP" //150-199
						+ "GSPGT" + "PGSRS" + "RTPSL" + "PTPPT" + "REPKK" + "VAVVR" + "TPPKS" + "PSSAK" + "SRLQT" + "APVPM" //200-249
						+ "PDLKN" + "VKSKI" + "GSTEN" + "LKHQP" + "GGGKV" + "QIINK" + "KLDLS" + "NVQSK" + "CGSKD" + "NIKHV" //250-299
						+ "PGGGS" + "VQIVY" + "KPVDL" + "SKVTS" + "KCGSL" + "GNIHH" + "KPGGG" + "QVEVK" + "SEKLD" + "FKDRV" //300-349
						+ "QSKIG" + "SLDNI" + "THVPG" + "GGNKK" + "IETHK" + "LTFRE" + "NAKAK" + "TDHGA" + "EIVYK" + "SPVVS" //350-399
						+ "GDTSP" + "RHLSN" + "VSSTG" + "SIDMV" + "DSPQL" + "ATLAD" + "EVSAS" + "LAKQG" + "L";


		//true = by aa
		//false = by fragment
		Boolean method = true;

		String fileName = "PRJ-382-LFQ-Tau All PSMs_20180629";
		String outputFileName = "phosSitesOutput";
		String type = "557622";

		System.out.println("Reading file...");

		ArrayList<ArrayList<String>> week5 = extractFile(fileName + ".xlsx"); //[frag, sight of mod, abundance]
        /*for (ArrayList<String> tempp : week5) {
        	System.out.println(tempp);
        }*/
		ArrayList<ArrayList<Integer>> fragmentLocation = findIndexes(tauIsoformF, week5); //[start, end] for each fragment
		ArrayList<ArrayList<Integer>> modPlacementInTau = findModSite(fragmentLocation, week5); //[placement of modification(s)] for each fragment
		LinkedHashMap<String, ArrayList<Double>> totalAbundances = calcAbundances(week5, modPlacementInTau, tauIsoformF, fragmentLocation, type, method);

		System.out.println("Finished compiling");
		Workbook workbook = new XSSFWorkbook(); // new HSSFWorkbook() for generating `.xls` file
		CreationHelper createHelper = workbook.getCreationHelper();

		// Create a Sheet
		Sheet sheet = workbook.createSheet("Data Set");

		// Create a Font for styling header cells
		Font headerFont = workbook.createFont();

		// Create a CellStyle with the font
		CellStyle headerCellStyle = workbook.createCellStyle();
		headerCellStyle.setFont(headerFont);

		// Create a Row
		Row headerRow = sheet.createRow(0);

		//if (fragOption == 0) {
		ArrayList<ArrayList<Object>> list = new ArrayList<ArrayList<Object>>();
		int counter = 0;
		for (ArrayList<Integer> ends : fragmentLocation) {
			if (!ends.contains(-1)) {
				ArrayList<Object> temp = new ArrayList<Object>();
				temp.add(week5.get(counter).get(0));
				temp.addAll(ends);
				list.add(temp);
			}

			counter++;
		}
		System.out.println("Writing into file: " + outputFileName);

		//generates tables for coverage maps

		//make "all fragments" coverage map
		Cell cellF = headerRow.createCell(0);
		Cell cellX = headerRow.createCell(1);
		Cell cellR = headerRow.createCell(2);
		Cell cellY = headerRow.createCell(3);
		cellX.setCellValue("Start");
		cellR.setCellValue("Range");
		cellY.setCellValue("End");
		cellF.setCellValue("All Fragments in file:");
		int rowNum = 1;
		Set<ArrayList<Object>> set = new HashSet<ArrayList<Object>>(list);
		for (ArrayList<Object> thing : set) {
			if (!thing.get(1).equals(-1)) {
				Row currRow = sheet.createRow(rowNum);
				currRow.createCell(0).setCellValue((String) thing.get(0));
				currRow.createCell(1).setCellValue((int) thing.get(1));
				currRow.createCell(2).setCellValue((int) thing.get(2) - (int) thing.get(1));
				currRow.createCell(3).setCellValue((int) thing.get(2));
				rowNum++;
			}
		}

		//make "modified fragments" coverage map
		Cell cellF2 = headerRow.createCell(5);
		Cell cellX2 = headerRow.createCell(6);
		Cell cellR2 = headerRow.createCell(7);
		Cell cellY2 = headerRow.createCell(8);
		cellX2.setCellValue("Start");
		cellR2.setCellValue("Range");
		cellY2.setCellValue("End");
		cellF2.setCellValue("Modified Fragments in file:");
		int rowNumber = 1;
		counter = 0;
		for (ArrayList<String> info : week5) {
			ArrayList<Object> toRemove = new ArrayList<Object>();
			toRemove.add(week5.get(counter).get(0));
			toRemove.add(fragmentLocation.get(counter).get(0));
			toRemove.add(fragmentLocation.get(counter).get(1));
			if (week5.get(counter).get(2) == null) {
				list.remove(toRemove);
			}
			counter++;
		}
		Set<ArrayList<Object>> newset = new HashSet<ArrayList<Object>>(list);
		for (ArrayList<Object> thing : newset) {
			Row currRow = sheet.createRow(rowNumber);
			currRow.createCell(5).setCellValue((String) thing.get(0));
			currRow.createCell(6).setCellValue((int) thing.get(1));
			currRow.createCell(7).setCellValue((int) thing.get(2) - (int) thing.get(1));
			currRow.createCell(8).setCellValue((int) thing.get(2));
			rowNumber++;
		}


		int timer = 0;
		for (Entry<String, ArrayList<Double>> entry : totalAbundances.entrySet()) {
			String aa = entry.getKey();
			ArrayList<Double> ab = entry.getValue();
			double count = 0.0;

			for (Double item : ab) {
				count += item;
			}
			if (!aa.contains("~") && count > 0) {
				String[] temp = aa.split("-");
				System.out.println(temp[0] + "\t" + temp[1] + "\t" + count);
			}
			if (aa.contains("~") && count > 0) {
				String[] temp = aa.split("-");
				System.out.println(temp[0] + "\t" + temp[1] + "\t" + count);
			}
		}

		// Resize all columns to fit the content size
		for (int i = 0; i < set.size() + 1; i++) {
			sheet.autoSizeColumn(i);
		}

		// Write the output to a file
		FileOutputStream fileOut = new FileOutputStream(outputFileName + ".xlsx");
		workbook.write(fileOut);
		fileOut.close();

		// Closing the workbook
		workbook.close();
		System.out.println("Done");
	}
}