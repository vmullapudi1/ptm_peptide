package utsw.joachimiak.vish;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class Peptide {

	private String fileID;
	private String annotatedSeq;
	private int phosphorylationCount;
	private String sequence;
	private int tauIndex;
	private String[] phosphorylations;
	private int[] tauPhosLocalization;
	private double abundance;
	private boolean modificationInPeptideFormat;
	Peptide(String fileID, String seq, String modifications, double abundance) {
		setFileID(fileID);
		setAnnotatedSeq(seq);
		modificationInPeptideFormat = isPeptideFormat(modifications);
		//phosphorylationCount=countPhosphorylations(modifications);
		setPhosphorylations(initPhosphorylations(modifications));
		setAbundance(abundance);
		tauPhosLocalization = new int[phosphorylations.length];
		Arrays.fill(tauPhosLocalization, -1);
		setSequence(seq.substring(seq.indexOf('.') + 1, seq.lastIndexOf('.')));
	}

	String getSequence() {
		return sequence;
	}

	private void setSequence(String sequence) {
		this.sequence = sequence;
	}

	double getAbundance() {
		return abundance;
	}

	private void setAbundance(double abundance) {
		this.abundance = abundance;
	}

	private String[] initPhosphorylations(String modifications) {
		if (modifications.equals("") || modifications == null) {
			return new String[0];
		}
		//Todo-use number from the phosphorylation count to go directly to array
		ArrayList<String> phosphorylations = new ArrayList<>();

		if (modificationInPeptideFormat) {
			//TODO bug in this regex-phosphorylation extracting algorithm. Non-regex version below works on original dataset,
			//	but can't parse out the new peptide file format. New version should hopefully be robust to both file formats
			modifications = modifications.split("Phospho \\[")[1];
			Pattern phosRegex = Pattern.compile("[STY][0-9]+?");
			Matcher m = phosRegex.matcher(modifications);
			while (m.find()) {
				phosphorylations.add(m.group());
			}

		} else {
			String[] split = modifications.split(";");
			for (String s : split) {
				if (s.contains("Phospho")) {
					phosphorylations.add(s.replaceAll("\\D", ""));
				}
			}
		}

		String[] s = new String[phosphorylations.size()];
		return phosphorylations.toArray(s);
	}

	//todo implement this method and make it robust to non-assigned phophorylations and interference from other
	//	modifications
	private int countPhosphorylations(String modifications) {
		if (modifications.matches("[STY]")) {
			modifications = modifications.substring(modifications.indexOf("["), modifications.indexOf("]"));
		}
		return modifications.split(";").length;
	}

	private boolean isPeptideFormat(String modifications) {
		return !(modifications.equals("") || modifications == null || modifications.matches("Phospho\\("));
	}

	String getFileID() {
		return fileID;
	}

	private void setFileID(String fileID) {
		this.fileID = fileID;
	}

	String getAnnotatedSeq() {
		return annotatedSeq;
	}

	private void setAnnotatedSeq(String annotatedSeq) {
		this.annotatedSeq = annotatedSeq;
	}

	String[] getPhosphorylations() {
		return phosphorylations;
	}

	private void setPhosphorylations(String[] phosphorylations) {
		this.phosphorylations = phosphorylations;
	}

	int[] getTauPhosLocalization() {
		return tauPhosLocalization;
	}

	void setTauPhosLocalization(int[] tauPhosLocalization) {
		this.tauPhosLocalization = tauPhosLocalization;
	}

	int getTauIndex() {
		return tauIndex;
	}

	void setTauIndex(int tauIndex) {
		this.tauIndex = tauIndex;
	}

	@Override
	public String toString() {
		return "FileID " + getFileID() + " Seq: " + getAnnotatedSeq();
	}

}
