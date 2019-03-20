package utsw.joachimiak.vish;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class Peptide {

	private String fileID;
	private String annotatedSeq;
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
		if (modifications.equals("") || modifications == null || !modifications.contains("Phospho")) {
			return new String[0];
		}
		//Todo-use number from the phosphorylation count to go directly to array
		ArrayList<String> phosphorylations = new ArrayList<>();

		//if (modificationInPeptideFormat) {
			//Assumes that the peptide will not have more than 0-2 digits specifying the localization
			//And that the only modifications to serine, threonine, and tyrosine will be phosphorylation
			//TODO bug in this regex-phosphorylation extracting algorithm. Non-regex version below works on original dataset,
			//	but can't parse out the new peptide file format. New version should hopefully be robust to both file formats

			Pattern phosRegex = Pattern.compile("[STY]([\\d]{0,2})");
			Matcher m = phosRegex.matcher(modifications);
			while (m.find()) {
				phosphorylations.add(m.group());
			}

		/*} else {

			String[] split = modifications.split(";");
			for (String s : split) {
				if (s.contains("Phospho")) {
					phosphorylations.add(s.replaceAll("\\D", ""));
				}
			}
		}*/

		String[] s = new String[phosphorylations.size()];
		return phosphorylations.toArray(s);
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

	void setProteinPhosLocalizations(int[] tauPhosLocalization) {
		this.tauPhosLocalization = tauPhosLocalization;
	}

	int getTauIndex() {
		return tauIndex;
	}

	void setPeptideProteinIndex(int tauIndex) {
		this.tauIndex = tauIndex;
	}

	@Override
	public String toString() {
		return "FileID " + getFileID() + " Seq: " + getAnnotatedSeq();
	}

}
