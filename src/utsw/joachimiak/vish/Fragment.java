package utsw.joachimiak.vish;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class Fragment {

	private String fileID;
	private String annotatedSeq;
	private String sequence;
	private int indexInProtein = -1;
	private String[] phosphorylations;
	private int[] proteinPhosLocalizations;
	private double abundance;
	boolean containsUnlocalizedPhosphorylation;

	Fragment(String fileID, String seq, String modifications, double abundance) {
		setFileID(fileID);
		setAnnotatedSeq(seq);
		setPhosphorylations(initPhosphorylations(modifications));
		setAbundance(abundance);
		proteinPhosLocalizations = new int[phosphorylations.length];
		Arrays.fill(proteinPhosLocalizations, -1);
		setSequence(seq.substring(seq.indexOf('.') + 1, seq.lastIndexOf('.')));
	}

	Fragment(String fileID, String seq, String[] phosphorylations, int[] localizedPhosSites, int protIndex) {
		setFileID(fileID);
		setSequence(seq);
		setPhosphorylations(phosphorylations);
		proteinPhosLocalizations = new int[phosphorylations.length];
		Arrays.fill(proteinPhosLocalizations, -1);
		setSequence(seq);
		indexInProtein = protIndex;
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
		ArrayList<String> phosphorylations = new ArrayList<>();

		//Assumes that the peptide will not have more than 0-2 digits specifying the localization
		//And that the only modifications to serine, threonine, and tyrosine will be phosphorylation
		Pattern phosRegex = Pattern.compile("[STY]([\\d]{0,2})");
		Matcher m = phosRegex.matcher(modifications);
		while (m.find()) {
			phosphorylations.add(m.group());
		}
		String[] s = new String[phosphorylations.size()];
		return phosphorylations.toArray(s);
	}

	String getFileID() {
		return fileID;
	}

	private void setFileID(String fileID) {
		this.fileID = fileID;
	}

	private String getAnnotatedSeq() {
		return annotatedSeq;
	}

	private void setAnnotatedSeq(String annotatedSeq) {
		this.annotatedSeq = annotatedSeq;
	}

	String[] getPhosphorylations() {
		return phosphorylations;
	}

	private void setPhosphorylations(String[] phosphorylations) {
		this.phosphorylations = phosphorylations.clone();
	}

	int[] getProteinPhosLocalizations() {
		return proteinPhosLocalizations;
	}

	void setProteinPhosLocalizations(int[] tauPhosLocalization) {
		this.proteinPhosLocalizations = tauPhosLocalization.clone();
	}

	int getIndexInProtein() {
		return indexInProtein;
	}

	void setPeptideProteinIndex(int index) {
		this.indexInProtein = index;
	}

	@Override
	public String toString() {
		return "FileID " + getFileID() + " Seq: " + getAnnotatedSeq();
	}

}
