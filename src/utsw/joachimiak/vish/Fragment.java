package utsw.joachimiak.vish;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class Fragment {

	boolean containsUnlocalizedPhosphorylation;
	private String fileID;
	private String annotatedSeq;
	private String sequence;
	private int indexInProtein = -1;
	private String[] phosphorylations;
	private int[] proteinPhosLocalizations;
	private double abundance;

	Fragment(final String fileID, final String seq, final String modifications, final double abundance) {
		this.setFileID(fileID);
		this.setAnnotatedSeq(seq);
		this.setPhosphorylations(this.initPhosphorylations(modifications));
		this.setAbundance(abundance);
		this.proteinPhosLocalizations = new int[this.phosphorylations.length];
		Arrays.fill(this.proteinPhosLocalizations, -1);
		this.setSequence(seq.substring(seq.indexOf('.') + 1, seq.lastIndexOf('.')));
	}

	Fragment(final String fileID, final String seq, final String[] phosphorylations, final int[] localizedPhosSites, final int protIndex, final boolean containsUnlocalizedPhosphorylation) {
		this.setFileID(fileID);
		this.setSequence(seq);
		this.setPhosphorylations(phosphorylations);
		this.proteinPhosLocalizations = localizedPhosSites;
		Arrays.fill(this.proteinPhosLocalizations, -1);
		this.setSequence(seq);
		this.indexInProtein = protIndex;
		this.containsUnlocalizedPhosphorylation = containsUnlocalizedPhosphorylation;
	}

	String getSequence() {
		return this.sequence;
	}

	private void setSequence(final String sequence) {
		this.sequence = sequence;
	}

	double getAbundance() {
		return this.abundance;
	}

	private void setAbundance(final double abundance) {
		this.abundance = abundance;
	}

	private String[] initPhosphorylations(final String modifications) {
		if ("".equals(modifications) || !modifications.contains("Phospho")) {
			return new String[0];
		}
		final ArrayList<String> phosphorylations = new ArrayList<>();

		//Assumes that the peptide will not have more than 0-2 digits specifying the localization
		//And that the only modifications to serine, threonine, and tyrosine will be phosphorylation
		final Pattern phosRegex = Pattern.compile("[STY]([\\d]{0,2})");
		final Matcher m = phosRegex.matcher(modifications);
		while (m.find()) {
			phosphorylations.add(m.group());
		}
		final String[] s = new String[phosphorylations.size()];
		return phosphorylations.toArray(s);
	}

	String getFileID() {
		return this.fileID;
	}

	private void setFileID(final String fileID) {
		this.fileID = fileID;
	}

	private String getAnnotatedSeq() {
		return this.annotatedSeq;
	}

	private void setAnnotatedSeq(final String annotatedSeq) {
		this.annotatedSeq = annotatedSeq;
	}

	String[] getPhosphorylations() {
		return this.phosphorylations;
	}

	private void setPhosphorylations(final String[] phosphorylations) {
		this.phosphorylations = phosphorylations.clone();
	}

	int[] getProteinPhosLocalizations() {
		return this.proteinPhosLocalizations;
	}

	void setProteinPhosLocalizations(final int[] tauPhosLocalization) {
		proteinPhosLocalizations = tauPhosLocalization.clone();
	}

	int getIndexInProtein() {
		return this.indexInProtein;
	}

	void setPeptideProteinIndex(final int index) {
		indexInProtein = index;
	}

	@Override
	public String toString() {
		return "FileID " + this.getFileID() + " Seq: " + this.getAnnotatedSeq();
	}

}