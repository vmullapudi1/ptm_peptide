package utsw.joachimiak.vish;

import java.util.ArrayList;
import java.util.Arrays;

class Peptide {

	private String fileID;
	private String annotatedSeq;

	private String sequence;
	private int tauIndex;

	Peptide(String fileID, String seq, String modifications, double abundance) {
		setFileID(fileID);
		setAnnotatedSeq(seq);
		setPhosphorylations(initphosphorylations(modifications));
		setAbundance(abundance);
		tauPhosLocalization= new int[phosphorylations.length];
		Arrays.fill(tauPhosLocalization,-1);
		sequence = seq.substring(seq.indexOf('.') + 1, seq.lastIndexOf('.'));
	}

	private String[] phosphorylations;

	public String getSequence() {
		return sequence;
	}

	private int[] tauPhosLocalization;

	double getAbundance() {
		return abundance;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	private double abundance;

	private void setAbundance(double abundance) {
		this.abundance = abundance;
	}

	private String[] initphosphorylations(String modifications) {
		ArrayList<String> phosphorylations = new ArrayList<>();
		for (String s : modifications.split(";")) {
			if (s.contains("Phospho")) {
				phosphorylations.add(s.replaceAll("\\D", ""));
			}
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

	public int[] getTauPhosLocalization() {
		return tauPhosLocalization;
	}

	public void setTauPhosLocalization(int[] tauPhosLocalization) {
		this.tauPhosLocalization = tauPhosLocalization;
	}

	public int getTauIndex() {
		return tauIndex;
	}

	public void setTauIndex(int tauIndex) {
		this.tauIndex = tauIndex;
	}

	@Override
	public String toString() {
		return "FileID " + getFileID() + " Seq: " + getAnnotatedSeq();
	}

}
