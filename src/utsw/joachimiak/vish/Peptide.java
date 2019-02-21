package utsw.joachimiak.vish;

import java.util.Arrays;

class Peptide {

	private String fileID;
	private String annotatedSeq;
	private String[] phosphorylations;

	public int[] getTauPhosLocalization() {
		return tauPhosLocalization;
	}

	public void setTauPhosLocalization(int[] tauPhosLocalization) {
		this.tauPhosLocalization = tauPhosLocalization;
	}

	private int[] tauPhosLocalization;

	double getAbundance() {
		return abundance;
	}

	void setAbundance(double abundance) {
		this.abundance = abundance;
	}

	private double abundance;
	Peptide(String fileID, String seq, String[] phosphorylations, double abundance) {
		setFileID(fileID);
		setAnnotatedSeq(seq);
		setPhosphorylations(phosphorylations);
		setAbundance(abundance);
		tauPhosLocalization= new int[phosphorylations.length];
		Arrays.fill(tauPhosLocalization,-1);
	}

	String getFileID() {
		return fileID;
	}

	void setFileID(String fileID) {
		this.fileID = fileID;
	}

	String getAnnotatedSeq() {
		return annotatedSeq;
	}

	void setAnnotatedSeq(String annotatedSeq) {
		this.annotatedSeq = annotatedSeq;
	}

	String[] getPhosphorylations() {
		return phosphorylations;
	}

	void setPhosphorylations(String[] phosphorylations) {
		this.phosphorylations = phosphorylations;
	}

	@Override
	public String toString() {
		return "FileID " + getFileID() + " Seq: " + getAnnotatedSeq();
	}

}
