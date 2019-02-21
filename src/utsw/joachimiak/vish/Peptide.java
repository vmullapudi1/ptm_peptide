package utsw.joachimiak.vish;

class Peptide {

	private String fileID;
	private String annotatedSeq;
	private String[] phosphorylations;

	Peptide(String fileID, String seq, String[] phosphorylations) {
		setFileID(fileID);
		setAnnotatedSeq(seq);
		setPhosphorylations(phosphorylations);
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
