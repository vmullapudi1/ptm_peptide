package utsw.joachimiak.vish;

class Abundance {
	double phosphorylated;
	double total;

	Abundance(final double phosphorylated, final double total) {
		this.phosphorylated = phosphorylated;
		this.total = total;
	}

	@Override
	public String toString() {
		return "phosphorylated: " + this.phosphorylated + " total: " + this.total;
	}
}