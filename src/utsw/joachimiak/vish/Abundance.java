package utsw.joachimiak.vish;

class Abundance {
	double phosphorylated;
	double total;

	Abundance(double phosphorylated, double total) {
		this.phosphorylated = phosphorylated;
		this.total = total;
	}

	@Override
	public String toString() {
		return "phosphorylated: " + phosphorylated + " total: " + total;
	}
}
