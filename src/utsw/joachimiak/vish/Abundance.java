package utsw.joachimiak.vish;

class Abundance {
	double phosphorylated;
	double total;

	Abundance(double phos, double tot) {
		this.phosphorylated = phos;
		this.total = tot;
	}

	@Override
	public String toString() {
		return "phosphorylated: " + this.phosphorylated + " total: " + this.total;
	}
}