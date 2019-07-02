# ptm_peptide

This project is a system to extract and analyze phosphorylation data produced by possttranslational modification analysis performed by the CPFP proteomics pipeline.
The code ingests .csv files containing fragment sequences, postranslational modification data, and abundance data and delivers a csv. file detailing the modification proportion of each residue or peptide fragment output by the PTM analysis.

Future improvements:
Include a config file to set the protein sequence and other rutnime options
Generalize the project to perform analysis of other posttranslational modifications, such as acetylation, GlcONAc, geranylation, farnesylation, etc.
Allow for piping I/O for integration into other analysis pipelines
