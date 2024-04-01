# Cappable-seq RNA-seq
This repository contains files and analysis for RNA-seq derived from cappable seq.
## Genome
Reads will align to both human and WSN genomes as we expect viral RNAs to be captured in cappable-seq. Therefore, I will create a custom combo of hg38 annotated from RefSeq and our lab's WSN genome with custom annotations.
### Procedure
This will be detailed in a jupyter notebook in the `/genome` folder. `STAR` and `bowtie2` are both able to create indices from this custom reference genome. Annotations were based on sequences obtained from GenBank. Annotation is done using the `gb2gtf.py` script and this procedure is detailed in the notebook.