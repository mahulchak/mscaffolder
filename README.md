# scaffolder
A comparative genome scaffolding tool 

Scaffolder scaffolds a genome using an existing high quality genome as the reference. It aligns the two genomes using nucmer from MUMmer and then orders and orients the contigs in the candidate genome guided by their alignments to the reference genome. 

<b>What is needed?</b>

Repeatmasker for masking the repeats, A reference genome fasta file, a genome that will be scaffolded, a working MUMmer installation.

<b>Workflow:</b>

1. Mask the repeats in your genome (required) and the reference genome (optional) such that the sequence lengths do not change after masking. Lets say that after masking the genomes, you have a 'reference.masked.fasta" and "your_genome.masked.fasta". The unmasked files are "reference.fasta" and "your_genome.fasta".

2. Run nucmer for 1-to-1 alignment between the reference genome and your genome.

  ```
   nucmer -mum -prefix mygenome reference.repmasked.fasta your_genome.masked.fasta
  ```
3. Run scaffolder (use the unmasked genome file as file input).

  ```
   scaffolder -D mygenome.delta -f your_genome.fasta > my_scaffold.fasta
  ```

<b>What else do you need to know?</b>

  a) You can use a subset of the chromosomes from the reference genome for scaffolding purposes.
  
  b) Note that the file input for scaffolder is not the masked file.
  
  c) You will notice a file called "ctgmap.txt" after scaffolder runs. That file has the information about the unscaffolded contigs.
  
  d) I will soon be adding another utility for manually adding an unscaffolded contig into an existing scaffold.
