# mScaffolder
A comparative genome scaffolding tool based on MUMmer 

mScaffolder scaffolds a genome using an existing high quality genome as the reference. It aligns the two genomes using nucmer utility from MUMmer and then orders and orients the contigs of the candidate genome guided by their alignments to the reference genome. Please send your questions and comments to mchakrab@uci.edu. 

<b>What is needed?</b>

<i>Repeatmasker</i> or a similar program for masking the repeats, A reference genome fasta file, a genome that will be scaffolded, a working <i>MUMmer</i> installation.

<b>Workflow:</b>

1. Compile <i>mscaffolder</i>. Please go inside the folder and type 

  ```
   make
  ```
2. Mask the repeats in your genome and the reference genome such that the sequence lengths do not change after masking. Lets say that after masking the genomes, you have a 'reference.masked.fasta" and "your_genome.masked.fasta". The unmasked files are "reference.fasta" and "your_genome.fasta".

3. Run nucmer and delta-filter for 1-to-1 alignment between the reference genome and your genome.

  ```
   nucmer -mum -prefix mygenome reference.repmasked.fasta your_genome.masked.fasta
   delta-filter -m mygenome.delta > mygenome.mdelta

  ```
4. Run mscaffolder (use the <b>unmasked</b> genome file and NOT the masked genome file as fasta file input).

  ```
   mscaffolder -md mygenome.mdelta -f your_genome.fasta > my_scaffold.fasta
  ```

<b>What else do you need to know?</b>

  a) You can use a subset of the chromosomes from the reference genome for scaffolding purposes.
  
  b) Note that the file input for scaffolder is not the masked file.
  
  c) The unscaffolded contigs are named with a 'U_' prefix. You will see that a file called "ctgmap.txt" is generated after scaffolder runs. That file has the information about the unscaffolded contigs and information about which contig maps to which reference chromosome.
  
  
<b>Citation</b>
Coming soon!



