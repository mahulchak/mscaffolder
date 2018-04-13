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

   delta-filter -1 mygenome.delta > mygenome.1delta

   delta-filter -m mygenome.delta > mygenome.mdelta

  ```
4. Run mscaffolder (use the <b>unmasked</b> genome file and NOT the masked genome file as fasta file input).

  ```
   mscaffolder -d1 mygenome.1delta -md mygenome.mdelta -f your_genome.fasta -ul n > my_scaffold.fasta

  ```

In some cases, you may want to order the contigs in a specifc way. E.g. you may like how mummerplot orders your contigs to a reference genome. To scaffold based on a user provided contig order, mscaffolder needs to be run as -

  ```
   mscaffolder -d1 mygenome.1delta -md mygenome.mdelta -f your_genome.fasta -ul y -l query_list > my_scaffold.fasta

  ```

i.e. the switch 'ul' needs to be turned on and a list needs to be provided. To obtain the list from mummer, run mummerplot on your 1 delta file as follows:
  
  ```
   mummerplot -fat -filter -postscript mygenome.1delta

   sed 's/["|,|\|)|(]//g' out.gp |tail -n +3 |head -n -30|awk '{if(NF >1)print $1"\t"$2}' > query_list

  ```

  If you want to create the list yourself, just make sure that it has "y tics" as the header (y and tics are separated by a tab) like -

  ```
   y    tics
   ctg1
   ctg2
   ....
   ctgn
 
  ```
<b>What else do you need to know?</b>

  a) You can use a subset of the chromosomes from the reference genome for scaffolding purposes.
  
  b) Note that the file input for scaffolder is not the masked file.
  
  c) The unscaffolded contigs are named with a "U_" prefix. You will see that a file called "ctgmap.txt" is generated after scaffolder runs. That file has the information about the unscaffolded contigs and information about which contig maps to which reference chromosome.
  
  
<b>Citation</b>
https://www.nature.com/articles/s41588-017-0010-y



