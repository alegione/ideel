# ideel
Indels are not ideal - quick test for interrupted ORFs in bacterial/microbial genomes

# Dependencies:
* Python
  * Snakemake
  * altair
  * pandas
  * selenium
* Prodigal
* Diamond

You will need a diamond index of UniProt TREMBL, called uniprot_trembl.diamond.dmnd and place it in the 'ideel' directory

# run

Clone the repo.  Make a directory called "genomes", put assemblies in there with .fa file extension, then:

```
snakemake -s ideel.smk
```
