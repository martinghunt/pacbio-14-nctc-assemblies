# pacbio-14-nctc-assemblies

14 [NCTC][nctc home] samples for testing PacBio assemblers.
They are the 14 samples used in the [Circlator paper][circlator paper].

Assemblies and files for analysis are all in this github repository.
Raw sequences reads are in the ENA.

The file `sample_data.tsv` lists accession IDs for the raw reads and the
reference assembly of each sample, and some basic stats (assembly
size, number of reads etc).

Each directory `NCTCxxxxx/` contains all the files
relating to that sample.
The FASTA files in each directory are:

* `ref.fa` - the reference sequence
* `miniasm.fa` - an assembly made with [miniasm][miniasm github]
  (preprint [here][miniasm arxiv]).
* `hgap.fa` - an assembly made with [HGAP][hgap github]
  (publication [here][hgap paper]).


## HGAP assemblies

Details TBC...

## miniasm assemblies

Made with version 0.2 using the filtered subreads output during a run
of HGAP. The three commands run were:

    minimap -Sw5 -L100 -m0 -t4 $reads $reads | gzip -1 > miniasm.paf.gz
    miniasm -f $reads miniasm.paf.gz > miniasm.gfa
    awk '$1=="S" {print ">"$2"\n"$3} ' miniasm.gfa > miniasm.fa

where `$reads` is the FASTQ file of reads, and the final output
FASTA file of contigs is called `miniasm.fa`.


[hgap github]: https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP
[hgap paper]: http://www.nature.com/nmeth/journal/v10/n6/abs/nmeth.2474.html
[miniasm arxiv]: http://arxiv.org/abs/1512.01801
[miniasm github]: http://www.genomebiology.com/2015/16/1/294
[circlator paper]: http://www.genomebiology.com/2015/16/1/294
[nctc home]: https://www.phe-culturecollections.org.uk/collections/nctc.aspx
