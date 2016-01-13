# pacbio-14-nctc-assemblies

14 [NCTC][nctc home] samples for testing PacBio assemblers.
They are the 14 samples used in the [Circlator paper][circlator paper].

Assemblies and files for analysis are all in this github repository.
Raw sequencing reads are in the [ENA][ena homepage]. The filtered
subreads FASTQ and corrected reads FASTA files made when running HGAP are available
from
[ftp://ngs.sanger.ac.uk/production/pathogens/mh12/pacbio-14-nctc-assemblies](ftp://ngs.sanger.ac.uk/production/pathogens/mh12/pacbio-14-nctc-assemblies).

The file `sample_data.tsv` lists accession IDs for the raw reads and the
reference assembly of each sample, and some basic stats (assembly
size, number of reads etc).

Each directory `NCTCxxxxx/` contains all the files
relating to that sample.
The FASTA files in each directory are:

* `ref.fa` - the reference sequence
* `canu.1.0.fa` - as assembly made with version 1.0 of [canu][canu github]
* `miniasm.fa` - an assembly made with [miniasm][miniasm github]
  (preprint [here][miniasm arxiv]).
* `hgap.fa` - an assembly made with [HGAP][hgap github]
  (publication [here][hgap paper]).


## Canu assemblies

Made with canu version 1.0. The filtered subreads were used as input with

    -pacbio-raw filtered_subreads.fq

and the genome size was set to the length of the reference genome for
each sample, using

    genomeSize=$length

where `$length` was taken from the file `sample_data.tsv`.

The only other options changed were cluster-specific:

    maxThreads=8 maxMemory=16 useGrid=0


## HGAP assemblies

Details TBC...


## miniasm assemblies

Made with version 0.2 (and [minimap][minimap github] version 0.2)
using the filtered subreads output during a run
of HGAP. The three commands run were:

    minimap -Sw5 -L100 -m0 -t4 $reads $reads | gzip -1 > miniasm.paf.gz
    miniasm -f $reads miniasm.paf.gz > miniasm.gfa
    awk '$1=="S" {print ">"$2"\n"$3} ' miniasm.gfa > miniasm.fa

where `$reads` is the FASTQ file of reads, and the final output
FASTA file of contigs is called `miniasm.fa`.


## To do

* Gather HGAP assembler version/options etc
* Rename miniasm files with version number in them
* Add PBcR assemblies
* Add SPRAI assemblies
* Run Quast on all assemblies/refs
* Try quiver on miniasm assemblies


[canu github]: https://github.com/marbl/canu
[ena homepage]: http://www.ebi.ac.uk/ena
[ftp reads]: ftp://ngs.sanger.ac.uk/production/pathogens/mh12/pacbio-14-nctc-assemblies/
[hgap github]: https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP
[hgap paper]: http://www.nature.com/nmeth/journal/v10/n6/abs/nmeth.2474.html
[miniasm arxiv]: http://arxiv.org/abs/1512.01801
[miniasm github]: https://github.com/lh3/miniasm
[minimap github]: https://github.com/lh3/minimap
[circlator paper]: http://www.genomebiology.com/2015/16/1/294
[nctc home]: https://www.phe-culturecollections.org.uk/collections/nctc.aspx
