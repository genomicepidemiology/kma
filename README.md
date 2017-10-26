/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 All rights reserved.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/


Introduction

KMA is a tool designed to map raw reads directly against databases, in an 
ultra-fast manner using seed and extend. KMA is particulary good at aligning 
short reads against highly redundant databases, where unique matches often does 
not exist for each read. Where non-unique matches are resolved using the 
"ConClave" sorting scheme.
For samples contaminated with host DNA, such as human for metagenomic samples,
KMA is able to remove this by modifying the database in a way that neither
increase the database size nor the computation time significantly.
KMA supports fastq and fasta format, with sequences of any length. For mapping
nanopore reads it is recommended to convert the reads to fasta before mapping, 
as the trimming is done from the ends.

If you use KMA for your published research, then please cite the KMA paper.


Compilation

Any c-compiler for C99 can be used, here GCC.
gcc -O3 -o kma KMA.c -lm
gcc -O3 -o kma_index KMA_index.c
gcc -O3 -o kma_shm KMA_SHM.c


Indexing
KMA supports two forms of indexing, standard and sparse. The standard makes a
full indexing of the input for alignment. The sparse format index only k-mers 
with a certain prefix, and the alignement is not performed whereas this relies
solely on k-mer counting.

# Options are:	Desc:						Default:		Requirements:
#
#	-i			Input/query file name		STDIN
#	-o			Output file					REQUIRED
#	-batch		Batch input file
#	-deCon		File name of contamination	None/False
#	-t_db		Add to existing DB			None/False
#	-k			Kmersize					16
#	-ML			Minimum length of templates	0
#	-CS			Chain size					8 MB
#	-MS			Max chain size				14 GB
#	-SM			Start at max chain size		False
#	-QS			Start query size			1M
#	-Sparse		Make Sparse DB
#				('-' for no prefix)			None/False
#	-ht			Homology template			1.0
#	-hq			Homology query				1.0
#	-and		Both homolgy thresholds
#				has to be reached			or
#	-h			Shows this help message
#

Resistance indexing example:
kma_index -i ResFinder.fsa -o ResFinder

Resistance indexing, with decontamination information on the human genome example:
kma_index -i ResFinder.fsa -o ResFinder -deCon human.fsa

cgMLST indexing example, "-CS" sets the reallocation size:
kma_index -i cgMLSTscheme.fsa -o cgMLSTscheme -CS 256

Sparse indexing of bacterial genomes, with prefix "ATG":
kma_index -i bac_o.fsa -o bac_o -Sparse ATG

Sparse indexing of bacterial genomes, with prefix "ATG" and homolgy reduction 
to 90% template and query coverage:
kma_index -i bac_o.fsa -o bac_o -Sparse ATG -ht 0.9 -hq 0.9 -and


Shared memory

KMA supports mapping against databases in shared memory, meaning that several
processes of KMA can run i parallel on the same computer without having to load
the database more than once.
WARNING:
To avoid segments of shared memory hanging around, do not kill the program
while running. Let it finish and then take it down when the database is not
needed anymore. The shared database will be available untill it is destroyed or
the computer is restarted.

# Options are:	Desc:						Default:	Requirements:
#
#	-t_db		Template DB					None		REQUIRED
#	-deCon		Setup contamination DB		False
#	-destroy	Destroy shared DB			False
#	-Sparse		Sparse DB					False
#	-IOD		Dump converted DB			False
#	-IOL		Load previous converted DB	False
#	-h			Shows this help message
#

Setup shared memory for standard indexing:
kma_shm -t_db ResFinder

Setup shared memory for decontaminated indexing:
kma_shm -t_db ResFinder -deCon

Setup shared memory for sparse indexing:
kma_shm -t_db ResFinder -Sparse

When the database should not be used anymore, when you can take it down again
with the same command used to put it up with the addition of the "-destroy" 
option.


Mapping with KMA.

# Options are:	Desc:							Default:	Requirements:
#
#	-i			Input/query file name			None		REQUIRED
#	-o			Output file						None		REQUIRED
#	-t_db		Template DB						None		REQUIRED
#	-k			Kmersize						16
#	-e			evalue							0.05
#	-delta		Allocation size for sequences	512
#	-mem_mode	Use kmers to choose best
#				template, and save memory		False
#	-ex_mode	Searh kmers exhaustively		False
#	-deCon		Remove contamination			False
#	-dense		Do not allow insertions
#				in assembly						False
#	-ref_fsa	Consensus sequnce will
#				have "n" instead of gaps		False
#	-matrix		Print assembly matrix			False
#	-mp			Minimum phred score				30
#	-5p			Cut a constant number of
#				nucleotides from the 5 prime.	0
#	-CS			Chain size						8 MB
#	-MS			Max chain size					14 GB
#	-Sparse		Run KmerFinder					False
#	-ID			Minimum ID						1.0%
#	-ss			Sparse sorting (q,c,d)			q
#	-NW			Use Needleman-Wunsch			False
#	-shm		Use shared DB made by kma_ssm	False
#	-1t1		Skip HMM						False
#	-mrs		Minimum alignment score score,
#				normalized to alignment length	0.0
#	-h			Shows this help message
#

Mapping short reads against resistance genes:
kma -i sample_1.fastq sample_2.fastq -o out_res -t_db ResFinder -NW -1t1

Mapping against resistance genes, with decontamination:
kma -i sample_1.fastq sample_2.fastq -o out_res -t_db ResFinder -NW -1t1 -deCon

Mapping against resistance genes, with shared database:
kma -i sample_1.fastq sample_2.fastq -o out_res -t_db ResFinder -NW -1t1 -shm

Mapping long reads against resisntance genes (potential for several matches per
read):
kma -i sample_1.fastq sample_2.fastq -o out_res -t_db ResFinder -NW

Mapping agains large database, where memory consumption is a concern.
kma -i sample_1.fastq sample_2.fastq -o out_cgMLST -t_db cgMLSTscheme -NW /
-mem_mode

Sparse mapping.
kma -i sample_1.fastq sample_2.fastq -o out_spec -t_db bac_o -Sparse

Sparse mapping, with shared database.
kma -i sample_1.fastq sample_2.fastq -o out_spec -t_db bac_o -Sparse -shm


Output:
*.res		Result file, containing summary of output.
*.aln		Consensus alignment.
*.fsa		Consensus sequences drawn from mappings.
*.frag		Information about each mapping read, containing:
			Read, #matches, aln score, start, end, template, read name
*.mat.gz	Count of each called nucleotide on each position in all mapped 
			templates, requires that the "-matrix" option is enabled when 
			mapping. The columns are:
			Ref. nucleotide, #A, #T, #C, #G, #N, #-.


Questions
If you have any questions regarding KMA, or wishes for changes in future
versions of KMA, then feel free to mail me at: plan@dtu.dk.

