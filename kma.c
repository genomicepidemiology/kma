/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ankers.h"
#include "align.h"
#include "alnfrags.h"
#include "assembly.h"
#include "chain.h"
#include "hashmapkma.h"
#include "kma.h"
#include "kmapipe.h"
#include "kmeranker.h"
#include "kmers.h"
#include "kmmap.h"
#include "mt1.h"
#include "penalties.h"
#include "pherror.h"
#include "qseqs.h"
#include "runinput.h"
#include "runkma.h"
#include "savekmers.h"
#include "sparse.h"
#include "spltdb.h"
#include "stdstat.h"
#include "tmp.h"
#include "vcf.h"
#include "version.h"

char * strjoin(char **strings, int len) {
	
	int i, new_len, escape;
	char *newStr, *stringPtr;
	
	new_len = len + 16;
	escape = 0;
	for(i = 0; i < len; ++i) {
		if(*strings[i] == '-') {
			escape = 0;
		} else if(escape) {
			new_len += 2;
		}
		new_len += strlen(strings[i]);
		if(strncmp(strings[i], "-i", 2) == 0) {
			escape = 1;
		}
	}
	
	newStr = smalloc(new_len);
	
	*newStr = 0;
	escape = 0;
	stringPtr = newStr;
	for(i = 0; i < len; ++i) {
		if(*strings[i] == '-') {
			escape = 0;
		}
		
		if(escape) {
			*stringPtr = '\"';
			++stringPtr;
		}
		new_len = strlen(strings[i]);
		strcpy(stringPtr, strings[i]);
		stringPtr += new_len;
		if(escape) {
			*stringPtr = '\"';
			++stringPtr;
		}
		*stringPtr = ' ';
		++stringPtr;
		
		if(*strings[i] == '-' && (strings[i][1] == 'i' || strings[i][1] == 'o')) {
			escape = 1;
		}
	}
	
	return newStr;
}

static void helpMessage(int exitStatus) {
	
	FILE *out = exitStatus ? stderr : stdout;
	
	fprintf(out, "# KMA-%s maps and/or aligns raw reads to a template database.\n", KMA_VERSION);
	
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options:", "Desc:", "Default:");
	
	fprintf(out, "#\n# Input:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Single end input(s)", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ipe", "Paired end input(s)", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-int", "Interleaved input(s)", "");
	
	fprintf(out, "#\n# Output:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output prefix", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ef", "Output additional features", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-vcf", "Output vcf file, 2 to apply FT", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-sam", "Output sam, 4/2096 for mapped/aligned", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nc", "No consensus file", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nc", "No aln file", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nf", "No frag file", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-matrix", "Output assembly matrix", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-a", "Output all template mappings", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-and", "Use both mrs and p-value on consensus", "or");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-oa", "Use neither mrs or p-value on consensus", "False");
	
	fprintf(out, "#\n# Consensus:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bc", "Minimum support to call bases", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bcNano", "Altered indel calling for ONT data", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bcd", "Minimum depth to cal bases", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bcg", "Maintain insignificant gaps", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ID", "Minimum consensus ID", "1.0%");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dense", "Skip insertion in consensus", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ref_fsa", "Use n's on indels", "False");
	
	fprintf(out, "#\n# General:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t_db", "Template DB", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-p", "P-value", "0.05");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-shm", "Use DB in shared memory", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mmap", "Memory map *.comp.b", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t", "Number of threads", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-status", "Extra status", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-verbose", "Extra verbose", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-c", "Citation", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-v", "Version", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this help message", "");
	
	fprintf(out, "#\n# Template mapping:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ConClave", "ConClave version", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mem_mode", "Base ConClave on template mappings", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-proxi", "Proximity scoring (negative for soft)", "False/1.0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ex_mode", "Searh kmers exhaustively", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-deCon", "Remove contamination", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-Sparse", "Only count kmers", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ss", "Sparse sorting (q,c,d)", "q");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-Mt1", "Map everything to one template", "False/0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-pm", "Pairing method (p,u,f)", "u");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-1t1", "One query to one template", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-hmm", "Use a HMM to assign template(s)", "True");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ck", "Count k-mers over pseudo alignment", "False");
	
	fprintf(out, "#\n# Chaining:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-k", "K-mersize", "DB defined");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ts", "Trim front of seeds", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ssa", "Seeds soround alignments", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ex_mode", "Searh kmers exhaustively", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fpm", "Pairing method (p,u,f)", "u");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mq", "Minimum mapping quality", "0");
	
	fprintf(out, "#\n# Alignment:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ca", "Circular alignments", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mrs", "Minimum relative alignment score", "0.5");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mrc", "Minimum query coverage", "0.0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum alignment length", "16");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-reward", "Score for match", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-penalty", "Penalty for mismatch", "2");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-gapopen", "Penalty for gap opening", "3");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-gapextend", "Penalty for gap extension", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-per", "Reward for pairing reads", "7");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-Npenalty", "Penalty matching N", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-transition", "Penalty for transition", "2");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-transversion", "Penalty for transversion", "2");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-sasm", "Skip alignment", "False");
	
	fprintf(out, "#\n# Trimming:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mp", "Minimum phred score", "20");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-eq", "Minimum avg. quality score", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-5p", "Trim 5 prime", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-3p", "Trim 3 prime", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum length", "16");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-xl", "Maximum length on se", "2147483647");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-boot", "Bootstrap sub-sequence", "False");
	
	fprintf(out, "#\n# Presets:\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-apm", "Sets both pm and fpm", "u");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-cge", "Set CGE penalties and rewards", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mint2", "Set 2ng gen Mintyper preset", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mint3", "Set 3rd gen Mintyper preset", "False");
	
	fprintf(out, "#\n");
	
	exit(exitStatus);
}

int kma_main(int argc, char *argv[]) {
	
	static const double prob[128] = {
		1.0000000000000000, 0.7943282347242815, 0.6309573444801932, 0.5011872336272722, 0.3981071705534972, 0.3162277660168379, 0.2511886431509580, 0.1995262314968880,
		0.1584893192461113, 0.1258925411794167, 0.1000000000000000, 0.0794328234724281, 0.0630957344480193, 0.0501187233627272, 0.0398107170553497, 0.0316227766016838,
		0.0251188643150958, 0.0199526231496888, 0.0158489319246111, 0.0125892541179417, 0.0100000000000000, 0.0079432823472428, 0.0063095734448019, 0.0050118723362727,
		0.0039810717055350, 0.0031622776601684, 0.0025118864315096, 0.0019952623149689, 0.0015848931924611, 0.0012589254117942, 0.0010000000000000, 0.0007943282347243,
		0.0006309573444802, 0.0005011872336273, 0.0003981071705535, 0.0003162277660168, 0.0002511886431510, 0.0001995262314969, 0.0001584893192461, 0.0001258925411794,
		0.0001000000000000, 0.0000794328234724, 0.0000630957344480, 0.0000501187233627, 0.0000398107170553, 0.0000316227766017, 0.0000251188643151, 0.0000199526231497,
		0.0000158489319246, 0.0000125892541179, 0.0000100000000000, 0.0000079432823472, 0.0000063095734448, 0.0000050118723363, 0.0000039810717055, 0.0000031622776602,
		0.0000025118864315, 0.0000019952623150, 0.0000015848931925, 0.0000012589254118, 0.0000010000000000, 0.0000007943282347, 0.0000006309573445, 0.0000005011872336,
		0.0000003981071706, 0.0000003162277660, 0.0000002511886432, 0.0000001995262315, 0.0000001584893192, 0.0000001258925412, 0.0000001000000000, 0.0000000794328235,
		0.0000000630957344, 0.0000000501187234, 0.0000000398107171, 0.0000000316227766, 0.0000000251188643, 0.0000000199526231, 0.0000000158489319, 0.0000000125892541,
		0.0000000100000000, 0.0000000079432823, 0.0000000063095734, 0.0000000050118723, 0.0000000039810717, 0.0000000031622777, 0.0000000025118864, 0.0000000019952623,
		0.0000000015848932, 0.0000000012589254, 0.0000000010000000, 0.0000000007943282, 0.0000000006309573, 0.0000000005011872, 0.0000000003981072, 0.0000000003162278,
		0.0000000002511886, 0.0000000001995262, 0.0000000001584893, 0.0000000001258925, 0.0000000001000000, 0.0000000000794328, 0.0000000000630957, 0.0000000000501187,
		0.0000000000398107, 0.0000000000316228, 0.0000000000251189, 0.0000000000199526, 0.0000000000158489, 0.0000000000125893, 0.0000000000100000, 0.0000000000079433,
		0.0000000000063096, 0.0000000000050119, 0.0000000000039811, 0.0000000000031623, 0.0000000000025119, 0.0000000000019953, 0.0000000000015849, 0.0000000000012589,
		0.0000000000010000, 0.0000000000007943, 0.0000000000006310, 0.0000000000005012, 0.0000000000003981, 0.0000000000003162, 0.0000000000002512, 0.0000000000001995};
	static int minPhred, minQ, fiveClip, threeClip, ConClave, mem_mode;
	static int fileCounter, fileCounter_PE, fileCounter_INT, Ts, Tv, minlen;
	static int extendedFeatures, spltDB, thread_num, kmersize, targetNum, mq;
	static int ref_fsa, print_matrix, print_all, sam, vcf, Mt1, bcd, one2one;
	static int sparse_run, ts, maxlen, **d, status = 0;
	static unsigned xml, nc, nf, shm, exhaustive, verbose;
	static char *outputfilename, *templatefilename, **templatefilenames;
	static char **inputfiles, **inputfiles_PE, **inputfiles_INT, ss;
	static double ID_t, scoreT, coverT, mrc, evalue, minFrac, support;
	static Penalties *rewards;
	int i, j, args, exe_len, fileCount, size, escape, tmp, step1, step2;
	unsigned totFrags;
	char *to2Bit, *exeBasic, *myTemplatefilename;
	FILE *templatefile, *ioStream;
	time_t t0, t1;
	Qseqs qseq;
	HashMapKMA *templates;
	
	step1 = 0;
	step2 = 0;
	
	if(argc) {
		if(sizeof(long unsigned) != 8) {
			fprintf(stderr, "Need a 64-bit system.\n");
			exit(1);
		}
		
		/* SET DEFAULTS */
		ConClave = 1;
		verbose = 0;
		vcf = 0;
		xml = 0;
		sam = 0;
		nc = 0;
		nf = 0;
		targetNum = 0;
		spltDB = 0;
		extendedFeatures = 0;
		minPhred = 20;
		minQ = 0;
		fiveClip = 0;
		threeClip = 0;
		sparse_run = 0;
		fileCounter = 0;
		fileCounter_PE = 0;
		fileCounter_INT = 0;
		outputfilename = 0;
		templatefilename = 0;
		print_matrix = 0;
		print_all = 0;
		ref_fsa = 0;
		kmersize = 0;
		ts = 0;
		minlen = 16;
		maxlen = 2147483647;
		evalue = 0.05;
		support = 0.0;
		minFrac = 1.0;
		exhaustive = 0;
		shm = 0;
		mq = 0;
		bcd = 1;
		scoreT = 0.5;
		coverT = 0.1;
		mrc = 0.1;
		ID_t = 1.0;
		one2one = 0;
		ss = 'q';
		mem_mode = 0;
		rewards = smalloc(sizeof(Penalties));
		rewards->M = 1;
		rewards->MM = -2;
		rewards->U = -1;
		rewards->W1 = -3;
		rewards->Wl = -6;
		rewards->Mn = 0;
		rewards->PE = 7;
		Tv = -2;
		Ts = -2;
		thread_num = 1;
		inputfiles_PE = 0;
		inputfiles_INT = 0;
		inputfiles = 0;
		templatefilenames = 0;
		tmp = 0;
		Mt1 = 0;
		inputfiles = 0;
		deConPrintPtr = printPtr;
		
		/* PARSE COMMAND LINE OPTIONS */
		args = 1;
		while(args < argc) {
			if(strcmp(argv[args], "-t_db") == 0) {
				if(++args < argc) {
					templatefilename = malloc(strlen(argv[args]) + 64);
					if(!templatefilename) {
						ERROR();
					}
					strcpy(templatefilename, argv[args]);
					++targetNum;
					templatefilenames = realloc(templatefilenames, targetNum * sizeof(char *));
					templatefilenames[targetNum - 1] = templatefilename;
				}
				while(++args < argc && *argv[args] != '-') {
					templatefilename = malloc(strlen(argv[args]) + 64);
					if(!templatefilename) {
						ERROR();
					}
					strcpy(templatefilename, argv[args]);
					++targetNum;
					templatefilenames = realloc(templatefilenames, targetNum * sizeof(char *));
					templatefilenames[targetNum - 1] = templatefilename;
				}
				--args;
			} else if(strcmp(argv[args], "-i") == 0) {
				++args;
				fileCount = fileCounter;
				for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
					++fileCounter;
				}
				if(fileCounter == 0) {
					fprintf(stderr, "No files were specified.\n");
					exit(1);
				} else {
					inputfiles = realloc(inputfiles, fileCounter * sizeof(char *));
					if(!inputfiles) {
						ERROR();
					}
				}
				
				for(i = fileCount; i < fileCounter; ++i, ++args) {
					inputfiles[i] = argv[args];
				}
				--args;
			} else if(strcmp(argv[args], "-ipe") == 0) {
				++args;
				fileCount = fileCounter_PE;
				for(i = args; i < argc && strncmp(argv[i], "-", 1) != 0; ++i) {
					++fileCounter_PE;
				}
				if(fileCounter_PE % 2) {
					fprintf(stderr, "Uneven number of paired end files.\n");
					exit(1);
				} else if(fileCounter_PE == 0) {
					fprintf(stderr, "No paired end files were specified.\n");
					exit(1);
				} else {
					inputfiles_PE = realloc(inputfiles_PE, fileCounter_PE * sizeof(char *));
					if(!inputfiles_PE) {
						ERROR();
					}
				}
				
				for(i = fileCount; i < fileCounter_PE; ++i, ++args) {
					inputfiles_PE[i] = argv[args];
				}
				--args;
			} else if(strcmp(argv[args], "-int") == 0) {
				++args;
				fileCount = fileCounter_INT;
				for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
					++fileCounter_INT;
				}
				if(fileCounter_INT == 0) {
					fprintf(stderr, "No interleaved files were specified.\n");
					exit(1);
				}
				inputfiles_INT = realloc(inputfiles_INT, fileCounter_INT * sizeof(char *));
				if(!inputfiles_INT) {
					ERROR();
				}
				for(i = fileCount; i < fileCounter_INT; ++i, ++args) {
					inputfiles_INT[i] = argv[args];
				}
				--args;
			} else if(strcmp(argv[args], "-pm") == 0) {
				++args;
				if(args < argc) {
					if(*(argv[args]) == 'p') {
						save_kmers_pair = &save_kmers_penaltyPair;
					} else if(*(argv[args]) == 'u') {
						save_kmers_pair = &save_kmers_unionPair;
					} else if(*(argv[args]) == 'f') {
						save_kmers_pair = &save_kmers_forcePair;
					} else {
						fprintf(stderr, "Invalid argument at pairing method: \"-pm\"\n");
						fprintf(stderr, "Options are:\n");
						fprintf(stderr, "p:\tReward for pairing.\n");
						fprintf(stderr, "u:\tUnion of best hits.\n");
						fprintf(stderr, "f:\tForce paring.\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-fpm") == 0) {
				++args;
				if(args < argc) {
					if(*(argv[args]) == 'p') {
						alnFragsPE = &alnFragsPenaltyPE;
					} else if(*(argv[args]) == 'u') {
						alnFragsPE = &alnFragsUnionPE;
					} else if(*(argv[args]) == 'f') {
						alnFragsPE = &alnFragsForcePE;
					} else {
						fprintf(stderr, "Invalid argument at fine pairing method: \"-fpm\"\n");
						fprintf(stderr, "Options are:\n");
						fprintf(stderr, "p:\tReward for pairing.\n");
						fprintf(stderr, "u:\tUnion of best hits.\n");
						fprintf(stderr, "f:\tForce paring.\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-apm") == 0) {
				++args;
				if(args < argc) {
					if(*(argv[args]) == 'p') {
						alnFragsPE = &alnFragsPenaltyPE;
						save_kmers_pair = &save_kmers_penaltyPair;
					} else if(*(argv[args]) == 'u') {
						alnFragsPE = &alnFragsUnionPE;
						save_kmers_pair = &save_kmers_unionPair;
					} else if(*(argv[args]) == 'f') {
						alnFragsPE = &alnFragsForcePE;
						save_kmers_pair = &save_kmers_forcePair;
					} else {
						fprintf(stderr, "Invalid argument at fine pairing method: \"-fpm\"\n");
						fprintf(stderr, "Options are:\n");
						fprintf(stderr, "p:\tReward for pairing.\n");
						fprintf(stderr, "u:\tUnion of best hits.\n");
						fprintf(stderr, "f:\tForce paring.\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-ConClave") == 0) {
				++args;
				if(args < argc) {
					ConClave = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0 || ConClave < 0 || 2 < ConClave) {
						fprintf(stderr, " Invalid ConClave version specified.\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-o") == 0) {
				++args;
				if(args < argc) {
					outputfilename = malloc(strlen(argv[args]) + 64);
					if(!outputfilename) {
						ERROR();
					}
					strcpy(outputfilename, argv[args]);
				}
			} else if(strcmp(argv[args], "-deCon") == 0) {
				deConPrintPtr = &deConPrint;
				printPairPtr = &deConPrintPair;
			} else if(strcmp(argv[args], "-shm") == 0) {
				++args;
				if(args < argc && argv[args][0] != '-') {
					shm = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid shm-lvl specified.\n");
						exit(1);
					}
				} else {
					--args;
					shm = 3;
				}
			} else if(strcmp(argv[args], "-mmap") == 0 || strcmp(argv[args], "-swap") == 0) {
				shm |= 32;
				hashMapKMA_destroy = &hashMapKMA_munmap;
			} else if(strcmp(argv[args], "-t") == 0) {
				++args;
				if(args < argc && argv[args][0] != '-') {
					thread_num = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid number of threads specified.\n");
						exit(1);
					}
				} else {
					--args;
				}
				if(thread_num < 1) {
					thread_num = 1;
				}
			} else if(strcmp(argv[args], "-s1") == 0) {
				step1 = 1;
			} else if(strcmp(argv[args], "-s2") == 0) {
				step2 = 1;
			} else if(strcmp(argv[args], "-mem_mode") == 0) {
				mem_mode = 1;
				alignLoadPtr = &alignLoad_fly_mem;
				ankerPtr = &ankerAndClean_MEM;
			} else if(strcmp(argv[args], "-ex_mode") == 0) {
				exhaustive = 1;
			} else if(strcmp(argv[args], "-k") == 0) {
				++args;
				if(args < argc) {
					kmersize = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "# Invalid kmersize parsed\n");
						exit(1);
					} else if(kmersize == 0) {
						fprintf(stderr, "# Invalid kmersize parsed, using default\n");
						kmersize = 16;
					} else if(kmersize > 31) {
						fprintf(stderr, "# Invalid kmersize parsed, max size is 31\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-ts") == 0) {
				++args;
				if(args < argc) {
					ts = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0 || ts < 0 || ts > 30) {
						fprintf(stderr, "# Invalid seed trim parsed\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-ssa") == 0) {
				trimSeedsPtr = &trimSeedsNoLead;
				leadTailAlnPtr = &skipLeadAln;
				trailTailAlnPtr = &skipTrailAln;
			} else if(strcmp(argv[args], "-ml") == 0) {
				++args;
				if(args < argc) {
					minlen = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "# Invalid kmersize parsed\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-xl") == 0) {
				++args;
				if(args < argc) {
					maxlen = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "# Invalid minimum length parsed\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-mp") == 0) {
				++args;
				if(args < argc) {
					minPhred = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "# Invalid minimum phred score parsed\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-eq") == 0) {
				++args;
				if(args < argc) {
					minQ = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "# Invalid average quality score parsed\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-mq") == 0) {
				++args;
				if(args < argc) {
					mq = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "# Invalid minimum mapping quality parsed\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-5p") == 0) {
				++args;
				if(args < argc) {
					fiveClip = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-5p\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-3p") == 0) {
				++args;
				if(args < argc) {
					threeClip = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-3p\".\n");
						exit(4);
					}
				}
			} else if(strcmp(argv[args], "-dense") == 0) {
				alnToMatPtr = alnToMatDense;
			} else if(strcmp(argv[args], "-sasm") == 0) {
				assembly_KMA_Ptr = &skip_assemble_KMA;
				ID_t = 0.0;
			} else if(strcmp(argv[args], "-matrix") == 0) {
				print_matrix = 1;
			} else if(strcmp(argv[args], "-a") == 0) {
				print_all = 1;
			} else if(strcmp(argv[args], "-ref_fsa") == 0) {
				ref_fsa = 1;
				if(++args < argc && *(argv[args]) != '-') {
					ref_fsa = strtoul(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-ref_fsa\".\n");
						exit(4);
					} else if(ref_fsa == 0) {
						ref_fsa = 2;
					}
				} else {
					--args;
				}
			} else if(strcmp(argv[args], "-Sparse") == 0) {
				sparse_run = 1;
			} else if(strcmp(argv[args], "-1t1") == 0) {
				kmerScan = &save_kmers;
				one2one = 1;
			} else if(strcmp(argv[args], "-ck") == 0) {
				get_kmers_for_pair_ptr = &get_kmers_for_pair_count;
			} else if(strcmp(argv[args], "-hmm") == 0) {
				kmerScan = &save_kmers_HMM;
				one2one = 0;
			} else if(strcmp(argv[args], "-proxi") == 0) {
				if(++args < argc) {
					minFrac = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0 || minFrac < 0 || 1 < minFrac) {
						fprintf(stderr, "Invalid argument at \"-proxi\".\n");
						exit(1);
					} else {
						/* set proximity parameter */
						getMatch = &getProxiMatch;
						getMatchSparse = &getProxiMatchSparse;
						getSecondForce = &getSecondProxiForce;
						getSecondPen = &getSecondProxiPen;
						getF = &getF_Proxi;
						getR = &getR_Proxi;
						getChainTemplates = &getProxiChainTemplates;
						getMatch((int *)(&minFrac), 0);
						getMatchSparse((int *)(&minFrac), 0, 0, 0, 0, 0);
						getSecondPen((int *)(&minFrac), 0, 0, 0, 0, 0, 0, 0);
						getF((int *)(&minFrac), 0, 0, 0, 0);
						ankerAndClean((int *)(&minFrac), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
						ankerAndClean_MEM((int *)(&minFrac), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
						getProxiChainTemplates(0, (const Penalties *)(&minFrac), 0, 0, 0, 0, 0);
					}
				} else {
					fprintf(stderr, "Need argument at: \"-proxi\".\n");
					exit(1);
				}
				
			} else if(strcmp(argv[args], "-ca") == 0) {
				chainSeedsPtr = &chainSeeds_circular;
			} else if(strcmp(argv[args], "-ss") == 0) {
				if(++args < argc) {
					if(argv[args][0] == 'q') {
						ss = 'q';
					} else if(argv[args][0] == 'c') {
						ss = 'c';
					} else if(argv[args][0] == 'd') {
						ss = 'd';
					} else {
						fprintf(stderr, "Invalid argument parsed to option: \"-ss\", using default.\n");
					}
				}
			} else if(strcmp(argv[args], "-p") == 0 || strcmp(argv[args], "-e") == 0) {
				if(++args < argc) {
					evalue = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0 || evalue < 0 || 1.0 < evalue) {
						fprintf(stderr, "Invalid argument at \"%s\".\n", argv[--args]);
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-bc") == 0) {
				if(++args < argc && argv[args][0] != '-') {
					significantBase = &significantAndSupport;
					support = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0 || support < 0 || 1 < support) {
						fprintf(stderr, "Invalid argument at \"-bc\".\n");
						exit(1);
					} else {
						significantAndSupport(0, 0, support);
					}
				} else {
					--args;
					significantBase = &significantNuc;
				}
			} else if(strcmp(argv[args], "-bc90") == 0) {
				significantBase = &significantAnd90Nuc;
			} else if(strcmp(argv[args], "-bcg") == 0) {
				baseCall = &orgBaseCaller;
			} else if(strcmp(argv[args], "-bcNano") == 0) {
				if(significantBase == &significantNuc) {
					significantBase = &significantAnd90Nuc;
				}
				baseCall = &nanoCaller;
			} else if(strcmp(argv[args], "-bcd") == 0) {
				++args;
				if(args < argc) {
					bcd = strtol(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-bcd\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-ID") == 0) {
				++args;
				if(args < argc) {
					ID_t = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-ID\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-mrs") == 0) {
				++args;
				if(args < argc) {
					scoreT = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-mrs\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-mrc") == 0) {
				++args;
				if(args < argc) {
					mrc = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-mrc\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-mct") == 0) {
				++args;
				if(args < argc) {
					coverT = strtod(argv[args], &exeBasic);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-mct\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-reward") == 0) {
				++args;
				if(args < argc) {
					rewards->M = strtol(argv[args], &exeBasic, 10);
					rewards->M = abs(rewards->M);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-reward\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-penalty") == 0) {
				++args;
				if(args < argc) {
					rewards->MM = strtol(argv[args], &exeBasic, 10);
					rewards->MM = MIN(-rewards->MM, rewards->MM);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-penalty\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-gapopen") == 0) {
				++args;
				if(args < argc) {
					rewards->W1 = strtol(argv[args], &exeBasic, 10);
					rewards->W1 = MIN(-rewards->W1, rewards->W1);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-gapopen\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-gapextend") == 0) {
				++args;
				if(args < argc) {
					rewards->U = strtol(argv[args], &exeBasic, 10);
					rewards->U = MIN(-rewards->U, rewards->U);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-gapextend\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-localopen") == 0) {
				/* add to help */
				++args;
				if(args < argc) {
					rewards->Wl = strtol(argv[args], &exeBasic, 10);
					rewards->Wl = MIN(-rewards->Wl, rewards->Wl);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-localopen\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-Npenalty") == 0) {
				/* add to help */
				++args;
				if(args < argc) {
					rewards->Mn = strtol(argv[args], &exeBasic, 10);
					rewards->Mn = MIN(-rewards->Mn, rewards->Mn);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-localopen\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-per") == 0) {
				++args;
				if(args < argc) {
					rewards->PE = strtol(argv[args], &exeBasic, 10);
					rewards->PE = abs(rewards->PE);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-per\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-transition") == 0) {
				/* add to help */
				++args;
				if(args < argc) {
					Ts = strtol(argv[args], &exeBasic, 10);
					Ts = MIN(-Ts, Ts);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-localopen\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-transversion") == 0) {
				/* add to help */
				++args;
				if(args < argc) {
					Tv = strtol(argv[args], &exeBasic, 10);
					Tv = MIN(-Tv, Tv);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-localopen\".\n");
						exit(1);
					}
				}
			} else if(strcmp(argv[args], "-and") == 0) {
				cmp = &cmp_and;
			} else if(strcmp(argv[args], "-oa") == 0) {
				cmp = &cmp_true;
				ID_t = 0.0;
			} else if(strcmp(argv[args], "-boot") == 0) {
				printFsa_ptr = &bootFsa;
			} else if(strcmp(argv[args], "-Mt1") == 0) {
				++args;
				if(args < argc) {
					Mt1 = strtol(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-Mt1\".\n");
						exit(1);
					}
				}
				if(Mt1 < 1) {
					fprintf(stderr, "Invalid template specified at \"-Mt1\"\n");
					exit(1);
				}
				printFsa_ptr = &printFsaMt1;
				printFsa_pair_ptr = &printFsa_pairMt1;
			} else if(strcmp(argv[args], "-ef") == 0) {
				if((args + 1) < argc && *(argv[args + 1]) != '-') {
					++args;
					extendedFeatures = strtol(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-ef\".\n");
						exit(1);
					}
				} else {
					extendedFeatures = 1;
				}
			} else if(strcmp(argv[args], "-vcf") == 0) {
				vcf = 1;
				if(++args < argc) {
					if(argv[args][0] != '-') {
						vcf = strtol(argv[args], &exeBasic, 10);
						if(*exeBasic != 0) {
							fprintf(stderr, "Invalid argument at \"-vcf\".\n");
							exit(1);
						}
					} else {
						--args;
					}
				}
			} else if(strcmp(argv[args], "-xml") == 0) {
				xml = 1;
				if(++args < argc) {
					if(argv[args][0] == '-' && argv[args][1] == '-' && argv[args][2] == 0) {
						xml = 2;
					} else {
						--args;
					}
				}
			} else if(strcmp(argv[args], "-sam") == 0) {
				sam = 1;
				if(++args < argc) {
					if(argv[args][0] != '-') {
						sam = strtol(argv[args], &exeBasic, 10);
						if(*exeBasic != 0) {
							fprintf(stderr, "Invalid argument at \"-sam\".\n");
							exit(1);
						}
					} else {
						--args;
					}
				}
			} else if(strcmp(argv[args], "-nc") == 0) {
				nc = 1;
			} else if(strcmp(argv[args], "-na") == 0) {
				nc = 2;
			} else if(strcmp(argv[args], "-nf") == 0) {
				nf = 1;
			} else if(strcmp(argv[args], "-cge") == 0) {
				scoreT = 0.75;
				rewards->M = 1;
				rewards->MM = -3;
				rewards->W1 = -5;
				rewards->U = -1;
				rewards->PE = 17;
			} else if(strcmp(argv[args], "-tmp") == 0) {
				tmp = 1;
				if(++args < argc) {
					if(argv[args][0] != '-') {
						if(argv[args][strlen(argv[args]) - 1] != '/') {
							fprintf(stderr, "Invalid output directory specified.\n");
							exit(1);
						}
						tmpF(argv[args]);
						tmp = 0;
					} else {
						--args;
					}
				}
			} else if(strcmp(argv[args], "-spltDB") == 0) {
				spltDB = 1;
			} else if(strcmp(argv[args], "-status") == 0) {
				kmaPipe = &kmaPipeFork;
			} else if(strcmp(argv[args], "-verbose") == 0) {
				if(++args < argc && argv[args][0] != '-') {
					verbose = strtol(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-verbose\".\n");
						exit(1);
					}
				} else {
					verbose = 1;
					--args;
				}
			} else if(strcmp(argv[args], "-mint2") == 0) {
				/* equivalent to:
				-1t1, -mem_mode, -ca, -cge, -mq 1, -ref_fsa 2, -dense, 
				-bcg, -bcd 10, -bc 0.9 -vcf -ef */
				kmerScan = &save_kmers;
				one2one = 1;
				mem_mode = 1;
				alignLoadPtr = &alignLoad_fly_mem;
				ankerPtr = &ankerAndClean_MEM;
				chainSeedsPtr = &chainSeeds_circular;
				scoreT = 0.75;
				rewards->M = 1;
				rewards->MM = -3;
				rewards->W1 = -5;
				rewards->U = -1;
				rewards->PE = 17;
				mq = 1;
				ref_fsa = 2;
				alnToMatPtr = alnToMatDense;
				baseCall = &orgBaseCaller;
				bcd = 10;
				significantBase = &significantAndSupport;
				significantAndSupport(0, 0, 0.9);
				vcf = 1;
				extendedFeatures = 1;
			} else if(strcmp(argv[args], "-mint3") == 0) {
				/* equivalent to:
				-1t1, -mem_mode, -ca, -mq 1, -ref_fsa, -dense, 
				-bcNano, -bcd 10, -bc 0.7, -vcf -ef */
				kmerScan = &save_kmers;
				one2one = 1;
				mem_mode = 1;
				alignLoadPtr = &alignLoad_fly_mem;
				ankerPtr = &ankerAndClean_MEM;
				chainSeedsPtr = &chainSeeds_circular;
				mq = 1;
				ref_fsa = 2;
				alnToMatPtr = alnToMatDense;
				baseCall = &nanoCaller;
				bcd = 10;
				significantBase = &significantAndSupport;
				significantAndSupport(0, 0, 0.7);
				vcf = 1;
				extendedFeatures = 1;
			} else if(strcmp(argv[args], "-v") == 0) {
				fprintf(stdout, "KMA-%s\n", KMA_VERSION);
				exit(0);
			} else if(strcmp(argv[args], "-c") == 0) {
				fprintf(stdout, "Philip T.L.C. Clausen, Frank M. Aarestrup & Ole Lund, \"Rapid and precise alignment of raw reads against redundant databases with KMA\", BMC Bioinformatics, 2018;19:307.\n");
				exit(0);
			} else if(strcmp(argv[args], "-h") == 0) {
				helpMessage(0);
			} else {
				fprintf(stderr, " Invalid option:\t%s\n", argv[args]);
				fprintf(stderr, " Printing help message:\n");
				helpMessage(1);
			}
			++args;
		}
		preseed(0, 0, exhaustive);
		trimSeedsPtr(0, ts);
		
		if(sam && kmaPipe != &kmaPipeThread) {
			fprintf(stderr, "\"-sam\" and \"-status\" cannot coincide.\n");
			kmaPipe = &kmaPipeThread;
		}
		
		if(spltDB || targetNum != 1) {
			printPtr = &print_ankers_spltDB;
			if(deConPrintPtr != &deConPrint) {
				deConPrintPtr = printPtr;
			}
			kmerScan = &save_kmers;
			one2one = 1;
		}
		
		if(get_kmers_for_pair_ptr == &get_kmers_for_pair_count) {
			if(one2one) {
				kmerScan = &save_kmers_count;
			}
		}
		
		if(ref_fsa == 1) {
			if(baseCall == nanoCaller) {
				baseCall = &refNanoCaller;
			} else {
				baseCall = &refCaller;
			}
		}
		
		if(outputfilename == 0 || templatefilename == 0) {
			fprintf(stderr, " Too few arguments handed\n");
			fprintf(stderr, " Printing help message:\n");
			helpMessage(1);
		} else if(tmp) {
			/* set tmp files */
			tmpF(outputfilename);
		}
		
		if(fileCounter == 0 && fileCounter_PE == 0 && fileCounter_INT == 0) {
			inputfiles = smalloc(sizeof(char*));
			inputfiles[0] = "--";
			fileCounter = 1;
		}
		
		/* set scoring matrix */
		rewards->MM = (Ts + Tv - 1) / 2; /* avg. of transition and transversion, rounded down */
		d = smalloc(5 * sizeof(int *) + 25 * sizeof(int));
		*d = (int *) (d + 5);
		i = 0;
		while(i < 4) {
			j = 4;
			d[i][j] = rewards->Mn;
			while(j--) {
				d[i][j] = Tv;
			}
			d[i][(i - 2) < 0 ? (i + 2) : (i - 2)] = Ts;
			d[i][i] = rewards->M;
			j = i++;
			d[i] = d[j] + 5;
		}
		i = 5;
		while(i--) {
			d[4][i] = rewards->Mn;
		}
		d[4][4] = 0;
		rewards->d = (int **) d;
		
		if(spltDB && targetNum != 1) {
			/* allocate space for commands */
			escape = 0;
			size = argc + strlen(outputfilename) + 32;
			for(args = 0; args < argc; ++args) {
				if(*argv[args] == '-') {
					escape = 0;
				} else if(escape) {
					size += 2;
				}
				size += strlen(argv[args]);
				if(strncmp(argv[args], "-i", 2) == 0) {
					escape = 1;
				}
			}
			exeBasic = smalloc(size);
			
			fprintf(stderr, "# Map\n");
			for(i = 0; i < targetNum; ++i) {
				to2Bit = exeBasic;
				*to2Bit = 0;
				args = -1;
				while(++args < argc) {
					if(strcmp(argv[args], "-t_db") == 0) {
						escape = 1;
						while(escape && ++args < argc) {
							if(*argv[args] == '-') {
								escape = 0;
							}
						}
						--args;
					} else {
						if(*argv[args] == '-') {
							escape = 0;
						}
						
						if(escape) {
							*to2Bit = '\"';
							++to2Bit;
						}
						
						exe_len = strlen(argv[args]);
						strcpy(to2Bit, argv[args]);
						to2Bit += exe_len;
						
						if(escape) {
							*to2Bit = '\"';
							++to2Bit;
						}
						*to2Bit = ' ';
						++to2Bit;
						
						if(strncmp(argv[args], "-i", 2) == 0) {
							escape = 1;
						}
					}
				}
				fprintf(stdout, "%s-t_db %s -s2 > %s.%d &\n", exeBasic, templatefilenames[i], outputfilename, i);
			}
			
			fprintf(stderr, "# Reduce:\n");
			to2Bit = exeBasic;
			*to2Bit = 0;
			args = -1;
			while(++args < argc) {
				if(strcmp(argv[args], "-spltDB") != 0) {
					if(*argv[args] == '-') {
						escape = 0;
					}
					
					if(escape) {
						*to2Bit = '\"';
						++to2Bit;
					}
					
					exe_len = strlen(argv[args]);
					strcpy(to2Bit, argv[args]);
					to2Bit += exe_len;
					
					if(escape) {
						*to2Bit = '\"';
						++to2Bit;
					}
					*to2Bit = ' ';
					++to2Bit;
					
					if(strncmp(argv[args], "-i", 2) == 0) {
						escape = 1;
					}
				}
			}
			fprintf(stdout, "%s\n", exeBasic);
			
			return 0;
		}
		
		templatefilename = *templatefilenames;
		ioStream = stdout;
		anker_rc(0, 0, one2one, 0, 0, 0);
		anker_rc_comp(0, 0, (unsigned char *)(&one2one), 0, 0, 0, 0, 0);
	} else {
		if(strcmp(*argv, "-s1") == 0) {
			step1 = 1;
		} else if(strcmp(*argv, "-s2") == 0) {
			step2 = 1;
		}
		ioStream = (FILE *) argv[1];
	}
	
	if(step1) {
		t0 = clock();
		/* set to2Bit conversion */
		to2Bit = smalloc(384); /* 128 * 3 = 384 -> OS independent */
		for(i = 0; i < 384; ++i) {
			to2Bit[i] = 8;
		}
		to2Bit += 128;
		to2Bit['\n'] = 16;
		to2Bit['A'] = 0;
		to2Bit['C'] = 1;
		to2Bit['G'] = 2;
		to2Bit['T'] = 3;
		to2Bit['N'] = 4;
		to2Bit['a'] = 0;
		to2Bit['c'] = 1;
		to2Bit['g'] = 2;
		to2Bit['t'] = 3;
		to2Bit['n'] = 4;
		to2Bit['R'] = 0;
		to2Bit['Y'] = 1;
		to2Bit['S'] = 2;
		to2Bit['W'] = 3;
		to2Bit['K'] = 2;
		to2Bit['M'] = 0;
		to2Bit['B'] = 1;
		to2Bit['D'] = 0;
		to2Bit['H'] = 3;
		to2Bit['V'] = 2;
		to2Bit['X'] = 4;
		to2Bit['r'] = 0;
		to2Bit['y'] = 1;
		to2Bit['s'] = 2;
		to2Bit['w'] = 3;
		to2Bit['k'] = 2;
		to2Bit['m'] = 0;
		to2Bit['b'] = 1;
		to2Bit['d'] = 0;
		to2Bit['h'] = 3;
		to2Bit['v'] = 2;
		to2Bit['x'] = 4;
		to2Bit['U'] = 3;
		to2Bit['u'] = 3;
		
		if(sparse_run) {
			templates = smalloc(sizeof(HashMapKMA));
			exe_len = strlen(templatefilename);
			myTemplatefilename = smalloc(exe_len + 64);
			strcpy(myTemplatefilename, templatefilename);
			
			if(deConPrintPtr == deConPrint) {
				strcat(myTemplatefilename, ".decon.comp.b");
			} else {
				strcat(myTemplatefilename, ".comp.b");
			}
			templatefile = sfopen(myTemplatefilename, "rb" );
			loadPrefix(templates, templatefile);
			fclose(templatefile);
			myTemplatefilename[exe_len] = 0;
			kmersize = templates->kmersize;
			if(templates->prefix_len) {
				templates->mask = 0;
				templates->mask = (~templates->mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (templates->prefix_len << 1));
			}
			
			/* merge reads */
			if(fileCounter_PE > 0) {
				inputfiles = realloc(inputfiles, (fileCounter + fileCounter_PE) * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
				for(i = 0; i < fileCounter_PE; ++i, ++fileCounter) {
					inputfiles[fileCounter] = inputfiles_PE[i];
				}
				free(inputfiles_PE);
				fprintf(stderr, "Paired end information is not considered in Sparse mode.\n");
			}
			if(fileCounter_INT > 0) {
				inputfiles = realloc(inputfiles, (fileCounter + fileCounter_INT) * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
				for(i = 0; i < fileCounter_INT; ++i, ++fileCounter) {
					inputfiles[fileCounter] = inputfiles_INT[i];
				}
				free(inputfiles_INT);
				fprintf(stderr, "Interleaved information is not considered in Sparse mode.\n");
			}
			
			run_input_sparse(templates, inputfiles, fileCounter, minPhred, minQ, fiveClip, threeClip, kmersize, to2Bit, prob, ioStream);
			hashMapKMA_destroy(templates);
			free(myTemplatefilename);
		} else {
			if(Mt1) {
				myTemplatefilename = smalloc(strlen(templatefilename) + 64);
				strcpy(myTemplatefilename, templatefilename);
				strcat(myTemplatefilename, ".length.b");
				templatefile = sfopen(myTemplatefilename, "rb");
				fseek(templatefile, (Mt1 + 1) * sizeof(int), SEEK_CUR);
				sfread(&qseq.len, sizeof(int), 1, templatefile);
				fclose(templatefile);
				printFsaMt1(0, &qseq, 0, ioStream);
				printFsa_pairMt1(0, &qseq, 0, 0, 0, ioStream);
			} else {
				myTemplatefilename = 0;
			}
			totFrags = 0;
			
			/* SE */
			if(fileCounter > 0) {
				totFrags += run_input(inputfiles, fileCounter, minPhred, minQ, fiveClip, threeClip, minlen, maxlen, to2Bit, prob, ioStream);
			}
			
			/* PE */
			if(fileCounter_PE > 0) {
				totFrags += run_input_PE(inputfiles_PE, fileCounter_PE, minPhred, minQ, fiveClip, threeClip, minlen, to2Bit, prob, ioStream);
			}
			
			/* INT */
			if(fileCounter_INT > 0) {
				totFrags += run_input_INT(inputfiles_INT, fileCounter_INT, minPhred, minQ, fiveClip, threeClip, minlen, to2Bit, prob, ioStream);
			}
			
			status |= errno;
			
			if(Mt1) {
				Mt1 = -1;
				sfwrite(&Mt1, sizeof(int), 1, ioStream);
				free(myTemplatefilename);
			}
		}
		free((to2Bit - 128));
		if(kmaPipe == &kmaPipeFork) {
			t1 = clock();
			fprintf(stderr, "#\n# Total time used for converting query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
		} else {
			fprintf(stderr, "#\n# Query converted\n#\n");
		}
	} else if(Mt1) {
		myTemplatefilename = smalloc(strlen(templatefilename) + 64);
		strcpy(myTemplatefilename, templatefilename);
		runKMA_Mt1(myTemplatefilename, outputfilename, strjoin(argv, argc), kmersize, minlen, rewards, ID_t, mq, scoreT, mrc, evalue, support, bcd, Mt1, ref_fsa, print_matrix, vcf, xml, sam, nc, nf, thread_num);
		free(myTemplatefilename);
		fprintf(stderr, "# Closing files\n");
	} else if(step2) {
		myTemplatefilename = smalloc(strlen(templatefilename) + 64);
		strcpy(myTemplatefilename, templatefilename);
		status |= save_kmers_batch(myTemplatefilename, "-s1", shm, thread_num, exhaustive, rewards, ioStream, sam, minlen, scoreT, coverT);
		free(myTemplatefilename);
	} else if(sparse_run) {
		myTemplatefilename = smalloc(strlen(templatefilename) + 64);
		strcpy(myTemplatefilename, templatefilename);
		status |= save_kmers_sparse_batch(myTemplatefilename, outputfilename, "-s1", ID_t, evalue, ss, shm);
		free(myTemplatefilename);
		fprintf(stderr, "# Closing files\n");
	} else {
		exeBasic = strjoin(argv, argc);
		myTemplatefilename = smalloc(strlen(templatefilename) + 64);
		strcpy(myTemplatefilename, templatefilename);
		if(spltDB == 0 && targetNum != 1) {
			status |= runKMA_spltDB(templatefilenames, targetNum, outputfilename, argc, argv, ConClave, kmersize, minlen, rewards, extendedFeatures, ID_t, mq, scoreT, mrc, evalue, support, bcd, ref_fsa, print_matrix, print_all, vcf, xml, sam, nc, nf, shm, thread_num, verbose);
		} else if(mem_mode) {
			status |= runKMA_MEM(myTemplatefilename, outputfilename, exeBasic, ConClave, kmersize, minlen, rewards, extendedFeatures, ID_t, mq, scoreT, mrc, evalue, support, bcd, ref_fsa, print_matrix, print_all, vcf, xml, sam, nc, nf, shm, thread_num, verbose);
		} else {
			status |= runKMA(myTemplatefilename, outputfilename, exeBasic, ConClave, kmersize, minlen, rewards, extendedFeatures, ID_t, mq, scoreT, mrc, evalue, support, bcd, ref_fsa, print_matrix, print_all, vcf, xml, sam, nc, nf, shm, thread_num, verbose);
		}
		free(myTemplatefilename);
		fprintf(stderr, "# Closing files\n");
	}
	
	fflush(stderr);
	fflush(stdout);
	status |= errno;
	if(status < 0) {
		status = 1;
	}
	
	return status;
}
