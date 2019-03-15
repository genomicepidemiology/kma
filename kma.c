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
#include "kmers.h"
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

static void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# KMA-%s mapps raw reads to a template database.\n", KMA_VERSION);
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-i\t\tInput file name(s)\t\tSTDIN\n");
	fprintf(helpOut, "#\t-ipe\t\tInput paired end file name(s)\n");
	fprintf(helpOut, "#\t-int\t\tInput interleaved file name(s)\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t%s\n", "DB defined");
	fprintf(helpOut, "#\t-e\t\tevalue\t\t\t\t0.05\n");
	fprintf(helpOut, "#\t-ConClave\tConClave version\t\t1\n");
	fprintf(helpOut, "#\t-mem_mode\tUse kmers to choose best\n#\t\t\ttemplate, and save memory\tFalse\n");
	fprintf(helpOut, "#\t-ex_mode\tSearh kmers exhaustively\tFalse\n");
	fprintf(helpOut, "#\t-ef\t\tPrint additional features\tFalse\n");
	fprintf(helpOut, "#\t-vcf\t\tMake vcf file, 2 to apply FT\tFalse/0\n");
	fprintf(helpOut, "#\t-deCon\t\tRemove contamination\t\tFalse\n");
	fprintf(helpOut, "#\t-dense\t\tDo not allow insertions\n#\t\t\tin assembly\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ref_fsa\tConsensus sequnce will\n#\t\t\thave \"n\" instead of gaps\tFalse\n");
	fprintf(helpOut, "#\t-matrix\t\tPrint assembly matrix\t\tFalse\n");
	fprintf(helpOut, "#\t-a\t\tPrint all best mappings\t\tFalse\n");
	fprintf(helpOut, "#\t-mp\t\tMinimum phred score\t\t20\n");
	fprintf(helpOut, "#\t-5p\t\tCut a constant number of\n#\t\t\tnucleotides from the 5 prime.\t0\n");
	fprintf(helpOut, "#\t-Sparse\t\tOnly count kmers\t\tFalse\n");
	fprintf(helpOut, "#\t-Mt1\t\tMap only to \"num\" template.\t0 / False\n");
	fprintf(helpOut, "#\t-ID\t\tMinimum ID\t\t\t1.0%%\n");
	fprintf(helpOut, "#\t-ss\t\tSparse sorting (q,c,d)\t\tq\n");
	fprintf(helpOut, "#\t-pm\t\tPairing method (p,u,f)\t\tu\n");
	fprintf(helpOut, "#\t-fpm\t\tFine Pairing method (p,u,f)\tu\n");
	fprintf(helpOut, "#\t-apm\t\tSets both pm and fpm\t\tu\n");
	fprintf(helpOut, "#\t-shm\t\tUse shared DB made by kma_shm\t0 (lvl)\n");
	//fprintf(helpOut, "#\t-swap\t\tSwap DB to disk\t\t\t0 (lvl)\n");
	fprintf(helpOut, "#\t-1t1\t\tSkip HMM\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ck\t\tCount kmers instead of\n#\t\t\tpseudo alignment\t\tFalse\n");
	fprintf(helpOut, "#\t-ca\t\tMake circular alignments\tFalse\n");
	fprintf(helpOut, "#\t-boot\t\tBootstrap sequence\t\tFalse\n");
	fprintf(helpOut, "#\t-bc\t\tBase calls should be\n#\t\t\tsignificantly overrepresented.\t[True]\n");
	fprintf(helpOut, "#\t-bc90\t\tBase calls should be both\n#\t\t\tsignificantly overrepresented,\n#\t\t\tand have 90%% agreement.\t\tFalse\n");
	fprintf(helpOut, "#\t-bcNano\t\tCall bases at suspicious\n#\t\t\tdeletions, made for nanopore.\tFalse\n");
	fprintf(helpOut, "#\t-bcd\t\tMinimum depth at base\t\t1\n");
	fprintf(helpOut, "#\t-bcg\t\tMaintain insignificant gaps\n");
	fprintf(helpOut, "#\t-and\t\tBoth mrs and p_value thresholds\n#\t\t\thas to reached to in order to\n#\t\t\treport a template hit.\t\tor\n");
	fprintf(helpOut, "#\t-mq\t\tMinimum mapping quality\t\t0\n");
	fprintf(helpOut, "#\t-mrs\t\tMinimum alignment score,\n#\t\t\tnormalized to alignment length\t0.50\n");
	fprintf(helpOut, "#\t-reward\t\tScore for match\t\t\t1\n");
	fprintf(helpOut, "#\t-penalty\tPenalty for mismatch\t\t-2\n");
	fprintf(helpOut, "#\t-gapopen\tPenalty for gap opening\t\t-3\n");
	fprintf(helpOut, "#\t-gapextend\tPenalty for gap extension\t-1\n");
	fprintf(helpOut, "#\t-per\t\tReward for pairing reads\t7\n");
	fprintf(helpOut, "#\t-cge\t\tSet CGE penalties and rewards\tFalse\n");
	fprintf(helpOut, "#\t-t\t\tNumber of threads\t\t1\n");
	fprintf(helpOut, "#\t-v\t\tVersion\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int kma_main(int argc, char *argv[]) {
	
	int i, j, args, exe_len, minPhred, fiveClip, sparse_run, mem_mode, Mt1;
	int step1, step2, fileCounter, fileCounter_PE, fileCounter_INT, status;
	int ConClave, extendedFeatures, vcf, targetNum, size, escape, spltDB, mq;
	int ref_fsa, print_matrix, print_all, one2one, thread_num, kmersize, bcd;
	int **d, W1, U, M, MM, PE;
	unsigned shm, exhaustive;
	long unsigned totFrags;
	char *exeBasic, *outputfilename, *templatefilename, **templatefilenames;
	char **inputfiles, **inputfiles_PE, **inputfiles_INT, *to2Bit;
	char Date[11], ss;
	double ID_t, scoreT, evalue, support;
	FILE *templatefile;
	time_t t0, t1;
	struct tm *tm;
	Qseqs qseq;
	HashMapKMA *templates;
	Penalties *rewards;
	
	if(sizeof(long unsigned) != 8) {
		fprintf(stderr, "Need a 64-bit system.\n");
		exit(3);
	}
	
	/* SET DEFAULTS */
	ConClave = 1;
	totFrags = 0;
	vcf = 0;
	targetNum = 0;
	spltDB = 0;
	extendedFeatures = 0;
	status = 0;
	assembly_KMA_Ptr = &assemble_KMA_threaded;
	cmp = &cmp_or;
	minPhred = 20;
	fiveClip = 0;
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
	evalue = 0.05;
	exhaustive = 0;
	shm = 0;
	mq = 0;
	bcd = 1;
	scoreT = 0.5;
	step1 = 0;
	step2 = 0;
	ID_t = 1.0;
	one2one = 0;
	ss = 'q';
	mem_mode = 0;
	M = 1;
	MM = -2;
	W1 = -3;
	U = -1;
	PE = 7;
	thread_num = 1;
	kmerScan = &save_kmers_HMM;
	get_kmers_for_pair_ptr = &get_kmers_for_pair;
	save_kmers_pair = &save_kmers_unionPair;
	alnFragsPE = &alnFragsUnionPE;
	printPairPtr = &printPair;
	printPtr = &print_ankers;
	printFsa_pair_ptr = &printFsa_pair;
	deConPrintPtr = printPtr;
	ankerPtr = &ankerAndClean;
	alignLoadPtr = &alignLoad_fly;
	destroyPtr = &alignClean;
	printFsa_ptr = &printFsa;
	inputfiles_PE = 0;
	inputfiles_INT = 0;
	inputfiles = 0;
	templatefilenames = 0;
	Mt1 = 0;
	significantBase = &significantNuc; //-bc
	baseCall = &baseCaller;
	chainSeedsPtr = &chainSeeds;
	inputfiles = 0;
	
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
			status = fileCounter;
			for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
				++fileCounter;
			}
			if(fileCounter == 0) {
				fprintf(stderr, "No files were specified.\n");
				exit(3);
			} else {
				inputfiles = realloc(inputfiles, fileCounter * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
			}
			
			for(i = status; i < fileCounter; ++i, ++args) {
				inputfiles[i] = argv[args];
			}
			--args;
		} else if(strcmp(argv[args], "-ipe") == 0) {
			++args;
			status = fileCounter_PE;
			for(i = args; i < argc && strncmp(argv[i], "-", 1) != 0; ++i) {
				++fileCounter_PE;
			}
			if(fileCounter_PE % 2) {
				fprintf(stderr, "Uneven number of paired end files.\n");
				exit(3);
			} else if(fileCounter_PE == 0) {
				fprintf(stderr, "No paired end files were specified.\n");
				exit(3);
			} else {
				inputfiles_PE = realloc(inputfiles_PE, fileCounter_PE * sizeof(char *));
				if(!inputfiles_PE) {
					ERROR();
				}
			}
			
			for(i = status; i < fileCounter_PE; ++i, ++args) {
				inputfiles_PE[i] = argv[args];
			}
			--args;
		} else if(strcmp(argv[args], "-int") == 0) {
			++args;
			status = fileCounter_INT;
			for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
				++fileCounter_INT;
			}
			if(fileCounter_INT == 0) {
				fprintf(stderr, "No interleaved files were specified.\n");
				exit(3);
			}
			inputfiles_INT = realloc(inputfiles_INT, fileCounter_INT * sizeof(char *));
			if(!inputfiles_INT) {
				ERROR();
			}
			for(i = status; i < fileCounter_INT; ++i, ++args) {
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
					exit(4);
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
					exit(4);
				}
			} else {
				--args;
				shm = 3;
			}
		} else if(strcmp(argv[args], "-t") == 0) {
			++args;
			if(args < argc && argv[args][0] != '-') {
				thread_num = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid number of threads specified.\n");
					exit(4);
				}
			} else {
				--args;
			}
			if(thread_num < 1) {
				thread_num = 1;
			}
		} else if(strcmp(argv[args], "-swap") == 0) {
			++args;
			if(!(args < argc && argv[args][0] != '-')) {
				--args;
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
					exit(4);
				} else if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 32) {
					kmersize = 32;
				}
			}
		} else if(strcmp(argv[args], "-mp") == 0) {
			++args;
			if(args < argc) {
				minPhred = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum phred score parsed\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-mq") == 0) {
			++args;
			if(args < argc) {
				mq = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum mapping quality parsed\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-5p") == 0) {
			++args;
			if(args < argc) {
				fiveClip = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-5p\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-dense") == 0) {
			assembly_KMA_Ptr = &assemble_KMA_dense_threaded;
		} else if(strcmp(argv[args], "-matrix") == 0) {
			print_matrix = 1;
		} else if(strcmp(argv[args], "-a") == 0) {
			print_all = 1;
			mem_mode = 1;
			alignLoadPtr = &alignLoad_fly_mem;
			ankerPtr = &ankerAndClean_MEM;
		} else if(strcmp(argv[args], "-ref_fsa") == 0) {
			ref_fsa = 1;
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
		} else if(strcmp(argv[args], "-1t1") == 0) {
			kmerScan = &save_kmers;
			one2one = 1;
		} else if(strcmp(argv[args], "-ck") == 0) {
			get_kmers_for_pair_ptr = &get_kmers_for_pair_count;
		} else if(strcmp(argv[args], "-ca") == 0) {
			chainSeedsPtr = &chainSeeds_circular;
		} else if(strcmp(argv[args], "-ss") == 0) {
			++args;
			if(args < argc) {
				if(argv[args][0] == 'q') {
					ss = 'q';
				} else if(argv[args][0] == 'c') {
					ss = 'c';
				} else if(argv[args][0] == 'd') {
					ss = 'd';
				} else {
					fprintf(stderr, "# Invalid argument parsed to option: \"-ss\", using default.\n");
				}
			}
		} else if(strcmp(argv[args], "-e") == 0) {
			++args;
			if(args < argc) {
				evalue = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-e\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-bc") == 0) {
			if(++args < argc && argv[args][0] != '-') {
				significantBase = &significantAndSupport;
				support = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0 || 1 < support) {
					fprintf(stderr, "Invalid argument at \"-bc\".\n");
					exit(4);
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
					fprintf(stderr, "Invalid argument at \"-ID\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-ID") == 0) {
			++args;
			if(args < argc) {
				ID_t = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-ID\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-mrs") == 0) {
			++args;
			if(args < argc) {
				scoreT = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-mrs\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-reward") == 0) {
			++args;
			if(args < argc) {
				M = strtol(argv[args], &exeBasic, 10);
				M = abs(M);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-reward\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-penalty") == 0) {
			++args;
			if(args < argc) {
				MM = strtol(argv[args], &exeBasic, 10);
				MM = MIN(-MM, MM);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-penalty\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-gapopen") == 0) {
			++args;
			if(args < argc) {
				W1 = strtol(argv[args], &exeBasic, 10);
				W1 = MIN(-W1, W1);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-gapopen\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-gapextend") == 0) {
			++args;
			if(args < argc) {
				U = strtol(argv[args], &exeBasic, 10);
				U = MIN(-U, U);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-gapextend\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-per") == 0) {
			++args;
			if(args < argc) {
				PE = strtol(argv[args], &exeBasic, 10);
				PE = abs(PE);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-per\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-and") == 0) {
			cmp = &cmp_and;
		} else if(strcmp(argv[args], "-boot") == 0) {
			printFsa_ptr = &bootFsa;
		} else if(strcmp(argv[args], "-Mt1") == 0) {
			++args;
			if(args < argc) {
				Mt1 = strtol(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-Mt1\".\n");
					exit(4);
				}
			}
			if(Mt1 < 1) {
				fprintf(stderr, "Invalid template specified at \"-Mt1\"\n");
				exit(3);
			}
			printFsa_ptr = &printFsaMt1;
			printFsa_pair_ptr = &printFsa_pairMt1;
		} else if(strcmp(argv[args], "-ef") == 0) {
			if((args + 1) < argc && *(argv[args + 1]) != '-') {
				++args;
				extendedFeatures = strtol(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-Mt1\".\n");
					exit(4);
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
						exit(4);
					}
				} else {
					--args;
				}
			}
		} else if(strcmp(argv[args], "-cge") == 0) {
			scoreT = 0.75;
			M = 1;
			MM = -3;
			W1 = -5;
			U = -1;
			PE = 17;
		} else if(strcmp(argv[args], "-spltDB") == 0) {
			spltDB = 1;
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA-%s\n", KMA_VERSION);
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
	
	if(ref_fsa) {
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
	}
	
	if(fileCounter == 0 && fileCounter_PE == 0 && fileCounter_INT == 0) {
		inputfiles = malloc(sizeof(char*));
		if(!inputfiles) {
			ERROR();
		}
		inputfiles[0] = "--";
		fileCounter = 1;
	}
	status = 0;
	
	/* set scoring matrix */
	d = smalloc(5 * sizeof(int *));
	for(i = 0; i < 4; ++i) {
		d[i] = smalloc(5 * sizeof(int));
		for(j = 0; j < 4; ++j) {
			d[i][j] = MM;
		}
		d[i][i] = M;
	}
	d[4] = smalloc(5 * sizeof(int));
	for(i = 0; i < 5; ++i) {
		d[4][i] = U;
		d[i][4] = U;
	}
	d[4][4] = 0;
	rewards = smalloc(sizeof(Penalties));
	rewards->d = (int **) d;
	rewards->W1 = W1;
	rewards->U = U;
	rewards->M = M;
	rewards->MM = MM;
	rewards->PE = PE;
	
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
			size += strlen(argv[i]);
			if(strncmp(argv[i], "-i", 2) == 0) {
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
	} else {
		setvbuf(stdout, NULL, _IOFBF, 1048576);
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
			
			if(deConPrintPtr == deConPrint) {
				strcat(templatefilename, ".decon.comp.b");
			} else {
				strcat(templatefilename, ".comp.b");
			}
			templatefile = sfopen(templatefilename, "rb" );
			loadPrefix(templates, templatefile);
			fclose(templatefile);
			templatefilename[exe_len] = 0;
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
			
			run_input_sparse(templates, inputfiles, fileCounter, minPhred, fiveClip, kmersize, to2Bit);
		} else {
			if(Mt1) {
				strcat(templatefilename, ".length.b");
				templatefile = sfopen(templatefilename, "rb");
				fseek(templatefile, (Mt1 + 1) * sizeof(int), SEEK_CUR);
				fread(&qseq.len, sizeof(int), 1, templatefile);
				fclose(templatefile);
				printFsaMt1(0, &qseq, 0);
			}
			kmersize = 16;
			
			/* SE */
			if(fileCounter > 0) {
				totFrags += run_input(inputfiles, fileCounter, minPhred, fiveClip, kmersize, to2Bit);
			}
			
			/* PE */
			if(fileCounter_PE > 0) {
				totFrags += run_input_PE(inputfiles_PE, fileCounter_PE, minPhred, fiveClip, kmersize, to2Bit);
			}
			
			/* INT */
			if(fileCounter_INT > 0) {
				totFrags += run_input_INT(inputfiles_INT, fileCounter_INT, minPhred, fiveClip, kmersize, to2Bit);
			}
			
			if(Mt1) {
				Mt1 = -1;
				sfwrite(&Mt1, sizeof(int), 1, stdout);
			}
			
			if(extendedFeatures && targetNum == 1) {
				strcat(outputfilename, ".mapstat");
				templatefile = sfopen(outputfilename, "wb");
				fprintf(templatefile, "## method\tKMA\n");
				fprintf(templatefile, "## version\t%s\n", KMA_VERSION);
				fprintf(templatefile, "## database %s\n", noFolder(templatefilename));
				fprintf(templatefile, "## fragmentCount\t%lu\n", totFrags);
				time(&t1);
				tm = localtime(&t1);
				strftime(Date, sizeof(Date), "%Y-%m-%d", tm);
				fprintf(templatefile, "## date\t%s\n", Date);
				//fprintf(templatefile, "## date\t%s", ctime(&t1));
				fprintf(templatefile, "## command\t%s", *argv);
				argc -= step1;
				argc -= step2;
				for(args = 1; args < argc; ++args) {
					fprintf(templatefile, " %s", argv[args]);
				}
				fprintf(templatefile, "\n");
				fclose(templatefile);
			}
		}
		fflush(stdout);
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for converting query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else if(Mt1) {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-s1");
		
		runKMA_Mt1(templatefilename, outputfilename, exeBasic, kmersize, rewards, ID_t, mq, scoreT, evalue, bcd, Mt1, ref_fsa, print_matrix, vcf, thread_num);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	} else if(step2) {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-s1");
		status = save_kmers_batch(templatefilename, exeBasic, shm, thread_num, exhaustive, rewards);
		fflush(stdout);
	} else if(sparse_run) {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-s1");
		status = save_kmers_sparse_batch(templatefilename, outputfilename, exeBasic, ID_t, evalue, ss, shm);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	} else {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-s2");
		
		if(spltDB == 0 && targetNum != 1) {
			status = runKMA_spltDB(templatefilenames, targetNum, outputfilename, argc, argv, ConClave, kmersize, rewards, extendedFeatures, ID_t, mq, scoreT, evalue, bcd, ref_fsa, print_matrix, print_all, vcf, shm, thread_num);
		} else if(mem_mode) {
			status = runKMA_MEM(templatefilename, outputfilename, exeBasic, ConClave, kmersize, rewards, extendedFeatures, ID_t, mq, scoreT, evalue, bcd, ref_fsa, print_matrix, print_all, vcf, shm, thread_num);
		} else {
			status = runKMA(templatefilename, outputfilename, exeBasic, ConClave, kmersize, rewards, extendedFeatures, ID_t, mq, scoreT, evalue, bcd, ref_fsa, print_matrix, print_all, vcf, shm, thread_num);
		}
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	}
	
	return status | errno;
}
