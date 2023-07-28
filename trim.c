/* Philip T.L.C. Clausen May 2021 plan@dtu.dk */

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
#include <stdlib.h>
#include <stdio.h>
#include "compdna.h"
#include "pherror.h"
#include "runinput.h"
#include "qseqs.h"
#include "trim.h"

void printTrimFsa(Qseqs *header, Qseqs *qseq, Qseqs *qual, CompDNA *compressor, FILE *out) {
	
	const char bases[6] = "ACGTN-";
	static FILE *OUT = 0;
	int i, error;
	unsigned char *ptr;
	
	/* init */
	if(!OUT) {
		OUT = out;
		return;
	}
	error = 0;
	
	/* translate to human readable bases */
	ptr = qseq->seq;
	i = qseq->len + 1;
	while(--i) {
		*ptr = bases[*ptr];
		++ptr;
	}
	
	/* print */
	if(qual) {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		sfwrite("\n+\n", sizeof(unsigned char), 3, OUT);
		sfwrite(qual->seq, sizeof(unsigned char), qual->len, OUT);
		error += (putc('\n', OUT) == EOF);
	} else {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		error += (putc('\n', OUT) == EOF);
	}
	
	if(error) {
		ERROR();
	}
}

void printTrimFsa_pair(Qseqs *header, Qseqs *qseq, Qseqs *qual, Qseqs *header_r, Qseqs *qseq_r, Qseqs *qual_r, CompDNA *compressor, FILE *out) {
	
	const char bases[6] = "ACGTN-";
	static FILE *OUT = 0;
	int i, error;
	unsigned char *ptr;
	
	/* init */
	if(!OUT) {
		OUT = out;
		return;
	}
	error = 0;
	
	/* translate to human readable bases */
	ptr = qseq->seq;
	i = qseq->len + 1;
	while(--i) {
		*ptr = bases[*ptr];
		++ptr;
	}
	ptr = qseq_r->seq;
	i = qseq_r->len + 1;
	while(--i) {
		*ptr = bases[*ptr];
		++ptr;
	}
	
	/* print */
	if(qual) {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		sfwrite("\n+\n", sizeof(unsigned char), 3, OUT);
		sfwrite(qual->seq, sizeof(unsigned char), qual->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(header_r->seq, sizeof(unsigned char), header_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq_r->seq, sizeof(unsigned char), qseq_r->len, OUT);
		sfwrite("\n+\n", sizeof(unsigned char), 3, OUT);
		sfwrite(qual_r->seq, sizeof(unsigned char), qual_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
	} else {
		sfwrite(header->seq, sizeof(unsigned char), header->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq->seq, sizeof(unsigned char), qseq->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(header_r->seq, sizeof(unsigned char), header_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
		sfwrite(qseq_r->seq, sizeof(unsigned char), qseq_r->len, OUT);
		error += (putc('\n', OUT) == EOF);
	}
	
	if(error) {
		ERROR();
	}
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#kma trim trims sequences\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s)", "STDIN");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ipe", "Paired input files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-int", "Inerleaved input file(s)", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "STDOUT");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-oi", "Interleaved output file", "STDOUT");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum length", "16");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-xl", "Maximum length", "2147483647");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mp", "Minimum phred", "20");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mi", "Minimum internal phred score", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-eq", "Minimum average quality", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-5p", "Trim 5 prime", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-3p", "Trim 3 prime", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	
	return (out == stderr);
}

int trim_main(int argc, char *argv[]) {
	
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
	int i, args, fileCounter, fileCounter_PE, fileCounter_INT, fileCount;
	int minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, maxlen;
	long unsigned totFrags;
	char **inputfiles, **inputfiles_PE, **inputfiles_INT;
	char *to2Bit, *outputfilename, *outputfilename_int, *exeBasic;
	FILE *out, *out_int;
	
	/* set defaults */
	minPhred = 20;
	minmaskQ = 0;
	minQ = 0;
	fiveClip = 0;
	threeClip = 0;
	minlen = 16;
	maxlen = 2147483647;
	fileCounter = 0;
	fileCounter_PE = 0;
	fileCounter_INT = 0;
	inputfiles = 0;
	inputfiles_PE = 0;
	inputfiles_INT = 0;
	outputfilename = 0;
	outputfilename_int = 0;
	out = stdout;
	out_int = stdout;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
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
		} else if(strcmp(argv[args], "-o") == 0) {
			++args;
			if(args < argc) {
				outputfilename = smalloc(strlen(argv[args]) + 64);
				strcpy(outputfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-oi") == 0) {
			++args;
			if(args < argc) {
				outputfilename_int = smalloc(strlen(argv[args]) + 64);
				strcpy(outputfilename_int, argv[args]);
			}
		} else if(strcmp(argv[args], "-ml") == 0) {
			++args;
			if(args < argc) {
				minlen = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum length parsed\n");
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
		} else if(strcmp(argv[args], "-mi") == 0) {
			++args;
			if(args < argc) {
				minmaskQ = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid internal minimum phred score parsed\n");
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
		} else if(strcmp(argv[args], "-h") == 0) {
			return helpMessage(stdout);
		} else {
			fprintf(stderr, " Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, " Printing help message:\n");
			return helpMessage(stderr);
		}
		++args;
	}
	
	/* set ptrs */
	if(outputfilename) {
		out = sfopen(outputfilename, "wb");
	}
	printFsa_ptr = &printTrimFsa;
	printFsa_ptr(0, 0, 0, 0, out);
	if(outputfilename_int) {
		out_int = sfopen(outputfilename_int, "wb");
	}
	printFsa_pair_ptr = &printTrimFsa_pair;
	printFsa_pair_ptr(0, 0, 0, 0, 0, 0, 0, out_int);
	
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
	
	/* trim sequences */
	totFrags = 0;
	if(!fileCounter && !fileCounter_PE && !fileCounter_INT) {
		fileCounter = 1;
		inputfiles = smalloc(sizeof(char *));
		*inputfiles = "--";
	}
	if(minPhred < minmaskQ) {
		minPhred = minmaskQ;
	}
	
	/* SE */
	if(fileCounter > 0) {
		totFrags += run_input(inputfiles, fileCounter, minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, maxlen, to2Bit, prob, out);
	}
	
	/* PE */
	if(fileCounter_PE > 0) {
		totFrags += run_input_PE(inputfiles_PE, fileCounter_PE, minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, to2Bit, prob, out);
	}
	
	/* INT */
	if(fileCounter_INT > 0) {
		totFrags += run_input_INT(inputfiles_INT, fileCounter_INT, minPhred, minmaskQ, minQ, fiveClip, threeClip, minlen, to2Bit, prob, out);
	}
	
	fprintf(stderr, "#\n# Total number of query fragment after trimming:\t%lu\n", totFrags);
	
	return 0;
}
