/* Philip T.L.C. Clausen Jan 2017 s123580@student.dtu.dk */

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

char *comp_base;

int chomp(char *string) {
	/* remove trailing spaces and newlines */
	int k = strlen(string) - 1;
	while(string[k] == '\n' || isspace(string[k]))
		k--;
	k++;
	string[k] = '\0';
	return k;
}

int strpos_last(const char* str1, const char* str2) {
	char* strp;
	int i, len1, len2;
	
	len1 = strlen(str1);
	len2 = strlen(str2);
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	strp = (char*)(str1 + len1 - len2);
	for(i = len1 - len2; i >= 0; i--) {
		if(*strp == *str2) {
			if(strncmp(strp,str2,len2)==0)
				return i;
		}
		strp--;
	}
	return -1;
}

void chop_seq(char *header, char *qseq, int seqlen, char *outputdirname, int seqcount, int size, int pair, FILE *statfile) {
	
	if(seqlen == 0) {
		return;
	}
	int i, j, end;
	char *outputfilename, *read;
	FILE *outputfile;
	
	/* 
	* outputfilename = "gzip -c > " + outputdirname + "chopDB" + seqcount + ".fastq.gz" 
	*                   10            strlen()         6         10          10
	*/
	outputfilename = malloc((strlen(outputdirname) + 36) * sizeof(char));
	if(outputfilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	sprintf(outputfilename, "gzip -c > %schopDB%d.fastq.gz", outputdirname, seqcount);
	outputfile = popen(outputfilename, "w");
	if(outputfile == NULL) {
		fprintf(stderr, "File corruption:\t%s\n", outputfilename);
		exit(1);
	}
	end = seqlen - size + 1;
	for(i = 0; i < end; i++) {
		read = qseq + i;
		fprintf(outputfile, "@%d\n", i);
		for(j = 0; j < size; j++) {
			fprintf(outputfile, "%c", read[j]);
		}
		fprintf(outputfile, "\n+\n");
		for(j = 0; j < size; j++) {
			fprintf(outputfile, "%c", 'I');
		}
		fprintf(outputfile, "\n");
	}
	
	pclose(outputfile);
	
	/* print stats */
	fprintf(statfile, "%d\t%d\t%f\t%s\n", seqcount, end, 1.0 * end * size / seqlen, header + 1);
	
	/* clean up */
	free(outputfilename);
	
}

void chop_seq_pair(char *header, char *qseq, int seqlen, char *outputdirname, int seqcount, int size, int pair, FILE *statfile) {
	
	if(seqlen == 0) {
		return;
	}
	int i, j, end;
	char **outputfilename, *read;
	FILE **outputfile;
	
	/* 
	* outputfilename = "gzip -c > " + outputdirname + "chopDB" + seqcount + "_x.fastq.gz" 
	*                   10            strlen()         6         10          12
	*/
	outputfilename = malloc(2 * sizeof(char*));
	if(outputfilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	outputfilename[0] = malloc((strlen(outputdirname) + 38) * sizeof(char));
	outputfilename[1] = malloc((strlen(outputdirname) + 38) * sizeof(char));
	if(outputfilename[0] == NULL || outputfilename[1] == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	sprintf(outputfilename[0], "gzip -c > %schopDB%d_1.fastq.gz", outputdirname, seqcount);
	sprintf(outputfilename[1], "gzip -c > %schopDB%d_2.fastq.gz", outputdirname, seqcount);
	outputfile = malloc(2 * sizeof(FILE*));
	if(outputfile == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	outputfile[0] = popen(outputfilename[0], "w");
	outputfile[1] = popen(outputfilename[1], "w");
	if(outputfile[0] == NULL || outputfile[1] == NULL) {
		fprintf(stderr, "File corruption:\t%s\n", outputfilename[0]);
		exit(1);
	}
	end = seqlen - (2 * size + pair) + 1;
	/* print forward */
	/*for(i = 0; i < size; i++) {
	fprintf(outputfile[0], "@-1/1\n");
	for(j = 0; j < size; j++) {
		fprintf(outputfile[0], "%c", qseq[j]);
	}
	fprintf(outputfile[0], "\n+\n");
	for(j = 0; j < size; j++) {
		fprintf(outputfile[0], "%c", 'I');
	}
	fprintf(outputfile[0], "\n");
	}*/
	for(i = 0; i < end; i++) {
		read = qseq + i;
		fprintf(outputfile[0], "@%d/1\n", i);
		for(j = 0; j < size; j++) {
			fprintf(outputfile[0], "%c", read[j]);
		}
		fprintf(outputfile[0], "\n+\n");
		for(j = 0; j < size; j++) {
			fprintf(outputfile[0], "%c", 'I');
		}
		fprintf(outputfile[0], "\n");
	}
	/* print reverse */
	/*for(i = 0; i < size; i++) {
	read = qseq + (size + pair);
	fprintf(outputfile[1], "@%d/2\n", 0);
	for(j = size - 1; j >= 0; j--) {
		fprintf(outputfile[1], "%c", comp_base[read[j]]);
	}
	fprintf(outputfile[1], "\n+\n");
	for(j = 0; j < size; j++) {
		fprintf(outputfile[1], "%c", 'I');
	}
	fprintf(outputfile[1], "\n");
	}*/
	for(i = 0; i < end; i++) {
		read = qseq + (i + size + pair);
		fprintf(outputfile[1], "@%d/2\n", i);
		/*for(j = 0; j < size; j++) {
			fprintf(outputfile[1], "%c", read[j]);
		}*/
		for(j = size - 1; j >= 0; j--) {
			fprintf(outputfile[1], "%c", comp_base[read[j]]);
		}
		fprintf(outputfile[1], "\n+\n");
		for(j = 0; j < size; j++) {
			fprintf(outputfile[1], "%c", 'I');
		}
		fprintf(outputfile[1], "\n");
		
	}
	
	pclose(outputfile[0]);
	pclose(outputfile[1]);
	
	/* print stats */
	fprintf(statfile, "%d\t%d\t%f\t%s\n", seqcount, 2 * end, 2.0 * end * size / seqlen, header + 1);
	
	/* clean up */
	free(outputfilename[0]);
	free(outputfilename[1]);
	free(outputfilename);
	free(outputfile);
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# chop_DB takes a fasta file and chops it into raw reads.\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "# Options:\t\tDesc:\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\t-i\t\tTemplate DB\t\tstdin\n");
	fprintf(helpOut, "#\t-o\t\tOutput dir\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-size\t\tSize of reads\t\t100\n");
	fprintf(helpOut, "#\t-pair\t\tSpace bestween reads\tFalse / 0\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int args, i, line_size, l_len, seqcount, seqlen, q_size, size, pair;
	volatile int c;
	char *inputfilename, *outputdirname, *statfilename, *line, *header, *qseq;
	FILE *inputfile, *statfile;
	void (*chopSeq_ptr)(char*, char*, int, char*, int, int, int, FILE*);
	
	chopSeq_ptr = &chop_seq;
	inputfilename = 0;
	outputdirname = 0;
	size = 100;
	pair = 0;
	inputfile = stdin;
	
	/* parse cmd-line */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			args++;
			if(args < argc) {
				inputfilename = strdup(argv[args]);
				if(inputfilename == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-o") == 0) {
			args++;
			if(args < argc) {
				outputdirname = strdup(argv[args]);
			}
		} else if(strcmp(argv[args], "-size") == 0) {
			args++;
			if(args < argc) {
				size = atoi(argv[args]);
				if(size <= 0) {
					fprintf(stderr, "Invalid size specified, using default\n");
					size = 100;
				}
			}
		} else if(strcmp(argv[args], "-pair") == 0) {
			args++;
			if(args < argc) {
				pair = atoi(argv[args]);
				if(pair < 0) {
					fprintf(stderr, "Invalid insert size specified, using default\n");
					pair = 0;
				}
			}
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(-1);
		}
		args++;
	}
	if(outputdirname == 0) {
		fprintf(stderr, "# Too few arguments handed\n");
		fprintf(stderr, "# Printing help message:\n");
		helpMessage(-1);
	} else {
		if(outputdirname[strlen(outputdirname) - 1] != '/') {
			outputdirname = realloc(outputdirname, strlen(outputdirname) + 1);
			if(outputdirname == NULL) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			outputdirname[strlen(outputdirname) - 1] = '/';
			outputdirname[strlen(outputdirname)] = '\0';
		}
		mkdir(outputdirname, 0775);
	}
	
	/* open input file */
	if(strcmp(inputfilename + (strlen(inputfilename) - 3), ".gz") == 0) {
		inputfilename = realloc(inputfilename, (strlen(inputfilename) + 11) * sizeof(char));
		if(inputfilename == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(i = strlen(inputfilename); i >= 0; i--) {
			inputfilename[i + 10] = inputfilename[i];
		}
		strncpy(inputfilename, "gunzip -c ", 10);
		inputfile = popen(inputfilename, "r");
	} else {
		inputfile = fopen(inputfilename, "r");
	}
	if(inputfile == NULL) {
		fprintf(stderr, "No such file.\n");
		exit(-1);
	}
	
	/* open statfile */
	if(inputfilename) {
		statfilename = malloc(strlen(outputdirname) + strlen(inputfilename + strpos_last(inputfilename, "/") + 1) + 10);
		if(statfilename == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		sprintf(statfilename, "%s%s_stat.txt", outputdirname, inputfilename + strpos_last(inputfilename, "/") + 1);
	} else {
		statfilename = malloc(strlen(outputdirname) + 9);
		if(statfilename == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		sprintf(statfilename, "%sstat.txt", outputdirname);
	}
	statfile = fopen(statfilename, "w");
	if(statfile == NULL) {
		fprintf(stderr, "File corruption:\t%s\n", statfilename);
		exit(1);
	}
	/* set paired end info */
	if(pair != 0) {
		chopSeq_ptr = &chop_seq_pair;
		comp_base = malloc(128 * sizeof(char));
		if(comp_base == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(i = 0; i < 128; i++) {
			comp_base[i] = 'N';
		}
		comp_base['A'] = 'T';
		comp_base['T'] = 'A';
		comp_base['C'] = 'G';
		comp_base['G'] = 'C';
		comp_base['a'] = 't';
		comp_base['t'] = 'a';
		comp_base['c'] = 'g';
		comp_base['g'] = 'c';
		comp_base['n'] = 'n';
	}
	
	/* parse file */
	seqcount = 0;
	line_size = 1024;
	line = malloc(line_size * sizeof(char));
	header = malloc(line_size * sizeof(char));
	seqlen = 0;
	q_size = 4096;
	qseq = malloc(q_size * sizeof(char));
	if(line == NULL || header == NULL || qseq == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	while(!feof(inputfile) && fgets(line, line_size, inputfile)) {
		l_len = chomp(line);
		if(*line == '>') { //new seq
			//chop_seq(header, qseq, seqlen, outputdirname, seqcount, size, pair, statfile);
			(*chopSeq_ptr)(header, qseq, seqlen, outputdirname, seqcount, size, pair, statfile);
			seqlen = 0;
			seqcount++;
			strcpy(header, line);
			if(l_len == (line_size - 1)) { //Skip rest header
				while((c = fgetc(inputfile)) != EOF && c != '\n');
			}
		} else { //seq
			if((l_len + seqlen) >= q_size) {
				q_size *= 2;
				qseq = realloc(qseq, q_size * sizeof(char));
				if(qseq == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
			}
			strcpy(qseq + seqlen, line);
			seqlen += l_len;
		}
	}
	/* chop last seq */
	//chop_seq(header, qseq, seqlen, outputdirname, seqcount, size, pair, statfile);
	(*chopSeq_ptr)(header, qseq, seqlen, outputdirname, seqcount, size, pair, statfile);
	seqlen = 0;
	
	/* close file */
	if(strcmp(inputfilename + (strlen(inputfilename) - 3), ".gz") == 0) {
		pclose(inputfile);
	} else {
		fclose(inputfile);
	}
	fclose(statfile);
	
	return 0;
}
