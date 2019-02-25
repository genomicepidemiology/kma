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
#include <stdio.h>
#include <stdlib.h>
#include "ankers.h"
#include "compdna.h"
#include "pherror.h"
#include "qseqs.h"

void (*printPtr)(int*, CompDNA*, int, const Qseqs*);
void (*printPairPtr)(int*, CompDNA*, int, const Qseqs*, CompDNA*, int, const Qseqs*);
void (*deConPrintPtr)(int*, CompDNA*, int, const Qseqs*);

void print_ankers(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header) {
	
	int infoSize[6];
	
	infoSize[0] = qseq->seqlen;
	infoSize[1] = qseq->complen;
	infoSize[2] = qseq->N[0];
	infoSize[3] = rc_flag;
	infoSize[4] = *out_Tem;
	infoSize[5] = header->len;
	sfwrite(infoSize, sizeof(int), 6, stdout);
	sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, stdout);
	if(qseq->N[0]) {
		sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], stdout);
	}
	sfwrite(out_Tem + 1, sizeof(int), *out_Tem, stdout);
	sfwrite(header->seq, 1, header->len, stdout);
}

void print_ankers_Sparse(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header) {
	
	int infoSize[6];
	
	infoSize[0] = qseq->seqlen;
	infoSize[1] = qseq->complen;
	infoSize[2] = qseq->N[0];
	infoSize[3] = -(abs(rc_flag));
	infoSize[4] = *out_Tem;
	infoSize[5] = header->len;
	sfwrite(infoSize, sizeof(int), 6, stdout);
	
	sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, stdout);
	sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], stdout);
	sfwrite(out_Tem + 1, sizeof(int), *out_Tem, stdout);
	sfwrite(header->seq, 1, header->len, stdout);
	
}


int find_contamination(int *out_Tem, const int contamination) {
	
	int i;
	
	i = *out_Tem + 1;
	out_Tem += i;
	while(--i) {
		if(*--out_Tem == contamination) {
			return i;
		}
	}
	
	return 0;
}

int find_contamination2(int *out_Tem, const int contamination) {
	
	int i;
	
	i = *out_Tem + 1;
	out_Tem += i;
	while(--i) {
		if(*--out_Tem == contamination) {
			return i;
		} else if(0 < *out_Tem) {
			return 0;
		}
	}
	
	return 0;
}

void deConPrint(int *out_Tem, CompDNA *qseq, int rc_flag, const Qseqs *header) {
	
	int contPos;
	
	if((contPos = find_contamination(out_Tem, out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	if((contPos = find_contamination2(out_Tem, -out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	
	if(0 < *out_Tem) {
		printPtr(out_Tem, qseq, rc_flag, header);
	}
}

void deConPrintPair(int *out_Tem, CompDNA *qseq, int bestScore, const Qseqs *header, CompDNA *qseq_r, int bestScore_r, const Qseqs *header_r) {
	
	int contPos;
	
	if((contPos = find_contamination(out_Tem, out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	if((contPos = find_contamination2(out_Tem, -out_Tem[-3])) != 0) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		--*out_Tem;
	}
	
	if(0 < *out_Tem) {
		contPos = *out_Tem;
		*out_Tem = 0;
		printPtr(out_Tem, qseq, bestScore, header);
		*out_Tem = contPos;
		printPtr(out_Tem, qseq_r, bestScore_r, header_r);
	}
}

void printPair(int *out_Tem, CompDNA *qseq, int bestScore, const Qseqs *header, CompDNA *qseq_r, int bestScore_r, const Qseqs *header_r) {
	
	int contPos;
	
	contPos = *out_Tem;
	*out_Tem = 0;
	printPtr(out_Tem, qseq, bestScore, header);
	*out_Tem = contPos;
	out_Tem[-1]++;
	printPtr(out_Tem, qseq_r, bestScore_r, header_r);
}

int get_ankers(int *out_Tem, CompDNA *qseq, Qseqs *header, FILE *inputfile) {
	
	int infoSize[6];
	
	if(fread(infoSize, sizeof(int), 6, inputfile)) {
		qseq->seqlen = infoSize[0];
		qseq->complen = infoSize[1];
		*out_Tem = infoSize[4];
		header->len = infoSize[5];
		
		/* reallocate */
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
		}
		
		qseq->N[0] = infoSize[2];
		if(header->size <= header->len) {
			free(header->seq);
			header->size = header->len << 1;
			header->seq = malloc(header->size);
			if(!header->seq) {
				ERROR();
			}
		}
		
		fread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		fread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		fread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		fread(header->seq, 1, header->len, inputfile);
	} else {
		infoSize[3] = 0;
	}
	
	/* return score */
	return infoSize[3];
}
