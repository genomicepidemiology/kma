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


#ifndef COMPDNA
typedef struct compDNA CompDNA;
struct compDNA {
	int seqlen;
	int size;
	int complen;
	long unsigned *seq;
	int *N;
};
#define COMPDNA 1
#endif

void allocComp(CompDNA *compressor, int size);
void freeComp(CompDNA *compressor);
void resetComp(CompDNA *compressor);
void compDNA(CompDNA *compressor, unsigned char *seq, int seqlen);
int compDNAref(CompDNA *compressor, unsigned char *qseq, int seqlen);
void unCompDNA(CompDNA *compressor, unsigned char *seq);
long unsigned binRev(long unsigned mer);
void rc_comp(CompDNA *compressor, CompDNA *compressor_rc);
void comp_rc(CompDNA *compressor);
void dumpComp(CompDNA *compressor, FILE* file);
int loadComp(CompDNA *compressor, FILE* file);
int getComp(CompDNA *compressor, FILE* file);
