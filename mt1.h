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
#include "compdna.h"
#include "penalties.h"
#include "qseqs.h"

void printFsaMt1(Qseqs *header, Qseqs *qseq, Qseqs *qual, CompDNA *compressor, FILE *out);
void printFsa_pairMt1(Qseqs *header, Qseqs *qseq, Qseqs *qual, Qseqs *header_r, Qseqs *qseq_r, Qseqs *qual_r, CompDNA *compressor, FILE *out);
void runKMA_Mt1(char *templatefilename, char *outputfilename, char *exePrev, int kmersize, int minlen, Penalties *rewards, double ID_t, double Depth_t, int mq, double scoreT, double mrc, double evalue, double support, int bcd, int Mt1, int ref_fsa, int print_matrix, long unsigned tsv, int vcf, int xml, int sam, int nc, int nf, int thread_num);
