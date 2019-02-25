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
#include <pthread.h>
#include <stdio.h>
#include "chain.h"
#include "compdna.h"
#include "hashmapindex.h"
#include "qseqs.h"

#ifndef ALNTHREAD
typedef struct aln_thread Aln_thread;
struct aln_thread {
	pthread_t id;
	int *matched_templates;
	int *bestTemplates;
	int *bestTemplates_r;
	int *best_start_pos;
	int *best_end_pos;
	int *template_lengths;
	long unsigned *alignment_scores;
	long unsigned *uniq_alignment_scores;
	long *index_indexes;
	long *seq_indexes;
	FILE *inputfile;
	FILE *frag_out_raw;
	int index_in;
	int seq_in;
	int kmersize;
	int mq;
	double scoreT;
	CompDNA *qseq_comp;
	CompDNA *qseq_r_comp;
	Qseqs *qseq;
	Qseqs *qseq_r;
	Qseqs *header;
	Qseqs *header_r;
	AlnPoints *points;
	NWmat *NWmatrices;
	HashMap_index **templates_index;
	struct aln_thread *next;
};
#define ALNTHREAD 1;
#endif

/* pointer defining how align paired end reds */
void (*alnFragsPE)(HashMap_index**, int*, int*, int, double, CompDNA*, CompDNA*, unsigned char*, unsigned char*, Qseqs*, Qseqs*, int, int*, int*, long unsigned*, long unsigned*, int*, int*, int, int, long*, long*, FILE*, AlnPoints *, NWmat *, volatile int *, volatile int *);
void alnFragsSE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, int rc_flag, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, int q_len, int kmersize, Qseqs *header, int *bestTemplates, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB);
void alnFragsUnionPE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB);
void alnFragsPenaltyPE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB);
void alnFragsForcePE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB);
void * alnFrags_threaded(void * arg);
