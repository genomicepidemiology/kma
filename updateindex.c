/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
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
#include "compdna.h"
#include "hashmap.h"
#include "hashmapindex.h"
#include "pherror.h"
#include "qualcheck.h"
#include "stdnuc.h"
#include "updateindex.h"

int updateDBs(HashMap *templates, CompDNA *qseq, unsigned template, int MinKlen, double homQ, double homT, unsigned *template_ulengths, unsigned *template_slengths) {
	
	int i, j, end, shifter;
	
	if(qseq->seqlen < templates->kmersize) {
		return 0;
	}
	
	/* set parameters */
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	
	/* iterate sequence */
	for(i = 1, j = 0; i <= qseq->N[0]; ++i) {
		end = qseq->N[i] - templates->kmersize + 1;
		for(;j < end; ++j) {
			/* update hashMap */
			hashMap_add(templates, getKmer(qseq->seq, j, shifter), template);
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	
	return 1;
}

int updateDBs_sparse(HashMap *templates, CompDNA *qseq, unsigned template, int MinKlen, double homQ, double homT, unsigned *template_ulengths, unsigned *template_slengths) {
	
	int i, j, end, rc, prefix_len, prefix_shifter, shifter;
	long unsigned prefix;
	
	if(qseq->seqlen < templates->kmersize) {
		return 0;
	}
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	
	/* test homology and length */
	if(QualCheck(templates, qseq, MinKlen, homQ, homT, template_ulengths)) {
		template_slengths[template] = 0;
		template_ulengths[template] = 0;
		for(rc = 0; rc < 2; ++rc) {
			/* revers complement */
			if(rc) {
				comp_rc(qseq);
			}
			
			/* iterate seq */
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			j = 0;
			if(prefix_len) {
				for(i = 1; i <= qseq->N[0]; ++i) {
					end = qseq->N[i] - prefix_len - templates->kmersize + 1;
					for(;j < end; ++j) {
						if(getKmer(qseq->seq, j, prefix_shifter) == prefix) {
							/* add kmer */
							if(hashMap_add(templates, getKmer(qseq->seq, j + prefix_len, shifter), template)) {
								template_ulengths[template]++;
							}
							template_slengths[template]++;
						}
					}
					j = qseq->N[i] + 1;
				}
				qseq->N[0]--;
			} else {
				for(i = 1; i <= qseq->N[0]; ++i) {
					end = qseq->N[i] - templates->kmersize + 1;
					for(;j < end; ++j) {
						/* add kmer */
						if(hashMap_add(templates, getKmer(qseq->seq, j, shifter), template)) {
							template_ulengths[template]++;
						}
						template_slengths[template]++;
					}
					j = qseq->N[i] + 1;
				}
				qseq->N[0]--;
			}
		}
		return 1;
	}
	
	return 0;
}

void updateAnnots(CompDNA *qseq, int DB_size, int kmerindex, FILE *seq_out, FILE *index_out, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths) {
	
	/* Dump annots */
	dumpIndex(qseq, kmerindex, seq_out, index_out);
	
	(*template_lengths)[DB_size] = qseq->seqlen;
	if((DB_size + 1) >= **template_lengths) {
		**template_lengths <<= 1;
		*template_lengths = realloc(*template_lengths, **template_lengths * sizeof(unsigned));
		if(!*template_lengths) {
			ERROR();
		}
	}
}

void updateAnnots_sparse(CompDNA *qseq, int DB_size, int kmerindex, FILE *seq_out, FILE *index_out, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths) {
	
	/* Dump annots */
	dumpIndex(qseq, kmerindex, seq_out, index_out);
	
	(*template_lengths)[DB_size] = qseq->seqlen;
	if((DB_size + 1) >= **template_ulengths) {
		**template_ulengths <<= 1;
		*template_slengths = realloc(*template_slengths, **template_ulengths * sizeof(unsigned));
		*template_ulengths = realloc(*template_ulengths, **template_ulengths * sizeof(unsigned));
		*template_lengths = realloc(*template_lengths, **template_ulengths * sizeof(unsigned));
		if(!template_lengths || !template_slengths || !template_ulengths) {
			ERROR();
		}
	}
}

void dumpSeq(CompDNA *qseq, int kmerindex, FILE *seq_out, FILE *index_out) {
	cfwrite(qseq->seq, sizeof(long unsigned), (qseq->seqlen >> 5) + 1, seq_out);
}

void makeIndexing(CompDNA *compressor, int kmerindex, FILE *seq_out, FILE *index_out) {
	
	int i, j, end, shifter;
	HashMap_index *template_index;
	
	/* allocate index */
	template_index = smalloc(sizeof(HashMap_index));
	template_index->len = compressor->seqlen;
	template_index->size = compressor->seqlen << 1;
	template_index->kmerindex = kmerindex;
	template_index->index = calloc(template_index->size, sizeof(int));
	if(!template_index->index) {
		ERROR();
	}
	
	/* load index */
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmerindex << 1);
	template_index->seq = compressor->seq;
	compressor->N[0]++;
	compressor->N[compressor->N[0]] = compressor->seqlen + 1;
	j = 0;
	for(i = 1; i <= compressor->N[0]; ++i) {
		end = compressor->N[i] - kmerindex;
		for(;j < end; ++j) {
			hashMapIndex_add(template_index, getKmer(compressor->seq, j, shifter), j);
		}
		j = compressor->N[i] + 1;
	}
	compressor->N[0]--;
	
	/* dump index */
	hashMap_index_dump(template_index, seq_out, index_out);
	
	free(template_index->index);
	free(template_index);
}
