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
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "align.h"
#include "assembly.h"
#include "chain.h"
#include "filebuff.h"
#include "hashmapindex.h"
#include "kmapipe.h"
#include "nw.h"
#include "pherror.h"
#include "qseqs.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"

void * (*assembly_KMA_Ptr)(void *);
int (*significantBase)(int, int, double);
unsigned char (*baseCall)(unsigned char, unsigned char, int, int, double, Assembly*);

void updateFrags(FileBuff *dest, Qseqs *qseq, Qseqs *header, char *template_name, int *stats) {
	
	int check, avail;
	char *update;
	
	avail = dest->bytes;
	check = 47 + qseq->len + header->len + strlen(template_name);
	
	/* flush buffer */
	if(avail < check) {
		writeGzFileBuff(dest);
		
		/* seq is too big, reallocate buffer */
		if(dest->bytes < check) {
			resetGzFileBuff(dest, check << 1);
		}
		avail = dest->bytes;
	}
	
	/* update buffer with fragment */
	memcpy(dest->next, qseq->seq, qseq->len);
	dest->next += qseq->len;
	avail -= qseq->len;
	/* stats */
	update = (char *) dest->next;
	check = sprintf(update, "\t%d\t%d\t%d\t%d\t%s\t", stats[0], stats[1], stats[2], stats[3], template_name);
	dest->next += check;
	avail -= check;
	/* header */
	header->seq[header->len - 1] = '\n';
	memcpy(dest->next, header->seq, header->len);
	dest->next += header->len;
	avail -= header->len;
	
	dest->bytes = avail;
	
	/* equivalent with:
	fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
	*/
}

void updateMatrix(FileBuff *dest, char *template_name, long unsigned *template_seq, AssemInfo *matrix, int t_len) {
	
	unsigned i, pos, check, avail, asm_len;
	char *update, bases[] = "ACGTN-";
	Assembly *assembly;
	
	/* check buffer capacity */
	check = strlen(template_name) + 2;
	if(dest->bytes < check) {
		writeGzFileBuff(dest);
	}
	update = (char *) dest->next;
	avail = dest->bytes - check;
	
	/* fill in header */
	check -= 2;
	*update++ = '#';
	memcpy(update, template_name, check);
	update += check;
	*update++ = '\n';
	
	/* fill in rows */
	asm_len = matrix->len;
	assembly = matrix->assmb;
	i = 0;
	for(pos = 0; asm_len != 0; --asm_len, pos = assembly[pos].next) {
		/* check buffer capacity */
		if(avail < 38) {
			dest->bytes = avail;
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		
		/* update with row */
		if(pos < t_len) {
			check = sprintf(update, "%c\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", bases[getNuc(template_seq, i)], assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
			++i;
		} else {
			check = sprintf(update, "-\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
		}
		avail -= check;
		update += check;
	}
	
	/* update with last newline */
	if(avail == 0) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	*update++ = '\n';
	dest->next = (unsigned char *) update;
	dest->bytes = avail - 1;
}


int significantNuc(int X, int Y, double evalue) {
	return (p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue && Y < X);
}

int significantAnd90Nuc(int X, int Y, double evalue) {
	return (p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue && (9 * (X + Y) <= 10 * X));
}

int significantAndSupport(int X, int Y, double evalue) {
	
	static double support = 0;
	
	if(support == 0) {
		support = evalue;
	}
	
	return (p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue && (support * (X + Y) <= X));
}

unsigned char baseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = '-';
	} else {
		/* Use MnNemars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
			if(bestNuc == '-' && tNuc != '-') {
				bestNuc = 'n';
			} else {
				bestNuc = tolower(bestNuc);
			}
		}
	}
	
	return bestNuc;
}

unsigned char orgBaseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0 || bestNuc == '-') {
		bestNuc = '-';
	} else if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) { /* McNemars test */
		bestNuc = tolower(bestNuc);
	}
	
	return bestNuc;
}

unsigned char refCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0 || (bestNuc == '-' && tNuc != '-')) {
		bestNuc = 'n';
	} else if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
		bestNuc = tolower(bestNuc);
	}
	
	return bestNuc;
}

unsigned char nanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	int j, bestBaseScore;
	const char bases[] = "ACGTN-";
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = '-';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
			if(bestNuc == '-' && tNuc != '-') {
				bestBaseScore = 0;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < calls->counts[j]) {
						bestBaseScore = calls->counts[j];
						bestNuc = j;
					}
				}
				if(bestBaseScore == 0) {
					bestNuc = '-';
				} else {
					bestNuc = tolower(bases[bestNuc]);
				}
			} else {
				bestNuc = tolower(bestNuc);
			}
		}
	}
	
	return bestNuc;
}

unsigned char refNanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, double evalue, Assembly *calls) {
	
	int j, bestBaseScore;
	const char bases[] = "ACGTN-";
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = 'n';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore, evalue) == 0) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < calls->counts[j]) {
						bestBaseScore = calls->counts[j];
						bestNuc = j;
					}
				}
				if(bestBaseScore == 0) {
					bestNuc = 'n';
				} else {
					bestNuc = tolower(bases[bestNuc]);
				}
			} else {
				bestNuc = tolower(bestNuc);
			}
		} else if(bestNuc == '-') {
			bestNuc = 'n';
		}
	}
	
	return bestNuc;
}

void * assemble_KMA_threaded(void *arg) {
	
	static volatile int excludeIn[1] = {0}, excludeOut[1] = {0}, excludeMatrix[1] = {0}, mainTemplate = -2, thread_wait = 0;
	static char *template_name;
	static HashMap_index *template_index;
	Assemble_thread *thread = arg;
	int i, j, t_len, aln_len, start, end, bias, myBias, gaps, pos, asm_len;
	int read_score, depthUpdate, bestBaseScore, bestScore, template, spin;
	int nextTemplate, file_i, file_count, delta, thread_num, mq, status, bcd;
	int stats[4], buffer[7];
	unsigned coverScore;
	long unsigned depth, depthVar;
	const char bases[] = "ACGTN-";
	double score, scoreT, evalue;
	unsigned char bestNuc;
	FILE **files, *file;
	AlnScore alnStat;
	Assembly *assembly;
	FileBuff *frag_out;
	Assem *aligned_assem;
	Aln *aligned, *gap_align;
	Qseqs *qseq, *header;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	mq = thread->mq;
	scoreT = thread->scoreT;
	evalue = thread->evalue;
	bcd = thread->bcd;
	spin = thread->spin;
	thread_num = thread->thread_num;
	
	if(template != -2) {
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(excludeMatrix);
			mainTemplate = template;
			unlock(excludeMatrix);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(excludeMatrix);
		template_name = thread->template_name;
		template_index = thread->template_index;
		t_len = template_index->len;
		matrix->len = t_len;
		if(matrix->size < (t_len << 1)) {
			matrix->size = (t_len << 1);
			free(matrix->assmb);
			matrix->assmb = malloc(matrix->size * sizeof(Assembly));
			if(!matrix->assmb) {
				ERROR();
			}
		}
		
		/* cpy template seq */
		assembly = matrix->assmb;
		for(i = 0, j = 1; i < t_len; ++i, ++j) {
			assembly[i].counts[0] = 0;
			assembly[i].counts[1] = 0;
			assembly[i].counts[2] = 0;
			assembly[i].counts[3] = 0;
			assembly[i].counts[4] = 0;
			assembly[i].counts[5] = 0;
			assembly[i].next = j;
		}
		/* circularize */
		assembly[t_len - 1].next = 0;
		
		/* start threads */
		aligned_assem->score = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		unlock(excludeMatrix);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(excludeMatrix);
		template = mainTemplate;
		if(template != -1) {
			t_len = template_index->len;
		}
		unlock(excludeMatrix);
		if(template == -1) {
			return NULL;
		}
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			lockTime(excludeIn, spin);
			file = files[file_i];
			if(file != 0) {
				fread(buffer, sizeof(int), 7, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = malloc(qseq->size);
						if(!qseq->seq) {
							ERROR();
						}
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = malloc(header->size);
						if(!header->seq) {
							ERROR();
						}
					}
					fread(qseq->seq, 1, qseq->len, file);
					fread(header->seq, 1, header->len, file);
					unlock(excludeIn);
					
					if(delta < qseq->len) {
						delta = qseq->len << 1;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = smalloc((delta + 1) << 1);
						aligned->s = smalloc((delta + 1) << 1);
						aligned->q = smalloc((delta + 1) << 1);
						gap_align->t = smalloc((delta + 1) << 1);
						gap_align->s = smalloc((delta + 1) << 1);
						gap_align->q = smalloc((delta + 1) << 1);
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template_index, qseq->seq, qseq->len, points)) {
						/* Start with alignment */
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						alnStat = KMA(template_index, qseq->seq, qseq->len, aligned, gap_align, stats[2], MIN(t_len, stats[3]), mq, scoreT, points, NWmatrices);
						
						/* get read score */
						aln_len = alnStat.len;
						start = alnStat.pos;
						end = start + aln_len - alnStat.gaps;
						
						/* Get normed score */
						read_score = alnStat.score;
						if(0 < aln_len) {
							score = 1.0 * read_score / aln_len;
						} else {
							score = 0;
						}
						
						if(0 < read_score && scoreT <= score) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							//lock(excludeMatrix);
							lockTime(excludeMatrix, 10)
							aligned_assem->score += read_score;
							
							/* diff */
							i = 0;
							pos = start;
							assembly = matrix->assmb;
							while(i < aln_len) {
								if(aligned->t[i] == 5) { // Template gap, insertion
									if(t_len <= pos) {
										assembly[pos].counts[aligned->q[i]]++;
										++i;
										pos = assembly[pos].next;
									} else {
										/* get estimate for non insertions */
										myBias = 0;
										for(j = 0; j < 6; ++j) {
											myBias += assembly[pos].counts[j];
										}
										if(myBias > 0) {
											--myBias;
										}
										/* find position of insertion */
										gaps = pos;
										if(pos != 0) {
											--pos;
										} else {
											pos = t_len - 1;
										}
										while(assembly[pos].next != gaps) {
											pos = assembly[pos].next;
										}
										while(i < aln_len && aligned->t[i] == 5) {
											assembly[pos].next = matrix->len++;
											if(matrix->len == matrix->size) {
												matrix->size <<= 1;
												matrix->assmb = realloc(assembly, matrix->size * sizeof(Assembly));
												if(!matrix->assmb) {
													matrix->size >>= 1;
													matrix->size += 1024;
													matrix->assmb = realloc(assembly, matrix->size * sizeof(Assembly));
													if(!matrix->assmb) {
														ERROR();
													}
												}
												assembly = matrix->assmb;
											}
											pos = assembly[pos].next;
											assembly[pos].next = gaps;
											assembly[pos].counts[0] = 0;
											assembly[pos].counts[1] = 0;
											assembly[pos].counts[2] = 0;
											assembly[pos].counts[3] = 0;
											assembly[pos].counts[4] = 0;
											assembly[pos].counts[5] = myBias;
											assembly[pos].counts[aligned->q[i]]++;
											
											++i;
										}
										pos = assembly[pos].next;
									}
								} else if(t_len <= pos) { // Old template gap, not present in this read
									assembly[pos].counts[5]++;
									pos = assembly[pos].next;
								} else {
									assembly[pos].counts[aligned->q[i]]++;
									++i;
									pos = assembly[pos].next;
								}
							}
							
							unlock(excludeMatrix);
							
							/* Convert fragment */
							for(i = 0; i < qseq->len; ++i) {
								 qseq->seq[i] = bases[qseq->seq[i]];
							}
							qseq->seq[qseq->len] = 0;
							
							/* Save fragment */
							//lock(excludeOut);
							lockTime(excludeOut, 10);
							updateFrags(frag_out, qseq, header, template_name, stats);
							unlock(excludeOut);
							//fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
						}
					}
				} else if(nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						kmaPipe(0, 0, file, &status);
						errno |= status;
					}
					files[file_i] = 0;
					unlock(excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					fseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(excludeIn);
				} else {
					/* Move pointer back */
					fseek(file, (-7) * sizeof(int), SEEK_CUR);
					unlock(excludeIn);
					++file_i;
				}
			} else {
				unlock(excludeIn);
				++file_i;
			}
		}
		lock(excludeIn);
		--thread_wait;
		unlock(excludeIn);
	} while(thread->num != 0);
	
	wait_atomic(thread_wait);
	
	if(aligned_assem->score == 0) {
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		
		return NULL;
	}
	
	/* diff */
	/* pre on dense */
	/* Pepare and make alignment on consensus */
	asm_len = matrix->len;
	assembly = matrix->assmb;
	if(aligned_assem->size <= asm_len) {
		aligned_assem->size = (asm_len + 1) << 1;
		free(aligned_assem->t);
		free(aligned_assem->s);
		free(aligned_assem->q);
		aligned_assem->t = smalloc(aligned_assem->size);
		aligned_assem->s = smalloc(aligned_assem->size);
		aligned_assem->q = smalloc(aligned_assem->size);
	}
	
	/* Call nucleotides for the consensus */
	/* diff */
	i = 0;
	pos = 0;
	depth = 0;
	depthVar = 0;
	aln_len = 0;
	while(i < asm_len) {
		/* call template */
		if(pos < t_len) {
			aligned_assem->t[i] = bases[getNuc(template_index->seq, pos)]; 
		} else {
			aligned_assem->t[i] = '-';
		}
		
		/* call query */
		bestNuc = 5;
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; ++j) {
			if(bestScore < assembly[pos].counts[j]) {
				bestScore = assembly[pos].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[pos].counts[j];
		}
		bestNuc = bases[bestNuc];
		
		/* check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < assembly[pos].counts[j]) {
						bestBaseScore = assembly[pos].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[pos].counts[5];
		}
		
		/* determine base at current position */
		if(bcd <= depthUpdate) {
			bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[pos]);
		} else {
			bestNuc = baseCall('-', aligned_assem->t[i], 0, 0, evalue, &assembly[pos]);
		}
		aligned_assem->q[i] = bestNuc;
		
		if(bestNuc != '-') {
			depth += depthUpdate;
			depthVar += (depthUpdate * depthUpdate);
			++aln_len;
		}
		
		++i;
		pos = assembly[pos].next;
	}
	
	/* Trim alignment on consensus */
	coverScore = 0;
	bias = 0;
	for(i = 0; i < asm_len; ++i) {
		if(aligned_assem->t[i] == '-' && aligned_assem->q[i] == '-') {
			++bias;
		} else {
			aligned_assem->t[i - bias] = aligned_assem->t[i];
			aligned_assem->q[i - bias] = aligned_assem->q[i];
			if(tolower(aligned_assem->t[i]) == tolower(aligned_assem->q[i])) {
				aligned_assem->s[i - bias] = '|';
				++coverScore;
			} else {
				aligned_assem->s[i - bias] = '_';
			}
		}
	}
	asm_len -= bias;
	aligned_assem->t[asm_len] = 0;
	aligned_assem->s[asm_len] = 0;
	aligned_assem->q[asm_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->depthVar = depthVar;
	aligned_assem->len = asm_len;
	aligned_assem->aln_len = aln_len;
	
	return NULL;
}

void * assemble_KMA_dense_threaded(void *arg) {
	
	static volatile int excludeIn[1] = {0}, excludeOut[1] = {0}, excludeMatrix[1] = {0}, mainTemplate = -2, thread_wait = 0;
	static char *template_name;
	static HashMap_index *template_index;
	Assemble_thread *thread = arg;
	int i, j, t_len, aln_len, start, end, file_i, file_count, template, spin;
	int pos, read_score, bestScore, depthUpdate, bestBaseScore, nextTemplate;
	int thread_num, mq, status, bcd, stats[4], buffer[7];
	unsigned coverScore, delta;
	long unsigned depth, depthVar;
	const char bases[] = "ACGTN-";
	double score, scoreT, evalue;
	unsigned char bestNuc;
	FILE **files, *file;
	AlnScore alnStat;
	Assembly *assembly;
	FileBuff *frag_out;
	Assem *aligned_assem;
	Aln *aligned, *gap_align;
	Qseqs *qseq, *header;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	mq = thread->mq;
	scoreT = thread->scoreT;
	evalue = thread->evalue;
	bcd = thread->bcd;
	spin = thread->spin;
	thread_num = thread->thread_num;
	
	if(template != -2) {
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(excludeOut);
			mainTemplate = template;
			unlock(excludeOut);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(excludeOut);
		template_name = thread->template_name;
		template_index = thread->template_index;
		t_len = template_index->len;
		matrix->len = t_len;
		
		/* diff */
		if(aligned_assem->size <= t_len) {
			aligned_assem->size = t_len + 1;
			free(aligned_assem->t);
			free(aligned_assem->s);
			free(aligned_assem->q);
			aligned_assem->t = malloc(t_len + 1);
			aligned_assem->s = malloc(t_len + 1);
			aligned_assem->q = malloc(t_len + 1);
			if(!aligned_assem->t || !aligned_assem->s || !aligned_assem->q) {
				ERROR();
			}
		}
		if(matrix->size <= t_len) {
			matrix->size = t_len + 1;
			free(matrix->assmb);
			matrix->assmb = malloc(matrix->size * sizeof(Assembly));
			if(!matrix->assmb) {
				ERROR();
			}
		}
		
		/* cpy template seq */
		assembly = matrix->assmb;
		for(i = 0, j = 1; i < t_len; ++i, ++j) {
			/* diff */
			aligned_assem->t[i] = getNuc(template_index->seq, i);
			assembly[i].counts[0] = 0;
			assembly[i].counts[1] = 0;
			assembly[i].counts[2] = 0;
			assembly[i].counts[3] = 0;
			assembly[i].counts[4] = 0;
			assembly[i].counts[5] = 0;
			assembly[i].next = j;
		}
		/* circularize */
		assembly[t_len - 1].next = 0;
		
		/* start threads */
		aligned_assem->score = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		unlock(excludeOut);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(excludeOut);
		template = mainTemplate;
		if(template != -1) {
			t_len = template_index->len;
			assembly = matrix->assmb;
		}
		unlock(excludeOut);
		if(template == -1) {
			return NULL;
		}
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			lockTime(excludeIn, spin);
			file = files[file_i];
			if(file != 0) {
				fread(buffer, sizeof(int), 7, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = malloc(qseq->size);
						if(!qseq->seq) {
							ERROR();
						}
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = malloc(header->size);
						if(!header->seq) {
							ERROR();
						}
					}
					fread(qseq->seq, 1, qseq->len, file);
					fread(header->seq, 1, header->len, file);
					unlock(excludeIn);
					
					if(delta < qseq->size) {
						delta = qseq->size;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = malloc((delta + 1) << 1);
						aligned->s = malloc((delta + 1) << 1);
						aligned->q = malloc((delta + 1) << 1);
						gap_align->t = malloc((delta + 1) << 1);
						gap_align->s = malloc((delta + 1) << 1);
						gap_align->q = malloc((delta + 1) << 1);
						if(!aligned->t || !aligned->s || !aligned->q || !gap_align->t || !gap_align->s || !gap_align->q) {
							ERROR();
						}
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template_index, qseq->seq, qseq->len, points)) {
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						/* Start with alignment */
						alnStat = KMA(template_index, qseq->seq, qseq->len, aligned, gap_align, stats[2], MIN(t_len, stats[3]), mq, scoreT, points, NWmatrices);
						
						/* get read score */
						aln_len = alnStat.len;
						start = alnStat.pos;
						end = start + aln_len - alnStat.gaps;
						
						/* Get normed score */
						read_score = alnStat.score;
						if(0 < aln_len) {
							score = 1.0 * read_score / aln_len;
						} else {
							score = 0;
							read_score = 0;
						}
						
						if(0 < read_score && scoreT <= score) {
							
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							//lock(excludeMatrix);
							lockTime(excludeMatrix, 10)
							aligned_assem->score += read_score;
							
							/* diff */
							for(i = 0, pos = start; i < aln_len; ++i) {
								if(aligned->t[i] == aligned_assem->t[pos]) {
									assembly[pos].counts[aligned->q[i]]++;
									pos = assembly[pos].next;
								}
							}
							unlock(excludeMatrix);
							
							/* Convert fragment */
							for(i = 0; i < qseq->len; ++i) {
								 qseq->seq[i] = bases[qseq->seq[i]];
							}
							qseq->seq[qseq->len] = 0;
							
							/* Save fragment */
							//lock(excludeOut);
							lockTime(excludeOut, 10);
							updateFrags(frag_out, qseq, header, template_name, stats);
							unlock(excludeOut);
							
							//fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
						}
					}
				} else if (nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						kmaPipe(0, 0, file, &status);
						errno |= status;
					}
					files[file_i] = 0;
					unlock(excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					fseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(excludeIn);
				} else {
					/* Move pointer back */
					fseek(file, (-7) * sizeof(int), SEEK_CUR);
					unlock(excludeIn);
					++file_i;
				}
			} else {
				unlock(excludeIn);
				++file_i;
			}
		}
		
		lock(excludeIn);
		--thread_wait;
		unlock(excludeIn);
		
	} while(thread->num != 0);
	
	wait_atomic(thread_wait);
	
	if(aligned_assem->score == 0) {
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		
		return NULL;
	}
	
	/* Make consensus assembly by majority voting */
	/* diff */
	depth = 0;
	depthVar = 0;
	coverScore = 0;
	aln_len = 0;
	for(i = 0; i < t_len; ++i) {
		/* call template */
		aligned_assem->t[i] = bases[aligned_assem->t[i]];
		
		/* call query */
		bestNuc = 5;
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; ++j) {
			if(bestScore < assembly[i].counts[j]) {
				bestScore = assembly[i].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[i].counts[j];
		}
		bestNuc = bases[bestNuc];
		
		/* Check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < assembly[i].counts[j]) {
						bestBaseScore = assembly[i].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[i].counts[5];
		}
		
		/* determine base at current position */
		if(bcd <= depthUpdate) {
			bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, evalue, &assembly[i]);
		} else {
			bestNuc = baseCall('-', aligned_assem->t[i], 0, 0, evalue, &assembly[i]);
		}
		aligned_assem->q[i] = bestNuc;
		
		if(bestNuc != '-') {
			depth += depthUpdate;
			depthVar += (depthUpdate * depthUpdate);
			++aln_len;
		}
		
		if(tolower(aligned_assem->q[i]) == tolower(aligned_assem->t[i])) {
			aligned_assem->s[i] = '|';
			++coverScore;
		} else {
			aligned_assem->s[i] = '_';
		}
	}
	aligned_assem->t[t_len] = 0;
	aligned_assem->s[t_len] = 0;
	aligned_assem->q[t_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->depthVar = depthVar;
	aligned_assem->len = t_len;
	aligned_assem->aln_len = aln_len;
	
	return NULL;
}
