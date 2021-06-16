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
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "align.h"
#include "alnfrags.h"
#include "ankers.h"
#include "assembly.h"
#include "chain.h"
#include "compdna.h"
#include "ef.h"
#include "filebuff.h"
#include "frags.h"
#include "hashmapcci.h"
#include "kmapipe.h"
#include "nw.h"
#include "penalties.h"
#include "pherror.h"
#include "printconsensus.h"
#include "qseqs.h"
#include "runkma.h"
#include "sam.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "tmp.h"
#include "updatescores.h"
#include "vcf.h"
#include "xml.h"
#ifndef _WIN32
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#else
typedef int key_t;
#define ftok(charPtr, integer) (0)
#define shmget(key, size, permission) ((size != 0) ? (-1) : (-key))
#define shmat(shmid, NULL_Ptr, integer) (NULL)
#define shmdt(dest) fprintf(stderr, "sysV not available on Windows.\n")
#define shmctl(shmid, cmd, buf) fprintf(stderr, "sysV not available on Windows.\n")
#endif
#define nameSkip(infile, c) while((c = fgetc(infile)) != '\n' && c != EOF)


unsigned char * ustrdup(unsigned char *src, size_t n) {
	
	unsigned char *dest;
	
	dest = smalloc(n);
	memcpy(dest, src, n);
	
	return dest;
}

int load_DBs_KMA(char *templatefilename, long unsigned **alignment_scores, long unsigned **uniq_alignment_scores, int **template_lengths, unsigned shm) {
	
	/* load DBs needed for KMA */
	int file_len, shmid, DB_size;
	FILE *DB_file;
	key_t key;
	
	/* allocate DBs */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	DB_file = sfopen(templatefilename, "rb");
	sfread(&DB_size, sizeof(int), 1, DB_file);
	if(shm & 4) {
		key = ftok(templatefilename, 'l');
		shmid = shmget(key, DB_size * sizeof(int), 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared length\n");
			exit(1);
		} else {
			*template_lengths = shmat(shmid, NULL, 0);
		}
	} else {
		*template_lengths = smalloc(DB_size * sizeof(int));
		
		/* load lengths */
		sfread(*template_lengths, sizeof(int), DB_size, DB_file);
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	/* allocate pointers */
	*alignment_scores = calloc(DB_size, sizeof(long unsigned));
	*uniq_alignment_scores = calloc(DB_size, sizeof(long unsigned));
	if(!*alignment_scores || !*uniq_alignment_scores) {
		ERROR();
	}
	
	return DB_size;
}

char * nameLoad(Qseqs *name, FILE *infile) {
	
	int c, size;
	unsigned char *ptr;
	
	ptr = name->seq;
	size = name->size - 1;
	while((c = fgetc(infile)) != '\n' && c != EOF) {
		*ptr++ = c;
		if(--size == 0) {
			size = name->size - 1;
			name->seq = realloc(name->seq, (name->size <<= 1));
			if(!name->seq) {
				ERROR();
			}
			ptr = name->seq + size;
		}
	}
	*ptr = 0;
	
	return (char *) name->seq;
}

int runKMA(char *templatefilename, char *outputfilename, char *exePrev, int ConClave, int kmersize, int minlen, Penalties *rewards, int extendedFeatures, double ID_t, int mq, double scoreT, double mrc, double evalue, double support, int bcd, int ref_fsa, int print_matrix, int print_all, int vcf, int xml, int sam, int nc, int nf, unsigned shm, int thread_num, int verbose) {
	
	int i, j, tmp_template, tmp_tmp_template, file_len, bestTemplate, tot;
	int template, bestHits, t_len, start, end, aln_len, status, rand, sparse;
	int fragCount, fileCount, maxFrag, coverScore, tmp_start, tmp_end, score;
	int seq_in_no, DB_size, flag, counter, stats[5];
	int *bestTemplates, *bestTemplates_r, *best_start_pos, *best_end_pos;
	int *matched_templates, *template_lengths;
	unsigned randScore, *fragmentCounts, *readCounts;
	long read_score, best_read_score, *seq_indexes;
	long unsigned Nhits, template_tot_ulen, bestNum, seqin_size;
	long unsigned *w_scores, *uniq_alignment_scores, *alignment_scores;
	double tmp_score, bestScore, id, q_id, cover, q_cover, p_value;
	long double depth, expected, q_value;
	FILE *inputfile, *frag_in_raw, *res_out, *name_file;
	FILE *alignment_out, *consensus_out, *frag_out_raw, **template_fragments;
	FILE *extendedFeatures_out, *xml_out;
	time_t t0, t1;
	FileBuff *frag_out, *frag_out_all, *matrix_out, *vcf_out;
	Aln *aligned, *gap_align;
	Assem *aligned_assem;
	Frag **alignFrags, *alignFrag;
	CompDNA *qseq_comp, *qseq_r_comp;
	Qseqs *qseq, *qseq_r, *header, *header_r, *template_name;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	Assemble_thread *threads, *thread;
	Aln_thread *alnThreads, *alnThread;
	HashMapCCI **templates_index;
	
	/* get lengths and names */
	file_len = strlen(templatefilename);
	DB_size = load_DBs_KMA(templatefilename, &alignment_scores, &uniq_alignment_scores, &template_lengths, shm);
	templatefilename[file_len] = 0;
	template_name = setQseqs(256);
	strcat(templatefilename, ".name");
	name_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	
	/* print sam-header */
	if(sam) {
		saminit(template_name, name_file, template_lengths, DB_size);
	}
	
	/* open pipe */
	status = 0;
	inputfile = kmaPipe("-s2", "rb", 0, 0);
	if(!inputfile) {
		ERROR();
	} else {
		setvbuf(inputfile, NULL, _IOFBF, CHUNK);
	}
	
	/* load databases */
	strcat(templatefilename, ".seq.b");
	seq_in_no = open(templatefilename, O_RDONLY);
	if(seq_in_no == -1) {
		ERROR();
	}
	seqin_size = 4 * lseek(seq_in_no, 0, SEEK_END);
	if(lseek(seq_in_no, 0, SEEK_SET) != 0) {
		ERROR();
	}
	templatefilename[file_len] = 0;
	templates_index = calloc(DB_size, sizeof(HashMapCCI*));
	if(!templates_index) {
		ERROR();
	}
	alignLoadPtr = &alignLoad_fly;
	if(!kmersize) {
		kmersize = *template_lengths;
	}
	if(kmersize < 4 || 32 < kmersize) {
		kmersize = 16;
	}
	
	/* allocate stuff */
	file_len = strlen(outputfilename);
	seq_indexes = smalloc((DB_size + 1) * sizeof(long));
	/* make file indexes of template indexing */
	*seq_indexes = 0;
	seq_indexes[1] = 0;
	for(i = 2; i < DB_size; ++i) {
		seq_indexes[i] = seq_indexes[i - 1] + ((template_lengths[i - 1] >> 5) + 1) * sizeof(long unsigned);
	}
	
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		if(nf == 0) {
			strcat(outputfilename, ".frag.gz");
			frag_out = gzInitFileBuff(CHUNK);
			openFileBuff(frag_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			frag_out = 0;
		}
		if(nc == 0) {
			strcat(outputfilename, ".aln");
			alignment_out = sfopen(outputfilename, "w");
			outputfilename[file_len] = 0;
			strcat(outputfilename, ".fsa");
			consensus_out = sfopen(outputfilename, "w");
			outputfilename[file_len] = 0;
		} else if(nc == 2) {
			alignment_out = 0;
			strcat(outputfilename, ".fsa");
			consensus_out = sfopen(outputfilename, "w");
			outputfilename[file_len] = 0;
		} else {
			alignment_out = 0;
			consensus_out = 0;
		}
		frag_out_raw = tmpF(0);
		if(!frag_out_raw) {
			ERROR();
		}
		if(print_matrix) {
			matrix_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".mat.gz");
			openFileBuff(matrix_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(print_all) {
			strcat(outputfilename, ".frag_raw.gz");
			frag_out_all = gzInitFileBuff(CHUNK);
			openFileBuff(frag_out_all, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			frag_out_all = 0;
		}
		if(vcf) {
			vcf_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".vcf.gz");
			openFileBuff(vcf_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			vcf_out = 0;
		}
	} else {
		fprintf(stderr, " No output file specified!\n");
		exit(1);
	}
	
	fprintf(stderr, "# Running KMA.\n");
	t0 = clock();
	
	
	/* allocate stuff */
	i = 1;
	alnThreads = 0;
	while(i < thread_num) {
		/* allocate stuff */
		matched_templates = smalloc(((DB_size + 1) << 1) * sizeof(int));
		bestTemplates = smalloc(((DB_size + 1) << 1) * sizeof(int));
		bestTemplates_r = smalloc(((DB_size + 1) << 1) * sizeof(int));
		best_start_pos = smalloc((DB_size << 1) * sizeof(int));
		best_end_pos = smalloc((DB_size << 1) * sizeof(int));
		qseq_comp = smalloc(sizeof(CompDNA));
		qseq_r_comp = smalloc(sizeof(CompDNA));
		allocComp(qseq_comp, 1024);
		allocComp(qseq_r_comp, 1024);
		
		/* allocate matrcies for NW */
		NWmatrices = smalloc(sizeof(NWmat));
		NWmatrices->NW_s = 1024 * 1024;
		NWmatrices->NW_q = 1024;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		NWmatrices->rewards = rewards;
		
		/* move it to the thread */
		alnThread = smalloc(sizeof(Aln_thread));
		alnThread->matched_templates = matched_templates;
		alnThread->bestTemplates = bestTemplates;
		alnThread->bestTemplates_r = bestTemplates_r;
		alnThread->best_start_pos = best_start_pos;
		alnThread->best_end_pos = best_end_pos;
		alnThread->alignment_scores = alignment_scores;
		alnThread->uniq_alignment_scores = uniq_alignment_scores;
		alnThread->seq_indexes = seq_indexes;
		alnThread->inputfile = inputfile;
		alnThread->frag_out_raw = frag_out_raw;
		alnThread->frag_out_all = frag_out_all;
		alnThread->seq_in = seq_in_no;
		alnThread->qseq_comp = qseq_comp;
		alnThread->qseq_r_comp = qseq_r_comp;
		alnThread->qseq = setQseqs(1024);
		alnThread->qseq_r = setQseqs(1024);
		alnThread->header = setQseqs(256);
		alnThread->header_r = setQseqs(256);
		alnThread->points = seedPoint_init(1024, rewards);
		alnThread->NWmatrices = NWmatrices;
		alnThread->kmersize = kmersize;
		alnThread->minlen = minlen;
		alnThread->template_lengths = template_lengths;
		alnThread->templates_index = templates_index;
		alnThread->mq = mq;
		alnThread->sam = sam;//(sam == 1) ? 1 : 0;
		alnThread->scoreT = scoreT;
		alnThread->next = alnThreads;
		alnThreads = alnThread;
		
		
		/* start thread */
		if((errno = pthread_create(&alnThread->id, NULL, &alnFrags_threaded, alnThread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			alnThreads = alnThread->next;
			free(alnThread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* allocate stuff */
	matched_templates = smalloc(((DB_size + 1) << 1) * sizeof(int));
	bestTemplates = smalloc(((DB_size + 1) << 1) * sizeof(int));
	bestTemplates_r = smalloc(((DB_size + 1) << 1) * sizeof(int));
	best_start_pos = smalloc((DB_size << 1) * sizeof(int));
	best_end_pos = smalloc((DB_size << 1) * sizeof(int));
	qseq = setQseqs(1024);
	qseq_r = setQseqs(1024);
	header = setQseqs(256);
	header_r = setQseqs(256);
	qseq_comp = smalloc(sizeof(CompDNA));
	qseq_r_comp = smalloc(sizeof(CompDNA));
	allocComp(qseq_comp, 1024);
	allocComp(qseq_r_comp, 1024);
	points = seedPoint_init(1024, rewards);
	
	/* allocate matrcies for NW */
	NWmatrices = smalloc(sizeof(NWmat));
	NWmatrices->NW_s = 1024 * 1024;
	NWmatrices->NW_q = 1024;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	NWmatrices->rewards = rewards;
	
	/* strat main thread */
	alnThread = smalloc(sizeof(Aln_thread));
	alnThread->matched_templates = matched_templates;
	alnThread->bestTemplates = bestTemplates;
	alnThread->bestTemplates_r = bestTemplates_r;
	alnThread->best_start_pos = best_start_pos;
	alnThread->best_end_pos = best_end_pos;
	alnThread->alignment_scores = alignment_scores;
	alnThread->uniq_alignment_scores = uniq_alignment_scores;
	alnThread->seq_indexes = seq_indexes;
	alnThread->inputfile = inputfile;
	alnThread->frag_out_raw = frag_out_raw;
	alnThread->frag_out_all = frag_out_all;
	alnThread->seq_in = seq_in_no;
	alnThread->qseq_comp = qseq_comp;
	alnThread->qseq_r_comp = qseq_r_comp;
	alnThread->qseq = qseq;
	alnThread->qseq_r = qseq_r;
	alnThread->header = header;
	alnThread->header_r = header_r;
	alnThread->points = points;
	alnThread->NWmatrices = NWmatrices;
	alnThread->kmersize = kmersize;
	alnThread->minlen = minlen;
	alnThread->mq = mq;
	alnThread->sam = sam;//(sam == 1) ? 1 : 0;
	alnThread->scoreT = scoreT;
	alnThread->template_lengths = template_lengths;
	alnThread->templates_index = templates_index;
	alnThread->next = 0;
	
	/* Get alignments */
	alnFrags_threaded(alnThread);
	free(alnThread);
	
	/* join threads */
	for(alnThread = alnThreads; alnThread != 0; alnThread = alnThread->next) {
		/* join thread */
		if((errno = pthread_join(alnThread->id, NULL))) {
			ERROR();
		}
	}
	kmaPipe(0, 0, inputfile, &i);
	status |= i;
	
	i = 0;
	sfwrite(&i, sizeof(int), 1, frag_out_raw);
	fflush(frag_out_raw);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(header->size < header_r->size) {
		free(header->seq);
		header->size = header_r->size;
		header->seq = smalloc(header->size);
	}
	if(qseq->size < qseq_r->size) {
		free(qseq->seq);
		qseq->size = qseq_r->size;
		qseq->seq = smalloc(qseq->size);
	}
	
	/* Patricks features */
	if(extendedFeatures) {
		strcat(outputfilename, ".mapstat");
		extendedFeatures_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		initExtendedFeatures(extendedFeatures_out, templatefilename, *matched_templates, exePrev);
	} else {
		extendedFeatures_out = 0;
	}
	
	if(extendedFeatures || xml) {
		fragmentCounts = calloc(DB_size, sizeof(unsigned));
		readCounts = calloc(DB_size, sizeof(unsigned));
		if(!fragmentCounts || !readCounts) {
			ERROR();
		}
	} else {
		fragmentCounts = 0;
		readCounts = 0;
	}
	
	if(xml) {
		if(xml == 2) {
			xml_out = openInitXML("--", templatefilename, *matched_templates, 1, &exePrev);
		} else {
			strcat(outputfilename, ".xml");
			xml_out = openInitXML(outputfilename, templatefilename, *matched_templates, 1, &exePrev);
			outputfilename[file_len] = 0;
		}
	} else {
		xml_out = 0;
	}
	
	/* clean threads */
	i = 1;
	threads = 0;
	for(alnThread = alnThreads; alnThread != 0; alnThread = alnThreads) {
		alnThreads = alnThread->next;
		
		/* reuse what is possible */
		thread = smalloc(sizeof(Assemble_thread));
		thread->num = i;
		thread->thread_num = thread_num;
		thread->template = -2;
		thread->mq = mq;
		thread->minlen = minlen;
		thread->scoreT = scoreT;
		thread->mrc = mrc;
		thread->evalue = evalue;
		thread->bcd = bcd;
		thread->sam = sam;
		thread->ef = extendedFeatures;
		thread->seq_in = 0;
		thread->kmersize = kmersize;
		thread->frag_out = frag_out;
		thread->xml_out = xml_out;
		thread->NWmatrices = alnThread->NWmatrices;
		thread->qseq = alnThread->qseq;
		thread->header = alnThread->header;
		thread->points = alnThread->points;
		thread->points->len = 0;
		thread->next = threads;
		threads = thread;
		++i;
		
		/* check header and qseq */
		if(header->size < thread->header->size) {
			header->size = thread->header->size;
			free(header->seq);
			header->seq = smalloc(header->size);
		}
		if(qseq->size < thread->qseq->size) {
			qseq->size = thread->qseq->size;
			free(qseq->seq);
			qseq->seq = smalloc(qseq->size);
		}
		
		/* free the rest */
		free(alnThread->matched_templates);
		free(alnThread->bestTemplates);
		free(alnThread->bestTemplates_r);
		free(alnThread->best_start_pos);
		free(alnThread->best_end_pos);
		freeComp(alnThread->qseq_comp);
		free(alnThread->qseq_comp);
		freeComp(alnThread->qseq_r_comp);
		free(alnThread->qseq_r_comp);
		destroyQseqs(alnThread->qseq_r);
		destroyQseqs(alnThread->header_r);
		free(alnThread);
	}
	
	if(kmaPipe == &kmaPipeFork) {
		t1 = clock();
		fprintf(stderr, "#\n# KMA mapping time\t%.2f s.\n", difftime(t1, t0) / 1000000);
	} else {
		fprintf(stderr, "# KMA mapping done\n");
	}
	fprintf(stderr, "#\n# Sort, output and select KMA alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped read
	Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(Frag*));
	w_scores = calloc(DB_size, sizeof(long unsigned));
	if(!alignFrags || !w_scores) {
		ERROR();
	}
	outputfilename[file_len] = 0;
	frag_in_raw = frag_out_raw;
	rewind(frag_in_raw);
	outputfilename[file_len] = 0;
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		ERROR();
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	
	/* Get expected values */
	sparse = 0;
	template_tot_ulen = 0;
	i = DB_size;
	while(--i) {
		template_tot_ulen += template_lengths[i];
	}
	
	/* ConClave */
	if(ConClave == 1) {
		while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			flag = stats[4];
			
			sfread(qseq->seq, 1, qseq->len, frag_in_raw);
			sfread(header->seq, 1, header->len, frag_in_raw);
			sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				start = 0;
				end = 0;
				/* iterate hits */
				for(i = 0; i != bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					tmp_start = best_start_pos[i];
					tmp_end = best_end_pos[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					//if(tmp_score > bestScore) {
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(alignment_scores[tmp_template] == best_read_score) {
					//} else if(tmp_score == bestScore) {
						if(tmp_score > bestScore) {
						//if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						//} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(fragmentCounts) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				sfread(stats, sizeof(int), 3, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				flag = stats[2];
				sfread(qseq->seq, 1, qseq->len, frag_in_raw);
				sfread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = smalloc(sizeof(Frag));
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->buffer[6] = flag;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags, DB_size);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags, DB_size);
		++fileCount;
	} else if(ConClave == 2) {
		/* find potential template candidates */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			/* best templates, skip rest */
			fseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				/* iterate hits */
				for(i = 0; i != bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					//if(tmp_score > bestScore) {
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
					} else if(alignment_scores[tmp_template] == best_read_score) {
					//} else if(tmp_score == bestScore) {
						if(tmp_score > bestScore) {
						//if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						//} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
			}
			w_scores[abs(bestTemplate)] += read_score;
			
			if(stats[2] < 0) {
				sfread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* discard insignifiacant templates */
		Nhits = 0;
		template = DB_size;
		while(--template) {
			Nhits += w_scores[template];
		}
		
		template = DB_size;
		while(--template) {
			if((read_score = w_scores[template])) {
				t_len = template_lengths[template];
				//expected = (Nhits - read_score) * (1.0 * t_len) / (template_tot_ulen - t_len + etta);
				expected = t_len;
				expected /= MAX(1,(template_tot_ulen - t_len));
				expected *= (Nhits - read_score);
				//q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= read_score - expected;
				p_value  = p_chisqr(q_value);
				if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
					w_scores[template] = 0;
				}
			}
		}
		
		/* identify sorting keys */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			if(bestHits != 1) {
				/* best templates, skip rest */
				fseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
				sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
				bestTemplate = 0;
				i = bestHits;
				while(i--) {
					template = abs(bestTemplates[i]);
					if(w_scores[template]) {
						if(bestTemplate) {
							bestTemplate = 0;
							break;
						} else {
							bestTemplate = template;
						}
					}
				}
				
				if(bestTemplate) {
					uniq_alignment_scores[bestTemplate] += read_score;
				}
			} else {
				/* skip rest */
				fseek(frag_in_raw, qseq->len + header->len + 4 * sizeof(int), SEEK_CUR);
			}
			
			if(stats[2] < 0) {
				sfread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* choose the templates */
		memset(w_scores, 0, DB_size * sizeof(long unsigned));
		while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			flag = stats[4];
			
			sfread(qseq->seq, 1, qseq->len, frag_in_raw);
			sfread(header->seq, 1, header->len, frag_in_raw);
			sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best according to sorting keys */
			if(bestHits != 1) {
				
				bestTemplate = 0;
				bestScore = 0;
				start = 0;
				end = 0;
				
				tot = 0;
				i = bestHits;
				while(i--) {
					tot += uniq_alignment_scores[abs(bestTemplates[i])];
				}
				
				if(tot && 16 <= qseq->len) {
					/* get seed */
					rand = qseq->seq[0];
					i = -1;
					j = qseq->len;
					while(++i < 7) {
						rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
					}
					/* minimal standard */
					rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
					if (rand <= 0) {
						rand += 0x7fffffff;
					}
					
					tmp_score = rand;
					tmp_score /= INT_MAX;
					randScore = tmp_score * tot;
					
					score = 0;
					i = 0;
					while(i != bestHits) {
						score += uniq_alignment_scores[abs(bestTemplates[i])];
						if(randScore < score) {
							bestTemplate = bestTemplates[i];
							start = best_start_pos[i];
							end = best_end_pos[i];
							i = bestHits;
						} else {
							++i;
						}
						
					}
					
					if(bestTemplate == 0) {
						tot = 0;
					}
				} else {
					tot = 0;
				}
				
				if(tot == 0) {
					bestTemplate = -1;
					best_read_score = 0;
					bestNum = 0;
					
					/* iterate hits */
					for(i = 0; i != bestHits; ++i) {
						tmp_tmp_template = bestTemplates[i];
						tmp_start = best_start_pos[i];
						tmp_end = best_end_pos[i];
						if(tmp_tmp_template < 0) {
							tmp_template = -tmp_tmp_template;
						} else {
							tmp_template = tmp_tmp_template;
						}
						tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(alignment_scores[tmp_template] == best_read_score) {
						//} else if(tmp_score == bestScore) {
							if(tmp_score > bestScore) {
							//if(alignment_scores[tmp_template] > best_read_score) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
							//} else if(alignment_scores[tmp_template] == best_read_score) {
								if(uniq_alignment_scores[tmp_template] > bestNum) {
									bestTemplate = tmp_tmp_template;
									best_read_score = alignment_scores[tmp_template];
									bestScore = tmp_score;
									bestNum = uniq_alignment_scores[tmp_template];
									start = tmp_start;
									end = tmp_end;
								} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
									bestTemplate = tmp_tmp_template;
									best_read_score = alignment_scores[tmp_template];
									bestScore = tmp_score;
									bestNum = uniq_alignment_scores[tmp_template];
									start = tmp_start;
									end = tmp_end;
								}
							}
						}
					}
				}
				
				
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(fragmentCounts) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				sfread(stats, sizeof(int), 3, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				flag = stats[2];
				sfread(qseq->seq, 1, qseq->len, frag_in_raw);
				sfread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = smalloc(sizeof(Frag));
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->buffer[6] = flag;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags, DB_size);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags, DB_size);
		++fileCount;
	}
	
	fragCount = 0;
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(matched_templates);
	free(bestTemplates);
	destroyQseqs(qseq_r);
	fclose(frag_out_raw);
	if(frag_out_all) {
		destroyGzFileBuff(frag_out_all);
	}
	
	/* Get expected values */
	Nhits = 0;
	i = DB_size;
	while(--i) {
		Nhits += w_scores[i];
	}
	Nhits = Nhits ? Nhits : 1;
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(vcf) {
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* Get expected values */
	points->len = 0;
	
	/* preallocate assembly matrices */
	aligned_assem = smalloc(sizeof(Assem));
	matrix = smalloc(sizeof(AssemInfo));
	matrix->size = qseq->size;
	for(i = 0; i < DB_size; ++i) {
		if(matrix->size < template_lengths[i]) {
			matrix->size = template_lengths[i];
		}
	}
	if(alnToMatPtr == &alnToMat) {
		matrix->size <<= 1;
	} else {
		matrix->size++;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
	/* allocate matrcies for NW */
	for(thread = threads; thread != 0; thread = thread->next) {
		/* get remaining thread info */
		aligned = smalloc(sizeof(Aln));
		gap_align = smalloc(sizeof(Aln));
		aligned->t = smalloc((qseq->size + 1) << 1);
		aligned->s = smalloc((qseq->size + 1) << 1);
		aligned->q = smalloc((qseq->size + 1) << 1);
		gap_align->t = smalloc((qseq->size + 1) << 1);
		gap_align->s = smalloc((qseq->size + 1) << 1);
		gap_align->q = smalloc((qseq->size + 1) << 1);
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->file_count = fileCount;
		thread->files = template_fragments;
		thread->aligned_assem = aligned_assem;
		thread->matrix = matrix;
		thread->spin = (sparse < 0) ? 10 : 100;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
		}
	}
	
	/* start main thread */
	aligned = smalloc(sizeof(Aln));
	gap_align = smalloc(sizeof(Aln));
	aligned->t = smalloc((qseq->size + 1) << 1);
	aligned->s = smalloc((qseq->size + 1) << 1);
	aligned->q = smalloc((qseq->size + 1) << 1);
	gap_align->t = smalloc((qseq->size + 1) << 1);
	gap_align->s = smalloc((qseq->size + 1) << 1);
	gap_align->q = smalloc((qseq->size + 1) << 1);
	thread = smalloc(sizeof(Assemble_thread));
	thread->num = 0;
	thread->thread_num = thread_num;
	thread->template = 0;
	thread->mq = mq;
	thread->minlen = minlen;
	thread->scoreT = scoreT;
	thread->mrc = mrc;
	thread->evalue = evalue;
	thread->bcd = bcd;
	thread->sam = sam;
	thread->ef = extendedFeatures;
	thread->seq_in = 0;
	thread->kmersize = kmersize;
	thread->file_count = fileCount;
	thread->files = template_fragments;
	thread->frag_out = frag_out;
	thread->xml_out = xml_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->next = 0;
	thread->spin = (sparse < 0) ? 10 : 100;
	if(assembly_KMA_Ptr == &skip_assemble_KMA) {
		alignLoadPtr = &alignLoad_skip;
	}
	if(verbose) {
		fprintf(stderr, "# Template\tScore\tProgress\n");
	}
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	counter = 0;
	for(template = 1; template < DB_size; ++template) {
		if(w_scores[template] > 0) {
			if(verbose) {
				counter += w_scores[template];
				fprintf(stderr, "# %d / %d\t%lu\t%3lu%%\n", template, DB_size, w_scores[template], 100 * counter / Nhits);
			}
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = t_len;
			expected /= MAX(1, (template_tot_ulen - t_len));
			expected *= (Nhits - read_score);
			if(0 < expected) {
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= (read_score - expected);
			} else {
				q_value = read_score;
			}
			p_value  = p_chisqr(q_value);
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len))) {
				thread->template_name = nameLoad(template_name, name_file);
				thread->template_index = templates_index[template];
				if(xml) {
					newIterXML(xml_out, template, t_len, thread->template_name);
				}
				/* Do assembly */
				//status |= assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
				thread->template = template;
				thread->t_len = t_len;
				assembly_KMA_Ptr(thread);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = aligned_assem->aln_len;
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					aln_len = 0;
					id = 0;
				}
				
				if(xml) {
					capIterXML(xml_out, DB_size, seqin_size, t_len, readCounts[template], p_value, read_score, aligned_assem->q, aln_len);
				}
				
				if(ID_t <= id && 0 < id) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						thread->template_name, read_score, (unsigned) expected, t_len, id, cover, q_id, q_cover, (double) depth, (double) q_value, p_value);
					if(nc != 1) {
						printConsensus(aligned_assem, thread->template_name, alignment_out, consensus_out, ref_fsa);
					}
					/* print matrix */
					if(matrix_out) {
						updateMatrix(matrix_out, thread->template_name, templates_index[template]->seq, matrix, t_len);
					}
					if(extendedFeatures) {
						printExtendedFeatures(thread->template_name, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
					}
					if(vcf) {
						updateVcf(thread->template_name, aligned_assem->t, evalue, support, bcd, t_len, matrix, vcf, vcf_out);
					}
				}
			} else {
				if((sam && !(sam & 2096)) || ID_t == 0.0) {
					thread->template_index = templates_index[template];
					thread->template_name = nameLoad(template_name, name_file);
					thread->template = template;
					skip_assemble_KMA(thread);
					//skip_assemble_KMA(template, sam, t_len, thread->template_name, fileCount, template_fragments, aligned_assem, qseq, header);
					if(ID_t == 0.0) {
						depth = aligned_assem->depth;
						depth /= t_len;
						aln_len = aligned_assem->aln_len;
						cover = 100.0 * aln_len / t_len;
						q_cover = 100.0 * t_len / aln_len;
						fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
							thread->template_name, read_score, (unsigned) expected, t_len, 0.0, cover, 0.0, q_cover, (double) depth, (double) q_value, p_value);
						if(extendedFeatures) {
							printExtendedFeatures(thread->template_name, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
						}
					}
				} else {
					nameSkip(name_file, end);
				}
			}
		} else {
			nameSkip(name_file, end);
		}
		hashMapCCI_destroy(templates_index[template]);
	}
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	close(seq_in_no);
	fclose(res_out);
	if(alignment_out) {
		fclose(alignment_out);
		fclose(consensus_out);
	}
	fclose(name_file);
	if(frag_out) {
		destroyGzFileBuff(frag_out);
	}
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(extendedFeatures) {
		fclose(extendedFeatures_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	if(xml) {
		closeCapXML(xml_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}

int runKMA_MEM(char *templatefilename, char *outputfilename, char *exePrev, int ConClave, int kmersize, int minlen, Penalties *rewards, int extendedFeatures, double ID_t, int mq, double scoreT, double mrc, double evalue, double support, int bcd, int ref_fsa, int print_matrix, int print_all, int vcf, int xml, int sam, int nc, int nf, unsigned shm, int thread_num, int verbose) {
	
	/* runKMA_MEM is a memory saving version of runKMA,
	   at the cost it chooses best templates based on kmers
	   instead of alignment score. */
	
	int i, j, tmp_template, tmp_tmp_template, file_len, score, rand, rc_flag;
	int template, bestHits, t_len, start, delta, aln_len, fragCount, maxFrag;
	int fileCount, coverScore, tmp_start, tmp_end, bestTemplate, status, tot;
	int sparse, progress, seq_in_no, DB_size, end, flag, flag_r;
	int *matched_templates, *bestTemplates, *best_start_pos, *best_end_pos;
	int *template_lengths, stats[5];
	unsigned randScore, *fragmentCounts, *readCounts;
	long best_read_score, read_score, seq_seeker;
	long unsigned Nhits, template_tot_ulen, bestNum, counter, seqin_size;
	long unsigned *w_scores, *uniq_alignment_scores, *alignment_scores;
	double tmp_score, bestScore, id, cover, q_id, q_cover, p_value;
	long double depth, q_value, expected;
	FILE *inputfile, *frag_in_raw, *res_out, *name_file;
	FILE *alignment_out, *consensus_out, *frag_out_raw, **template_fragments;
	FILE *extendedFeatures_out, *xml_out;
	time_t t0, t1;
	FileBuff *frag_out, *frag_out_all, *matrix_out, *vcf_out;
	Aln *aligned, *gap_align;
	Assem *aligned_assem;
	Frag **alignFrags, *alignFrag;
	CompDNA *qseq_comp, *qseq_r_comp;
	Qseqs *qseq, *qseq_r, *header, *header_r, *template_name;
	AssemInfo *matrix;
	AlnPoints *points;
	NWmat *NWmatrices;
	Assemble_thread *threads, *thread;
	HashMapCCI *template_index;
	
	/* get lengths and names */
	file_len = strlen(templatefilename);
	DB_size = load_DBs_KMA(templatefilename, &alignment_scores, &uniq_alignment_scores, &template_lengths, shm);
	templatefilename[file_len] = 0;
	template_name = setQseqs(256);
	strcat(templatefilename, ".name");
	name_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	if((verbose & 2)) {
		progress = 1;
		verbose = 0;
	} else {
		progress = 0;
	}
	
	/* print sam-header */
	if(sam) {
		saminit(template_name, name_file, template_lengths, DB_size);
	}
	
	/* open pipe */
	status = 0;
	inputfile = kmaPipe("-s2", "rb", 0, 0);
	if(!inputfile) {
		ERROR();
	} else {
		setvbuf(inputfile, NULL, _IOFBF, CHUNK);
	}
	
	/* load databases */
	alignLoadPtr = &alignLoad_fly_mem;
	if(!kmersize) {
		kmersize = *template_lengths;
	}
	if(kmersize < 4 || 32 < kmersize) {
		kmersize = 16;
	}
	strcat(templatefilename, ".seq.b");
	seq_in_no = open(templatefilename, O_RDONLY);
	if(seq_in_no == -1) {
		ERROR();
	}
	seqin_size = 4 * lseek(seq_in_no, 0, SEEK_END);
	if(lseek(seq_in_no, 0, SEEK_SET) != 0) {
		ERROR();
	}
	templatefilename[file_len] = 0;
	
	/* allocate stuff */
	file_len = strlen(outputfilename);
	qseq_comp = malloc(sizeof(CompDNA));
	qseq_r_comp = malloc(sizeof(CompDNA));
	if(!qseq_comp || !qseq_r_comp) {
		ERROR();
	}
	delta = 1024;
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	qseq = setQseqs(delta);
	qseq_r = setQseqs(delta);
	header = setQseqs(256);
	header_r = setQseqs(256);
	points = seedPoint_init(delta, rewards);
	
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		if(nf == 0) {
			strcat(outputfilename, ".frag.gz");
			frag_out = gzInitFileBuff(CHUNK);
			openFileBuff(frag_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			frag_out = 0;
		}
		if(nc == 0) {
			strcat(outputfilename, ".aln");
			alignment_out = sfopen(outputfilename, "w");
			outputfilename[file_len] = 0;
			strcat(outputfilename, ".fsa");
			consensus_out = sfopen(outputfilename, "w");
			outputfilename[file_len] = 0;
		} else if(nc == 2) {
			alignment_out = 0;
			strcat(outputfilename, ".fsa");
			consensus_out = sfopen(outputfilename, "w");
			outputfilename[file_len] = 0;
		} else {
			alignment_out = 0;
			consensus_out = 0;
		}
		frag_out_raw = tmpF(0);
		if(!frag_out_raw) {
			ERROR();
		}
		if(print_matrix) {
			matrix_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".mat.gz");
			openFileBuff(matrix_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(print_all) {
			strcat(outputfilename, ".frag_raw.gz");
			frag_out_all = gzInitFileBuff(CHUNK);
			openFileBuff(frag_out_all, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			frag_out_all = 0;
		}
		if(vcf) {
			vcf_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".vcf.gz");
			openFileBuff(vcf_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			vcf_out = 0;
		}
	} else {
		fprintf(stderr, " No output file specified!\n");
		exit(1);
	}
	
	fprintf(stderr, "# Collecting k-mer scores.\n");
	t0 = clock();
	
	/* Get alignments */
	matched_templates = malloc(((DB_size + 1) << 1) * sizeof(int));
	best_start_pos = calloc((DB_size << 1), sizeof(int));
	best_end_pos = malloc((DB_size << 1) * sizeof(int));
	if(!matched_templates || !best_start_pos || !best_end_pos) {
		ERROR();
	}
	bestTemplates = (matched_templates + 1);
	
	/* consider printPair */
	Nhits = 0;
	t_len = 0;
	read_score = 0;
	while((rc_flag = get_ankers(matched_templates, qseq_comp, header, &flag, inputfile)) != 0) {
		if(*matched_templates) { // SE
			read_score = 0;
			qseq_r->len = 0;
		} else { // PE
			read_score = get_ankers(matched_templates, qseq_r_comp, header_r, &flag_r, inputfile);
			read_score = labs(read_score);
			qseq_r->len = qseq_r_comp->seqlen;
		}
		qseq->len = qseq_comp->seqlen;
		
		if(kmersize <= qseq->len) {
			if(delta <= MAX(qseq->len, qseq_r->len)) {
				delta = MAX(qseq->len, qseq_r->len);
				delta <<= 1;
				qseq->size = delta;
				qseq_r->size = delta;
				free(qseq->seq);
				free(qseq_r->seq);
				qseq->seq = malloc(delta);
				qseq_r->seq = malloc(delta);
				if(!qseq->seq || !qseq_r->seq) {
					ERROR();
				}
			}
			unCompDNA(qseq_comp, qseq->seq);
			
			
			/* reverse complement seq */
			best_read_score = abs(rc_flag);
			
			for(i = 1, bestHits = 0; i <= *matched_templates; ++i, ++bestHits) {
				best_end_pos[bestHits] = template_lengths[abs(matched_templates[i])];
			}
			
			if(rc_flag < 0 && 0 < matched_templates[*matched_templates]) {
				bestHits = -bestHits;
			}
			if(read_score && kmersize <= qseq_r->len) {
				unCompDNA(qseq_r_comp, qseq_r->seq);
				update_Scores_pe(qseq->seq, qseq->len, qseq_r->seq, qseq_r->len, bestHits, best_read_score + read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, flag, flag_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			} else {
				update_Scores(qseq->seq, qseq->len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, flag, alignment_scores, uniq_alignment_scores, frag_out_raw);
			}
			
			/* dump seq to all */
			if(frag_out_all) {
				updateAllFrag(qseq->seq, qseq->len, abs(bestHits), best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_all);
				if(read_score) {
					updateAllFrag(qseq_r->seq, qseq_r->len, abs(bestHits), read_score, best_start_pos, best_end_pos, bestTemplates, header_r, frag_out_all);
				}
			}
			
			/* verbose */
			if(verbose && verbose++ == 1024) {
				Nhits += verbose - 1;
				fprintf(stderr, "# Scored %ld query sequences.\n", Nhits);
				verbose = 1;
			}
		}
	}
	
	/* verbose */
	if(verbose) {
		Nhits += verbose - 1;
		fprintf(stderr, "# Scored %ld query sequences.\n", Nhits);
		verbose = 1;
	}
	kmaPipe(0, 0, inputfile, &i);
	status |= i;
	i = 0;
	sfwrite(&i, sizeof(int), 1, frag_out_raw);
	fflush(frag_out_raw);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(header->size < header_r->size) {
		destroyQseqs(header);
		header = header_r;
	} else {
		destroyQseqs(header_r);
	}
	header_r = 0;
	if(qseq->size < qseq_r->size) {
		destroyQseqs(qseq);
		qseq = qseq_r;
	} else {
		destroyQseqs(qseq_r);
	}
	qseq_r = 0;
	
	if(frag_out_all) {
		destroyGzFileBuff(frag_out_all);
	}
	if(kmaPipe == &kmaPipeFork) {
		t1 = clock();
		fprintf(stderr, "#\n# Time for score collecting:\t%.2f s.\n", difftime(t1, t0) / 1000000);
	} else {
		fprintf(stderr, "# Score collection done\n");
	}
	fprintf(stderr, "#\n# Sort, output and select k-mer alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(Frag*));
	w_scores = calloc(DB_size, sizeof(long unsigned));
	if(!alignFrags || !w_scores) {
		ERROR();
	}
	outputfilename[file_len] = 0;
	frag_in_raw = frag_out_raw;
	rewind(frag_in_raw);
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		ERROR();
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	
	/* Patricks features */
	if(extendedFeatures) {
		strcat(outputfilename, ".mapstat");
		extendedFeatures_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		initExtendedFeatures(extendedFeatures_out, templatefilename, *matched_templates, exePrev);
	} else {
		extendedFeatures_out = 0;
	}
	if(extendedFeatures || xml) {
		fragmentCounts = calloc(DB_size, sizeof(unsigned));
		readCounts = calloc(DB_size, sizeof(unsigned));
		if(!fragmentCounts || !readCounts) {
			ERROR();
		}
	} else {
		fragmentCounts = 0;
		readCounts = 0;
	}
	
	
	if(xml) {
		if(xml == 2) {
			xml_out = openInitXML("--", templatefilename, *matched_templates, 1, &exePrev);
		} else {
			strcat(outputfilename, ".xml");
			xml_out = openInitXML(outputfilename, templatefilename, *matched_templates, 1, &exePrev);
			outputfilename[file_len] = 0;
		}
	} else {
		xml_out = 0;
	}
	
	/* Get expected values */
	sparse = 0;
	template_tot_ulen = 0;
	i = DB_size;
	while(--i) {
		template_tot_ulen += template_lengths[i];
	}
	
	/* ConClave */
	if(ConClave == 1) {
		while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			flag = stats[4];
			
			sfread(qseq->seq, 1, qseq->len, frag_in_raw);
			sfread(header->seq, 1, header->len, frag_in_raw);
			sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				start = 0;
				end = 0;
				/* iterate hits */
				for(i = 0; i != bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					tmp_start = best_start_pos[i];
					tmp_end = best_end_pos[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					if(tmp_score > bestScore) {
					//if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					//} else if(alignment_scores[tmp_template] == best_read_score) {
					} else if(tmp_score == bestScore) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(fragmentCounts) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				sfread(stats, sizeof(int), 3, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				flag = stats[2];
				sfread(qseq->seq, 1, qseq->len, frag_in_raw);
				sfread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = smalloc(sizeof(Frag));
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->buffer[6] = flag;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags, DB_size);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags, DB_size);
		++fileCount;
	} else if(ConClave == 2) {
		/* find potential template candidates */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			/* best templates, skip rest */
			fseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				/* iterate hits */
				for(i = 0; i != bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					if(tmp_score > bestScore) {
					//if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
					//} else if(alignment_scores[tmp_template] == best_read_score) {
					} else if(tmp_score == bestScore) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
			}
			w_scores[abs(bestTemplate)] += read_score;
			
			if(stats[2] < 0) {
				sfread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* discard insignifiacant templates */
		Nhits = 0;
		template = DB_size;
		while(--template) {
			Nhits += w_scores[template];
		}
		
		template = DB_size;
		while(--template) {
			if((read_score = w_scores[template])) {
				t_len = template_lengths[template];
				//expected = (Nhits - read_score) * (t_len / (template_tot_ulen - t_len + etta));
				expected = t_len;
				expected /= MAX(1, (template_tot_ulen - t_len));
				expected *= (Nhits - read_score);
				//q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= read_score - expected;
				p_value  = p_chisqr(q_value);
				if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
					w_scores[template] = 0;
				}
			}
		}
		
		/* identify sorting keys */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			if(bestHits != 1) {
				/* best templates, skip rest */
				fseek(frag_in_raw, qseq->len + header->len + (2 * bestHits + 1) * sizeof(int), SEEK_CUR);
				sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
				bestTemplate = 0;
				i = bestHits;
				while(i--) {
					template = abs(bestTemplates[i]);
					if(w_scores[template]) {
						if(bestTemplate) {
							bestTemplate = 0;
							break;
						} else {
							bestTemplate = template;
						}
					}
				}
				
				if(bestTemplate) {
					uniq_alignment_scores[bestTemplate] += read_score;
				}
			} else {
				/* skip rest */
				fseek(frag_in_raw, qseq->len + header->len + 4 * sizeof(int), SEEK_CUR);
			}
			
			if(stats[2] < 0) {
				sfread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1] + sizeof(int), SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* choose the templates */
		memset(w_scores, 0, DB_size * sizeof(long unsigned));
		while(fread(stats, sizeof(int), 5, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			flag = stats[4];
			
			sfread(qseq->seq, 1, qseq->len, frag_in_raw);
			sfread(header->seq, 1, header->len, frag_in_raw);
			sfread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			sfread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best according to sorting keys */
			if(bestHits != 1) {
				bestTemplate = 0;
				bestScore = 0;
				start = 0;
				end = 0;
				
				tot = 0;
				i = bestHits;
				while(i--) {
					tot += uniq_alignment_scores[abs(bestTemplates[i])];
				}
				
				if(tot && 16 <= qseq->len) {
					/* get seed */
					rand = qseq->seq[0];
					i = -1;
					j = qseq->len;
					while(++i < 7) {
						rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
					}
					/* minimal standard */
					rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
					if (rand <= 0) {
						rand += 0x7fffffff;
					}
					
					tmp_score = rand;
					tmp_score /= INT_MAX;
					randScore = tmp_score * tot;
					
					score = 0;
					i = 0;
					while(i != bestHits) {
						score += uniq_alignment_scores[abs(bestTemplates[i])];
						if(randScore < score) {
							bestTemplate = bestTemplates[i];
							start = best_start_pos[i];
							end = best_end_pos[i];
							i = bestHits;
						} else {
							++i;
						}
					}
					
					if(bestTemplate == 0) {
						tot = 0;
					}
				} else {
					tot = 0;
				}
				
				if(tot == 0) {
					bestTemplate = -1;
					best_read_score = 0;
					bestNum = 0;
					
					/* iterate hits */
					for(i = 0; i != bestHits; ++i) {
						tmp_tmp_template = bestTemplates[i];
						tmp_start = best_start_pos[i];
						tmp_end = best_end_pos[i];
						if(tmp_tmp_template < 0) {
							tmp_template = -tmp_tmp_template;
						} else {
							tmp_template = tmp_tmp_template;
						}
						tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
						if(tmp_score > bestScore) {
						//if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						//} else if(alignment_scores[tmp_template] == best_read_score) {
						} else if(tmp_score == bestScore) {
							//if(tmp_score > bestScore) {
							if(alignment_scores[tmp_template] > best_read_score) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							//} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
							} else if(alignment_scores[tmp_template] == best_read_score) {
								if(uniq_alignment_scores[tmp_template] > bestNum) {
									bestTemplate = tmp_tmp_template;
									best_read_score = alignment_scores[tmp_template];
									bestScore = tmp_score;
									bestNum = uniq_alignment_scores[tmp_template];
									start = tmp_start;
									end = tmp_end;
								} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
									bestTemplate = tmp_tmp_template;
									best_read_score = alignment_scores[tmp_template];
									bestScore = tmp_score;
									bestNum = uniq_alignment_scores[tmp_template];
									start = tmp_start;
									end = tmp_end;
								}
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(fragmentCounts) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = smalloc(sizeof(Frag));
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->buffer[6] = flag;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				sfread(stats, sizeof(int), 3, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				flag = stats[2];
				sfread(qseq->seq, 1, qseq->len, frag_in_raw);
				sfread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = smalloc(sizeof(Frag));
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->buffer[6] = flag;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags, DB_size);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags, DB_size);
		++fileCount;
	}
	
	fragCount = 0;
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(matched_templates);
	fclose(frag_out_raw);
	
	/* Get expected values */
	Nhits = 0;
	i = DB_size;
	while(--i) {
		Nhits += w_scores[i];
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(vcf) {
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* preallocate assembly matrices */
	matrix = smalloc(sizeof(AssemInfo));
	aligned_assem = smalloc(sizeof(Assem));
	matrix->size = delta;
	for(i = 0; i < DB_size; ++i) {
		if(matrix->size < template_lengths[i]) {
			matrix->size = template_lengths[i];
		}
	}
	
	template_index = smalloc(sizeof(HashMapCCI));
	template_index->size = 0;
	hashMapCCI_initialize(template_index, matrix->size, kmersize);
	/*
	if(alignLoadPtr != alignLoad_fly_shm) {
		hashMapCCI_initialize(template_index, matrix->size, kmersize);
	} else {
		template_index->seq = 0;
		template_index->index = 0;
	}
	*/
	if(alnToMatPtr == &alnToMat) {
		matrix->size <<= 1;
	} else {
		matrix->size++;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(Assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
	/* allocate matrcies for NW */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		/* allocate matrices */
		NWmatrices = smalloc(sizeof(NWmat));
		NWmatrices->NW_s = 1024 * 1024;
		NWmatrices->NW_q = 1024;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		NWmatrices->rewards = rewards;
		
		aligned = smalloc(sizeof(Aln));
		gap_align = smalloc(sizeof(Aln));
		aligned->t = smalloc((delta + 1) << 1);
		aligned->s = smalloc((delta + 1) << 1);
		aligned->q = smalloc((delta + 1) << 1);
		gap_align->t = smalloc((delta + 1) << 1);
		gap_align->s = smalloc((delta + 1) << 1);
		gap_align->q = smalloc((delta + 1) << 1);
		
		/* move it to the thread */
		thread = smalloc(sizeof(Assemble_thread));
		thread->num = i;
		thread->thread_num = thread_num;
		thread->mq = mq;
		thread->minlen = minlen;
		thread->scoreT = scoreT;
		thread->mrc = mrc;
		thread->evalue = evalue;
		thread->bcd = bcd;
		thread->sam = sam;
		thread->ef = extendedFeatures;
		thread->seq_in = seq_in_no;
		thread->kmersize = kmersize;
		thread->template = -2;
		thread->file_count = fileCount;
		thread->files = template_fragments;
		thread->frag_out = frag_out;
		thread->xml_out = xml_out;
		thread->aligned_assem = aligned_assem;
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->NWmatrices = NWmatrices;
		thread->matrix = matrix;
		thread->qseq = setQseqs(qseq->size);
		thread->header = setQseqs(header->size);
		thread->points = seedPoint_init(delta, rewards);
		thread->points->len = 0;
		thread->spin = (sparse < 0) ? 10 : 100;
		thread->template_index = template_index;
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* start main thread */
	NWmatrices = smalloc(sizeof(NWmat));
	NWmatrices->NW_s = 1024 * 1024;
	NWmatrices->NW_q = 1024;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	NWmatrices->rewards = rewards;
	
	aligned = smalloc(sizeof(Aln));
	gap_align = smalloc(sizeof(Aln));
	aligned->t = smalloc((delta + 1) << 1);
	aligned->s = smalloc((delta + 1) << 1);
	aligned->q = smalloc((delta + 1) << 1);
	gap_align->t = smalloc((delta + 1) << 1);
	gap_align->s = smalloc((delta + 1) << 1);
	gap_align->q = smalloc((delta + 1) << 1);
	
	/* move it to the thread */
	thread = smalloc(sizeof(Assemble_thread));
	thread->num = 0;
	thread->thread_num = thread_num;
	thread->mq = mq;
	thread->minlen = minlen;
	thread->scoreT = scoreT;
	thread->mrc = mrc;
	thread->evalue = evalue;
	thread->bcd = bcd;
	thread->sam = sam;
	thread->ef = extendedFeatures;
	thread->seq_in = seq_in_no;
	thread->kmersize = kmersize;
	thread->template = 0;
	thread->file_count = fileCount;
	thread->files = template_fragments;
	thread->frag_out = frag_out;
	thread->xml_out = xml_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->template_index = template_index;
	thread->next = 0;
	thread->spin = (sparse < 0) ? 10 : 100;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	seq_seeker = 0;
	counter = 0;
	if(progress) {
		fprintf(stderr, "# Progress:\t%3d%%\r", 0);
		fflush(stderr);
	}
	if(assembly_KMA_Ptr == &skip_assemble_KMA) {
		alignLoadPtr = &alignLoad_skip;
	}
	if(verbose) {
		fprintf(stderr, "# Template\tScore\tProgress\n");
	}
	for(template = 1; template < DB_size; ++template) {
		if(w_scores[template] > 0) {
			if(progress) {
				counter += w_scores[template];
				fprintf(stderr, "# Progress:\t%3lu%%\r", 100 * counter / Nhits);
				fflush(stderr);
			} else if(verbose) {
				counter += w_scores[template];
				fprintf(stderr, "# %d / %d\t%lu\t%3lu%%\n", template, DB_size, w_scores[template], 100 * counter / Nhits);
			}
			
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = t_len;
			expected /= MAX(1, (template_tot_ulen - t_len));
			expected *= (Nhits - read_score);
			if(0 < expected) {
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= (read_score - expected);
			} else {
				q_value = read_score;
			}
			p_value  = p_chisqr(q_value);
			
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len))) {
				/* load DB */
				seq_seeker *= sizeof(long unsigned);
				lseek(seq_in_no, seq_seeker, SEEK_CUR);
				seq_seeker = 0;
				thread->template_index->len = 0;
				thread->template_name = nameLoad(template_name, name_file);
				
				if(xml) {
					newIterXML(xml_out, template, t_len, thread->template_name);
				}
				
				/* Do assembly */
				//status |= assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
				thread->template = template;
				thread->t_len = t_len;
				assembly_KMA_Ptr(thread);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = aligned_assem->aln_len;
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
					q_id = 0;
					depth = aligned_assem->depth;
					depth /= t_len;
					aln_len = aligned_assem->aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 0;
				}
				
				if(xml) {
					capIterXML(xml_out, DB_size, seqin_size, t_len, readCounts[template], p_value, read_score, aligned_assem->q, aln_len);
				}
				if(ID_t <= id) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						thread->template_name, read_score, (unsigned) expected, t_len, id, cover, q_id, q_cover, (double) depth, (double) q_value, p_value);
					if(nc != 1) {
						printConsensus(aligned_assem, thread->template_name, alignment_out, consensus_out, ref_fsa);
					}
					/* print matrix */
					if(matrix_out) {
						updateMatrix(matrix_out, thread->template_name, thread->template_index->seq, matrix, t_len);
					}
					if(extendedFeatures) {
						printExtendedFeatures(thread->template_name, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
					}
					if(vcf) {
						updateVcf(thread->template_name, aligned_assem->t, evalue, support, bcd, t_len, matrix, vcf, vcf_out);
					}
				}
			} else {
				if((sam && !(sam & 2096)) || ID_t == 0.0) {
					/* load DB */
					seq_seeker *= sizeof(long unsigned);
					lseek(seq_in_no, seq_seeker, SEEK_CUR);
					seq_seeker = 0;
					thread->template_index = alignLoad_skip(thread->template_index, seq_in_no, template_lengths[template], kmersize, 0);
					thread->template_name = nameLoad(template_name, name_file);
					thread->template = template;
					skip_assemble_KMA(thread);
					//skip_assemble_KMA(template, sam, t_len, thread->template_name, fileCount, template_fragments, aligned_assem, qseq, header);
					
					if(ID_t == 0.0) {
						depth = aligned_assem->depth;
						depth /= t_len;
						aln_len = aligned_assem->aln_len;
						cover = 100.0 * aln_len / t_len;
						q_cover = 0;
						fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
							thread->template_name, read_score, (unsigned) expected, t_len, 0.0, cover, 0.0, q_cover, (double) depth, (double) q_value, p_value);
						if(extendedFeatures) {
							printExtendedFeatures(thread->template_name, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
						}
					}
				} else {
					nameSkip(name_file, end);
				}
				seq_seeker += ((template_lengths[template] >> 5) + 1);
			}
		} else {
			nameSkip(name_file, end);
			seq_seeker += ((template_lengths[template] >> 5) + 1);
		}
	}
	if(progress) {
		fprintf(stderr, "\n");
	}
	
	/* clear index */
	hashMapCCI_destroy(thread->template_index);
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	close(seq_in_no);
	fclose(res_out);
	if(alignment_out) {
		fclose(alignment_out);
		fclose(consensus_out);
	}
	fclose(name_file);
	if(frag_out) {
		destroyGzFileBuff(frag_out);
	}
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(extendedFeatures) {
		fclose(extendedFeatures_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	if(xml) {
		closeCapXML(xml_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}
