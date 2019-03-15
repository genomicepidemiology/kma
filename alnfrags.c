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
#include "align.h"
#include "alnfrags.h"
#include "ankers.h"
#include "chain.h"
#include "compdna.h"
#include "hashmapindex.h"
#include "pherror.h"
#include "qseqs.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"
#include "updatescores.h"

void alnFragsSE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, int rc_flag, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, int q_len, int kmersize, Qseqs *header, int *bestTemplates, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, bestHits, aln_len;
	int start, end, W1;
	double score, bestScore;
	AlnScore alnStat;
	
	/* reverse complement seq */
	if(rc_flag < 0) {
		if(qseq_r_comp->size < qseq_comp->size) {
			freeComp(qseq_r_comp);
			allocComp(qseq_r_comp, qseq_comp->size);
		}
		rc_comp(qseq_comp, qseq_r_comp);
		unCompDNA(qseq_r_comp, qseq_r);
		qseq_r_comp->N[0]++;
		qseq_r_comp->N[qseq_r_comp->N[0]] = q_len;
	}
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = q_len;
	
	bestScore = 0;
	best_read_score = 0;
	bestHits = 0;
	W1 = NWmatrices->rewards->W1;
	
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		
		/* check if index DB is loaded */
		lock(excludeDB);
		if(template >= 0 && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, index_in, template_lengths[template], kmersize, seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, index_in, template_lengths[-template], kmersize, seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseq */
		if(template < 0) {
			alnStat = KMA_score(templates_index[-template], qseq_r, q_len, qseq_r_comp, mq, scoreT, points, NWmatrices);
		} else {
			alnStat = KMA_score(templates_index[template], qseq, q_len, qseq_comp, mq, scoreT, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		start = alnStat.pos;
		end = start + aln_len - alnStat.gaps;
		if(template_lengths[abs(template)] < end) {
			end -= template_lengths[abs(template)];
		}
		
		/* penalty for non complete mapping */
		read_score = alnStat.score;
		/* full gene award */
		if((start == 0) && (end == template_lengths[abs(template)])) {
			read_score += abs(W1);
		}
		//read_score += (((start != 0) + (end != template_lengths[abs(template)])) * W1);
		
		/* Get normed score */
		if(aln_len > 0) {
			score = 1.0 * read_score / aln_len;
		} else {
			score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			if(score > bestScore) { // save as best match
				bestScore = score;
				best_read_score = read_score;
				*bestTemplates = template;
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score > best_read_score) { // save as best match
				bestScore = score;
				best_read_score = read_score;
				*bestTemplates = template;
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score == best_read_score) { // update best match
				bestTemplates[bestHits] = template;
				best_start_pos[bestHits] = start;
				best_end_pos[bestHits] = end;
				++bestHits;
			}
		}
	}
	if(best_read_score > kmersize) {
		lock(excludeOut);
		update_Scores(qseq, q_len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
	}
}

void alnFragsUnionPE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, best_read_score_r;
	int compScore, bestHits, bestHits_r, aln_len, start, end, rc, W1;
	double score;
	AlnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	W1 = NWmatrices->rewards->W1;
	start = 0;
	end = 0;
	score = 0;
	best_read_score = 0;
	best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, index_in, template_lengths[template], kmersize, seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, index_in, template_lengths[-template], kmersize, seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		template = abs(template);
		
		/* align qseqs */
		alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score <= read_score) {
			best_read_score = read_score;
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score_r <= read_score) {
			best_read_score_r = read_score;
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
		} else {
			bestTemplates_r[t_i] = 0;
			if(bestTemplates[t_i] != 0) {
				best_start_pos[t_i] = -1;
				best_end_pos[t_i] = -1;
			}
		}
		
		read_score += bestTemplates[t_i];
		if(best_start_pos[t_i] == 0 && best_end_pos[t_i] == template_lengths[abs(template)]) {
			read_score += abs(W1);
		}
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	if(best_read_score && best_read_score_r) {
		/* both matched */
		if(compScore && compScore == (best_read_score + best_read_score_r)) {
			/* proper pair */
			bestHits = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(compScore == (bestTemplates[t_i] + bestTemplates_r[t_i])) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header_r, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
				}
				lock(excludeOut);
				update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			}
		} else {
			/* unmaided pair */
			bestHits = 0;
			bestHits_r = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] == best_read_score) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				} else if(bestTemplates_r[t_i] == best_read_score_r) {
					bestTemplates_r[bestHits_r] = matched_templates[t_i];
					best_start_pos[bestHits_r] = best_start_pos[t_i];
					best_end_pos[bestHits_r] = best_end_pos[t_i];
					++bestHits_r;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
			}
			lock(excludeOut);
			update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
			update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
	} else if(best_read_score) {
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates[t_i] == best_read_score) {
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = best_start_pos[t_i];
				best_end_pos[bestHits] = best_end_pos[t_i];
				++bestHits;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
		} else if(!rc) {
			strrc(qseq, qseq_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
	} else if(best_read_score_r) {
		bestHits_r = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates_r[t_i] == best_read_score_r) {
				bestTemplates_r[bestHits_r] = matched_templates[t_i];
				best_start_pos[bestHits_r] = best_start_pos[t_i];
				best_end_pos[bestHits_r] = best_end_pos[t_i];
				++bestHits_r;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates_r < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates_r[t_i] = -bestTemplates_r[t_i];
			}
		} else if(!rc) {
			strrc(qseq_r, qseq_r_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
	}
}

void alnFragsPenaltyPE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, best_read_score_r;
	int compScore, bestHits, bestHits_r, aln_len, start, end, rc, W1, PE;
	double score;
	AlnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	W1 = NWmatrices->rewards->W1;
	PE = NWmatrices->rewards->PE;
	start = 0;
	end = 0;
	score = 0;
	best_read_score = 0;
	best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, index_in, template_lengths[template], kmersize, seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, index_in, template_lengths[-template], kmersize, seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		template = abs(template);
		
		/* align qseqs */
		alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score <= read_score) {
			best_read_score = read_score;
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		alnStat = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score_r <= read_score) {
			best_read_score_r = read_score;
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
		} else {
			bestTemplates_r[t_i] = 0;
			if(bestTemplates[t_i] != 0) {
				best_start_pos[t_i] = -1;
				best_end_pos[t_i] = -1;
			}
		}
		
		read_score += (bestTemplates[t_i] + PE);
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	if(best_read_score && best_read_score_r) {
		/* both matched */
		if(compScore && (best_read_score + best_read_score_r) <= compScore) {
			/* proper pair */
			compScore -= PE;
			bestHits = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(compScore == (bestTemplates[t_i] + bestTemplates_r[t_i])) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header_r, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
				}
				lock(excludeOut);
				update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
				unlock(excludeOut);
			}
		} else {
			/* unmaided pair */
			bestHits = 0;
			bestHits_r = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] == best_read_score) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				} else if(bestTemplates_r[t_i] == best_read_score_r) {
					bestTemplates_r[bestHits_r] = matched_templates[t_i];
					best_start_pos[bestHits_r] = best_start_pos[t_i];
					best_end_pos[bestHits_r] = best_end_pos[t_i];
					++bestHits_r;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
			}
			lock(excludeOut);
			update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
			update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
	} else if(best_read_score) {
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates[t_i] == best_read_score) {
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = best_start_pos[t_i];
				best_end_pos[bestHits] = best_end_pos[t_i];
				++bestHits;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
		} else if(!rc) {
			strrc(qseq, qseq_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
	} else if(best_read_score_r) {
		bestHits_r = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates_r[t_i] == best_read_score_r) {
				bestTemplates_r[bestHits_r] = matched_templates[t_i];
				best_start_pos[bestHits_r] = best_start_pos[t_i];
				best_end_pos[bestHits_r] = best_end_pos[t_i];
				++bestHits_r;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates_r < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates_r[t_i] = -bestTemplates_r[t_i];
			}
		} else if(!rc) {
			strrc(qseq_r, qseq_r_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
		unlock(excludeOut);
	}
}

void alnFragsForcePE(HashMap_index **templates_index, int *matched_templates, int *template_lengths, int mq, double scoreT, CompDNA *qseq_comp, CompDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, Qseqs *header, Qseqs *header_r, int kmersize, int *bestTemplates, int *bestTemplates_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, int *best_start_pos, int *best_end_pos, int seq_in, int index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, AlnPoints *points, NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, bestHits, aln_len, W1;
	int start, end, rc;
	double score, bestScore;
	AlnScore alnStat, alnStat_r;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	W1 = NWmatrices->rewards->W1;
	start = 0;
	end = 0;
	score = 0;
	bestScore = 0;
	best_read_score = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(templates_index[template], seq_in, index_in, template_lengths[template], kmersize, seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(templates_index[-template], seq_in, index_in, template_lengths[-template], kmersize, seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		template = abs(template);
		
		/* align qseq */
		alnStat = KMA_score(templates_index[template], qseq, qseq_comp->seqlen, qseq_comp, mq, scoreT, points, NWmatrices);
		if(0 < alnStat.score) {
			alnStat_r = KMA_score(templates_index[template], qseq_r, qseq_r_comp->seqlen, qseq_r_comp, mq, scoreT, points, NWmatrices);
			
			/* get read score */
			if(0 < alnStat_r.score) {
				aln_len = alnStat.len + alnStat_r.len;
				
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(alnStat.pos < alnStat_r.pos) {
					start = alnStat.pos;
					end = alnStat_r.pos + alnStat_r.len - alnStat_r.gaps;
				} else {
					start = alnStat_r.pos;
					end = alnStat.pos + alnStat.len - alnStat.gaps;
				}
				
				read_score = alnStat.score + alnStat_r.score;
				if(start == 0 && end == template_lengths[abs(template)]) {
					read_score += abs(W1);
				}
				score = 1.0 * read_score / aln_len;
			} else {
				read_score = 0;
			}
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			if(score > bestScore) { // save as best match
				bestScore = score;
				best_read_score = read_score;
				*bestTemplates = matched_templates[t_i];
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score > best_read_score) { // save as best match
				bestScore = score;
				best_read_score = read_score;
				*bestTemplates = matched_templates[t_i];
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score == best_read_score) { // update best match
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = start;
				best_end_pos[bestHits] = end;
				++bestHits;
			}
		}
	}
	
	if(best_read_score) {
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
			lock(excludeOut);
			update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header_r, header, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		} else {
			if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
			}
			lock(excludeOut);
			update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, alignment_scores, uniq_alignment_scores, frag_out_raw);
			unlock(excludeOut);
		}
	}
}

void * alnFrags_threaded(void * arg) {
	
	static volatile int excludeIn[1] = {0}, excludeOut[1] = {0}, excludeDB[1] = {0};
	Aln_thread *thread = arg;
	int rc_flag, read_score, delta, index_in, seq_in, kmersize, mq;
	int *matched_templates, *bestTemplates, *bestTemplates_r;
	int *best_start_pos, *best_end_pos, *template_lengths;
	long *index_indexes, *seq_indexes;
	long unsigned *alignment_scores, *uniq_alignment_scores;
	double scoreT;
	FILE *inputfile, *frag_out_raw;
	CompDNA *qseq_comp, *qseq_r_comp;
	Qseqs *qseq, *qseq_r, *header, *header_r;
	AlnPoints *points;
	NWmat *NWmatrices;
	HashMap_index **templates_index;
	
	/* get input */
	matched_templates = thread->matched_templates;
	bestTemplates = thread->bestTemplates;
	bestTemplates_r = thread->bestTemplates_r;
	best_start_pos = thread->best_start_pos;
	best_end_pos = thread->best_end_pos;
	alignment_scores = thread->alignment_scores;
	uniq_alignment_scores = thread->uniq_alignment_scores;
	index_indexes = thread->index_indexes;
	seq_indexes = thread->seq_indexes;
	inputfile = thread->inputfile;
	frag_out_raw = thread->frag_out_raw;
	index_in = thread->index_in;
	seq_in = thread->seq_in;
	qseq_comp = thread->qseq_comp;
	qseq_r_comp = thread->qseq_r_comp;
	qseq = thread->qseq;
	qseq_r = thread->qseq_r;
	header = thread->header;
	header_r = thread->header_r;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	kmersize = thread->kmersize;
	mq = thread->mq;
	scoreT = thread->scoreT;
	template_lengths = thread->template_lengths;
	templates_index = thread->templates_index;
	
	delta = qseq->size;
	read_score = 0;
	//lock(excludeIn);
	lockTime(excludeIn, 65536);
	while((rc_flag = get_ankers(matched_templates, qseq_comp, header, inputfile)) != 0) {
		if(*matched_templates) { // SE
			read_score = 0;
		} else { // PE
			read_score = get_ankers(matched_templates, qseq_r_comp, header_r, inputfile);
			read_score = labs(read_score);
			qseq_r->len = qseq_r_comp->seqlen;
		}
		unlock(excludeIn);
		qseq->len = qseq_comp->seqlen;
		
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
		
		if(kmersize <= qseq->len) {
			if(read_score && kmersize <= qseq_r->len) { // PE
				alnFragsPE(templates_index, matched_templates, template_lengths, mq, scoreT, qseq_comp, qseq_r_comp, qseq->seq, qseq_r->seq, header, header_r, kmersize, bestTemplates, bestTemplates_r, alignment_scores, uniq_alignment_scores, best_start_pos, best_end_pos, seq_in, index_in, seq_indexes, index_indexes, frag_out_raw, points, NWmatrices, excludeOut, excludeDB);
			} else { // SE
				alnFragsSE(templates_index, matched_templates, template_lengths, mq, scoreT, rc_flag, qseq_comp, qseq_r_comp, qseq->seq, qseq_r->seq, qseq->len, kmersize, header, bestTemplates, alignment_scores, uniq_alignment_scores, best_start_pos, best_end_pos, seq_in, index_in, seq_indexes, index_indexes, frag_out_raw, points, NWmatrices, excludeOut, excludeDB);
			}
		}
		lock(excludeIn);
	}
	unlock(excludeIn);
	
	return NULL;
}
