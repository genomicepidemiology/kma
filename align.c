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
#include <string.h>
#include "align.h"
#include "chain.h"
#include "compdna.h"
#include "hashmapindex.h"
#include "nw.h"
#include "stdnuc.h"
#include "stdstat.h"

AlnScore KMA(const HashMap_index *template_index, const unsigned char *qseq, int q_len, Aln *aligned, Aln *Frag_align, int min, int max, int mq, double scoreT, AlnPoints *points, NWmat *matrices) {
	
	int i, j, k, bias, prev, start, stop, t_len, value, end, band, mem_count;
	int t_l, t_s, t_e, q_s, q_e, score, shifter, kmersize, U, M, **d;
	long unsigned key, mask;
	unsigned char nuc;
	AlnScore Stat, NWstat;
	Penalties *rewards;
	
	/* Extract indexes and template sequence */
	rewards = points->rewards;
	U = rewards->U;
	M = rewards->M;
	d = rewards->d;
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	key = 0;
	/* circular, skip boundaries */
	if(min < max) {
		min = 0;
		max = t_len;
	}
	/* find seeds */
	if(points->len) {
		mem_count = points->len;
	} else {
		mem_count = 0;
		i = 0;
		while(i < q_len) {
			end = charpos(qseq, 4, i, q_len);
			if(end == -1) {
				end = q_len;
			}
			
			if(i < end - kmersize) {
				key = makeKmer(qseq, i, kmersize - 1);
				i += (kmersize - 1);
			} else {
				i = end + 1;
			}
			
			while(i < end) {
				key = ((key << 2) | qseq[i]) & mask;
				value = hashMap_index_get_bound(template_index, key, min, max, shifter);
				
				if(value == 0) {
					++i;
				} else if(0 < value) {
					i -= (kmersize - 1);
					
					/* backseed for overlapping seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
					}
					
					/* get start positions */
					points->qStart[mem_count] = j + 1;
					points->tStart[mem_count] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
					}
					
					/* get end positions */
					points->qEnd[mem_count] = i;
					points->tEnd[mem_count] = value + 1;
					
					/* calculate weight */
					points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
					++mem_count;
					
					/* realloc seeding points */
					if(mem_count == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				} else {
					/* move counter back */
					i -= (kmersize - 1);
					
					/* get position in hashmap */
					stop = hashMap_index_getDubPos(template_index, key, value, shifter);
					
					/* get all mems */
					bias = i;
					while(0 <= stop) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(template_index->index[stop]);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[mem_count] = j + 1;
						points->tStart[mem_count] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[mem_count] = k;
						points->tEnd[mem_count] = value + 1;
						
						/* calculate weight */
						points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
						++mem_count;
						
						/* realloc seeding points */
						if(mem_count == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						stop = hashMap_index_getNextDubPos(template_index, key, min, max, stop, shifter);
					}
					i = bias + 1;
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
					
				}
			}
			i = end + 1;
		}
	}
	
	if(mem_count) {
		points->len = mem_count;
	} else {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		aligned->s[0] = 0;
		aligned->len = 0;
		points->len = 0;
		return Stat;
	}
	
	/* get best seed chain, returns best starting point */
	start = chainSeedsPtr(points, q_len, t_len, kmersize, &aligned->mapQ);
	score = points->score[start];
	
	//fprintf(stderr, "%d\t%d\t%d\t%d\n", t_len, mem_count, aligned->mapQ, score);
	if(aligned->mapQ < mq || score < kmersize || (score << 1) < scoreT * (points->qEnd[points->len - 1] - points->qStart[0])) {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		aligned->s[0] = 0;
		aligned->len = 0;
		points->len = 0;
		return Stat;
	}
	
	/* initialize */
	Stat.len = 0;
	Stat.score = 0;
	Stat.gaps = 0;
	value = points->tStart[start] - 1;
	Stat.pos = value;
	i = points->qStart[start];
	
	/* align leading tail */
	if(i != 0) {
		/* get boundaries */
		t_s = 0;
		t_e = value;
		q_s = 0;
		q_e = i;
		if((q_e << 1) < t_e || (q_e + 64) < t_e) { // big leading template gap, cut down
			t_s = t_e - MIN(64, (q_e << 1));
		} else if((t_e << 1) < q_e || (t_e + 64) < q_e) { // big leading query gap, cut down
			q_s = q_e - MIN(64, (t_e << 1));
		}
		
		/* align */
		if(t_e - t_s > 0 && q_e - q_s > 0) {
			band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
			if(q_e - q_s <= band || t_e - t_s <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
				NWstat = NW(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, matrices);
			} else {
				NWstat = NW_band(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, band, matrices);
			}
			
			/* trim leading gaps */
			bias = 0;
			if(t_s == 0) {
				while(bias < NWstat.len && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
					if(Frag_align->t[bias] == 5) {
						--NWstat.gaps;
					}
					++bias;
				}
				NWstat.len -= bias;
				/*if(bias) {
					NWstat.score -= (W1 + (bias - 1) * U);
				}*/
			}
			
			memcpy(aligned->t, Frag_align->t + bias, NWstat.len);
			memcpy(aligned->s, Frag_align->s + bias, NWstat.len);
			memcpy(aligned->q, Frag_align->q + bias, NWstat.len);
			Stat.pos -= (NWstat.len - NWstat.gaps);
			Stat.score = NWstat.score;
			Stat.len = NWstat.len;
			Stat.gaps = NWstat.gaps;
		}
	}
	
	/* piece seeds together */
	stop = 1;
	while(stop) {
		/* MEM */
		q_s = points->qStart[start];
		end = points->qEnd[start] - q_s;
		memcpy(aligned->t + Stat.len, qseq + q_s, end);
		memset(aligned->s + Stat.len, '|', end);
		memcpy(aligned->q + Stat.len, qseq + q_s, end);
		Stat.len += end;
		end = points->qEnd[start];
		for(i = points->qStart[start]; i < end; ++i) {
			nuc = qseq[i];
			Stat.score += d[nuc][nuc];
		}
		
		/* join MEMs */
		if(points->next[start]) {
			/* get positions between seed-extends */
			q_s = points->qEnd[start];
			t_s = points->tEnd[start] - 1;
			start = points->next[start];
			
			/* check if next MEM is a semi match */
			if(points->qStart[start] < q_s) {
				points->tStart[start] += (q_s - points->qStart[start]);
				points->qStart[start] = q_s;
			}
			t_e = points->tStart[start] - 1;
			
			if(t_e < t_s) {
				if(t_s <= points->tEnd[start]) {
					points->qStart[start] += (t_s - t_e);
					t_e = t_s;
					t_l = t_e - t_s;
				} else {
					/* circular joining */
					Frag_align->pos = t_len;
					t_l = t_len - t_s + t_e;
				}
			} else {
				t_l = t_e - t_s;
			}
			q_e = points->qStart[start];
			
			/* piece seed-extends together */
			if(abs(t_l - q_e + q_s) * U > q_len * M || t_l > q_len || q_e - q_s > (q_len >> 1)) {
				/* gap is too big to give a positive score */
				Stat.score = 0;
				Stat.len = 1;
				Stat.gaps = 0;
				aligned->s[0] = 0;
				aligned->len = 0;
				points->len = 0;
				return Stat;
			}
			if((t_l > 0 || q_e - q_s > 0)) {
				band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
				if(q_e - q_s <= band || t_l <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
					NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, matrices);
				} else {
					NWstat = NW_band(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, band, matrices);
					//NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align);
				}
				memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
				memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
				memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
				Stat.score += NWstat.score;
				Stat.len += NWstat.len;
				Stat.gaps += NWstat.gaps;
			}
		} else {
			stop = 0;
		}
	}
	
	/* align trailing tail */
	/* Get intervals in query and template to align */
	q_s = points->qEnd[start];
	t_s = points->tEnd[start] - 1;
	q_e = q_len;
	t_e = t_len;
	if(((q_len - q_s) << 1) < (t_len - t_s) || (q_len - q_s + 64) < (t_len - t_s)) { // big trailing template gap, cut down
		t_e = t_s + MIN(64, ((q_len - q_s) << 1));
	} else if(((t_len - t_s) << 1) < (q_len - q_s) || (t_len - t_s + 64) < (q_len - q_s)) { // big leading query gap, cut down
		q_e = q_s + MIN(64, ((t_len - t_s) << 1));
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
		if(q_e - q_s <= band || t_e - t_s <= band) {//|| abs(t_e - t_s - q_e - q_s) >= 32) {
			NWstat = NW(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, matrices);
		} else {
			NWstat = NW_band(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, band, matrices);
			//NWstat = NW(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align);
		}
		/* trim trailing gaps */
		if(t_e == t_len) {
			bias = NWstat.len - 1;
			while(bias && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
				if(Frag_align->t[bias] == 5) {
					--NWstat.gaps;
				}
				--bias;
			}
			++bias;
			/*
			if(bias != NWstat.len) {
				NWstat.score -= (W1 + (NWstat.len - bias) * U);
				NWstat.len = bias;
			}
			*/
		}
		
		memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
		memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
		memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
		Stat.score += NWstat.score;
		Stat.len += NWstat.len;
		Stat.gaps += NWstat.gaps;
	}
	aligned->s[Stat.len] = 0;
	aligned->len = Stat.len;
	points->len = 0;
	
	return Stat;
}

AlnScore KMA_score(const HashMap_index *template_index, const unsigned char *qseq, int q_len, const CompDNA *qseq_comp, int mq, double scoreT, AlnPoints *points, NWmat *matrices) {
	
	int i, j, k, l, bias, prev, start, stop, t_len, value, end, band;
	int t_l, t_s, t_e, q_s, q_e, mem_count, score, kmersize;
	int U, M, **d;
	unsigned mapQ, shifter;
	long unsigned key;
	unsigned char nuc;
	AlnScore Stat, NWstat;
	Penalties *rewards;
	
	/* Extract indexes and template sequence */
	rewards = points->rewards;
	U = rewards->U;
	M = rewards->M;
	d = rewards->d;
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (template_index->kmerindex << 1);
	
	/* find seeds */
	mem_count = 0;
	for(i = 1, j = 0; i <= qseq_comp->N[0]; ++i) {
		end = qseq_comp->N[i] - kmersize + 1;
		while(j < end) {
			value = hashMap_index_get(template_index, getKmer(qseq_comp->seq, j, shifter), shifter);
			
			if(value == 0) {
				++j;
			} else if(0 < value) {
				/* backseed for ambiguos seeds */
				prev = value - 2;
				for(k = j - 1; 0 <= k && 0 <= prev && qseq[k] == getNuc(template_index->seq, prev); --k) {
					--prev;
				}
				
				/* get start positions */
				points->qStart[mem_count] = k + 1;
				points->tStart[mem_count] = prev + 2;
				
				
				/* skip k-mer bases */
				value += (kmersize - 1);
				j += kmersize;
				
				/* extend */
				end += (kmersize - 1);
				while(j < end && value < t_len && qseq[j] == getNuc(template_index->seq, value)) {
					++j;
					++value;
				}
				end -= (kmersize - 1);
				
				/* get end positions */
				points->qEnd[mem_count] = j;
				points->tEnd[mem_count] = value + 1;
				
				/* calculate weight */
				points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
				++mem_count;
				
				/* realloc seeding points */
				if(mem_count == points->size) {
					seedPoint_realloc(points, points->size << 1);
				}
			} else {
				/* get position in hashmap */
				key = getKmer(qseq_comp->seq, j, shifter);
				stop = hashMap_index_getDubPos(template_index, key, value, shifter);
				
				/* get all mems */
				bias = j;
				while(0 <= stop) {
					/* get mem info */
					l = j;
					/* backseed for overlapping seeds */
					value = abs(template_index->index[stop]);
					prev = value - 2;
					for(k = l - 1; 0 <= k && 0 <= prev && qseq[k] == getNuc(template_index->seq, prev); --k) {
						--prev;
					}
					
					/* get start positions */
					points->qStart[mem_count] = k + 1;
					points->tStart[mem_count] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					l += kmersize;
					
					/* extend */
					end += (kmersize - 1);
					while(l < end && value < t_len && qseq[l] == getNuc(template_index->seq, value)) {
						++l;
						++value;
					}
					end -= (kmersize - 1);
					
					/* get end positions */
					points->qEnd[mem_count] = l;
					points->tEnd[mem_count] = value + 1;
					
					/* calculate weight */
					points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
					++mem_count;
					
					/* realloc seeding points */
					if(mem_count == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					if(bias < l) {
						bias = l;
					}
					
					stop = hashMap_index_getNextDubPos(template_index, key, 0, t_len, stop, shifter);
				}
				j = bias + 1;
			}
		}
		j = qseq_comp->N[i] + 1;
	}
	
	if(mem_count) {
		points->len = mem_count;
	} else {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		return Stat;
	}
	
	/* get best seed chain, returns best starting point */
	start = chainSeedsPtr(points, q_len, t_len, kmersize, &mapQ);
	score = points->score[start];
	
	//if(score < (q_len >> 5)) {
	if(mapQ < mq || score < kmersize || (score << 1) < scoreT * (points->qEnd[points->len - 1] - points->qStart[0])) {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		return Stat;
	}
	
	/* initialize */
	Stat.len = 0;
	Stat.score = 0;
	Stat.gaps = 0;
	value = points->tStart[start] - 1;
	Stat.pos = value;
	i = points->qStart[start];
	
	/* align leading tail */
	if(i != 0) {
		/* get boundaries */
		t_s = 0;
		t_e = value;
		q_s = 0;
		q_e = i;
		if((q_e << 1) < t_e || (q_e + 64) < t_e) { // big leading template gap, cut down
			t_s = t_e - MIN(64, (q_e << 1));
		} else if((t_e << 1) < q_e || (t_e + 64) < q_e) { // big leading query gap, cut down
			q_s = q_e - MIN(64, (t_e << 1));
		}
		
		/* align */
		if(t_e - t_s > 0 && q_e - q_s > 0) {
			band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
			if(q_e - q_s <= band || t_e - t_s <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
				NWstat = NW_score(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, matrices, t_len);
			} else {
				NWstat = NW_band_score(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, band, matrices, t_len);
			}
			Stat.pos -= (NWstat.len - NWstat.gaps);
			Stat.score = NWstat.score;
			Stat.len = NWstat.len;
			Stat.gaps = NWstat.gaps;
		}
	}
	
	/* piece seeds together */
	stop = 1;
	while(stop) {
		/* MEM */
		q_s = points->qStart[start];
		end = points->qEnd[start] - q_s;
		Stat.len += end;
		end = points->qEnd[start];
		for(i = points->qStart[start]; i < end; ++i) {
			nuc = qseq[i];
			Stat.score += d[nuc][nuc];
		}
		
		/* join MEMs */
		if(points->next[start]) {
			/* get positions between seed-extends */
			q_s = points->qEnd[start];
			t_s = points->tEnd[start] - 1;
			start = points->next[start];
			
			/* check if next MEM is a semi match, or a circular joining */
			if(points->qStart[start] < q_s) {
				points->tStart[start] += (q_s - points->qStart[start]);
				points->qStart[start] = q_s;
			}
			t_e = points->tStart[start] - 1;
			
			if(t_e < t_s) {
				if(t_s <= points->tEnd[start]) {
					points->qStart[start] += (t_s - t_e);
					t_e = t_s;
					t_l = t_e - t_s;
				} else {
					/* circular joining */
					t_l = t_len - t_s + t_e;
				}
			} else {
				t_l = t_e - t_s;
			}
			q_e = points->qStart[start];
			
			/* piece seed-extends together */
			if(abs(t_l - q_e + q_s) * U > q_len * M || t_l > q_len || q_e - q_s > (q_len >> 1)) {
				/* gap is too big to give a positive score */
				Stat.score = 0;
				Stat.len = 1;
				Stat.gaps = 0;
				return Stat;
			}
			if((t_l > 0 || q_e - q_s > 0)) {
				band = 4 * abs(t_l - q_e + q_s) + 64;
				if(q_e - q_s <= band || t_l <= band) {
					NWstat = NW_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, matrices, t_len);
				} else {
					NWstat = NW_band_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, band, matrices, t_len);
				}
				Stat.score += NWstat.score;
				Stat.len += NWstat.len;
				Stat.gaps += NWstat.gaps;
			}
		} else {
			stop = 0;
		}
	}
	
	/* align trailing tail */
	/* Get intervals in query and template to align */
	q_s = points->qEnd[start];
	t_s = points->tEnd[start] - 1;
	q_e = q_len;
	t_e = t_len;
	if((t_len - t_s) > (q_len - q_s + 64) || (t_len - t_s) > ((q_len - q_s) << 1)) { // big trailing template gap, cut down
		t_e = t_s + MIN(64, ((q_len - q_s) << 1));
	} else if ((q_len - q_s) > (t_len - t_s + 64) || (q_len - q_s) > ((t_len - t_s) << 1)) { // big leading query gap, cut down
		q_e = q_s + MIN(64, ((t_len - t_s) << 1));
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
		if(q_e - q_s <= band || t_e - t_s <= band) {//|| abs(t_e - t_s - q_e - q_s) >= 32) {
			NWstat = NW_score(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, matrices, t_len);
		} else {
			NWstat = NW_band_score(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, band, matrices, t_len);
		}
		Stat.score += NWstat.score;
		Stat.len += NWstat.len;
		Stat.gaps += NWstat.gaps;
	}
	
	return Stat;
}

int preseed(const HashMap_index *template_index, unsigned char *qseq, int q_len) {
	
	static int exhaustive = 1;
	int i, shifter, kmersize;
	
	if(exhaustive) {
		exhaustive = q_len;
		return 0;
	}
	
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	for(i = 0; i < q_len; i += kmersize) {
		if(hashMap_index_get_bound(template_index, makeKmer(qseq, i, kmersize), 0, template_index->len, shifter) != 0) {
			return 0;
		}
	}
	
	return i;
}

void intcpy(int *dest, int *src, int size) {
	
	while(size--) {
		*dest++ = *src++;
	}
}

int anker_rc(const HashMap_index *template_index, unsigned char *qseq, int q_len, AlnPoints *points) {
	
	int i, j, k, rc, end, stop, score, score_r, value, t_len, prev, bias;
	int bestScore, mem_count, totMems, shifter, kmersize;
	long unsigned key, mask;
	
	t_len = template_index->len;
	kmersize = template_index->kmerindex;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	key = 0;
	
	/* find seeds */
	bestScore = 0;
	score = 0;
	score_r = 0;
	mem_count = 0;
	totMems = 0;
	points->len = 0;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			strrc(qseq, q_len);
			score = score_r;
			points->len = mem_count;
		}
		score_r = 0;
		mem_count = 0;
		i = preseed(template_index, qseq, q_len);
		while(i < q_len) {
			end = charpos(qseq, 4, i, q_len);
			if(end == -1) {
				end = q_len;
			}
			if(i < end - kmersize) {
				key = makeKmer(qseq, i, kmersize - 1);
				i += (kmersize - 1);
			} else {
				i = end + 1;
			}
			
			while(i < end) {
				key = ((key << 2) | qseq[i]) & mask;
				value = hashMap_index_get_bound(template_index, key, 0, t_len, shifter);
				if(value == 0) {
					++i;
				} else if(0 < value) {
					i -= (kmersize - 1);
					
					/* backseed for ambiguos seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
						++score_r;
					}
					
					/* get start positions */
					points->qStart[totMems] = j + 1;
					points->tStart[totMems] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					score_r += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
						++score_r;
					}
					
					/* get end positions */
					points->qEnd[totMems] = i;
					points->tEnd[totMems] = value + 1;
					
					/* calculate weight */
					points->weight[totMems] = (points->tEnd[totMems] - points->tStart[totMems]);
					++mem_count;
					++totMems;
					
					/* realloc seeding points */
					if(totMems == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				} else {
					/* move counter back */
					i -= (kmersize - 1);
					score_r += kmersize;
					
					/* get position in hashmap */
					stop = hashMap_index_getDubPos(template_index, key, value, shifter);
					
					/* get all mems */
					bias = i;
					while(0 <= stop) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(template_index->index[stop]);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[totMems] = j + 1;
						points->tStart[totMems] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[totMems] = k;
						points->tEnd[totMems] = value + 1;
						
						/* calculate weight */
						points->weight[totMems] = (points->qEnd[totMems] - points->qStart[totMems]);
						++mem_count;
						++totMems;
						
						/* realloc seeding points */
						if(totMems == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						
						stop = hashMap_index_getNextDubPos(template_index, key, 0, t_len, stop, shifter);
					}
					i = bias + 1;
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				}
			}
			i = end + 1;
		}
		if(bestScore < score_r) {
			bestScore = score_r;
		}
	}
	
	if(bestScore * kmersize < (q_len - kmersize - bestScore)) {
		bestScore = 0;
		points->len = 0;
	} else if(bestScore == score) {
		strrc(qseq, q_len);
	} else {
		/* move mems down */
		if(points->len) {
			intcpy(points->tStart, points->tStart + points->len, mem_count);
			intcpy(points->tEnd, points->tEnd + points->len, mem_count);
			intcpy(points->qStart, points->qStart + points->len, mem_count);
			intcpy(points->qEnd, points->qEnd + points->len, mem_count);
			intcpy(points->weight, points->weight + points->len, mem_count);
		}
		points->len = mem_count;
	}
	
	return bestScore;
}
