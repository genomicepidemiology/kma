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
#include <stdlib.h>
#include "kmeranker.h"
#include "penalties.h"

KmerAnker * getBestChainTemplates(KmerAnker *src, const Penalties *rewards, int *bests, int *Score, int *extendScore, char *include) {
	
	/* get set of best templates and silences the chain, except for the initial anker */
	int template, score, bestScore, Wl, W1, U, M, MM, gaps;
	unsigned i, j, start, end, SU;
	unsigned *values;
	short unsigned *values_s;
	KmerAnker *node, *prev;
	
	/* mark template candidates */
	if(!src) {
		return 0;
	} else if((SU = *include)) {
		values_s = (short unsigned *) src->values;
		i = *values_s + 1;
		*bests += i;
		values_s += i;
		bests += i;
		while(--i) {
			include[(*--bests = *--values_s)] ^= 1;
		}
	} else {
		values = src->values;
		i = *values + 1;
		*bests += i;
		values += i;
		bests += i;
		while(--i) {
			include[(*--bests = *--values)] ^= 1;
		}
	}
	M = rewards->M;
	MM = rewards->MM;
	U = rewards->U;
	W1 = rewards->W1;
	Wl = rewards->Wl;
	bestScore = src->score;
	
	/* get chaining scores */
	for(node = src; node != 0; node = node->next) {
		if(SU) {
			values_s = (short unsigned *) node->values;
			i = *values_s + 1;
			values_s += i;
		} else {
			values = node->values;
			i = *values + 1;
			values += i;
		}
		
		start = node->start;
		end = node->end;
		
		while(--i) {
			template = SU ? *--values_s : *--values;
			if(!include[template]) {
				score = Score[template];
				gaps = start - extendScore[template];
				if(gaps == 0) {
					/* continued match */
					score += M;
				} else if(gaps < 0) {
					/* deletion */
					score += (W1 + abs(gaps + 1) * U);
				} else if(gaps == 1) {
					/* one unmatched between */
					score += MM;
				} else {
					/* unmatched bases between */
					if((gaps = (gaps - 2) * M) < 0) {
						score += gaps;
					}
					score += (MM << 1);
				}
				
				/* check local chain is needed */
				if(score < Wl) {
					score = Wl + node->weight;
					/* mark as terminating */
					extendScore[template] = 0;
				} else {
					score += node->weight;
					/* mark descendant */
					extendScore[template] = end;
				}
				
				/* update Scores */
				Score[template] = score;
			}
		}
		
		/* mark as used */
		node->score = 0;
		prev = node;
	}
	src->score = bestScore;
	
	/* get best templates */
	j = 0;
	for(i = 1; i <= *bests; ++i) {
		if(Score[(template = bests[i])] == bestScore) {
			bests[++j] = template;
		}
		/* clear Score arrays */
		Score[template] = 0;
		extendScore[template] = 0;
		include[template] = 0;
	}
	*bests = j;
	
	return prev;
}

KmerAnker * pruneAnkers(KmerAnker *V_score, int kmersize) {
	
	KmerAnker *node, *prev;
	
	while(V_score->score < kmersize) {
		++V_score;
	}
	
	prev = V_score;
	node = V_score->descend;
	while(node) {
		if(kmersize <= node->score) {
			prev->descend = node;
			prev = node;
		}
		node = node->descend;
	}
	prev->descend = 0;
	
	return V_score;
}

unsigned getStartAnker(KmerAnker *src) {
	
	while(src->next) {
		src = src->next;
	}
	return src->start;
}

KmerAnker * getBestAnker(KmerAnker **src, unsigned *ties) {
	
	KmerAnker *node, *prev, *best;
	
	*ties = 0;
	prev = *src;
	while(prev && prev->score == 0) {
		prev = prev->descend;
	}
	*src = prev;
	best = prev;
	node = prev->descend;
	
	while(node) {
		if(best->score < node->score) {
			best = node;
			*ties = 0;
		} else if(best->score == node->score) {
			if((node->end - node->start) < (best->end - best->start)) {
				best = node;
				*ties = 0;
			} else if((node->end - node->start) == (best->end - best->start)) {
				best = node;
				++*ties;
			}
		} else if(node->score != 0) {
			prev->descend = node;
			prev = node;
		}
		node = node->descend;
	}
	prev->descend = 0;
	
	return best;
}

KmerAnker * getTieAnker(KmerAnker *stop, KmerAnker *src, int score, int len) {
	
	/* search downwards */
	while(--src != stop) {
		if(src->score == score && (src->end - src->start) == len) {
			return src;
		}
	}
	
	return 0;
}
