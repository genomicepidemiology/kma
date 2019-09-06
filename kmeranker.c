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
#include <math.h>
#include "kmeranker.h"
#include "penalties.h"

KmerAnker * getBestChainTemplates(KmerAnker *src, const Penalties *rewards, int kmersize, int *bests, int *Score, int *extendScore, char *include) {
	
	/* get set of best templates and silences the chain, except for the initial anker */
	int i, j, template, score, bestScore, Wl, W1, U, M, MM, gaps;
	int start, end, pos, SU;
	unsigned *values;
	short unsigned *values_s;
	KmerAnker *node, *prev;
	
	/* mark template candidates */
	if(!src) {
		return 0;
	} else if((SU = *include)) {
		values_s = (short unsigned *) src->values;
		i = *values_s;
		*bests = i;
		++i;
		values_s += i;
		bests += i;
		while(--i) {
			include[(*--bests = *--values_s)] ^= 1;
		}
	} else {
		values = src->values;
		i = *values;
		*bests = i;
		++i;
		values += i;
		bests += i;
		while(--i) {
			include[(*--bests = *--values)] ^= 1;
		}
	}
	--bests;
	M = rewards->M;
	MM = rewards->MM;
	U = rewards->U;
	W1 = rewards->W1;
	Wl = rewards->Wl;
	bestScore = src->score;
	values_s = 0;
	values = 0;
	
	/* get chaining scores */
	for(node = src; node; node = node->next) {
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
			if(include[template]) {
				score = Score[template];
				pos = extendScore[template];
				gaps = pos - end;
				
				/* extend chain */
				if(pos == 0) {
					if(node->cStart) {
						score = W1 + (node->cStart - 1) * U;
						score = node->weight + (Wl < score ? score : Wl);
					} else {
						score = node->weight;
					}
				} else {
					if(gaps < 0) {
						if(gaps == -kmersize) {
							score += node->weight - (kmersize - 1) * M;
						} else {
							score += (W1 + (-gaps - 1) * U) + node->weight + gaps * M;
						}
					} else if(gaps == 0) {
						score += node->weight + W1;
					} else if(gaps <= 2) {
						score += node->weight + gaps * MM;
					} else if((MM * 2 + (gaps - 2) * M) < 0) {
						score += node->weight + (MM * 2 + (gaps - 2) * M);
					} else {
						score += node->weight - log(gaps * M);
					}
					
					/* verify extension */
					if(score < 0 || score == bestScore) {
						include[template] = 0;
					}
				}
				extendScore[template] = start;
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

KmerAnker * pruneAnkers(KmerAnker *V_score, int kmersize, double mrs) {
	
	KmerAnker *node, *prev;
	
	if(!V_score || !V_score->score) {
		return 0;
	}
	while((V_score->score < kmersize || V_score->score < mrs * (V_score->end - V_score->cStart)) && (V_score = V_score->descend));
	
	if(!V_score) {
		return 0;
	}
	prev = V_score;
	node = V_score->descend;
	while(node) {
		if(kmersize <= node->score && mrs * (node->end - node->cStart) <= node->score) {
			prev->descend = node;
			prev = node;
		}
		node = node->descend;
	}
	prev->descend = 0;
	
	return V_score;
}

KmerAnker * getBestAnker(KmerAnker **src, unsigned *ties) {
	
	KmerAnker *node, *prev, *best;
	
	*ties = 0;
	prev = *src;
	while(prev && prev->score == 0) {
		prev = prev->descend;
	}
	*src = prev;
	if(!(best = prev)) {
		return 0;
	}
	node = prev->descend;
	while(node) {
		if(best->score < node->score) {
			best = node;
			*ties = 0;
		} else if(best->score == node->score) {
			if((node->end - node->start) < (best->end - best->cStart)) {
				best = node;
				*ties = 0;
			} else if((node->end - node->start) == (best->end - best->cStart)) {
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

KmerAnker * getTieAnker(int stop, KmerAnker *src, int score, int len) {
	
	if(!src || src->start == stop) {
		return 0;
	}
	
	/* search downwards */
	while((--src)->start != stop) {
		if(src->score == score && (src->end - src->cStart) == len) {
			return src;
		}
	}
	
	return 0;
}
