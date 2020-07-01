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
#include "penalties.h"

#ifndef KMERANKER
#define KMERANKER 1;
typedef struct kmerAnker KmerAnker;
struct kmerAnker {
	int score;
	int weight;
	unsigned start;
	unsigned end;
	unsigned ties;
	unsigned *values;
	struct kmerAnker *descend; /* descending anker */
};
#endif

extern KmerAnker * (*getChainTemplates)(KmerAnker*, const Penalties*, int, int*, int*, int*, char*);
KmerAnker * getBestChainTemplates(KmerAnker *src, const Penalties *rewards, int kmersize, int *bests, int *Score, int *extendScore, char *include);
KmerAnker * getProxiChainTemplates(KmerAnker *src, const Penalties *rewards, int kmersize, int *bests, int *Score, int *extendScore, char *include);
KmerAnker * pruneAnkers(KmerAnker *V_score, int kmersize);
KmerAnker * getBestAnker(KmerAnker **src, unsigned *ties);
KmerAnker * getTieAnker(int stop, KmerAnker *src, int score);
