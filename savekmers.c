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
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "ankers.h"
#include "compdna.h"
#include "hashmapkma.h"
#include "penalties.h"
#include "pherror.h"
#include "qseqs.h"
#include "savekmers.h"
#include "stdnuc.h"
#include "stdstat.h"
#include "threader.h"

int loadFsa(CompDNA *qseq, Qseqs *header, FILE *inputfile) {
	
	int buffer[4];
	
	if(fread(buffer, sizeof(int), 4, inputfile)) {
		qseq->seqlen = buffer[0];
		qseq->complen = buffer[1];
		/* if pair, header->len < 0 */
		header->len = abs(buffer[3]);
		
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
		}
		qseq->N[0] = buffer[2];
		
		if(header->size <= header->len) {
			header->size = header->len << 1;
			free(header->seq);
			header->seq = smalloc(header->size);
		}
		fread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		fread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		fread(header->seq, 1, header->len, inputfile);
	} else {
		qseq->seqlen = 0;
		return 0;
	}
	
	return buffer[3];
}

void * save_kmers_threaded(void *arg) {
	
	static volatile int excludeIn[1] = {0}, excludeOut[1] = {0};
	static unsigned readNum = 0;
	KmerScan_thread *thread = arg;
	int *Score, *Score_r, *bestTemplates, *bestTemplates_r, *regionTemplates;
	int *regionScores, *extendScore, go, spltDB, exhaustive;
	FILE *inputfile;
	HashMapKMA *templates;
	CompDNA *qseq, *qseq_r;
	Qseqs *header, *header_r;
	Penalties *rewards;
	
	templates = thread->templates;
	exhaustive = thread->exhaustive;
	qseq = thread->qseq;
	qseq_r = thread->qseq_r;
	header = thread->header;
	bestTemplates = thread->bestTemplates;
	bestTemplates_r = thread->bestTemplates_r;
	rewards = thread->rewards;
	if(save_kmers_pair != &save_kmers_unionPair) {
		regionScores = calloc(templates->DB_size << 1, sizeof(int));
		if(!regionScores) {
			ERROR();
		}
	} else {
		regionScores = 0;
	}
	extendScore = calloc((templates->DB_size + 1), sizeof(int));
	header_r = malloc(sizeof(Qseqs));
	if(!extendScore || !header_r) {
		ERROR();
	}
	header_r = setQseqs(256);
	regionTemplates = smalloc(((templates->DB_size << 1) + 4) * sizeof(int));
	Score = calloc(templates->DB_size, sizeof(int));
	Score_r = calloc(templates->DB_size, sizeof(int));
	if(!Score || !Score_r) {
		ERROR();
	}
	inputfile = thread->inputfile;
	*Score = thread->num;
	*bestTemplates++ = templates->DB_size;
	*bestTemplates_r++ = templates->DB_size;
	*regionTemplates++ = templates->DB_size;
	*bestTemplates++ = thread->num;
	*bestTemplates_r++ = thread->num;
	*regionTemplates++ = thread->num;
	*bestTemplates++ = 0;
	*bestTemplates_r++ = 0;
	*regionTemplates++ = 0;
	spltDB = thread->spltDB;
	
	go = 1;
	while(go != 0) {
		
		/* load qseqs */
		lock(excludeIn);
		if((go = loadFsa(qseq, header, inputfile)) < 0) {
			/* PE */
			loadFsa(qseq_r, header_r, inputfile);
		}
		if(spltDB) {
			bestTemplates[-1] = ++readNum;
			bestTemplates_r[-1] = readNum;
			regionTemplates[-1] = readNum;
		}
		unlock(excludeIn);
		
		/* allocate memory */
		if(qseq_r->size < qseq->size && 0 < go) {
			freeComp(qseq_r);
			allocComp(qseq_r, qseq->size);
		}
		
		/* find ankers */
		if(0 < go) {
			kmerScan(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, qseq_r, header, extendScore, exhaustive, excludeOut);
		} else if(go < 0) {
			save_kmers_pair(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores, qseq, qseq_r, header, header_r, extendScore, exhaustive, excludeOut);
		}
	}
	
	return NULL;
}

int get_kmers_for_pair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, l, rc, end, HIT, gaps, score, Ms, MMs, Us, W1s, template, SU;
	int hitCounter, bestSeqCount, kmersize, shifter, W1, U, M, MM;
	int *bests, *Scores;
	unsigned *values, *last, n;
	short unsigned *values_s;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	bests = bestTemplates;
	Scores = Score;
	bestSeqCount = 0;
	kmersize = templates->kmersize;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			bests = bestTemplates_r;
			Scores = Score_r;
			comp_rc(qseq);
		}
		/* Make quick check of the qseq */
		HIT = exhaustive;
		hitCounter = 0;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		/* If deltamer qseq hits, then continue */
		if(HIT) {
			/* Scan the deltamer exhaustively, and collect scores in Score*/
			last = 0;
			gaps = 0;
			HIT = 0;
			Ms = 0;
			MMs = 0;
			Us = 0;
			W1s = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
						if(values == last) {
							if(kmersize < gaps) {
								Ms += kmersize;
								gaps -= kmersize;
								if(gaps) {
									/* go for best scenario */
									if(gaps == 1) {
										MMs += 2;
									} else {
										gaps -= 2;
										if((MM << 1) + gaps * M < 0) {
											Ms += gaps;
											MMs += 2;
										}
									}
								} else {
									++MMs;
								}
							} else if (gaps) {
								--gaps;
								++W1s;
								Us += gaps;
							} else {
								++Ms;
							}
							HIT = j;
							gaps = 0;
						} else {
							if(last) {
								if(HIT) {
									HIT += kmersize;
								} else {
									HIT = j + kmersize;
								}
								score = Ms * M + MMs * MM + Us * U + W1s * W1;
								if(SU) {
									values_s = (short unsigned *) last;
									l = (*values_s) + 1;
									while(--l) {
										Scores[(template = values_s[l])] += score;
										extendScore[template] = HIT;
									}
									
									score = kmersize * M;
									MMs = MM << 1;
									values_s = (short unsigned *) values;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										if(j < extendScore[(template = values_s[l])]) {
											if(extendScore[template] == HIT) {
												Scores[template] += M;
											} else {
												gaps = extendScore[template] - j - 1;
												Scores[template] += (W1 + gaps * U);
											}
										} else if(Scores[template] != 0) {
											Scores[template] += score;
											if((gaps = extendScore[template] - j)) {
												if(gaps == 1) {
													Scores[template] += MMs;
												} else {
													gaps -= 2;
													if((Ms = MMs + gaps * M) < 0) {
														Scores[template] += Ms;
													}
												}
											} else {
												Scores[template] += MM;
											}
										} else {
											Scores[template] = score;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								} else {
									l = (*last) + 1;
									while(--l) {
										Scores[(template = last[l])] += score;
										extendScore[template] = HIT;
									}
									
									score = kmersize * M;
									MMs = MM << 1;
									n = *values;
									for(l = 1; l <= n; ++l) {
										if(j < extendScore[(template = values[l])]) {
											if(extendScore[template] == HIT) {
												Scores[template] += M;
											} else {
												gaps = extendScore[template] - j - 1;
												Scores[template] += (W1 + gaps * U);
											}
										} else if(Scores[template] != 0) {
											Scores[template] += score;
											if((gaps = extendScore[template] - j)) {
												if(gaps == 1) {
													Scores[template] += MMs;
												} else {
													gaps -= 2;
													if((Ms = MMs + gaps * M) < 0) {
														Scores[template] += Ms;
													}
												}
											} else {
												Scores[template] += MM;
											}
										} else {
											Scores[template] = score;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								}
							} else if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								Ms = kmersize * M;
								for(l = 1; l <= n; ++l) {
									Scores[(template = values_s[l])] = Ms;
									bests[l] = template;
								}
								*bests = n;
							} else {
								n = *values;
								Ms = kmersize * M;
								for(l = 1; l <= n; ++l) {
									Scores[(template = values[l])] = Ms;
									bests[l] = template;
								}
								*bests = n;
							}
							
							HIT = 0;
							gaps = 0;
							Ms = 0;
							MMs = 0;
							Us = 0;
							W1s = 0;
							last = values;
						}
						++hitCounter;
					} else {
						++gaps;
					}
				}
				j = qseq->N[i] + 1;
			}
			
			if(last) {
				score = Ms * M + MMs * MM + Us * U + W1s * W1;
				if(SU) {
					values_s = (short unsigned *) last;
					l = (*values_s) + 1;
					while(--l) {
						Scores[values_s[l]] += score;
					}
				} else {
					l = (*last) + 1;
					while(--l) {
						Scores[last[l]] += score;
					}
				}
				for(l = *bests; l != 0; --l) {
					extendScore[(template = bests[l])] = 0;
					if(Scores[template] < 0) {
						Scores[template] = 0;
					}
				}
			}
			
			if(bestSeqCount < hitCounter) {
				bestSeqCount = hitCounter;
			}
		}
		qseq->N[0]--;
	}
	
	return bestSeqCount;
}

int get_kmers_for_pair_count(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, rc, end, HIT, hitCounter, bestSeqCount, reps, SU, kmersize;
	int shifter, *bests, *Scores;
	unsigned *values, *last, n;
	short unsigned *values_s;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	bests = bestTemplates;
	Scores = Score;
	bestSeqCount = 0;
	kmersize = templates->kmersize;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			bests = bestTemplates_r;
			Scores = Score_r;
			comp_rc(qseq);
		}
		/* Make quick check of the qseq */
		HIT = exhaustive;
		hitCounter = 0;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		/* If deltamer qseq hits, then continue */
		if(HIT) {
			/* Scan the deltamer exhaustively, and collect scores in Score*/
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) last;
									n = *values_s + 1;
									while(--n) {
										if((Scores[*++values_s] += reps) == reps) {
											bests[++*bests] = *values_s;
										}
									}
								} else {
									n = *last + 1;
									while(--n) {
										if((Scores[*++last] += reps) == reps) {
											bests[++*bests] = *last;
										}
									}
								}
								hitCounter += reps;
							}
							reps = 1;
							last = values;
						}
					}
				}
				j = qseq->N[i] + 1;
			}
			
			if(last) {
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s + 1;
					while(--n) {
						if((Scores[*++values_s] += reps) == reps) {
							bests[++*bests] = *values_s;
						}
					}
				} else {
					n = *last + 1;
					while(--n) {
						if((Scores[*++last] += reps) == reps) {
							bests[++*bests] = *last;
						}
					}
				}
				hitCounter += reps;
			}
			
			reps = 0;
			if(bestSeqCount < hitCounter) {
				bestSeqCount = hitCounter;
			}
		}
		qseq->N[0]--;
	}
	
	return bestSeqCount;
}

int get_kmers_for_pair_Sparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	int i, j, n, end, rc, prefix_len, hitCounter, reps, n_kmers, kmersize;
	int HIT, SU, *bests, *Scores;
	unsigned shifter, prefix_shifter, *values, *last;
	short unsigned *values_s;
	long unsigned prefix;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->prefix_len << 1);
	
	if(*extendScore) {
		bests = bestTemplates_r;
		Scores = Score_r;
	} else {
		bests = bestTemplates;
		Scores = Score;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	hitCounter = 0;
	n_kmers = 0;
	end = qseq->seqlen;
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			if(rc) {
				comp_rc(qseq);
			}
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize - prefix_len + 1;
				while(j < end) {
					if(getKmer(qseq->seq, j, prefix_shifter) == prefix) {
						if((values = hashMap_get(templates, getKmer(qseq->seq, j + prefix_len, shifter)))) {
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s + 1;
								while(--n) {
									if(Scores[*++values_s]++ == 0) {
										bests[++*bests] = *values_s;
									}
								}
							} else {
								n = *values + 1;
								while(--n) {
									if(Scores[*++values]++ == 0) {
										bests[++*bests] = *values;
									}
								}
							}
							++hitCounter;
						}
						++n_kmers;
					}
					++j;
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		if(hitCounter) {
			hitCounter *= (((qseq->seqlen - kmersize + 1) << 1) / n_kmers);
		}
	} else {
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s + 1;
									while(--n) {
										if((Scores[*++values_s] += reps) == reps) {
											bests[++*bests] = *values_s;
										}
									}
								} else {
									n = *values + 1;
									while(--n) {
										if((Scores[*++last] += reps) == reps) {
											bests[++*bests] = *last;
										}
									}
								}
								hitCounter += reps;
							}
							reps = 1;
							last = values;
						}
					}
				}
				j = qseq->N[i] + 1;
			}
			if(last) {
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s + 1;
					while(--n) {
						if((Scores[*++values_s] += reps) == reps) {
							bests[++*bests] = *values_s;
						}
					}
				} else {
					n = *last + 1;
					while(--n) {
						if((Scores[*++last] += reps) == reps) {
							bests[++*bests] = *last;
						}
					}
				}
				hitCounter += reps;
			}
		}
		qseq->N[0]--;
	}
	
	return hitCounter;
}

int get_kmers_for_pair_pseoudoSparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, int *extendScore, const int exhaustive) {
	
	int i, j, l, n, end, template, hitCounter, gaps, Ms, MMs, Us, W1s;
	int W1, U, M, MM, HIT, SU, kmersize, score, *bests, *Scores;
	unsigned shifter, *values, *last;
	short unsigned *values_s;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return 0;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	if(*extendScore) {
		bests = bestTemplates_r;
		Scores = Score_r;
	} else {
		bests = bestTemplates;
		Scores = Score;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	hitCounter = 0;
	end = qseq->seqlen;
	
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	if(HIT) {
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		for(i = 1; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
					if(values == last) {
						if(kmersize < gaps) {
							Ms += kmersize;
							gaps -= kmersize;
							if(gaps) {
								/* go for best scenario */
								if(gaps == 1) {
									MMs += 2;
								} else {
									gaps -= 2;
									if((MM << 1) + gaps * M < 0) {
										Ms += gaps;
										MMs += 2;
									}
								}
							} else {
								++MMs;
							}
						} else if (gaps) {
							--gaps;
							++W1s;
							Us += gaps;
						} else {
							++Ms;
						}
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							if(HIT) {
								HIT += kmersize;
							} else {
								HIT = j + kmersize;
							}
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							if(SU) {
								values_s = (short unsigned *) last;
								l = *values_s + 1;
								while(--l) {
									Scores[*++values_s] += score;
									extendScore[*values_s] = HIT;
								}
							} else {
								l = *last + 1;
								while(--l) {
									Scores[*++last] += score;
									extendScore[*last] = HIT;
								}
							}
							score = kmersize * M;
							MMs = MM << 1;
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s + 1;
								while(--n) {
									template = *++values_s;
									if(j < extendScore[template]) {
										if(extendScore[template] == HIT) {
											Scores[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Scores[template] += (W1 + gaps * U);
										}
									} else if(Scores[template] != 0) {
										Scores[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Scores[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Scores[template] += Ms;
												}
											}
										} else {
											Scores[template] += MM;
										}
									} else {
										Scores[template] = score;
										bests[++*bests] = template;
									}
								}
							} else {
								n = *(last = values) + 1;
								while(--n) {
									template = *++values;
									if(j < extendScore[template]) {
										if(extendScore[template] == HIT) {
											Scores[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Scores[template] += (W1 + gaps * U);
										}
									} else if(Scores[template] != 0) {
										Scores[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Scores[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Scores[template] += Ms;
												}
											}
										} else {
											Scores[template] += MM;
										}
									} else {
										Scores[template] = score;
										bests[++*bests] = template;
									}
								}
							}
						} else if(SU) {
							values_s = (short unsigned *) values;
							n = *values_s + 1;
							Ms = kmersize * M;
							*bests = 0;
							while(--n) {
								Scores[*++values_s] = Ms;
								bests[++*bests] = *values_s;
							}
						} else {
							n = *(last = values) + 1;
							Ms = kmersize * M;
							*bests = 0;
							while(--n) {
								Scores[*++last] = Ms;
								bests[++*bests] = *last;
							}
						}
						HIT = 0;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
						last = values;
					}
					++hitCounter;
				} else {
					++gaps;
				}
			}
			j = qseq->N[i] + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = *values_s + 1;
				while(--l) {
					Scores[*++values_s] += score;
				}
			} else {
				l = *last + 1;
				while(--l) {
					Scores[*++last] += score;
				}
			}
			for(l = *bests; l != 0; --l) {
				extendScore[(template = bests[l])] = 0;
				if(Scores[template] < 0) {
					Scores[template] = 0;
				}
			}
		}
	}
	qseq->N[0]--;
	
	return hitCounter;
}

void getFirstForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, bestHits;
	
	bestHits = 0;
	
	for(i = 1; i <= *bestTemplates; ++i) {
		++bestHits;
		regionTemplates[bestHits] = bestTemplates[i];
		regionScores[bestHits] = Score[bestTemplates[i]];
		Score[bestTemplates[i]] = 0;
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		++bestHits;
		regionTemplates[bestHits] = -bestTemplates_r[i];
		regionScores[bestHits] = Score_r[bestTemplates_r[i]];
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
}

int getSecondForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, score, bestHits, bestScore;
	
	bestHits = 0;
	bestScore = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			if((score = Score[regionTemplates[i]])) {
				score += regionScores[i];
				if(bestScore < score) {
					bestScore = score;
					bestHits = 1;
					regionTemplates[bestHits] = regionTemplates[i];
				} else if(bestScore == score) {
					++bestHits;
					regionTemplates[bestHits] = regionTemplates[i];
				}
			}
		} else {
			if((score = Score_r[-regionTemplates[i]])) {
				score += regionScores[i];
				if(bestScore < score) {
					bestScore = score;
					bestHits = 1;
					regionTemplates[bestHits] = regionTemplates[i];
				} else if(bestScore == score) {
					++bestHits;
					regionTemplates[bestHits] = regionTemplates[i];
				}
			}
		}
	}
	*regionTemplates = bestHits;
	
	for(i = *bestTemplates; i != 0; --i) {
		Score[bestTemplates[i]] = 0;
	}
	for(i = *bestTemplates_r; i != 0; --i) {
		Score_r[bestTemplates_r[i]] = 0;
	}
	
	return bestScore;
}

int getFirstPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, score, bestScore, bestHits;
	
	bestScore = 0;
	bestHits = 0;
	/* save forward matches */
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore < (score = Score[bestTemplates[i]])) {
			bestScore = score;
		}
		++bestHits;
		regionTemplates[bestHits] = bestTemplates[i];
		regionScores[bestHits] = score;
		Score[bestTemplates[i]] = 0;
	}
	
	/* save reverse matches */
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore < (score = Score_r[bestTemplates_r[i]])) {
			bestScore = score;
		}
		++bestHits;
		regionTemplates[bestHits] = -bestTemplates_r[i];
		regionScores[bestHits] = score;
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getSecondPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, int bestScore, int PE) {
	
	int i, score, bestScore_r, compScore, bestHits;
	
	/* get best scoring tempates */
	bestScore_r = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < Score[bestTemplates[i]]) {
			bestScore_r = Score[bestTemplates[i]];
		}
	}
	bestHits = *bestTemplates;
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < Score_r[bestTemplates_r[i]]) {
			bestScore_r = Score_r[bestTemplates_r[i]];
		}
		++bestHits;
		bestTemplates[bestHits] = -bestTemplates_r[i];
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	if(bestScore_r) {
		compScore = bestScore + bestScore_r - PE;
		compScore = MAX(0, compScore);
		for(i = 1; i <= *regionTemplates; ++i) {
			if(0 < regionTemplates[i]) {
				/* we got one */
				if(0 < (score = Score_r[regionTemplates[i]])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
						bestHits = 1;
						regionTemplates[bestHits] = regionTemplates[i];
					} else if(compScore == score) {
						++bestHits;
						regionTemplates[bestHits] = regionTemplates[i];
					}
				}
			} else {
				/* we got one */
				if(0 < (score = Score[-regionTemplates[i]])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
						bestHits = 1;
						regionTemplates[bestHits] = regionTemplates[i];
					} else if(compScore == score) {
						++bestHits;
						regionTemplates[bestHits] = regionTemplates[i];
					}
				}
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
		/* clear scores */
		for(i = *bestTemplates; i != 0; --i) {
			if(0 < bestTemplates[i]) {
				Score[bestTemplates[i]] = 0;
			} else {
				Score_r[-bestTemplates[i]] = 0;
			}
		}
	} else {
		/* get bestHits from each as SE */
		for(i = 1; i <= *regionTemplates; ++i) {
			if(bestScore == regionScores[i]) {
				++bestHits;
				regionTemplates[bestHits] = regionTemplates[i];
			}
		}
		*regionTemplates = bestHits;
		
		bestHits = 0;
		for(i = 1; i <= *bestTemplates; ++i) {
			if(0 < bestTemplates[i]) {
				if(bestScore_r == Score[bestTemplates[i]]) {
					++bestHits;
					bestTemplates[bestHits] = bestTemplates[i];
				}
				Score[bestTemplates[i]] = 0;
			} else {
				if(bestScore_r == Score_r[-bestTemplates[i]]) {
					++bestHits;
					bestTemplates[bestHits] = bestTemplates[i];
				}
				Score_r[-bestTemplates[i]] = 0;
			}
		}
		*bestTemplates = bestHits;
	}
	
	return bestScore_r;
}

int getF_Best(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	int i, score, bestScore, bestHits;
	
	bestScore = 0;
	bestHits = 0;
	
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore < (score = Score[bestTemplates[i]])) {
			bestScore = score;
			bestHits = 1;
			regionTemplates[bestHits] = bestTemplates[i];
		} else if(bestScore == score) {
			++bestHits;
			regionTemplates[bestHits] = bestTemplates[i];
		}
		Score[bestTemplates[i]] = 0;
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore < (score = Score_r[bestTemplates_r[i]])) {
			bestScore = score;
			bestHits = 1;
			regionTemplates[bestHits] = -bestTemplates_r[i];
		} else if(bestScore == score) {
			++bestHits;
			regionTemplates[bestHits] = -bestTemplates_r[i];
		}
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getR_Best(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	int i, j, score, bestScore_r, bestHits;
	
	/* get best scoring tempates */
	bestScore_r = 0;
	bestHits = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < (score = Score[bestTemplates[i]])) {
			for(j = bestHits; j != 0; --j) {
				Score[bestTemplates[j]] = 0;
			}
			bestScore_r = score;
			bestHits = 1;
			bestTemplates[bestHits] = bestTemplates[i];
		} else if(bestScore_r == score) {
			++bestHits;
			bestTemplates[bestHits] = bestTemplates[i];
		} else {
			Score[bestTemplates[i]] = 0;
		}
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < (score = Score_r[bestTemplates_r[i]])) {
			for(j = bestHits; j != 0; --j) {
				if(0 < bestTemplates[j]) {
					Score[bestTemplates[j]] = 0;
				} else {
					Score_r[-bestTemplates[j]] = 0;
				}
			}
			bestScore_r = score;
			bestHits = 1;
			bestTemplates[bestHits] = -bestTemplates_r[i];
		} else if(bestScore_r == score) {
			++bestHits;
			bestTemplates[bestHits] = -bestTemplates_r[i];
		} else {
			Score_r[bestTemplates_r[i]] = 0;
		}
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			/* we got one */
			if(Score_r[regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		} else {
			/* we got one */
			if(Score[-regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
	}
	
	/* clear scores */
	for(i = *bestTemplates; i != 0; --i) {
		if(0 < bestTemplates[i]) {
			Score[bestTemplates[i]] = 0;
		} else {
			Score_r[-bestTemplates[i]] = 0;
		}
	}
	
	return bestScore_r;
}

void save_kmers_Sparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, j, k, l, n, end, rc, prefix_len, template, hitCounter, HIT, SU;
	int M, MM, n_kmers, score, bestScore, bestHits, reps, kmersize;
	unsigned shifter, prefix_shifter, *values, *last;
	short unsigned *values_s;
	long unsigned prefix;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	M = rewards->M;
	MM = rewards->MM;
	*bestTemplates = 0;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->prefix_len << 1);
	hitCounter = 0;
	n_kmers = 0;
	bestScore = 0;
	end = qseq->seqlen;
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			if(rc) {
				comp_rc(qseq);
			}
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize - prefix_len + 1;
				while(j < end) {
					if(getKmer(qseq->seq, j, prefix_shifter) == prefix) {
						if((values = hashMap_get(templates, getKmer(qseq->seq, j + prefix_len, shifter)))) {
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								for(k = 1; k <= n; ++k) {
									if(Score[(template = values_s[k])] == 0) {
										(*bestTemplates)++;
										bestTemplates[*bestTemplates] = template;
									}
									Score[template]++;
								}
							} else {
								n = *values;
								for(k = 1; k <= n; ++k) {
									if(Score[(template = values[k])] == 0) {
										(*bestTemplates)++;
										bestTemplates[*bestTemplates] = template;
									}
									Score[template]++;
								}
							}
							++hitCounter;
						}
						++n_kmers;
					}
					++j;
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		
		/* get best match(es) */
		if(hitCounter * kmersize > (n_kmers - hitCounter)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				score = Score[(template = bestTemplates[l])];
				score = score * kmersize * M + (n_kmers - score) * MM;
				if(score > bestScore) {
					bestScore = score;
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(score == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				Score[bestTemplates[l]] = 0;
			}
			*bestTemplates = 0;
		}
		end = n_kmers - hitCounter - bestScore;
	} else {
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										if(Score[(template = values_s[l])]) {
											Score[template] += reps;
										} else {
											Score[template] = reps;
											bestTemplates[0]++;
											bestTemplates[*bestTemplates] = template;
										}
									}
								} else {
									n = *values;
									for(l = 1; l <= n; ++l) {
										if(Score[(template = last[l])]) {
											Score[template] += reps;
										} else {
											Score[template] = reps;
											bestTemplates[0]++;
											bestTemplates[*bestTemplates] = template;
										}
									}
								}
								hitCounter += reps;
							}
							reps = 1;
							last = values;
						}
					}
				}
				j = qseq->N[i] + 1;
			}
			if(last) {
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						if(Score[(template = values_s[l])]) {
							Score[template] += reps;
						} else {
							Score[template] = reps;
							bestTemplates[0]++;
							bestTemplates[*bestTemplates] = template;
						}
					}
				} else {
					n = *last;
					for(l = 1; l <= n; ++l) {
						if(Score[(template = last[l])]) {
							Score[template] += reps;
						} else {
							Score[template] = reps;
							bestTemplates[0]++;
							bestTemplates[*bestTemplates] = template;
						}
					}
				}
				hitCounter += reps;
			}
		}
		qseq->N[0]--;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				Score[template] = 0;
			}
			*bestTemplates = 0;
		}
		end = qseq->seqlen + 1 - bestScore;
	}
	
	if(bestScore) {
		if(bestScore * kmersize > end) {
			lock(excludeOut);
			deConPrintPtr(bestTemplates, qseq, bestScore, header);
			unlock(excludeOut);
		}
	}
}

void save_kmers_pseuodeSparse(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, j, l, n, end, template, hitCounter, gaps, Ms, MMs, Us, W1s;
	int HIT, SU, score, bestScore, bestHits, kmersize;
	int W1, U, M, MM;
	unsigned shifter, *values, *last;
	short unsigned *values_s;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	*bestTemplates = 0;
	hitCounter = 0;
	bestScore = 0;
	end = qseq->seqlen;
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	if(HIT) {
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		for(i = 1; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
					if(values == last) {
						if(kmersize < gaps) {
							Ms += kmersize;
							gaps -= kmersize;
							if(gaps) {
								/* go for best scenario */
								if(gaps == 1) {
									MMs += 2;
								} else {
									gaps -= 2;
									if((MM << 1) + gaps * M < 0) {
										Ms += gaps;
										MMs += 2;
									}
								}
							} else {
								++MMs;
							}
						} else if (gaps) {
							--gaps;
							++W1s;
							Us += gaps;
						} else {
							++Ms;
						}
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							if(HIT) {
								HIT += kmersize;
							} else {
								HIT = j + kmersize;
							}
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							if(SU) {
								values_s = (short unsigned *) last;
								l = (*values_s) + 1;
								while(--l) {
									Score[(template = values_s[l])] += score;
									extendScore[template] = HIT;
								}
							} else {
								l = (*last) + 1;
								while(--l) {
									Score[(template = last[l])] += score;
									extendScore[template] = HIT;
								}
							}
							
							score = kmersize * M;
							MMs = MM << 1;
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values_s[l])]) {
										if(extendScore[template] == HIT) {
											Score[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score[template] += (W1 + gaps * U);
										}
									} else if(Score[template] != 0) {
										Score[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score[template] += Ms;
												}
											}
										} else {
											Score[template] += MM;
										}
									} else {
										Score[template] = score;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							} else {
								n = *values;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values[l])]) {
										if(extendScore[template] == HIT) {
											Score[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score[template] += (W1 + gaps * U);
										}
									} else if(Score[template] != 0) {
										Score[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score[template] += Ms;
												}
											}
										} else {
											Score[template] += MM;
										}
									} else {
										Score[template] = score;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							}
						} else if(SU) {
							values_s = (short unsigned *) values;
							n = *values_s;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score[(template = values_s[l])] = Ms;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						} else {
							n = *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score[(template = values[l])] = Ms;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						}
						HIT = 0;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
						last = values;
					}
					++hitCounter;
				} else {
					++gaps;
				}
			}
			j = qseq->N[i] + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score[last[l]] += score;
				}
			}
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				if(Score[template] < 0) {
					Score[template] = 0;
				}
			}
		}
	}
	qseq->N[0]--;
	
	/* get best match(es) */
	if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
		bestHits = 0;
		for(l = 1; l <= *bestTemplates; ++l) {
			if(Score[(template = bestTemplates[l])] > bestScore) {
				bestScore = Score[template];
				bestHits = 1;
				bestTemplates[bestHits] = template;
			} else if(Score[template] == bestScore) {
				++bestHits;
				bestTemplates[bestHits] = template;
			}
			Score[template] = 0;
			extendScore[template] = 0;
		}
		*bestTemplates = bestHits;
	} else {
		for(l = *bestTemplates; l != 0; --l) {
			extendScore[(template = bestTemplates[l])] = 0;
			Score[template] = 0;
		}
		*bestTemplates = 0;
	}
	end = qseq->seqlen + 1 - bestScore;
	
	if(bestScore) {
		if(bestScore * kmersize > end) {
			lock(excludeOut);
			deConPrintPtr(bestTemplates, qseq, bestScore, header);
			unlock(excludeOut);
		}
	}
}

void save_kmers(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, j, l, end, HIT, gaps, score, Ms, MMs, Us, W1s, W1, U, M, MM;
	int template, bestHits, hitCounter, bestScore, bestScore_r, kmersize;
	unsigned *values, *last, n, SU, shifter;
	short unsigned *values_s;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	bestScore = 0;
	bestScore_r = 0;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	W1 = rewards->W1;
	U = rewards->U;
	M = rewards->M;
	MM = rewards->MM;
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	
	/* Search forward strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates = 0;
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		end = qseq->seqlen;
		for(i = 1; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
					if(values == last) {
						if(kmersize < gaps) {
							Ms += kmersize;
							gaps -= kmersize;
							if(gaps) {
								/* go for best scenario */
								if(gaps == 1) {
									MMs += 2;
								} else {
									gaps -= 2;
									if((MM << 1) + gaps * M < 0) {
										Ms += gaps;
										MMs += 2;
									}
								}
							} else {
								++MMs;
							}
						} else if (gaps) {
							--gaps;
							++W1s;
							Us += gaps;
						} else {
							++Ms;
						}
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							if(HIT) {
								HIT += kmersize;
							} else {
								HIT = j + kmersize;
							}
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							
							if(SU) {
								values_s = (short unsigned *) last;
								l = (*values_s) + 1;
								while(--l) {
									Score[(template = values_s[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								values_s = (short unsigned *) values;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values_s[l])]) {
										if(extendScore[template] == HIT) {
											Score[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score[template] += (W1 + gaps * U);
										}
									} else if(Score[template] != 0) {
										Score[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score[template] += Ms;
												}
											}
										} else {
											Score[template] += MM;
										}
									} else {
										Score[template] = score;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							} else {
								l = (*last) + 1;
								while(--l) {
									Score[(template = last[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								n = *values;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values[l])]) {
										if(extendScore[template] == HIT) {
											Score[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score[template] += (W1 + gaps * U);
										}
									} else if(Score[template] != 0) {
										Score[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score[template] += Ms;
												}
											}
										} else {
											Score[template] += MM;
										}
									} else {
										Score[template] = score;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							}
						} else if(SU) {
							values_s = (short unsigned *) values;
							n = *values_s;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score[(template = values_s[l])] = Ms;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						} else {
							n = *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score[(template = values[l])] = Ms;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						}
						HIT = 0;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
						last = values;
					}
					++hitCounter;
				} else {
					++gaps;
				}
			}
			j = qseq->N[i] + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score[last[l]] += score;
				}
			}
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				if(Score[template] < 0) {
					Score[template] = 0;
				}
			}
		}
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				Score[template] = 0;
			}
			*bestTemplates = 0;
		}
	}
	qseq->N[0]--;
	
	/* search rc strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	j = 0;
	qseq_r->N[0]++;
	qseq_r->N[qseq_r->N[0]] = qseq_r->seqlen;
	for(i = 1; i <= qseq_r->N[0] && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq_r->seq, j, shifter))) {
				HIT = 1;
			}
		}
		j = qseq_r->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates_r = 0;
		last = 0;
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		end = qseq_r->seqlen;
		for(i = 1; i <= qseq_r->N[0]; ++i) {
			end = qseq_r->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq_r->seq, j, shifter)))) {
					if(values == last) {
						if(kmersize < gaps) {
							Ms += kmersize;
							gaps -= kmersize;
							if(gaps) {
								/* go for best scenario */
								if(gaps == 1) {
									MMs += 2;
								} else {
									gaps -= 2;
									if((MM << 1) + gaps * M < 0) {
										Ms += gaps;
										MMs += 2;
									}
								}
							} else {
								++MMs;
							}
						} else if (gaps) {
							--gaps;
							++W1s;
							Us += gaps;
						} else {
							++Ms;
						}
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							if(HIT) {
								HIT += kmersize;
							} else {
								HIT = j + kmersize;
							}
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							if(SU) {
								values_s = (short unsigned *) last;
								l = (*values_s) + 1;
								while(--l) {
									Score_r[(template = values_s[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								values_s = (short unsigned *) values;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values_s[l])]) {
										if(extendScore[template] == HIT) {
											Score_r[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score_r[template] += (W1 + gaps * U);
										}
									} else if(Score_r[template] != 0) {
										Score_r[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score_r[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score_r[template] += Ms;
												}
											}
										} else {
											Score_r[template] += MM;
										}
									} else {
										Score_r[template] = score;
										bestTemplates_r[0]++;
										bestTemplates_r[*bestTemplates_r] = template;
									}
								}
							} else {
								l = (*last) + 1;
								while(--l) {
									Score_r[(template = last[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								n = *values;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values[l])]) {
										if(extendScore[template] == HIT) {
											Score_r[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score_r[template] += (W1 + gaps * U);
										}
									} else if(Score_r[template] != 0) {
										Score_r[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score_r[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score_r[template] += Ms;
												}
											}
										} else {
											Score_r[template] += MM;
										}
									} else {
										Score_r[template] = score;
										bestTemplates_r[0]++;
										bestTemplates_r[*bestTemplates_r] = template;
									}
								}
							}
						} else if(SU) {
							values_s = (short unsigned *) values;
							n = *values_s;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score_r[(template = values_s[l])] = Ms;
								bestTemplates_r[l] = template;
							}
							*bestTemplates_r = n;
						} else {
							n = *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score_r[(template = values[l])] = Ms;
								bestTemplates_r[l] = template;
							}
							*bestTemplates_r = n;
						}
						HIT = 0;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
						last = values;
					}
					++hitCounter;
				} else {
					++gaps;
				}
			}
			j = qseq_r->N[i] + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score_r[last[l]] += score;
				}
			}
		}
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; ++l) {
				if(Score_r[(template = bestTemplates_r[l])] > bestScore_r) {
					bestScore_r = Score_r[template];
					bestHits = 1;
					bestTemplates_r[bestHits] = template;
				} else if(Score_r[template] == bestScore_r) {
					++bestHits;
					bestTemplates_r[bestHits] = template;
				}
				Score_r[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates_r = bestHits;
		} else {
			for(l = *bestTemplates_r; l != 0; --l) {
				extendScore[(template = bestTemplates_r[l])] = 0;
				Score_r[template] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	if(bestScore > 0 || bestScore_r > 0) {
		end = qseq->seqlen + 1;
		if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r))) {
			if(bestScore > bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, bestScore, header);
				unlock(excludeOut);
			} else if(bestScore < bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header);
				unlock(excludeOut);
			} else {
				/* merge */
				for(i = 1; i <= *bestTemplates_r; ++i) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, -bestScore, header);
				unlock(excludeOut);
			}
		}
	}
}

int save_kmers_intCount(const HashMapKMA *templates, int *bestTemplates, int *Score, CompDNA *qseq, unsigned *values, unsigned pos, const unsigned shifter) {
	
	unsigned n, reps, stepSize, *next;
	short unsigned *values_s;
	
	if(templates->DB_size < USHRT_MAX) {
		n = 1;
	} else {
		n = 0;
	}
	
	pos += (stepSize = templates->kmersize - 1);
	reps = 1;
	while(stepSize) {
		next = hashMap_get(templates, getKmer(qseq->seq, pos, shifter));
		if(next == values) {
			reps += stepSize;
		} else {
			if(next) {
				save_kmers_intCount(templates, bestTemplates, Score, qseq, next, pos, shifter);
			}
			stepSize >>= 1;
		}
		pos += stepSize;
	}
	
	/* score templates */
	if(n) {
		values_s = (short unsigned *) values;
		n = *values_s + 1;
		while(--n) {
			if((Score[*++values_s] += reps) == reps) {
				bestTemplates[++*bestTemplates] = *values_s;
			}
		}
	} else {
		n = *values + 1;
		while(--n) {
			if((Score[*++values] += reps) == reps) {
				bestTemplates[++*bestTemplates] = *values;
			}
		}
	}
	
	return pos + 1;
}

void save_kmers_count(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, j, l, end, HIT, bestHits, hitCounter, bestScore, bestScore_r, reps;
	int n, template, SU, kmersize;
	unsigned shifter, *values, *last;
	short unsigned *values_s;
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		return;
	} else if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	bestScore = 0;
	bestScore_r = 0;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	
	/* Search forward strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq->seq, j, shifter))) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates = 0;
		last = 0;
		reps = 0;
		j = 0;
		end = qseq->seqlen;
		for(i = 1; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq->seq, j, shifter)))) {
					if(values == last) {
						++reps;
					} else {
						if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s + 1;
								while(--n) {
									if((Score[*++values_s] += reps) == reps) {
										bestTemplates[++*bestTemplates] = *values_s;
									}
								}
							} else {
								n = *last + 1;
								while(--n) {
									if((Score[*++last] += reps) == reps) {
										bestTemplates[++*bestTemplates] = *last;
									}
								}
							}
							hitCounter += reps;
						}
						reps = 1;
						last = values;
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		if(last) {
			if(SU) {
				values_s = (short unsigned *) last;
				n = *values_s + 1;
				while(--n) {
					if((Score[*++values_s] += reps) == reps) {
						bestTemplates[++*bestTemplates] = *values_s;
					}
				}
			} else {
				n = *last + 1;
				while(--n) {
					if((Score[*++last] += reps) == reps) {
						bestTemplates[++*bestTemplates] = *last;
					}
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				Score[bestTemplates[l]] = 0;
			}
			*bestTemplates = 0;
		}
	}
	qseq->N[0]--;
	
	/* search rc strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	j = 0;
	qseq_r->N[0]++;
	qseq_r->N[qseq_r->N[0]] = qseq_r->seqlen;
	for(i = 1; i <= qseq_r->N[0] && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq_r->seq, j, shifter))) {
				HIT = 1;
			}
		}
		j = qseq_r->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates_r = 0;
		last = 0;
		reps = 0;
		j = 0;
		end = qseq_r->seqlen;
		for(i = 1; i <= qseq_r->N[0]; ++i) {
			end = qseq_r->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq_r->seq, j, shifter)))) {
					if(values == last) {
						++reps;
					} else {
						if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s + 1;
								while(--n) {
									if((Score_r[*++values_s] += reps) == reps) {
										bestTemplates_r[++*bestTemplates_r] = *values_s;
									}
								}
							} else {
								n = *last + 1;
								while(--n) {
									if((Score_r[*++last] += reps) == reps) {
										bestTemplates_r[++*bestTemplates_r] = *last;
									}
								}
							}
							hitCounter += reps;
						}
						last = values;
						reps = 1;
					}
				}
			}
			j = qseq_r->N[i] + 1;
		}
		if(last) {
			if(SU) {
				values_s = (short unsigned *) last;
				n = *values_s + 1;
				while(--n) {
					if((Score_r[*++values_s] += reps) == reps) {
						bestTemplates_r[++*bestTemplates_r] = *values_s;
					}
				}
			} else {
				n = *last + 1;
				while(--n) {
					if((Score_r[*++last] += reps) == reps) {
						bestTemplates_r[++*bestTemplates_r] = *last;
					}
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; ++l) {
				if(Score_r[(template = bestTemplates_r[l])] > bestScore_r) {
					bestScore_r = Score_r[template];
					bestHits = 1;
					bestTemplates_r[bestHits] = template;
				} else if(Score_r[template] == bestScore_r) {
					++bestHits;
					bestTemplates_r[bestHits] = template;
				}
				Score_r[template] = 0;
			}
			*bestTemplates_r = bestHits;
		} else {
			for(l = *bestTemplates_r; l != 0; --l) {
				Score_r[bestTemplates_r[l]] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	if(bestScore > 0 || bestScore_r > 0) {
		end = qseq->seqlen + 1;
		if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r))) {
			if(bestScore > bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, bestScore, header);
				unlock(excludeOut);
			} else if(bestScore < bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header);
				unlock(excludeOut);
			} else {
				/* merge */
				for(i = 1; i <= *bestTemplates_r; ++i) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, -bestScore, header);
				unlock(excludeOut);
			}
		}
	}
}

void save_kmers_unionPair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, const Qseqs *header_r, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, bestScore, bestScore_r, hitCounter, kmersize;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	kmersize = templates->kmersize;
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive)) && (qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
		/* got hits */
		bestScore = getF_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		
		if(bestScore * kmersize < (qseq->seqlen - bestScore)) {
			bestScore = 0;
		}
	} else {
		bestScore = 0;
		if(hitCounter) {
			i = *bestTemplates + 1;
			while(--i) {
				Score[bestTemplates[i]] = 0;
			}
			*bestTemplates = 0;
			i = *bestTemplates_r + 1;
			while(--i) {
				Score_r[bestTemplates_r[i]] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore, exhaustive)) && (qseq_r->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
		if(bestScore) {
			bestScore_r = getR_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		} else {
			bestScore_r = getF_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		}
		if(bestScore_r * kmersize < (qseq_r->seqlen - bestScore_r)) {
			bestScore_r = 0;
			*regionTemplates = abs(*regionTemplates);
		}
	} else {
		bestScore_r = 0;
		i = *bestTemplates + 1;
		while(--i) {
			Score[bestTemplates[i]] = 0;
		}
		*bestTemplates = 0;
		i = *bestTemplates_r + 1;
		while(--i) {
			Score_r[bestTemplates_r[i]] = 0;
		}
		*bestTemplates_r = 0;
	}
	
	if(0 < bestScore && 0 < bestScore_r) {
		if(*regionTemplates < 0) {
			*regionTemplates = -(*regionTemplates);
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
					bestScore_r = -bestScore_r;
				}
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore_r, header_r);
				unlock(excludeOut);
			} else {
				comp_rc(qseq_r);
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq_r, bestScore_r, header_r, qseq, bestScore, header);
				unlock(excludeOut);
			}
		} else {
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
				}
			} else {
				for(i = 1; i <= *regionTemplates; ++i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(regionTemplates, qseq, bestScore, header);
			unlock(excludeOut);
			if(0 < bestTemplates[1]) {
				comp_rc(qseq_r);
				if(bestTemplates[*bestTemplates] < 0) {
					bestScore_r = -bestScore_r;
				}
			} else {
				for(i = 1; i <= *bestTemplates; ++i) {
					bestTemplates[i] = -bestTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(bestTemplates, qseq_r, bestScore_r, header_r);
			unlock(excludeOut);
		}
	} else if(bestScore) {
		if(0 < regionTemplates[1]) {
			comp_rc(qseq);
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore = -bestScore;
			}
		} else {
			for(i = 1; i <= *regionTemplates; ++i) {
				regionTemplates[i] = -regionTemplates[i];
			}
		}
		lock(excludeOut);
		deConPrintPtr(regionTemplates, qseq, bestScore, header);
		unlock(excludeOut);
	} else if(bestScore_r) {
		if(0 < regionTemplates[1]) {
			comp_rc(qseq_r);
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore_r = -bestScore_r;
			}
		} else {
			for(i = 1; i <= *regionTemplates; ++i) {
				regionTemplates[i] = -regionTemplates[i];
			}
		}
		lock(excludeOut);
		deConPrintPtr(regionTemplates, qseq_r, bestScore_r, header_r);
		unlock(excludeOut);
	}
}

void save_kmers_penaltyPair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, const Qseqs *header_r, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, bestScore, bestScore_r, compScore, hitCounter, hitCounter_r;
	int kmersize;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	kmersize = templates->kmersize;
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive))) {
		/* got hits */
		bestScore = getFirstPen(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores);
	} else {
		bestScore = 0;
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter_r = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore, exhaustive))) {
		if(0 < bestScore) {
			bestScore_r = getSecondPen(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores, bestScore, rewards->PE);
		} else {
			bestScore_r = getF_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		}
	} else {
		bestScore_r = 0;
	}
	
	if(0 < bestScore && 0 < bestScore_r) {
		if(*regionTemplates < 0) {
			compScore = MIN((hitCounter + hitCounter_r), (bestScore + bestScore_r));
			if((qseq->seqlen + qseq_r->seqlen - compScore - (kmersize << 1)) < compScore * kmersize) {
				*regionTemplates = -(*regionTemplates);
				if(0 < regionTemplates[1]) {
					comp_rc(qseq);
					if(regionTemplates[*regionTemplates] < 0) {
						bestScore = -bestScore;
						bestScore_r = -bestScore_r;
					}
					lock(excludeOut);
					printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore_r, header_r);
					unlock(excludeOut);
				} else {
					comp_rc(qseq_r);
					for(i = *regionTemplates; i != 0; --i) {
						regionTemplates[i] = -regionTemplates[i];
					}
					lock(excludeOut);
					printPairPtr(regionTemplates, qseq_r, bestScore_r, header_r, qseq, bestScore, header);
					unlock(excludeOut);
				}
			}
		} else {
			hitCounter = MIN(hitCounter, bestScore);
			if((qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
				if(0 < regionTemplates[1]) {
					comp_rc(qseq);
					if(regionTemplates[*regionTemplates] < 0) {
						bestScore = -bestScore;
					}
				} else {
					for(i = *regionTemplates; i != 0; --i) {
						regionTemplates[i] = -regionTemplates[i];
					}
				}
				lock(excludeOut);
				deConPrintPtr(regionTemplates, qseq, bestScore, header);
				unlock(excludeOut);
			}
			hitCounter_r = MIN(hitCounter_r, bestScore_r);
			if((qseq_r->seqlen - hitCounter_r - kmersize) < hitCounter_r * kmersize) {
				if(0 < bestTemplates[1]) {
					comp_rc(qseq_r);
					if(bestTemplates[*bestTemplates] < 0) {
						bestScore_r = -bestScore_r;
					}
				} else {
					for(i = *bestTemplates; i != 0; --i) {
						bestTemplates[i] = -bestTemplates[i];
					}
				}
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq_r, bestScore_r, header_r);
				unlock(excludeOut);
			}
		}
	} else if(0 < bestScore) {
		hitCounter = MIN(hitCounter, bestScore);
		if((qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
				}
			} else {
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(regionTemplates, qseq, bestScore, header);
			unlock(excludeOut);
		}
	} else if(0 < bestScore_r) {
		hitCounter_r = MIN(hitCounter_r, bestScore_r);
		if((qseq_r->seqlen - hitCounter_r - kmersize) < hitCounter_r * kmersize) {
			if(0 < regionTemplates[1]) {
				comp_rc(qseq_r);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore_r = -bestScore_r;
				}
			} else {
				for(i = 1; i <= *regionTemplates; ++i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(regionTemplates, qseq_r, bestScore_r, header_r);
			unlock(excludeOut);
		}
	}
}

void save_kmers_forcePair(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, const Qseqs *header_r, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	int i, bestScore, hitCounter, hitCounter_r, kmersize;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	kmersize = templates->kmersize;
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(templates, rewards, bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore, exhaustive))) {
		/* got hits */
		getFirstForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores);
	} else {
		return;
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter_r = get_kmers_for_pair_ptr(templates, rewards, bestTemplates_r, bestTemplates, Score_r, Score, qseq_r, extendScore, exhaustive)) && 
		(qseq->seqlen + qseq_r->seqlen - hitCounter - hitCounter_r - (kmersize << 1)) < (hitCounter + hitCounter_r) * kmersize && 
	(bestScore = getSecondForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores))) {
		
		if((qseq->seqlen + qseq_r->seqlen - bestScore) < bestScore * kmersize) {
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore = -bestScore;
			}
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore, header_r);
				unlock(excludeOut);
			} else {
				comp_rc(qseq_r);
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq_r, bestScore, header_r, qseq, bestScore, header);
				unlock(excludeOut);
			}
		}
	} else if(hitCounter || hitCounter_r) {
		i = *bestTemplates + 1;
		while(--i) {
			Score[bestTemplates[i]] = 0;
		}
		*bestTemplates = 0;
		i = *bestTemplates_r + 1;
		while(--i) {
			Score_r[bestTemplates_r[i]] = 0;
		}
		*bestTemplates_r = 0;
	}
}

void save_kmers_HMM(const HashMapKMA *templates, const Penalties *rewards, int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, CompDNA *qseq, CompDNA *qseq_r, const Qseqs *header, int *extendScore, const int exhaustive, volatile int *excludeOut) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	and is the time determining step */
	static int **RegionTemplates, **TmpNs, *Sizes, *template_lengths, minLen = 0;
	static unsigned ***tVF_scores, ***tVR_scores;
	int i, j, k, l, N, n, i_r, j_r, seqlen, seqend, end, HIT, Ncheck, SU;
	int hitCounter, template, bestScore, bestHits, start, stop, kmersize;
	int start_cut, end_cut, DB_size, deCon, *regionTemplates;
	unsigned *values, *last, *rlast, **VF_scores, **VR_scores, num, shifter;
	short unsigned *values_s;
	int *tmpNs, reps, rreps;
	double Ms, Ns, Ms_prev, Ns_prev, HMM_param[8];
	
	if(qseq->seqlen < (kmersize = templates->kmersize)) {
		if(minLen == 0) {
			/* initial allocate */
			template_lengths = bestTemplates_r;
			n = *bestTemplates;
			Sizes = smalloc(n * sizeof(int));
			tVF_scores = smalloc(n * sizeof(int **));
			tVR_scores = smalloc(n * sizeof(int **));
			RegionTemplates = smalloc(n * sizeof(int*));
			TmpNs = smalloc(n * sizeof(int *));
			i = n;
			while(i--) {
				RegionTemplates[i] = smalloc(((templates->DB_size << 1) + 4) * sizeof(int));
				Sizes[i] = 256;
				TmpNs[i] = smalloc(256 * sizeof(int));
				tVF_scores[i] = calloc(256, sizeof(unsigned *));
				tVR_scores[i] = calloc(256, sizeof(unsigned *));
				if(!tVF_scores[i] || !tVR_scores[i]) {
					ERROR();
				}
			}
			
			/* set minLen */
			minLen = template_lengths[1] / 2;
			for(i = 1; i < templates->DB_size; ++i) {
				if(minLen > (template_lengths[i] / 2)) {
					minLen = template_lengths[i] / 2;
				}
			}
			minLen = MIN(minLen, templates->kmersize);
		}
		return;
	} else if((DB_size = templates->DB_size) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	deCon = deConPrintPtr == &deConPrint;
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1);
	
	/* calculate HMM parameters */
	HMM_param[0] = log(1 - pow(0.25, kmersize));
	HMM_param[1] = log(pow(0.25, kmersize));
	HMM_param[2] = log(1 - pow(0.25, kmersize - 1) * 0.75);
	HMM_param[3] = log(pow(0.25, kmersize - 1) * 0.75);
	HMM_param[4] = log(1 - 1.0 / kmersize * 0.75 * 0.25);
	HMM_param[5] = log(1.0 / kmersize * 0.75 * 0.25);
	HMM_param[6] = log(0.75);
	HMM_param[7] = log(0.25);
	
	if(Sizes[*Score] < qseq->size) {
		Sizes[*Score] = qseq->size;
		free(tVF_scores[*Score]);
		free(tVR_scores[*Score]);
		free(TmpNs[*Score]);
		TmpNs[*Score] = smalloc(qseq->size * sizeof(int));
		tVF_scores[*Score] = calloc(qseq->size, sizeof(unsigned *));
		tVR_scores[*Score] = calloc(qseq->size, sizeof(unsigned *));
		if(!tVF_scores[*Score] || !tVR_scores[*Score]) {
			ERROR();
		}
	}
	regionTemplates = RegionTemplates[*Score];
	*regionTemplates++ = templates->DB_size;
	*regionTemplates++ = *Score;
	*regionTemplates++ = 0;
	VF_scores = tVF_scores[*Score];
	VR_scores = tVR_scores[*Score];
	tmpNs = TmpNs[*Score];
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	seqlen = qseq->seqlen;
	seqend = seqlen - kmersize + 1;
	i = 0;
	i_r = seqlen - kmersize;
	N = 1;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = seqlen;
	while(N <= qseq->N[0]) {
		/* find a seed */
		HIT = 0;
		end = qseq->N[N] - kmersize + 1;
		if(exhaustive) {
			while(i < end && !HIT) {
				if(hashMap_get(templates, getKmer(qseq->seq, i, shifter)) || hashMap_get(templates, getKmer(qseq_r->seq, i_r, shifter))) {
					HIT = 1;
				} else {
					++i;
					--i_r;
				}
			}
		} else {
			while(i < end && !HIT) {
				if(hashMap_get(templates, getKmer(qseq->seq, i, shifter)) || hashMap_get(templates, getKmer(qseq_r->seq, i_r, shifter))) {
					HIT = 1;
				} else {
					i += kmersize;
					i_r -= kmersize;
				}
			}
		}
		
		/* evaluate seed */
		if(HIT) {
			/* set scores attr */
			bestScore = 0;
			*bestTemplates = 0;
			hitCounter = 1;
			
			/* save seed */
			VF_scores[i] = hashMap_get(templates, getKmer(qseq->seq, i, shifter));
			VR_scores[i] = hashMap_get(templates, getKmer(qseq_r->seq, i_r, shifter));
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			Ms = 0;
			Ns = 0;
			
			/* extend backward */
			j = i - 1;
			j_r = i_r + 1;
			n = N - 1;
			Ncheck = (n > 0) ? -1 : qseq->N[n];
			while(j >= 0) {
				if(j == Ncheck) {
					
					k = j;
					while(k >= kmersize && k < (j - kmersize)) {
						/* update next N check */
						if(k == Ncheck) {
							j = Ncheck;
							--n;
							Ncheck = (n > 0) ? -1 : qseq->N[n];
						}
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							break;
						}
						--k;
						Ns_prev = Ns;
						Ms_prev = Ms;
					}
					
					if(k >= kmersize && k < (j - kmersize)) {
						j = k - 1;
						break;
					} else {
						j = k;
						j_r = seqlen - kmersize - k;
					}
				} else {
					VF_scores[j] = hashMap_get(templates, getKmer(qseq->seq, j, shifter));
					VR_scores[j] = hashMap_get(templates, getKmer(qseq_r->seq, j_r, shifter));
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						++hitCounter;
						
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[5] >= Ms_prev + HMM_param[3] + HMM_param[5]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[5];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[5];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[7] >= Ns_prev + HMM_param[1] + HMM_param[7]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[7];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[7];
							--j;
							break;
						}
					} else {
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							--j;
							break;
						}
					}
				}
				--j;
				++j_r;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			start = j + 1;
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			
			/* extend forward */
			j = i + 1;
			j_r = i_r - 1;
			Ncheck = qseq->N[N] - kmersize + 1;
			while(j < seqend) {
				if(j == Ncheck) {
					
					k = j;
					while(k < seqend && k < (j + kmersize)) {
						/* update next N check */
						if(k == Ncheck) {
							j = Ncheck;
							++N;
							Ncheck = (N == qseq->N[0]) ? seqlen : qseq->N[N] - kmersize + 1;
						}
						
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							break;
						}
						++k;
						Ns_prev = Ns;
						Ms_prev = Ms;
					}
					
					if(k < seqend && k < (j + kmersize)) {
						j = k;
						break;
					} else {
						j = k;
						j_r = seqlen - kmersize - k;
					}
				} else {
					VF_scores[j] = hashMap_get(templates, getKmer(qseq->seq, j, shifter));
					VR_scores[j] = hashMap_get(templates, getKmer(qseq_r->seq, j_r, shifter));
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						++hitCounter;
						
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[5] >= Ms_prev + HMM_param[3] + HMM_param[5]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[5];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[5];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[7] >= Ns_prev + HMM_param[1] + HMM_param[7]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[7];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[7];
							++j;
							break;
						}
					} else {
						/* Non match seq */
						if(Ns_prev + HMM_param[0] + HMM_param[4] >= Ms_prev + HMM_param[3] + HMM_param[4]) {
							Ns = Ns_prev + HMM_param[0] + HMM_param[4];
						} else {
							Ns = Ms_prev + HMM_param[3] + HMM_param[4];
						}
						
						/* Match seq */
						if(Ms_prev + HMM_param[2] + HMM_param[6] >= Ns_prev + HMM_param[1] + HMM_param[6]) {
							Ms = Ms_prev + HMM_param[2] + HMM_param[6];
						} else {
							Ms = Ns_prev + HMM_param[1] + HMM_param[6];
							++j;
							break;
						}
					}
				}
				++j;
				--j_r;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			stop = j + kmersize - 1;
			
			/* evaluate hit */
			if(hitCounter > 0 && (hitCounter * kmersize > (stop - start - hitCounter + kmersize)) 
				&& ((stop - start) > minLen || start == 0 || stop == seqlen)) {
				if(deCon) {
					if(SU) {
						for(k = start; k < j; ++k) {
							if(((values_s = (short unsigned *) VF_scores[k]) && values_s[*values_s] == DB_size) 
								|| ((values_s = (short unsigned *) VR_scores[k]) && values_s[*values_s] == DB_size)) {
								--hitCounter;
							}
						}
					} else {
						for(k = start; k < j; ++k) {
							if(((values = VF_scores[k]) && values[*values] == DB_size)
								|| ((values = VR_scores[k]) && values[*values] == DB_size)) {
								--hitCounter;
							}
						}
					}
				}
				
				/* accept hit */
				if(hitCounter > 0) {
					/* gain total scores and mapping templates for this region */
					*bestTemplates = 0;
					*bestTemplates_r = 0;
					
					last = 0;
					reps = 0;
					rlast = 0;
					rreps = 0;
					for(k = start; k < j; ++k) {
						/* forward */
						if(VF_scores[k]) {
							if(VF_scores[k] == last) {
								++reps;
							} else {
								if(last) {
									if(SU) {
										values_s = (short unsigned *) last;
										num = *values_s + 1;
										while(--num) {
											if((Score[*++values_s] += reps) == reps) {
												bestTemplates[++*bestTemplates] = *values_s;
											}
										}
									} else {
										n = *last + 1;
										while(--n) {
											if((Score[*++last] += reps) == reps) {
												bestTemplates[++*bestTemplates] = *last;
											}
										}
									}
								}
								reps = 1;
								last = VF_scores[k];
							}
						}
						
						/* rc */
						if(VR_scores[k]) {
							if(VR_scores[k] == rlast) {
								++rreps;
							} else {
								if(rlast) {
									if(SU) {
										values_s = (short unsigned *) rlast;
										num = *values_s + 1;
										while(--num) {
											if((Score_r[*++values_s] += rreps) == rreps) {
												bestTemplates_r[++*bestTemplates_r] = *values_s;
											}
										}
									} else {
										num = *rlast + 1;
										while(--num) {
											if((Score_r[*++rlast] += rreps) == rreps) {
												bestTemplates_r[++*bestTemplates_r] = *rlast;
											}
										}
									}
								}
								rreps = 1;
								rlast = VR_scores[k];
							}
						}
						
					}
					if(last) {
						if(SU) {
							values_s = (short unsigned *) last;
							num = *values_s + 1;
							while(--num) {
								if((Score[*++values_s] += reps) == reps) {
									bestTemplates[++*bestTemplates] = *values_s;
								}
							}
						} else {
							num = *last + 1;
							while(--num) {
								if((Score[*++last] += reps) == reps) {
									bestTemplates[++*bestTemplates] = *last;
								}
							}
						}
					}
					if(rlast) {
						if(SU) {
							values_s = (short unsigned *) rlast;
							num = *values_s + 1;
							while(--num) {
								if((Score_r[*++values_s] += rreps) == rreps) {
									bestTemplates_r[++*bestTemplates_r] = *values_s;
								}
							}
						} else {
							num = *rlast + 1;
							while(--num) {
								if((Score_r[*++rlast] += rreps) == rreps) {
									bestTemplates_r[++*bestTemplates_r] = *rlast;
								}
							}
						}
					}
					
					/* cut out template hits */
					while(HIT != 0) {
						/* get best score */
						bestScore = 0;
						bestHits = 0;
						
						/* forward */
						for(k = 1; k <= *bestTemplates; ++k) {
							template = bestTemplates[k];
							if(Score[template] > bestScore) {
								bestScore = Score[template];
								bestHits = 1;
								regionTemplates[bestHits] = template;
							} else if(Score[template] == bestScore) {
								if(Score[template]) {
									++bestHits;
									regionTemplates[bestHits] = template;
								} else {
									bestTemplates[k] = bestTemplates[*bestTemplates];
									bestTemplates[0]--;
									--k;
								}
							}
						}
						
						/* rc */
						for(k = 1; k <= *bestTemplates_r; ++k) {
							template = bestTemplates_r[k];
							if(Score_r[template] > bestScore) {
								bestScore = Score_r[template];
								bestHits = 1;
								regionTemplates[bestHits] = -template;
							} else if(Score_r[template] == bestScore) {
								if(Score_r[template]) {
									++bestHits;
									regionTemplates[bestHits] = -template;
								} else {
									bestTemplates_r[k] = bestTemplates_r[*bestTemplates_r];
									bestTemplates_r[0]--;
									--k;
								}
							}
						}
						
						*regionTemplates = bestHits;
						if(bestScore > 0) {
							/* find limits of match */
							start_cut = j;
							for(k = 1; k <= bestHits; ++k) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = start; l < start_cut; ++l) {
									if(VR_scores[l] && intpos_bin_contaminationPtr(VR_scores[l], template) != -1) {
										start_cut = l;
									}
									if(VF_scores[l] && intpos_bin_contaminationPtr(VF_scores[l], template) != -1) {
										start_cut = l;
									}
								}
							}
							end_cut = start_cut;
							for(k = 1; k <= bestHits; ++k) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = j; l > end_cut; --l) {
									if(VR_scores[l] && intpos_bin_contaminationPtr(VR_scores[l], template) != -1) {
										end_cut = l;
									}
									if(VF_scores[l] && intpos_bin_contaminationPtr(VF_scores[l], template) != -1) {
										end_cut = l;
									}
								}
							}
							
							/* evaluate best hit */
							if(bestScore * kmersize > (end_cut - start_cut - bestScore + kmersize)) {
								/* check for hits on rc */
								HIT = (regionTemplates[*regionTemplates] > 0) ? 1 : -1;
								/* print */
								if(start != 0 && j != qseq->seqlen) {
									ankerAndClean(regionTemplates, Score, Score_r, template_lengths, VF_scores, VR_scores, tmpNs, qseq, HIT, bestScore, start_cut, end_cut, header, excludeOut);
								} else {
									ankerPtr(regionTemplates, Score, Score_r, template_lengths, VF_scores, VR_scores, tmpNs, qseq, HIT, bestScore, start_cut, end_cut, header, excludeOut);
								}
							} else {
								/* clear scores */
								for(k = 1; k <= *bestTemplates; ++k) {
									Score[bestTemplates[k]] = 0;
								}
								for(k = 1; k <= *bestTemplates_r; ++k) {
									Score_r[bestTemplates_r[k]] = 0;
								}
								HIT = 0;
							}
						} else {
							/* clear scores */
							for(k = 1; k <= *bestTemplates; ++k) {
								Score[bestTemplates[k]] = 0;
							}
							for(k = 1; k <= *bestTemplates_r; ++k) {
								Score_r[bestTemplates_r[k]] = 0;
							}
							
							HIT = 0;
						}
					}
				}
			}
			
			/* clear scores */
			for(k = start; k < j; ++k) {
				VF_scores[k] = 0;
				VR_scores[k] = 0;
			}
			
			i = stop + 1;
			i_r = seqlen - kmersize - i;
		} else {
			++N;
		}
	}
}

void ankerAndClean(int *regionTemplates, int *Score, int *Score_r, int *template_lengths, unsigned **VF_scores, unsigned **VR_scores, int *tmpNs, CompDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, const Qseqs *header, volatile int *excludeOut) {
	
	int k, l, bestHitsCov, template, DB_size;
	unsigned *values, n, SU;
	short unsigned *values_s;
	double thisCov, bestCov;
	CompDNA tmpQseq;
	
	if((DB_size = regionTemplates[-3]) < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	/* here */
	// make sure cuts isn't random seeds
	
	bestHitsCov = *regionTemplates;
	bestCov = 0;
	for(k = 1; k <= *regionTemplates; ++k) {
		template = regionTemplates[k];
		if(template < 0) {
			template = -template;
			thisCov = 1.0 * Score_r[template] / template_lengths[template];
		} else {
			thisCov = 1.0 * Score[template] / template_lengths[template];
		}
		if(thisCov > bestCov) {
			bestCov = thisCov;
		}
	}
	
	for(k = start_cut; k <= end_cut; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				n = *values_s;
				for(l = 1; l <= n; ++l) {
					template = values_s[l];
					if(Score[template] != bestScore && template != DB_size) {
						thisCov = 1.0 * Score[template] / template_lengths[template];
						if(thisCov > bestCov) {
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = template;
						} else if(thisCov == bestCov) {
							++bestHitsCov;
							regionTemplates[bestHitsCov] = template;
						}
					}
					Score[template]--;
				}
			} else {
				values = VF_scores[k];
				n = *values;
				for(l = 1; l <= n; ++l) {
					template = values[l];
					if(Score[template] != bestScore && template != DB_size) {
						thisCov = 1.0 * Score[template] / template_lengths[template];
						if(thisCov > bestCov) {
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = template;
						} else if(thisCov == bestCov) {
							++bestHitsCov;
							regionTemplates[bestHitsCov] = template;
						}
					}
					Score[template]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				n = *values_s;
				for(l = 1; l <= n; ++l) {
					template = values_s[l];
					if(Score_r[template] != bestScore && template != DB_size) {
						thisCov = 1.0 * Score_r[template] / template_lengths[template];
						if(thisCov > bestCov) {
							HIT = -1;
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = -template;
						} else if(thisCov == bestCov) {
							HIT = -1;
							++bestHitsCov;
							regionTemplates[bestHitsCov] = -template;
						}
					}
					Score_r[template]--;
				}
			} else {
				values = VR_scores[k];
				n = *values;
				for(l = 1; l <= n; ++l) {
					template = values[l];
					if(Score_r[template] != bestScore && template != DB_size) {
						thisCov = 1.0 * Score_r[template] / template_lengths[template];
						if(thisCov > bestCov) {
							HIT = -1;
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = -template;
						} else if(thisCov == bestCov) {
							HIT = -1;
							++bestHitsCov;
							regionTemplates[bestHitsCov] = -template;
						}
					}
					Score_r[template]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	*regionTemplates = bestHitsCov;
	
	/* clear nearest templates on both sides of match */
	for(k = ((start_cut - 92) < 0) ? 0 : (start_cut - 92); k < start_cut; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	for(k = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92); k > end_cut; --k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	
	
	
	/* modify limits of match seq */
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	start_cut = (start_cut >> 5) << 5;
	end_cut = ((end_cut >> 5) << 5) + 32;
	end_cut = (end_cut < qseq->seqlen) ? end_cut : qseq->seqlen;
	tmpQseq.seqlen = (end_cut - start_cut);
	tmpQseq.seq = qseq->seq + (start_cut >> 5);
	tmpQseq.N = tmpNs;
	
	for(k = 1, l = 0; k < qseq->N[0]; ++k) {
		if(start_cut <= qseq->N[k]) {
			++l;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				--l;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	--tmpQseq.seqlen;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		--tmpQseq.seqlen;
		--l;
	}
	++tmpQseq.seqlen;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	lock(excludeOut);
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header);
	unlock(excludeOut);
}

void ankerAndClean_MEM(int *regionTemplates, int *Score, int *Score_r, int *template_lengths, unsigned **VF_scores, unsigned **VR_scores, int *tmpNs, CompDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, const Qseqs *header, volatile int *excludeOut) {
	
	int k, l, SU;
	unsigned *values;
	short unsigned *values_s;
	CompDNA tmpQseq;
	
	if(regionTemplates[-3] < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	/* clean up scores */
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	for(k = start_cut; k < end_cut; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	
	/* modify limits of match seq */
	start_cut = (start_cut >> 5) << 5;
	end_cut = ((end_cut >> 5) << 5) + 32;
	end_cut = (end_cut < qseq->seqlen) ? end_cut : qseq->seqlen;
	tmpQseq.seqlen = (end_cut - start_cut);
	tmpQseq.seq = qseq->seq + (start_cut >> 5);
	tmpQseq.N = tmpNs;
	
	for(k = 1, l = 0; k < qseq->N[0]; ++k) {
		if(start_cut <= qseq->N[k]) {
			++l;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				--l;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	--tmpQseq.seqlen;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		--tmpQseq.seqlen;
		--l;
	}
	++tmpQseq.seqlen;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	lock(excludeOut);
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header);
	unlock(excludeOut);
}
