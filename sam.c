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

#include <stdio.h>
#include <stdlib.h>
#include "nw.h"
#include "pherror.h"
#include "qseqs.h"
#include "runkma.h"

int makeCigar(Qseqs *Cigar, const Aln *aligned) {
	
	int len, cLen, rep, consume;
	char op, pop, *s, *cigar;
	unsigned char *t, *q;
	
	if(Cigar->size < aligned->len << 1) {
		Cigar->size <<= 1;
		free(Cigar->seq);
		Cigar->seq = smalloc(Cigar->size);
	} else if(aligned->len == 0) {
		return 0;
	}
	
	len = aligned->len + 1;
	t = aligned->t;
	s = aligned->s;
	q = aligned->q;
	cigar = (char *) Cigar->seq;
	
	cLen = 0;
	rep = 1;
	if(*s == '|') {
		pop = '=';
	} else if(*t == '-') {
		pop = 'I';
	} else if(*q == '-') {
		pop = 'D';
	} else {
		pop = 'X';
	}
	if(pop != 'I') {
		consume = 1;
	} else {
		consume = -1;
	}
	while(--len) {
		if(*s == '|') {
			op = '=';
		} else if(*t == '-') {
			op = 'I';
		} else if(*q == '-') {
			op = 'D';
		} else {
			op = 'X';
		}
		if(0 < consume) {
			--consume;
			if(op != 'I') {
				consume = -consume;
			}
		}
		if(op == pop) {
			++rep;
		} else {
			cLen += sprintf(cigar + cLen, "%d%c", rep, pop);
			rep = 1;
		}
	}
	Cigar->len = cLen + sprintf(cigar + cLen, "%d%c", rep, pop);
	consume = consume < 0 ? 0 : consume;
	
	return consume;
}

void saminit(Qseqs *template_name, FILE *name_file, int *template_lengths, int DB_size) {
	
	while(--DB_size) {
		fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n", nameLoad(template_name, name_file), *++template_lengths);
	}
	fseek(name_file, 0, SEEK_SET);
}

int samwrite(const Qseqs *qseq, const Qseqs *header, const Qseqs *Qual, char *rname, const Aln *aligned, const int *stats, Qseqs *Cigar) {
	
	int flag, pos, mapQ, pnext, tlen;
	char *qname, *cigar, *rnext, *seq, *qual;
	
	
	qname = (char *) header->seq;
	seq = (char *) qseq->seq;
	if(Qual) {
		qual = (char *) Qual->seq;
	} else {
		qual = "*";
	}
	
	if(rname) {
		mapQ = 254 < aligned->mapQ ? 254 : aligned->mapQ;
		pos = stats[2] + makeCigar(Cigar, aligned);
		tlen = stats[3] - pos;
		flag = stats[4];
		cigar = (char *) Cigar->seq;
	} else {
		rname = "*";
		mapQ = 0;
		pos = 0;
		tlen = 0;
		flag = *stats;
		cigar = "*";
	}
	rnext = "*";
	pnext = 0;
	
	return fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n", qname, flag, rname, pos, mapQ, cigar, rnext, pnext, tlen, seq, qual);
}
