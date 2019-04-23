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

#include "stdnuc.h"

long unsigned getKmer(long unsigned *compressor, unsigned cPos, const unsigned shifter) {
	
	unsigned iPos = (cPos & 31) << 1;
	cPos >>= 5;
	
	return (iPos <= shifter) ? ((compressor[cPos] << iPos) >> shifter) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifter);
}

long unsigned makeKmer(const unsigned char *qseq, unsigned pos, unsigned size) {
	
	long unsigned key = qseq[pos];
	
	size += pos;
	for(++pos; pos < size; ++pos) {
		key = (key << 2) | qseq[pos];
	}
	
	return key;
}

int charpos(const unsigned char *src, unsigned char target, int start, int len) {
	
	unsigned char *ptr;
	
	ptr = (unsigned char *) src + --start;
	while(++start < len) {
		if(*++ptr == target) {
			return start;
		}
	}
	
	return -1;
}

void strrc(unsigned char *qseq, int q_len) {
	
	int i, j, seqlen;
	unsigned char carry, comp[6] = {3, 2, 1, 0, 4, 5};
	
	seqlen = q_len >> 1;
	
	for(i = 0, j = q_len - 1; i < seqlen; ++i, --j) {
		carry = comp[qseq[i]];
		qseq[i] = comp[qseq[j]];
		qseq[j] = carry;
	}
	if(q_len & 1) {
		qseq[seqlen] = comp[qseq[seqlen]];
	}
	
}

void strtranslate(unsigned char *strp, char *trans) {
	--strp;
	while(*++strp) {
		*strp = trans[*strp];
	}
}

void nibble2base(unsigned char *seq, int len) {
	
	const char bases[6] = "ACGTN-";
	
	seq += len;
	*seq-- = 0;
	++len;
	while(--len) {
		*seq = bases[*seq];
		--seq;
	}	
}
