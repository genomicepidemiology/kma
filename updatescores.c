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
#include "pherror.h"
#include "qseqs.h"
#include "updatescores.h"

void update_Scores(unsigned char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, int flag, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = score;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { //Only one best match
		if(*template < 0) {
			template[0] = -template[0];
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter;
		while(i--) {
			alignment_scores[abs(template[i])] += score;
		}
	}
}

void update_Scores_pe(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, int counter, int score, int *start, int *end, int *template, Qseqs *header, Qseqs *header_r, int flag, int flag_r, long unsigned *alignment_scores, long unsigned *uniq_alignment_scores, FILE *frag_out_raw) {
	
	int i, buffer[5];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = -score;
	buffer[3] = header->len;
	buffer[4] = flag;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 5, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	buffer[0] = qr_len;
	buffer[1] = header_r->len;
	buffer[2] = flag_r;
	sfwrite(buffer, sizeof(int), 3, frag_out_raw);
	sfwrite(qseq_r, 1, qr_len, frag_out_raw);
	sfwrite(header_r->seq, 1, header_r->len, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { //Only one best match
		if(*template < 0) {
			template[0] = -template[0];
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter;
		while(i--) {
			alignment_scores[abs(template[i])] += score;
		}
	}
}
