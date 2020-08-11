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
#include "assembly.h"
#include "printconsensus.h"

void printConsensus(Assem *aligned_assem, char *header, FILE *alignment_out, FILE *consensus_out, int ref_fsa) {
	
	int i, aln_len, bias;
	
	/* print alignment */
	aln_len = aligned_assem->len;
	if(alignment_out) {
		fprintf(alignment_out, "# %s\n", header);
		for(i = 0; i < aln_len; i += 60) {
			fprintf(alignment_out, "%-10s\t%.60s\n", "template:", aligned_assem->t + i);
			fprintf(alignment_out, "%-10s\t%.60s\n", "", aligned_assem->s + i);
			fprintf(alignment_out, "%-10s\t%.60s\n\n", "query:", aligned_assem->q + i);
		}
	}
	/* Prepare consensus */
	if(ref_fsa == 0) {
		for(i = 0, bias = 0; i < aln_len; ++i, ++bias) {
			aligned_assem->q[bias] = aligned_assem->q[i];
			if(aligned_assem->q[i] == '-') {
				--bias;
			}
		}
		aln_len = bias;
		aligned_assem->q[aln_len] = 0;
	} else if(ref_fsa == 1) {
		for(i = 0; i < aln_len; ++i) {
			if(aligned_assem->q[i] == '-') {
				aligned_assem->q[i] = 'n';
			}
		}
	}
	
	/* Print consensus */
	fprintf(consensus_out, ">%s\n", header);
	for(i = 0; i < aln_len; i += 60) {
		fprintf(consensus_out, "%.60s\n", aligned_assem->q + i);
	}
}
