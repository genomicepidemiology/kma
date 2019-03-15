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
#include <math.h>
#include <stdio.h>
#include "assembly.h"
#include "ef.h"
#include "stdnuc.h"
#include "vcf.h"
#include "version.h"

void getExtendedFeatures(char *template_name, AssemInfo *matrix, long unsigned *template_seq, int t_len, Assem *aligned_assem, unsigned fragmentCount, unsigned readCount, FILE *outfile) {
	
	unsigned pos, depthUpdate, maxDepth, nucHighVarSum;
	long unsigned snpSum, insertSum, deletionSum;
	long double var, nucHighVar;
	Assembly *assembly;
	
	if(matrix) {
		/* iterate matrix to get:
			Nuc_high_depth_variance
			Depth_max
			Snp_sum
			Inserts_sum
			Deletions_sum
		*/
		nucHighVar = aligned_assem->depth;
		nucHighVar /= t_len;
		var = aligned_assem->depthVar;
		var /= t_len;
		var -= (nucHighVar * nucHighVar);
		nucHighVar += (3 * sqrt(var));
		
		nucHighVarSum = 0;
		maxDepth = 0;
		snpSum = 0;
		insertSum = 0;
		deletionSum = 0;
		
		assembly = matrix->assmb;
		pos = 0;
		do {
			depthUpdate = assembly[pos].counts[0] + assembly[pos].counts[1] + assembly[pos].counts[2] + assembly[pos].counts[3] + assembly[pos].counts[4];
			
			if(pos < t_len) {
				deletionSum += assembly[pos].counts[5];
				snpSum += (depthUpdate - assembly[pos].counts[getNuc(template_seq, pos)]);
			} else {
				insertSum += depthUpdate;
			}
			
			depthUpdate += assembly[pos].counts[5];
			
			if(maxDepth < depthUpdate) {
				maxDepth = depthUpdate;
			}
			if(nucHighVar < depthUpdate) {
				++nucHighVarSum;
			}
		} while((pos = assembly[pos].next) != 0);
		
		
		fprintf(outfile, "%s\t%u\t%u\t%lu\t%u\t%u\t%lu\t%f\t%u\t%u\t%lu\t%lu\t%lu\n", template_name, readCount, fragmentCount, aligned_assem->score, aligned_assem->aln_len, aligned_assem->cover, aligned_assem->depth, (double) var, nucHighVarSum, maxDepth, snpSum, insertSum, deletionSum);
	} else {
		fprintf(outfile, "%s\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%u\t%u\t%u\t%u\t%u\n", template_name, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0);
	}
}
