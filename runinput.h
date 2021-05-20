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

#include "compdna.h"
#include "qseqs.h"

/* pointers determining how to deliver the input */
extern void (*printFsa_ptr)(Qseqs*, Qseqs*, CompDNA*, FILE*);
extern void (*printFsa_pair_ptr)(Qseqs*, Qseqs*, Qseqs*, Qseqs*, CompDNA*, FILE*);
long unsigned run_input(char **inputfiles, int fileCount, int minPhred, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, char *trans, const double *prob, FILE *out);
long unsigned run_input_PE(char **inputfiles, int fileCount, int minPhred, int minQ, int fiveClip, int threeClip, int minlen, char *trans, const double *prob, FILE *out);
long unsigned run_input_INT(char **inputfiles, int fileCount, int minPhred, int minQ, int fiveClip, int threeClip, int minlen, char *trans, const double *prob, FILE *out);
void bootFsa(Qseqs *header, Qseqs *qseq, CompDNA *compressor, FILE *out);
void printFsa(Qseqs *header, Qseqs *qseq, CompDNA *compressor, FILE *out);
void printFsa_pair(Qseqs *header, Qseqs *qseq, Qseqs *header_r, Qseqs *qseq_r, CompDNA *compressor, FILE *out);
