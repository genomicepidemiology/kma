/* Philip T.L.C. Clausen Jan 2020 plan@dtu.dk */

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashmapkma.h"
#include "matrix.h"
#include "pherror.h"
#include "runkma.h"

HashMapKMA * loadValues(const char *filename);
void destroyValues(HashMapKMA *src);
void kmerSimilarity(HashMapKMA *DB, Matrix *Dist, int *N);
void kmerDist(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format);
void kmerShared(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format);
void kmerQuery(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format);
void kmerTemplate(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format);
void kmerAvg(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format);
void runDist(char *templatefilename, char *outputfilename, int flag, int format);
int dist_main(int argc, char *argv[]);
