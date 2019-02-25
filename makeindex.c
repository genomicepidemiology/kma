/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compdna.h"
#include "filebuff.h"
#include "hashmap.h"
#include "makeindex.h"
#include "pherror.h"
#include "qseqs.h"
#include "seqparse.h"
#include "updateindex.h"

void makeDB(HashMap *templates, int kmerindex, char **inputfiles, int fileCount, char *outputfilename, int appender, char *trans, int MinLen, int MinKlen, double homQ, double homT, unsigned **template_lengths, unsigned **template_ulengths, unsigned **template_slengths) {
	
	int fileCounter, file_len, bias, FASTQ;
	char *filename;
	unsigned char *seq;
	FILE *index_out, *seq_out, *length_out, *name_out;
	Qseqs *header, *qseq;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	/* allocate */
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* open files */
	file_len = strlen(outputfilename);
	strcat(outputfilename, 	".length.b");
	length_out = sfopen(outputfilename, "wb");
	outputfilename[file_len] = 0;
	if(appender) {
		strcat(outputfilename, 	".name");
		name_out = sfopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
	} else {
		strcat(outputfilename, 	".name");
		name_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
	}
	if(appender) {
		strcat(outputfilename, ".seq.b");
		seq_out = sfopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
	} else {
		strcat(outputfilename, ".seq.b");
		seq_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
	}
	
	if(dumpIndex == &makeIndexing) {
		if(appender) {
			strcat(outputfilename, ".index.b");
			index_out = sfopen(outputfilename, "ab");
			outputfilename[file_len] = 0;
		} else {
			strcat(outputfilename, ".index.b");
			index_out = sfopen(outputfilename, "wb");
			outputfilename[file_len] = 0;
			cfwrite(&kmerindex, sizeof(int), 1, index_out);
		}
	} else {
		index_out = 0;
	}
	
	fprintf(stderr, "# Updating DBs\n");
	/* iterate inputfiles */
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		/* open file */
		filename = inputfiles[fileCounter];
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
			
			/* parse the file */
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
				if(qseq->len >= compressor->size) {
					freeComp(compressor);
					allocComp(compressor, qseq->len << 1);
				}
				bias = compDNAref(compressor, qseq->seq, qseq->len);
				if(qseq->len > MinLen && update_DB(templates, compressor, templates->DB_size, MinKlen, homQ, homT, *template_ulengths, *template_slengths)) {
					/* Update annots */
					seq = header->seq + header->len;
					while(isspace(*--seq)) {
						*seq = 0;
					}
					
					if(bias > 0) {
						fprintf(name_out, "%s B%d\n", header->seq + 1, bias);
					} else {
						fprintf(name_out, "%s\n", header->seq + 1);
					}
					updateAnnotsPtr(compressor, templates->DB_size, kmerindex, seq_out, index_out, template_lengths, template_ulengths, template_slengths);
					fprintf(stderr, "# Added:\t%s\n", header->seq + 1);
					
					if(++templates->DB_size == USHRT_MAX) {
						/* convert values to unsigned */
						convertToU(templates);
					}
				} else {
					fprintf(stderr, "# Skipped:\t%s\n", header->seq + 1);
				}
			}
			
			/* close file buffer */
			if(FASTQ & 4) {
				gzcloseFileBuff(inputfile);
			} else {
				closeFileBuff(inputfile);
			}
		}
	}
	
	/* Dump annots */
	cfwrite(&templates->DB_size, sizeof(int), 1, length_out);
	if(*template_ulengths != 0) {
		**template_ulengths = 0;
		**template_slengths = 0;
		cfwrite(*template_lengths, sizeof(unsigned), templates->DB_size, length_out);
		cfwrite(*template_slengths, sizeof(unsigned), templates->DB_size, length_out);
		cfwrite(*template_ulengths, sizeof(unsigned), templates->DB_size, length_out);
	} else {
		**template_lengths = kmerindex;
		cfwrite(*template_lengths, sizeof(unsigned), templates->DB_size, length_out);
	}
	if(index_out) {
		fclose(index_out);
	}
	fclose(seq_out);
	fclose(length_out);
	fclose(name_out);
	
	if(templates->n) {
		fprintf(stderr, "# Templates key-value pairs:\t%lu.\n", templates->n);// / 1048576);
	} else {
		fprintf(stderr, "DB is empty!!!\n");
		exit(1);
	}
	
	/* clean */
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyFileBuff(inputfile);
}
