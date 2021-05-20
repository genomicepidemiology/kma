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
#include <stdlib.h>
#include <stdio.h>
#include "compdna.h"
#include "filebuff.h"
#include "pherror.h"
#include "runinput.h"
#include "qseqs.h"
#include "seqparse.h"
#include "stdstat.h"

void (*printFsa_ptr)(Qseqs*, Qseqs*, CompDNA*, FILE*) = &printFsa;
void (*printFsa_pair_ptr)(Qseqs*, Qseqs*, Qseqs*, Qseqs*, CompDNA*, FILE*) = &printFsa_pair;

long unsigned run_input(char **inputfiles, int fileCount, int minPhred, int minQ, int fiveClip, int threeClip, int minlen, int maxlen, char *trans, const double *prob, FILE *out) {
	
	int fileCounter, phredScale, phredCut, start, end;
	unsigned FASTQ;
	long unsigned count;
	char *filename;
	unsigned char *seq;
	Qseqs *header, *qseq, *qual;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	compressor = smalloc(sizeof(CompDNA));
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	count = 0;
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		}
		
		/* parse the file */
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			phredCut = phredScale + minPhred;
			
			/* parse reads */
			while(FileBuffgetFq(inputfile, header, qseq, qual, trans)) {
				if(qseq->len <= maxlen) {
					/* trim */
					seq = qual->seq;
					start = fiveClip;
					end = qseq->len - 1 - threeClip;
					end = end < 0 ? 0 : end;
					while(end >= 0 && seq[end] < phredCut) {
						--end;
					}
					++end;
					while(start < end && seq[start] < phredCut) {
						++start;
					}
					/*
					for(i = start; i < end; ++i) {
						if(seq[i] < phredCut) {
							seq[i] = 4;
						}
					}
					*/
					qseq->len = end - start;
					/* print */
					if(qseq->len > minlen && minQ <= eQual(seq + start, qseq->len, minQ, prob - phredScale)) {
						/* dump seq */
						qseq->seq += start;
						printFsa_ptr(header, qseq, compressor, out);
						qseq->seq -= start;
						++count;
					}
				}
			}
		} else if(FASTQ & 2) {
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
				if(qseq->len <= maxlen) {
					/* remove leading and trailing N's */
					start = 0;
					end = qseq->len - 1;
					seq = qseq->seq;
					while(end >= 0 && seq[end] == 4) {
						--end;
					}
					++end;
					while(start < end && seq[start] == 4) {
						++start;
					}
					qseq->len = end - start;
					if(qseq->len > minlen) {
						/* dump seq */
						qseq->seq += start;
						printFsa_ptr(header, qseq, compressor, out);
						qseq->seq -= start;
						++count;
					}
				}
			}
		}
		
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
	}
	
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
	
	return count;
}

long unsigned run_input_PE(char **inputfiles, int fileCount, int minPhred, int minQ, int fiveClip, int threeClip, int minlen, char *trans, const double *prob, FILE *out) {
	
	int fileCounter, phredScale, phredCut, start, start2, end;
	unsigned FASTQ, FASTQ2;
	long unsigned count;
	char *filename;
	unsigned char *seq;
	double eq1, eq2;
	Qseqs *header, *qseq, *qual, *header2, *qseq2, *qual2;
	FileBuff *inputfile, *inputfile2;
	CompDNA *compressor;
	
	
	compressor = malloc(sizeof(CompDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	header2 = setQseqs(256);
	qseq2 = setQseqs(1024);
	qual2 = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	inputfile2 = setFileBuff(CHUNK);
	count = 0;
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		
		filename = inputfiles[fileCounter];
		/* determine filetype and open it */
		FASTQ = openAndDetermine(inputfile, filename);
		++fileCounter;
		filename = inputfiles[fileCounter];
		FASTQ2 = openAndDetermine(inputfile2, filename);
		if(FASTQ == FASTQ2) {
			fprintf(stderr, "# Reading inputfile:\t%s %s\n", inputfiles[fileCounter-1], filename);
		} else {
			fprintf(stderr, "Inputfiles:\t%s %s\nAre in different format.\n", inputfiles[fileCounter-1], filename);
			FASTQ = 0;
		}
		
		/* parse the file */
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			if(phredScale == 0) {
				phredScale = getPhredFileBuff(inputfile2);
			}
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			phredCut = phredScale + minPhred;
			
			/* parse reads */
			//while(FileBuffgetFq(inputfile, header, qseq, qual) && FileBuffgetFq(inputfile2, header2, qseq2, qual2)) {
			/* this ensures reading of truncated files */
			while((FileBuffgetFq(inputfile, header, qseq, qual, trans) | FileBuffgetFq(inputfile2, header2, qseq2, qual2, trans))) {
				/* trim forward */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1 - threeClip;
				end = end < 0 ? 0 : end;
				while(end >= 0 && seq[end] < phredCut) {
					--end;
				}
				++end;
				while(start < end && seq[start] < phredCut) {
					++start;
				}
				/*
				for(i = start; i < end; ++i) {
					if(seq[i] < phredCut) {
						seq[i] = 4;
					}
				}
				*/
				qseq->len = end - start;
				eq1 = eQual(seq + start, qseq->len, minQ, prob - phredScale);
				
				/* trim reverse */
				seq = qual2->seq;
				start2 = fiveClip;
				end = qseq2->len - 1;
				while(end >= 0 && seq[end] < phredCut) {
					--end;
				}
				++end;
				while(start2 < end && seq[start2] < phredCut) {
					++start2;
				}
				/*
				for(i = start; i < end; ++i) {
					if(seq[i] < phredCut) {
						seq[i] = 4;
					}
				}
				*/
				qseq2->len = end - start2;
				eq2 = eQual(seq + start2, qseq2->len, minQ, prob - phredScale);
				
				/* print */
				if(qseq->len > minlen && qseq2->len > minlen && minQ <= eq1 && minQ <= eq2) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor, out);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > minlen && minQ <= eq1) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor, out);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > minlen && minQ <= eq2) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor, out);
					qseq2->seq -= start2;
					++count;
				}
			}
		} else if(FASTQ & 2) {
			while((FileBuffgetFsa(inputfile, header, qseq, trans) | FileBuffgetFsa(inputfile2, header2, qseq2, trans))) {
				/* remove leading and trailing N's */
				start = 0;
				end = qseq->len - 1;
				seq = qseq->seq;
				while(end >= 0 && seq[end] == 4) {
					--end;
				}
				++end;
				while(start < end && seq[start] == 4) {
					++start;
				}
				qseq->len = end - start;
				start2 = 0;
				end = qseq2->len - 1;
				seq = qseq2->seq;
				while(end >= 0 && seq[end] == 4) {
					--end;
				}
				++end;
				while(start2 < end && seq[start2] == 4) {
					++start2;
				}
				qseq2->len = end - start2;
				
				/* print */
				if(qseq->len > minlen && qseq2->len > minlen) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor, out);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > minlen) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor, out);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > minlen) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor, out);
					qseq2->seq -= start2;
					++count;
				}
			}
		}
		
		--fileCounter;
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		++fileCounter;
		if(FASTQ2 & 4) {
			gzcloseFileBuff(inputfile2);
		} else {
			closeFileBuff(inputfile2);
		}
		
	}
	
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyQseqs(header2);
	destroyQseqs(qseq2);
	destroyQseqs(qual2);
	destroyFileBuff(inputfile);
	destroyFileBuff(inputfile2);
	
	return count;
}

long unsigned run_input_INT(char **inputfiles, int fileCount, int minPhred, int minQ, int fiveClip, int threeClip, int minlen, char *trans, const double *prob, FILE *out) {
	
	int fileCounter, phredScale, phredCut, start, start2, end;
	unsigned FASTQ;
	long unsigned count;
	char *filename;
	unsigned char *seq;
	double eq1, eq2;
	Qseqs *header, *qseq, *qual, *header2, *qseq2, *qual2;
	FileBuff *inputfile;
	CompDNA *compressor;
	
	
	compressor = malloc(sizeof(CompDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	header2 = setQseqs(256);
	qseq2 = setQseqs(1024);
	qual2 = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	count = 0;
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		}
		
		/* parse the file */
		if(FASTQ & 1) {
			/* get phred scale */
			phredScale = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredScale);
			phredCut = phredScale + minPhred;
			
			/* parse reads */
			while((FileBuffgetFq(inputfile, header, qseq, qual, trans) | FileBuffgetFq(inputfile, header2, qseq2, qual2, trans))) {
				/* trim forward */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1 - threeClip;
				end = end < 0 ? 0 : end;
				while(end >= 0 && seq[end] < phredCut) {
					--end;
				}
				++end;
				while(start < end && seq[start] < phredCut) {
					++start;
				}
				/*
				for(i = start; i < end; ++i) {
					if(seq[i] < phredCut) {
						seq[i] = 4;
					}
				}
				*/
				qseq->len = end - start;
				eq1 = eQual(seq + start, qseq->len, minQ, prob - phredScale);
				
				/* trim reverse */
				seq = qual2->seq;
				start2 = fiveClip;
				end = qseq2->len - 1;
				while(end >= 0 && seq[end] < phredCut) {
					--end;
				}
				++end;
				while(start2 < end && seq[start2] < phredCut) {
					++start2;
				}
				/*
				for(i = start; i < end; ++i) {
					if(seq[i] < phredCut) {
						seq[i] = 4;
					}
				}
				*/
				qseq2->len = end - start2;
				eq2 = eQual(seq + start2, qseq2->len, minQ, prob - phredScale);
				
				/* print */
				if(qseq->len > minlen && qseq2->len > minlen && minQ <= eq1 && minQ <= eq2) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor, out);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > minlen && minQ <= eq1) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor, out);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > minlen && minQ <= eq2) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor, out);
					qseq2->seq -= start2;
					++count;
				}
			}
		} else if(FASTQ & 2) {
			while((FileBuffgetFsa(inputfile, header, qseq, trans) | FileBuffgetFsa(inputfile, header2, qseq2, trans))) {
				/* remove leading and trailing N's */
				start = 0;
				end = qseq->len - 1;
				seq = qseq->seq;
				while(end >= 0 && seq[end] == 4) {
					--end;
				}
				++end;
				while(start < end && seq[start] == 4) {
					++start;
				}
				qseq->len = end - start;
				
				start2 = 0;
				end = qseq2->len - 1;
				seq = qseq2->seq;
				while(end >= 0 && seq[end] == 4) {
					--end;
				}
				++end;
				while(start2 < end && seq[start2] == 4) {
					++start2;
				}
				qseq2->len = end - start2;
				
				/* print */
				if(qseq->len > minlen && qseq2->len > minlen) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor, out);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > minlen) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor, out);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > minlen) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor, out);
					qseq2->seq -= start2;
					++count;
				}
			}
		}
		
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyQseqs(header2);
	destroyQseqs(qseq2);
	destroyQseqs(qual2);
	destroyFileBuff(inputfile);
	
	return count;
}

void bootFsa(Qseqs *header, Qseqs *qseq, CompDNA *compressor, FILE *out) {
	
	int i, end, buffer[4];
	
	/* bootstrap in pieces of 1024 */
	buffer[3] = header->len;
	end = qseq->len - 1024;
	for(i = 0; i < end; i += 512) {
		compDNA(compressor, qseq->seq + i, 1024);
		/* print */
		buffer[0] = compressor->seqlen;
		buffer[1] = compressor->complen;
		buffer[2] = compressor->N[0];
		
		sfwrite(buffer, sizeof(int), 4, out);
		sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
		sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
		sfwrite((header->seq + 1), 1, header->len, out);
		resetComp(compressor);
	}
	
	compDNA(compressor, qseq->seq + i, qseq->len - i);
	/* print */
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header->seq + 1), 1, header->len, out);
	resetComp(compressor);
}

void printFsa(Qseqs *header, Qseqs *qseq, CompDNA *compressor, FILE *out) {
	
	int buffer[4];
	
	/* translate to 2bit */
	if(qseq->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq->len);
	}
	compDNA(compressor, qseq->seq, qseq->len);
	
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	buffer[3] = header->len;
	
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header->seq + 1), 1, header->len, out);
	
	resetComp(compressor);
}

void printFsa_pair(Qseqs *header, Qseqs *qseq, Qseqs *header_r, Qseqs *qseq_r, CompDNA *compressor, FILE *out) {
	
	int buffer[4];
	
	/* translate to 2bit */
	if(qseq->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq->len);
	}
	compDNA(compressor, qseq->seq, qseq->len);
	
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	buffer[3] = -header->len;
	
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header->seq + 1), 1, header->len, out);
	resetComp(compressor);
	
	/* translate to 2bit */
	if(qseq_r->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq_r->len);
	}
	compDNA(compressor, qseq_r->seq, qseq_r->len);
	
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	buffer[3] = header_r->len;
	
	sfwrite(buffer, sizeof(int), 4, out);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, out);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], out);
	sfwrite((header_r->seq + 1), 1, header_r->len, out);
	resetComp(compressor);
}
