/* Philip T.L.C. Clausen Jan 2017 s123580@student.dtu.dk */

/*
 Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 All rights reserved.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>

/*
 STRUCTURES
*/
struct aln {
	char *t;  /* template */
	char *s;  /* score */
	char *q;  /* query */
	unsigned pos; /* start of aln, relative to template */
	int score; /* aln score */
	/* start of aln, relative to query */
};

struct assem {
	char *t;  /* template */
	char *s;  /* score */
	char *q;  /* query */
	unsigned cover;
	unsigned depth;
	unsigned size;
};

struct frag {
	char *qseq;
	int score;
	int start;
	int end;
	char *header;
	struct frag *next;
};

struct frag_stats {
	char *seq;
	int *seq_int;
	int stats[4];
};

struct hashTable {
	long unsigned key;
	unsigned *value;
	struct hashTable *next;
};

struct hashMap {
	long unsigned size;
	long unsigned n;
	struct hashTable **table;
};

struct hashTable_index {
	long unsigned key;
	int value;
	struct hashTable_index *next;
};

struct hashMap_index {
	long unsigned size;
	long unsigned n;
	struct hashTable_index **table;
};

union DNA_tree {
	int value;
	union DNA_tree *next;
};

struct Hit {
	unsigned n;
	unsigned tot;
};

struct qseqs {
	int size;
	int len;
	char *seq;
};

struct FileBuff {
	int pos;
	int bytes;
	int buffSize;
	char *buffer;
	FILE *file;
};

struct hashTable_shm {
	long unsigned key;
	unsigned value;
	unsigned next;
};

struct hashMap_shm {
	long unsigned size;
	long unsigned n;
	unsigned *exist; //size long
	struct hashTable_shm *table; //n long
	unsigned *values; // undef long
};

/*
  GLOBAL VARIABLES
*/
struct hashMap *templates;
struct hashMap_shm *templates_shm;
union DNA_tree *templates_align; // hashTree for linking k-mers to posistions
struct hashMap_index **templates_index; // hashMap for linking k-mers to posistions
char **template_seqs; //Whole sequences
char **template_names; //Convertion from int to string name
int *template_lengths, *template_ulengths; // lengths of sequences
double **alignment_scores, evalue, ID_t; // global scores for each tmeplate
int mincoverage, delta, deCon, contamination, print_matrix, ref_fsa, exhaustive, mem_mode, minLen;
unsigned **VF_scores, **VR_scores;
unsigned kmersize, DB_size, shm;
char bases[] = "ATCGN-";
int com_bases[] = {1, 0, 3, 2, 4};
char rcBases[] = "10324";
int *bestTemplates, *bestTemplates_r, *regionTemplates, *Score, *Score_r;
long unsigned *convertNum;
unsigned INITIAL_SIZE = 1048576;
long unsigned MAX_SIZE;
int *input_convertion;
double HMM_param[8];
int **D, **Q, **P, **E;
int W1, U, d[5][5], SW, M, MM, SW_s;
double scoreT = 0;

/* 
	Set function pointers 
*/
struct assem (*assemblyPtr)(int, FILE**, int, FILE*, FILE*, char*, int*, struct aln, struct aln);
int (*ankerPtr)(int*, int*, int, int, int, int, int, int, char*, int);
int (*deConPtr)(int*, int*, int, int, int, int, int, int, char*, int);
void (*aligner)(const int, const char*, int*, struct aln*, struct aln*);
void (*destroyPtr)(int);
void (*alignLoadPtr)(int, FILE*, long unsigned, int);
unsigned * (*hashMap_get)(long unsigned);
void (*kmerScan)(int*, int*, int, char*, int);

/*
 FUNCTIONS
*/

/* BASIC FUNCTIONS */
long unsigned quinaryToDecimal(int *seq, int offset, int key_size) {
	
	int i;
	long unsigned result = 0;
	
	for(i = 0; i < key_size; i++) {
		result += seq[i + offset] * convertNum[i];
	}
	
	return result;
	
}

int chrpos(const char* str1, const char str2) {
	int i, len1;
	
	len1 = strlen(str1);
	if(len1 == 0) {
		return -1;
	}
	
	for(i = 0; i < len1; i++) {
		if(str1[i] == str2)
			return i;
	}
	return -1;
}

int strpos(const char* str1, const char* str2) {
	char* strp;
	int i, len1, len2;
	
	len1 = strlen(str1);
	len2 = strlen(str2);
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	strp = (char*)(str1);
	for(i = 0; i <= len1 - len2; i++) {
		if(*strp == *str2) {
			if(strncmp(strp,str2,len2)==0)
				return i;
		}
		strp++;
	}
	return -1;
}

int strpos_last(const char* str1, const char* str2) {
	char* strp;
	int i, len1, len2;
	
	len1 = strlen(str1);
	len2 = strlen(str2);
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	strp = (char*)(str1 + len1 - len2);
	for(i = len1 - len2; i >= 0; i--) {
		if(*strp == *str2) {
			if(strncmp(strp,str2,len2)==0)
				return i;
		}
		strp--;
	}
	return -1;
}

int uintpos(const unsigned* str1, const unsigned str2) {
	int i, len1;
	
	len1 = str1[0];
	if(len1 == 0) {
		return -1;
	}
	
	for(i = 1; i <= len1; i++) {
		if(str1[i] == str2)
			return i;
	}
	return -1;
}

int uintpos_bin(const unsigned* str1, const unsigned str2) {
	int i, pos, upLim, downLim;
	
	upLim = str1[0];
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) / 2;
	while(upLim - downLim > 0) {
		if(str1[pos] == str2) {
			return pos;
		} else if(str1[pos] < str2) {
			downLim = pos + 1;
		} else {
			upLim = pos - 1;
		}
		pos = (upLim + downLim) / 2;
	}
	if(str1[pos] == str2) {
		return pos;
	}
	return -1;
}

int intpos(const int* str1, const int str2) {
	int i, len1;
	
	len1 = str1[0];
	if(len1 == 0) {
		return -1;
	}
	
	for(i = 1; i <= len1; i++) {
		if(str1[i] == str2)
			return i;
	}
	return -1;
}

int chomp(char *string) {
	/* remove trailing newlines */
	int k = strlen(string) - 1;
	while(string[k] == '\n')
		k--;
	k++;
	string[k] = '\0';
	return k;
}

void insert(char *dest, char src, int location, int dest_len) {
	int i;
	dest[dest_len + 1] = '\0';
	for(i = dest_len; i > location; i--) {
		dest[i] = dest[i - 1];
	}
	dest[location] = src;
}

int replace_chars(char *dest, char src) {
	int i, j, bias, len;
	while(dest[0] == src) {
		dest++;
	}
	len = strlen(dest);
	if(len == 0)
		return len;
	
	bias = 0;
	for(i = 1; i < len && dest[i]; i++) {
		if(dest[i] == src) {
			bias++;
		} else {
			dest[i - bias] = dest[i];
		}
	}
	len -= bias;
	return len;
}

int int_eq(int *s1, int *s2, int len) {
	if(len == 0) {
		return 1;
	}
	int i;
	for(i = 0; i < len; i++) {
		if(s1[i] != s2[i]) {
			return 0;
		}
	}
	return 1;
}

char * fget_line(char *line, int *line_size, int *l_len, FILE *file) {
	
	int i, grow;
	i = 0;
	*l_len = 0;
	grow = *line_size;
	
	while(!feof(file) && fgets((line + i), grow, file) != NULL) {
		*l_len = chomp(line);
		if(*l_len == (*line_size - 1)) { //realloc
			i = *l_len;
			*line_size += grow;
			line = realloc(line, *line_size);
			if(line == NULL) {
				fprintf(stderr, "fget_line error.\n");
				exit(-1);
			}
		} else {
			return line;
		}
	}
	line[0] = '\0';
	return line;
}

char * fget_to(char *line, int *line_size, FILE *file, char stopChar) {
	int c, i, grow;
	
	grow = *line_size;
	i = 0;
	while((c = fgetc(file)) != EOF && (line[i] = c) != stopChar) {
		if((i++) == *line_size) {
			*line_size += grow;
			line = realloc(line, *line_size);
			if(line == NULL) {
				fprintf(stderr, "Memory error\n");
				exit(-1);
			}
		}
	}
	line[i] = '\0';
	
	return line;
}

struct qseqs * setQseqs(int size) {
	
	struct qseqs * dest;
	
	dest = malloc(sizeof(struct qseqs));
	if(!dest) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	dest->len = 0;
	dest->size = size;
	dest->seq = malloc(size);
	if(!dest->seq) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	return dest;
}

void destroyQseqs(struct qseqs *dest) {
	free(dest->seq);
	free(dest);
}

struct FileBuff * setFileBuff(int buffSize) {
	
	struct FileBuff *dest;
	
	dest = malloc(sizeof(struct FileBuff));
	if(!dest) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	dest->pos = 0;
	dest->file = 0;
	dest->buffSize = buffSize;
	dest->buffer = malloc(buffSize);
	if(!dest->buffer) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	return dest;
}

void openFileBuff(struct FileBuff *dest, char *filename, char *mode) {
	
	dest->file = fopen(filename, mode);
	if(!dest->file) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
}

void popenFileBuff(struct FileBuff *dest, char *filename, char *mode) {
	
	dest->file = popen(filename, mode);
	if(!dest->file) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
}

void closeFileBuff(struct FileBuff *dest) {
	fclose(dest->file);
	dest->file = 0;
}

void pcloseFileBuff(struct FileBuff *dest) {
	pclose(dest->file);
	dest->file = 0;
}

void destroyFileBuff(struct FileBuff *dest) {
	free(dest->buffer);
	free(dest);
}

int buffFileBuff(struct FileBuff *dest) {
	dest->pos = 0;
	dest->bytes = fread(dest->buffer, 1, dest->buffSize, dest->file);
	return dest->bytes;
}

int chunkPos(char* seq, int start, int end) {
	
	int i;
	
	for(i = start; i < end; i++) {
		if(seq[i] == '\n') {
			return i;
		}
	}
	
	return end;
}

int FileBuffgetFsa(struct FileBuff *dest, struct qseqs *header, struct qseqs *seq) {
	
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	
	/* test2.b human: 3101810364 */
	int chunk, increase;
	
	/* get header */
	header->len = 0;
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		increase = chunk - dest->pos;
		if(increase + header->len > header->size) {
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(header->seq + header->len, dest->buffer + dest->pos, increase);
		header->len += increase;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	increase = chunk - dest->pos;
	if(increase + header->len > header->size) {
		header->size <<= 1;
		header->seq = realloc(header->seq, header->size);
		if(!header->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	strncpy(header->seq + header->len, dest->buffer + dest->pos, increase);
	header->len += increase;
	header->seq[header->len] = 0;
	dest->pos = chunk + 1;
	
	/* get seq */
	seq->len = 0;
	while(dest->pos < dest->bytes && dest->buffer[dest->pos] != '>' && dest->bytes) {
		while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
			increase = chunk - dest->pos;
			if(increase + seq->len > seq->size) {
				seq->size <<= 1;
				seq->seq = realloc(seq->seq, seq->size);
				if(!seq->seq) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
			}
			strncpy(seq->seq + seq->len, dest->buffer + dest->pos, increase);
			seq->len += increase;
			if(!buffFileBuff(dest)) {
				return 1;
			}
		}
		increase = chunk - dest->pos;
		if(increase + seq->len > seq->size) {
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		
		strncpy(seq->seq + seq->len, dest->buffer + dest->pos, chunk - dest->pos);
		seq->len += (chunk - dest->pos);
		dest->pos = chunk + 1;
	}
	seq->seq[seq->len] = 0;
	
	return 1;
}

int FileBuffgetFsaSeq(struct FileBuff *dest, struct qseqs *seq) {
	
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	
	int chunk, increase;
	
	/* skip header */
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos = chunk + 1;
	
	/* get seq */
	seq->len = 0;
	while(dest->buffer[dest->pos] != '>' && dest->bytes) {
		while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
			increase = chunk - dest->pos;
			if(increase + seq->len > seq->size) {
				seq->size <<= 1;
				seq->seq = realloc(seq->seq, seq->size);
				if(!seq->seq) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
			}
			strncpy(seq->seq + seq->len, dest->buffer + dest->pos, increase);
			seq->len += increase;
			if(!buffFileBuff(dest)) {
				return 1;
			}
		}
		increase = chunk - dest->pos;
		if(increase + seq->len > seq->size) {
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(seq->seq + seq->len, dest->buffer + dest->pos, chunk - dest->pos);
		seq->len += (chunk - dest->pos);
		dest->pos = chunk + 1;
	}
	
	return 1;
}

int FileBuffgetFq(struct FileBuff *dest, struct qseqs *header, struct qseqs *seq, struct qseqs *qual) {
	
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	
	int chunk, increase;
	
	/* get header */
	header->len = 0;
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		increase = chunk - dest->pos;
		if(increase + header->len > header->size) {
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(header->seq + header->len, dest->buffer + dest->pos, increase);
		header->len += increase;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	increase = chunk - dest->pos;
	if(increase + header->len > header->size) {
		header->size <<= 1;
		header->seq = realloc(header->seq, header->size);
		if(!header->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	strncpy(header->seq + header->len, dest->buffer + dest->pos, increase);
	header->len += increase;
	dest->pos = chunk + 1;
	
	/* get seq */
	seq->len = 0;
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		increase = chunk - dest->pos;
		if(increase + seq->len > seq->size) {
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(seq->seq + seq->len, dest->buffer + dest->pos, increase);
		seq->len += increase;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	increase = chunk - dest->pos;
	if(increase + seq->len > seq->size) {
		seq->size <<= 1;
		seq->seq = realloc(seq->seq, seq->size);
		if(!seq->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	strncpy(seq->seq + seq->len, dest->buffer + dest->pos, chunk - dest->pos);
	seq->len += (chunk - dest->pos);
	dest->pos = chunk + 1;
	
	/* skip info */
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos = chunk + 1;
	
	/* get qual */
	if(qual->size < seq->size) {
		free(qual->seq);
		qual->size = seq->size;
		qual->seq = malloc(qual->size);
		if(!qual->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	qual->len = 0;
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		increase = chunk - dest->pos;
		strncpy(qual->seq + qual->len, dest->buffer + dest->pos, increase);
		qual->len += increase;
		if(!buffFileBuff(dest)) {
			return 1;
		}
	}
	increase = chunk - dest->pos;
	strncpy(qual->seq + qual->len, dest->buffer + dest->pos, increase);
	qual->len += increase;
	dest->pos = chunk + 1;
	
	return 1;
}

int FileBuffgetFqSeq(struct FileBuff *dest, struct qseqs *seq, struct qseqs *qual) {
	
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	
	int chunk, increase;
	
	/* skip header */
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos = chunk + 1;
	
	/* get seq */
	seq->len = 0;
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		increase = chunk - dest->pos;
		if(increase + seq->len > seq->size) {
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(seq->seq + seq->len, dest->buffer + dest->pos, increase);
		seq->len += increase;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	increase = chunk - dest->pos;
	if(increase + seq->len > seq->size) {
		seq->size <<= 1;
		seq->seq = realloc(seq->seq, seq->size);
		if(!seq->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	strncpy(seq->seq + seq->len, dest->buffer + dest->pos, chunk - dest->pos);
	seq->len += (chunk - dest->pos);
	dest->pos = chunk + 1;
	
	/* skip info */
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos = chunk + 1;
	
	/* get qual */
	if(qual->size < seq->size) {
		free(qual->seq);
		qual->size = seq->size;
		qual->seq = malloc(qual->size);
		if(!qual->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	qual->len = 0;
	while((chunk = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		increase = chunk - dest->pos;
		strncpy(qual->seq + qual->len, dest->buffer + dest->pos, increase);
		qual->len += increase;
		if(!buffFileBuff(dest)) {
			return 1;
		}
	}
	increase = chunk - dest->pos;
	strncpy(qual->seq + qual->len, dest->buffer + dest->pos, increase);
	qual->len += increase;
	dest->pos = chunk + 1;
	
	return 1;
}

int getPhredFileBuff(struct FileBuff *dest) {
	
	int i, seek;
	
	while(dest->pos < dest->bytes) {
		/* skip header, seq and info */
		for(i = 0; i < 3; i++) {
			seek = 1;
			while(seek) {
				if(dest->pos < dest->bytes) {
					if(dest->buffer[dest->pos] == '\n') {
						seek = 0;
					}
					dest->pos++;
				} else {
					dest->pos = 0;
					return 33;
				}
			}
		}
		/* get Phred scale */
		seek = 1;
		while(seek) {
			if(dest->pos < dest->bytes) {
				if(dest->buffer[dest->pos] == '\n') {
					seek = 0;
				} else if(dest->buffer[dest->pos] < 59) {
					dest->pos = 0;
					return 33;
				} else if(dest->buffer[dest->pos] > 84) {
					dest->pos = 0;
					return 64;
				}
				dest->pos++;
			} else {
				dest->pos = 0;
				return 33;
			}
		}
	}
	
	dest->pos = 0;
	return 33;
}

double fastp(double q) {
	/* P-value from quantile in a chi-square distribution */
	double p = 1.0;
	if(q > 114.5242) {
		p = 1e-26;
	} else if(q > 109.9604) {
		p = 1e-25;
	} else if(q > 105.3969) {
		p = 1e-24;
	} else if(q > 100.8337) {
		p = 1e-23;
	} else if(q > 96.27476) {
		p = 1e-22;
	} else if(q > 91.71701) {
		p = 1e-21;
	} else if(q > 87.16164) {
		p = 1e-20;
	} else if(q > 82.60901) {
		p = 1e-19;
	} else if(q > 78.05917) {
		p = 1e-18;
	} else if(q > 73.51245) {
		p = 1e-17;
	} else if(q > 68.96954) {
		p = 1e-16;
	} else if(q > 64.43048) {
		p = 1e-15;
	} else if(q > 59.89615) {
		p = 1e-14;
	} else if(q > 55.36699) {
		p = 1e-13;
	} else if(q > 50.84417) {
		p = 1e-12;
	} else if(q > 46.32844) {
		p = 1e-11;
	} else if(q > 41.82144) {
		p = 1e-10;
	} else if(q > 37.32489) {
		p = 1e-9;
	} else if(q > 32.84127) {
		p = 1e-8;
	} else if(q > 28.37395) {
		p = 1e-7;
	} else if(q > 23.92814) {
		p = 1e-6;
	} else if(q > 19.51139) {
		p = 1e-5;
	} else if(q > 15.13671) {
		p = 1e-4;
	} else if(q > 10.82759) {
		p = 1e-3;
	} else if(q > 6.634897) {
		p = 0.01;
	} else if(q > 3.841443) {
		p = 0.05;
	} else if(q > 2.705532) {
		p = 0.1;
	} else if(q > 2.072251) {
		p = 0.15;
	} else if(q > 1.642374) {
		p = 0.2;
	} else if(q > 1.323304) {
		p = 0.25;
	} else if(q > 1.074194) {
		p = 0.3;
	} else if(q > 0.8734571) {
		p = 0.35;
	} else if(q > 0.7083263) {
		p = 0.4;
	} else if(q > 0.5706519) {
		p = 0.45;
	} else if(q > 0.4549364) {
		p = 0.5;
	} else if(q > 0.3573172) {
		p = 0.55;
	} else if(q > 0.2749959) {
		p = 0.6;
	} else if(q > 0.2059001) {
		p = 0.65;
	} else if(q > 0.1484719) {
		p = 0.7;
	} else if(q > 0.1015310) {
		p = 0.75;
	} else if(q > 0.06418475) {
		p = 0.8;
	} else if(q > 0.03576578) {
		p = 0.85;
	} else if(q > 0.01579077) {
		p = 0.9;
	} else if(q > 0.00393214) {
		p = 0.95;
	} else if(q >= 0.0) {
		p = 1.0;
	} else {
		p = 1.00 - fastp(-1 * q);
	}
	return p;
}

double p_chisqr(double q) {
	if(q > 49) {
		/* Get p-val from table */
		return fastp(q);
	}
	return 1 - 1.772453850 * erf(sqrt(0.5 * q)) / tgamma(0.5);
}

/* DNA SPECIFIC FUNCTIONS */
int find_contamination(int *out_Tem) {
	int i;
	for(i = out_Tem[0]; i > 0; i--) {
		if(out_Tem[i] == contamination) {
			return i;
		}
	}
	return -1;
}

int find_contamination2(int *out_Tem, int contamination_s) {
	int i;
	for(i = out_Tem[0]; i > 0; i--) {
		if(out_Tem[i] == contamination_s) {
			return i;
		}
	}
	return -1;
}

void convert_DNA(char *qseq, int *qseq_int, int seqlen) {
	int i;
	for(i = 0; i < seqlen; i++) {
		qseq_int[i] = chrpos(bases, qseq[i]);
	}
}

void decode_DNA(char *qseq, int *qseq_int, int seqlen) {
	int i;
	for(i = 0; i < seqlen; i++) {
		qseq_int[i] = input_convertion[qseq[i]];
	}
}

void strrc(char *qseq, int q_len) {
        
        int i, seqlen;
        char carry;
        
        seqlen = q_len >> 1;
        
        for(i = 0; i < seqlen; i++) {
                carry = rcBases[qseq[i] - '0'];
                qseq[i] = rcBases[qseq[q_len - i - 1] - '0'];
                qseq[q_len - i - 1] = carry;
        }
        if(q_len & 1) {
                qseq[seqlen] = rcBases[qseq[seqlen] - '0'];
        }
        
}

/* DNA TREE */
void set_node(union DNA_tree *dest) {
	dest->next = NULL;
}

void tree_addKey(union DNA_tree *dest, int *key, int offset, int key_size, int value) {
	
	int i, j;
	union DNA_tree *node = dest;
	for(i = offset; i < offset + key_size; i++) {
		/* Create node if not there */
		if(node->next == NULL) {
			node->next = malloc(5 * sizeof(union DNA_tree));
			if(!node->next) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			for(j = 0; j < 5; j++) {
				set_node(&(node->next[j]));
			}
		}
		node = &(node->next[key[i]]);
	}
	
	if(node->next == NULL) {
		node->value = value;
	} else {
		node->value = -2;
	}
	
}

int tree_addValue(union DNA_tree *dest, int *key, int offset, int key_size, int value) {
	
	int i, j;
	union DNA_tree *node = dest;
	for(i = offset; i < offset + key_size; i++) {
		/* Create node if not there */
		if(node->next == NULL) {
			node->next = malloc(5 * sizeof(union DNA_tree));
			if(!node->next) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			for(j = 0; j < 5; j++) {
				set_node(&(node->next[j]));
			}
		}
		node = &(node->next[key[i]]);
	}
	
	if(node->next == NULL) {
		node->value = value;
		return 1;
	} else {
		node->value = -2;
		return 0;
	}
	
}

void tree_addPath(union DNA_tree *dest, int *key, int offset, int key_size, int value) {
	
	int i, j;
	union DNA_tree *node = dest;
	for(i = offset; i < offset + key_size; i++) {
		/* Create node if not there */
		if(node->next == NULL) {
			node->next = malloc(5 * sizeof(union DNA_tree));
			if(!node->next) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			for(j = 0; j < 5; j++) {
				set_node(&(node->next[j]));
			}
		}
		node = &(node->next[key[i]]);
	}
}

int tree_getValue(union DNA_tree *dest, int* key, int offset, int key_size) {
	
	int i;
	union DNA_tree *node = (union DNA_tree*)(dest);
	
	for(i = offset; i < offset + key_size; i++) {
		if(node->next == NULL) {
			return (-1) * (i - offset);
		}
		node = (union DNA_tree*)(&(node->next[key[i]]));
	}
	
	if(node->next == NULL) {
		return (-1) * (i - offset);
	} else if(node->value == -2) {
		return (-1) * (key_size + 1);
	} else {
		return node->value;
	}
}

int tree_hasValue(union DNA_tree *dest, char* key, int offset, int key_size) {
	
	int i;
	union DNA_tree *node = (union DNA_tree*)(dest);
	
	for(i = offset; i < offset + key_size; i++) {
		if(node->next == NULL) {
			return 0;
		}
		node = (union DNA_tree*)(&(node->next[key[i] - '0']));
	}
	return 1;
}

void dump_tree(union DNA_tree *dest, FILE *file, int value) {
	/* Save the hash to a file */
	int i;
	union DNA_tree *node;
	node = (union DNA_tree*)(dest);
	
	/* Dump node */
	fwrite(node, sizeof(union DNA_tree), 1, file);
	
	/* Recursive dump for next nodes, pointed to from this node */
	if(value < kmersize && node->next) {
		value++;
		for(i = 0; i < 5; i++) {
			dump_tree((union DNA_tree*)(&(node->next[i])), file, value);
		}
	}
}

void load_tree(union DNA_tree *dest, FILE *file, int value) {
	/* Load hash from file, opposite of dump_hash */
	int i;
	union DNA_tree *node;
	node = (union DNA_tree*)(dest);
	
	/* Load node */
	fread(node, sizeof(union DNA_tree), 1, file);
	
	/* Recursive load next nodes, pointed to from this node */
	if(value < kmersize && node->next) {
		value++;
		node->next = calloc(5, sizeof(union DNA_tree));
		if(!node->next) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(i = 0; i < 5; i++) {
			load_tree((union DNA_tree*)(&(node->next[i])), file, value);
		}
	}
}

void dump_values(int **values, FILE *file) {
	int i;
	fwrite(&(int){values[0][0]}, sizeof(int), 1, file);
	for(i = 1; i <= values[0][0]; i++) {
		if(values[i] != NULL) {
			fwrite(&(int){values[i][0]}, sizeof(int), 1, file);
			fwrite(values[i], (values[i][0] + 1) * sizeof(int), 1, file);
		} else {
			fprintf(stderr, "# NULL at %d\n", i);
			exit(-1);
		}
	}
}

void destroy_tree(union DNA_tree *dest, int value) {
	/* Destroy all nodes and leaves of a tree */
	int i;
	union DNA_tree *node;
	node = (union DNA_tree*)(dest);
	/* Recursive destroy for next nodes, pointed to from this node */
	if(value < kmersize && node->next) {
		value++;
		for(i = 0; i < 5; i++) {
			destroy_tree((union DNA_tree*)(&(node->next[i])), value);
		}
		/* destroy this node */
		free(node->next);
	}
	
}

int ** load_values(FILE *file) {
	int **values;
	int i, size;
	fread(&size, sizeof(int), 1, file);
	values = malloc((size + 1) * sizeof(int*));
	if(!values) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	values[0] = malloc(2 * sizeof(int));
	if(!values[0]) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	values[0][0] = size;
	for(i = 1; i <= values[0][0]; i++) {
		fread(&size, sizeof(int), 1, file);
		values[i] = malloc((size + 1) * sizeof(int));
		if(!values[i]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		fread(values[i], (size + 1) * sizeof(int), 1, file);
	}
	return values;
}

void alignDump(char *qseq, int q_len, union DNA_tree *align, int template, const char *align_path) {
	
	int out_len;
	char *out_name;
	FILE *align_out;
	/* open outputfile */
	out_len = strlen(align_path);
	out_name = calloc(out_len + 256, sizeof(char));
	if(!out_name) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	sprintf(out_name, "%s%d.b", align_path, template);
	align_out = fopen(out_name, "wb");
	if(!align_out) {
		fprintf(stderr, "File coruption: %s\n", out_name);
		exit(1);
	}
	/* dump seq */
	fwrite(qseq, q_len * sizeof(char), 1, align_out);
	/* dump index */
	dump_tree(align, align_out, 0);
	
	/* close up */
	fclose(align_out);
	destroy_tree(&templates_align[DB_size], 0);
}

/* HASHMAP FUNCTIONS */
void initialize_hashMap(struct hashMap *dest, unsigned newSize) {
	/* set hashMap */
	dest->size = newSize;
	dest->n = 0;
	/* set hashTable */
	dest->table = calloc(newSize, sizeof(struct hashTable*));
	if(!dest->table) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
}

void set_hashTable(struct hashTable *dest) {
	dest->value = 0;
	dest->next = 0;
}

void hashMap_destroy(struct hashMap *dest) {
	
	int i;
	struct hashTable *node, *next;
	
	/* destroy hashtable elements*/
	for(i = 0; i < dest->size; i++) {
		if(dest->table[i] != 0) {
			for(node = dest->table[i]; node != 0; node = next) {
				next = node->next;
				free(node->value);
				free(node);
			}
			dest->table[i] = 0;
		}
	}
	free(dest->table);
	dest->table = 0;
}

void hashMap_addKey(struct hashMap *dest, long unsigned key, unsigned value) {
	unsigned i, j, index;
	struct hashTable *node;
	
	/* check if allocation is needed */
	if(dest->n + 1 >= dest->size && dest->size < MAX_SIZE) {
		/* copy content of the hashTable */
		struct hashTable *carryOver;
		struct hashTable *next;
		unsigned old_n = dest->n;
		carryOver = 0;
		for(i = 0; i < dest->size; i++) {
			for(node = dest->table[i]; node != 0; node = next) {
				next = node->next;
				node->next = carryOver;
				carryOver = node;
			}
		}
		
		free(dest->table);
		
		/* allocate hash size */
		if(dest->size + INITIAL_SIZE > MAX_SIZE) {
			initialize_hashMap(dest, MAX_SIZE);
		} else {
			initialize_hashMap(dest, dest->size + INITIAL_SIZE);
		}
		dest->n = old_n;
		
		/* refill hash */
		for(node = carryOver; node != 0; node = next) {
			next = node->next;
			index = node->key % dest->size;
			node->next = dest->table[index];
			dest->table[index] = node;
			
		}
	}
	
	/* get index */
	index = key % dest->size;
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		dest->n++;
		dest->table[index] = malloc(sizeof(struct hashTable));
		if(!dest->table[index]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node = dest->table[index];
		node->value = calloc(2, sizeof(unsigned));
		if(!node->value) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node->value[0] = 1;
		node->value[1] = value;
		node->key = key;
		node->next = 0;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(node->value == 0) { // New value
				dest->n++;
				node->value = calloc(2, sizeof(unsigned));
				if(node->value == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				node->value[0] = 1;
				node->value[1] = value;
				node->key = key;
				return;
			} else if(key == node->key) { // Keys match change value
				node->value[0]++;
				node->value = realloc(node->value, (node->value[0] + 1) * sizeof(unsigned));
				if(node->value == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				node->value[node->value[0]] = value;
				return;
			}
			if(node->next == 0) { // This chain is filled, create next
				node->next = malloc(sizeof(struct hashTable));
				if(node->next == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				set_hashTable(node->next);
			}
		}
	}
}

void hashMap_CountUp(struct hashMap *dest, long unsigned key, unsigned offset) {
	unsigned i, j, index;
	struct hashTable *node;
	
	/* check if allocation is needed */
	if(dest->n + 1 >= dest->size && dest->size < MAX_SIZE) {
		/* copy content of the hashTable */
		struct hashTable *carryOver;
		struct hashTable *next;
		unsigned old_n = dest->n;
		carryOver = 0;
		for(i = 0; i < dest->size; i++) {
			for(node = dest->table[i]; node != 0; node = next) {
				next = node->next;
				node->next = carryOver;
				carryOver = node;
			}
		}
		
		free(dest->table);
		
		/* allocate hash size */
		if(dest->size + INITIAL_SIZE > MAX_SIZE) {
			initialize_hashMap(dest, MAX_SIZE);
		} else {
			initialize_hashMap(dest, dest->size + INITIAL_SIZE);
		}
		dest->n = old_n;
		
		/* refill hash */
		for(node = carryOver; node != 0; node = next) {
			next = node->next;
			index = node->key % dest->size;
			node->next = dest->table[index];
			dest->table[index] = node;
			
		}
	}
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		dest->n++;
		dest->table[index] = malloc(sizeof(struct hashTable));
		if(!dest->table[index]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node = dest->table[index];
		node->value = malloc(sizeof(unsigned));
		if(!node->value) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node->value[0] = offset;
		node->key = key;
		node->next = 0;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(key == node->key) { // Keys match change value
				node->value[0]++;
				return;
			} else if(node->next == 0) { // This chain is filled, create next
				dest->n++;
				node->next = malloc(sizeof(struct hashTable));
				if(!node->next) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				node = node->next;
				node->next = 0;
				node->key = key;
				node->value = malloc(sizeof(unsigned));
				if(!node->value) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				node->value[0] = offset;
				return;
			}
		}
	}
}

unsigned * hashMap_getValue(struct hashMap *dest, long unsigned key) {
	unsigned index;
	struct hashTable *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	for(node = dest->table[index]; node != 0; node = node->next) {
		if(node->key == key) {
			/* Keys match, we found it */
			return node->value;
		}
	}
	return 0;
}

unsigned * hashMap_getGlobal(long unsigned key) {
	unsigned index;
	struct hashTable *node;
	
	/* get index */
	index = key % templates->size;
	
	/* find pos */
	for(node = templates->table[index]; node != 0; node = node->next) {
		if(node->key == key) {
			/* Keys match, we found it */
			return node->value;
		}
	}
	return 0;
}

void hashMap_del(struct hashMap *dest, long unsigned key) {
	unsigned index;
	struct hashTable *node, *prev;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	prev = 0;
	node = dest->table[index];
	while(node != 0) {
		if(node->key == key) {
			/* Keys match, we found it */
			if(prev == 0) {
				dest->table[index] = node->next;
				free(node->value);
				free(node);
				return;
			} else {
				prev->next = node->next;
				free(node->value);
				free(node);
				return;
			}
		}
		prev = node;
		node = node->next;
	}
}

void hashMap_dump(struct hashMap *dest, FILE *file) {
	
	int i;
	struct hashTable *node;
	
	/* dump hashMap */
	fwrite(&kmersize, sizeof(unsigned), 1, file);
	fwrite(&dest->size, sizeof(long unsigned), 1, file);
	fwrite(&dest->n, sizeof(long unsigned), 1, file);
	
	/* dump hashtable */
	for(i = 0; i < dest->size; i++) {
		if(dest->table[i] != 0) {
			for(node = dest->table[i]; node != 0; node = node->next) {
				fwrite(&(node->key), sizeof(long unsigned), 1, file);
				fwrite(&(unsigned){node->value[0] + 1}, sizeof(unsigned), 1, file);
				fwrite(node->value, (node->value[0] + 1) * sizeof(unsigned), 1, file);
			}
		}
	}
}

void hashMap_load(struct hashMap *dest, FILE *file) {
	
	unsigned i, index, tmp_size;
	long unsigned key;
	struct hashTable *node;
	
	/* load hashMap */
	fread(&kmersize, sizeof(unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	if(dest->size > MAX_SIZE) {
		dest->size = MAX_SIZE;
	}
	dest->table = calloc(dest->size, sizeof(struct hashTable*));
	if(dest->table == NULL) {
		fprintf(stderr, "# Failed to allocate k-mer DB\n");
		fprintf(stderr, "# This will hurt performance, but i will try make the hashing list smaller.\n");
		while(dest->table == NULL) {
			dest->size -= INITIAL_SIZE;
			if(dest->size < INITIAL_SIZE) {
				fprintf(stderr, "# It is impossible to store the DB!!!\n");
				exit(-1);
			}
			dest->table = calloc(dest->size, sizeof(struct hashTable*));
		}
	}
	/* load hashTable */
	for(i = 0; i < dest->n; i++) {
		fread(&key, sizeof(long unsigned), 1, file);
		index = key % dest->size;
		if(dest->table[index] == 0) {
			dest->table[index] = malloc(sizeof(struct hashTable));
			if(!dest->table[index]) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}			
			node = dest->table[index];
			node->key = key;
			fread(&tmp_size, sizeof(unsigned), 1, file);
			node->value = malloc(tmp_size * sizeof(unsigned));
			fread(node->value, tmp_size * sizeof(unsigned), 1, file);
			node->next = 0;
		} else {
			node = malloc(sizeof(struct hashTable));
			if(!node) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}			
			node->key = key;
			fread(&tmp_size, sizeof(unsigned), 1, file);
			node->value = malloc(tmp_size * sizeof(unsigned));
			if(!node->value) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}			
			fread(node->value, tmp_size * sizeof(unsigned), 1, file);
			node->next = dest->table[index];
			dest->table[index] = node;
		}
	}
}

/* hashMap indexes */

void initialize_hashMap_index(struct hashMap_index *dest, unsigned newSize) {
	/* set hashMap */
	dest->size = newSize;
	dest->n = 0;
	/* set hashTable */
	dest->table = calloc(newSize, sizeof(struct hashTable_index*));
	if(!dest->table) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
}

void hashMap_destroy_index(struct hashMap_index *dest) {
	
	int i;
	struct hashTable_index *node, *next;
	
	/* destroy hashtable elements*/
	for(i = 0; i < dest->size; i++) {
		if(dest->table[i] != 0) {
			for(node = dest->table[i]; node != 0; node = next) {
				next = node->next;
				free(node);
			}
			dest->table[i] = 0;
		}
	}
	free(dest->table);
	dest->table = 0;
}

void hashMap_addIndex(struct hashMap_index *dest, long unsigned key, int indexNum) {
	unsigned i, j, index;
	struct hashTable_index *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		dest->n++;
		dest->table[index] = malloc(sizeof(struct hashTable_index));
		if(!dest->table[index]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node = dest->table[index];
		node->value = indexNum;
		node->key = key;
		node->next = 0;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(key == node->key) { // Keys match change value
				node->value = -2;
				return;
			} else if(node->next == 0) { // This chain is filled, create next
				dest->n++;
				node->next = malloc(sizeof(struct hashTable));
				if(!node->next) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}			
				node = node->next;
				node->next = 0;
				node->key = key;
				node->value = indexNum;
				return;
			}
		}
	}
}

void hashMap_CountIndex(struct hashMap_index *dest, long unsigned key) {
	unsigned i, j, index;
	struct hashTable_index *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		dest->n++;
		dest->table[index] = malloc(sizeof(struct hashTable_index));
		if(!dest->table[index]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}					
		node = dest->table[index];
		node->value = 1;
		node->key = key;
		node->next = 0;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(key == node->key) { // Keys match change value
				node->value++;
				return;
			} else if(node->next == 0) { // This chain is filled, create next
				dest->n++;
				node->next = malloc(sizeof(struct hashTable));
				if(!node->next) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}			
				node = node->next;
				node->next = 0;
				node->key = key;
				node->value = 1;
				return;
			}
		}
	}
}

int hashMap_getIndex(struct hashMap_index *dest, long unsigned key) {
	unsigned index;
	struct hashTable_index *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	for(node = dest->table[index]; node != 0; node = node->next) {
		if(node->key == key) {
			/* Keys match, we found it */
			return node->value;
		}
	}
	return -1;
}

/* HASHMAP_SHM FUNCTIONS */

void hashMap_shm_load(struct hashMap_shm *dest, FILE* file, char *filename) {
	
	key_t key;
	int id, i;
	unsigned valueSize;
	
	if(file == NULL) {
		fprintf(stderr, "No such shared DB.\n");
		exit(-1);
	}
	
	fread(&contamination, sizeof(unsigned), 1, file);
	fread(&kmersize, sizeof(unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	fread(&valueSize, sizeof(unsigned), 1, file);
	
	/* get tables */
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, dest->size * sizeof(unsigned), 0666);
	if(id < 0) {
		key = ftok(filename, 'e');
		id = shmget(key, dest->size * sizeof(unsigned), 0666);
		if(id < 0) {
			fprintf(stderr, "Sharing not available for this DB.\n");
			exit(-1);
		}
	}
	
	dest->exist = (unsigned *) shmat(id, NULL, 0);
	
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, dest->n * sizeof(struct hashTable_shm), 0666);
	if(id < 0) {
		key = ftok(filename, 't');
		id = shmget(key, dest->n * sizeof(struct hashTable_shm), 0666);
		if(id < 0) {
			fprintf(stderr, "Sharing not available for this DB.\n");
			exit(-1);
		}
	}
	dest->table = (struct hashTable_shm *) shmat(id, NULL, 0);
	
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, valueSize * sizeof(unsigned), 0666);
	if(id < 0) {
		key = ftok(filename, 'v');
		id = shmget(key, valueSize * sizeof(unsigned), 0666);
		if(id < 0) {
			fprintf(stderr, "Sharing not available for this DB.\n");
			exit(-1);
		}
	}
	dest->values = (unsigned *) shmat(id, NULL, 0);
	
}

unsigned * hashMap_shm_getValue(struct hashMap_shm *dest, long unsigned key) {
	
	unsigned index, node_index;
	
	/* get index */
	index = key % dest->size;
	
	if(dest->exist[index] == dest->n + 1) {
		return 0;
	}
	
	/* link through the list */
	struct hashTable_shm *table = dest->table;
	node_index = dest->exist[index];
	while(table[node_index].next != 0) {
		if(table[node_index].key == key) {
			return (dest->values + table[node_index].value);
		}
		node_index = table[node_index].next;
	}
	
	if(table[node_index].key == key) {
		return (dest->values + table[node_index].value);
	}
	
	return 0;
}

unsigned * hashMap_shm_getGlobal(long unsigned key) {
	
	unsigned index, node_index;
	struct hashMap_shm *dest = templates_shm;
	
	/* get index */
	index = key % dest->size;
	
	if(dest->exist[index] == dest->n + 1) {
		return 0;
	}
	
	/* link through the list */
	struct hashTable_shm *table = dest->table;
	node_index = dest->exist[index];
	while(table[node_index].next != 0) {
		if(table[node_index].key == key) {
			return (dest->values + table[node_index].value);
		}
		node_index = table[node_index].next;
	}
	
	if(table[node_index].key == key) {
		return (dest->values + table[node_index].value);
	}
	
	return 0;
}

void hashMap_shm_detach(struct hashMap_shm *dest) {
	shmdt(dest->exist);
	shmdt(dest->table);
	shmdt(dest->values);
}

/* DB LOADING */
void alignLoad_fly(int template, FILE *align_in, long unsigned file_index, int SET) {
	
	int i, end, file_len, *qseq_int;
	/* set file pointer */
	fseek(align_in, file_index, SET);
	
	/* allocate seq */
	template_seqs[template] = malloc((template_lengths[template] + 1) * sizeof(char));
	qseq_int = malloc(template_lengths[template] * sizeof(int));
	if(!qseq_int || !template_seqs[template]) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* read seq */
	fread(template_seqs[template], template_lengths[template] * sizeof(char), 1, align_in);
	template_seqs[template][template_lengths[template]] = '\0';
	
	/* create index */
	for(i = 0; i < template_lengths[template]; i++) {
		qseq_int[i] = template_seqs[template][i] - '0';
	}
	end = template_lengths[template] - kmersize + 1;
	for(i = 0; i < end; i++) {
		tree_addValue(&templates_align[template], qseq_int, i, kmersize, i);
	}
	for(i = end; i < template_lengths[template]; i++) {
		tree_addPath(&templates_align[template], qseq_int, i, template_lengths[template] - i, i);
	}
	
	free(qseq_int);
}

void alignLoad_fly_SW(int template, FILE *align_in, long unsigned file_index, int SET) {
	
	int i, end, file_len, *qseq_int;
	long unsigned key;
	/* set file pointer */
	fseek(align_in, file_index, SET);
	
	/* allocate seq */
	template_seqs[template] = malloc((template_lengths[template] + 1) * sizeof(char));
	qseq_int = malloc(template_lengths[template] * sizeof(int));
	if(!qseq_int || !template_seqs[template]) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* read seq */
	fread(template_seqs[template], template_lengths[template] * sizeof(char), 1, align_in);
	template_seqs[template][template_lengths[template]] = '\0';
	for(i = 0; i < template_lengths[template]; i++) {
		qseq_int[i] = template_seqs[template][i] - '0';
	}
	
	/* allocate hashMap */
	end = template_lengths[template] - kmersize + 1;
	templates_index[template] = malloc(sizeof(struct hashMap_index));
	if(!templates_index[template]) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	initialize_hashMap_index(templates_index[template], end);
	
	/* create indexes */
	key = quinaryToDecimal(qseq_int, 0, kmersize);
	hashMap_addIndex(templates_index[template], key, 0);
	for(i = 1; i < end; i++) {
		key = (key - qseq_int[i - 1] * convertNum[0]) * 5 + qseq_int[i + kmersize - 1];
		hashMap_addIndex(templates_index[template], key, i);
	}
	
	free(qseq_int);
}

void alignClean(int template) {
	destroy_tree(&templates_align[template], 0);
}

void alignClean_SW(int template) {
	hashMap_destroy_index(templates_index[template]);
	free(templates_index[template]);
}

void load_DBs_KMA(char *filename) {
	/* load DBs needed for KMA */
	int i, file_len;
	char *templatefilename;
	FILE *DB_file;
	
	/* reallocate filename, to get all DBs */
	templatefilename = strdup(filename);
	if(templatefilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	file_len = strlen(templatefilename);
	templatefilename = realloc(templatefilename, (file_len + 10) * sizeof(char));
	if(templatefilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	/* allocate DBs */
	strcat(templatefilename, ".length.b");
	templatefilename[file_len + 9] = '\0';
	DB_file = fopen(templatefilename, "rb");
	if(DB_file == NULL) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
		exit(-1);
	}
	templatefilename[file_len] = '\0';
	fread(&DB_size, sizeof(int), 1, DB_file);
	template_lengths = malloc(DB_size * sizeof(int));
	template_seqs = calloc(DB_size, sizeof(char*));
	template_names = malloc(DB_size * sizeof(char*));
	if(SW) {
		templates_index = calloc(DB_size, sizeof(struct hashMap_index*));
		if(!templates_index) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
	} else {
		templates_align = malloc(DB_size * sizeof(union DNA_tree));
		if(!templates_align) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
	}
	alignment_scores = malloc(DB_size * sizeof(double*));
	if(!template_lengths || !template_seqs || !template_names || !alignment_scores) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	contamination = DB_size + 1;
	for(i = 0; i < DB_size; i++) {
		alignment_scores[i] = calloc(2, sizeof(double));
		if(!alignment_scores[i]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
	}
	/* load lengths */
	fread(template_lengths, DB_size * sizeof(int), 1, DB_file);
	fclose(DB_file);
	
	/* load names */
	strcat(templatefilename, ".name.b");
	templatefilename[file_len + 7] = '\0';
	DB_file = fopen(templatefilename, "rb");
	if(DB_file == NULL) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
		exit(-1);
	}
	templatefilename[file_len] = '\0';
	for(i = 0; i < DB_size; i++) {
		template_names[i] = malloc(256 * sizeof(char));
		if(!template_names[i]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		fread(template_names[i], 256 * sizeof(char), 1, DB_file);
		template_names[i] = realloc(template_names[i], (strlen(template_names[i]) + 1) * sizeof(char));
		if(!template_names[i]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		template_names[i][strlen(template_names[i])] = '\0';
	}
	fclose(DB_file);
	
	free(templatefilename);
}

void load_DBs_Sparse(const char *filename) {
	/* load DBs needed for KMA */
	int i, file_len;
	char *templatefilename;
	FILE *DB_file;
	
	/* reallocate filename, to get all DBs */
	templatefilename = strdup(filename);
	if(templatefilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	file_len = strlen(templatefilename);
	templatefilename = realloc(templatefilename, (file_len + 10) * sizeof(char));
	if(templatefilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* Open DB */
	strcat(templatefilename, ".spaLen.b");
	templatefilename[file_len + 9] = '\0';
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "File coruption: %s\n", templatefilename);
		exit(1);
	}
	templatefilename[file_len] = '\0';
	
	/* allocate DBs */
	fread(&DB_size, sizeof(int), 1, DB_file);
	template_lengths = malloc(DB_size * sizeof(int));
	template_ulengths = malloc(DB_size * sizeof(int));
	template_names = malloc(DB_size * sizeof(char*));
	if(!template_lengths || !template_ulengths || !template_names) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	/* load lengths */
	fread(template_lengths, DB_size * sizeof(int), 1, DB_file);
	fread(template_ulengths, DB_size * sizeof(int), 1, DB_file);
	fclose(DB_file);
	
	/* load names */
	strcat(templatefilename, ".name.b");
	templatefilename[file_len + 7] = '\0';
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "File coruption: %s\n", templatefilename);
		exit(1);
	}
	templatefilename[file_len] = '\0';
	for(i = 0; i < DB_size; i++) {
		template_names[i] = malloc(256 * sizeof(char));
		if(template_names[i] == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		fread(template_names[i], 256 * sizeof(char), 1, DB_file);
		template_names[i] = realloc(template_names[i], (strlen(template_names[i]) + 1) * sizeof(char));
		if(template_names[i] == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		template_names[i][strlen(template_names[i])] = '\0';
	}
	
	fclose(DB_file);
	free(templatefilename);
}

/* METHOD SPECIFIC METHODS */
void print_ankers(int *out_Tem, int *qseq_int, int seqlen, int rc_flag, char *header, int header_len) {
	
	int i, contPos;
	/* check for contamination */
	fwrite(&(int){seqlen}, sizeof(int), 1, stdout);
	fwrite(qseq_int, seqlen * sizeof(int), 1, stdout);
	fwrite(&(int){rc_flag}, sizeof(int), 1, stdout);
	fwrite(&(int){out_Tem[0]}, sizeof(int), 1, stdout);
	fwrite(out_Tem, (out_Tem[0] + 1) * sizeof(int), 1, stdout);
	fwrite(&(int){header_len}, sizeof(int), 1, stdout);
	fwrite(header, header_len, 1, stdout);
}

int ankerAndClean(int *regionTemplates, int *qseq_int, int start, int stop, int HIT, int bestScore, int start_cut, int end_cut, char *header, int header_len) {
	
	int k, l, bestHitsCov, template;
	double thisCov, bestCov;
	unsigned *values;
	
	bestHitsCov = regionTemplates[0];
	bestCov = 0;
	for(k = 1; k <= regionTemplates[0]; k++) {
		template = regionTemplates[k];
		if(template < 0) {
			template *= (-1);
			thisCov = 1.0 * Score_r[template] / template_lengths[template];
		} else {
			thisCov = 1.0 * Score[template] / template_lengths[template];
		}
		if(thisCov > bestCov) {
			bestCov = thisCov;
		}
	}
	for(k = start_cut; k <= end_cut; k++) {
		if(VF_scores[k] != 0) {
			values = VF_scores[k];
			for(l = 1; l <= values[0]; l++) {
				if(Score[values[l]] != bestScore && values[l] != contamination) {
					thisCov = 1.0 * Score[values[l]] / template_lengths[values[l]];
					if(thisCov > bestCov) {
						bestCov = thisCov;
						bestHitsCov = regionTemplates[0] + 1;
						regionTemplates[bestHitsCov] = values[l];
					} else if(thisCov == bestCov) {
						bestHitsCov++;
						regionTemplates[bestHitsCov] = values[l];
					}
				}
				Score[values[l]]--;
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k] != 0) {
			values = VR_scores[k];
			for(l = 1; l <= values[0]; l++) {
				if(Score_r[values[l]] != bestScore && values[l] != contamination) {
					thisCov = 1.0 * Score_r[values[l]] / template_lengths[values[l]];
					if(thisCov > bestCov) {
						HIT = -1;
						bestCov = thisCov;
						bestHitsCov = regionTemplates[0] + 1;
						regionTemplates[bestHitsCov] = (-1) * values[l];
					} else if(thisCov == bestCov) {
						bestHitsCov++;
						regionTemplates[bestHitsCov] = (-1) * values[l];
					}
				}
				Score_r[values[l]]--;
			}
			VR_scores[k] = 0;
		}
	}
	regionTemplates[0] = bestHitsCov;
	print_ankers(regionTemplates, qseq_int + start, stop - start, HIT * bestScore, header, header_len);
	
	return HIT;
}

int ankerAndClean_MEM(int *regionTemplates, int *qseq_int, int start, int stop, int HIT, int bestScore, int start_cut, int end_cut, char *header, int header_len) {
	
	int k, l;
	unsigned *values;
	
	print_ankers(regionTemplates, qseq_int + start, stop - start, HIT * bestScore, header, header_len);
	
	/* clean up scores */
	for(k = start_cut; k <= end_cut; k++) {
		if(VF_scores[k] != 0) {
			values = VF_scores[k];
			for(l = 1; l <= values[0]; l++) {
				Score[values[l]]--;
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k] != 0) {
			values = VR_scores[k];
			for(l = 1; l <= values[0]; l++) {
				Score_r[values[l]]--;
			}
			VR_scores[k] = 0;
		}
	}
	
	return HIT;
}

int ankerAndClean_1t1(int *regionTemplates, int *qseq_int, int start, int stop, int HIT, int bestScore, int start_cut, int end_cut, char *header, int header_len) {
	print_ankers(regionTemplates, qseq_int, stop, HIT * bestScore, header, header_len);
	return 1;
}

int deConAnkers(int *regionTemplates, int *qseq_int, int start, int stop, int HIT, int bestScore, int start_cut, int end_cut, char *header, int header_len) {
	
	int contPos;
	
	if((contPos = find_contamination(regionTemplates)) != -1) {
		regionTemplates[contPos] = regionTemplates[regionTemplates[0]];
		regionTemplates[0]--;
	}
	if((contPos = find_contamination2(regionTemplates, (-1) * contamination)) != -1) {
		regionTemplates[contPos] = regionTemplates[regionTemplates[0]];
		regionTemplates[0]--;
	}
	if(regionTemplates[0] > 0) {
		return (*ankerPtr)(regionTemplates, qseq_int, start, stop, HIT, bestScore, start_cut, end_cut, header, header_len);
	} else {
		return 0;
	}
}

FILE * printFrags(char *filename, struct frag **alignFrags) {
	
	int i, seqlen;
	FILE *OUT;
	struct frag *alignFrag, *next;
	
	OUT = fopen(filename, "wb");
	if(!OUT) {
		fprintf(stderr, "File coruption: %s\n", filename);
		exit(1);
	}
	for(i = 0; i < DB_size; i++) {
		if(alignFrags[i] != 0) {
			for(alignFrag = alignFrags[i]; alignFrag != 0; alignFrag = next) {
				next = alignFrag->next;
				
				seqlen = strlen(alignFrag->qseq);
				fwrite(&(int){i}, sizeof(int), 1, OUT);
				fwrite(&(int){seqlen}, sizeof(int), 1, OUT);
				fwrite(alignFrag->qseq, seqlen * sizeof(char), 1, OUT);
				fwrite(&(int){alignFrag->score}, sizeof(int), 1, OUT);
				fwrite(&(int){alignFrag->start}, sizeof(int), 1, OUT);
				fwrite(&(int){alignFrag->end}, sizeof(int), 1, OUT);
				seqlen = strlen(alignFrag->header)+1;
				fwrite(&(int){seqlen}, sizeof(int), 1, OUT);
				fwrite(alignFrag->header, seqlen * sizeof(char), 1, OUT);
				
				free(alignFrag->qseq);
				free(alignFrag->header);
				free(alignFrag);
			}
			alignFrags[i] = 0;
		}
	}
	fwrite(&(int){-1}, sizeof(int), 1, OUT);
	fclose(OUT);
	
	return fopen(filename, "rb");
}

void save_kmers(int *qseq_int, int *qseq_int_r, int seqlen, char *header, int header_len) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, k, l, end, HIT, bestHits, hitCounter, bestScore, bestScore_r, contPos, maxScore;
	long unsigned key;
	unsigned *values;
	
	if(seqlen < kmersize) {
		return;
	}
	
	bestScore = 0;
	bestScore_r = 0;
	end = seqlen - kmersize + 1;
	/* Search forward strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	for(k = 0; k < end && !HIT; k += kmersize) {
		if((*hashMap_get)(quinaryToDecimal(qseq_int, k, kmersize))) {
			HIT = 1;
		}
	}
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		bestTemplates[0] = 0;
		key = quinaryToDecimal(qseq_int, 0, kmersize);
		values = (*hashMap_get)(key);
		if(values) {
			for(l = 1; l <= values[0]; l++) {
				if(Score[values[l]] > 0) {
					Score[values[l]]++;
				} else {
					Score[values[l]] = 1;
					bestTemplates[0]++;
					bestTemplates[bestTemplates[0]] = values[l];
				}	
			}
			hitCounter++;
		}
		for(k = 1; k < end; k++) {
			key = (key - qseq_int[k - 1] * convertNum[0]) * 5 + qseq_int[k + kmersize - 1];
			values = (*hashMap_get)(key);
			/* update scores */
			if(values) {
				for(l = 1; l <= values[0]; l++) {
					if(Score[values[l]] > 0) {
						Score[values[l]]++;
					} else {
						Score[values[l]] = 1;
						bestTemplates[0]++;
						bestTemplates[bestTemplates[0]] = values[l];
					}	
				}
				hitCounter++;
			}
		}
		/* get best match(es) */
		if(hitCounter - (1.0 * end - hitCounter) / kmersize > 1) {
			bestHits = 0;
			for(l = 1; l <= bestTemplates[0]; l++) {
				if(Score[bestTemplates[l]] > bestScore) {
					bestScore = Score[bestTemplates[l]];
					bestHits = 1;
					bestTemplates[bestHits] = bestTemplates[l];
				} else if(Score[bestTemplates[l]] == bestScore) {
					bestHits++;
					bestTemplates[bestHits] = bestTemplates[l];
				}
				Score[bestTemplates[l]] = 0;
			}
			bestTemplates[0] = bestHits;
		} else {
			for(l = 1; l <= bestTemplates[0]; l++) {
				Score[bestTemplates[l]] = 0;
			}
			bestTemplates[0] = 0;
		}
	}
	
	/* search rc strand
	   reverse complement qseq */
	for(i = 0; i < seqlen; i++) {
		qseq_int_r[i] = com_bases[qseq_int[seqlen - 1 - i]];
	}
	/* Make quick check of the qseq */
	HIT = exhaustive;
	for(k = 0; k < end && !HIT; k += kmersize) {
		if((*hashMap_get)(quinaryToDecimal(qseq_int_r, k, kmersize))) {
			HIT = 1;
		}
	}
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		bestTemplates_r[0] = 0;
		key = quinaryToDecimal(qseq_int_r, 0, kmersize);
		values = (*hashMap_get)(key);
		if(values) {
			for(l = 1; l <= values[0]; l++) {
				if(Score[values[l]] > 0) {
					Score[values[l]]++;
				} else {
					Score[values[l]] = 1;
					bestTemplates_r[0]++;
					bestTemplates_r[bestTemplates_r[0]] = values[l];
				}	
			}
			hitCounter++;
		}
		for(k = 1; k < end; k++) {
			key = (key - qseq_int_r[k - 1] * convertNum[0]) * 5 + qseq_int_r[k + kmersize - 1];
			values = (*hashMap_get)(key);
			/* update scores */
			if(values) {
				for(l = 1; l <= values[0]; l++) {
					if(Score[values[l]] > 0) {
						Score[values[l]]++;
					} else {
						Score[values[l]] = 1;
						bestTemplates_r[0]++;
						bestTemplates_r[bestTemplates_r[0]] = values[l];
					}	
				}
				hitCounter++;
			}
		}
		
		/* get best match(es) */
		if(hitCounter - (1.0 * end - hitCounter) / kmersize > 1) {
			bestHits = 0;
			for(l = 1; l <= bestTemplates_r[0]; l++) {
				if(Score[bestTemplates_r[l]] > bestScore_r) {
					bestScore_r = Score[bestTemplates_r[l]];
					bestHits = 1;
					bestTemplates_r[bestHits] = bestTemplates_r[l];
				} else if(Score[bestTemplates_r[l]] == bestScore_r) {
					bestHits++;
					bestTemplates_r[bestHits] = bestTemplates_r[l];
				}
				Score[bestTemplates_r[l]] = 0;
			}
			bestTemplates_r[0] = bestHits;
		} else {
			for(l = 1; l <= bestTemplates_r[0]; l++) {
				Score[bestTemplates_r[l]] = 0;
			}
			bestTemplates_r[0] = 0;
		}
	}
	
	/* Validate best match */
	if((bestScore >= bestScore_r && bestScore - (1.0 * end - bestScore) / kmersize > 1) || (bestScore < bestScore_r && bestScore_r - (1.0 * end - bestScore_r) / kmersize > 1)) {
		/* merge results */
		if(bestScore == bestScore_r) {
			for(i = 1; i <= bestTemplates_r[0]; i++) {
				bestTemplates[0]++;
				bestTemplates[bestTemplates[0]] = (-1) * bestTemplates_r[i];
			}
			deConAnkers(bestTemplates, qseq_int, 0, seqlen, -1, bestScore, 0, 0, header, header_len);
			//print_ankers(bestTemplates, qseq_int, seqlen, (-1) * bestScore, header, header_len);
		} else if(bestScore_r > bestScore) {
			deConAnkers(bestTemplates_r, qseq_int_r, 0, seqlen, 1, bestScore, 0, 0, header, header_len);
			//print_ankers(bestTemplates_r, qseq_int_r, seqlen, bestScore, header, header_len);
		} else {
			deConAnkers(bestTemplates, qseq_int, 0, seqlen, 1, bestScore, 0, 0, header, header_len);
			//print_ankers(bestTemplates, qseq_int, seqlen, bestScore, header, header_len);
		}
	}
}

void save_kmers_HMM(int *qseq_int, int *qseq_int_r, int seqlen, char *header, int header_len) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, k, l, m, end, HIT, bestHits, bestHitsCov, hitCounter, bestScore, contPos, maxScore, start, stop, start_cut, end_cut, part_seqlen, template;
	long unsigned key, key_r;
	unsigned *values, *values_r;
	double Ms, Ns, Ms_prev, Ns_prev, bestCov, thisCov;
	
	if(seqlen < kmersize) {
		return;
	}
	end = seqlen - kmersize + 1;
	
	for(i = 0; i < seqlen; i++) {
		qseq_int_r[i] = com_bases[qseq_int[seqlen - 1 - i]];
	}
	
	i = 0;
	while(i < end) {
		HIT = 0;
		if(exhaustive) {
			key = quinaryToDecimal(qseq_int, i, kmersize);
			key_r = quinaryToDecimal(qseq_int_r, end - i, kmersize);
			if((*hashMap_get)(key) || (*hashMap_get)(key_r)) {
				HIT = 1;
			}
			i++;
			while(i < end && !HIT) {
				key = (key - qseq_int[i - 1] * convertNum[0]) * 5 + qseq_int[i + kmersize - 1];
				key_r = (key_r - qseq_int_r[end - i + kmersize]) / 5 + qseq_int_r[end - i] * convertNum[0];
				if((*hashMap_get)(key) || (*hashMap_get)(key_r)) {
					HIT = 1;
				} else {
					i++;
				}
			}
		} else {
			while(i < end && !HIT) {
				if((*hashMap_get)(quinaryToDecimal(qseq_int, i, kmersize)) || (*hashMap_get)(quinaryToDecimal(qseq_int_r, end - i, kmersize))) {
					HIT = 1;
				} else {
					i += kmersize;
				}
			}
		}
		
		if(HIT) {
			bestScore = 0;
			bestTemplates[0] = 0;
			hitCounter = 1;
			key = quinaryToDecimal(qseq_int, i, kmersize);
			VF_scores[i] = (*hashMap_get)(key);
			
			key_r = quinaryToDecimal(qseq_int_r, end - i, kmersize);
			VR_scores[i] = (*hashMap_get)(key_r);
			
			j = i - 1;
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			while(j >= 0) { // && Ms >= Ns) {
				key = (key - qseq_int[j + kmersize]) / 5 + qseq_int[j] * convertNum[0];
				VF_scores[j] = (*hashMap_get)(key);
				key_r = (key_r - qseq_int_r[end - j - 1] * convertNum[0]) * 5 + qseq_int_r[end - j + kmersize - 1];
				//key_r = quinaryToDecimal(qseq_int_r, end - j, kmersize);
				VR_scores[j] = (*hashMap_get)(key_r);
				/* HMM */
				if(VF_scores[j] || VR_scores[j]) {
					hitCounter++;
					
					/* HMM */
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
						j--;
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
						j--;
						break;
					}
				}
				j--;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			start = j + 1;
			
			key = quinaryToDecimal(qseq_int, i, kmersize);
			key_r = quinaryToDecimal(qseq_int_r, end - i, kmersize);
			j = i + 1;
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			while(j < end) { // && Ms >= Ns) {
				key = (key - qseq_int[j - 1] * convertNum[0]) * 5 + qseq_int[j + kmersize - 1];
				VF_scores[j] = (*hashMap_get)(key);
				key_r = (key_r - qseq_int_r[end - j + kmersize]) / 5 + qseq_int_r[end - j] * convertNum[0];
				//key_r = quinaryToDecimal(qseq_int_r, end - j, kmersize);
				VR_scores[j] = (*hashMap_get)(key_r);
				/* HMM */
				if(VF_scores[j] || VR_scores[j]) {
					/* update scores */
					hitCounter++;
					/* HMM */
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
						j++;
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
						j++;
						break;
					}
				}
				j++;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			stop = j + kmersize - 1;
			
			if((stop - start) > kmersize && (start == 0 || stop == seqlen || (stop - start) > minLen)) {
				if(deCon) {
					for(k = start; k < j; k++) {
						if(VF_scores[k]) {
							values = VF_scores[k];
							if(values[values[0]] == contamination) {
								hitCounter--;
							}
						}
						if(VR_scores[k]) {
							values = VR_scores[k];
							if(values[values[0]] == contamination) {
								hitCounter--;
							}
						}
					}
				}
				
				if(hitCounter - (1.0 * (stop - start) - hitCounter) / kmersize > 1) {
					/* gain total scores and mapping templates for this region */
					bestTemplates[0] = 0;
					bestTemplates_r[0] = 0;
					for(k = start; k < j; k++) {
						if(VF_scores[k]) {
							values = VF_scores[k];
							for(l = 1; l <= values[0]; l++) {
								if(Score[values[l]] > 0) {
									Score[values[l]]++;
								} else {
									Score[values[l]] = 1;
									bestTemplates[0]++;
									bestTemplates[bestTemplates[0]] = values[l];
								}
							}
						}
						if(VR_scores[k]) {
							values_r = VR_scores[k];
							for(l = 1; l <= values_r[0]; l++) {
								if(Score_r[values_r[l]] > 0) {
									Score_r[values_r[l]]++;
								} else {
									Score_r[values_r[l]] = 1;
									bestTemplates_r[0]++;
									bestTemplates_r[bestTemplates_r[0]] = values_r[l];
								}	
							}
						}
					}
					
					/* clear scores */
					HIT = 1;
					while(HIT != 0) {
						/* get best score */
						bestScore = 0;
						HIT = 1;
						bestHits = 0;
						for(k = 1; k <= bestTemplates[0]; k++) {
							template = bestTemplates[k];
							if(Score[template] > bestScore) {
								bestScore = Score[template];
								bestHits = 1;
								regionTemplates[bestHits] = template;
							} else if(Score[template] == bestScore) {
								if(Score[template] == 0) {
									bestTemplates[k] = bestTemplates[bestTemplates[0]];
									bestTemplates[0]--;
									k--;
								} else {
									bestHits++;
									regionTemplates[bestHits] = template;
								}
							}
						}
						for(k = 1; k <= bestTemplates_r[0]; k++) {
							template = bestTemplates_r[k];
							if(Score_r[template] > bestScore) {
								bestScore = Score_r[template];
								bestHits = 1;
								regionTemplates[bestHits] = (-1) * template;
								HIT = -1;
							} else if(Score_r[template] == bestScore) {
								if(Score_r[template] == 0) {
									bestTemplates_r[k] = bestTemplates_r[bestTemplates_r[0]];
									bestTemplates_r[0]--;
									k--;
								} else {
									bestHits++;
									regionTemplates[bestHits] = (-1) * template;
									HIT = -1;
								}
							}
						}
						regionTemplates[0] = bestHits;
						
						if(bestScore > 0) {
							start_cut = stop;
							for(k = 1; k <= regionTemplates[0]; k++) {
								template = (regionTemplates[k] < 0) ? (-1) * regionTemplates[k] : regionTemplates[k];
								for(l = start; l < start_cut; l++) {
									if(VR_scores[l] && uintpos_bin(VR_scores[l], template) != -1) {
										start_cut = l;
									}
									if(VF_scores[l] && uintpos_bin(VF_scores[l], template) != -1) {
										start_cut = l;
									}
								}
							}
							end_cut = start_cut;
							for(k = 1; k <= regionTemplates[0]; k++) {
								template = (regionTemplates[k] < 0) ? (-1) * regionTemplates[k] : regionTemplates[k];
								//for(l = end_cut; l < stop; l++) {
								for(l = stop - kmersize + 1; l > end_cut; l--) {
									if(VR_scores[l] && uintpos_bin(VR_scores[l], template) != -1) {
										end_cut = l;
									}
									if(VF_scores[l] && uintpos_bin(VF_scores[l], template) != -1) {
										end_cut = l;
									}
								}
							}
							if(bestScore - (1.0 * (end_cut - start_cut) - bestScore) / kmersize > 1) {
								HIT = (*deConPtr)(regionTemplates, qseq_int, start, stop, HIT, bestScore, start_cut, end_cut, header, header_len);
								
								start_cut += kmersize - 1;
								end_cut -= kmersize + 1;
								
								for(l = start_cut; l < end_cut; l++) {
									qseq_int[l] = 4;
								}
								
							} else {
								for(k = 1; k <= bestTemplates[0]; k++) {
									Score[bestTemplates[k]] = 0;
								}
								for(k = 1; k <= bestTemplates_r[0]; k++) {
									Score_r[bestTemplates_r[k]] = 0;
								}
								for(l = start; l < stop; l++) {
									VR_scores[l] = 0;
									VF_scores[l] = 0;
								}
								HIT = 0;
							}
						} else {
							for(l = start; l < stop; l++) {
								VR_scores[l] = 0;
								VF_scores[l] = 0;
							}
							HIT = 0;
						}
					}
					for(l = start; l < stop; l++) {
						VR_scores[l] = 0;
						VF_scores[l] = 0;
					}
					for(k = 1; k <= bestTemplates[0]; k++) {
						Score[bestTemplates[k]] = 0;
					}
					for(k = 1; k <= bestTemplates_r[0]; k++) {
						Score_r[bestTemplates_r[k]] = 0;
					}
				} else {
					/* clear scores */
					for(l = start; l < stop; l++) {
						VR_scores[l] = 0;
						VF_scores[l] = 0;
					}
					
				}
			} else if(hitCounter > 0) {
				/* clear scores */
				for(l = start; l < stop; l++) {
					VR_scores[l] = 0;
					VF_scores[l] = 0;
				}
			}
			i = stop + 1;
		}	
	}
}

unsigned save_kmers_sparse(int seqlen, int *prefix, int prefix_len, struct hashMap_index *foundKmers, int *qseq_int) {
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, k, end;
	long unsigned key;
	unsigned hits = 0;
	
	/* get rc */
	for(i = seqlen; i > 0; i--) {
		qseq_int[2 * seqlen - i] = com_bases[qseq_int[i - 1]];
	}
	end = seqlen - kmersize + 1 - prefix_len;
	for(i = 0; i < 2 * seqlen; i += seqlen) {
		for(j = i; j < i + end; j++) {
			/* check if prefix matches */
			if(int_eq((qseq_int + j), prefix, prefix_len)) {
				/* check if k-mer exists */
				key = quinaryToDecimal(qseq_int, j + prefix_len, kmersize);
				if((*hashMap_get)(key)) {
					hashMap_CountIndex(foundKmers, key);
				}
				hits++;
			}
		}
	}
	return hits;
}

struct hashTable * collect_Kmers(int *Scores, int *Scores_tot, struct hashMap_index *foundKmers, struct Hit *hits) {
	unsigned i, j, *value;
	struct hashTable_index *node, *node_next;
	struct hashTable *kmerNode, *kmerList;
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	for(i = 0; i < foundKmers->size; i++) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			value = (*hashMap_get)(node->key);
			if(value) {
				hits->n++;
				hits->tot += node->value;
				
				kmerNode = malloc(sizeof(struct hashTable));
				if(kmerNode == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				kmerNode->value = malloc((value[0] + 1) * sizeof(unsigned));
				if(kmerNode->value == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				kmerNode->value[0] = value[0];
				for(j = 1; j <= value[0]; j++) {
					kmerNode->value[j] = value[j];
					Scores[value[j]]++;
					Scores_tot[value[j]] += node->value;
				}
				kmerNode->key = node->value;
				
				kmerNode->next = kmerList;
				kmerList = kmerNode;
			}
			free(node);
		}
	}
	free(foundKmers->table);
	
	return kmerList;
}

struct hashTable ** collect_Kmers_deCon(int *Scores, int *Scores_tot, struct hashMap_index *foundKmers, struct Hit *hits) {
	unsigned i, j, *value;
	struct hashTable_index *node, *node_next;
	struct hashTable *kmerNode, *kmerList;
	struct hashTable *decon_node, *deconList;
	struct hashTable **Returner;
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	deconList = 0;
	decon_node = 0;
	
	Returner = malloc(2 * sizeof(struct hashTable *));
	if(!Returner) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	for(i = 0; i < foundKmers->size; i++) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			value = (*hashMap_get)(node->key);
			if(value) {
				/* check for contamination */
				hits->n++;
				hits->tot += node->value;
				
				if(value[value[0]] == contamination) {
					decon_node = malloc(sizeof(struct hashTable));
					if(decon_node == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					decon_node->value = malloc((value[0] + 1) * sizeof(unsigned));
					if(decon_node->value == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					decon_node->value[0] = value[0];
					for(j = 1; j <= value[0]; j++) {
						decon_node->value[j] = value[j];
					}
					decon_node->key = node->value;
					
					decon_node->next = deconList;
					deconList = decon_node;
				} else {
					kmerNode = malloc(sizeof(struct hashTable));
					if(kmerNode == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					kmerNode->value = malloc((value[0] + 1) * sizeof(unsigned));
					if(kmerNode->value == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					kmerNode->value[0] = value[0];
					for(j = 1; j <= value[0]; j++) {
						kmerNode->value[j] = value[j];
						Scores[value[j]]++;
						Scores_tot[value[j]] += node->value;
					}
					kmerNode->key = node->value;
					
					kmerNode->next = kmerList;
					kmerList = kmerNode;
				}
			}
			free(node);
		}
	}
	free(foundKmers->table);
	
	Returner[0] = kmerList;
	Returner[1] = deconList;
	
	return Returner;
}

struct hashTable * withDraw_Kmers(int *Scores, int *Scores_tot, struct hashTable *kmerList, int template, struct Hit *hits) {
	unsigned i;
	struct hashTable *node, *node_next, *prev;;
	prev = 0;
	
	if(kmerList == 0) {
		return 0;
	}
	
	node = kmerList;
	while(node != 0) {
		if(uintpos_bin(node->value, template) != -1) {
			hits->n--;
			hits->tot -= node->key;
			for(i = 1; i <= node->value[0]; i++) {
				Scores[node->value[i]]--;
				Scores_tot[node->value[i]] -= node->key;
			}
			if(prev == 0) {
				kmerList = node->next;
				free(node->value);
				free(node);
				node = kmerList;
			} else {
				prev->next = node->next;
				free(node->value);
				free(node);
				node = prev->next;
			}
			/* early stopping */
			if(Scores[template] == 0 && Scores_tot[template] == 0) {
				return kmerList;
			}
			
		} else {
			prev = node;
			node = node->next;
		}
	}
	
	if(prev != 0) {
		prev->next = 0;
	}
	
	return kmerList;
}

struct Hit withDraw_Contamination(int *Scores, int *Scores_tot, struct hashTable *kmerList, struct hashTable *deConTable, int template, struct Hit hits) {
	unsigned i, j, belong;
	struct hashTable *node, *prev;
	struct hashTable *cont_node, *cont_prev;
	if(kmerList != 0) {
		prev = 0;
		node = kmerList;
		cont_prev = 0;
		cont_node = deConTable;
		while(node != 0) {
			/* check if k-mer belongs to template */
			belong = 0;
			for(j = node->value[0]; j > 0 && !belong; j--) {
				if(node->value[j] == template) {
					belong = 1;
				}
			}
			if(belong) { //withdraw score
				hits.n--;
				hits.tot -= node->key;
				for(j = node->value[0]; j > 0; j--) {
					Scores[node->value[j]]--;
					Scores_tot[node->value[j]] -= node->key;
				}
				cont_node->value = node->value;
				cont_node->key = node->key;
				cont_node->next = malloc(sizeof(struct hashTable));
				if(cont_node->next == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				cont_prev = cont_node;
				cont_node = cont_node->next;
				
				/* delete node */
				if(prev == 0) {
					kmerList = node->next;
					free(node);
					node = kmerList;
				} else {
					prev->next = node->next;
					free(node);
					node = prev->next;
				}
			} else {
				prev = node;
				node = node->next;
			}
		}
		if(cont_prev != 0) {
			free(cont_prev->next);
			cont_prev->next = 0;
		} else {
			free(deConTable);
			deConTable = 0;
		}
	} else {
		free(deConTable);
		deConTable = 0;
	}
	return hits;
}

void run_input(char **inputfiles, int fileCount, int minPhred, int fiveClip) {
	
	freopen(NULL, "wb", stdout);
	
	int i, j, k, FASTA, FASTQ, fileCounter, phredCut, start, end, *qseq_int;
	char *filename, *cmd, *zipped;
	struct qseqs *header, *qseq, *qual;
	struct FileBuff *inputfile;
	
	qseq_int = malloc(delta * sizeof(int));
	cmd = malloc(1);
	zipped = strdup(".gz");
	if(!cmd || !zipped || !qseq_int) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	header = setQseqs(128);
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* set input conversion */
	input_convertion = malloc(128 * sizeof(int));
	if(!input_convertion) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	for(i = 0; i < 128; i++) {
		input_convertion[i] = 4;
	}
	input_convertion['A'] = 0;
	input_convertion['T'] = 1;
	input_convertion['C'] = 2;
	input_convertion['G'] = 3;
	input_convertion['a'] = 0;
	input_convertion['t'] = 1;
	input_convertion['c'] = 2;
	input_convertion['g'] = 3;
	
	for(fileCounter = 0; fileCounter < fileCount; fileCounter++) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			cmd = realloc(cmd, (strlen(filename) + strlen("gunzip -c ") + 1));
			if(!cmd) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			sprintf(cmd, "gunzip -c %s", filename);
			popenFileBuff(inputfile, cmd, "r");
		} else if(strncmp(filename, "--", 2) == 0) {
			inputfile->file = stdin;
		} else {
			openFileBuff(inputfile, filename, "r");
		}
		fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		
		/* Get first char and determine the format */
		buffFileBuff(inputfile);
		FASTQ = 0;
		FASTA = 0;
		if(inputfile->buffer[0] == '@') { //FASTQ
			FASTQ = 1;
		} else if(inputfile->buffer[0] == '>') { //FASTA
			FASTA = 1;
		} else {
			fprintf(stderr, "Cannot determine format of file:\t%s\n", filename);
		}
		
		/* parse the file */
		if(FASTQ) {
			/* get phred scale */
			phredCut = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			while(FileBuffgetFq(inputfile, header, qseq, qual)) {
				/* trim */
				start = fiveClip;
				end = qseq->len - 1;
				while(end >= 0 && qual->seq[end] < phredCut) {
					end--;
				}
				end++;
				while(start < end && qual->seq[start] < phredCut) {
					start++;
				}
				qseq->len = end - start;
				/* print */
				if(qseq->len > kmersize) {
					/* translate to 2bit */
					if(qseq->len >= delta) {
						free(qseq_int);
						delta = qseq->len * 2;
						qseq_int = malloc(delta * sizeof(int));
						if(!qseq_int) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
					}
					
					decode_DNA((qseq->seq + start), qseq_int, qseq->len);
					fwrite(&qseq->len, sizeof(int), 1, stdout);
					fwrite(qseq_int, qseq->len * sizeof(int), 1, stdout);
					fwrite(&header->len, sizeof(int), 1, stdout);
					fwrite((header->seq + 1), header->len, 1, stdout);
				}
			}
		} else if(FASTA) {
			while(FileBuffgetFsa(inputfile, header, qseq)) {
				if(qseq->len > kmersize) {
					/* translate to 2bit */
					if(qseq->len >= delta) {
						free(qseq_int);
						delta = qseq->len * 2;
						qseq_int = malloc(delta * sizeof(int));
						if(!qseq_int) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
					}
					decode_DNA(qseq->seq, qseq_int, qseq->len);
					fwrite(&qseq->len, sizeof(int), 1, stdout);
					fwrite(qseq_int, qseq->len * sizeof(int), 1, stdout);
					fwrite(&header->len, sizeof(int), 1, stdout);
					fwrite((header->seq + 1), header->len, 1, stdout);
				}
			}
		}
		
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			pcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	
	free(cmd);
	free(zipped);
	free(qseq_int);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
}

void run_input_sparse(char **inputfiles, int fileCount, int minPhred, int fiveClip) {
	
	freopen(NULL, "wb", stdout);
	
	int i, j, k, FASTA, FASTQ, fileCounter, phredCut, start, end, kmerCount, *qseq_int;
	char *filename, *cmd, *zipped;
	struct qseqs *qseq, *qual;
	struct FileBuff *inputfile;
	
	qseq_int = malloc(delta * sizeof(int));
	cmd = malloc(1);
	zipped = strdup(".gz");
	if(!cmd || !zipped || !qseq_int) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	qseq = setQseqs(1024);
	qual = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* set input conversion */
	input_convertion = malloc(128 * sizeof(int));
	if(!input_convertion) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	for(i = 0; i < 128; i++) {
		input_convertion[i] = 4;
	}
	input_convertion['A'] = 0;
	input_convertion['T'] = 1;
	input_convertion['C'] = 2;
	input_convertion['G'] = 3;
	input_convertion['a'] = 0;
	input_convertion['t'] = 1;
	input_convertion['c'] = 2;
	input_convertion['g'] = 3;
	
	for(fileCounter = 0; fileCounter < fileCount; fileCounter++) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			cmd = realloc(cmd, (strlen(filename) + strlen("gunzip -c ") + 1));
			if(!cmd) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			sprintf(cmd, "gunzip -c %s", filename);
			popenFileBuff(inputfile, cmd, "r");
		} else if(strncmp(filename, "--", 2) == 0) {
			inputfile->file = stdin;
		} else {
			openFileBuff(inputfile, filename, "r");
		}
		fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		
		/* Get first char and determine the format */
		buffFileBuff(inputfile);
		FASTQ = 0;
		FASTA = 0;
		if(inputfile->buffer[0] == '@') { //FASTQ
			FASTQ = 1;
		} else if(inputfile->buffer[0] == '>') { //FASTA
			FASTA = 1;
		} else {
			fprintf(stderr, "Cannot determine format of file:\t%s\n", filename);
		}
		
		/* parse the file */
		if(FASTQ) {
			/* get phred scale */
			phredCut = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			while(FileBuffgetFqSeq(inputfile, qseq, qual)) {
				/* trim */
				start = fiveClip;
				end = qseq->len - 1;
				while(end >= 0 && qual->seq[end] < phredCut) {
					end--;
				}
				end++;
				while(start < end && qual->seq[start] < phredCut) {
					start++;
				}
				qseq->len = end - start;
				/* print */
				if(qseq->len > kmersize) {
					if(qseq->len >= delta) {
						free(qseq_int);
						delta = qseq->len * 2;
						qseq_int = malloc(delta * sizeof(int));
						if(!qseq_int) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
					}
					
					/* print */
					decode_DNA(qseq->seq, qseq_int, qseq->len);
					fwrite(&qseq->len, sizeof(int), 1, stdout);
					fwrite(qseq_int, qseq->len * sizeof(int), 1, stdout);
					
				}
			}
		} else if(FASTA) {
			while(FileBuffgetFsaSeq(inputfile, qseq)) {
				if(qseq->len > kmersize) {
					if(qseq->len >= delta) {
						free(qseq_int);
						delta = qseq->len * 2;
						qseq_int = malloc(delta * sizeof(int));
						if(!qseq_int) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
					}
					
					/* print */
					decode_DNA(qseq->seq, qseq_int, qseq->len);
					fwrite(&qseq->len, sizeof(int), 1, stdout);
					fwrite(qseq_int, qseq->len * sizeof(int), 1, stdout);
				}
			}
		}
		
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			pcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	
	free(cmd);
	free(zipped);
	free(qseq_int);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
}

void save_kmers_batch(char *templatefilename, char *exePrev) {
	
	/* open pipe */
	FILE *inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "File coruption\n");
		exit(1);
	}
	freopen(NULL, "wb", stdout);
	FILE *templatefile;
	int i, file_len, *qseq_int, *qseq_int_r, l_len, header_size;
	char *header;
	time_t t0, t1;
	
	t0 = clock();
	
	/* initialize seqs */
	qseq_int = malloc(delta * sizeof(int));
	qseq_int_r = malloc(delta * sizeof(int));
	if(!qseq_int || !qseq_int_r) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	header_size = 0;
	header = NULL;
	
	/* load hash tree */
	file_len = strlen(templatefilename);
	templatefilename = realloc(templatefilename, (file_len + 17) * sizeof(char));
	if(!templatefilename) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	if(shm) {
		templates_shm = malloc(sizeof(struct hashMap_shm));
		if(!templates_shm) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		strcat(templatefilename, ".m");
		templatefilename[file_len + 2] = '\0';
		templatefile = fopen(templatefilename, "rb" );
		if(templatefile == NULL) {
			fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
			exit(-1);
		}
		hashMap_shm_load(templates_shm, templatefile, templatefilename);
		fclose(templatefile);
		hashMap_get = &hashMap_shm_getGlobal;
	} else {
		if(deCon) {
			strcat(templatefilename, ".decontaminated.b");
			templatefilename[file_len + strlen(".decontaminated.b")] = '\0';
			templatefile = fopen(templatefilename, "rb" );
		} else {
			strcat(templatefilename, ".b");
			templatefilename[file_len + strlen(".b")] = '\0';
			templatefile = fopen(templatefilename, "rb" );
		}
		if(templatefile == NULL) {
			fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
			exit(-1);
		}
		fread(&contamination, sizeof(int), 1, templatefile);
		templates = malloc(sizeof(struct hashMap));
		if(!templates) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		hashMap_load(templates, templatefile);
		fclose(templatefile);
		hashMap_get = &hashMap_getGlobal;
	}
	templatefilename[file_len] = '\0';
	/* initialize convertNum */
	convertNum = malloc(kmersize * sizeof(long unsigned));
	if(!convertNum) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	convertNum[kmersize - 1] = 1;
	for(i = kmersize - 2; i >= 0; i--) {
		convertNum[i] = convertNum[i+1] * 5;
	}
	/* calculate HMM parameters */
	HMM_param[0] = log(1 - pow(0.25, kmersize));
	HMM_param[1] = log(pow(0.25, kmersize));
	HMM_param[2] = log(1 - pow(0.25, kmersize - 1) * 0.75);
	HMM_param[3] = log(pow(0.25, kmersize - 1) * 0.75);
	HMM_param[4] = log(1 - 1.0 / kmersize * 0.75 * 0.25);
	HMM_param[5] = log(1.0 / kmersize * 0.75 * 0.25);
	HMM_param[6] = log(0.75);
	HMM_param[7] = log(0.25);
	
	
	DB_size = contamination + 1;
	bestTemplates = calloc(DB_size + 1, sizeof(int));
	bestTemplates_r = calloc(DB_size + 1, sizeof(int));
	regionTemplates = calloc(2 * DB_size + 1, sizeof(int));
	Score = calloc(DB_size, sizeof(int));
	Score_r = calloc(DB_size, sizeof(int));
	VF_scores = calloc(delta, sizeof(int*));
	VR_scores = calloc(delta, sizeof(int*));
	/* load lengths */
	template_lengths = malloc(DB_size * sizeof(int));
	if(!bestTemplates || !bestTemplates_r || !regionTemplates || !Score || !Score_r || !VF_scores || !VR_scores || !template_lengths) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	strcat(templatefilename, ".length.b");
	templatefilename[file_len + 9] = '\0';
	templatefile = fopen(templatefilename, "rb");
	if(!templatefile) {
		fprintf(stderr, "File coruption: %s\n", templatefilename);
		exit(1);
	}
	templatefilename[file_len] = '\0';
	fread(&DB_size, sizeof(int), 1, templatefile);
	fread(template_lengths, DB_size * sizeof(int), 1, templatefile);
	fclose(templatefile);
	minLen = template_lengths[0] / 2;
	for(i = 0; i < DB_size; i++) {
		if(minLen > (template_lengths[i] / 2)) {
			minLen = template_lengths[i] / 2;
		}
		/* convert to number of kmers */
		template_lengths[i] -= (kmersize - 1);
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mer ankers\n");
	
	while(fread(&l_len, sizeof(int), 1, inputfile)) {
		if(l_len >= delta) {
			free(qseq_int);
			free(qseq_int_r);
			free(VF_scores);
			free(VR_scores);
			delta = 2 * l_len;
			qseq_int = malloc(delta * sizeof(int));
			qseq_int_r = malloc(delta * sizeof(int));
			VF_scores = calloc(delta, sizeof(int*));
			VR_scores = calloc(delta, sizeof(int*));
			if(!qseq_int || !qseq_int_r || !VF_scores || !VR_scores) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		fread(qseq_int, l_len * sizeof(int), 1, inputfile);
		fread(&i, sizeof(int), 1, inputfile);
		if(i > header_size) {
			header_size = i;
			free(header);
			header = malloc(header_size);
			if(!header) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		fread(header, i, 1, inputfile);
		(*kmerScan)(qseq_int, qseq_int_r, l_len, header, i);
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used ankering query: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
	pclose(inputfile);
	
	/* detach DB */
	if(shm) {
		hashMap_shm_detach(templates_shm);
	}
}

void save_kmers_sparse_batch(char *templatefilename, char *outputfilename, char *exePrev, char ss) {
	
	/* open pipe */
	FILE *inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "File corruption.\n");
		exit(1);
	}
	FILE *templatefile;
	int i, file_len, prefix_len;
	int *prefix;
	time_t t0, t1;
	
	t0 = clock();
	/* load hash tree */
	file_len = strlen(templatefilename);
	templatefilename = realloc(templatefilename, (file_len + 29) * sizeof(char));
	if(!templatefilename) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	if(deCon) {
		strcat(templatefilename, ".sparse.decontaminated.b");
		templatefilename[file_len + strlen(".sparse.decontaminated.b")] = '\0';
		templatefile = fopen(templatefilename, "rb" );	
	} else {
		strcat(templatefilename, ".sparse.b");
		templatefilename[file_len + strlen(".sparse.b")] = '\0';
		templatefile = fopen(templatefilename, "rb" );
	}
	if(templatefile == NULL) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
		exit(-1);
	}
	fread(&contamination, sizeof(int), 1, templatefile);
	DB_size = contamination + 1;
	/* load prefix */
	fread(&prefix_len, sizeof(int), 1, templatefile);
	if(prefix_len != 0) {
		prefix = malloc(prefix_len * sizeof(int));
		if(!prefix) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		fread(prefix, prefix_len * sizeof(int), 1, templatefile);
	} else {
		prefix = malloc(sizeof(int));
		if(!prefix) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		fread(prefix, sizeof(int), 1, templatefile);
	}
	
	/* Load hashMap */
	if(shm) {
		fclose(templatefile);
		if(deCon) {
			templatefilename[file_len + strlen(".sparse.decontaminated.")] = 'm';
		} else {
			templatefilename[file_len + strlen(".sparse.")] = 'm';
		}
		templates_shm = malloc(sizeof(struct hashMap_shm));
		if(!templates_shm) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		templatefile = fopen(templatefilename, "rb" );
		if(!templatefile) {
			fprintf(stderr, "File coruption: %s\n", templatefilename);
			exit(1);
		}		
		hashMap_shm_load(templates_shm, templatefile, templatefilename);
		fclose(templatefile);
		hashMap_get = &hashMap_shm_getGlobal;
		
		templates = malloc(sizeof(struct hashMap));
		if(!templates) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		templates->n = templates_shm->n;
		templates->size = templates_shm->size;
	} else {
		templates = malloc(sizeof(struct hashMap));
		if(!templates) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		hashMap_load(templates, templatefile);
		fclose(templatefile);
		hashMap_get = &hashMap_getGlobal;
	}
	templatefilename[file_len] = '\0';
	/* initialize convertNum */
	convertNum = malloc(kmersize * sizeof(long unsigned));
	if(!convertNum) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	convertNum[kmersize - 1] = 1;
	for(i = kmersize - 2; i >= 0; i--) {
		convertNum[i] = convertNum[i+1] * 5;
	}
	/* load template attributes */
	load_DBs_Sparse(templatefilename);
	
	/* open output file */
	file_len = strlen(outputfilename);
	/*outputfilename = realloc(outputfilename, (file_len + 5) * sizeof(char));
	if(!outputfilename) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}*/
	strcat(outputfilename, ".spa");
	outputfilename[file_len + 4] = '\0';
	FILE *sparse_out = fopen(outputfilename, "w");
	if(!sparse_out) {
		fprintf(stderr, "File coruption: %s\n", outputfilename);
		exit(1);
	}
	
	struct hashMap_index *foundKmers;
	struct hashTable *kmerList;
	foundKmers = malloc(sizeof(struct hashMap_index));
	if(!foundKmers) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	foundKmers->size = 0;
	for(i = 0; i < DB_size; i++) {
		if(template_ulengths[i] > foundKmers->size) {
			foundKmers->size = template_ulengths[i];
		}
	}
	foundKmers->size *= 2;
	initialize_hashMap_index(foundKmers, foundKmers->size);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mer\n");
	unsigned Ntot;//, Nhits, w_Nhits, Nhits_tot, w_Nhits_tot;
	struct Hit Nhits, w_Nhits;
	int l_len, stop, *qseq_int;
	delta *= 2;
	qseq_int = malloc(delta * sizeof(int));
	if(!qseq_int) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	Ntot = 0;
	
	while(fread(&l_len, sizeof(int), 1, inputfile)) {
		if(2 * l_len >= delta) {
			free(qseq_int);
			delta = 4 * l_len;
			qseq_int = malloc(delta * sizeof(int));
			if(!qseq_int) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		fread(qseq_int, l_len * sizeof(int), 1, inputfile);
		Ntot += save_kmers_sparse(l_len, prefix, prefix_len, foundKmers, qseq_int);
	}
	pclose(inputfile);
	t1 = clock();
	fprintf(stderr, "#\n# Total time used to identify k-mers in query: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding best matches and output results.\n");
	int *Scores, *Scores_tot, *w_Scores, *w_Scores_tot, *SearchList;
	int lenTresh;
	Scores = calloc(DB_size, sizeof(int));
	Scores_tot = calloc(DB_size, sizeof(int));
	SearchList = malloc(DB_size * sizeof(int));
	if(!Scores || !Scores_tot || !SearchList) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	int template, score, tmp_score;
	double depth, cover, expected, q_value, p_value, tot_depth, tot_cover, query_cover, tot_query_cover;
	double tmp_depth, tmp_cover, tmp_expected, tmp_q, tmp_p;
	double etta = 1.0e-6;
	
	fprintf(sparse_out, "#Template\tNum\tScore\tExpected\tTemplate length\tquery_coverage\tCoverage\tDepth\ttot_query_coverage\ttot_coverage\ttot_depth\tq_value\tp_value\n");
	stop = 0;
	if(deCon) {
		/* start by removing contamination */
		unsigned element, *score_templates, belong;
		struct hashTable *deConTable;
		int score_add, score_tot_add;
		struct hashTable *node, *prev;
		
		/* collect scores */
		struct hashTable **Collecter = collect_Kmers_deCon(Scores, Scores_tot, foundKmers, &Nhits);
		kmerList = Collecter[0];
		deConTable = Collecter[1];
		free(Collecter);
		
		fprintf(stderr, "# total number of matches: %d of %d kmers\n", Nhits.tot, Ntot);
		/* copy scores */
		w_Scores = malloc(DB_size * sizeof(int));
		w_Scores_tot = malloc(DB_size * sizeof(int));
		if(!w_Scores || !w_Scores_tot) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		
		for(i = 0; i < DB_size; i++) {
			w_Scores[i] = Scores[i];
			w_Scores_tot[i] = Scores_tot[i];
			if(Scores[i] == 0) {
				SearchList[i] = 0;
			} else {
				SearchList[i] = 1;
			}
		}
		w_Nhits.n = Nhits.n;
		w_Nhits.tot = Nhits.tot;
		
		if(w_Scores[contamination] > 0 || w_Scores_tot[contamination] > 0) {
			fprintf(stderr, "# Failed at removing contamination\n");
			exit(-1);
		} else {
			SearchList[contamination] = 0;
		}
		
		/* get best matches */
		while(! stop) {
			/* get best match, depth first then coverage */
			depth = 0;
			cover = 0;
			score = 0;
			if(ss == 'q') {
				for(i = 0; i < DB_size; i++) {
					if(SearchList[i] && w_Scores_tot[i] >= score) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						if(tmp_cover >= ID_t) {
							tmp_score = w_Scores_tot[i];
							tmp_depth = 1.0 * w_Scores_tot[i] / template_lengths[i];
							if(tmp_score > score || (tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && template_ulengths[i] > template_ulengths[template]))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else if (ss == 'd') {
				for(i = 0; i < DB_size; i++) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						if(tmp_cover >= ID_t) {
							tmp_score = w_Scores_tot[i];
							tmp_depth = 1.0 * w_Scores_tot[i] / template_lengths[i];
							if(tmp_depth > depth || (tmp_depth == depth && (tmp_cover > cover || (tmp_cover == cover && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								// calculate p_value
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else {
				for(i = 0; i < DB_size; i++) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						if(tmp_cover >= ID_t) {
							tmp_score = w_Scores_tot[i];
							tmp_depth = 1.0 * w_Scores_tot[i] / template_lengths[i];
							if(tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								/* calculate p_value */
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			}
			
			/* validate best match */
			if(cover >= ID_t) {
				/* with draw contamination k-mers matching this template */
				score_add = 0;
				score_tot_add = 0;
				if(deConTable != 0) {
					prev = 0;
					node = deConTable;
					while(node != 0) {
						if(uintpos_bin(node->value, template) != -1) {
							score_add++;
							score_tot_add += node->key;
							if(prev == 0) {
								deConTable = node->next;
								free(node);
								free(node->value);
								node = deConTable;
							} else {
								prev->next = node->next;
								free(node);
								free(node->value);
								node = prev->next;
							}
						} else {
							prev = node;
							node = node->next;
						}
					}
				}
				
				/* Calculate new attributes */
				query_cover = 100.0 * (w_Scores_tot[template] + score_tot_add) / Ntot;
				cover = 100.0 * (w_Scores[template] + score_add) / template_ulengths[template];
				depth = 1.0 * (w_Scores_tot[template] + score_tot_add) / template_lengths[template];
				/* calc tot values */
				tot_cover = 100.0 * (Scores[template] + score_add) / template_ulengths[template];
				tot_depth = 1.0 * (Scores_tot[template] + score_tot_add) / template_lengths[template];
				tot_query_cover = 100.0 * (Scores_tot[template] + score_tot_add) / Ntot;
				
				/* output results */
				fprintf(sparse_out, "%s\t%d\t%d\t%d\t%d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%f\t%e\n", 
					template_names[template], template, score, (int) expected, template_lengths[template], query_cover, cover, depth, tot_query_cover, tot_cover, tot_depth, q_value, p_value);
				
				/* update scores */
				kmerList = withDraw_Kmers(w_Scores, w_Scores_tot, kmerList, template, &w_Nhits);
				
				if(w_Scores[template] != 0 || w_Scores_tot[template] != 0) {
					fprintf(stderr, "# Failed updating the scores\n");
					//exit(-1);
					//stop = 1;
					SearchList[template] = 0;
				} else {
					SearchList[template] = 0;
				}
				if(kmerList == 0) {
					stop = 1;
				}
			} else {
				stop = 1;
			}
		}
	} else {
		/* collect scores */
		kmerList = collect_Kmers(Scores, Scores_tot, foundKmers, &Nhits);
		
		fprintf(stderr, "# total number of matches: %d of %d kmers\n", Nhits.tot, Ntot);
		/* copy scores */
		w_Scores = malloc(DB_size * sizeof(int));
		w_Scores_tot = malloc(DB_size * sizeof(int));
		if(!w_Scores || !w_Scores_tot) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		
		for(i = 0; i < DB_size; i++) {
			w_Scores[i] = Scores[i];
			w_Scores_tot[i] = Scores_tot[i];
			if(Scores[i] == 0) {
				SearchList[i] = 0;
			} else {
				SearchList[i] = 1;
			}
		}
		w_Nhits.n = Nhits.n;
		w_Nhits.tot = Nhits.tot;
		
		if(kmerList == 0) {
			stop = 1;
		}
		
		while(! stop) {
			/* get best match, depth first then coverage */
			depth = 0;
			cover = 0;
			score = 0;
			if(ss == 'q') {
				for(i = 0; i < DB_size; i++) {
					if(SearchList[i] && w_Scores_tot[i] >= score) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						if(tmp_cover >= ID_t) {
							tmp_score = w_Scores_tot[i];
							tmp_depth = 1.0 * w_Scores_tot[i] / template_lengths[i];
							if(tmp_score > score || (tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && template_ulengths[i] > template_ulengths[template]))))) {
								// calculate p_value
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else if (ss == 'd') {
				for(i = 0; i < DB_size; i++) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						if(tmp_cover >= ID_t) {
							tmp_score = w_Scores_tot[i];
							tmp_depth = 1.0 * w_Scores_tot[i] / template_lengths[i];
							if(tmp_depth > depth || (tmp_depth == depth && (tmp_cover > cover || (tmp_cover == cover && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								// calculate p_value
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			} else {
				for(i = 0; i < DB_size; i++) {
					if(SearchList[i]) {
						tmp_cover = 100.0 * w_Scores[i] / template_ulengths[i];
						if(tmp_cover >= ID_t) {
							tmp_score = w_Scores_tot[i];
							tmp_depth = 1.0 * w_Scores_tot[i] / template_lengths[i];
							if(tmp_cover > cover || (tmp_cover == cover && (tmp_depth > depth || (tmp_depth == depth && (tmp_score > score || (tmp_score == score && template_ulengths[i] > template_ulengths[template])))))) {
								// calculate p_value
								tmp_expected = 1.0 * (Nhits.tot - w_Scores_tot[i]) * template_ulengths[i] / (templates->n - template_ulengths[i] + etta);
								tmp_q = (tmp_score - tmp_expected) * (tmp_score - tmp_expected) / (tmp_score + tmp_expected);
								tmp_p = p_chisqr(tmp_q);
								if(tmp_p <= evalue && tmp_score > tmp_expected) {
									score = tmp_score;
									cover = tmp_cover;
									depth = tmp_depth;
									template = i;
									expected = tmp_expected;
									p_value = tmp_p;
									q_value = tmp_q;
								} else {
									SearchList[i] = 0;
								}
							}
						} else {
							SearchList[i] = 0;
						}
					}
				}
			}
			
			/* validate best match */
			if(cover >= ID_t) {
				/* output results */
				query_cover = 100.0 * w_Scores_tot[template] / Ntot;
				/* calc tot values */
				tot_cover = 100.0 * Scores[template] / template_ulengths[template];
				tot_depth = 1.0 * Scores_tot[template] / template_lengths[template];
				tot_query_cover = 100.0 * Scores_tot[template] / Ntot;
				fprintf(sparse_out, "%s\t%d\t%d\t%d\t%d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n", 
					template_names[template], template, score, (int) expected, template_lengths[template], query_cover, cover, depth, tot_query_cover, tot_cover, tot_depth, q_value, p_value);
				
				/* update scores */
				kmerList = withDraw_Kmers(w_Scores, w_Scores_tot, kmerList, template, &w_Nhits);
				if(w_Scores[template] != 0 || w_Scores_tot[template] != 0) {
					fprintf(stderr, "# Failed updating the scores\n");
					SearchList[template] = 0;
				} else {
					SearchList[template] = 0;
				}
				if(kmerList == 0) {
					stop = 1;
				}
			} else {
				stop = 1;
			}
		}
	}
	fclose(sparse_out);
	t1 = clock();
	fprintf(stderr, "# Total for finding and outputting best matches: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
	/* detach DB */
	if(shm) {
		hashMap_shm_detach(templates_shm);
	}
}

void SW_global(char *template, char *query, int k, struct aln *aligned) {
	
	int m, n, thisScore, t_len, q_len, aln_len, pos[2];
	
	t_len = strlen(template);
	q_len = strlen(query);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			aligned->score = 0;
			aligned->t[0] = '\0';
			aligned->s[0] = '\0';
			aligned->q[0] = '\0';
		} else if(t_len == 0) {
			memset(aligned->s, '_', q_len);
			aligned->s[q_len] = '\0';
			for(n = 0; n < q_len; n++) {
				aligned->t[n] = '5';	
			}
			aligned->t[q_len] = '\0';
			strncpy(aligned->q, query, q_len);
			aligned->q[q_len] = '\0';
			aligned->score = W1 + (q_len - 1) * U;
		} else {
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = '\0';
			for(m = 0; m < t_len; m++) {
				aligned->q[m] = '5';
			}
			aligned->q[t_len] = '\0';
			strncpy(aligned->t, template, t_len);
			aligned->t[t_len] = '\0';
			aligned->score = W1 + (t_len - 1) * U;
		}
		return;
	}
	
	/* check matrix size */
	if(2 * t_len >= SW_s || 2 * q_len >= SW_s) {
		for(m = 0; m < SW_s; m++) {
			free(D[m]);
			free(Q[m]);
			free(P[m]);
			free(E[m]);
		}
		free(D);
		free(Q);
		free(P);
		free(E);
		if(t_len > q_len) {
			SW_s = 2 * t_len + 1;
		} else {
			SW_s = 2 * q_len + 1;
		}
		D = malloc(SW_s * sizeof(int*));
		Q = malloc(SW_s * sizeof(int*));
		P = malloc(SW_s * sizeof(int*));
		E = malloc(SW_s * sizeof(int*));
		if(!D || !Q || !P || !E) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(m = 0; m < SW_s; m++) {
			D[m] = malloc(SW_s * sizeof(int));
			Q[m] = malloc(SW_s * sizeof(int));
			P[m] = malloc(SW_s * sizeof(int));
			E[m] = malloc(SW_s * sizeof(int));
			if(!D[m] || !Q[m] || !P[m] || !E[m]) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
	}
	
	/* fill in start penalties */
	if(k == 1) { // end alignment
		if(t_len > q_len) {
			for(m = t_len - 1; m > q_len; m--) {
				D[m][q_len] = W1 + (t_len - 1 - m) * U;
				Q[m][q_len] = D[m][q_len];
				E[m][q_len] = 5;
			}
			E[t_len - 1][q_len] = 4;
			for(m = q_len; m >= 0; m--) {
				D[m][q_len] = 0;
				Q[m][q_len] = 0;
				E[m][q_len] = 0;
			}
			for(n = q_len - 1; n >= 0; n--) {
				D[t_len][n] = W1 + (q_len - 1 - n) * U;
				P[t_len][n] = D[t_len][n];
				E[t_len][n] = 3;
			}
			E[t_len][q_len - 1] = 2;
		} else {
			for(m = t_len - 1; m >= 0; m--) {
				D[m][q_len] = W1 + (t_len - 1 - m) * U;
				Q[m][q_len] = D[m][q_len];
				E[m][q_len] = 5;
			}
			E[t_len - 1][q_len] = 4;
			for(n = q_len - 1; n > t_len; n--) {
				D[t_len][n] = W1 + (q_len - 1 - n) * U;
				P[t_len][n] = D[t_len][n];
				E[t_len][n] = 3;
			}
			for(n = t_len; n >= 0; n--) {
				D[t_len][n] = 0;
				P[t_len][n] = 0;
				E[t_len][n] = 0;
			}
			E[t_len][q_len - 1] = 2;
		}
		D[t_len][q_len] = 0;
		E[t_len][q_len] = 0;
	} else {
		for(m = t_len - 1; m >= 0; m--) {
			D[m][q_len] = W1 + (t_len - 1 - m) * U;
			Q[m][q_len] = D[m][q_len];
			E[m][q_len] = 5;
		}
		E[t_len - 1][q_len] = 4;
		for(n = q_len - 1; n >= 0; n--) {
			D[t_len][n] = W1 + (q_len - 1 - n) * U;
			P[t_len][n] = D[t_len][n];
			E[t_len][n] = 3;
		}
		E[t_len][q_len - 1] = 2;
		D[t_len][q_len] = 0;
		E[t_len][q_len] = 0;
	}
	
	for(m = t_len - 1; m >= 0; m--) {
		for(n = q_len - 1; n >= 0; n--) {
			/* Update D */
			D[m][n] = D[m + 1][n + 1] + d[template[m] - '0'][query[n] - '0'];
			E[m][n] = 1;
			
			/* update Q gap */
			Q[m][n] = D[m][n + 1] + W1;
			thisScore = Q[m][n + 1] + U;
			if(thisScore > Q[m][n]) {
				Q[m][n] = thisScore;
			}
			
			/* update P gap */
			P[m][n] = D[m + 1][n] + W1;
			thisScore = P[m + 1][n] + U;
			if(thisScore > P[m][n]) {
				P[m][n] = thisScore;
			}
			
			/* Update E trace */
			if(Q[m][n] > D[m][n]) {
				D[m][n] = Q[m][n];
				E[m][n] = 3;
			} 
			if(P[m][n] > D[m][n]) {
				D[m][n] = P[m][n];
				E[m][n] = 5;
			}
		}
	}
	
	/* get start position of alignment */
	aligned->score = D[0][0];
	pos[0] = 0;
	pos[1] = 0;
	if(k == -1) { // front alignment
		if(t_len > q_len) {
			for(m = 1; m < t_len; m++) {
				if(D[m][0] > aligned->score) {
					aligned->score = D[m][0];
					pos[0] = m;
				}
			}
		} else {
			for(n = 1; n < q_len; n++) {
				if(D[0][n] > aligned->score) {
					aligned->score = D[0][n];
					pos[1] = n;
				}
			}
		}
	}
	
	/* make back tracing */
	m = pos[0];
	n = pos[1];
	aln_len = 0;
	while(E[m][n] != 0) {
		if(E[m][n] == 1) {
			aligned->t[aln_len] = template[m];
			aligned->q[aln_len] = query[n];
			if(template[m] == query[n]) {
				aligned->s[aln_len] = '|';
			} else {
				aligned->s[aln_len] = '_';
			}
			m++;
			n++;
		} else if(E[m][n] >= 4) {
			aligned->t[aln_len] = template[m];
			aligned->q[aln_len] = '5';
			aligned->s[aln_len] = '_';
			m++;
		} else {
			aligned->t[aln_len] = '5';
			aligned->q[aln_len] = query[n];
			aligned->s[aln_len] = '_';
			n++;
		}
		aln_len++;
	}
	aligned->t[aln_len] = '\0';
	aligned->s[aln_len] = '\0';
	aligned->q[aln_len] = '\0';
	
}

void KMA_SW(const int template_name, const char *queryseq, int *qseq_int, struct aln *aligned, struct aln *Frag_align) {
	int i, j, prev, prev_index, stop, t_len, q_len, value, aln_len, frag_len, query_stop, template_stop, end;
	char base;
	long unsigned key;
	
	/* Extract indexes and template sequence */
	struct hashMap_index *template_index;
	const char *template; // const to avoid copy and changing the DB
	template = template_seqs[template_name];
	template_index = templates_index[template_name];
	
	/* allocate space for alignment */
	t_len = template_lengths[template_name];
	q_len = strlen(queryseq);
	aligned->pos = t_len;
	
	char *frag_q, *frag_t;
	frag_q = calloc((q_len + 2), sizeof(char));
	frag_t = calloc(2 * (q_len + 1) + 1, sizeof(char));
	if(!frag_q || !frag_t) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	aligned->s[0] = '\0';
	aligned->t[0] = '\0';
	aligned->q[0] = '\0';
	aligned->score = 0;
	/* Align the query and the template */
	prev = 0;
	prev_index = 0;
	i = 0;
	aln_len = 0;
	end = (q_len - kmersize + 1);
	key = quinaryToDecimal(qseq_int, i, kmersize);
	value = hashMap_getIndex(template_index, key);
	while(i < end && prev_index < t_len) {
		if(value > -1 && value >= prev_index) { // match, start prolonging the alignment
			if(value < aligned->pos) { // Start of alignment
				if(aligned->pos != t_len) { //control consistent increasing path
					aligned->s[0] = '\0';
					aligned->t[0] = '\0';
					aligned->q[0] = '\0';
					aligned->score = 0;
					free(frag_q);
					frag_q = NULL;
					free(frag_t);
					frag_t = NULL;
					return;
				}
				
				if(value - prev_index > 2 * (i - prev)) { // big leading gap, cut down
					strncpy(frag_t, template + (value - 2 * (i - prev)), 2 * (i - prev));
					frag_t[2 * (i - prev)] = '\0';
					strncpy(frag_q, queryseq + prev, i - prev);
					frag_q[i - prev] = '\0';
				} else if (i - prev > 2 * (value - prev_index)) {
					strncpy(frag_t, template + prev_index, value - prev_index);
					frag_t[value - prev_index] = '\0';
					strncpy(frag_q, queryseq + (i - 2 * (value - prev_index)), 2 * (value - prev_index));
					frag_q[2 * (value - prev_index)] = '\0';
				} else { // small leading gap
					strncpy(frag_t, template + prev_index, value - prev_index);
					frag_t[value - prev_index] = '\0';
					strncpy(frag_q, queryseq + prev, i - prev);
					frag_q[i - prev] = '\0';
				}
				
				SW_global(frag_t, frag_q, -1, Frag_align);
				frag_len = strlen(Frag_align->s);
				
				/* get start pos, relative to the template */
				aligned->pos = value;
				for(j = 0; j < frag_len; j++) {
					if(Frag_align->t[j] != '5') {
						aligned->pos--;
					}
				}
				
				strncpy(aligned->s, Frag_align->s, frag_len);
				strncpy(aligned->t, Frag_align->t, frag_len);
				strncpy(aligned->q, Frag_align->q, frag_len);
				aln_len = frag_len;
				aligned->s[aln_len] = '\0';
				aligned->t[aln_len] = '\0';
				aligned->q[aln_len] = '\0';
				aligned->score = Frag_align->score;
			} else { // Extend alignment
				strncpy(frag_q, queryseq + prev, i - prev);
				frag_q[i - prev] = '\0';
				if(value - prev_index > 2 * q_len) { // Gap is too big to give a positive score
					aligned->s[0] = '\0';
					aligned->t[0] = '\0';
					aligned->q[0] = '\0';
					
					free(frag_q);
					frag_q = NULL;
					free(frag_t);
					frag_t = NULL;
					aligned->score = 0;
					return;
				}
				strncpy(frag_t, template + prev_index, value - prev_index);
				frag_t[value - prev_index] = '\0';
				
				SW_global(frag_t, frag_q, 0, Frag_align);
				frag_len = strlen(Frag_align->s);
				
				strncpy((aligned->t + aln_len), Frag_align->t, frag_len);
				strncpy((aligned->s + aln_len), Frag_align->s, frag_len);
				strncpy((aligned->q + aln_len), Frag_align->q, frag_len);
				aln_len += frag_len;
				aligned->score += Frag_align->score;
			}
			
			/* Expand from current matchin k-mer */
			for(j = aln_len; j < aln_len + kmersize; j++) {
				aligned->s[j] = '|';
				base = queryseq[i + j - aln_len];
				aligned->t[j] = base;
				aligned->q[j] = base;
				base -= '0';
				aligned->score += d[base][base];
			}
			aln_len += kmersize;
			aligned->s[aln_len] = '\0';
			aligned->t[aln_len] = '\0';
			aligned->q[aln_len] = '\0';
			
			i += kmersize;
			prev_index = value + kmersize;
			
			/* prolong the sequence as long as possible */
			stop = 0;
			while(i < q_len && prev_index < t_len && !stop) {
				if(template[prev_index] == queryseq[i]) {
					aligned->s[aln_len] = '|';
					aligned->t[aln_len] = template[prev_index];
					aligned->q[aln_len] = queryseq[i];
					aligned->score += d[queryseq[i] - '0'][template[prev_index] - '0'];
					aln_len++;
					prev_index++;
					i++;
				} else {
					stop = 1;
				}
			}
			
			/* Update i and previous index */
			prev = i;
			if(i < end) {
				key = quinaryToDecimal(qseq_int, i, kmersize);
				value = hashMap_getIndex(template_index, key);
			}
		} else {
			i++;
			if(i < end) {
				key = (key - qseq_int[i - 1] * convertNum[0]) * 5 + qseq_int[i + kmersize - 1];
				value = hashMap_getIndex(template_index, key);
			}
		}
	}
	if(prev_index == 0 || prev == 0) { // No valid ankers were found
		aligned->s[0] = '\0';
		aligned->t[0] = '\0';
		aligned->q[0] = '\0';
		aligned->score = 0;
	} else {
		
		/* align last / trailing gap */
		if(t_len - prev_index > 2 * (q_len - prev + 1)) {
			strncpy(frag_t, template + prev_index, 2 * (q_len - prev + 1));
			frag_t[2 * (q_len - prev + 1)] = '\0';
			strncpy(frag_q, queryseq + prev, q_len - prev + 1);
			frag_q[q_len - prev + 1] = '\0';
		} else if(q_len - prev > 2 * (t_len - prev_index + 1)) {
			strncpy(frag_t, template + prev_index, t_len - prev_index);
			frag_t[t_len - prev_index] = '\0';
			strncpy(frag_q, queryseq + prev, 2 * (t_len - prev_index + 1));
			frag_q[2 * (t_len - prev_index + 1)] = '\0';
		} else {
			strncpy(frag_t, template + prev_index, t_len - prev_index);
			frag_t[t_len - prev_index] = '\0';
			strncpy(frag_q, queryseq + prev, q_len - prev + 1);
			frag_q[q_len - prev + 1] = '\0';
		}
		
		if(frag_t[0] != '\0' && frag_q[0] != '\0') {
			SW_global(frag_t, frag_q, 1, Frag_align);
			frag_len = strlen(Frag_align->s);
			strncpy((aligned->t + aln_len), Frag_align->t, frag_len);
			strncpy((aligned->s + aln_len), Frag_align->s, frag_len);
			strncpy((aligned->q + aln_len), Frag_align->q, frag_len);
			aln_len += frag_len;
			aligned->score += Frag_align->score;
			
			/* Cut alignment */
			template_stop = 0;
			query_stop = 0;
			aln_len--;
			while(aln_len >= 0 && !(query_stop && template_stop)) {
				if(aligned->t[aln_len] < '5') {
					template_stop = 1;
				} else {
					aligned->score -= U;
				}
				if(aligned->q[aln_len] < '5') {
					query_stop = 1;
				} else {
					aligned->score -= U;
				}
				aln_len--;
				frag_len--;
			}
			frag_len++;
			aln_len += 2;
		}
		aligned->s[aln_len] = '\0';
		aligned->t[aln_len] = '\0';
		aligned->q[aln_len] = '\0';
	}
	
	free(frag_q);
	frag_q = NULL;
	free(frag_t);
	frag_t = NULL;
}

void AlignFrag(char *frag_template, char *frag_query, int k, union DNA_tree *template_index, struct aln *aligned) {
	/* Initialize values */
	int i, j, index, pieces, aln_len;
	struct aln gap;
	int t_len = strlen(frag_template);
	int q_len = strlen(frag_query);
	
	k--;
	if(t_len < k && t_len <= q_len) {
		k = t_len;
	} else if(q_len < k && q_len <= t_len) {
		k = q_len;
	}
	
	/* Alignment done, exception */
	if(k <= 0) {
		if(q_len < t_len) {
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = '\0';
			for(i = 0; i < t_len - q_len; i++)
				aligned->q[i] = '5';
			for(i = 1; i <= q_len; i++)
				aligned->q[t_len-i] = frag_query[q_len-i];
			aligned->q[t_len] = '\0';
			strncpy(aligned->t, frag_template, t_len);
			aligned->t[t_len] = '\0';
		} else if(t_len < q_len) {
			for(i = 0; i < q_len; i++)
				aligned->s[i] = '_';
			aligned->s[q_len] = '\0';
			for(i = 0; i < q_len - t_len; i++)
				aligned->t[i] = '5';
			for(i = 1; i <= t_len; i++)
				aligned->t[q_len-i] = frag_template[t_len-i];
			aligned->t[q_len] = '\0';
			strncpy(aligned->q, frag_query, q_len);
			aligned->q[q_len] = '\0';
		} else {
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = '\0';
			strncpy(aligned->q, frag_query, q_len);
			aligned->q[q_len] = '\0';
			strncpy(aligned->t, frag_template, t_len);
			aligned->t[t_len] = '\0';
		}
		return;
	} else { //Allocate memory for alignment
		if(t_len < q_len) {
			pieces = q_len;
		} else {
			pieces = t_len;
		}	
	}
	
	char *submer;
	submer = malloc((k + 1) * sizeof(char));
	if(!submer) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	submer[k] = '\0';
	
	char **aln_pieces_t, **aln_pieces_s, **aln_pieces_q;
	aln_pieces_t = malloc((t_len + q_len + 1) * sizeof(char*));
	aln_pieces_s = malloc((t_len + q_len + 1) * sizeof(char*));
	aln_pieces_q = malloc((t_len + q_len + 1) * sizeof(char*));
	gap.t = calloc((t_len + q_len + 1), sizeof(char));
	gap.s = calloc((t_len + q_len + 1), sizeof(char));
	gap.q = calloc((t_len + q_len + 1), sizeof(char));
	if(!aln_pieces_t || !aln_pieces_s || !aln_pieces_q || !gap.t || !gap.s || !gap.q) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* Align using decreasing k-mers */
	pieces = 0;
	i = q_len;
	while(i >= k) {
		index = -1;
		if(tree_hasValue(template_index, frag_query, i - k, k)) {
			strncpy(submer, frag_query+(i-k), k);
			index = strpos_last(frag_template, submer);
		}
		/* Match found */
		if(index != -1) {
			/* Align previous gap */
			AlignFrag((frag_template + (index + k)), (frag_query + i), k, template_index, &gap);
			aln_len = strlen(gap.s);
			aln_pieces_t[pieces] = malloc((aln_len + 1) * sizeof(char));
			aln_pieces_s[pieces] = malloc((aln_len + 1) * sizeof(char));
			aln_pieces_q[pieces] = malloc((aln_len + 1) * sizeof(char));
			if(!aln_pieces_t[pieces] || !aln_pieces_s[pieces] || !aln_pieces_q[pieces]) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			strncpy(aln_pieces_t[pieces], gap.t, aln_len);
			strncpy(aln_pieces_s[pieces], gap.s, aln_len);
			strncpy(aln_pieces_q[pieces], gap.q, aln_len);
			aln_pieces_t[pieces][aln_len] = '\0';
			aln_pieces_s[pieces][aln_len] = '\0';
			aln_pieces_q[pieces][aln_len] = '\0';
			pieces++;
			
			/* add the match */
			aln_len = k;
			aln_pieces_t[pieces] = malloc((aln_len + 1) * sizeof(char));
			aln_pieces_s[pieces] = malloc((aln_len + 1) * sizeof(char));
			aln_pieces_q[pieces] = malloc((aln_len + 1) * sizeof(char));
			if(!aln_pieces_t[pieces] || !aln_pieces_s[pieces] || !aln_pieces_q[pieces]) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			strncpy(aln_pieces_t[pieces], submer, k);
			memset(aln_pieces_s[pieces], '|', k);
			strncpy(aln_pieces_q[pieces], submer, k);
			aln_pieces_t[pieces][k] = '\0';
			aln_pieces_s[pieces][k] = '\0';
			aln_pieces_q[pieces][k] = '\0';
			pieces++;
			
			/* discard mapped sequences */
			frag_template[index] = '\0';
			frag_query[i-k] = '\0';
			i -= (k + 1);
			index--;
			
			while(i >= 0 && index >= 0) {
				if(frag_query[i] == frag_template[index]) {
					aln_pieces_t[pieces] = malloc(2 * sizeof(char));
					aln_pieces_s[pieces] = malloc(2 * sizeof(char));
					aln_pieces_q[pieces] = malloc(2 * sizeof(char));
					if(!aln_pieces_t[pieces] || !aln_pieces_s[pieces] || !aln_pieces_q[pieces]) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					aln_pieces_t[pieces][0] = frag_template[index];
					aln_pieces_s[pieces][0] = '|';
					aln_pieces_q[pieces][0] = frag_query[i];
					aln_pieces_t[pieces][1] = '\0';
					aln_pieces_s[pieces][1] = '\0';
					aln_pieces_q[pieces][1] = '\0';
					i--;
					index--;
				} else {
					index = -1;
				}
			}
		} else {
			i--;
		}
	}
	/* align final gap */
	AlignFrag(frag_template, frag_query, k, template_index, &gap);
	aln_len = strlen(gap.s);
	aln_pieces_t[pieces] = malloc((aln_len + 1) * sizeof(char));
	aln_pieces_s[pieces] = malloc((aln_len + 1) * sizeof(char));
	aln_pieces_q[pieces] = malloc((aln_len + 1) * sizeof(char));
	if(!aln_pieces_t[pieces] || !aln_pieces_s[pieces] || !aln_pieces_q[pieces]) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	strncpy(aln_pieces_t[pieces], gap.t, aln_len);
	strncpy(aln_pieces_s[pieces], gap.s, aln_len);
	strncpy(aln_pieces_q[pieces], gap.q, aln_len);
	aln_pieces_t[pieces][aln_len] = '\0';
	aln_pieces_s[pieces][aln_len] = '\0';
	aln_pieces_q[pieces][aln_len] = '\0';
	
	aln_len = 0;
	for(i = pieces; i >= 0; i--) {
		for(j = 0; aln_pieces_t[i][j] != '\0'; j++) {
			aligned->t[aln_len] = aln_pieces_t[i][j];
			aligned->s[aln_len] = aln_pieces_s[i][j];
			aligned->q[aln_len] = aln_pieces_q[i][j];
			aln_len++;
		}
		free(aln_pieces_t[i]);
		aln_pieces_t[i] = NULL;
		free(aln_pieces_s[i]);
		aln_pieces_s[i] = NULL;
		free(aln_pieces_q[i]);
		aln_pieces_q[i] = NULL;
	}
	aligned->t[aln_len] = '\0';
	aligned->s[aln_len] = '\0';
	aligned->q[aln_len] = '\0';
	
	
	free(submer);
	submer = NULL;
	free(gap.t);
	free(gap.s);
	free(gap.q);
	gap.t = NULL;
	gap.s = NULL;
	gap.q = NULL;
	free(aln_pieces_t);
	aln_pieces_t = 0;
	free(aln_pieces_s);
	aln_pieces_s = 0;
	free(aln_pieces_q);
	aln_pieces_q = 0;
}

void AlignFragLast(const char *frag_t, const char *frag_q, int k, union DNA_tree *template_index, struct aln *aligned) {
	/* Initialize values */
	int i, j, index, aln_len, gap_len;
	struct aln gap;
	int t_len = strlen(frag_t);
	int q_len = strlen(frag_q);
	
	k--;
	if(t_len < k && t_len <= q_len) {
		k = t_len;
	} else if(q_len < k && q_len <= t_len) {
		k = q_len;
	}
	
	/* Alignment done, exception */
	if(k <= 0) {
		if(q_len < t_len) {
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = '\0';
			strncpy(aligned->q, frag_q, q_len);
			for(i = q_len; i < t_len; i++) {
				aligned->q[i] = '5';
			}
			aligned->q[t_len] = '\0';
			strncpy(aligned->t, frag_t, t_len);
			aligned->t[t_len] = '\0';
		} else if(t_len < q_len) {
			memset(aligned->s, '_', q_len);
			aligned->s[q_len] = '\0';
			strncpy(aligned->t, frag_t, t_len);
			for(i = 0; i < q_len; i++) {
				aligned->t[i] = '5';	
			}
			aligned->t[q_len] = '\0';
			strncpy(aligned->q, frag_q, q_len);
			aligned->q[q_len] = '\0';
		} else {
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = '\0';
			strncpy(aligned->q, frag_q, q_len);
			aligned->q[q_len] = '\0';
			strncpy(aligned->t, frag_t, t_len);
			aligned->t[t_len] = '\0';
		}
		return;
	}
	
	char *submer, *tmp_q, *tmp_t, *frag_template, *frag_query;
	submer = malloc((k + 1) * sizeof(char));
	if(!submer) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	submer[k] = '\0';
	frag_template = (char*)(frag_t);
	frag_query = (char*)(frag_q);
	tmp_q = malloc((q_len + 1) * sizeof(char));
	tmp_t = malloc((t_len + 1) * sizeof(char));
	
	gap.t = malloc((t_len + q_len + 1) * sizeof(char));
	gap.s = malloc((t_len + q_len + 1) * sizeof(char));
	gap.q = malloc((t_len + q_len + 1) * sizeof(char));
	if(!gap.t || !gap.s || !gap.q || !tmp_q || !tmp_t) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	aligned->t[0] = '\0';
	aligned->s[0] = '\0';
	aligned->q[0] = '\0';
	
	/* Align using decreasing k-mers */
	i = 0;
	aln_len = 0;
	while(i < (q_len - k + 1)) {
		index = -1;
		if(tree_hasValue(template_index, frag_query, i, k)) {
			strncpy(submer, frag_query+i, k);
			index = strpos(frag_template, submer);
		}
		/* Match found */
		if(index != -1) {
			/* Align previous gap */
			strncpy(tmp_t, frag_template, index);
			tmp_t[index] = '\0';
			strncpy(tmp_q, frag_query, i);
			tmp_q[i] = '\0';
			AlignFragLast(tmp_t, tmp_q, k, template_index, &gap);
			gap_len = strlen(gap.s);
			strncpy((aligned->t + aln_len), gap.t, gap_len);
			strncpy((aligned->s + aln_len), gap.s, gap_len);
			strncpy((aligned->q + aln_len), gap.q, gap_len);
			aln_len += gap_len;
			
			/* add the match */
			for(j = aln_len; j < aln_len + k; j++) {
				aligned->s[j] = '|';
			}
			strncpy((aligned->t + aln_len), submer, k);
			strncpy((aligned->q + aln_len), submer, k);
			aln_len += k;
			
			/* discard mapped sequences */
			frag_template = (frag_template + index + k);
			frag_query = (frag_query + i + k);
			t_len -= (index + k);
			q_len -= (i + k);
			i = 0;
			index = 0;
			while(i < q_len && index < t_len) {
				if(frag_query[0] == frag_template[0]) {
					aligned->t[aln_len] = frag_template[0];
					aligned->s[aln_len] = '|';
					aligned->q[aln_len] = frag_query[0];
					frag_template = (frag_template + 1);
					frag_query = (frag_query + 1);
					q_len--;
					t_len--;
				} else {
					index = t_len;
				}
			}
			
		} else {
			i++;
		}
	}
	/* align final gap */
	AlignFragLast(frag_template, frag_query, k, template_index, &gap);
	gap_len = strlen(gap.s);
	strncpy((aligned->t + aln_len), gap.t, gap_len);
	strncpy((aligned->s + aln_len), gap.s, gap_len);
	strncpy((aligned->q + aln_len), gap.q, gap_len);
	aln_len += gap_len;
	
	aligned->t[aln_len] = '\0';
	aligned->s[aln_len] = '\0';
	aligned->q[aln_len] = '\0';
	
	
	free(tmp_q);
	tmp_q = NULL;
	free(tmp_t);
	tmp_t = NULL;
	free(submer);
	submer = NULL;
	free(gap.t);
	free(gap.s);
	free(gap.q);
	gap.t = NULL;
	gap.s = NULL;
	gap.q = NULL;
}

void KMA(const int template_name, const char *queryseq, int *queryseq_int, struct aln *aligned, struct aln *Frag_align) {
	int i, j, prev, prev_index, stop, t_len, q_len, value, aln_len, frag_len, query_stop, template_stop, nextK, end;
	struct aln Frag_align_tmp;
	
	/* Extract indexes and template sequence */
	union DNA_tree *template_index;
	const char *template; // const to avoid copy and changing the DB
	template = template_seqs[template_name];
	template_index = &templates_align[template_name];
	
	/* allocate space for alignment */
	t_len = template_lengths[template_name];
	q_len = strlen(queryseq);
	aligned->pos = t_len;
	
	char *frag_q, *frag_t;
	frag_q = calloc((q_len + 2), sizeof(char));
	frag_t = calloc(2 * (q_len + 1) + 1, sizeof(char));
	if(!frag_q || !frag_t) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	aligned->s[0] = '\0';
	aligned->t[0] = '\0';
	aligned->q[0] = '\0';
	aligned->score = 0;
	/* Align the query and the template */
	prev = 0;
	prev_index = 0;
	i = 0;
	aln_len = 0;
	nextK = -1;
	end = q_len - kmersize + 1;
	while(i < end && prev_index < t_len) {
		value = tree_getValue(template_index, queryseq_int, i, kmersize);
		if(value > -1 && value >= prev_index) { // match, start prolonging the alignment
			/* align the previous gap */
			if(i - prev < 0 || prev < 0) { // Deletion exception
				prev = i;
			}
			
			if(value < aligned->pos) { // Start of alignment
				if(aligned->pos != t_len) { //control consistent increasing path
					aligned->s[0] = '\0';
					aligned->t[0] = '\0';
					aligned->q[0] = '\0';
					aligned->score = 0;
					free(frag_q);
					frag_q = NULL;
					free(frag_t);
					frag_t = NULL;
					return;
				}
				
				if(value - prev_index > 2 * (i - prev)) { // big leading gap, cut down
					strncpy(frag_t, template + (value - 2 * (i - prev)), 2 * (i - prev));
					frag_t[2 * (i - prev)] = '\0';
					strncpy(frag_q, queryseq + prev, i - prev);
					frag_q[i - prev] = '\0';
				} else if (i - prev > 2 * (value - prev_index)) {
					strncpy(frag_t, template + prev_index, value - prev_index);
					frag_t[value - prev_index] = '\0';
					strncpy(frag_q, queryseq + (i - 2 * (value - prev_index)), 2 * (value - prev_index));
					frag_q[2 * (value - prev_index)] = '\0';
				} else { // small leading gap
					strncpy(frag_t, template + prev_index, value - prev_index);
					frag_t[value - prev_index] = '\0';
					strncpy(frag_q, queryseq + prev, i - prev);
					frag_q[i - prev] = '\0';
				}
				
				AlignFrag(frag_t, frag_q, (-1) * nextK, template_index, Frag_align);
				nextK = -1;
				frag_len = strlen(Frag_align->s);
				
				/* Cut alignment */
				template_stop = 0;
				query_stop = 0;
				j = 0;
				while(j < frag_len && !(query_stop && template_stop)) {
					if(Frag_align->t[j] != '5') {
						template_stop = 1;
					}
					if(Frag_align->q[j] != '5') {
						query_stop = 1;
					}
					j++;
				}
				if(j != 0) {
					j--;
				}
				Frag_align_tmp.t = (Frag_align->t + j);
				Frag_align_tmp.s = (Frag_align->s + j);
				Frag_align_tmp.q = (Frag_align->q + j);
				frag_len -= j;
				
				/* get start pos, relative to the template */
				aligned->pos = value;
				for(j = 0; j < frag_len; j++) {
					if(Frag_align_tmp.t[j] != '5') {
						aligned->pos--;
					}
					if(Frag_align_tmp.s[j] == '|') {
						aligned->score += M;
					} else {
						aligned->score += MM;
					}
				}
				strncpy(aligned->s, Frag_align_tmp.s, frag_len);
				strncpy(aligned->t, Frag_align_tmp.t, frag_len);
				strncpy(aligned->q, Frag_align_tmp.q, frag_len);
				aln_len = frag_len;
				aligned->s[aln_len] = '\0';
				aligned->t[aln_len] = '\0';
				aligned->q[aln_len] = '\0';
				
			} else { // Extend alignment
				strncpy(frag_q, queryseq + prev, i - prev);
				frag_q[i - prev] = '\0';
				if(value - prev_index > 2 * q_len) { // Gap is too big to give a positive score
					aligned->s[0] = '\0';
					aligned->t[0] = '\0';
					aligned->q[0] = '\0';
					
					free(frag_q);
					frag_q = NULL;
					free(frag_t);
					frag_t = NULL;
					return;
				}
				strncpy(frag_t, template + prev_index, value - prev_index);
				frag_t[value - prev_index] = '\0';
				AlignFrag(frag_t, frag_q, (-1) * nextK, template_index, Frag_align);
				nextK = -1;
				frag_len = strlen(Frag_align->s);
				
				strncpy((aligned->t + aln_len), Frag_align->t, frag_len);
				strncpy((aligned->s + aln_len), Frag_align->s, frag_len);
				strncpy((aligned->q + aln_len), Frag_align->q, frag_len);
				aln_len += frag_len;
				for(j = 0; j < frag_len; j++) {
					if(Frag_align_tmp.s[j] == '|') {
						aligned->score += M;
					} else {
						aligned->score += MM;
					}
				}
			}
			
			/* Expand from current matchin k-mer */
			for(j = aln_len; j < aln_len + kmersize; j++) {
				aligned->s[j] = '|';
			}
			strncpy((aligned->t + aln_len), queryseq + i, kmersize);
			strncpy((aligned->q + aln_len), queryseq + i, kmersize);
			aln_len += kmersize;
			aligned->score += kmersize * M;
			
			prev = i + kmersize;
			prev_index = value + kmersize;
			
			/* prolong the sequence as long as possible */
			stop = 0;
			while(prev < q_len && prev_index < t_len && !stop) {
				if(template[prev_index] == queryseq[prev]) {
					aligned->s[aln_len] = '|';
					aligned->t[aln_len] = template[prev_index];
					aligned->q[aln_len] = queryseq[prev];
					aligned->score += M;
					aln_len++;
					prev_index++;
					prev++;
					i++;
				} else {
					stop = 1;
				}
			}
			
			/* Update i and previous index */
			i += kmersize;
		} else if(value < nextK) {
			nextK = value;
			i++;
		} else {
			i++;
		}
	}
	
	if(prev_index == 0 || prev == 0) { // No valid ankers were found
		aligned->s[0] = '\0';
		aligned->t[0] = '\0';
		aligned->q[0] = '\0';
	} else {
		
		/* align last / trailing gap */
		if(t_len - prev_index > 2 * (q_len - prev + 1)) {
			strncpy(frag_t, template + prev_index, 2 * (q_len - prev + 1));
			frag_t[2 * (q_len - prev + 1)] = '\0';
			strncpy(frag_q, queryseq + prev, q_len - prev + 1);
			frag_q[q_len - prev + 1] = '\0';
		} else if(q_len - prev > 2 * (t_len - prev_index + 1)) {
			strncpy(frag_t, template + prev_index, t_len - prev_index);
			frag_t[t_len - prev_index] = '\0';
			strncpy(frag_q, queryseq + prev, 2 * (t_len - prev_index + 1));
			frag_q[2 * (t_len - prev_index + 1)] = '\0';
		} else {
			strncpy(frag_t, template + prev_index, t_len - prev_index);
			frag_t[t_len - prev_index] = '\0';
			strncpy(frag_q, queryseq + prev, q_len - prev + 1);
			frag_q[q_len - prev + 1] = '\0';
		}
		AlignFragLast(frag_t, frag_q, (-1) * nextK, template_index, Frag_align);
		frag_len = strlen(Frag_align->s);
		strncpy((aligned->t + aln_len), Frag_align->t, frag_len);
		strncpy((aligned->s + aln_len), Frag_align->s, frag_len);
		strncpy((aligned->q + aln_len), Frag_align->q, frag_len);
		aln_len += frag_len;
		for(j = 0; j < frag_len; j++) {
			if(Frag_align_tmp.s[j] == '|') {
				aligned->score += M;
			} else {
				aligned->score += MM;
			}
		}
		/* Cut alignment */
		template_stop = 0;
		query_stop = 0;
		aln_len--;
		while(aln_len >= 0 && !(query_stop && template_stop)) {
			if(aligned->t[aln_len] < '5') {
				template_stop = 1;
			} else {
				aligned->score -= MM;
			}
			if(aligned->q[aln_len] < '5') {
				query_stop = 1;
			} else {
				aligned->score -= MM;
			}
			aln_len--;
		}
		aln_len += 2;
		aligned->s[aln_len] = '\0';
		aligned->t[aln_len] = '\0';
		aligned->q[aln_len] = '\0';
	}
	
	free(frag_q);
	frag_q = NULL;
	free(frag_t);
	frag_t = NULL;
}

struct assem assemble_KMA(int template, FILE **files, int file_count, FILE *frag_out, FILE *matrix_out, char *outputfilename, int *qseq_int, struct aln aligned, struct aln gap_align) {
	
	int i, j, k, t_len, q_len, aln_len, asm_len, s_len, start, end, template_stop, query_stop, bestScore, depthUpdate, bestBaseScore, bias, myBias, pos, file_len, matches, gaps, read_score, nextTemplate, file_i;
	double E_depth, q_value, p_value, tmp_Edepth, score;
	unsigned depth, coverScore;
	char *qseq, bestNuc, *header, *backBone, *consensus;
	int stats[4], **assembly, max_asmLen, *assemBias;
	struct assem aligned_assem;
	FILE *file;
	qseq = malloc(delta);
	header = malloc(1024);
	file_len = strlen(outputfilename);
	if(!qseq || !header) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* Allocate assembly arrays */
	t_len = template_lengths[template];
	backBone = calloc(3 * t_len + 1, sizeof(char));
	if(!backBone) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	strncpy(backBone, template_seqs[template], t_len);
	assemBias = calloc(t_len + 1, sizeof(int));
	assembly = calloc(3 * t_len + 1, sizeof(int*));
	if(!assembly || !assemBias) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	for(i = 0; i < t_len; i++) {
		assembly[i] = calloc(6, sizeof(int));
		if(!assembly[i]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
	}
	
	asm_len = t_len;
	max_asmLen = 3 * t_len;
	
	/* load reads of this template */
	file_i = 0;
	while(file_i < file_count) {
		file = files[file_i];
		if(file != 0) {
			fread(&nextTemplate, sizeof(int), 1, file);
			if(nextTemplate == template) {
				/* load frag */
				fread(&q_len, sizeof(int), 1, file);
				fread(qseq, q_len * sizeof(char), 1, file);
				qseq[q_len] = '\0';
				for(i = 0; i < q_len; i++) {
					qseq_int[i] = qseq[i] - '0';
				}
				stats[0] = 1;
				fread(&stats[1], sizeof(int), 1, file);
				fread(&stats[2], sizeof(int), 1, file);
				fread(&stats[3], sizeof(int), 1, file);
				
				fread(&i, sizeof(int), 1, file);
				fread(header, i, 1, file);
				
				/* Update assembly with read */
				/* Start with alignment */
				(*aligner)(template, qseq, qseq_int, &aligned, &gap_align);
				
				/* get read score */
				aln_len = strlen(aligned.s);
				start = aligned.pos;
				end = aligned.pos + aln_len;
				bias = 0;
				for(i = 0; i < aln_len; i++) {
					if(aligned.t[i] == '5') {
						bias++;
					}
				}
				
				/* Get score */
				end -= bias;
				read_score = aligned.score;
				if(aln_len > 0) {
					score = 1.0 * read_score / aln_len;
				} else {
					score = 0;
				}
				
				if(score > scoreT) {
					stats[0] = 1;
					stats[1] = read_score;
					stats[2] = start;
					stats[3] = end;
					
					/* Update backbone and counts */
					bias = assemBias[start];
					pos = start + bias;
					myBias = 0;
					for(i = 0; i < aln_len; i++) {
						if(aligned.t[i] == '5') {
							myBias++;
						}
						if(aligned.t[i] == backBone[pos]) {
							assembly[pos][aligned.q[i] - '0']++;
							pos++;
						} else if(backBone[pos] == '5') {
							assembly[pos][5]++;
							pos++;
							i--;
						} else if(aligned.t[i] == '5' && (start + i - myBias) >= 0) {
							/* realloc assembly, if needed */
							if(asm_len + 1 >= max_asmLen) {
								backBone = realloc(backBone, 2 * max_asmLen * sizeof(char));
								assembly = realloc(assembly, 2 * max_asmLen * sizeof(int*));
								if(!assembly || !backBone) {
									fprintf(stderr, "OOM\n");
									exit(1);
								}
								for(j = asm_len; j < 2 * max_asmLen; j++) {
									backBone[j] = '\0';
									assembly[j] = NULL;
								}
								max_asmLen *= 2;
							}
							/* insert into assembly */
							insert(backBone, '5', pos, asm_len);
							for(j = asm_len; j > pos; j--) {
								assembly[j] = assembly[j - 1];
							}
							assembly[pos] = calloc(6, sizeof(int));
							if(!assembly[pos]) {
								fprintf(stderr, "OOM\n");
								exit(1);
							}
							if((pos - 1) >= 0 && pos + 1 < asm_len) {
								for(j = 0; j < 6; j++) {
									assembly[pos][5] += assembly[pos+1][j];
									assembly[pos][5] += assembly[pos-1][j];
								}
								assembly[pos][5] /= 2;
							} else {
								for(j = 0; j < 6; j++) {
									assembly[pos][5] += assembly[pos+1][j];
								}
							}
							assembly[pos][aligned.q[i] - '0']++;
							asm_len++;
							
							for(j = start + i - myBias; j < t_len; j++) {
								assemBias[j]++;
							}
							pos++;
						}
					}
					
					/* Save fragment */
					for(i = 0; i < q_len; i++) {
						fprintf(frag_out, "%c", bases[qseq_int[i]]);
					}
					fprintf(frag_out, "\t%d\t%d\t%d\t%d\t%s\t%s\n", stats[0], stats[1], stats[2], stats[3], template_names[template], header);
				}
			} else if(nextTemplate == -1) {
				fclose(file);
				files[file_i] = 0;
				file_i++;
			} else if(nextTemplate < template) {
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len * sizeof(char) + 3 * sizeof(int), SEEK_CUR);
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len, SEEK_CUR);
			} else {
				/* Move pointer back */
				fseek(file, ((-1) * sizeof(int)), SEEK_CUR);
				file_i++;
			}
		} else {
			file_i++;
		}
	}
	free(header);
	free(assemBias);
	/* Make consensus assembly by majority voting */
	/* Get total abundance, for later use as expected depth */
	E_depth = 0.0;
	
	/* Convert backBone */
	for(i = 0; i < asm_len; i++) {
		backBone[i] = bases[backBone[i] - '0'];
	}
	
	/* Call nucleotides for the consensus */
	consensus = NULL;
	consensus = calloc(asm_len + 1, sizeof(char));
	if(!consensus) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	depth = 0;
	for(i = 0; i < asm_len; i++) {
		bestNuc = '-';
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; j++) {
			if(assembly[i][j] > bestScore) {
				bestScore = assembly[i][j];
				bestNuc = bases[j];
			}
			depthUpdate += assembly[i][j];
		}
		/* Check for minor base call */
		if(2 * bestScore < depthUpdate) {
			bestBaseScore = 0;
			bestNuc = 'n';
			for(j = 0; j < 5; j++) {
				if(assembly[i][j] > bestBaseScore) {
					bestBaseScore = assembly[i][j];
					bestNuc = tolower(bases[j]);
				}
			}
			bestScore = depthUpdate - assembly[i][5];
		}
		/* Determine base at current position */
		if(depthUpdate == 0) {
			consensus[i] = '-';
		} else {
			/* Use MC Neymars test to test significance of the base call */
			q_value = pow(bestScore - (depthUpdate - bestScore), 2) / depthUpdate;
			p_value = p_chisqr(q_value);
			if(p_value <= evalue && bestScore > (depthUpdate - bestScore)) {
				consensus[i] = bestNuc;
			} else if(bestNuc == '-') {
				consensus[i] = bestNuc;
			} else {
				/* Test significance of the base depth */
				/*tmp_Edepth = (E_depth - depthUpdate) / (t_len - 1);
				q_value = pow(depthUpdate - tmp_Edepth, 2) / (depthUpdate + tmp_Edepth);
				p_value = p_chisqr(q_value);
				if(p_value <= evalue && depthUpdate < tmp_Edepth) {
					consensus[i] = '-';
				} else {
					consensus[i] = tolower(bestNuc);
				}*/
				consensus[i] = tolower(bestNuc);
			}
		}
		if(bestNuc != '-')
			depth += depthUpdate;
	}
	
	/* Pepare and make alignment on consensus */
	aligned_assem.t = NULL;
	aligned_assem.s = NULL;
	aligned_assem.q = NULL;
	aligned_assem.s = calloc((asm_len + 1), sizeof(char));
	aligned_assem.t = calloc((asm_len + 1), sizeof(char));
	aligned_assem.q = calloc((asm_len + 1), sizeof(char));
	if(!aligned_assem.t || !aligned_assem.s || !aligned_assem.q) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	coverScore = 0;
	bias = 0;
	for(i = 0; i < asm_len; i++) {
		if(!(consensus[i] == '-' && backBone[i] == '-')) {
			aligned_assem.q[i - bias] = consensus[i];
			aligned_assem.t[i - bias] = backBone[i];
			if(tolower(consensus[i]) == tolower(backBone[i])) {
				aligned_assem.s[i - bias] = '|';
				coverScore++;
			} else {
				aligned_assem.s[i - bias] = '_';
			}
		} else {
			bias++;
		}
	}
	aligned_assem.cover = coverScore;
	aligned_assem.depth = depth;
	aligned_assem.t[asm_len - bias] = '\0';
	aligned_assem.s[asm_len - bias] = '\0';
	aligned_assem.q[asm_len - bias] = '\0';
	
	/* print matrix */
	if(print_matrix && coverScore > 0) {
		fprintf(matrix_out, "#%s\n", template_names[template]);
		for(i = 0; i < asm_len; i++) {
			fprintf(matrix_out, "%c", backBone[i]);
			for(j = 0; j < 6; j++) {
				fprintf(matrix_out, "\t%d", assembly[i][j]);
			}
			fprintf(matrix_out, "\n");
		}
		fprintf(matrix_out, "\n");
	}
	for(i = 0; i < asm_len; i++) {
		free(assembly[i]);
		assembly[i] = 0;
	}
	free(assembly);
	free(consensus);
	free(backBone);
	
	return aligned_assem;
}

struct assem assemble_KMA_dense(int template, FILE **files, int file_count, FILE *frag_out, FILE *matrix_out, char *outputfilename, int *qseq_int, struct aln aligned, struct aln gap_align) {
	
	int **assembly, i, j, k, t_len, q_len, aln_len, asm_len, s_len, start, end, template_stop, query_stop, bestScore, depthUpdate, bestBaseScore, bias, file_len, matches, gaps, read_score, nextTemplate, file_i;
	double E_depth, q_value, p_value, tmp_Edepth, score;
	unsigned depth, coverScore, pos;
	char *qseq, bestNuc, *backBone, *consensus, *header;
	int stats[4];
	struct assem aligned_assem;
	FILE *file;
	qseq = malloc(delta);
	header = malloc(1024);
	file_len = strlen(outputfilename);
	if(!qseq || !header) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* Allocate assembly arrays */
	t_len = template_lengths[template];
	asm_len = t_len;
	consensus = calloc(asm_len + 1, sizeof(char));
	backBone = calloc((t_len + 1), sizeof(char));
	if(!consensus || !backBone) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	strncpy(backBone, template_seqs[template], t_len);
	assembly = malloc((t_len + 1) * sizeof(int*));
	if(!assembly) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	for(i = 0; i < t_len + 1; i++) {
		assembly[i] = calloc(6, sizeof(int));
		if(!assembly[i]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
	}
	
	/* load reads of this template */
	file_i = 0;
	while(file_i < file_count) {
		file = files[file_i];
		if(file != 0) {
			fread(&nextTemplate, sizeof(int), 1, file);
			
			if(nextTemplate == template) {
				/* load frag */
				fread(&q_len, sizeof(int), 1, file);
				fread(qseq, q_len * sizeof(char), 1, file);
				qseq[q_len] = '\0';
				for(i = 0; i < q_len; i++) {
					qseq_int[i] = qseq[i] - '0';
				}
				stats[0] = 1;
				fread(&stats[1], sizeof(int), 1, file);
				fread(&stats[2], sizeof(int), 1, file);
				fread(&stats[3], sizeof(int), 1, file);
				
				fread(&i, sizeof(int), 1, file);
				fread(header, i, 1, file);
				
				/* Update assembly with read */
				/* Start with alignment */
				(*aligner)(template, qseq, qseq_int, &aligned, &gap_align);
				
				/* get read score */
				aln_len = strlen(aligned.s);
				start = aligned.pos;
				end = aligned.pos + aln_len;
				bias = 0;
				for(i = 0; i < aln_len; i++) {
					if(aligned.t[i] == '5') {
						bias++;
					}
				}
				
				/* Get score */
				end -= bias;
				read_score = aligned.score;
				if(aln_len > 0) {
					score = 1.0 * read_score / aln_len;
				} else {
					score = 0;
				}
				
				if(score > scoreT) {
					stats[0] = 1;
					stats[1] = read_score;
					stats[2] = start;
					stats[3] = end;
					
					/* Update backbone and counts */
					bias = 0;
					pos = start;
					for(i = 0; i < aln_len; i++) {
						if(aligned.t[i] == backBone[pos]) {
							assembly[pos][aligned.q[i] - '0']++;
							pos++;
						}
					}
					/* Save fragment */
					for(i = 0; i < q_len; i++) {
						fprintf(frag_out, "%c", bases[qseq_int[i]]);
					}
					fprintf(frag_out, "\t%d\t%d\t%d\t%d\t%s\t%s\n", stats[0], stats[1], stats[2], stats[3], template_names[template], header);
				}
			} else if (nextTemplate == -1) {
				fclose(file);
				files[file_i] = 0;
				file_i++;
			} else if(nextTemplate < template) {
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len * sizeof(char) + 3 * sizeof(int), SEEK_CUR);
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len, SEEK_CUR);
			} else {
				/* Move pointer back */
				fseek(file, ((-1) * sizeof(int)), SEEK_CUR);
				
				file_i++;
			}
		} else {
			file_i++;
		}
	}
	free(header);
	/* Make consensus assembly by majority voting */
	/* Get total abundance, for later use as expected delth */
	E_depth = 0.0;
	for(i = 0; i < asm_len; i++) {
		for(j = 0; j < 5; j++) {
			E_depth += assembly[i][j];
		}
	}
	
	/* Convert backBone */
	for(i = 0; i < asm_len; i++) {
		backBone[i] = bases[backBone[i] - '0'];
	}
	
	/* Call nucleotides for the consensus */
	depth = 0;
	for(i = 0; i < asm_len; i++) {
		bestNuc = '-';
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; j++) {
			if(assembly[i][j] > bestScore) {
				bestScore = assembly[i][j];
				bestNuc = bases[j];
			}
			depthUpdate += assembly[i][j];
		}
		/* Check for minor base call */
		if(2 * bestScore < depthUpdate) {
			bestBaseScore = 0;
			for(j = 0; j < 5; j++) {
				if(assembly[i][j] > bestBaseScore) {
					bestBaseScore = assembly[i][j];
					bestNuc = tolower(bases[j]);
				}
			}
			bestScore = depthUpdate - assembly[i][5];
		}
		/* Determine base at current position */
		if(depthUpdate == 0) {
			consensus[i] = '-';
		} else {
			/* Use MC Neymars test to test significance of the base call */
			q_value = pow(bestScore - (depthUpdate - bestScore), 2) / depthUpdate;
			p_value = p_chisqr(q_value);
			if(p_value <= evalue && bestScore > (depthUpdate - bestScore)) {
				consensus[i] = bestNuc;
			} else if(bestNuc == '-') {
				consensus[i] = bestNuc;
			} else {
				/* Test significance of the base depth */
				/*tmp_Edepth = (E_depth - depthUpdate) / (t_len - 1);
				q_value = pow(depthUpdate - tmp_Edepth, 2) / (depthUpdate + tmp_Edepth);
				p_value = p_chisqr(q_value);
				if(p_value <= evalue && depthUpdate < tmp_Edepth) {
					consensus[i] = '-';
				} else {
					consensus[i] = tolower(bestNuc);
				}*/
				consensus[i] = tolower(bestNuc);
			}
		}
		if(bestNuc != '-')
			depth += depthUpdate;
	}
	
	/* Pepare and make alignment on consensus */
	aligned_assem.s = calloc((asm_len + 1), sizeof(char));
	aligned_assem.t = calloc((asm_len + 1), sizeof(char));
	aligned_assem.q = calloc((asm_len + 1), sizeof(char));
	if(!aligned_assem.t || !aligned_assem.s || !aligned_assem.q) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	coverScore = 0;
	bias = 0;
	for(i = 0; i < asm_len; i++) {
		aligned_assem.q[i] = consensus[i];
		aligned_assem.t[i] = backBone[i];
		if(tolower(consensus[i]) == tolower(backBone[i])) {
			aligned_assem.s[i] = '|';
			coverScore++;
		} else {
			aligned_assem.s[i] = '_';
		}
	}
	aligned_assem.cover = coverScore;
	aligned_assem.depth = depth;
	
	/* print matrix */
	if(print_matrix && coverScore > 0) {
		fprintf(matrix_out, "#%s\n", template_names[template]);
		for(i = 0; i < asm_len; i++) {
			fprintf(matrix_out, "%c", backBone[i]);
			for(j = 0; j < 6; j++) {
				fprintf(matrix_out, "\t%d", assembly[i][j]);
			}
			fprintf(matrix_out, "\n");
		}
		fprintf(matrix_out, "\n");
	}
	
	
	free(consensus);
	consensus = NULL;
	free(backBone);
	backBone = NULL;
	for(i = 0; i < asm_len; i++) {
		free(assembly[i]);
		assembly[i] = 0;
	}
	free(assembly);
	assembly = NULL;
	
	return aligned_assem;
}

void update_Scores(char *qseq, int counter, int score, int *start, int *end, int *template, char *header, FILE *frag_out_raw) {
	
	int i;
	/* print read:  qseq, num, score, start, end, template */
	fprintf(frag_out_raw, "%s\t%d\t%d\t%d", qseq, counter, score, start[0]);
	for(i = 1; i < counter; i++) {
		fprintf(frag_out_raw, ",%d", start[i]);
	}
	fprintf(frag_out_raw, "\t%d", end[0]);
	for(i = 1; i < counter; i++) {
		fprintf(frag_out_raw, ",%d", end[i]);
	}
	fprintf(frag_out_raw, "\t%d", template[0]);
	for(i = 1; i < counter; i++) {
		fprintf(frag_out_raw, ",%d", template[i]);
	}
	fprintf(frag_out_raw, "\t%s\n", header);
	/* update scores */
	if(counter == 1) { //Only one best match
		if(template[0] < 0) {
			template[0] *= (-1);
		}
		alignment_scores[template[0]][0] += (1.0 * score) / (end[0] - start[0]);
		alignment_scores[template[0]][1]++;
	} else {
		for(i = 0; i < counter; i++) {
			if(template[i] < 0) {
				template[i] *= (-1);
			}
			alignment_scores[template[i]][0] += (1.0 * score) / (end[i] - start[i]);
		}
	}
}

void runKMA(char *templatefilename, char *outputfilename, char *exePrev) {
	
	/* open pipe */
	FILE *inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "File corruption.\n");
		exit(1);
	}
	int i, j, k, t_i, tmp_template, tmp_tmp_template, file_len, best_read_score, template, template_stop, query_stop, s_len, gaps, bestHits, read_counter, start, end, Nhits, pos, asm_len, q_len, frags, aln_len, depthUpdate, coverScore, start_pos, end_pos, bias, matches, read_score, bestTemplate, template_tot_ulen, t_len, rc_flag;
	int start_s, end_s, tmp_start, tmp_end, *matched_templates, *bestTemplates, *best_start_pos, *best_end_pos, *qseq_int, *qseq_int_r, stats[4], *w_scores, goOn, start_cut, end_cut, qbias, fragCount, fileCount, maxFrag;
	double etta, tmp_score, score, bestScore, bestBaseScore, E_depth, tmp_Edepth, depth, q_value, p_value, expected, cover, bestNum;
	char bestNuc, *outZipped, *qseq, *qseq_r, *tmp_string, *line, *line_sep, *template_s, *template_string, *template_acc, *header;
	struct aln aligned, gap_align;
	struct assem aligned_assem;
	struct frag **alignFrags, *alignFrag;
	FILE *align_in, *res_out, *frag_out, *alignment_out, *consensus_out, *frag_out_raw, *frag_in_raw, *fragmentIN, *matrix_out, **template_fragments;
	long unsigned *file_indexes;
	int line_size = 5000 * sizeof(char);
	time_t t0, t1;
	
	load_DBs_KMA(templatefilename);
	
	strcat(templatefilename, ".align.b");
	templatefilename[strlen(templatefilename) + strlen(".align.b")] = '\0';
	align_in = fopen(templatefilename, "rb");
	if(align_in == NULL) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
		exit(-1);
	}
	
	file_len = strlen(outputfilename);
	
	outputfilename[file_len] = '\0';
	outZipped = malloc((strlen("gunzip -c .frag_raw.gz") + file_len + 8) * sizeof(char));
	if(!outZipped) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	line = malloc(line_size);
	header = malloc(1024 * sizeof(char));
	if(!line || !header) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	
	
	
	/* load indexes */
	file_indexes = malloc(DB_size * sizeof(long unsigned));
	if(!file_indexes) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	file_indexes[0] = 0;
	for(i = 1; i < DB_size; i++) {
		file_indexes[i] = file_indexes[i - 1] + template_lengths[i - 1];
	}
	
	aligned.t = malloc(2*(delta + 1) * sizeof(char));
	aligned.s = malloc(2*(delta + 1) * sizeof(char));
	aligned.q = malloc(2*(delta + 1) * sizeof(char));
	gap_align.t = malloc(2*(delta + 1) * sizeof(char));
	gap_align.s = malloc(2*(delta + 1) * sizeof(char));
	gap_align.q = malloc(2*(delta + 1) * sizeof(char));
	qseq_int = calloc(delta + 1, sizeof(int));
	qseq_int_r = calloc(delta + 1, sizeof(int));
	qseq = malloc((delta + 1) * sizeof(char));
	qseq_r = malloc((delta + 1) * sizeof(char));
	if(!aligned.t || !aligned.s || !aligned.q || !gap_align.t || !gap_align.s || !gap_align.q || !qseq || !qseq_r || !qseq_int || !qseq_int_r) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	if(SW) {
		SW_s = delta;
		D = malloc(SW_s * sizeof(int*));
		Q = malloc(SW_s * sizeof(int*));
		P = malloc(SW_s * sizeof(int*));
		E = malloc(SW_s * sizeof(int*));
		if(!D || !Q || !P || !E) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(i = 0; i < SW_s; i++) {
			D[i] = malloc(SW_s * sizeof(int));
			Q[i] = malloc(SW_s * sizeof(int));
			P[i] = malloc(SW_s * sizeof(int));
			E[i] = malloc(SW_s * sizeof(int));
			if(!D[i] || !Q[i] || !P[i] || !E[i]) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		convertNum = malloc(kmersize * sizeof(long unsigned));
		if(!convertNum) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		convertNum[kmersize - 1] = 1;
		for(i = kmersize - 2; i >= 0; i--) {
			convertNum[i] = convertNum[i+1] * 5;
		}
	}
	/* etta = small value to avoid zero-divisio */
	etta = 1.0e-6;
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		outputfilename[file_len + 4] = '\0';
		res_out = fopen(outputfilename, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".frag.gz");
		outputfilename[file_len + 8] = '\0';
		sprintf(outZipped, "gzip -c > %s", outputfilename);
		frag_out = popen(outZipped, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".aln");
		outputfilename[file_len + 4] = '\0';
		alignment_out = fopen(outputfilename, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".fsa");
		outputfilename[file_len + 4] = '\0';
		consensus_out = fopen(outputfilename, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".frag_raw.gz");
		outputfilename[file_len + 12] = '\0';
		sprintf(outZipped, "gzip -c > %s", outputfilename);
		frag_out_raw = popen(outZipped, "w");
		outputfilename[file_len] = '\0';
		if(print_matrix) {
			strcat(outputfilename, ".mat.gz");
			outputfilename[file_len + 7] = '\0';
			sprintf(outZipped, "gzip -c > %s", outputfilename);
			matrix_out = popen(outZipped, "w");
			if(!matrix_out) {
				fprintf(stderr, "File coruption: %s\n", outputfilename);
				exit(1);
			}
			
			outputfilename[file_len] = '\0';
		} else {
			matrix_out = 0;
		}
	} else {
		fprintf(stderr, "# No output file specified!\n");
		exit( 1 );
	}
	if(!res_out || !frag_out || !consensus_out || !frag_out_raw || !alignment_out) {
		fprintf(stderr, "Could not create outputfiles\n");
		exit(1);
	}
	fprintf(stderr, "# Running KMA.\n");
	t0 = clock();
	
	/* Get alignments */
	read_counter = 0;
	matched_templates = malloc(2 * (DB_size + 1) * sizeof(int));
	best_start_pos = malloc(2 * DB_size * sizeof(int));
	best_end_pos = malloc(2 * DB_size * sizeof(int));
	bestTemplates = malloc(2 * DB_size * sizeof(int));
	if(!matched_templates || !best_start_pos || !best_end_pos || !bestTemplates) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	while(fread(&q_len, sizeof(int), 1, inputfile)) {
		if(q_len >= delta) {
			delta = 2 * q_len;
			free(qseq_int);
			free(qseq_int_r);
			qseq_int = malloc(delta * sizeof(int));
			qseq_int_r = malloc(delta * sizeof(int));
			free(aligned.t);
			free(aligned.s);
			free(aligned.q);
			aligned.t = malloc(2*(delta + 1) * sizeof(char));
			aligned.s = malloc(2*(delta + 1) * sizeof(char));
			aligned.q = malloc(2*(delta + 1) * sizeof(char));
			free(gap_align.t);
			free(gap_align.s);
			free(gap_align.q);
			gap_align.t = malloc(2*(delta + 1) * sizeof(char));
			gap_align.s = malloc(2*(delta + 1) * sizeof(char));
			gap_align.q = malloc(2*(delta + 1) * sizeof(char));
			free(qseq);
			free(qseq_r);
			qseq = malloc((delta + 1) * sizeof(char));
			qseq_r = malloc((delta + 1) * sizeof(char));
			if(!aligned.t || !aligned.s || !aligned.q || !gap_align.t || !gap_align.s || !gap_align.q || !qseq || !qseq_r || !qseq_int || !qseq_int_r) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		
		fread(qseq_int, q_len * sizeof(int), 1, inputfile);
		for(i = 0; i < q_len; i++) {
			qseq[i] = qseq_int[i] + '0';
		}
		qseq[q_len] = '\0';
		
		/* reverse complement seq */
		fread(&rc_flag, sizeof(int), 1, inputfile);
		if(rc_flag < 0) {
			for(i = 0; i < q_len; i++) {
				qseq_int_r[i] = com_bases[qseq_int[q_len - 1 - i]];
				qseq_r[i] = qseq_int_r[i] + '0';
			}
			qseq_r[q_len] = '\0';
		}
		
		/* Get number of matched templates */
		fread(&bestHits, sizeof(int), 1, inputfile);
		fread(matched_templates, (bestHits + 1) * sizeof(int), 1, inputfile);
		
		/* get header */
		fread(&t_i, sizeof(int), 1, inputfile);
		fread(header, t_i, 1, inputfile);
		
		bestScore = 0;
		best_read_score = 0;
		bestHits = 0;
		for(t_i = 1; t_i <= matched_templates[0]; t_i++) {
			template = matched_templates[t_i];
			/* check if index DB is loaded */
			if(template >= 0 && template_seqs[template] == 0) {
				(*alignLoadPtr)(template, align_in, file_indexes[template], SEEK_SET);
			} else if(template < 0 && template_seqs[(-1) * template] == 0) {
				(*alignLoadPtr)((-1) * template, align_in, file_indexes[(-1) * template], SEEK_SET);
			}
			
			/* align qseq */
			if(template < 0) {
				(*aligner)((-1) * template, qseq_r, qseq_int_r, &aligned, &gap_align);
			} else {
				(*aligner)(template, qseq, qseq_int, &aligned, &gap_align);
			}
			
			/* get read score */
			aln_len = strlen(aligned.s);
			start = aligned.pos;
			end = aligned.pos + aln_len;
			bias = 0;
			for(i = 0; i < aln_len; i++) {
				if(aligned.t[i] == '5') {
					bias++;
				}
			}
			
			/* Get score */
			end -= bias;
			read_score = aligned.score;
			if(aln_len > 0) {
				score = 1.0 * read_score / aln_len;
			} else {
				score = 0;
			}
			
			/* save best match(es) */
			if(score > scoreT) {
				if(score > bestScore) { // save as best match
					bestScore = score;
					best_read_score = read_score;
					bestTemplates[0] = template;
					best_start_pos[0] = start;
					best_end_pos[0] = end;
					bestHits = 1;
				} else if(score == bestScore && read_score > best_read_score) { // save as best match
					bestScore = score;
					best_read_score = read_score;
					bestTemplates[0] = template;
					best_start_pos[0] = start;
					best_end_pos[0] = end;
					bestHits = 1;
				} else if(score == bestScore && read_score == best_read_score) { // update best match
					bestTemplates[bestHits] = template;
					best_start_pos[bestHits] = start;
					best_end_pos[bestHits] = end;
					bestHits++;
				}
			}
		}
		
		if(best_read_score > kmersize) {
			update_Scores(qseq, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		}
	}
	pclose(inputfile);
	pclose(frag_out_raw);
	free(matched_templates);
	t1 = clock();
	fprintf(stderr, "#\n# KMA mapping time\t%d s.\n", (int) difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select KMA alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(int));
	if(!alignFrags || !w_scores) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	outputfilename[file_len] = '\0';
	sprintf(outZipped, "gunzip -c %s.frag_raw.gz", outputfilename);
	frag_in_raw = popen(outZipped, "r");
	if(!frag_in_raw) {
		fprintf(stderr, "File coruption: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[file_len] = '\0';
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	while( ! feof(frag_in_raw) && frag_in_raw) {
		qseq = fget_to(qseq, &delta, frag_in_raw, '\t');
		if(qseq[0] != '\0') {
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			bestHits = atoi(line);
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			score = atoi(line);
			for(i = 0; i < bestHits - 1; i++) {
				line = fget_to(line, &line_size, frag_in_raw, ',');
				best_start_pos[i] = atoi(line);
			}
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			best_start_pos[bestHits - 1] = atoi(line);
			for(i = 0; i < bestHits - 1; i++) {
				line = fget_to(line, &line_size, frag_in_raw, ',');
				best_end_pos[i] = atoi(line);
			}
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			best_end_pos[bestHits - 1] = atoi(line);
			for(i = 0; i < bestHits - 1; i++) {
				line = fget_to(line, &line_size, frag_in_raw, ',');
				bestTemplates[i] = atoi(line);
			}
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			bestTemplates[bestHits - 1] = atoi(line);
			
			fgets(header, 1024, frag_in_raw);
			header[strlen(header) - 1] = '\0';
			
			if(strlen(qseq) > kmersize) {
				/* Several mapped templates, choose best */
				if(bestHits > 1) {
					bestTemplate = 0;
					bestScore = 0;
					best_read_score = 0;
					bestNum = 0;
					start = 0;
					end = 0;
					/* iterate hits */
					for(i = 0; i < bestHits; i++) {
						tmp_tmp_template = bestTemplates[i];
						tmp_start = best_start_pos[i];
						tmp_end = best_end_pos[i];
						if(tmp_tmp_template < 0) {
							tmp_template = (-1) * tmp_tmp_template;
						} else {
							tmp_template = tmp_tmp_template;
						}
						tmp_score = alignment_scores[tmp_template][0] / template_lengths[tmp_template];
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template][0] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template][0];
							bestScore = tmp_score;
							bestNum = alignment_scores[tmp_template][1];
							start = tmp_start;
							end = tmp_end;
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template][0] > best_read_score) {
						} else if(alignment_scores[tmp_template][0] == best_read_score) {
							if(tmp_score > bestScore) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template][0];
								bestScore = tmp_score;
								bestNum = alignment_scores[tmp_template][1];
								start = tmp_start;
								end = tmp_end;
							} else if(tmp_score == bestScore && alignment_scores[tmp_template][1] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template][0];
								bestScore = tmp_score;
								bestNum = alignment_scores[tmp_template][1];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
					/* reverse complement seq */
					if(bestTemplate < 0) {
						bestTemplate *= (-1);
						strrc(qseq, strlen(qseq));
					}
					w_scores[bestTemplate] += score;
					
					/* load frag info */
					if(alignFrags[bestTemplate] == 0) {
						alignFrags[bestTemplate] = malloc(sizeof(struct frag));
						alignFrag = alignFrags[bestTemplate];
						if(!alignFrag) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
						alignFrag->qseq = strdup(qseq);
						alignFrag->header = strdup(header);
						alignFrag->score = (int) score;
						alignFrag->start = start;
						alignFrag->end = end;
						alignFrag->next = 0;
						if(!alignFrag->qseq || !alignFrag->header) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					} else {
						alignFrag = malloc(sizeof(struct frag));
						if(!alignFrag) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
						alignFrag->qseq = strdup(qseq);
						alignFrag->header = strdup(header);
						alignFrag->score = (int) score;
						alignFrag->start = start;
						alignFrag->end = end;
						alignFrag->next = alignFrags[bestTemplate];
						alignFrags[bestTemplate] = alignFrag;
						if(!alignFrag->qseq || !alignFrag->header) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					}
					
					fragCount++;
					if(fragCount >= maxFrag) {
						sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
						template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
						if(!template_fragments[fileCount]) {
							fprintf(stderr, "File coruption: %s\n", outputfilename);
							exit(1);
						}
						outputfilename[file_len] = '\0';
						fileCount++;
						fragCount = 0;
						/* control fileamount */
						if(fileCount >= DB_size) {
							template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
							if(!template_fragments) {
								fprintf(stderr, "OOM\n");
								exit(1);
							}
						}
					}
				} else {
					bestTemplate = *bestTemplates;
					/* reverse complement seq */
					if(bestTemplate < 0) {
						bestTemplate *= (-1);
						strrc(qseq, strlen(qseq));
					}
					w_scores[bestTemplate] += score;
					
					/* load frag info */
					if(alignFrags[bestTemplate] == 0) {
						alignFrags[bestTemplate] = malloc(sizeof(struct frag));
						alignFrag = alignFrags[bestTemplate];
						if(!alignFrag) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
						alignFrag->qseq = strdup(qseq);
						alignFrag->header = strdup(header);
						alignFrag->score = (int) score;
						alignFrag->start = start;
						alignFrag->end = end;
						alignFrag->next = 0;
						if(!alignFrag->qseq || !alignFrag->header) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					} else {
						alignFrag = malloc(sizeof(struct frag));
						if(!alignFrag) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
						alignFrag->qseq = strdup(qseq);
						alignFrag->header = strdup(header);
						alignFrag->score = (int) score;
						alignFrag->start = start;
						alignFrag->end = end;
						alignFrag->next = alignFrags[bestTemplate];
						alignFrags[bestTemplate] = alignFrag;
						if(!alignFrag->qseq || !alignFrag->header) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					}
					
					fragCount++;
					if(fragCount >= maxFrag) {
						sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
						template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
						if(!template_fragments[fileCount]) {
							fprintf(stderr, "File coruption: %s\n", outputfilename);
							exit(1);
						}
						outputfilename[file_len] = '\0';
						fileCount++;
						fragCount = 0;
						/* control fileamount */
						if(fileCount >= DB_size) {
							template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
							if(!template_fragments) {
								fprintf(stderr, "OOM\n");
								exit(1);
							}
							template_fragments[fileCount] = 0;
						}
					}
				}
			}
		}
	}
	/*
	free(best_start_pos);
	free(best_end_pos);
	free(bestTemplates);
	free(line);
	free(qseq);
	free(qseq_r);
	free(qseq_int_r);
	free(header);
	*/
	pclose(frag_in_raw);
	sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
	template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
	if(!template_fragments[fileCount]) {
		fprintf(stderr, "File coruption: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[file_len] = '\0';
	fileCount++;
	fragCount = 0;
	free(alignFrags);
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%d s.\n", (int) difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate length\tCoverage\tDepth\tq_value\tp_value\n");
	
	/* Get expected values */
	template_tot_ulen = 0;
	Nhits = 0;
	for(i = 0; i < DB_size; i++) {
		Nhits += w_scores[i];
		template_tot_ulen += template_lengths[i];
	}
	/* Do local assemblies of fragments mapping to the same template */
	aligned_assem.t = NULL;
	aligned_assem.s = NULL;
	aligned_assem.q = NULL;
	for(template = 0; template < DB_size; template++) {
		if(w_scores[template] > 0) {
			/* make p_value to see whether assembly is feasable */
			score = w_scores[template];
			t_len = template_lengths[template];
			expected = (Nhits - score) * (1.0 * t_len) / (template_tot_ulen - t_len + etta);
			q_value = pow(score - expected, 2) / (expected + score + etta);
			p_value  = p_chisqr(q_value);
			
			if(p_value <= evalue && (score > expected || evalue == 1.0)) {
				/* Do assembly */
				free(aligned_assem.t);
				free(aligned_assem.s);
				free(aligned_assem.q);
				aligned_assem.t = NULL;
				aligned_assem.s = NULL;
				aligned_assem.q = NULL;
				aligned_assem = assemblyPtr(template, template_fragments, fileCount, frag_out, matrix_out, outputfilename, qseq_int, aligned, gap_align);
				
				coverScore = aligned_assem.cover;
				depth = aligned_assem.depth;
				
				/* Create output statistics */
				score = w_scores[template];
				t_len = template_lengths[template];
				expected = (Nhits - score) * (1.0 * t_len) / (template_tot_ulen - t_len);
				q_value = pow(score - expected, 2) / (expected + score + etta);
				p_value  = p_chisqr(q_value);
				
				/* Depth and depth_corr */
				depth /= t_len;
				cover = 100.0 * coverScore / t_len;
				
				if(cover >= ID_t) {
					/* Seperate template and accession number */
					template_acc = (char*)(template_names[template]);
					
					/* Output result */
					fprintf(res_out, "%-12s\t%8d\t%8d\t%8d\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_acc, (int) score, (int) expected, t_len, cover, depth, q_value, p_value);//, templates_descriptions[template]);
					
					/* print alignment */
					aln_len = strlen(aligned_assem.s);
					fprintf(alignment_out, "# %s\n", template_names[template]);
					for(i = 0; i < aln_len; i += 60) {
						fprintf(alignment_out, "%-10s\t", "template:");
						for(j = i; j < i + 60 && aligned_assem.t[j] != '\0'; j++) {
							fprintf(alignment_out, "%c", aligned_assem.t[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "");
						for(j = i; j < i + 60 && aligned_assem.s[j] != '\0'; j++) {
							fprintf(alignment_out, "%c", aligned_assem.s[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "query:");
						for(j = i; j < i + 60 && aligned_assem.q[j] != '\0'; j++) {
							fprintf(alignment_out, "%c", aligned_assem.q[j]);
						}
						fprintf(alignment_out, "\n\n");
					}
					
					/* Print consensus */
					fprintf(consensus_out, ">%s\n", template_names[template]);
					if(ref_fsa) {
						for(i = 0; i < aln_len; i += 60) {
							for(j = i; j < i + 60 && aligned_assem.q[j]; j++) {
								if(aligned_assem.q[j]) {
									if(aligned_assem.q[j] != '-') {
										fprintf(consensus_out, "%c", aligned_assem.q[j]);
									} else {
										fprintf(consensus_out, "%c", 'n');
									}
								}
							}
							fprintf(consensus_out, "\n");
						}
					} else {
						bias = 0;
						for(i = 0; i < aln_len; i += 60 + bias) {
							bias = 0;
							for(j = i; j < i + 60 + bias && aligned_assem.q[j]; j++) {
								if(aligned_assem.q[j] && aligned_assem.q[j] != '-') {
									fprintf(consensus_out, "%c", aligned_assem.q[j]);
								} else {
									bias++;
								}
							}
							fprintf(consensus_out, "\n");
						}
					}
				}
			}
		}
		
	}
	
	/* Close files */
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	pclose(frag_out);
	if(matrix_out) {
		pclose(matrix_out);
	}
	for(i = 0; i < fileCount; i++) {
		sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", i);
		remove(outputfilename);
		outputfilename[file_len] = '\0';
	}
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
}

void runKMA_MEM(char *templatefilename, char *outputfilename, char *exePrev) {
	/* runKMA_MEM is a memory saving version of runKMA, */
	/* at the cost it chooses best templates based on kmers instead of alignment score. */
	
	/* open pipe */
	FILE *inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "File corruption\n");
		exit(1);
	}
	int i, j, k, t_i, tmp_template, tmp_tmp_template, file_len, best_read_score, template, template_stop, query_stop, s_len, gaps, bestHits, read_counter, start, end, Nhits, pos, asm_len, q_len, frags, aln_len, depthUpdate, coverScore, start_pos, end_pos, bias, matches, read_score, bestTemplate, template_tot_ulen, t_len, rc_flag;
	int tmp_start, tmp_end, *matched_templates, *best_start_pos, *best_end_pos, *bestTemplates, *qseq_int, *qseq_int_r, stats[4], *w_scores, fragCount, fileCount, maxFrag;
	double etta, tmp_score, score, bestScore, bestBaseScore, E_depth, tmp_Edepth, depth, q_value, p_value, expected, cover, bestNum;
	char bestNuc;
	char *outZipped, *qseq, *qseq_r, *tmp_string, *line, *line_sep, *start_s, *end_s, *template_s, *template_string, *template_acc, *header;
	struct aln aligned, gap_align;
	struct assem aligned_assem;
	struct frag **alignFrags, *alignFrag;
	FILE *align_in, *res_out, *frag_out, *alignment_out, *consensus_out, *frag_out_raw, *frag_in_raw, *fragmentIN, *matrix_out, **template_fragments;
	int line_size = 5000 * sizeof(char);
	load_DBs_KMA(templatefilename);
	/*templatefilename = realloc(templatefilename, (strlen(templatefilename) + strlen(".align.b") + 1) * sizeof(char));
	if(!templatefilename) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}*/
	strcat(templatefilename, ".align.b");
	templatefilename[strlen(templatefilename) + strlen(".align.b")] = '\0';
	align_in = fopen(templatefilename, "rb");
	if(align_in == NULL) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist.\n");
		exit(-1);
	}
	
	file_len = strlen(outputfilename);
	/*outputfilename = realloc(outputfilename, (file_len + 13) * sizeof(char));
	if(!outputfilename) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}*/
	outputfilename[file_len] = '\0';
	outZipped = malloc((strlen("gunzip -c .frag_raw.gz") + file_len + 8) * sizeof(char));
	if(!outZipped) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	time_t t0, t1;
	line = malloc(line_size);
	header = malloc(1024 * sizeof(char));
	if(!line || !header) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	template_fragments = calloc(DB_size, sizeof(FILE*));
	qseq_int = calloc(delta, sizeof(int));
	qseq_int_r = calloc(delta, sizeof(int));
	qseq = malloc((delta + 1) * sizeof(char));
	qseq_r = malloc((delta + 1) * sizeof(char));
	aligned.t = malloc(2*(delta + 1) * sizeof(char));
	aligned.s = malloc(2*(delta + 1) * sizeof(char));
	aligned.q = malloc(2*(delta + 1) * sizeof(char));
	gap_align.t = malloc(2*(delta + 1) * sizeof(char));
	gap_align.s = malloc(2*(delta + 1) * sizeof(char));
	gap_align.q = malloc(2*(delta + 1) * sizeof(char));
	if(!aligned.t || !aligned.s || !aligned.q || !gap_align.t || !gap_align.s || !gap_align.q || !qseq || !qseq_r || !qseq_int || !qseq_int_r) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	if(SW) {
		SW_s = delta + 1;
		D = malloc(SW_s * sizeof(int*));
		Q = malloc(SW_s * sizeof(int*));
		P = malloc(SW_s * sizeof(int*));
		E = malloc(SW_s * sizeof(int*));
		if(!D || !Q || !P || !E) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(i = 0; i < SW_s; i++) {
			D[i] = malloc(SW_s * sizeof(int));
			Q[i] = malloc(SW_s * sizeof(int));
			P[i] = malloc(SW_s * sizeof(int));
			E[i] = malloc(SW_s * sizeof(int));
			if(!D[i] || !Q[i] || !P[i] || !E[i]) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		convertNum = malloc(kmersize * sizeof(long unsigned));
		if(!convertNum) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		convertNum[kmersize - 1] = 1;
		for(i = kmersize - 2; i >= 0; i--) {
			convertNum[i] = convertNum[i+1] * 5;
		}
	}
	/* etta = small value to avoid zero-divisio */
	etta = 1.0e-6;
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		outputfilename[file_len + 4] = '\0';
		res_out = fopen(outputfilename, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".frag.gz");
		outputfilename[file_len + 8] = '\0';
		sprintf(outZipped, "gzip -c > %s", outputfilename);
		frag_out = popen(outZipped, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".aln");
		outputfilename[file_len + 4] = '\0';
		alignment_out = fopen(outputfilename, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".fsa");
		outputfilename[file_len + 4] = '\0';
		consensus_out = fopen(outputfilename, "w");
		outputfilename[file_len] = '\0';
		strcat(outputfilename, ".frag_raw.gz");
		outputfilename[file_len + 12] = '\0';
		sprintf(outZipped, "gzip -c > %s", outputfilename);
		frag_out_raw = popen(outZipped, "w");
		outputfilename[file_len] = '\0';
		if(print_matrix) {
			strcat(outputfilename, ".mat.gz");
			outputfilename[file_len + 7] = '\0';
			sprintf(outZipped, "gzip -c > %s", outputfilename);
			matrix_out = popen(outZipped, "w");
			if(!matrix_out) {
				fprintf(stderr, "File coruption: %s\n", outputfilename);
				exit(1);
			}
			outputfilename[file_len] = '\0';
		} else {
			matrix_out = 0;
		}
	} else {
		fprintf(stderr, "# No output file specified!\n");
		exit( 1 );
	}
	if(!res_out || !frag_out || !consensus_out || !frag_out_raw || !alignment_out) {
		fprintf(stderr, "Could not create outputfiles\n");
		exit(1);
	}
	
	fprintf(stderr, "# Collecting k-mer scores.\n");
	t0 = clock();
	
	/* Get alignments */
	read_counter = 0;
	matched_templates = malloc(2 * (DB_size + 1) * sizeof(int));
	best_start_pos = malloc(2 * DB_size * sizeof(int));
	best_end_pos = malloc(2 * DB_size * sizeof(int));
	bestTemplates = malloc(2 * DB_size * sizeof(int));
	if(!matched_templates || !best_start_pos || !best_end_pos || !bestTemplates) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	while(fread(&q_len, sizeof(int), 1, inputfile)) {
		if(q_len >= delta) {
			free(qseq_int);
			free(qseq_int_r);
			delta = 2 * q_len;
			qseq_int = malloc(delta * sizeof(int));
			qseq_int_r = malloc(delta * sizeof(int));
			free(aligned.t);
			free(aligned.s);
			free(aligned.q);
			aligned.t = malloc(2*(delta + 1) * sizeof(char));
			aligned.s = malloc(2*(delta + 1) * sizeof(char));
			aligned.q = malloc(2*(delta + 1) * sizeof(char));
			free(gap_align.t);
			free(gap_align.s);
			free(gap_align.q);
			gap_align.t = malloc(2*(delta + 1) * sizeof(char));
			gap_align.s = malloc(2*(delta + 1) * sizeof(char));
			gap_align.q = malloc(2*(delta + 1) * sizeof(char));
			free(qseq);
			free(qseq_r);
			qseq = malloc((delta + 1) * sizeof(char));
			qseq_r = malloc((delta + 1) * sizeof(char));
			if(!aligned.t || !aligned.s || !aligned.q || !gap_align.t || !gap_align.s || !gap_align.q || !qseq || !qseq_r || !qseq_int || !qseq_int_r) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
		}
		fread(qseq_int, q_len * sizeof(int), 1, inputfile);
		for(i = 0; i < q_len; i++) {
			qseq[i] = qseq_int[i] + '0';
		}
		qseq[q_len] = '\0';
		
		/* score */
		fread(&rc_flag, sizeof(int), 1, inputfile);
		if(rc_flag < 0) {
			bestScore = (-1) * rc_flag;
		} else {
			bestScore = rc_flag;
		}
		
		/* Get number of matched templates */
		fread(&bestHits, sizeof(int), 1, inputfile);
		fread(matched_templates, (bestHits + 1) * sizeof(int), 1, inputfile);
		
		/* get header */
		fread(&i, sizeof(int), 1, inputfile);
		fread(header, i, 1, inputfile);
		
		for(i = 0; i < bestHits; i++) {
			best_end_pos[i] = q_len;
			bestTemplates[i] = matched_templates[i + 1];
		}
		update_Scores(qseq, bestHits, bestScore, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		
	}
	pclose(inputfile);
	pclose(frag_out_raw);
	t1 = clock();
	fprintf(stderr, "#\n# Time for score collecting:\t%d s.\n", (int) difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select k-mer alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(int));
	if(!alignFrags || !w_scores) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	outputfilename[file_len] = '\0';
	sprintf(outZipped, "gunzip -c %s.frag_raw.gz", outputfilename);
	frag_in_raw = popen(outZipped, "r");
	if(!frag_in_raw) {
		fprintf(stderr, "File coruption: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[file_len] = '\0';
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	while( ! feof(frag_in_raw) && frag_in_raw) {
		qseq = fget_to(qseq, &delta, frag_in_raw, '\t');
		if(qseq[0] != '\0') {
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			bestHits = atoi(line);
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			score = atoi(line);
			
			for(i = 0; i < bestHits - 1; i++) {
				line = fget_to(line, &line_size, frag_in_raw, ',');
				best_start_pos[i] = atoi(line);
			}
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			best_start_pos[bestHits - 1] = atoi(line);
			for(i = 0; i < bestHits - 1; i++) {
				line = fget_to(line, &line_size, frag_in_raw, ',');
				best_end_pos[i] = atoi(line);
			}
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			best_end_pos[bestHits - 1] = atoi(line);
			for(i = 0; i < bestHits - 1; i++) {
				line = fget_to(line, &line_size, frag_in_raw, ',');
				bestTemplates[i] = atoi(line);
			}
			line = fget_to(line, &line_size, frag_in_raw, '\t');
			bestTemplates[bestHits - 1] = atoi(line);
			
			fgets(header, 1024, frag_in_raw);
			header[strlen(header) - 1] = '\0';
			
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				start = 0;
				end = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; i++) {
					tmp_tmp_template = bestTemplates[i];
					tmp_start = best_start_pos[i];
					tmp_end = best_end_pos[i];
					if(tmp_tmp_template < 0) {
						tmp_template = (-1) * tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = alignment_scores[tmp_template][0] / (template_lengths[tmp_template] - kmersize + 1);
					if(tmp_score > bestScore) {
					//if(alignment_scores[tmp_template][0] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template][0];
						bestScore = tmp_score;
						bestNum = alignment_scores[tmp_template][1];
						start = tmp_start;
						end = tmp_end;
					//} else if(alignment_scores[tmp_template][0] == best_read_score) {
					} else if(tmp_score == bestScore) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template][0] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template][0];
							bestScore = tmp_score;
							bestNum = alignment_scores[tmp_template][1];
							start = tmp_start;
							end = tmp_end;
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template][1] > bestNum) {
						} else if(alignment_scores[tmp_template][0] == best_read_score && alignment_scores[tmp_template][1] > bestNum) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template][0];
							bestScore = tmp_score;
							bestNum = alignment_scores[tmp_template][1];
							start = tmp_start;
							end = tmp_end;
						}
					}
				}
				/* reverse complement seq */
				if(bestTemplate < 0) {
					bestTemplate *= (-1);
					strrc(qseq, strlen(qseq));
				}
				w_scores[bestTemplate] += score;
				
				/* load frag info */
				if(alignFrags[bestTemplate] == 0) {
					alignFrags[bestTemplate] = malloc(sizeof(struct frag));
					alignFrag = alignFrags[bestTemplate];
					if(!alignFrag) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					alignFrag->qseq = strdup(qseq);
					alignFrag->header = strdup(header);
					alignFrag->score = (int) score;
					alignFrag->start = start;
					alignFrag->end = end;
					alignFrag->next = 0;
					if(!alignFrag->qseq || !alignFrag->header) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
				} else {
					alignFrag = malloc(sizeof(struct frag));
					if(!alignFrag) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					alignFrag->qseq = strdup(qseq);
					alignFrag->header = strdup(header);
					alignFrag->score = (int) score;
					alignFrag->start = start;
					alignFrag->end = end;
					alignFrag->next = alignFrags[bestTemplate];
					alignFrags[bestTemplate] = alignFrag;
					if(!alignFrag->qseq || !alignFrag->header) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
				}
				
				fragCount++;
				if(fragCount >= maxFrag) {
					sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
					template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
					if(!template_fragments[fileCount]) {
						fprintf(stderr, "File coruption: %s\n", outputfilename);
						exit(1);
					}
					outputfilename[file_len] = '\0';
					fileCount++;
					fragCount = 0;
					/* control fileamount */
					if(fileCount >= DB_size) {
						template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
						if(!template_fragments) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
				/* reverse complement seq */
				if(bestTemplate < 0) {
					bestTemplate *= (-1);
					strrc(qseq, strlen(qseq));
				}
				w_scores[bestTemplate] += score;
				
				/* load frag info */
				if(alignFrags[bestTemplate] == 0) {
					alignFrags[bestTemplate] = malloc(sizeof(struct frag));
					alignFrag = alignFrags[bestTemplate];
					if(!alignFrag) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					alignFrag->qseq = strdup(qseq);
					alignFrag->header = strdup(header);
					alignFrag->score = (int) score;
					alignFrag->start = start;
					alignFrag->end = end;
					alignFrag->next = 0;
					if(!alignFrag->qseq || !alignFrag->header) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
				} else {
					alignFrag = malloc(sizeof(struct frag));
					if(!alignFrag) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					alignFrag->qseq = strdup(qseq);
					alignFrag->header = strdup(header);
					alignFrag->score = (int) score;
					alignFrag->start = start;
					alignFrag->end = end;
					alignFrag->next = alignFrags[bestTemplate];
					alignFrags[bestTemplate] = alignFrag;
					if(!alignFrag->qseq || !alignFrag->header) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
				}
				
				fragCount++;
				if(fragCount >= maxFrag) {
					sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
					template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
					if(!template_fragments[fileCount]) {
						fprintf(stderr, "File coruption: %s\n", outputfilename);
						exit(1);
					}
					outputfilename[file_len] = '\0';
					fileCount++;
					fragCount = 0;
					/* control fileamount */
					if(fileCount >= DB_size) {
						template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
						if(!template_fragments) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
						template_fragments[fileCount] = 0;
					}
				}
			}
		}
	}
	/*
	free(matched_templates);
	free(best_start_pos);
	free(best_end_pos);
	free(bestTemplates);
	free(header);
	*/
	pclose(frag_in_raw);
	
	sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
	template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
	if(!template_fragments[fileCount]) {
		fprintf(stderr, "File coruption: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[file_len] = '\0';
	fileCount++;
	fragCount = 0;
	free(alignFrags);
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%d s.\n", (int) difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate length\tCoverage\tDepth\tq_value\tp_value\tClass\tPhenotype\tPubMed\n");
	
	/* Get expected values */
	template_tot_ulen = 0;
	Nhits = 0;
	for(i = 0; i < DB_size; i++) {
		Nhits += w_scores[i];
		template_tot_ulen += template_lengths[i];
	}
	/* Do local assemblies of fragments mapping to the same template */
	aligned_assem.t = NULL;
	aligned_assem.s = NULL;
	aligned_assem.q = NULL;
	for(template = 0; template < DB_size; template++) {
		if(w_scores[template] > 0) {
			// Check p_value
			score = w_scores[template];
			t_len = template_lengths[template];
			expected = (Nhits - score) * (1.0 * t_len) / (template_tot_ulen - t_len + etta);
			q_value = pow(score - expected, 2) / (expected + score + etta);
			p_value  = p_chisqr(q_value);
			
			if(p_value <= evalue && (score > expected || evalue == 1.0)) {
				/* load DB */
				(*alignLoadPtr)(template, align_in, 0, SEEK_CUR);
				
				/* Do assembly */
				free(aligned_assem.t);
				free(aligned_assem.s);
				free(aligned_assem.q);
				aligned_assem.t = NULL;
				aligned_assem.s = NULL;
				aligned_assem.q = NULL;
				aligned_assem = (*assemblyPtr)(template, template_fragments, fileCount, frag_out, matrix_out, outputfilename, qseq_int, aligned, gap_align);
				
				coverScore = aligned_assem.cover;
				depth = aligned_assem.depth;
				
				/* Create output statistics */
				depth /= t_len;
				cover = 100.0 * coverScore / t_len;
				
				/* Seperate template and accession number */
				template_acc = (char*)(template_names[template]);
				
				/* Output significant result */
				if(cover >= ID_t) {
					fprintf(res_out, "%-12s\t%8d\t%8d\t%8d\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_acc, (int) score, (int) expected, t_len, cover, depth, q_value, p_value);//, templates_descriptions[template]);
					
					/* print alignment */
					aln_len = strlen(aligned_assem.s);
					fprintf(alignment_out, "# %s\n", template_names[template]);
					for(i = 0; i < aln_len; i += 60) {
						fprintf(alignment_out, "%-10s\t", "template:");
						for(j = i; j < i + 60 && aligned_assem.t[j] != '\0'; j++) {
							fprintf(alignment_out, "%c", aligned_assem.t[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "");
						for(j = i; j < i + 60 && aligned_assem.s[j] != '\0'; j++) {
							fprintf(alignment_out, "%c", aligned_assem.s[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "query:");
						for(j = i; j < i + 60 && aligned_assem.q[j] != '\0'; j++) {
							fprintf(alignment_out, "%c", aligned_assem.q[j]);
						}
						fprintf(alignment_out, "\n\n");
					}
					
					/* Print consensus */
					fprintf(consensus_out, ">%s\n", template_names[template]);
					if(ref_fsa) {
						for(i = 0; i < aln_len; i += 60) {
							for(j = i; j < i + 60 && aligned_assem.q[j]; j++) {
								if(aligned_assem.q[j]) {
									if(aligned_assem.q[j] != '-') {
										fprintf(consensus_out, "%c", aligned_assem.q[j]);
									} else {
										fprintf(consensus_out, "%c", 'n');
									}
								}
							}
							fprintf(consensus_out, "\n");
						}
					} else {
						bias = 0;
						for(i = 0; i < aln_len; i += 60 + bias) {
							bias = 0;
							for(j = i; j < i + 60 + bias && aligned_assem.q[j]; j++) {
								if(aligned_assem.q[j] && aligned_assem.q[j] != '-') {
									fprintf(consensus_out, "%c", aligned_assem.q[j]);
								} else {
									bias++;
								}
							}
							fprintf(consensus_out, "\n");
						}
					}
				}
				/* destroy this DB index */
				free(template_seqs[template]);
				template_seqs[template] = NULL;
				/* here */
				(*destroyPtr)(template);
			} else {
				fseek(align_in, template_lengths[template], SEEK_CUR);
			}
		} else {
			fseek(align_in, template_lengths[template], SEEK_CUR);
		}	
	}
	
	/* Close files */
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	pclose(frag_out);
	if(matrix_out) {
		pclose(matrix_out);
	}
	for(i = 0; i < fileCount; i++) {
		sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", i);
		remove(outputfilename);
		outputfilename[file_len] = '\0';
	}
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# KMA-1.0 mapps raw reads to a template database, for optimal performance it is designed to use 3 threads.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-i\t\tInput/query file name\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t16\n");
	fprintf(helpOut, "#\t-e\t\tevalue\t\t\t\t0.05\n");
	fprintf(helpOut, "#\t-delta\t\tAllocation size for sequences\t512\n");
	fprintf(helpOut, "#\t-mem_mode\tUse kmers to choose best\n#\t\t\ttemplate, and save memory\tFalse\n");
	fprintf(helpOut, "#\t-ex_mode\tSearh kmers exhaustively\tFalse\n");
	fprintf(helpOut, "#\t-deCon\t\tRemove contamination\t\tFalse\n");
	fprintf(helpOut, "#\t-dense\t\tDo not allow insertions\n#\t\t\tin assembly\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ref_fsa\tConsensus sequnce will\n#\t\t\thave \"n\" instead of gaps\tFalse\n");
	fprintf(helpOut, "#\t-matrix\t\tPrint assembly matrix\t\tFalse\n");
	fprintf(helpOut, "#\t-mp\t\tMinimum phred score\t\t30\n");
	fprintf(helpOut, "#\t-5p\t\tCut a constant number of\n#\t\t\tnucleotides from the 5 prime.\t0\n");
	fprintf(helpOut, "#\t-CS\t\tChain size\t\t\t8 MB\n");
	fprintf(helpOut, "#\t-MS\t\tMax chain size\t\t\t14 GB\n");
	fprintf(helpOut, "#\t-Sparse\t\tRun KmerFinder\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ID\t\tMinimum ID\t\t\t1.0%%\n");
	fprintf(helpOut, "#\t-ss\t\tSparse sorting (q,c,d)\t\tq\n");
	fprintf(helpOut, "#\t-NW\t\tUse Needleman-Wunsch\t\tFalse\n");
	fprintf(helpOut, "#\t-shm\t\tUse shared DB made by kma_shm\tFalse\n");
	fprintf(helpOut, "#\t-1t1\t\tSkip HMM\t\t\tFalse\n");
	fprintf(helpOut, "#\t-mrs\t\tMinimum alignment score score,\n#\t\t\tnormalized to alignment length\t0.0\n");
	fprintf(helpOut, "#\t-reward\t\tScore for match\t\t\t1\n");
	fprintf(helpOut, "#\t-penalty\tPenalty for mismatch\t\t-2\n");
	fprintf(helpOut, "#\t-gapopen\tPenalty for gap opening\t\t-3\n");
	fprintf(helpOut, "#\t-gapextend\tPenalty for gap extension\t-1\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int i, j, args, stop, exe_len, minPhred, fiveClip, sparse_run, fileCounter, step1, step2, step3;
	time_t t0, t1;
	char *exeBasic, **inputfiles, *outputfilename, *templatefilename, ss;
	/* SET DEFAULTS */
	assemblyPtr = &assemble_KMA;
	minPhred = 30;
	fiveClip = 0;
	sparse_run = 0;
	fileCounter = 1;
	outputfilename = 0;
	templatefilename = 0;
	print_matrix = 0;
	ref_fsa = 0;
	contamination = -1;
	kmersize = 16;
	evalue = 0.05;
	delta = 2048;
	exhaustive = 0;
	MAX_SIZE = 14 * (long unsigned)(pow(2, 30) + 0.5) / 8; //14 = #GB, 2^30 = GB, 8 = sizeof pointer
	step1 = 0;
	step2 = 0;
	step3 = 0;
	ID_t = 1.0;
	ss = 'q';
	deCon = 0;
	mem_mode = 0;
	M = 1;
	MM = -2;
	W1 = -3;
	U = -1;
	kmerScan = &save_kmers_HMM;
	ankerPtr = &ankerAndClean;
	alignLoadPtr = &alignLoad_fly;
	destroyPtr = &alignClean;
	deConPtr = ankerPtr;
	inputfiles = malloc(sizeof(char*));
	inputfiles[0] = malloc(3 * sizeof(char));
	if(!inputfiles || !inputfiles[0]) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	strncpy(inputfiles[0], "--", 2);
	inputfiles[0][2] = '\0';
	SW = 0;
	aligner = &KMA;
	shm = 0;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			stop = 0;
			args++;
			if(args < argc && (strncmp(argv[args], "-", 1) != 0 || strncmp(argv[args], "--", 2) == 0)) {
				inputfiles[0] = realloc(inputfiles[0], (strlen(argv[args]) + 1) * sizeof(char));
				if(!inputfiles[0]) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(inputfiles[0], argv[args]);
				inputfiles[0][strlen(argv[args])] = '\0';
				args++;
			}
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strncmp(argv[args], "--", 2) == 0) {
					fileCounter++;
					inputfiles = realloc(inputfiles, fileCounter * sizeof(char*));
					if(!inputfiles) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					inputfiles[fileCounter - 1] = malloc((strlen(argv[args]) + 1) * sizeof(char));
					if(!inputfiles[fileCounter - 1]) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					strcpy(inputfiles[fileCounter - 1], argv[args]);
					inputfiles[fileCounter - 1][strlen(argv[args])] = '\0';
					args++;
				} else {
					stop = 1;
				}
			}
			args--;
		} else if(strcmp(argv[args], "-delta") == 0) {
			args++;
			if(args < argc) {
				delta = atoi(argv[args]);
				if(delta == 0) {
					fprintf(stderr, "# Invalid delta specified.\n");
					exit(-1);
				}
			}
		} else if(strcmp(argv[args], "-o") == 0) {
			args++;
			if(args < argc) {
				outputfilename = malloc((strlen(argv[args]) + 64) * sizeof(char));
				if(!outputfilename) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(outputfilename, argv[args]);
				outputfilename[strlen(outputfilename)] = '\0';
			}
		} else if(strcmp(argv[args], "-deCon") == 0) {
			deCon = 1;
			deConPtr = &deConAnkers;
		} else if(strcmp(argv[args], "-NW") == 0 || strcmp(argv[args], "-SW") == 0) {
			SW = 1;
			aligner = &KMA_SW;
			alignLoadPtr = &alignLoad_fly_SW;
			destroyPtr = &alignClean_SW;
		} else if(strcmp(argv[args], "-shm") == 0) {
			shm = 1;
		} else if(strcmp(argv[args], "-step1") == 0) {
			step1 = 1;
		} else if(strcmp(argv[args], "-step2") == 0) {
			step2 = 1;
		} else if(strcmp(argv[args], "-step3") == 0) {
			step3 = 1;
		} else if(strcmp(argv[args], "-mem_mode") == 0) {
			mem_mode = 1;
			ankerPtr = &ankerAndClean_MEM;
		} else if(strcmp(argv[args], "-ex_mode") == 0) {
			exhaustive = 1;
		} else if(strcmp(argv[args], "-t_db") == 0) {
			args++;
			if(args < argc) {
				templatefilename = malloc((strlen(argv[args]) + 64) * sizeof(char));
				if(!templatefilename) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(templatefilename, argv[args]);
				templatefilename[strlen(argv[args])] = '\0';
			}
		} else if(strcmp(argv[args], "-k") == 0) {
			args++;
			if(args < argc) {
				kmersize = atoi(argv[args]);
				if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				}
			}
		} else if(strcmp(argv[args], "-mp") == 0) {
			args++;
			if(args < argc) {
				minPhred = atoi(argv[args]);	
			}
		} else if(strcmp(argv[args], "-5p") == 0) {
			args++;
			if(args < argc) {
				fiveClip = atoi(argv[args]);	
			}
		} else if(strcmp(argv[args], "-dense") == 0) {
			assemblyPtr = &assemble_KMA_dense;
		} else if(strcmp(argv[args], "-matrix") == 0) {
			print_matrix = 1;
		} else if(strcmp(argv[args], "-ref_fsa") == 0) {
			ref_fsa = 1;
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
		} else if(strcmp(argv[args], "-1t1") == 0) {
			kmerScan = &save_kmers;
			ankerPtr = &ankerAndClean_1t1;
		} else if(strcmp(argv[args], "-ss") == 0) {
			args++;
			if(args < argc) {
				if(argv[args][0] == 'q') {
					ss = 'q';
				} else if(argv[args][0] == 'c') {
					ss = 'c';
				} else if(argv[args][0] == 'd') {
					ss = 'd';
				} else {
					fprintf(stderr, "# Invalid argument parsed to option: \"-ss\", using default.\n");
				}
			}
		} else if(strcmp(argv[args], "-e") == 0) {
			args++;
			if(args < argc) {
				evalue = atof(argv[args]);
			}
		} else if(strcmp(argv[args], "-ID") == 0) {
			args++;
			if(args < argc) {
				ID_t = atof(argv[args]);
			}
		} else if(strcmp(argv[args], "-CS") == 0) {
			args++;
			if(args < argc) {
				INITIAL_SIZE = atoi(argv[args]) * 1048576 / 8;
				if(INITIAL_SIZE == 0) {
					fprintf(stderr, "# Invalid Chain Size parsed, using default\n");
					INITIAL_SIZE = 1048576;
				}
			}
		} else if(strcmp(argv[args], "-MS") == 0) {
			args++;
			if(args < argc) {
				MAX_SIZE = atoi(argv[args]) * (long unsigned)(pow(2, 30) + 0.5) / 8;
				if(MAX_SIZE == 0) {
					fprintf(stderr, "# Invalid Max Size parsed, using default\n");
					MAX_SIZE = 14 * (long unsigned)(pow(2, 30) + 0.5) / 8;
				}	
			}
		} else if(strcmp(argv[args], "-mrs") == 0) {
			args++;
			if(args < argc) {
				scoreT = atof(argv[args]);
			}
		} else if(strcmp(argv[args], "-reward") == 0) {
			args++;
			if(args < argc) {
				M = atoi(argv[args]);
			}
		} else if(strcmp(argv[args], "-penalty") == 0) {
			args++;
			if(args < argc) {
				MM = atoi(argv[args]);
			}
		} else if(strcmp(argv[args], "-gapopen") == 0) {
			args++;
			if(args < argc) {
				W1 = atoi(argv[args]);
			}
		} else if(strcmp(argv[args], "-gapextend") == 0) {
			args++;
			if(args < argc) {
				U = atoi(argv[args]);
			}
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(-1);
		}
		args++;
	}
	if(outputfilename == 0 || templatefilename == 0 || (step1 && fileCounter == 0)) {
		fprintf(stderr, "# Too few arguments handed\n");
		fprintf(stderr, "# Printing help message:\n");
		helpMessage(-1);
	}
	
	/* set scoring matrix */
	for(i = 0; i < 4; i++) {
		for(j = 0; j < 4; j++) {
			d[i][j] = MM;
		}
		d[i][i] = M;
	}
	for(i = 0; i < 5; i++) {
		d[4][i] = U;
		d[i][4] = U;
	}
	d[4][4] = 0;
	
	
	if(step1) {
		t0 = clock();
		if(sparse_run) {
			run_input_sparse(inputfiles, fileCounter, minPhred, fiveClip);
		} else {
			run_input(inputfiles, fileCounter, minPhred, fiveClip);
		}
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for converting query: %d s.\n#\n", (int) difftime(t1, t0) / 1000000);
	} else if(step2) {
		exe_len = strlen("-step1");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step1");
		exeBasic[exe_len] = '\0';
		
		save_kmers_batch(templatefilename, exeBasic);
	} else if(sparse_run) {
		exe_len = strlen("-step1");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step1");
		exeBasic[exe_len] = '\0';
		
		save_kmers_sparse_batch(templatefilename, outputfilename, exeBasic, ss);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
		fprintf(stdout, "DONE");
	} else if(mem_mode) {
		exe_len = strlen("-step2");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step2");
		exeBasic[exe_len] = '\0';
		runKMA_MEM(templatefilename, outputfilename, exeBasic);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
		fprintf(stdout, "DONE");
	} else {
		exe_len = strlen("-step2");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step2");
		exeBasic[exe_len] = '\0';
		runKMA(templatefilename, outputfilename, exeBasic);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
		fprintf(stdout, "DONE");
	}
	
	exit(0);
}
