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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <unistd.h>
#include <errno.h>
#define getNuc(Comp,pos)((Comp[pos >> 5] << ((pos & 31) << 1)) >> 62)

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
	int q_len;
	int bestHits;
	int score;
	int start;
	int end;
	char *header;
	struct frag *next;
};

struct hashTable {
	long unsigned key;
	int *value;
	struct hashTable *next;
};

struct hashMap_index {
	unsigned len; // seqlen
	unsigned size; // size of index
	int *index; // k-mer posititions in seq
	long unsigned *seq; // 2-bit sequence
};

struct hashTable_kmers {
	long unsigned key;
	int value;
	struct hashTable_kmers *next;
};

struct hashMap_kmers {
	unsigned size;
	unsigned n;
	struct hashTable_kmers **table;
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

struct compDNA {
	int seqlen;
	int size;
	int complen;
	long unsigned *seq;
	int *N;
};

struct compKmers {
	int n;
	int size;
	long unsigned *kmers;
};

struct hashMapKMA {
	unsigned kmersize;		// k
	long unsigned size;			// size of DB
	unsigned n;				// k-mers stored
	unsigned null_index;	// null value
	unsigned seqsize;		// size of seq
	unsigned v_index;		// size of values
	unsigned prefix_len;	// prefix length
	long unsigned prefix;	// prefix
	unsigned *exist;		// size long
	long unsigned *seq;		// compressed sequence of k-mers
	int *values;			// compressed values
	unsigned *key_index	;	// Relative
	unsigned *value_index;	// Relative
};

struct alnScore {
	int score;
	int len;
	int pos;
	int gaps;
};

struct diskOffsets {
	long exist;			//8 * sizeof(unsigned) + sizeof(long unsigned)
	long seq;			//exist + templates->size * sizeof(unsigned)
	long values;		//seq + templates->seqsize * sizeof(long unsigned)
	long key_index;		//values + templates->v_index * sizeof(unsigned)
	long value_index;	//key_index + templates->null_index * sizeof(unsigned)
	FILE *file;			// filedes of templates
};

/*
 	GLOBAL VARIABLES
*/
struct hashMapKMA *templates;
struct hashMap_index **templates_index;
struct diskOffsets *templates_offsets;
char **template_names, *to2Bit;
int *template_lengths, *template_ulengths, *bestTemplates, *bestTemplates_r;
int *regionTemplates, *Score, *Score_r, *tmpNs, exhaustive, minLen, ref_fsa;
int mincoverage, delta, deCon, contamination, print_matrix, print_all, shm;
int DB_size, shifter, r_shifter, **VF_scores, **VR_scores, *valuesFile, diskDB;
unsigned kmersize, *alignment_scores, *uniq_alignment_scores;
long unsigned mask;
double evalue, ID_t, scoreT, HMM_param[8];
int **D, **Q, **P, **E, W1, U, d[5][5], M, MM, NW_s;
int com_bases[] = {3, 2, 1, 0, 4};
char bases[] = "ACGTN-";

/*
	Set function pointers 
*/
void (*printPtr)(int*, struct compDNA *, int, char*, int);
void (*deConPrintPtr)(int*, struct compDNA *, int, char*, int);
void (*ankerPtr)(struct compDNA *, int, int, int, int, char*, int);
void (*assemblyPtr)(struct assem *, int, FILE**, int, FILE*, FILE*, char*, struct aln *, struct aln *, char *, char *);
void (*destroyPtr)(int);
struct hashMap_index * (*alignLoadPtr)(FILE*, FILE*, int, long unsigned, long unsigned);
int * (*hashMap_get)(long unsigned);
void (*kmerScan)(struct compDNA*, struct compDNA*, char*, int);
void (*printFsa_ptr)(struct qseqs *, struct qseqs *, struct compDNA *);

/*
 FUNCTIONS
*/

/* BASIC FUNCTIONS */
long unsigned getKmer(long unsigned *compressor, unsigned pos) {
	
	unsigned cPos, iPos;
	cPos = pos >> 5;
	iPos = (pos & 31) << 1;
	
	return (iPos <= shifter) ? ((compressor[cPos] << iPos) >> shifter) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifter);
}

long unsigned makeKmer(const char *qseq, unsigned pos, unsigned size) {
	
	long unsigned key = qseq[pos];
	
	size += pos;
	for(pos++; pos < size; pos++) {
		key = (key << 2) | qseq[pos];
	}
	
	return key;
}

int charpos(const char *src, char target, int start, int len) {
	
	int i;
	
	for(i = start; i < len; i++) {
		if(src[i] == target) {
			return i;
		}
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

int strrpos(const char* str1, const char* str2) {
	char *strp;
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

unsigned countChar(const char* str, char target) {
	
	char *strp;
	unsigned count;
	
	count = 0;
	for(strp = (char*)(str); *strp; strp++) {
		if(*strp == target) {
			count++;
		}
	}
	
	return count;
}

int intpos_bin(const int *str1, const int str2) {
	int i, pos, upLim, downLim;
	
	upLim = *str1;
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) / 2;
	while(0 < (upLim - downLim)) {
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

int chomp(char *string) {
	/* remove trailing whitespaces, and return length of string */
	int k = strlen(string) - 1;
	while(isspace(string[k])) {
		k--;
	}
	k++;
	string[k] = 0;
	return k;
}

void insert(char *dest, char src, int location, int dest_len) {
	int i;
	dest[dest_len + 1] = 0;
	for(i = dest_len; i > location; i--) {
		dest[i] = dest[i - 1];
	}
	dest[location] = src;
}

int replace_chars(char *dest, char src) {
	int i, j, bias, len;
	while(*dest == src) {
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

char * fget_line(char *line, int *line_size, int *l_len, FILE *file) {
	
	int i, grow;
	i = 0;
	*l_len = 0;
	grow = *line_size;
	
	while(!feof(file) && fgets((line + i), grow, file)) {
		*l_len = chomp(line);
		if(*l_len == (*line_size - 1)) { //realloc
			i = *l_len;
			*line_size += grow;
			line = realloc(line, *line_size);
			if(!line) {
				fprintf(stderr, "fget_line error.\n");
				exit(-1);
			}
		} else {
			return line;
		}
	}
	*line = 0;
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
	
	char *seq_ptr, *buff_ptr;
	int seqlen, seqsize, destpos, destbytes;
	
	/* get header */
	header->len = 0;
	while((seqlen = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		seqsize = seqlen - dest->pos;
		if(header->size < (seqsize + header->len)) {
			header->size = seqsize + header->len;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(header->seq + header->len, dest->buffer + dest->pos, seqsize);
		header->len += seqsize;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	seqsize = seqlen - dest->pos;
	if(header->size < (seqsize + header->len)) {
		header->size = seqsize + header->len;
		header->seq = realloc(header->seq, header->size);
		if(!header->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	strncpy(header->seq + header->len, dest->buffer + dest->pos, seqsize);
	header->len += seqsize;
	header->seq[header->len] = 0;
	dest->pos = seqlen + 1;
	
	/* get seq */
	seqlen = 0;
	seqsize = seq->size;
	seq_ptr = seq->seq;
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	destpos = dest->pos;
	destbytes = dest->bytes;
	buff_ptr = dest->buffer;
	
	while(buff_ptr[destpos] != '>') {
		/* accept char */
		seq_ptr[seqlen] = to2Bit[buff_ptr[destpos]];
		if(seq_ptr[seqlen] < 5) {
			seqlen++;
			if(seqlen == seqsize) {
				seq->size <<= 1;
				seq->seq = realloc(seq->seq, seq->size);
				if(!seq->seq) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				seqsize = seq->size;
				seq_ptr = seq->seq;
			}
		}
		destpos++;
		
		if(destpos == destbytes) {
			destpos = 0;
			if(!buffFileBuff(dest)) {
				break;
			}
			destbytes = dest->bytes;
		}
		
	}
	dest->pos = destpos;
	seq->len = seqlen;
	
	return 1;
}

int FileBuffgetFsaSeq(struct FileBuff *dest, struct qseqs *seq) {
	
	char *seq_ptr, *buff_ptr;
	int seqlen, seqsize, destpos, destbytes;
	
	/* skip header */
	while((dest->pos = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos++;
	
	/* get seq */
	seqlen = 0;
	seqsize = seq->size;
	seq_ptr = seq->seq;
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	destpos = dest->pos;
	destbytes = dest->bytes;
	buff_ptr = dest->buffer;
	
	while(buff_ptr[destpos] != '>') {
		/* accept char */
		seq_ptr[seqlen] = to2Bit[buff_ptr[destpos]];
		if(seq_ptr[seqlen] < 5) {
			seqlen++;
			if(seqlen == seqsize) {
				seq->size <<= 1;
				seq->seq = realloc(seq->seq, seq->size);
				if(!seq->seq) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				seqsize = seq->size;
				seq_ptr = seq->seq;
			}
		}
		destpos++;
		
		if(destpos == destbytes) {
			destpos = 0;
			if(!buffFileBuff(dest)) {
				break;
			}
			destbytes = dest->bytes;
		}
	}
	dest->pos = destpos;
	seq->len = seqlen;
	
	return 1;
}

int FileBuffgetFq(struct FileBuff *dest, struct qseqs *header, struct qseqs *seq, struct qseqs *qual) {
	
	
	int i, buff_end, seq_end;
	char *seq_ptr, *buff_ptr;
	
	/* get header */
	header->len = 0;
	while((seq_end = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		buff_end = seq_end - dest->pos;
		if(header->size < (buff_end + header->len)) {
			header->size = buff_end + header->len;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		strncpy(header->seq + header->len, dest->buffer + dest->pos, buff_end);
		header->len += buff_end;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	buff_end = seq_end - dest->pos;
	if(header->size < (buff_end + header->len)) {
		header->size = buff_end + header->len;
		header->seq = realloc(header->seq, header->size);
		if(!header->seq) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	strncpy(header->seq + header->len, dest->buffer + dest->pos, buff_end);
	header->len += buff_end;
	dest->pos = seq_end + 1;
	
	
	
	/* get seq */
	buff_ptr = (dest->buffer + dest->pos);
	seq_ptr = seq->seq;
	
	buff_end = dest->bytes - dest->pos;
	seq_end = seq->size;
	i = 0;
	seq->len = 0;
	if(!buff_end) {
		/* buff buffer and sync with seq*/
		if(!(buff_end = buffFileBuff(dest))) {
			return 0;
		}
		buff_ptr = dest->buffer;
	}
	
	while(buff_ptr[i] != '\n') {
		/* accept char */
		seq_ptr[i] = to2Bit[buff_ptr[i]];
		i++;
		
		/* check seq */
		if(i == seq_end) {
			seq_end += seq->size;
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			seq_ptr = seq->seq + seq->len;
		}
		
		/* check buffer */
		if(i == buff_end) {
			/* move seq offset */
			seq->len += i;
			seq_ptr += i;
			seq_end -= i;
			
			/* buff buffer and sync with seq*/
			if(!(buff_end = buffFileBuff(dest))) {
				return 0;
			}
			buff_ptr = dest->buffer;
			i = 0;
		}
	}
	seq->len += i;
	dest->pos += (i + 1);
	
	/* skip info */
	while((dest->pos = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos++;
	
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
	seq_end = dest->pos + seq->len;
	while(seq_end >= dest->bytes) {
		buff_end = dest->bytes - dest->pos;
		strncpy(qual->seq + qual->len, dest->buffer + dest->pos, buff_end);
		qual->len += buff_end;
		seq_end -= dest->bytes;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	buff_end = seq_end - dest->pos;
	strncpy(qual->seq + qual->len, dest->buffer + dest->pos, buff_end);
	qual->len = seq->len;
	dest->pos = seq_end + 1;
	
	return 1;
}

int FileBuffgetFqSeq2(struct FileBuff *dest, struct qseqs *seq, struct qseqs *qual) {
	
	int seqlen, seqsize, destpos, destbytes, chunk, increase;
	char *seq_ptr, *buff_ptr;
	
	/* skip header */
	while((dest->pos = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos++;
	
	/* get seq */
	seqlen = 0;
	seqsize = seq->size;
	seq_ptr = seq->seq;
	if(dest->pos == dest->bytes && !buffFileBuff(dest)) {
		return 0;
	}
	destpos = dest->pos;
	destbytes = dest->bytes;
	buff_ptr = dest->buffer;
	
	while(buff_ptr[destpos] != '\n') {
		/* accept char */
		seq_ptr[seqlen] = to2Bit[buff_ptr[destpos]];
		seqlen++;
		destpos++;
		
		if(seqlen == seqsize) {
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			seqsize = seq->size;
			seq_ptr = seq->seq;
		}
		
		if(destpos == destbytes) {
			if(!buffFileBuff(dest)) {
				return 0;
			}
			destpos = 0;
			destbytes = dest->bytes;
		}
		
	}
	dest->pos = destpos + 1;
	seq->len = seqlen;
	
	/* skip info */
	while((dest->pos = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos++;
	
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
	chunk = dest->pos + seq->len;
	while(chunk >= dest->bytes) {
		increase = dest->bytes - dest->pos;
		strncpy(qual->seq + qual->len, dest->buffer + dest->pos, increase);
		qual->len += increase;
		chunk -= dest->bytes;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	increase = chunk - dest->pos;
	strncpy(qual->seq + qual->len, dest->buffer + dest->pos, increase);
	qual->len = seq->len;
	dest->pos = chunk + 1;
	
	return 1;
}


int FileBuffgetFqSeq(struct FileBuff *dest, struct qseqs *seq, struct qseqs *qual) {
	
	int i, buff_end, seq_end;
	char *seq_ptr, *buff_ptr;
	
	/* skip header */
	while((dest->pos = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos++;
	
	
	/* get seq */
	buff_ptr = (dest->buffer + dest->pos);
	seq_ptr = seq->seq;
	
	buff_end = dest->bytes - dest->pos;
	seq_end = seq->size;
	i = 0;
	seq->len = 0;
	if(!buff_end) {
		/* buff buffer and sync with seq*/
		if(!(buff_end = buffFileBuff(dest))) {
			return 0;
		}
		buff_ptr = dest->buffer;
	}
	
	while(buff_ptr[i] != '\n') {
		/* accept char */
		seq_ptr[i] = to2Bit[buff_ptr[i]];
		i++;
		
		/* check seq */
		if(i == seq_end) {
			seq_end += seq->size;
			seq->size <<= 1;
			seq->seq = realloc(seq->seq, seq->size);
			if(!seq->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			seq_ptr = seq->seq + seq->len;
		}
		
		/* check buffer */
		if(i == buff_end) {
			/* move seq offset */
			seq->len += i;
			seq_ptr += i;
			seq_end -= i;
			
			/* buff buffer and sync with seq*/
			if(!(buff_end = buffFileBuff(dest))) {
				return 0;
			}
			buff_ptr = dest->buffer;
			i = 0;
		}
	}
	seq->len += i;
	dest->pos += (i + 1);
	
	/* skip info */
	while((dest->pos = chunkPos(dest->buffer, dest->pos, dest->bytes)) == dest->bytes) {
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	dest->pos++;
	
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
	seq_end = dest->pos + seq->len;
	while(seq_end >= dest->bytes) {
		strncpy(qual->seq + qual->len, dest->buffer + dest->pos, (buff_end = dest->bytes - dest->pos));
		qual->len += buff_end;
		seq_end -= dest->bytes;
		if(!buffFileBuff(dest)) {
			return 0;
		}
	}
	strncpy(qual->seq + qual->len, dest->buffer + dest->pos, seq_end - dest->pos);
	qual->len = seq->len;
	dest->pos = seq_end + 1;
	
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
	for(i = *out_Tem; i > 0; i--) {
		if(out_Tem[i] == contamination) {
			return i;
		}
	}
	return -1;
}

int find_contamination2(int *out_Tem, int contamination_s) {
	int i;
	for(i = *out_Tem; i > 0; i--) {
		if(out_Tem[i] == contamination_s) {
			return i;
		}
	}
	return -1;
}

void strrc(char *qseq, int q_len) {
	
	static char comp[6] = {3, 2, 1, 0, 4, 5};
	char carry;
	int i, j, seqlen;
	
	seqlen = q_len >> 1;
	
	for(i = 0, j = q_len - 1; i < seqlen; i++, j--) {
		carry = comp[qseq[i]];
		qseq[i] = comp[qseq[j]];
		qseq[j] = carry;
	}
	if(q_len & 1) {
		qseq[seqlen] = comp[qseq[seqlen]];
	}
	
}

void convertToNum(char *qseq, int q_len) {
	int i;
	for(i = 0; i < q_len; i++) {
		qseq[i] = to2Bit[qseq[i]];
	}
}

/*
	COMPRESSION FUNCTIONS
*/
void allocCompKmers(struct compKmers *compressor, int size) {
	
	compressor->n = 0;
	compressor->size = size;
	compressor->kmers = malloc(size * sizeof(long unsigned));
	if(!compressor->kmers) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
}

void reallocCompKmers(struct compKmers *compressor, int size) {
	
	compressor->kmers = realloc(compressor->kmers, size * sizeof(long unsigned));
	if(!compressor->kmers) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	compressor->size = size;
	
}

void allocComp(struct compDNA *compressor, int size) {
	
	compressor->seqlen = 0;
	if(size & 31) {
		compressor->size = (size >> 5) + 1;
		compressor->size <<= 5;
	} else {
		compressor->size = size;
	}
	
	compressor->seq = calloc(compressor->size >> 5, sizeof(long unsigned));
	compressor->N = malloc((compressor->size + 1) * sizeof(int));
	
	if(!compressor->seq || !compressor->N) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	compressor->N[0] = 0;
	
}

int pushCompKmers(struct compKmers *compressor, char *qseq, int kmersize) {
	
	int i;
	long unsigned key = 0;
	
	for(i = 0; i < kmersize; i++) {
		if(qseq[i] == 4) {
			return 0;
		} else {
			key = (key << 2) | qseq[i];
		}
	}
	
	compressor->kmers[compressor->n] = key;
	return 1;
}

void freeComp(struct compDNA *compressor) {
	
	compressor->seqlen = 0;
	compressor->complen = 0;
	compressor->size = 0;
	
	free(compressor->seq);
	free(compressor->N);
	
}

void resetComp(struct compDNA *compressor) {
	compressor->N[0] = 0;
	compressor->seqlen = 0;
	compressor->complen = 0;
}

void strtranslate(const char *qseq, char *trans) {
	char *strp;
	
	for(strp = (char*) qseq; *strp; strp++) {
		*strp = trans[*strp];
	}
}

int translateToKmersAndDump(long unsigned *Kmers, int n, int max, char *qseq, int seqlen, long unsigned prefix, int prefix_len) {
	
	int i, end, rc;
	long unsigned key;
	
	if(prefix_len) {
		for(rc = 0; rc < 2; rc++) {
			
			if(rc) {
				strrc(qseq, seqlen);
			}
			
			i = 0;
			while(i < seqlen) {
				end = charpos(qseq, 4, i, seqlen);
				if(end == -1) {
					end = seqlen;
				}
				if(i < end - kmersize - prefix_len) {
					key = makeKmer(qseq, i, prefix_len - 1);
					i += (prefix_len - 1);
					end -= kmersize;
				} else {
					i = end + 1;
				}
				
				while(i < end) {
					key = ((key << 2) | qseq[i]) & mask;
					i++;
					if(key == prefix) {
						Kmers[n] = makeKmer(qseq, i, kmersize);
						n++;
						if(n == max) {
							fwrite(Kmers, sizeof(long unsigned), n, stdout);
							n = 0;
						}
					}
				}
				i = end + kmersize + 1;
			}
		}
	} else {
		for(rc = 0; rc < 2; rc++) {
			if(rc) {
				strrc(qseq, seqlen);
			}
			
			i = 0;
			while(i < seqlen) {
				end = charpos(qseq, 4, i, seqlen);
				if(end == -1) {
					end = seqlen;
				}
				key = makeKmer(qseq, i, kmersize - 1);
				while(i < end) {
					key = ((key << 2) | qseq[i]) & mask;
					i++;
					Kmers[n] = key;
					n++;
					if(n == max) {
						fwrite(Kmers, sizeof(long unsigned), n, stdout);
						n = 0;
					}
				}
				i = end + kmersize + 1;
			}
		}
	}
	return n;
}

void compDNA(struct compDNA *compressor, char *seq, int seqlen) {
	
	int i, j, pos, end;
	
	compressor->seqlen = seqlen;
	if(seqlen & 31) {
		compressor->complen = (seqlen >> 5) + 1;
	} else {
		compressor->complen = seqlen >> 5;
	}
	
	
	for(i = 0; i < seqlen; i += 32) {
		end = (i + 32 < seqlen) ? i + 32 : seqlen;
		pos = i >> 5;
		for(j = i; j < end; j++) {
			if(seq[j] == 4) {
				compressor->seq[pos] <<= 2;
				compressor->N[0]++;
				compressor->N[compressor->N[0]] = j;
			} else {
				compressor->seq[pos] = (compressor->seq[pos] << 2) | seq[j];
			}
		}
	}
	if(seqlen & 31) {
		compressor->seq[pos] <<= (64 - ((seqlen & 31) << 1));
	}
}

void unCompDNA(struct compDNA *compressor, char *seq) {
	
	int i;
	
	/* get nucs */
	for(i = 0; i < compressor->seqlen; i++) {
		seq[i] = getNuc(compressor->seq, i);
	}
	
	/* get N's */
	for(i = 1; i <= compressor->N[0]; i++) {
		seq[compressor->N[i]] = 4;
	}
	
}

long unsigned binRev(long unsigned mer) {
	
	/* swap consecutive pairs */
	mer = ((mer >> 2) & 0x3333333333333333) | ((mer & 0x3333333333333333) << 2);
	/* swap nibbles */
	mer = ((mer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((mer & 0x0F0F0F0F0F0F0F0F) << 4);
	/* swap bytes */
	mer = ((mer >> 8) & 0x00FF00FF00FF00FF) | ((mer & 0x00FF00FF00FF00FF) << 8);
	/* swap 2-bytes */
	mer = ((mer >> 16) & 0x0000FFFF0000FFFF) | ((mer & 0x0000FFFF0000FFFF) << 16);
	/* swap 4-bytes */
	return ((mer >> 32) | (mer << 32));
}

void rc_comp(struct compDNA *compressor, struct compDNA *compressor_rc) {
	
	int i, j, shift, r_shift;
	
	compressor_rc->seqlen = compressor->seqlen;
	compressor_rc->complen = compressor->complen;
	
	/* reverse and complement*/
	for(i = 0, j = compressor->complen - 1; i < compressor->complen; i++, j--) {
		compressor_rc->seq[j] = binRev(~compressor->seq[i]);
	}
	
	/* shift */
	if((compressor->seqlen & 31)) {
		shift = (((compressor->complen << 5) - compressor->seqlen) << 1);
		r_shift = 64 - shift;
		for(i = 0, j = 1; j < compressor->complen; i++, j++) {
			compressor_rc->seq[i] = (compressor_rc->seq[i] << shift) | (compressor_rc->seq[j] >> r_shift);
		}
		compressor_rc->seq[i] <<= shift;
	}
	
	/* add N's */
	shift = compressor->seqlen - 1;
	compressor_rc->N[0] = compressor->N[0];
	for(i = 1, j = compressor->N[0]; i <= compressor->N[0]; i++, j--) {
		compressor_rc->N[i] = shift - compressor->N[j];
	}
	
}

void dumpComp(struct compDNA *compressor, FILE* file) {
	
	fwrite(&compressor->seqlen, sizeof(int), 1, file);
	fwrite(&compressor->complen, sizeof(int), 1, file);
	fwrite(compressor->seq, sizeof(long unsigned), compressor->complen, file);
	fwrite(compressor->N, sizeof(int), compressor->N[0] + 1, file);
	
}

int loadComp(struct compDNA *compressor, FILE* file) {
	
	if(fread(&compressor->seqlen, sizeof(int), 1, file)) {
		fread(&compressor->complen, sizeof(int), 1, file);
		fread(compressor->seq, sizeof(long unsigned), compressor->complen, file);
		fread(compressor->N, sizeof(int), 1, file);
		fread(compressor->N + 1, sizeof(int), compressor->N[0], file);
		return 1;
	}
	return 0;
}

int getComp(struct compDNA *compressor, FILE* file) {
	
	compressor->seqlen = 0;
	fread(&compressor->seqlen, sizeof(int), 1, file);
	if(!compressor->seqlen) {
		return 0;
	}
	
	fread(&compressor->complen, sizeof(int), 1, file);
	
	/* realloc */
	if(compressor->seqlen >= compressor->size) {
		free(compressor->N);
		free(compressor->seq);
		if(compressor->seqlen & 31) {
			compressor->size = (compressor->seqlen >> 5) + 1;
			compressor->size <<= 6;
		} else {
			compressor->size = compressor->seqlen << 1;
		}
		
		compressor->seq = calloc(compressor->size >> 5, sizeof(long unsigned));
		compressor->N = malloc((compressor->size + 1) * sizeof(int));
		if(!compressor->seq || !compressor->N) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		compressor->N[0] = 0;
	}
	
	fread(compressor->seq, sizeof(long unsigned), compressor->complen, file);
	fread(compressor->N, sizeof(int), 1, file);
	fread(compressor->N + 1, sizeof(int), compressor->N[0], file);
	
	return 1;
}

void dumpCompKmers(struct compKmers *compressor, FILE *file) {
	
	fwrite(&compressor->n, sizeof(int), 1, file);
	fwrite(compressor->kmers, sizeof(long unsigned), compressor->n, file);
	
}

int getCompKmers(struct compKmers *compressor, FILE *file) {
	fread(&compressor->n, sizeof(int), 1, file);
	
	if(!compressor->n) {
		return 0;
	}
	
	if(compressor->n >= compressor->size) {
		free(compressor->kmers);
		allocCompKmers(compressor, compressor->n << 1);
	}
	
	fread(compressor->kmers, sizeof(long unsigned), compressor->n, file);
	
	return 1;
}

/* HASHMAP FUNCTIONS */

int * hashMap_getGlobal_old(long unsigned key) {
	
	unsigned pos;
	
	pos = templates->exist[key & templates->size];
	
	if(pos != templates->null_index) {
		for(; templates->key_index[pos] != templates->null_index; pos++) {
			if(key == getKmer(templates->seq, templates->key_index[pos])) {
				return templates->values + (templates->value_index[pos]);
			}
		}
	}
	
	return 0;
}

int * hashMap_getGlobal(long unsigned key) {
	
	unsigned pos, kpos;
	long unsigned kmer;
	
	kpos = key & templates->size;
	pos = templates->exist[kpos];
	
	if(pos != templates->null_index) {
		kmer = getKmer(templates->seq, templates->key_index[pos]);
		while(key != kmer) {
			//if(kpos != (kmer & templates->size) || pos == templates->n) {
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			pos++;
			kmer = getKmer(templates->seq, templates->key_index[pos]);
		}
		return (templates->values + templates->value_index[pos]);
	}
	
	return 0;
}

struct diskOffsets * getOffstets(char *filename) {
	
	long unsigned mask;
	struct diskOffsets *src;
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (templates->kmersize << 1));
	
	src = malloc(sizeof(struct diskOffsets));
	if(!src) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	if((templates->size - 1) == mask) {
		src->exist = 7 * sizeof(unsigned) + 2 * sizeof(long unsigned);
		src->seq = src->exist;
		src->values = src->seq + templates->size * sizeof(unsigned);
		src->key_index = src->values;
		src->value_index = src->key_index;
	} else {
		src->exist = 7 * sizeof(unsigned) + 2 * sizeof(long unsigned);
		src->seq = src->exist + templates->size * sizeof(unsigned);
		src->values = src->seq + templates->seqsize * sizeof(long unsigned);
		src->key_index = src->values + templates->v_index * sizeof(unsigned);
		src->value_index = src->key_index + (templates->n + 1) * sizeof(unsigned);
	}
	src->file = fopen(filename, "rb");
	if(!src->file) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist. 2\n");
		exit(-1);
	}
	
	return src;
}

unsigned getUnsignedFromFile(long pos, FILE *file) {
	
	unsigned dest;
	
	fseek(file, pos, SEEK_SET);
	fread(&dest, sizeof(unsigned), 1, file);
	
	return dest;
}

long unsigned getKmerFromFile(long offset, unsigned pos, FILE *file) {
	
	int cPos, iPos;
	long unsigned kmer[2];
	
	cPos = pos >> 5;
	iPos = (pos & 31) << 1;
	fseek(file, offset + cPos * sizeof(long unsigned), SEEK_SET);
	
	if(iPos <= shifter) {
		fread(kmer, sizeof(long unsigned), 1, file);
		return (*kmer << iPos) >> shifter;
	} else {
		fread(kmer, sizeof(long unsigned), 2, file);
		return ((*kmer << iPos) | (kmer[1] >> (64 - iPos))) >> shifter;
	}
}

int * hashMap_getGlobal_disk(long unsigned key) {
	
	unsigned i, pos, kpos, key_index[4];
	long unsigned kmer;
	
	kpos = key & templates->size;
	pos = getUnsignedFromFile(templates_offsets->exist + kpos * sizeof(unsigned), templates_offsets->file);
	
	if(pos != templates->null_index) {
		/* get key indexes */
		fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
		fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
		i = 0;
		kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
		while(key != kmer) {
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			i++;
			pos++;
			if(i == 4) {
				fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
				fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
				i = 0;
			}
			kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
		}
		pos = getUnsignedFromFile(templates_offsets->value_index + pos * sizeof(unsigned), templates_offsets->file);
		fseek(templates_offsets->file, templates_offsets->values + pos * sizeof(unsigned), SEEK_SET);
		fread(valuesFile, sizeof(unsigned), 1, templates_offsets->file);
		fread(valuesFile + 1, sizeof(unsigned), *valuesFile, templates_offsets->file);
		
		return valuesFile;
	}
	
	return 0;
}

int * megaMap_getGlobal_disk(long unsigned key) {
	
	unsigned pos;
	
	pos = getUnsignedFromFile(templates_offsets->exist + (key & templates->size) * sizeof(unsigned), templates_offsets->file);
	
	if(pos != 1) {
		fseek(templates_offsets->file, templates_offsets->values + pos * sizeof(unsigned), SEEK_SET);
		fread(valuesFile, sizeof(unsigned), 1, templates_offsets->file);
		fread(valuesFile + 1, sizeof(unsigned), *valuesFile, templates_offsets->file);
		
		return valuesFile;
	}
	
	return 0;
}

int * hashMap_getGlobal_semDisk(long unsigned key) {
	
	unsigned i, pos, kpos, key_index[4];
	long unsigned kmer;
	
	kpos = key & templates->size;
	pos = getUnsignedFromFile(templates_offsets->exist + kpos * sizeof(unsigned), templates_offsets->file);
	
	if(pos != templates->null_index) {
		/* get key indices */
		fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
		fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
		i = 0;
		kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
		while(key != kmer) {
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			i++;
			pos++;
			if(i == 4) {
				fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
				fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
				i = 0;
			}
			kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
		}
		pos = getUnsignedFromFile(templates_offsets->value_index + pos * sizeof(unsigned), templates_offsets->file);
		return templates->values + pos;
	}
	
	return 0;
}

int * megaMap_getGlobal_semDisk(long unsigned key) {
	
	unsigned pos;
	
	pos = getUnsignedFromFile(templates_offsets->exist + (key & templates->size) * sizeof(unsigned), templates_offsets->file);
	
	if(pos != 1) {
		return templates->values + pos;
	}
	
	return 0;
}

void getPrefix(struct hashMapKMA *dest, FILE *file) {
	
	/* load sizes */
	fread(&contamination, sizeof(unsigned), 1, file);
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	kmersize = dest->kmersize;
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	dest->exist = 0;
	dest->seq = 0;
	dest->values = 0;
	dest->key_index = 0;
	dest->value_index = 0;
}

int * megaMap_getGlobal(long unsigned key) {
	
	if(templates->exist[key] != 1) {
		return templates->values + templates->exist[key];
	}
	return 0;
}

void hashMapKMA_load(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	key_t key;
	int shmid;
	long unsigned seekSize;
	
	/* load sizes */
	fread(&contamination, sizeof(unsigned), 1, file);
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	kmersize = dest->kmersize;
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	/* check shared memory, else load */
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), 0666);
	if(shmid < 0) {
		/* not shared, load */
		dest->exist = malloc(dest->size * sizeof(unsigned));
		if(!dest->exist) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		seekSize = 0;
		fread(dest->exist, sizeof(unsigned), dest->size, file);
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
		seekSize = dest->size * sizeof(unsigned);
	}
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	if((dest->size - 1) == mask) {
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(int), 0666);
		if(shmid < 0) {
			/* not shared, load */
			dest->values = malloc(dest->v_index * sizeof(int));
			if(!dest->values) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			fseek(file, seekSize, SEEK_CUR);
			seekSize = 0;
			fread(dest->values, sizeof(int), dest->v_index, file);
		} else {
			/* found */
			dest->values = shmat(shmid, NULL, 0);
			seekSize += dest->v_index * sizeof(int);
		}
		hashMap_get = &megaMap_getGlobal;
	} else {
		key = ftok(filename, 's');
		shmid = shmget(key, dest->seqsize * sizeof(long unsigned), 0666);
		if(shmid < 0) {
			/* not shared, load */
			dest->seq = malloc(dest->seqsize * sizeof(long unsigned));
			if(!dest->seq) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			fseek(file, seekSize, SEEK_CUR);
			seekSize = 0;
			fread(dest->seq, sizeof(long unsigned), dest->seqsize, file);
		} else {
			/* found */
			dest->seq = shmat(shmid, NULL, 0);
			seekSize += dest->seqsize * sizeof(long unsigned);
		}
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(int), 0666);
		if(shmid < 0) {
			/* not shared, load */
			dest->values = malloc(dest->v_index * sizeof(int));
			if(!dest->values) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			fseek(file, seekSize, SEEK_CUR);
			seekSize = 0;
			fread(dest->values, sizeof(int), dest->v_index, file);
		} else {
			/* found */
			dest->values = shmat(shmid, NULL, 0);
			seekSize += dest->v_index * sizeof(int);
		}
		key = ftok(filename, 'k');
		shmid = shmget(key, (dest->n + 1) * sizeof(unsigned), 0666);
		if(shmid < 0) {
			/* not shared, load */
			dest->key_index = malloc((dest->n + 1) * sizeof(unsigned));
			if(!dest->key_index) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			fseek(file, seekSize, SEEK_CUR);
			seekSize = 0;
			fread(dest->key_index, sizeof(unsigned), dest->n + 1, file);
		} else {
			/* found */
			dest->key_index = shmat(shmid, NULL, 0);
			seekSize += (dest->n + 1) * sizeof(unsigned);
		}
		key = ftok(filename, 'i');
		shmid = shmget(key, dest->n * sizeof(unsigned), 0666);
		if(shmid < 0) {
			/* not shared, load */
			dest->value_index = malloc(dest->n * sizeof(unsigned));
			if(!dest->value_index) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			fseek(file, seekSize, SEEK_CUR);
			seekSize = 0;
			fread(dest->value_index, sizeof(unsigned), dest->n, file);
		} else {
			/* found */
			dest->value_index = shmat(shmid, NULL, 0);
		}
	}
}

void hashMapKMA_load_shm(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	key_t key;
	int shmid;
	
	/* load sizes */
	fread(&contamination, sizeof(unsigned), 1, file);
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	kmersize = dest->kmersize;
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	/* check shared memory */
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), 0666);
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB not shared, see kma_shm\n");
		exit(-2);
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
	}
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	
	if((dest->size - 1) == mask) {
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(int), 0666);
		if(shmid < 0) {
			/* not shared */
			fprintf(stderr, "DB not shared, see kma_shm\n");
			exit(-2);
		} else {
			/* found */
			dest->values = shmat(shmid, NULL, 0);
		}
	} else {
		key = ftok(filename, 's');
		shmid = shmget(key, dest->seqsize * sizeof(long unsigned), 0666);
		if(shmid < 0) {
			/* not shared, load */
			fprintf(stderr, "DB not shared, see kma_shm\n");
			exit(-2);
		} else {
			/* found */
			dest->seq = shmat(shmid, NULL, 0);
		}
		
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(int), 0666);
		if(shmid < 0) {
			/* not shared */
			fprintf(stderr, "DB not shared, see kma_shm\n");
			exit(-2);
		} else {
			/* found */
			dest->values = shmat(shmid, NULL, 0);
		}
		
		key = ftok(filename, 'k');
		shmid = shmget(key, (dest->n + 1) * sizeof(unsigned), 0666);
		if(shmid < 0) {
			/* not shared */
			fprintf(stderr, "DB not shared, see kma_shm\n");
			exit(-2);
		} else {
			/* found */
			dest->key_index = shmat(shmid, NULL, 0);
		}
		
		key = ftok(filename, 'i');
		shmid = shmget(key, dest->n * sizeof(unsigned), 0666);
		if(shmid < 0) {
			/* not shared */
			fprintf(stderr, "DB not shared, see kma_shm\n");
			exit(-2);
		} else {
			/* found */
			dest->value_index = shmat(shmid, NULL, 0);
		}
	}
}

/* hashMap indexes */
void hashMap_index_initialize(struct hashMap_index *dest, int len) {
	
	dest->len = len;
	dest->size = len << 1;
	
	dest->index = malloc(dest->size * sizeof(int));
	dest->seq = malloc(((len >> 5) + 1) * sizeof(long unsigned));
	if(!dest->index || !dest->seq) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
}

void hashMap_index_set(struct hashMap_index *dest) {
	
	memset(dest->index, -1, dest->size * sizeof(int));
	memset(dest->seq, 0, ((dest->len >> 5) + 1) * sizeof(long unsigned));
}

void hashMap_index_destroy(struct hashMap_index *dest) {
	
	free(dest->index);
	free(dest->seq);
	free(dest);
}

int hashMap_index_get(struct hashMap_index *dest, long unsigned key) {
	
	unsigned index;
	int pos;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; index++) {
		if(getKmer(dest->seq, abs(pos) - 1) == key) {
			return pos;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; index++) {
			if(getKmer(dest->seq, abs(pos) - 1) == key) {
				return pos;
			}
		}
	}
	
	return 0;
}

int hashMap_index_getDub(struct hashMap_index *dest, long unsigned key, const char *qseq, int q_len, int *next) {
	
	int i, index, pos, max, score, mPos;
	max = 0;
	mPos = 0;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; index++) {
		if(pos < 0 && getKmer(dest->seq, (-1) - pos) == key) {
			pos = kmersize - pos + 1;
			score = 0;
			for(i = kmersize; i < q_len && pos < dest->len && getNuc(dest->seq, pos) == qseq[i]; i++, pos++) {
				score++;
			}
			if(score > max) {
				max = score;
				mPos = -dest->index[index];
			} else if(score == max) {
				mPos = 0;
				*next = score;
			}
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; index++) {
			if(pos < 0 && getKmer(dest->seq, (-1) - pos) == key) {
				pos = kmersize - pos + 1;
				score = 0;
				for(i = kmersize; i < q_len && pos < dest->len && getNuc(dest->seq, pos) == qseq[i]; i++, pos++) {
					score++;
				}
				if(score > max) {
					max = score;
					mPos = -dest->index[index];
				} else if(score == max) {
					mPos = 0;
					*next = score;
				}
			}
		}
	}
	
	return mPos;
}

int hashMap_index_get_bound(struct hashMap_index *dest, long unsigned key, int min, int max) {
	
	unsigned index;
	int pos;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; index++) {
		if(min < abs(pos) && abs(pos) < max && getKmer(dest->seq, abs(pos) - 1) == key) {
			return pos;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; index++) {
			if(min < abs(pos) && abs(pos) < max && getKmer(dest->seq, abs(pos) - 1) == key) {
				return pos;
			}
		}
	}
	
	return 0;
}

int hashMap_index_getDub_bound(struct hashMap_index *dest, long unsigned key, const char *qseq, int q_len, int *next, int min, int max) {
	
	int i, index, pos, maxScore, score, mPos;
	max = 0;
	mPos = 0;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; index++) {
		if(pos < 0 && max < pos && pos < min && getKmer(dest->seq, (-1) - pos) == key) {
			pos = kmersize - pos + 1;
			score = 0;
			for(i = kmersize; i < q_len && pos < max && getNuc(dest->seq, pos) == qseq[i]; i++, pos++) {
				score++;
			}
			if(score > maxScore) {
				maxScore = score;
				mPos = -dest->index[index];
			} else if(score == maxScore) {
				mPos = 0;
				*next = score;
			}
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; index++) {
			if(pos < 0 && max < pos && pos < min && getKmer(dest->seq, (-1) - pos) == key) {
				pos = kmersize - pos + 1;
				score = 0;
				for(i = kmersize; i < q_len && pos < max && getNuc(dest->seq, pos) == qseq[i]; i++, pos++) {
					score++;
				}
				if(score > maxScore) {
					maxScore = score;
					mPos = -dest->index[index];
				} else if(score == maxScore) {
					mPos = 0;
					*next = score;
				}
			}
		}
	}
	
	return mPos;
}

void hashMap_index_add(struct hashMap_index *dest, long unsigned key, int newpos) {
	
	int index, pos, neg;
	neg = 1;
	newpos++;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; index++) {
		if(pos > 0) {
			if(getKmer(dest->seq, pos - 1) == key) {
				dest->index[index] = 1 - dest->index[index];
				neg = -1;
			}
		} else {
			if(getKmer(dest->seq, 1 - pos) == key) {
				neg = -1;
			}
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; index++) {
			if(pos > 0) {
				if(getKmer(dest->seq, pos - 1) == key) {
					dest->index[index] = 1 - dest->index[index];
					neg = -1;
				}
			} else {
				if(getKmer(dest->seq, 1 - pos) == key) {
					neg = -1;
				}
			}
		}
	}
	
	if(index < dest->size) {
		dest->index[index] = neg * newpos;
	}
	
}

struct hashMap_index *hashMap_index_load(FILE *seq, FILE *index, int len) {
	
	struct hashMap_index *src;
	
	src = malloc(sizeof(struct hashMap_index));
	if(!src) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	hashMap_index_initialize(src, len);
	
	fread(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
	fread(src->index, sizeof(int), src->size, index);
	
	return src;
}

/* hashMap kmers for Sparse mapping */
void hashMap_kmers_initialize(struct hashMap_kmers *dest, unsigned newSize) {
	/* set hashMap */
	dest->size = newSize;
	dest->n = 0;
	/* set hashTable */
	dest->table = calloc(newSize, sizeof(struct hashTable_index*));
	if(!dest->table) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
}

void hashMap_kmers_CountIndex(struct hashMap_kmers *dest, long unsigned key) {
	unsigned i, j, index;
	struct hashTable_kmers *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		dest->n++;
		dest->table[index] = malloc(sizeof(struct hashTable_kmers));
		if(!dest->table[index]) {
			fprintf(stderr, "OOM\n");
			exit(-1);
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
				node->next = malloc(sizeof(struct hashTable_kmers));
				if(!node->next) {
					fprintf(stderr, "OOM\n");
					exit(-1);
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

/* DB LOADING */
struct hashMap_index * alignLoad_fly(FILE *seq_in, FILE *index_in, int len, long unsigned seq_index, long unsigned index_index) {
	
	/* move file pointer */
	fseek(index_in, index_index, SEEK_SET);
	fseek(seq_in, seq_index, SEEK_SET);
	
	return hashMap_index_load(seq_in, index_in, len);
}

struct hashMap_index * alignLoad_fly_mem(FILE *seq_in, FILE *index_in, int len, long unsigned seq_index, long unsigned index_index) {
	
	return hashMap_index_load(seq_in, index_in, len);
}

struct hashMap_index * alignLoad_fly_shm(FILE *seq_in, FILE *index_in, int len, long unsigned seq_index, long unsigned index_index) {
	
	struct hashMap_index *dest;
	
	dest = malloc(sizeof(struct hashMap_index));
	if(!dest) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	dest->len = len;
	dest->size = len << 1;
	dest->seq = templates_index[0]->seq + (seq_index / sizeof(long unsigned));
	dest->index = templates_index[0]->index + ((index_index - sizeof(int)) / sizeof(int));
	
	return dest;
}

struct hashMap_index * alignLoad_shm_initial(char *templatefilename, int file_len, FILE *seq_in, FILE *index_in) {
	
	key_t key;
	int shmid;
	struct hashMap_index *dest;
	long unsigned size;
	
	dest = malloc(sizeof(struct hashMap_index));
	if(!dest) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	dest->len = 0;
	dest->size = 0;
	
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".index.b");
	key = ftok(templatefilename, 'i');
	fread(&kmersize, sizeof(int), 1, index_in);
	fseek(index_in, 0, SEEK_END);
	size = ftell(index_in) - sizeof(int);
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		fprintf(stderr, "SHM error i.\n");
		exit(-2);
	} else {
		dest->index = shmat(shmid, NULL, 0);
	}
	
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".seq.b");
	key = ftok(templatefilename, 's');
	fseek(seq_in, 0, SEEK_END);
	size = ftell(seq_in);
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		fprintf(stderr, "SHM error s.\n");
		exit(-2);
	} else {
		dest->seq = shmat(shmid, NULL, 0);
	}
	rewind(seq_in);
	templatefilename[file_len] = 0;
	
	return dest;
}
void alignClean(int template) {
	
	hashMap_index_destroy(templates_index[template]);
	templates_index[template] = 0;
}

void alignClean_shm(int template) {
	free(templates_index[template]);
	templates_index[template] = 0;
}

void load_DBs_KMA(char *templatefilename) {
	/* load DBs needed for KMA */
	int i, j, file_len, file_size, shmid;
	FILE *DB_file;
	key_t key;
	
	/* allocate DBs */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist. 1\n");
		exit(-1);
	}
	fread(&DB_size, sizeof(int), 1, DB_file);
	contamination = DB_size + 1;
	if(shm & 4) {
		key = ftok(templatefilename, 'l');
		shmid = shmget(key, DB_size * sizeof(int), 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared length\n");
			exit(-2);
		} else {
			template_lengths = shmat(shmid, NULL, 0);
		}
	} else {
		template_lengths = malloc(DB_size * sizeof(int));
		if(!template_lengths) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		/* load lengths */
		fread(template_lengths, sizeof(int), DB_size, DB_file);
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	
	/* allocate pointers */
	template_names = malloc(DB_size * sizeof(char*));
	templates_index = calloc(DB_size, sizeof(struct hashMap_index*));
	alignment_scores = calloc(DB_size, sizeof(unsigned));
	uniq_alignment_scores = calloc(DB_size, sizeof(unsigned));
	if(!templates_index || !template_names || !alignment_scores || !uniq_alignment_scores) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	/* load names */
	strcat(templatefilename, ".name");
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "File coruption: %s\n", templatefilename);
		exit(-1);
	}
	
	/* get size of file */
	fseek(DB_file, 0, SEEK_END);
	file_size = ftell(DB_file);
	
	/* load file */
	if(shm & 16) {
		key = ftok(templatefilename, 'n');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared name\n");
			exit(-2);
		} else {
			template_names[0] = shmat(shmid, NULL, 0);
		}
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; i++) {
			if(template_names[0][i] == 0) {
				template_names[j] = template_names[0] + i + 1;
				j++;
			}
		}
	} else {
		rewind(DB_file);
		template_names[0] = malloc(file_size);
		if(!template_names[0]) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		fread(template_names[0], 1, file_size, DB_file);
		template_names[0][file_size - 1] = 0;
		
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; i++) {
			if(template_names[0][i] == '\n') {
				template_names[0][i] = 0;
				template_names[j] = template_names[0] + i + 1;
				j++;
			}
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
}

void load_DBs_Sparse(char *templatefilename) {
	/* load DBs needed for KMA */
	int i, j, file_len, file_size, shmid;
	FILE *DB_file;
	key_t key;
	
	
	/* Open DB */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	templatefilename[file_len + 9] = 0;
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "File coruption: %s\n", templatefilename);
		exit(-1);
	}
	
	/* allocate DBs */
	fread(&DB_size, sizeof(int), 1, DB_file);
	contamination = DB_size + 1;
	if(shm & 4) {
		fseek(DB_file, 0, SEEK_END);
		file_size = ftell(DB_file) - sizeof(int);
		key = ftok(templatefilename, 'l');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared length\n");
			exit(-2);
		} else {
			template_lengths = shmat(shmid, NULL, 0);
		}
		template_ulengths = template_lengths + DB_size;
	} else {
		template_lengths = malloc(DB_size * sizeof(int));
		template_ulengths = malloc(DB_size * sizeof(int));
		if(!template_lengths || !template_ulengths) {
			fprintf(stderr, "Out of memory\n");
			exit(-1);
		}
		/* load lengths */
		fread(template_lengths, sizeof(int), DB_size, DB_file);
		fread(template_ulengths, sizeof(int), DB_size, DB_file);
		
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	/* load names */
	strcat(templatefilename, ".name");
	templatefilename[file_len + 5] = 0;
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "File coruption: %s\n", templatefilename);
		exit(-1);
	}
	/* get size of file */
	fseek(DB_file, 0, SEEK_END);
	file_size = ftell(DB_file);
	template_names = malloc(DB_size * sizeof(char*));
	if(!template_names) {
		fprintf(stderr, "Out of memory\n");
		exit(-1);
	}
	
	/* load file */
	if(shm & 16) {
		key = ftok(templatefilename, 'n');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared name\n");
			exit(-2);
		} else {
			template_names[0] = shmat(shmid, NULL, 0);
		}
		for(i = 0, j = 2; j < DB_size; i++) {
			if(template_names[0][i] == 0) {
				template_names[j] = template_names[0] + i + 1;
				j++;
			}
		}
	} else {
		rewind(DB_file);
		template_names[0] = malloc(file_size);
		if(!template_names[0]) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		fread(template_names[0], 1, file_size, DB_file);
		template_names[0][file_size - 1] = 0;
		
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; i++) {
			if(template_names[0][i] == '\n') {
				template_names[0][i] = 0;
				template_names[j] = template_names[0] + i + 1;
				j++;
			}
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
}

/* METHOD SPECIFIC METHODS */
void print_ankers(int *out_Tem, struct compDNA *qseq, int rc_flag, char *header, int header_len) {
	
	dumpComp(qseq, stdout);
	fwrite(&(int){rc_flag}, sizeof(int), 1, stdout);
	fwrite(out_Tem, sizeof(int), *out_Tem + 1, stdout);
	fwrite(&(int){header_len}, sizeof(int), 1, stdout);
	fwrite(header, 1, header_len, stdout);
	
}

void ankerAndClean(struct compDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, char *header, int header_len) {
	
	int k, l, bestHitsCov, template, *values;
	double thisCov, bestCov;
	struct compDNA tmpQseq;
	
	/* here */
	// make sure cuts isn't random seeds
	
	bestHitsCov = *regionTemplates;
	bestCov = 0;
	for(k = 1; k <= *regionTemplates; k++) {
		template = regionTemplates[k];
		if(template < 0) {
			template = -template;
			thisCov = 1.0 * Score_r[template] / template_lengths[template];
		} else {
			thisCov = 1.0 * Score[template] / template_lengths[template];
		}
		if(thisCov > bestCov) {
			bestCov = thisCov;
		}
	}
	
	for(k = start_cut; k <= end_cut; k++) {
		if(VF_scores[k]) {
			values = VF_scores[k];
			for(l = 1; l <= *values; l++) {
				if(Score[values[l]] != bestScore && values[l] != contamination) {
					thisCov = 1.0 * Score[values[l]] / template_lengths[values[l]];
					if(thisCov > bestCov) {
						bestCov = thisCov;
						bestHitsCov = *regionTemplates + 1;
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
		if(VR_scores[k]) {
			values = VR_scores[k];
			for(l = 1; l <= *values; l++) {
				if(Score_r[values[l]] != bestScore && values[l] != contamination) {
					thisCov = 1.0 * Score_r[values[l]] / template_lengths[values[l]];
					if(thisCov > bestCov) {
						HIT = -1;
						bestCov = thisCov;
						bestHitsCov = *regionTemplates + 1;
						regionTemplates[bestHitsCov] = -values[l];
					} else if(thisCov == bestCov) {
						HIT = -1;
						bestHitsCov++;
						regionTemplates[bestHitsCov] = -values[l];
					}
				}
				Score_r[values[l]]--;
			}
			VR_scores[k] = 0;
		}
	}
	*regionTemplates = bestHitsCov;
	
	/* clear nearest templates on both sides of match */
	for(k = ((start_cut - 92) < 0) ? 0 : (start_cut - 92); k < start_cut; k++) {
		if(VF_scores[k]) {
			values = VF_scores[k];
			for(l = 1; l <= *values; l++) {
				Score[values[l]]--;
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			values = VR_scores[k];
			for(l = 1; l <= *values; l++) {
				Score_r[values[l]]--;
			}
			VR_scores[k] = 0;
		}
	}
	for(k = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92); k > end_cut; k--) {
		if(VF_scores[k]) {
			values = VF_scores[k];
			for(l = 1; l <= *values; l++) {
				Score[values[l]]--;
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			values = VR_scores[k];
			for(l = 1; l <= *values; l++) {
				Score_r[values[l]]--;
			}
			VR_scores[k] = 0;
		}
	}
	
	
	
	/* modify limits of match seq */
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	start_cut = (start_cut >> 5) << 5;
	end_cut = ((end_cut >> 5) << 5) + 32;
	end_cut = (end_cut < qseq->seqlen) ? end_cut : qseq->seqlen;
	tmpQseq.seqlen = (end_cut - start_cut);
	tmpQseq.seq = qseq->seq + (start_cut >> 5);
	tmpQseq.N = tmpNs;
	
	for(k = 1, l = 0; k < qseq->N[0]; k++) {
		if(start_cut <= qseq->N[k]) {
			l++;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				l--;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	tmpQseq.seqlen--;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		tmpQseq.seqlen--;
		l--;
	}
	tmpQseq.seqlen++;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header, header_len);
}

void ankerAndClean_MEM(struct compDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, char *header, int header_len) {
	
	int k, l, *values;
	struct compDNA tmpQseq;
	/* clean up scores */
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	for(k = start_cut; k < end_cut; k++) {
		if(VF_scores[k]) {
			values = VF_scores[k];
			for(l = 1; l <= *values; l++) {
				Score[values[l]]--;
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			values = VR_scores[k];
			for(l = 1; l <= *values; l++) {
				Score_r[values[l]]--;
			}
			VR_scores[k] = 0;
		}
	}
	
	/* modify limits of match seq */
	start_cut = (start_cut >> 5) << 5;
	end_cut = ((end_cut >> 5) << 5) + 32;
	end_cut = (end_cut < qseq->seqlen) ? end_cut : qseq->seqlen;
	tmpQseq.seqlen = (end_cut - start_cut);
	tmpQseq.seq = qseq->seq + (start_cut >> 5);
	tmpQseq.N = tmpNs;
	
	for(k = 1, l = 0; k < qseq->N[0]; k++) {
		if(start_cut <= qseq->N[k]) {
			l++;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				l--;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	tmpQseq.seqlen--;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		tmpQseq.seqlen--;
		l--;
	}
	tmpQseq.seqlen++;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header, header_len);
}

void deConPrint(int *out_Tem, struct compDNA *qseq, int rc_flag, char *header, int header_len) {
	
	int contPos;
	
	if((contPos = find_contamination(out_Tem)) != -1) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		out_Tem[0]--;
	}
	if((contPos = find_contamination2(out_Tem, -contamination)) != -1) {
		out_Tem[contPos] = out_Tem[*out_Tem];
		out_Tem[0]--;
	}
	
	if(*out_Tem > 0) {
		printPtr(out_Tem, qseq, rc_flag, header, header_len);
	}
}

FILE * printFrags(char *filename, struct frag **alignFrags) {
	
	int i, seqlen;
	FILE *OUT;
	struct frag *alignFrag, *next;
	
	OUT = fopen(filename, "wb");
	if(!OUT) {
		fprintf(stderr, "File coruption: %s\n", filename);
		exit(-1);
	}
	for(i = 0; i < DB_size; i++) {
		if(alignFrags[i] != 0) {
			for(alignFrag = alignFrags[i]; alignFrag != 0; alignFrag = next) {
				next = alignFrag->next;
				
				fwrite(&(int){i}, sizeof(int), 1, OUT);
				fwrite(&(int){alignFrag->q_len}, sizeof(int), 1, OUT);
				fwrite(alignFrag->qseq, 1, alignFrag->q_len, OUT);
				fwrite(&(int){alignFrag->bestHits}, sizeof(int), 1, OUT);
				fwrite(&(int){alignFrag->score}, sizeof(int), 1, OUT);
				fwrite(&(int){alignFrag->start}, sizeof(int), 1, OUT);
				fwrite(&(int){alignFrag->end}, sizeof(int), 1, OUT);
				seqlen = strlen(alignFrag->header)+1;
				fwrite(&(int){seqlen}, sizeof(int), 1, OUT);
				fwrite(alignFrag->header, 1, seqlen, OUT);
				
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

void save_kmers(struct compDNA *qseq, struct compDNA *qseq_r, char *header, int header_len) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, k, l, end, HIT, bestHits, hitCounter, bestScore, bestScore_r;
	int contPos, maxScore, reps, *values, *last;
	long unsigned key;
	
	if(qseq->seqlen < kmersize) {
		return;
	}
	
	bestScore = 0;
	bestScore_r = 0;
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	
	/* Search forward strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	j = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	for(i = 1; i <= qseq->N[0] && !HIT; i++) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(getKmer(qseq->seq, j))) {
				HIT = 1;
			}
		}
		j = qseq->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates = 0;
		last = 0;
		reps = 0;
		j = 0;
		for(i = 1; i <= qseq->N[0]; i++) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; j++) {
				if((values = hashMap_get(getKmer(qseq->seq, j)))) {
					if(values == last) {
						reps++;
					} else {
						if(last) {
							for(l = 1; l <= *last; l++) {
								if(Score[last[l]] > 0) {
									Score[last[l]] += reps;
								} else {
									Score[last[l]] = reps;
									bestTemplates[0]++;
									bestTemplates[*bestTemplates] = last[l];
								}
							}
							hitCounter += reps;
						}
						reps = 1;
						last = values;
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		if(last) {
			for(l = 1; l <= *last; l++) {
				if(Score[last[l]] > 0) {
					Score[last[l]] += reps;
				} else {
					Score[last[l]] = reps;
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = last[l];
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; l++) {
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
			*bestTemplates = bestHits;
		} else {
			for(l = 1; l <= *bestTemplates; l++) {
				Score[bestTemplates[l]] = 0;
			}
			*bestTemplates = 0;
		}
	}
	qseq->N[0]--;
	
	/* search rc strand */
	/* Make quick check of the qseq */
	HIT = exhaustive;
	j = 0;
	qseq_r->N[0]++;
	qseq_r->N[qseq_r->N[0]] = qseq_r->seqlen;
	for(i = 1; i <= qseq_r->N[0] && !HIT; i++) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(getKmer(qseq_r->seq, j))) {
				HIT = 1;
			}
		}
		j = qseq_r->N[i] + 1;
	}
	
	/* If deltamer qseq hits, then continue */
	if(HIT) {
		/* Scan the deltamer exhaustively, and collect scores in Score*/
		hitCounter = 0;
		*bestTemplates_r = 0;
		last = 0;
		reps = 0;
		j = 0;
		for(i = 1; i <= qseq_r->N[0]; i++) {
			end = qseq_r->N[i] - kmersize + 1;
			for(;j < end; j++) {
				if((values = hashMap_get(getKmer(qseq_r->seq, j)))) {
					if(values == last) {
						reps++;
					} else {
						if(last) {
							for(l = 1; l <= *last; l++) {
								if(Score_r[last[l]] > 0) {
									Score_r[last[l]] += reps;
								} else {
									Score_r[last[l]] = reps;
									bestTemplates_r[0]++;
									bestTemplates_r[*bestTemplates_r] = last[l];
								}	
							}
							hitCounter += reps;
						}
						reps = 1;
						last = values;
					}
				}
			}
			j = qseq_r->N[i] + 1;
		}
		if(last) {
			for(l = 1; l <= *last; l++) {
				if(Score_r[last[l]] > 0) {
					Score_r[last[l]] += reps;
				} else {
					Score_r[last[l]] = reps;
					bestTemplates_r[0]++;
					bestTemplates_r[*bestTemplates_r] = last[l];
				}	
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; l++) {
				if(Score_r[bestTemplates_r[l]] > bestScore_r) {
					bestScore_r = Score_r[bestTemplates_r[l]];
					bestHits = 1;
					bestTemplates_r[bestHits] = bestTemplates_r[l];
				} else if(Score_r[bestTemplates_r[l]] == bestScore_r) {
					bestHits++;
					bestTemplates_r[bestHits] = bestTemplates_r[l];
				}
				Score_r[bestTemplates_r[l]] = 0;
			}
			*bestTemplates_r = bestHits;
		} else {
			for(l = 1; l <= *bestTemplates_r; l++) {
				Score_r[bestTemplates_r[l]] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	if(bestScore > 0 || bestScore_r > 0) {
		if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore + kmersize)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r + kmersize))) {
			/* merge results */
			if(bestScore == bestScore_r) {
				for(i = 1; i <= *bestTemplates_r; i++) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				deConPrintPtr(bestTemplates, qseq, -bestScore, header, header_len);
			} else if(bestScore_r > bestScore) {
				deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header, header_len);
			} else {
				deConPrintPtr(bestTemplates, qseq, bestScore, header, header_len);
			}
		}
	}
}

void save_kmers_HMM(struct compDNA *qseq, struct compDNA *qseq_r, char *header, int header_len) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, k, l, N, n, i_r, j_r, seqlen, seqend, end, HIT, Ncheck, template;
	int hitCounter, bestScore, bestHits, start, stop, start_cut, end_cut;
	int reps, rreps, *values, *last, *rlast;
	double Ms, Ns, Ms_prev, Ns_prev;
	
	if(qseq->seqlen < kmersize) {
		return;
	}
	
	/* reverse complement qseq */
	rc_comp(qseq, qseq_r);
	
	seqlen = qseq->seqlen;
	seqend = seqlen - kmersize + 1;
	i = 0;
	i_r = seqlen - kmersize;
	N = 1;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = seqlen;
	while(N <= qseq->N[0]) {
		/* find a seed */
		HIT = 0;
		end = qseq->N[N] - kmersize + 1;
		if(exhaustive) {
			while(i < end && !HIT) {
				if(hashMap_get(getKmer(qseq->seq, i)) || hashMap_get(getKmer(qseq_r->seq, i_r))) {
					HIT = 1;
				} else {
					i++;
					i_r--;
				}
			}
		} else {
			while(i < end && !HIT) {
				if(hashMap_get(getKmer(qseq->seq, i)) || hashMap_get(getKmer(qseq_r->seq, i_r))) {
					HIT = 1;
				} else {
					i += kmersize;
					i_r -= kmersize;
				}
			}
		}
		
		/* evaluate seed */
		if(HIT) {
			/* set scores attr */
			bestScore = 0;
			*bestTemplates = 0;
			hitCounter = 1;
			
			/* save seed */
			VF_scores[i] = hashMap_get(getKmer(qseq->seq, i));
			VR_scores[i] = hashMap_get(getKmer(qseq_r->seq, i_r));
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			
			/* extend backward */
			j = i - 1;
			j_r = i_r + 1;
			n = N - 1;
			Ncheck = (n > 0) ? -1 : qseq->N[n];
			while(j >= 0) {
				if(j == Ncheck) {
					
					k = j;
					while(k >= kmersize && k < (j - kmersize)) {
						/* update next N check */
						if(k == Ncheck) {
							j = Ncheck;
							n--;
							Ncheck = (n > 0) ? -1 : qseq->N[n];
						}
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
							break;
						}
						k--;
						Ns_prev = Ns;
						Ms_prev = Ms;
					}
					
					if(k >= kmersize && k < (j - kmersize)) {
						j = k - 1;
						break;
					} else {
						j = k;
						j_r = seqlen - kmersize - k;
					}
				} else {
					VF_scores[j] = hashMap_get(getKmer(qseq->seq, j));
					VR_scores[j] = hashMap_get(getKmer(qseq_r->seq, j_r));
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						hitCounter++;
						
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
				}
				j--;
				j_r++;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			start = j + 1;
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			
			/* extend forward */
			j = i + 1;
			j_r = i_r - 1;
			Ncheck = qseq->N[N] - kmersize + 1;
			while(j < seqend) {
				if(j == Ncheck) {
					
					k = j;
					while(k < seqend && k < (j + kmersize)) {
						/* update next N check */
						if(k == Ncheck) {
							j = Ncheck;
							N++;
							Ncheck = (N == qseq->N[0]) ? seqlen : qseq->N[N] - kmersize + 1;
						}
						
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
							break;
						}
						k++;
						Ns_prev = Ns;
						Ms_prev = Ms;
					}
					
					if(k < seqend && k < (j + kmersize)) {
						j = k;
						break;
					} else {
						j = k;
						j_r = seqlen - kmersize - k;
					}
				} else {
					VF_scores[j] = hashMap_get(getKmer(qseq->seq, j));
					VR_scores[j] = hashMap_get(getKmer(qseq_r->seq, j_r));
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						hitCounter++;
						
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
				}
				j++;
				j_r--;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			stop = j + kmersize - 1;
			
			/* evaluate hit */
			if(hitCounter > 0 && (hitCounter * kmersize > (stop - start - hitCounter + kmersize)) && ((stop - start) > minLen || start == 0 || stop == seqlen)) {
				if(deCon) {
					for(k = start; k < j; k++) {
						if(((values =VF_scores[k]) && values[*values] == contamination) || ((values =VR_scores[k]) && values[*values] == contamination)) {
							hitCounter--;
						}
					}
				}
				
				/* accept hit */
				if(hitCounter > 0) {
					/* gain total scores and mapping templates for this region */
					*bestTemplates = 0;
					*bestTemplates_r = 0;
					
					last = 0;
					reps = 0;
					rlast = 0;
					rreps = 0;
					for(k = start; k < j; k++) {
						/* forward */
						if(VF_scores[k]) {
							if(VF_scores[k] == last) {
								reps++;
							} else {
								if(last) {
									for(l = 1; l <= *last; l++) {
										if(Score[last[l]] > 0) {
											Score[last[l]] += reps;
										} else {
											Score[last[l]] = reps;
											bestTemplates[0]++;
											bestTemplates[*bestTemplates] = last[l];
										}
									}
								}
								reps = 1;
								last = VF_scores[k];
							}
						}
						
						/* rc */
						if(VR_scores[k]) {
							if(VR_scores[k] == rlast) {
								rreps++;
							} else {
								if(rlast) {
									for(l = 1; l <= *rlast; l++) {
										if(Score_r[rlast[l]] > 0) {
											Score_r[rlast[l]] += rreps;
										} else {
											Score_r[rlast[l]] = rreps;
											bestTemplates_r[0]++;
											bestTemplates_r[*bestTemplates_r] = rlast[l];
										}
									}
								}
								rreps = 1;
								rlast = VR_scores[k];
							}
						}
						
					}
					if(last) {
						for(l = 1; l <= *last; l++) {
							if(Score[last[l]] > 0) {
								Score[last[l]] += reps;
							} else {
								Score[last[l]] = reps;
								bestTemplates[0]++;
								bestTemplates[*bestTemplates] = last[l];
							}
						}
					}
					if(rlast) {
						for(l = 1; l <= *rlast; l++) {
							if(Score_r[rlast[l]] > 0) {
								Score_r[rlast[l]] += rreps;
							} else {
								Score_r[rlast[l]] = rreps;
								bestTemplates_r[0]++;
								bestTemplates_r[*bestTemplates_r] = rlast[l];
							}
						}
					}
					
					/* cut out template hits */
					while(HIT != 0) {
						/* get best score */
						bestScore = 0;
						bestHits = 0;
						
						/* forward */
						for(k = 1; k <= *bestTemplates; k++) {
							template = bestTemplates[k];
							if(Score[template] > bestScore) {
								bestScore = Score[template];
								bestHits = 1;
								regionTemplates[bestHits] = template;
							} else if(Score[template] == bestScore) {
								if(Score[template]) {
									bestHits++;
									regionTemplates[bestHits] = template;
								} else {
									bestTemplates[k] = bestTemplates[*bestTemplates];
									bestTemplates[0]--;
									k--;
								}
							}
						}
						
						/* rc */
						for(k = 1; k <= *bestTemplates_r; k++) {
							template = bestTemplates_r[k];
							if(Score_r[template] > bestScore) {
								bestScore = Score_r[template];
								bestHits = 1;
								regionTemplates[bestHits] = -template;
							} else if(Score_r[template] == bestScore) {
								if(Score_r[template]) {
									bestHits++;
									regionTemplates[bestHits] = -template;
								} else {
									bestTemplates_r[k] = bestTemplates_r[*bestTemplates_r];
									bestTemplates_r[0]--;
									k--;
								}
							}
						}
						
						*regionTemplates = bestHits;
						
						if(bestScore > 0) {
							/* find limits of match */
							start_cut = j;
							for(k = 1; k <= bestHits; k++) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = start; l < start_cut; l++) {
									if(VR_scores[l] && intpos_bin(VR_scores[l], template) != -1) {
										start_cut = l;
									}
									if(VF_scores[l] && intpos_bin(VF_scores[l], template) != -1) {
										start_cut = l;
									}
								}
							}
							end_cut = start_cut;
							for(k = 1; k <= bestHits; k++) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = j; l > end_cut; l--) {
									if(VR_scores[l] && intpos_bin(VR_scores[l], template) != -1) {
										end_cut = l;
									}
									if(VF_scores[l] && intpos_bin(VF_scores[l], template) != -1) {
										end_cut = l;
									}
								}
							}
							
							/* evaluate best hit */
							if(bestScore * kmersize > (end_cut - start_cut - bestScore + kmersize)) {
								/* check for hits on rc */
								HIT = (regionTemplates[*regionTemplates] > 0) ? 1 : -1;
								/* print */
								ankerPtr(qseq, HIT, bestScore, start_cut, end_cut, header, header_len);
							} else {
								/* clear scores */
								for(k = 1; k <= *bestTemplates; k++) {
									Score[bestTemplates[k]] = 0;
								}
								for(k = 1; k <= *bestTemplates_r; k++) {
									Score_r[bestTemplates_r[k]] = 0;
								}
								HIT = 0;
							}
						} else {
							/* clear scores */
							for(k = 1; k <= *bestTemplates; k++) {
								Score[bestTemplates[k]] = 0;
							}
							for(k = 1; k <= *bestTemplates_r; k++) {
								Score_r[bestTemplates_r[k]] = 0;
							}
							
							HIT = 0;
						}
					}
				}
			}
			
			/* clear scores */
			for(k = start; k < j; k++) {
				VF_scores[k] = 0;
				VR_scores[k] = 0;
			}
			
			i = stop + 1;
			i_r = seqlen - kmersize - i;
		} else {
			N++;
		}
	}
}

void save_kmers_sparse(struct hashMap_kmers *foundKmers, struct compKmers *compressor) {
	/* save_kmers find ankering k-mers the in query sequence,
	and is the time determining step */
	
	int i;
	
	for(i = 0; i < compressor->n; i++) {
		if(hashMap_get(compressor->kmers[i])) {
			hashMap_kmers_CountIndex(foundKmers, compressor->kmers[i]);
		}
	}
	
}

struct hashTable * collect_Kmers(int *Scores, int *Scores_tot, struct hashMap_kmers *foundKmers, struct Hit *hits) {
	unsigned i, j;
	int *value;
	struct hashTable_kmers *node, *node_next;
	struct hashTable *kmerNode, *kmerList;
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	for(i = 0; i < foundKmers->size; i++) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			value = hashMap_get(node->key);
			if(value) {
				hits->n++;
				hits->tot += node->value;
				
				kmerNode = malloc(sizeof(struct hashTable));
				if(!kmerNode) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				kmerNode->value = malloc((*value + 1) * sizeof(unsigned));
				if(!kmerNode->value) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				kmerNode->value[0] = *value;
				for(j = 1; j <= *value; j++) {
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

struct hashTable ** collect_Kmers_deCon(int *Scores, int *Scores_tot, struct hashMap_kmers *foundKmers, struct Hit *hits) {
	unsigned i, j;
	int *value;
	struct hashTable_kmers *node, *node_next;
	struct hashTable *kmerNode, *kmerList;
	struct hashTable *decon_node, *deconList;
	struct hashTable **Returner;
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	deconList = 0;
	decon_node = 0;
	
	Returner = malloc(sizeof(struct hashTable *) << 1);
	if(!Returner) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	for(i = 0; i < foundKmers->size; i++) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			value = hashMap_get(node->key);
			if(value) {
				/* check for contamination */
				hits->n++;
				hits->tot += node->value;
				
				if(value[*value] == contamination) {
					decon_node = malloc(sizeof(struct hashTable));
					if(!decon_node) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					decon_node->value = malloc((*value + 1) * sizeof(unsigned));
					if(!decon_node->value) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					decon_node->value[0] = *value;
					for(j = 1; j <= *value; j++) {
						decon_node->value[j] = value[j];
					}
					decon_node->key = node->value;
					
					decon_node->next = deconList;
					deconList = decon_node;
				} else {
					kmerNode = malloc(sizeof(struct hashTable));
					if(!kmerNode) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					kmerNode->value = malloc((*value + 1) * sizeof(unsigned));
					if(!kmerNode->value) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					kmerNode->value[0] = *value;
					for(j = 1; j <= *value; j++) {
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
		if(intpos_bin(node->value, template) != -1) {
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
				if(!cont_node->next) {
					fprintf(stderr, "OOM\n");
					exit(-1);
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

void bootFsa(struct qseqs *header, struct qseqs *qseq, struct compDNA *compressor) {
	
	int i, end;
	
	/* bootstrap in pieces of 1024 */
	end = qseq->len - 1024;
	for(i = 0; i < end; i += 512) {
		compDNA(compressor, qseq->seq + i, 1024);
		/* print */
		dumpComp(compressor, stdout);
		fwrite(&header->len, sizeof(int), 1, stdout);
		fwrite((header->seq + 1), 1, header->len, stdout);
		resetComp(compressor);
	}
	
	compDNA(compressor, qseq->seq + i, qseq->len - i);
	/* print */
	dumpComp(compressor, stdout);
	fwrite(&header->len, sizeof(int), 1, stdout);
	fwrite((header->seq + 1), 1, header->len, stdout);
	resetComp(compressor);
}

void printFsa(struct qseqs *header, struct qseqs *qseq, struct compDNA *compressor) {
	
	/* translate to 2bit */
	if(qseq->len >= compressor->size) {
		freeComp(compressor);
		allocComp(compressor, qseq->len);
	}
	compDNA(compressor, qseq->seq, qseq->len);
	
	/* print */
	dumpComp(compressor, stdout);
	fwrite(&header->len, sizeof(int), 1, stdout);
	fwrite((header->seq + 1), 1, header->len, stdout);
	resetComp(compressor);
}

void run_input(char **inputfiles, int fileCount, int minPhred, int fiveClip) {
	
	int i, j, k, FASTA, FASTQ, fileCounter, phredCut, start, end;
	char *filename, *cmd, *zipped, *seq;
	struct qseqs *header, *qseq, *qual;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	freopen(NULL, "wb", stdout);
	
	compressor = malloc(sizeof(struct compDNA));
	cmd = malloc(1);
	zipped = strdup(".gz");
	if(!cmd || !zipped || !compressor) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(delta);
	qual = setQseqs(1024);
	inputfile = setFileBuff(1048576);
	
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
			openFileBuff(inputfile, filename, "rb");
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
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1;
				while(end >= 0 && seq[end] < phredCut) {
					end--;
				}
				end++;
				while(start < end && seq[start] < phredCut) {
					start++;
				}
				qseq->len = end - start;
				/* print */
				if(qseq->len > kmersize) {
					/* dump seq */
					printFsa_ptr(header, qseq, compressor);
					
					/* translate to 2bit */
					/*if(qseq->len >= compressor->size) {
						freeComp(compressor);
						allocComp(compressor, qseq->len);
					}
					compDNA(compressor, qseq->seq + start, qseq->len);
					*/
					/* print */
					/*dumpComp(compressor, stdout);
					fwrite(&header->len, sizeof(int), 1, stdout);
					fwrite((header->seq + 1), 1, header->len, stdout);
					resetComp(compressor);*/
				}
			}
		} else if(FASTA) {
			while(FileBuffgetFsa(inputfile, header, qseq)) {
				if(qseq->len > kmersize) {
					/* translate to 2bit */
					/*if(qseq->len >= compressor->size) {
						freeComp(compressor);
						allocComp(compressor, qseq->len);
					}
					compDNA(compressor, qseq->seq, qseq->len);*/
					
					printFsa_ptr(header, qseq, compressor);
					
					/* print */
					/*dumpComp(compressor, stdout);
					fwrite(&header->len, sizeof(int), 1, stdout);
					fwrite((header->seq + 1), 1, header->len, stdout);
					resetComp(compressor);*/
				}
			}
		}
		
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			pcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	fwrite(&(int){0}, sizeof(int), 1, stdout);
	
	free(cmd);
	free(zipped);
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
}

void run_input_sparse(char **inputfiles, int fileCount, int minPhred, int fiveClip) {
	
	int i, j, k, FASTA, FASTQ, fileCounter, phredCut, start, end, kmerCount;
	char *filename, *cmd, *zipped, *seq;
	struct qseqs *qseq, *qual;
	struct FileBuff *inputfile;
	struct compKmers *Kmers;
	freopen(NULL, "wb", stdout);
	
	Kmers = malloc(sizeof(struct compKmers));
	cmd = malloc(1);
	zipped = strdup(".gz");
	if(!cmd || !zipped || !Kmers) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocCompKmers(Kmers, delta);
	qseq = setQseqs(delta);
	qual = setQseqs(1024);
	inputfile = setFileBuff(1048576);
	
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
			openFileBuff(inputfile, filename, "rb");
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
		Kmers->n = 0;
		if(FASTQ) {
			/* get phred scale */
			phredCut = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			while(FileBuffgetFqSeq(inputfile, qseq, qual)) {
				/* trim */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1;
				while(end >= 0 && seq[end] < phredCut) {
					end--;
				}
				end++;
				while(start < end && seq[start] < phredCut) {
					start++;
				}
				qseq->len = end - start;
				
				/* print */
				if(qseq->len > kmersize) {
					/* translate to kmers */
					Kmers->n = translateToKmersAndDump(Kmers->kmers, Kmers->n, Kmers->size, qseq->seq + start, qseq->len, templates->prefix, templates->prefix_len);
				}
			}
			if(Kmers->n) {
				fwrite(Kmers->kmers, sizeof(long unsigned), Kmers->n, stdout);
				Kmers->n = 0;
			}
		} else if(FASTA) {
			while(FileBuffgetFsaSeq(inputfile, qseq)) {
				if(qseq->len > kmersize) {
					/* translate to kmers */
					Kmers->n = translateToKmersAndDump(Kmers->kmers, Kmers->n, Kmers->size, qseq->seq + start, qseq->len, templates->prefix, templates->prefix_len);
				}
			}
			if(Kmers->n) {
				fwrite(Kmers->kmers, sizeof(long unsigned), Kmers->n, stdout);
				Kmers->n = 0;
			}
		}
		
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			pcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	fwrite(&(int){0}, sizeof(int), 1, stdout);
	
	free(cmd);
	free(zipped);
	free(Kmers->kmers);
	free(Kmers);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
}

void save_kmers_batch(char *templatefilename, char *exePrev, int one2one) {
	
	int i, file_len, header_size, shmid;
	char *header;
	struct compDNA *qseq, *qseq_r;
	FILE *inputfile, *templatefile;
	time_t t0, t1;
	key_t key;
	
	/* open pipe */
	inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "File coruption\n");
		exit(-1);
	}
	freopen(NULL, "wb", stdout);
	t0 = clock();
	
	/* initialize seqs */
	qseq = malloc(sizeof(struct compDNA));
	qseq_r = malloc(sizeof(struct compDNA));
	if(!qseq || !qseq_r) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocComp(qseq, delta);
	allocComp(qseq_r, delta);
	header_size = 0;
	header = 0;
	
	/* load hashMap */
	file_len = strlen(templatefilename);
	if(deCon) {
		strcat(templatefilename, ".decon.comp.b");
	} else {
		strcat(templatefilename, ".comp.b");
	}
	templatefile = fopen(templatefilename, "rb" );
	if(!templatefile) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist. 2\n");
		exit(-1);
	}
	
	templates = malloc(sizeof(struct hashMapKMA));
	if(!templates) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	if(diskDB) {
		getPrefix(templates, templatefile);
		fclose(templatefile);
		templates_offsets = getOffstets(templatefilename);
		if(!one2one || diskDB & 2) {
			if(templates_offsets->key_index == templates_offsets->value_index) {
				hashMap_get = &megaMap_getGlobal_semDisk;
			} else {
				hashMap_get = &hashMap_getGlobal_semDisk;
			}
			key = ftok(templatefilename, 'v');
			shmid = shmget(key, templates->v_index * sizeof(unsigned), 0666);
			if(shmid < 0) {
				templates->values = malloc(templates->v_index * sizeof(unsigned));
				if(!templates->values) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				fseek(templates_offsets->file, templates_offsets->values, SEEK_SET);
				fread(templates->values, sizeof(unsigned), templates->v_index, templates_offsets->file);
			} else {
				templates->values = shmat(shmid, NULL, 0);
			}
		} else {
			if(templates_offsets->key_index == templates_offsets->value_index) {
				hashMap_get = &megaMap_getGlobal_disk;
			} else {
				hashMap_get = &hashMap_getGlobal_disk;
			}
			valuesFile = malloc((contamination + 1) * sizeof(unsigned));
			if(!valuesFile) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
	} else {
		hashMap_get = &hashMap_getGlobal;
		if((shm & 1) || (deCon && (shm & 2))) {
			hashMapKMA_load_shm(templates, templatefile, templatefilename);
		} else {
			hashMapKMA_load(templates, templatefile, templatefilename);
		}
		fclose(templatefile);
	}
	templatefilename[file_len] = 0;
	/* make indexing a masking problem */
	templates->size--;
	
	kmersize = templates->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	/* calculate HMM parameters */
	HMM_param[0] = log(1 - pow(0.25, kmersize));
	HMM_param[1] = log(pow(0.25, kmersize));
	HMM_param[2] = log(1 - pow(0.25, kmersize - 1) * 0.75);
	HMM_param[3] = log(pow(0.25, kmersize - 1) * 0.75);
	HMM_param[4] = log(1 - 1.0 / kmersize * 0.75 * 0.25);
	HMM_param[5] = log(1.0 / kmersize * 0.75 * 0.25);
	HMM_param[6] = log(0.75);
	HMM_param[7] = log(0.25);
	
	/* allocate scoring arrays */
	DB_size = contamination + 1;
	Score = calloc(DB_size, sizeof(int));
	Score_r = calloc(DB_size, sizeof(int));
	if(one2one) {
		bestTemplates = calloc((DB_size << 1) + 1, sizeof(int));
		bestTemplates_r = calloc(DB_size + 1, sizeof(int));
	} else {
		bestTemplates = calloc(DB_size + 1, sizeof(int));
		bestTemplates_r = calloc(DB_size + 1, sizeof(int));
		regionTemplates = calloc((DB_size << 1) + 1, sizeof(int));
		VF_scores = calloc(delta, sizeof(int*));
		VR_scores = calloc(delta, sizeof(int*));
		tmpNs = malloc(delta * sizeof(int));
		if(!regionTemplates || !VF_scores || !VR_scores || !tmpNs) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		/* load lengths */
		strcat(templatefilename, ".length.b");
		templatefile = fopen(templatefilename, "rb");
		if(!templatefile) {
			fprintf(stderr, "File coruption: %s\n", templatefilename);
			exit(-1);
		}
		fread(&DB_size, sizeof(int), 1, templatefile);
		if(shm & 4) {
			key = ftok(templatefilename, 'l');
			shmid = shmget(key, DB_size * sizeof(int), 0666);
			if(shmid < 0) {
				fprintf(stderr, "No shared length\n");
				exit(-2);
			} else {
				template_lengths = shmat(shmid, NULL, 0);
			}
		} else {
			template_lengths = malloc(DB_size * sizeof(int));
			if(!template_lengths) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			fread(template_lengths, sizeof(int), DB_size, templatefile);
		}
		templatefilename[file_len] = 0;
		fclose(templatefile);
		minLen = template_lengths[1] / 2;
		for(i = 1; i < DB_size; i++) {
			if(minLen > (template_lengths[i] / 2)) {
				minLen = template_lengths[i] / 2;
			}
		}
		minLen = (minLen < kmersize) ? kmersize : minLen;
	}
	if(!bestTemplates || !bestTemplates_r || !Score || !Score_r) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mer ankers\n");
	
	while(getComp(qseq, inputfile)) {
		if(qseq->seqlen >= delta) {
			delta = qseq->seqlen << 1;
			freeComp(qseq_r);
			allocComp(qseq_r, qseq->size);
			if(!one2one) {
				free(tmpNs);
				free(VF_scores);
				free(VR_scores);
				VF_scores = calloc(delta, sizeof(int*));
				VR_scores = calloc(delta, sizeof(int*));
				tmpNs = malloc(delta * sizeof(int));
				if(!VF_scores || !VR_scores || !tmpNs) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
			}
		}
		
		fread(&i, sizeof(int), 1, inputfile);
		if(i > header_size) {
			header_size = i << 1;
			free(header);
			header = malloc(header_size);
			if(!header) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		fread(header, i, 1, inputfile);
		
		kmerScan(qseq, qseq_r, header, i);
	}
	fwrite(&(int){0}, sizeof(int), 1, stdout);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used ankering query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	pclose(inputfile);
}

void save_kmers_sparse_batch(char *templatefilename, char *outputfilename, char *exePrev, char ss) {
	
	int i, file_len, prefix_len, stop, lenTresh, template, score, tmp_score;
	int score_add, score_tot_add, shmid;
	int *Scores, *Scores_tot, *w_Scores, *w_Scores_tot, *SearchList;
	unsigned Ntot, element, belong, *score_templates;
	double expected, q_value, p_value, etta, depth, cover, query_cover;
	double tot_depth, tot_cover, tot_query_cover, tmp_depth, tmp_cover;
	double tmp_expected, tmp_q, tmp_p;
	FILE *inputfile, *templatefile, *sparse_out;
	time_t t0, t1;
	key_t key;
	struct hashMap_kmers *foundKmers;
	struct hashTable *kmerList, *deConTable, *node, *prev, **Collecter;
	struct Hit Nhits, w_Nhits;
	struct compKmers *Kmers;
	
	/* here */
	// split input up
	
	/* open pipe */
	inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
	/* Load hashMap */
	t0 = clock();
	file_len = strlen(templatefilename);
	if(deCon) {
		strcat(templatefilename, ".decon.comp.b");
	} else {
		strcat(templatefilename, ".comp.b");
	}
	templatefile = fopen(templatefilename, "rb" );
	if(!templatefile) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist. 3\n");
		exit(-1);
	}
	templates = malloc(sizeof(struct hashMapKMA));
	if(!templates) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	if(diskDB) {
		getPrefix(templates, templatefile);
		fclose(templatefile);
		templates_offsets = getOffstets(templatefilename);
		if(diskDB & 2) {
			hashMap_get = &hashMap_getGlobal_semDisk;
			key = ftok(templatefilename, 'v');
			shmid = shmget(key, templates->v_index * sizeof(unsigned), 0666);
			if(shmid < 0) {
				templates->values = malloc(templates->v_index * sizeof(unsigned));
				if(!templates->values) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				fseek(templates_offsets->file, templates_offsets->values, SEEK_SET);
				fread(templates->values, sizeof(unsigned), templates->v_index, templates_offsets->file);
			} else {
				templates->values = shmat(shmid, NULL, 0);
			}
		} else {
			hashMap_get = &hashMap_getGlobal_disk;
			valuesFile = malloc((contamination + 1) * sizeof(unsigned));
			if(!valuesFile) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
	} else {
		hashMap_get = &hashMap_getGlobal;
		if((shm & 1) || (deCon && (shm & 2))) {
			hashMapKMA_load_shm(templates, templatefile, templatefilename);
		} else {
			hashMapKMA_load(templates, templatefile, templatefilename);
		}
		fclose(templatefile);
	}
	templatefilename[file_len] = 0;
	DB_size = contamination + 1;
	/* make indexing a masking problem */
	templates->size--;
	
	/* set dependent variables */
	kmersize = templates->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	/* load template attributes */
	load_DBs_Sparse(templatefilename);
	
	/* open output file */
	if(strcmp(outputfilename, "--") == 0) {
		sparse_out = stdout;
	} else {
		file_len = strlen(outputfilename);
		strcat(outputfilename, ".spa");
		sparse_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	}
	if(!sparse_out) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
	/* set hashMap for found kmers */
	foundKmers = malloc(sizeof(struct hashMap_kmers));
	if(!foundKmers) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	foundKmers->size = templates->n;
	hashMap_kmers_initialize(foundKmers, foundKmers->size);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mers\n");
	
	Kmers = malloc(sizeof(struct compKmers));
	if(!Kmers) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocCompKmers(Kmers, delta);
	Ntot = 0;
	
	/* count kmers */
	while((Kmers->n = fread(Kmers->kmers, sizeof(long unsigned), Kmers->size, inputfile))) {
		Ntot += Kmers->n;
		save_kmers_sparse(foundKmers, Kmers);
	}
	pclose(inputfile);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used to identify k-mers in query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding best matches and output results.\n");
	
	Scores = calloc(DB_size, sizeof(int));
	Scores_tot = calloc(DB_size, sizeof(int));
	SearchList = malloc(DB_size * sizeof(int));
	if(!Scores || !Scores_tot || !SearchList) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	etta = 1.0e-6;
	
	fprintf(sparse_out, "#Template\tNum\tScore\tExpected\tTemplate length\tquery_coverage\tCoverage\tDepth\ttot_query_coverage\ttot_coverage\ttot_depth\tq_value\tp_value\n");
	stop = 0;
	if(deCon) {
		/* start by removing contamination and collect scores */
		Collecter = collect_Kmers_deCon(Scores, Scores_tot, foundKmers, &Nhits);
		kmerList = Collecter[0];
		deConTable = Collecter[1];
		free(Collecter);
		
		fprintf(stderr, "# total number of matches: %d of %d kmers\n", Nhits.tot, Ntot);
		/* copy scores */
		w_Scores = malloc(DB_size * sizeof(int));
		w_Scores_tot = malloc(DB_size * sizeof(int));
		if(!w_Scores || !w_Scores_tot) {
			fprintf(stderr, "OOM\n");
			exit(-1);
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
						if(intpos_bin(node->value, template) != -1) {
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
			exit(-1);
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
		
		while(!stop) {
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
	fprintf(stderr, "# Total for finding and outputting best matches: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
}

struct alnScore NW(const long unsigned *template, const char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, struct aln *aligned) {
	
	int m, n, thisScore, t_len, q_len, pos[2];
	char *query, t_nuc;
	struct alnScore Stat;
	
	t_len = t_e - t_s;
	q_len = q_e - q_s;
	query = (char*)(queryOrg + q_s);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.gaps = 0;
			Stat.score = 0;
			aligned->s[0] = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.gaps = 0;
			Stat.score = W1 + (q_len - 1) * U;
			memset(aligned->s, '_', q_len);
			aligned->s[q_len] = 0;
			memset(aligned->t, 5, q_len);
			memcpy(aligned->q, query, q_len);
		} else {
			Stat.len = t_len;
			Stat.gaps = 0;
			Stat.score = W1 + (t_len - 1) * U;
			memset(aligned->s, '_', t_len);
			aligned->s[t_len] = 0;
			memset(aligned->q, 5, t_len);
			for(m = 0; m < t_len; m++) {
				aligned->t[m] = getNuc(template, (m + t_s));
			}
		}
		return Stat;
	}
	
	/* check matrix size */
	if(t_len >= NW_s || q_len >= NW_s) {
		for(m = 0; m < NW_s; m++) {
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
			NW_s = t_len << 1;
		} else {
			NW_s = q_len << 1;
		}
		D = malloc(NW_s * sizeof(int*));
		Q = malloc(NW_s * sizeof(int*));
		P = malloc(NW_s * sizeof(int*));
		E = malloc(NW_s * sizeof(int*));
		if(!D || !Q || !P || !E) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		for(m = 0; m < NW_s; m++) {
			D[m] = malloc(NW_s * sizeof(int));
			Q[m] = malloc(NW_s * sizeof(int));
			P[m] = malloc(NW_s * sizeof(int));
			E[m] = malloc(NW_s * sizeof(int));
			if(!D[m] || !Q[m] || !P[m] || !E[m]) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
	}
	
	/* fill in start penalties */
	if(k == 1) { // end alignment
		for(m = t_len; m >= 0; m--) {
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
		D[t_len][q_len] = 0;
		P[t_len][q_len] = 0;
		E[t_len][q_len] = 0;
	} else if (k == 2) { // extreme end alignment
		for(m = t_len; m >= 0; m--) {
			D[m][q_len] = 0;
			Q[m][q_len] = 0;
			E[m][q_len] = 0;
		}
		for(n = q_len; n >= 0; n--) {
			D[t_len][n] = 0;
			P[t_len][n] = 0;
			E[t_len][n] = 0;
		}
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
		Q[t_len][q_len] = 0;
		P[t_len][q_len] = 0;
		E[t_len][q_len] = 0;
	}
	
	
	for(m = t_len - 1; m >= 0; m--) {
		t_nuc = getNuc(template, (m + t_s));
		for(n = q_len - 1; n >= 0; n--) {
			/* Update D */
			D[m][n] = D[m + 1][n + 1] + d[t_nuc][query[n]];
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
	Stat.score = D[0][0];
	pos[0] = 0;
	pos[1] = 0;
	if(k == -1) { // front alignment
		for(m = 1; m < t_len; m++) {
			if(D[m][0] > Stat.score) {
				Stat.score = D[m][0];
				pos[0] = m;
			}
		}
	} else if(k == -2) { // extreme front alignment
		for(m = 1; m < t_len; m++) {
			if(D[m][0] > Stat.score) {
				Stat.score = D[m][0];
				pos[0] = m;
			}
		}
		for(n = 1; n < q_len; n++) {
			if(D[0][n] > Stat.score) {
				Stat.score = D[0][n];
				pos[1] = n;
			}
		}
		if(pos[1]) {
			pos[0] = 0;
		}
	}
	
	/* make back tracking */
	m = pos[0];
	n = pos[1];
	Stat.len = 0;
	Stat.gaps = 0;
	while(E[m][n] != 0) {
		if(E[m][n] == 1) {
			aligned->t[Stat.len] = getNuc(template, (m + t_s));
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = (aligned->t[Stat.len] == aligned->q[Stat.len]) ? '|' : '_';
			m++;
			n++;
		} else if(E[m][n] >= 4) {
			aligned->t[Stat.len] = getNuc(template, (m + t_s));
			aligned->q[Stat.len] = 5;
			aligned->s[Stat.len] = '_';
			m++;
		} else {
			aligned->t[Stat.len] = 5;
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = '_';
			Stat.gaps++;
			n++;
		}
		Stat.len++;
	}
	aligned->s[Stat.len] = 0;
	
	return Stat;
}

struct alnScore NW_score(const long unsigned *template, const char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e) {
	
	int m, n, thisScore, t_len, q_len, pos[2];
	char *query, t_nuc;
	struct alnScore Stat;
	
	t_len = t_e - t_s;
	q_len = q_e - q_s;
	query = (char*)(queryOrg + q_s);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.gaps = 0;
			Stat.score = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.gaps = 0;
			Stat.score = W1 + (q_len - 1) * U;
		} else {
			Stat.len = t_len;
			Stat.gaps = t_len;
			Stat.score = W1 + (t_len - 1) * U;
		}
		return Stat;
	}
	
	/* check matrix size */
	if(t_len >= NW_s || q_len >= NW_s) {
		for(m = 0; m < NW_s; m++) {
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
			NW_s = t_len << 1;
		} else {
			NW_s = q_len << 1;
		}
		D = malloc(NW_s * sizeof(int*));
		Q = malloc(NW_s * sizeof(int*));
		P = malloc(NW_s * sizeof(int*));
		E = malloc(NW_s * sizeof(int*));
		if(!D || !Q || !P || !E) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		for(m = 0; m < NW_s; m++) {
			D[m] = malloc(NW_s * sizeof(int));
			Q[m] = malloc(NW_s * sizeof(int));
			P[m] = malloc(NW_s * sizeof(int));
			E[m] = malloc(NW_s * sizeof(int));
			if(!D[m] || !Q[m] || !P[m] || !E[m]) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
	}
	
	/* fill in start penalties */
	if(k == 1) { // end alignment
		for(m = t_len; m >= 0; m--) {
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
		D[t_len][q_len] = 0;
		P[t_len][q_len] = 0;
		E[t_len][q_len] = 0;
	} else if(k == 2) { // extreme end alignment
		for(m = t_len; m >= 0; m--) {
			D[m][q_len] = 0;
			Q[m][q_len] = 0;
			E[m][q_len] = 0;
		}
		for(n = q_len; n >= 0; n--) {
			D[t_len][n] = 0;
			P[t_len][n] = 0;
			E[t_len][n] = 0;
		}
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
		Q[t_len][q_len] = 0;
		P[t_len][q_len] = 0;
		E[t_len][q_len] = 0;
	}
	
	
	for(m = t_len - 1; m >= 0; m--) {
		t_nuc = getNuc(template, (m + t_s));
		for(n = q_len - 1; n >= 0; n--) {
			/* Update D */
			D[m][n] = D[m + 1][n + 1] + d[t_nuc][query[n]];
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
	Stat.score = D[0][0];
	pos[0] = 0;
	pos[1] = 0;
	if(k == -1) { // front alignment
		for(m = 1; m < t_len; m++) {
			if(D[m][0] > Stat.score) {
				Stat.score = D[m][0];
				pos[0] = m;
			}
		}
	} else if(k == -2) { // extreme front alignment
		for(m = 1; m < t_len; m++) {
			if(D[m][0] > Stat.score) {
				Stat.score = D[m][0];
				pos[0] = m;
			}
		}
		for(n = 1; n < q_len; n++) {
			if(D[0][n] > Stat.score) {
				Stat.score = D[0][n];
				pos[1] = n;
			}
		}
		if(pos[1]) {
			pos[0] = 0;
		}
	}
	
	/* make back tracking */
	m = pos[0];
	n = pos[1];
	Stat.len = 0;
	Stat.gaps = 0;
	while(E[m][n] != 0) {
		if(E[m][n] == 1) {
			m++;
			n++;
		} else if(E[m][n] >= 4) {
			m++;
		} else {
			Stat.gaps++;
			n++;
		}
		Stat.len++;
	}
	
	return Stat;
}

struct alnScore KMA(const int template_name, const char *qseq, int q_len, struct aln *aligned, struct aln *Frag_align, int qBias, int min, int max) {
	
	int i, j, bias, prev, prev_index, stop, t_len, value, end, t_s, t_e, q_s, q_e;
	char nuc;
	long unsigned key;
	struct alnScore Stat, NWstat, recStat;
	struct aln *reAligned;
	struct hashMap_index *template_index;
	
	/* Extract indexes and template sequence */
	template_index = templates_index[template_name];
	t_len = template_lengths[template_name];
	
	/* Align the query and the template */
	recStat.score = 0;
	reAligned = 0;
	prev = 0;
	prev_index = 0;
	Stat.len = 0;
	Stat.score = 0;
	Stat.pos = t_len;
	Stat.gaps = 0;
	i = qBias;
	while(i < q_len) {
		end = charpos(qseq, 4, i, q_len);
		if(end == -1) {
			end = q_len;
		}
		
		if(i < end - kmersize) {
			key = makeKmer(qseq, i, kmersize - 1);
			i += (kmersize - 1);
		} else {
			i = end + 1;
		}
		
		while(i < end && prev_index < t_len) {
			key = ((key << 2) | qseq[i]) & mask;
			/* capture dublet */
			if((value = hashMap_index_get_bound(template_index, key, min, max)) < 0) {
				value = hashMap_index_getDub_bound(template_index, key, qseq + i, q_len - i, &j, -min, -max);
				if(value == 0) {
					i += j;
				}
			}
			
			if(prev_index < value) {
				value--;
				/* move position to start of kmer */
				i -= (kmersize - 1);
				if(value < Stat.pos) {
					if(Stat.pos != t_len) {
						/* should not happen */
						Stat.len = 1;
						Stat.score = 0;
						Stat.gaps = 0;
						/* check if any alternative alignment were better */
						if(0 < recStat.score) {
							/* move pointers */
							free(aligned->t);
							free(aligned->s);
							free(aligned->q);
							aligned->t = reAligned->t;
							aligned->s = reAligned->s;
							aligned->q = reAligned->q;
							free(reAligned);
							return recStat;
						} else if(recStat.score != 0) {
							free(reAligned->t);
							free(reAligned->s);
							free(reAligned->q);
							free(reAligned);
						}
						return Stat;
					}
					Stat.pos = value;
					
					/* Get intervals in query and template to align */
					if(value - prev_index > ((i - prev) << 1)) { // big leading template gap, cut down
						t_s = (value - ((i - prev) << 1));
						t_e = value;
						q_s = prev;
						q_e = i;
					} else if (i - prev > ((value - prev_index) << 1)) { // big leading query gap, cut down
						t_s = prev_index;
						t_e = value;
						q_s = i - ((value - prev_index) << 1);
						q_e = i;
					} else { // small leading gap
						t_s = prev_index;
						t_e = value;
						q_s = prev;
						q_e = i;
					}
					
					/* align leading gap */
					if(t_e - t_s > 0 && q_e - q_s > 0) {
						
						NWstat = NW(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align);
						
						/* trim leading gaps */
						bias = 0;
						if(t_s == 0) {
							while(bias < NWstat.len && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
								if(Frag_align->t[bias] == 5) {
									NWstat.gaps--;
								}
								bias++;
							}
							NWstat.len -= bias;
							/*if(bias) {
								NWstat.score -= (W1 + (bias - 1) * U);
							}*/
						}
						
						memcpy(aligned->t, Frag_align->t + bias, NWstat.len);
						memcpy(aligned->s, Frag_align->s + bias, NWstat.len);
						memcpy(aligned->q, Frag_align->q + bias, NWstat.len);
						Stat.pos -= (NWstat.len - NWstat.gaps);
						Stat.score = NWstat.score;
						Stat.len = NWstat.len;
						Stat.gaps = NWstat.gaps;
					} else {
						Stat.score = 0;
					}
				} else {
					/* piece seed-extends together */
					if(value - prev_index > (q_len << 1)) {
						/* gap is too big to give a positive score */
						Stat.score = 0;
						Stat.len = 1;
						Stat.gaps = 0;
						aligned->s[0] = 0;
						/* check if any alternative alignment were better */
						if(0 < recStat.score) {
							/* move pointers */
							free(aligned->t);
							free(aligned->s);
							free(aligned->q);
							aligned->t = reAligned->t;
							aligned->s = reAligned->s;
							aligned->q = reAligned->q;
							free(reAligned);
							return recStat;
						} else if(recStat.score != 0) {
							free(reAligned->t);
							free(reAligned->s);
							free(reAligned->q);
							free(reAligned);
						}
						return Stat;
					}
					
					/* get positions between seed-extends */
					t_s = prev_index;
					t_e = value;
					q_s = prev;
					q_e = i;
					
					NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align);
					memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
					memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
					memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
					Stat.score += NWstat.score;
					Stat.len += NWstat.len;
					Stat.gaps += NWstat.gaps;
				}
				
				/* Expand from matching k-mer */
				memcpy(aligned->t + Stat.len, qseq + i, kmersize);
				memset(aligned->s + Stat.len, '|', kmersize);
				memcpy(aligned->q + Stat.len, qseq + i, kmersize);
				
				Stat.len += kmersize;
				for(j = i + kmersize - 1; j >= i; j--) {
					Stat.score += d[qseq[j]][qseq[j]];
				}
				
				/* update positions */
				i += kmersize;
				prev_index = value + kmersize;
				
				/* extend */
				stop = 0;
				while(i < end && prev_index < t_len && !stop) {
					nuc = getNuc(template_index->seq, prev_index);
					if(qseq[i] == nuc) {
						Stat.score += d[nuc][nuc];
						aligned->t[Stat.len] = nuc;
						aligned->s[Stat.len] = '|';
						aligned->q[Stat.len] = nuc;
						Stat.len++;
						prev_index++;
						i++;
					} else {
						stop = 1;
					}
				}
				
				/* update positions */
				prev = i;
				
				if(i < end - kmersize) {
					key = makeKmer(qseq, i, kmersize - 1);
					i += (kmersize - 1);
				} else {
					i = end + 1;
				}
				
			} else {
				/* alternative alignment place */
				if(0 < value && value < prev_index && recStat.score == 0) {
					reAligned = malloc(sizeof(struct aln));
					if(!reAligned) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					reAligned->t = malloc((delta + 1) << 1);
					reAligned->s = malloc((delta + 1) << 1);
					reAligned->q = malloc((delta + 1) << 1);
					if(!reAligned->t || !reAligned->s || !reAligned->q) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					recStat = KMA(template_name, qseq, q_len, reAligned, Frag_align, i, min, max);
					if(recStat.score <= 0) {
						free(reAligned->t);
						free(reAligned->s);
						free(reAligned->q);
						free(reAligned);
						recStat.score = 0;
					}
				}
				i++;
			}
		}
		i = end + 1;
	}
	
	/* No valid ankers were found */
	if(prev_index == 0 || prev == 0) {
		/* No best mapping position found */
		Stat.score = 0;
		Stat.len = 1;
		aligned->s[0] = 0;
		/* check if any alternative alignment were better */
		if(0 < recStat.score) {
			/* move pointers */
			free(aligned->t);
			free(aligned->s);
			free(aligned->q);
			aligned->t = reAligned->t;
			aligned->s = reAligned->s;
			aligned->q = reAligned->q;
			free(reAligned);
			return recStat;
		} else if(recStat.score != 0) {
			free(reAligned->t);
			free(reAligned->s);
			free(reAligned->q);
			free(reAligned);
		}
		return Stat;
	}
	
	/* Get intervals in query and template to align */
	if(t_len - prev_index > ((q_len - prev) << 1)) { // big trailing template gap, cut down
		t_s = prev_index;
		t_e = t_s + ((q_len - prev) << 1);
		q_s = prev;
		q_e = q_len;
	} else if (q_len - prev > ((t_len - prev_index) << 1)) { // big leading query gap, cut down
		t_s = prev_index;
		t_e = t_len;
		q_s = prev;
		q_e = q_s + ((t_len - prev_index) << 1);
	} else { // small leading gap
		t_s = prev_index;
		t_e = t_len;
		q_s = prev;
		q_e = q_len;
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		NWstat = NW(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align);
		
		/* trim trailing gaps */
		/*
		if(t_e == t_len) {
			bias = NWstat.len - 1;
			while(bias && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
				if(Frag_align->t[bias] == 5) {
					NWstat.gaps--;
				}
				bias--;
			}
			bias++;
			if(bias != NWstat.len) {
				NWstat.score -= (W1 + (NWstat.len - bias) * U);
				NWstat.len = bias;
			}
		}
		*/
		memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
		memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
		memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
		Stat.score += NWstat.score;
		Stat.len += NWstat.len;
		Stat.gaps += NWstat.gaps;
	}
	aligned->s[Stat.len] = 0;
	
	/* check if any alternative alignment were better */
	if(Stat.score < recStat.score && 0 < recStat.score) {
		/* move pointers */
		free(aligned->t);
		free(aligned->s);
		free(aligned->q);
		aligned->t = reAligned->t;
		aligned->s = reAligned->s;
		aligned->q = reAligned->q;
		free(reAligned);
		return recStat;
	} else if(recStat.score != 0) {
		free(reAligned->t);
		free(reAligned->s);
		free(reAligned->q);
		free(reAligned);
	}
	return Stat;
}

struct alnScore KMA_score(const int template_name, const char *qseq, int q_len, const struct compDNA *qseq_comp, int qBias) {
	
	int i, j, k, prev, prev_index, stop, t_len, value, end, t_s, t_e, q_s, q_e;
	struct alnScore Stat, NWstat, recStat;
	struct hashMap_index *template_index;
	
	/* Extract indexes and template sequence */
	template_index = templates_index[template_name];
	t_len = template_lengths[template_name];
	
	/* Align the query and the template */
	recStat.score = 0;
	prev = 0;
	prev_index = 0;
	Stat.len = 0;
	Stat.score = 0;
	Stat.pos = t_len;
	Stat.gaps = 0;
	j = qBias;
	for(i = 1; i <= qseq_comp->N[0]; i++) {
		end = qseq_comp->N[i] - kmersize + 1;
		while(j < end && prev_index < t_len) {
			/* capture dublet */
			if((value = hashMap_index_get(template_index, getKmer(qseq_comp->seq, j))) < 0) {
				value = hashMap_index_getDub(template_index, getKmer(qseq_comp->seq, j), qseq + j, q_len - j, &k);
				if(value == 0) {
					j += k;
				}
			}
			
			if(prev_index < value) {
				value--;
				if(value < Stat.pos) {
					if(Stat.pos != t_len) {
						/* should not happen */
						Stat.len = 1;
						Stat.score = 0;
						Stat.gaps = 0;
						/* check if any alternative alignment were better */
						if(0 < recStat.score) {
							return recStat;
						}
						return Stat;
					}
					Stat.pos = value;
					
					/* Get intervals in query and template to align */
					if(value - prev_index > ((j - prev) << 1)) { // big leading template gap, cut down
						t_s = (value - ((j - prev) << 1));
						t_e = value;
						q_s = prev;
						q_e = j;
					} else if (j - prev > ((value - prev_index) << 1)) { // big leading query gap, cut down
						t_s = prev_index;
						t_e = value;
						q_s = j - ((value - prev_index) << 1);
						q_e = j;
					} else { // small leading gap
						t_s = prev_index;
						t_e = value;
						q_s = prev;
						q_e = j;
					}
					
					/* align leading gap */
					if(t_e - t_s > 0 && q_e - q_s > 0) {
						
						NWstat = NW_score(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e);
						
						Stat.pos -= (NWstat.len - NWstat.gaps); // NWstat.gaps: t_gap++, q_gap--
						Stat.score = NWstat.score;
						Stat.len = NWstat.len;
						Stat.gaps = NWstat.gaps;
					} else {
						Stat.score = 0;
					}
				} else {
					/* piece seed-extends together */
					if(value - prev_index > (q_len << 1)) {
						/* gap is too big to give a positive score */
						Stat.score = 0;
						Stat.len = 1;
						Stat.gaps = 0;
						/* check if any alternative alignment were better */
						if(0 < recStat.score) {
							return recStat;
						}
						return Stat;
					}
					
					/* get positions between seed-extends */
					t_s = prev_index;
					t_e = value;
					q_s = prev;
					q_e = j;
					
					NWstat = NW_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e);
					Stat.score += NWstat.score;
					Stat.len += NWstat.len;
					Stat.gaps += NWstat.gaps;
				}
				/* Expand from matching k-mer */
				Stat.len += kmersize;
				for(k = j + kmersize - 1; k >= j; k--) {
					Stat.score += d[qseq[k]][qseq[k]];
				}
				
				/* update positions */
				j += kmersize;
				prev_index = value + kmersize;
				
				/* extend */
				stop = 0;
				while(j < qseq_comp->N[i] && prev_index < t_len && !stop) {
					if(qseq[j] == getNuc(template_index->seq, prev_index)) {
						Stat.score += d[qseq[j]][qseq[j]];
						Stat.len++;
						prev_index++;
						j++;
					} else {
						stop = 1;
					}
				}
				/* update positions */
				prev = j;
			} else {
				/* alternative alignment */
				if(0 < value && value <= prev_index && recStat.score == 0) {
					recStat = KMA_score(template_name, qseq, q_len, qseq_comp, j);
				}
				j++;
			}
		}
		if(j < (qseq_comp->N[i] + 1)) {
			j = qseq_comp->N[i] + 1;
		}
	}
	
	/* No valid ankers were found */
	if(prev_index == 0 || prev == 0) {
		/* No best mapping position found */
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		/* check if any alternative alignment were better */
		if(0 < recStat.score) {
			return recStat;
		}
		return Stat;
	}
	
	/* Get intervals in query and template to align */
	if(t_len - prev_index > ((q_len - prev) << 1)) { // big trailing template gap, cut down
		t_s = prev_index;
		t_e = t_s + ((q_len - prev) << 1);
		q_s = prev;
		q_e = q_len;
	} else if (q_len - prev > ((t_len - prev_index) << 1)) { // big leading query gap, cut down
		t_s = prev_index;
		t_e = t_len;
		q_s = prev;
		q_e = q_s + ((t_len - prev_index) << 1);
	} else { // small leading gap
		t_s = prev_index;
		t_e = t_len;
		q_s = prev;
		q_e = q_len;
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		NWstat = NW_score(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e);
		Stat.score += NWstat.score;
		Stat.len += NWstat.len;
		Stat.gaps += NWstat.gaps;
	}
	
	/* check if any alternative alignment were better */
	if(Stat.score < recStat.score) {
		return recStat;
	}
	
	return Stat;
}

void assemble_KMA(struct assem *aligned_assem, int template, FILE **files, int file_count, FILE *frag_out, FILE *matrix_out, char *outputfilename, struct aln *aligned, struct aln *gap_align, char *qseq, char *header) {
	
	int i, j, t_len, q_len, aln_len, asm_len, start, end, bias, myBias, gaps;
	int read_score, depthUpdate, bestBaseScore, pos, file_len, bestScore;
	int nextTemplate, file_i, max_asmlen, nextGap, stats[4], *assemNext;
	unsigned depth, coverScore;
	short unsigned (*assembly)[6];
	double q_value, p_value, score;
	char bestNuc;
	FILE *file;
	struct alnScore alnStat;
	
	file_len = strlen(outputfilename);
	
	/* Allocate assembly arrays */
	t_len = template_lengths[template];
	asm_len = t_len;
	max_asmlen = t_len << 1;
	nextGap = t_len;
	assemNext = malloc(max_asmlen * sizeof(int));
	assembly = calloc(max_asmlen, 6 * sizeof(short unsigned));
	if(!assembly || !assemNext) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	/* init relative next */
	for(i = 0, j = 1; i < t_len; i++, j++) {
		assemNext[i] = j;
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
				fread(qseq, 1, q_len, file);
				
				/*
				stats[0] = 1;
				fread(&stats[1], sizeof(int), 1, file);
				fread(&stats[2], sizeof(int), 1, file);
				fread(&stats[3], sizeof(int), 1, file);
				*/
				fread(stats, sizeof(int), 4, file);
				fread(&i, sizeof(int), 1, file);
				fread(header, i, 1, file);
				
				/* Update assembly with read */
				/* Start with alignment */
				alnStat = KMA(template, qseq, q_len, aligned, gap_align, 0, stats[2], stats[3]);
				//alnStat = KMA(template, qseq, q_len, aligned, gap_align, 0, 0, t_len);
				
				/* get read score */
				aln_len = alnStat.len;
				start = alnStat.pos;
				end = start + aln_len - alnStat.gaps;
				
				/* Get normed score */
				read_score = alnStat.score;
				if(aln_len > 0) {
					score = 1.0 * read_score / aln_len;
				} else {
					score = 0;
				}
				
				if(read_score > kmersize && score > scoreT) {
					
					stats[1] = read_score;
					stats[2] = start;
					stats[3] = end;
					
					/* Update backbone and counts */
					i = 0;
					pos = start;
					while(i < aln_len) {
						if(aligned->t[i] == 5) { // Template gap, insertion
							if(pos >= t_len) {
								assembly[pos][aligned->q[i]]++;
								i++;
								pos = assemNext[pos];
							} else {
								/* get estimate for non insertions */
								myBias = 0;
								for(j = 0; j < 6; j++) {
									myBias += assembly[pos][j];
								}
								if(myBias > 0) {
									myBias--;
								}
								
								/* find position of insertion */
								gaps = pos;
								pos--;
								while(assemNext[pos] != gaps) {
									pos = assemNext[pos];
								}
								
								while(i < aln_len && aligned->t[i] == 5) {
									assemNext[pos] = nextGap;
									pos = assemNext[pos];
									assemNext[pos] = gaps;
									nextGap++;
									for(j = 0; j < 5; j++) {
										assembly[pos][j] = 0;
									}
									assembly[pos][5] = myBias;
									assembly[pos][aligned->q[i]]++;
									
									i++;
									if(nextGap == max_asmlen) {
										max_asmlen += t_len;
										assembly = realloc(assembly, max_asmlen * 6 * sizeof(short unsigned));
										assemNext = realloc(assemNext, max_asmlen * sizeof(int));
										if(!assembly || !assemNext) {
											fprintf(stderr, "OOM\n");
											exit(-1);
										}
									}
									asm_len++;
								}
								pos = assemNext[pos];
							}
						} else if(pos >= t_len) { // Old template gap, not present in this read
							assembly[pos][5]++;
							pos = assemNext[pos];
						} else { // (Match, mismatch) and (Query gap, deletion)
							assembly[pos][aligned->q[i]]++;
							i++;
							pos = assemNext[pos];
						}
					}
					
					/* Convert fragment */
					for(i = 0; i < q_len; i++) {
						 qseq[i] = bases[qseq[i]];
					}
					qseq[q_len] = 0;
					
					/* Save fragment */
					fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq, stats[0], stats[1], stats[2], stats[3], template_names[template], header);
				}
			} else if(nextTemplate == -1) {
				fclose(file);
				files[file_i] = 0;
				file_i++;
			} else if(nextTemplate < template) {
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len + 4 * sizeof(int), SEEK_CUR);
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len, SEEK_CUR);
			} else {
				/* Move pointer back */
				fseek(file, (-1) * sizeof(int), SEEK_CUR);
				file_i++;
			}
		} else {
			file_i++;
		}
	}
	
	/* Make consensus assembly by majority voting */
	/* Pepare and make alignment on consensus */
	if(aligned_assem->size <= asm_len) {
		aligned_assem->size = asm_len + 1;
		aligned_assem->t = malloc(asm_len + 1);
		aligned_assem->s = malloc(asm_len + 1);
		aligned_assem->q = malloc(asm_len + 1);
		if(!aligned_assem->t || !aligned_assem->s || !aligned_assem->q) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	/* Call nucleotides for the consensus */
	i = 0;
	pos = 0;
	depth = 0;
	while(i < asm_len) {
		/* call template */
		if(pos < t_len) {
			aligned_assem->t[i] = bases[getNuc(templates_index[template]->seq, pos)]; 
		} else {
			aligned_assem->t[i] = '-';
		}
		
		/* call query */
		bestNuc = 5;
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; j++) {
			if(bestScore < assembly[pos][j]) {
				bestScore = assembly[pos][j];
				bestNuc = j;
			}
			depthUpdate += assembly[pos][j];
		}
		bestNuc = bases[bestNuc];
		
		/* check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; j++) {
					if(bestBaseScore < assembly[pos][j]) {
						bestBaseScore = assembly[pos][j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[pos][5];
		}
		
		/* determine base at current position */
		if(depthUpdate == 0 || bestNuc == '-') {
			aligned_assem->q[i] = '-';
		} else {
			/* Use MC Neymars test to test significance of the base call */
			q_value = pow(bestScore - (depthUpdate - bestScore), 2) / depthUpdate;
			p_value = p_chisqr(q_value);
			if(p_value <= evalue && bestScore > (depthUpdate - bestScore)) {
				aligned_assem->q[i] = bestNuc;
			} else {
				aligned_assem->q[i] = tolower(bestNuc);
			}
		}
		if(bestNuc != '-') {
			depth += depthUpdate;
		}
		
		i++;
		pos = assemNext[pos];
	}
	
	/* print matrix */
	if(print_matrix) {
		fprintf(matrix_out, "#%s\n", template_names[template]);
		for(i = 0; i < asm_len; i++) {
			fprintf(matrix_out, "%c", aligned_assem->t[i]);
			for(j = 0; j < 6; j++) {
				fprintf(matrix_out, "\t%hu", assembly[i][j]);
			}
			fprintf(matrix_out, "\n");
		}
		fprintf(matrix_out, "\n");
	}
	
	/* Trim alignment on consensus */
	coverScore = 0;
	bias = 0;
	for(i = 0; i < asm_len; i++) {
		if(aligned_assem->t[i] == '-' && aligned_assem->q[i] == '-') {
			bias++;
		} else {
			aligned_assem->t[i - bias] = aligned_assem->t[i];
			aligned_assem->q[i - bias] = aligned_assem->q[i];
			if(tolower(aligned_assem->t[i]) == tolower(aligned_assem->q[i])) {
				aligned_assem->s[i - bias] = '|';
				coverScore++;
			} else {
				aligned_assem->s[i - bias] = '_';
			}
		}
	}
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->t[asm_len - bias] = 0;
	aligned_assem->s[asm_len - bias] = 0;
	aligned_assem->q[asm_len - bias] = 0;
	
	/* clean */
	free(assembly);
	free(assemNext);
}

void assemble_KMA_dense(struct assem *aligned_assem, int template, FILE **files, int file_count, FILE *frag_out, FILE *matrix_out, char *outputfilename, struct aln *aligned, struct aln *gap_align, char *qseq, char *header) {
	
	int i, j, t_len, q_len, aln_len, start, end, file_len, file_i, stats[4];
	int pos, read_score, bestScore, depthUpdate, bestBaseScore, nextTemplate;
	unsigned depth, coverScore;
	short unsigned (*assembly)[6];
	double q_value, p_value, score;
	char bestNuc;
	FILE *file;
	struct alnScore alnStat;
	
	file_len = strlen(outputfilename);
	/* Allocate assembly arrays */
	t_len = template_lengths[template];
	assembly = calloc(t_len + 1, 6 * sizeof(short unsigned));
	if(aligned_assem->size <= t_len) {
		aligned_assem->size = t_len + 1;
		aligned_assem->t = malloc(t_len + 1);
		aligned_assem->s = malloc(t_len + 1);
		aligned_assem->q = malloc(t_len + 1);
	}
	if(!assembly || !aligned_assem->t || !aligned_assem->s || !aligned_assem->q) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	/* cpy template seq */
	for(i = 0; i < t_len; i++) {
		aligned_assem->t[i] = getNuc(templates_index[template]->seq, i);
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
				fread(qseq, 1, q_len, file);
				
				fread(stats, sizeof(int), 4, file);
				
				fread(&i, sizeof(int), 1, file);
				fread(header, i, 1, file);
				
				/* Update assembly with read */
				/* Start with alignment */
				alnStat = KMA(template, qseq, q_len, aligned, gap_align, 0, stats[2], stats[3]);
				
				/* get read score */
				aln_len = alnStat.len;
				start = alnStat.pos;
				end = start + aln_len - alnStat.gaps;
				
				/* Get normed score */
				read_score = alnStat.score;
				if(aln_len > 0) {
					score = 1.0 * read_score / aln_len;
				} else {
					score = 0;
				}
				
				if(read_score > kmersize && score > scoreT) {
					
					stats[1] = read_score;
					stats[2] = start;
					stats[3] = end;
					
					/* Update backbone and counts */
					for(i = 0, pos = start; i < aln_len; i++) {
						if(aligned->t[i] == aligned_assem->t[pos]) {
							assembly[pos][aligned->q[i]]++;
							pos++;
						}
					}
					
					/* Convert fragment */
					for(i = 0; i < q_len; i++) {
						 qseq[i] = bases[qseq[i]];
					}
					qseq[q_len] = 0;
					
					/* Save fragment */
					/*for(i = 0; i < q_len; i++) {
						fprintf(frag_out, "%c", bases[qseq[i]]);
					}*/
					fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq, stats[0], stats[1], stats[2], stats[3], template_names[template], header);
				}
			} else if (nextTemplate == -1) {
				fclose(file);
				files[file_i] = 0;
				file_i++;
			} else if(nextTemplate < template) {
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len + 4 * sizeof(int), SEEK_CUR);
				fread(&q_len, sizeof(int), 1, file);
				fseek(file, q_len, SEEK_CUR);
			} else {
				/* Move pointer back */
				fseek(file, (-1) * sizeof(int), SEEK_CUR);
				file_i++;
			}
		} else {
			file_i++;
		}
	}
	
	/* Make consensus assembly by majority voting */
	depth = 0;
	coverScore = 0;
	for(i = 0; i < t_len; i++) {
		/* call template */
		aligned_assem->t[i] = bases[aligned_assem->t[i]];
		
		/* call query */
		bestNuc = 5;
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; j++) {
			if(bestScore < assembly[i][j]) {
				bestScore = assembly[i][j];
				bestNuc = j;
			}
			depthUpdate += assembly[i][j];
		}
		bestNuc = bases[bestNuc];
		
		/* Check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; j++) {
					if(bestBaseScore < assembly[i][j]) {
						bestBaseScore = assembly[i][j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[i][5];
		}
		
		/* determine base at current position */
		if(depthUpdate == 0 || bestNuc == '-') {
			aligned_assem->q[i] = '-';
		} else {
			/* Use MC Neymars test to test significance of the base call */
			q_value = pow(bestScore - (depthUpdate - bestScore), 2) / depthUpdate;
			p_value = p_chisqr(q_value);
			if(p_value <= evalue && bestScore > (depthUpdate - bestScore)) {
				aligned_assem->q[i] = bestNuc;
			} else {
				aligned_assem->q[i] = tolower(bestNuc);
			}
		}
		if(bestNuc != '-') {
			depth += depthUpdate;
		}
		
		if(tolower(aligned_assem->q[i]) == tolower(aligned_assem->t[i])) {
			aligned_assem->s[i] = '|';
			coverScore++;
		} else {
			aligned_assem->s[i] = '_';
		}
	}
	aligned_assem->t[t_len] = 0;
	aligned_assem->s[t_len] = 0;
	aligned_assem->q[t_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	
	/* print matrix */
	if(print_matrix && coverScore > 0) {
		fprintf(matrix_out, "#%s\n", template_names[template]);
		for(i = 0; i < t_len; i++) {
			fprintf(matrix_out, "%c", aligned_assem->t[i]);
			for(j = 0; j < 6; j++) {
				fprintf(matrix_out, "\t%hu", assembly[i][j]);
			}
			fprintf(matrix_out, "\n");
		}
		fprintf(matrix_out, "\n");
	}
	
	/* clean */
	free(assembly);
}

void update_Scores(char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, char *header, FILE *frag_out_raw) {
	
	int i, len;
	
	len = strlen(header) + 1;
	
	/* print read:  q_len, qseq, num, score, start, end, template */
	fwrite(&q_len, sizeof(int), 1, frag_out_raw);
	fwrite(qseq, 1, q_len, frag_out_raw);
	fwrite(&counter, sizeof(int), 1, frag_out_raw);
	fwrite(&score, sizeof(int), 1, frag_out_raw);
	fwrite(start, sizeof(int), counter, frag_out_raw);
	fwrite(end, sizeof(int), counter, frag_out_raw);
	fwrite(template, sizeof(int), counter, frag_out_raw);
	fwrite(&len, sizeof(int), 1, frag_out_raw);
	fwrite(header, 1, len, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { //Only one best match
		if(*template < 0) {
			template[0] = -template[0];
		}
		//alignment_scores[template[0]] += (1.0 * score) / (end[0] - start[0]);
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
		//uniq_alignment_scores[template[0]]++;
	} else {
		for(i = 0; i < counter; i++) {
			if(template[i] < 0) {
				template[i] = -template[i];
			}
			alignment_scores[template[i]] += score;
			//alignment_scores[template[i]] += (1.0 * score) / (end[i] - start[i]);
		}
	}
}

void runKMA(char *templatefilename, char *outputfilename, char *exePrev) {
	
	int i, j, t_i, tmp_template, tmp_tmp_template, file_len, best_read_score;
	int template, template_stop, query_stop, s_len, gaps, bestHits, t_len;
	int read_counter, start, end, pos, asm_len, q_len, frags, aln_len, rc_flag;
	int depthUpdate, coverScore, bias, matches, tmp_start, tmp_end;
	int bestTemplate, fragCount, fileCount, maxFrag, read_score, header_size;
	int *matched_templates, *bestTemplates, *best_start_pos, *best_end_pos;
	unsigned Nhits, template_tot_ulen, bestNum, *w_scores;
	long *index_indexes, *seq_indexes;
	char bestNuc, *outZipped, *qseq, *qseq_r, *header;
	double etta, tmp_score, score, bestScore, bestBaseScore, depth, q_value;
	double id, p_value, expected, cover, q_id, q_cover;
	FILE *index_in, *seq_in, *res_out, **template_fragments, *matrix_out;
	FILE *consensus_out, *alignment_out, *inputfile;
	FILE *frag_out_raw, *frag_in_raw, *frag_out, *frag_out_all;
	time_t t0, t1;
	struct aln *aligned, *gap_align;
	struct assem *aligned_assem;
	struct frag **alignFrags, *alignFrag;
	struct compDNA *qseq_comp, *qseq_r_comp;
	struct alnScore alnStat;
	
	/* open pipe */
	inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
	/* load databases */
	file_len = strlen(templatefilename);
	load_DBs_KMA(templatefilename);
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".index.b");
	index_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".seq.b");
	seq_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	if(!index_in || !seq_in) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist. 4\n");
		exit(-1);
	} else if(shm & 8) {
		templates_index[0] = alignLoad_shm_initial(templatefilename, file_len, seq_in, index_in);
		alignLoadPtr = &alignLoad_fly_shm;
		destroyPtr = &alignClean_shm;
	} else {
		fread(&kmersize, sizeof(int), 1, index_in);
	}
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	/* allocate stuff */
	file_len = strlen(outputfilename);
	outZipped = malloc(strlen("gunzip -c .frag_raw.gz") + file_len);
	qseq = malloc(delta + 1);
	qseq_r = malloc(delta + 1);
	qseq_comp = malloc(sizeof(struct compDNA));
	qseq_r_comp = malloc(sizeof(struct compDNA));
	header_size = 1024;
	header = malloc(header_size);
	index_indexes = malloc(DB_size * sizeof(long));
	seq_indexes = malloc(DB_size * sizeof(long));
	if(!qseq || !qseq_r || !header || !qseq_comp || !qseq_r_comp || !outZipped || !index_indexes || !seq_indexes) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	
	/* make file indexes of template indexing */
	*index_indexes = 0;
	*seq_indexes = 0;
	index_indexes[1] = sizeof(int);
	seq_indexes[1] = 0;
	for(i = 2; i < DB_size; i++) {
		index_indexes[i] = index_indexes[i - 1] + (template_lengths[i - 1] << 1) * sizeof(int);
		seq_indexes[i] = seq_indexes[i - 1] + ((template_lengths[i - 1] >> 5) + 1) * sizeof(long unsigned);
	}
	
	/* allocate matrcies for NW */
	NW_s = delta;
	D = malloc(NW_s * sizeof(int*));
	Q = malloc(NW_s * sizeof(int*));
	P = malloc(NW_s * sizeof(int*));
	E = malloc(NW_s * sizeof(int*));
	if(!D || !Q || !P || !E) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	for(i = 0; i < NW_s; i++) {
		D[i] = malloc(NW_s * sizeof(int));
		Q[i] = malloc(NW_s * sizeof(int));
		P[i] = malloc(NW_s * sizeof(int));
		E[i] = malloc(NW_s * sizeof(int));
		if(!D[i] || !Q[i] || !P[i] || !E[i]) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	
	/* etta = small value to avoid zero-divisio */
	etta = 1.0e-6;
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag.gz");
		sprintf(outZipped, "gzip -c > %s", outputfilename);
		frag_out = popen(outZipped, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".aln");
		alignment_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".fsa");
		consensus_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag_raw.b");
		frag_out_raw = fopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		if(print_matrix) {
			strcat(outputfilename, ".mat.gz");
			sprintf(outZipped, "gzip -c > %s", outputfilename);
			matrix_out = popen(outZipped, "w");
			if(!matrix_out) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				exit(-1);
			}
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(print_all) {
			strcat(outputfilename, ".frag_raw.gz");
			sprintf(outZipped, "gzip -c > %s", outputfilename);
			frag_out_all = popen(outZipped, "w");
			if(!frag_out_all) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				exit(-1);
			}
			outputfilename[file_len] = 0;
		} else {
			frag_out_all = 0;
		}
	} else {
		fprintf(stderr, "# No output file specified!\n");
		exit(-1);
	}
	if(!res_out || !frag_out || !consensus_out || !frag_out_raw || !alignment_out) {
		fprintf(stderr, "Could not create outputfiles\n");
		exit(-1);
	}
	fprintf(stderr, "# Running KMA.\n");
	t0 = clock();
	
	/* Get alignments */
	read_counter = 0;
	matched_templates = malloc(((DB_size + 1) << 1) * sizeof(int));
	best_start_pos = malloc((DB_size << 1) * sizeof(int));
	best_end_pos = malloc((DB_size << 1) * sizeof(int));
	bestTemplates = malloc((DB_size << 1) * sizeof(int));
	if(!matched_templates || !best_start_pos || !best_end_pos || !bestTemplates) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	while(getComp(qseq_comp, inputfile)) {
		q_len = qseq_comp->seqlen;
		
		if(q_len >= delta) {
			delta = q_len << 1;
			free(qseq);
			free(qseq_r);
			qseq = malloc(delta);
			qseq_r = malloc(delta);
			if(!qseq || !qseq_r) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		
		/* reverse complement seq */
		fread(&rc_flag, sizeof(int), 1, inputfile);
		if(rc_flag < 0) {
			if(qseq_comp->size > qseq_r_comp->size) {
				freeComp(qseq_r_comp);
				allocComp(qseq_r_comp, qseq_comp->size);
			}
			rc_comp(qseq_comp, qseq_r_comp);
			unCompDNA(qseq_r_comp, qseq_r);
			qseq_r_comp->N[0]++;
			qseq_r_comp->N[qseq_r_comp->N[0]] = q_len;
		}
		unCompDNA(qseq_comp, qseq);
		qseq_comp->N[0]++;
		qseq_comp->N[qseq_comp->N[0]] = q_len;
		
		/* Get number of matched templates */
		fread(&bestHits, sizeof(int), 1, inputfile);
		*matched_templates = bestHits;
		fread(matched_templates + 1, sizeof(int), bestHits, inputfile);
		
		/* get header */
		fread(&t_i, sizeof(int), 1, inputfile);
		if(t_i > header_size) {
			free(header);
			header_size = t_i << 1;
			header = malloc(header_size);
			if(!header) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		fread(header, t_i, 1, inputfile);
		header[t_i] = 0;
		
		bestScore = 0;
		best_read_score = 0;
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; t_i++) {
			template = matched_templates[t_i];
			/* check if index DB is loaded */
			if(template >= 0 && templates_index[template] == 0) {
				templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], seq_indexes[template], index_indexes[template]);
			} else if(template < 0 && templates_index[-template] == 0) {
				templates_index[-template] = alignLoadPtr(seq_in, index_in, template_lengths[-template], seq_indexes[-template], index_indexes[-template]);
			}
			
			/* align qseq */
			if(template < 0) {
				alnStat = KMA_score(-template, qseq_r, q_len, qseq_r_comp, 0);
			} else {
				alnStat = KMA_score(template, qseq, q_len, qseq_comp, 0);
			}
			
			/* get read score */
			aln_len = alnStat.len;
			start = alnStat.pos;
			end = start + aln_len - alnStat.gaps;
			
			/* penalty for non complete mapping */
			read_score = alnStat.score;
			read_score += (((start != 0) + (end != template_lengths[abs(template)])) * W1);
			
			/* Get normed score */
			if(aln_len > 0) {
				score = 1.0 * read_score / aln_len;
			} else {
				score = 0;
			}
			
			/* save best match(es) */
			if(read_score > kmersize && score > scoreT) {
				if(score > bestScore) { // save as best match
					bestScore = score;
					best_read_score = read_score;
					*bestTemplates = template;
					*best_start_pos = start;
					*best_end_pos = end;
					bestHits = 1;
				} else if(score == bestScore && read_score > best_read_score) { // save as best match
					bestScore = score;
					best_read_score = read_score;
					*bestTemplates = template;
					*best_start_pos = start;
					*best_end_pos = end;
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
			/*
			for(i = 0; i < q_len; i++) {
				fprintf(stdout, "%c", bases[qseq[i]]);
			}
			fprintf(stdout, "\n%d\t%d\t%d\t%s\n", best_read_score, *best_start_pos, *best_end_pos, header);
			*/
			update_Scores(qseq, q_len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		}
	}
	pclose(inputfile);
	fclose(frag_out_raw);
	free(matched_templates);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	
	t1 = clock();
	fprintf(stderr, "#\n# KMA mapping time\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select KMA alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(unsigned));
	if(!alignFrags || !w_scores) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	outputfilename[file_len] = 0;
	strcat(outputfilename, ".frag_raw.b");
	frag_in_raw = fopen(outputfilename, "rb");
	if(!frag_in_raw) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	outputfilename[file_len] = 0;
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	while(fread(&q_len, sizeof(int), 1, frag_in_raw)) {
		fread(qseq, 1, q_len, frag_in_raw);
		fread(&bestHits, sizeof(int), 1, frag_in_raw);
		fread(&read_score, sizeof(int), 1, frag_in_raw);
		
		fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
		fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
		fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		fread(&t_i, sizeof(int), 1, frag_in_raw);
		fread(header, 1, t_i, frag_in_raw);
		
		if(q_len > kmersize) {
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
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / template_lengths[tmp_template];
					if(tmp_score > bestScore) {
					//if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(tmp_score == bestScore) {
					//} else if(alignment_scores[tmp_template] == best_read_score) {
						//if(uniq_alignment_scores[tmp_template] > bestNum) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						//} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_score > bestScore) {
						} else if(alignment_scores[tmp_template] == best_read_score) {
						//} else if(tmp_score == bestScore) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if (uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq, q_len);
			}
			w_scores[bestTemplate] += read_score;
			
			/* dump frag info */
			if(alignFrags[bestTemplate] == 0) {
				alignFrags[bestTemplate] = malloc(sizeof(struct frag));
				alignFrag = alignFrags[bestTemplate];
				if(!alignFrag) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				alignFrag->q_len = q_len;
				alignFrag->qseq = malloc(q_len);
				alignFrag->header = strdup(header);
				if(!alignFrag->qseq || !alignFrag->header) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				memcpy(alignFrag->qseq, qseq, q_len);
				alignFrag->bestHits = bestHits;
				alignFrag->score = score;
				alignFrag->start = start;
				alignFrag->end = end;
				alignFrag->next = 0;
			} else {
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				alignFrag->q_len = q_len;
				alignFrag->qseq = malloc(q_len);
				alignFrag->header = strdup(header);
				if(!alignFrag->qseq || !alignFrag->header) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				memcpy(alignFrag->qseq, qseq, q_len);
				alignFrag->bestHits = bestHits;
				alignFrag->score = score;
				alignFrag->start = start;
				alignFrag->end = end;
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
			}
			
			fragCount++;
			if(fragCount >= maxFrag) {
				sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
				template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
				if(!template_fragments[fileCount]) {
					fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
					exit(-1);
				}
				outputfilename[file_len] = 0;
				fileCount++;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
				}
			}
			
			/* dump seq to all */
			if(frag_out_all) {
				for(i = 0; i < q_len; i++) {
					qseq[i] = bases[qseq[i]];
				}
				qseq[i] = 0;
				fprintf(frag_out_all, "%s\t%d\t%d\t%d", qseq, bestHits, read_score, *best_start_pos);
				for(i = 1; i < bestHits; i++) {
					fprintf(frag_out_all, ",%d", best_start_pos[i]);
				}
				fprintf(frag_out_all, "\t%d", *best_end_pos);
				for(i = 1; i < bestHits; i++) {
					fprintf(frag_out_all, ",%d", best_end_pos[i]);
				}
				fprintf(frag_out_all, "\t%d", *bestTemplates);
				for(i = 1; i < bestHits; i++) {
					fprintf(frag_out_all, ",%d", bestTemplates[i]);
				}
				fprintf(frag_out_all, "\t%s\n", header);
			}
		}
	}
	sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
	template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
	if(!template_fragments[fileCount]) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	outputfilename[file_len] = 0;
	fileCount++;
	fragCount = 0;
	fclose(frag_in_raw);
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(bestTemplates);
	free(qseq_r);
	if(frag_out_all) {
		pclose(frag_out_all);
	}
	/* remove frag_raw */
	strcat(outputfilename, ".frag_raw.b");
	remove(outputfilename);
	outputfilename[file_len] = 0;
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	
	/* Get expected values */
	template_tot_ulen = 0;
	Nhits = 0;
	for(i = 0; i < DB_size; i++) {
		Nhits += w_scores[i];
		template_tot_ulen += template_lengths[i];
	}
	
	/* Do local assemblies of fragments mapping to the same template */
	aligned = malloc(sizeof(struct aln));
	gap_align = malloc(sizeof(struct aln));
	if(!aligned || !gap_align) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	aligned_assem = malloc(sizeof(struct assem));
	aligned->t = malloc((delta + 1) << 1);
	aligned->s = malloc((delta + 1) << 1);
	aligned->q = malloc((delta + 1) << 1);
	gap_align->t = malloc((delta + 1) << 1);
	gap_align->s = malloc((delta + 1) << 1);
	gap_align->q = malloc((delta + 1) << 1);
	if(!aligned_assem || !aligned->t || !aligned->s || !aligned->q || !gap_align->t || !gap_align->s || !gap_align->q) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	aligned_assem->size = 0;
	aligned_assem->t = NULL;
	aligned_assem->s = NULL;
	aligned_assem->q = NULL;
	for(template = 1; template < DB_size; template++) {
		if(w_scores[template] > 0) {
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = (Nhits - read_score) * (1.0 * t_len) / (template_tot_ulen - t_len + etta);
			q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
			p_value  = p_chisqr(q_value);
			
			
			if((p_value <= evalue && score > expected) || ((1.0 * read_score / t_len) > scoreT)) {
				/* Do assembly */
				assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, matrix_out, outputfilename, aligned, gap_align, qseq, header);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = strlen(aligned_assem->q) - countChar(aligned_assem->q, '-');
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
				}
				
				if(id >= ID_t) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8u\t%8d\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_names[template], read_score, (int) expected, t_len, id, cover, q_id, q_cover, depth, q_value, p_value);
					
					/* print alignment */
					aln_len = strlen(aligned_assem->s);
					fprintf(alignment_out, "# %s\n", template_names[template]);
					for(i = 0; i < aln_len; i += 60) {
						fprintf(alignment_out, "%-10s\t", "template:");
						for(j = i; j < i + 60 && aligned_assem->t[j]; j++) {
							fprintf(alignment_out, "%c", aligned_assem->t[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "");
						for(j = i; j < i + 60 && aligned_assem->s[j]; j++) {
							fprintf(alignment_out, "%c", aligned_assem->s[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "query:");
						for(j = i; j < i + 60 && aligned_assem->q[j]; j++) {
							fprintf(alignment_out, "%c", aligned_assem->q[j]);
						}
						fprintf(alignment_out, "\n\n");
					}
					
					/* Print consensus */
					fprintf(consensus_out, ">%s\n", template_names[template]);
					if(ref_fsa) {
						for(i = 0; i < aln_len; i += 60) {
							for(j = i; j < i + 60 && aligned_assem->q[j]; j++) {
								if(aligned_assem->q[j]) {
									if(aligned_assem->q[j] != '-') {
										fprintf(consensus_out, "%c", aligned_assem->q[j]);
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
							for(j = i; j < i + 60 + bias && aligned_assem->q[j]; j++) {
								if(aligned_assem->q[j] && aligned_assem->q[j] != '-') {
									fprintf(consensus_out, "%c", aligned_assem->q[j]);
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
	free(qseq);
	free(header);
	if(aligned_assem->size != 0) {
		free(aligned_assem->t);
		free(aligned_assem->s);
		free(aligned_assem->q);
	}
	free(aligned_assem);
	
	/* Close files */
	fclose(index_in);
	fclose(seq_in);
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	pclose(frag_out);
	if(matrix_out) {
		pclose(matrix_out);
	}
	
	/* remove tmp files */
	for(i = 0; i < fileCount; i++) {
		sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", i);
		remove(outputfilename);
		outputfilename[file_len] = 0;
	}
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
}

void runKMA_MEM(char *templatefilename, char *outputfilename, char *exePrev) {
	/* runKMA_MEM is a memory saving version of runKMA,
	   at the cost it chooses best templates based on kmers
	   instead of alignment score. */
	
	int i, j, t_i, tmp_template, tmp_tmp_template, file_len, best_read_score;
	int template, template_stop, query_stop, s_len, gaps, bestHits, t_len;
	int read_counter, start, end, pos, asm_len, q_len, frags, aln_len, rc_flag;
	int depthUpdate, coverScore, bias, matches, tmp_start, tmp_end;
	int bestTemplate, fragCount, fileCount, maxFrag, read_score, header_size;
	int *matched_templates, *bestTemplates, *best_start_pos, *best_end_pos;
	unsigned Nhits, template_tot_ulen, bestNum, *w_scores;
	double etta, tmp_score, score, bestScore, bestBaseScore, depth, q_value;
	double p_value, expected, id, cover, q_id, q_cover;
	char bestNuc, *outZipped, *qseq, *qseq_r, *header;
	FILE *index_in, *seq_in, *res_out, **template_fragments, *matrix_out;
	FILE *consensus_out, *alignment_out, *inputfile;
	FILE *frag_out_raw, *frag_in_raw, *frag_out, *frag_out_all;
	time_t t0, t1;
	struct aln *aligned, *gap_align;
	struct assem *aligned_assem;
	struct frag **alignFrags, *alignFrag;
	struct compDNA *qseq_comp, *qseq_r_comp;
	struct alnScore alnStat;
	
	/* open pipe */
	inputfile = popen(exePrev, "r");
	if(!inputfile) {
		fprintf(stderr, "File corruption\n");
		exit(-1);
	}
	
	/* load databases */
	file_len = strlen(templatefilename);
	load_DBs_KMA(templatefilename);
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".index.b");
	index_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".seq.b");
	seq_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	if(!index_in || !seq_in) {
		fprintf(stderr, "Wrong format of DB, or DB does not exist. 5\n");
		exit(-1);
	}
	fread(&kmersize, sizeof(int), 1, index_in);
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	/* allocate stuff */
	file_len = strlen(outputfilename);
	outZipped = malloc(strlen("gunzip -c .frag_raw.gz") + file_len);
	qseq = malloc(delta + 1);
	qseq_r = malloc(delta + 1);
	qseq_comp = malloc(sizeof(struct compDNA));
	qseq_r_comp = malloc(sizeof(struct compDNA));
	header_size = 1024;
	header = malloc(header_size);
	if(!qseq || !qseq_r || !header || !qseq_comp || !qseq_r_comp || !outZipped) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	
	/* allocate matrcies for NW */
	NW_s = delta;
	D = malloc(NW_s * sizeof(int*));
	Q = malloc(NW_s * sizeof(int*));
	P = malloc(NW_s * sizeof(int*));
	E = malloc(NW_s * sizeof(int*));
	if(!D || !Q || !P || !E) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	for(i = 0; i < NW_s; i++) {
		D[i] = malloc(NW_s * sizeof(int));
		Q[i] = malloc(NW_s * sizeof(int));
		P[i] = malloc(NW_s * sizeof(int));
		E[i] = malloc(NW_s * sizeof(int));
		if(!D[i] || !Q[i] || !P[i] || !E[i]) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	
	/* etta = small value to avoid zero-divisio */
	etta = 1.0e-6;
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag.gz");
		sprintf(outZipped, "gzip -c > %s", outputfilename);
		frag_out = popen(outZipped, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".aln");
		alignment_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".fsa");
		consensus_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag_raw.b");
		frag_out_raw = fopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		if(print_matrix) {
			strcat(outputfilename, ".mat.gz");
			sprintf(outZipped, "gzip -c > %s", outputfilename);
			matrix_out = popen(outZipped, "w");
			if(!matrix_out) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				exit(-1);
			}
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(print_all) {
			strcat(outputfilename, ".frag_raw.gz");
			sprintf(outZipped, "gzip -c > %s", outputfilename);
			frag_out_all = popen(outZipped, "w");
			if(!frag_out_all) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				exit(-1);
			}
			outputfilename[file_len] = 0;
		} else {
			frag_out_all = 0;
		}
	} else {
		fprintf(stderr, "# No output file specified!\n");
		exit(-1);
	}
	if(!res_out || !frag_out || !consensus_out || !frag_out_raw || !alignment_out) {
		fprintf(stderr, "Could not create outputfiles\n");
		exit(-1);
	}
	
	fprintf(stderr, "# Collecting k-mer scores.\n");
	t0 = clock();
	
	/* Get alignments */
	read_counter = 0;
	matched_templates = malloc(((DB_size + 1) << 1) * sizeof(int));
	best_start_pos = calloc((DB_size << 1), sizeof(int));
	best_end_pos = malloc((DB_size << 1) * sizeof(int));
	bestTemplates = malloc((DB_size << 1) * sizeof(int));
	if(!matched_templates || !best_start_pos || !best_end_pos || !bestTemplates) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	while(getComp(qseq_comp, inputfile)) {
		q_len = qseq_comp->seqlen;
		
		if(q_len >= delta) {
			delta = q_len << 1;
			free(qseq);
			free(qseq_r);
			qseq = malloc(delta);
			qseq_r = malloc(delta);
			if(!qseq || !qseq_r) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		
		/* reverse complement seq */
		fread(&rc_flag, sizeof(int), 1, inputfile);
		if(rc_flag < 0) {
			best_read_score = -rc_flag;
			if(qseq_comp->size > qseq_r_comp->size) {
				freeComp(qseq_r_comp);
				allocComp(qseq_r_comp, qseq_comp->size);
			}
			rc_comp(qseq_comp, qseq_r_comp);
			unCompDNA(qseq_r_comp, qseq_r);
		} else {
			best_read_score = rc_flag;
		}
		unCompDNA(qseq_comp, qseq);
		
		/* Get number of matched templates */
		fread(&bestHits, sizeof(int), 1, inputfile);
		*matched_templates = bestHits;
		fread(matched_templates + 1, sizeof(int), bestHits, inputfile);
		
		/* get header */
		fread(&i, sizeof(int), 1, inputfile);
		if(i > header_size) {
			free(header);
			header_size = i << 1;
			header = malloc(header_size);
			if(!header) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		fread(header, i, 1, inputfile);
		header[i] = 0;
		
		for(i = 1, bestHits = 0; i <= *matched_templates; i++, bestHits++) {
			bestTemplates[bestHits] = matched_templates[i];
			//best_start_pos[bestHits] = 0;
			best_end_pos[bestHits] = template_lengths[abs(matched_templates[i])];
		}
		update_Scores(qseq, q_len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		
		/* dump seq to all */
		if(frag_out_all) {
			for(i = 0; i < q_len; i++) {
				qseq[i] = bases[qseq[i]];
			}
			qseq[i] = 0;
			fprintf(frag_out_all, "%s\t%d\t%d\t%d", qseq, bestHits, best_read_score, *best_start_pos);
			for(i = 1; i < bestHits; i++) {
				fprintf(frag_out_all, ",%d", best_start_pos[i]);
			}
			fprintf(frag_out_all, "\t%d", *best_end_pos);
			for(i = 1; i < bestHits; i++) {
				fprintf(frag_out_all, ",%d", best_end_pos[i]);
			}
			fprintf(frag_out_all, "\t%d", *bestTemplates);
			for(i = 1; i < bestHits; i++) {
				fprintf(frag_out_all, ",%d", bestTemplates[i]);
			}
			fprintf(frag_out_all, "\t%s\n", header);
		}
	}
	pclose(inputfile);
	fclose(frag_out_raw);
	free(matched_templates);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(frag_out_all) {
		pclose(frag_out_all);
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Time for score collecting:\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select k-mer alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(unsigned));
	if(!alignFrags || !w_scores) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	outputfilename[file_len] = 0;
	strcat(outputfilename, ".frag_raw.b");
	frag_in_raw = fopen(outputfilename, "rb");
	if(!frag_in_raw) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	outputfilename[file_len] = 0;
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	while(fread(&q_len, sizeof(int), 1, frag_in_raw)) {
		fread(qseq, 1, q_len, frag_in_raw);
		fread(&bestHits, sizeof(int), 1, frag_in_raw);
		fread(&read_score, sizeof(int), 1, frag_in_raw);
		
		fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
		fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
		fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
		
		fread(&t_i, sizeof(int), 1, frag_in_raw);
		fread(header, t_i, 1, frag_in_raw);
		
		if(q_len > kmersize) {
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
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					if(tmp_score > bestScore) {
					//if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					//} else if(alignment_scores[tmp_template] == best_read_score) {
					} else if(tmp_score == bestScore) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq, q_len);
			}
			w_scores[bestTemplate] += read_score;
			
			/* dump frag info */
			if(alignFrags[bestTemplate] == 0) {
				alignFrags[bestTemplate] = malloc(sizeof(struct frag));
				alignFrag = alignFrags[bestTemplate];
				if(!alignFrag) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				alignFrag->q_len = q_len;
				alignFrag->qseq = malloc(q_len);
				alignFrag->header = strdup(header);
				if(!alignFrag->qseq || !alignFrag->header) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				memcpy(alignFrag->qseq, qseq, q_len);
				alignFrag->bestHits = bestHits;
				alignFrag->score = (int) score;
				alignFrag->start = start;
				alignFrag->end = end;
				alignFrag->next = 0;
			} else {
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				alignFrag->q_len = q_len;
				alignFrag->qseq = malloc(q_len);
				alignFrag->header = strdup(header);
				if(!alignFrag->qseq || !alignFrag->header) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				memcpy(alignFrag->qseq, qseq, q_len);
				alignFrag->bestHits = bestHits;
				alignFrag->score = (int) score;
				alignFrag->start = start;
				alignFrag->end = end;
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
			}
			
			fragCount++;
			if(fragCount >= maxFrag) {
				sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
				template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
				if(!template_fragments[fileCount]) {
					fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
					exit(-1);
				}
				outputfilename[file_len] = 0;
				fileCount++;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
				}
			}
		}
	}
	sprintf(outputfilename, "%s%s%d", outputfilename, ".tmp_", fileCount);
	template_fragments[fileCount] = printFrags(outputfilename, alignFrags);
	if(!template_fragments[fileCount]) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	outputfilename[file_len] = 0;
	fileCount++;
	fragCount = 0;
	fclose(frag_in_raw);
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(bestTemplates);
	/* remove frag_raw */
	/*strcat(outputfilename, ".frag_raw.b");
	remove(outputfilename);
	outputfilename[file_len] = 0;*/
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	
	/* Get expected values */
	template_tot_ulen = 0;
	Nhits = 0;
	for(i = 0; i < DB_size; i++) {
		Nhits += w_scores[i];
		template_tot_ulen += template_lengths[i];
	}
	
	/* Do local assemblies of fragments mapping to the same template */
	aligned = malloc(sizeof(struct aln));
	gap_align = malloc(sizeof(struct aln));
	if(!aligned || !gap_align) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	aligned_assem = malloc(sizeof(struct assem));
	aligned->t = malloc((delta + 1) << 1);
	aligned->s = malloc((delta + 1) << 1);
	aligned->q = malloc((delta + 1) << 1);
	gap_align->t = malloc((delta + 1) << 1);
	gap_align->s = malloc((delta + 1) << 1);
	gap_align->q = malloc((delta + 1) << 1);
	if(!aligned_assem || !aligned->t || !aligned->s || !aligned->q || !gap_align->t || !gap_align->s || !gap_align->q) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	aligned_assem->size = 0;
	aligned_assem->t = NULL;
	aligned_assem->s = NULL;
	aligned_assem->q = NULL;
	for(template = 1; template < DB_size; template++) {
		if(w_scores[template] > 0) {
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = (Nhits - read_score) * (1.0 * t_len) / (template_tot_ulen - t_len + etta);
			q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
			p_value  = p_chisqr(q_value);
			
			if((p_value <= evalue && score > expected) || ((1.0 * read_score / t_len) > scoreT)) {
				/* load DB */
				templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], 0, 0);
				
				/* Do assembly */
				assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, matrix_out, outputfilename, aligned, gap_align, qseq, header);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = strlen(aligned_assem->q) - countChar(aligned_assem->q, '-');
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
				}
				if(id >= ID_t) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8u\t%8d\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_names[template], read_score, (int) expected, t_len, id, cover, q_id, q_cover, depth, q_value, p_value);
					
					/* print alignment */
					aln_len = strlen(aligned_assem->s);
					fprintf(alignment_out, "# %s\n", template_names[template]);
					for(i = 0; i < aln_len; i += 60) {
						fprintf(alignment_out, "%-10s\t", "template:");
						for(j = i; j < i + 60 && aligned_assem->t[j]; j++) {
							fprintf(alignment_out, "%c", aligned_assem->t[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "");
						for(j = i; j < i + 60 && aligned_assem->s[j]; j++) {
							fprintf(alignment_out, "%c", aligned_assem->s[j]);
						}
						fprintf(alignment_out, "\n%-10s\t", "query:");
						for(j = i; j < i + 60 && aligned_assem->q[j]; j++) {
							fprintf(alignment_out, "%c", aligned_assem->q[j]);
						}
						fprintf(alignment_out, "\n\n");
					}
					
					/* Print consensus */
					fprintf(consensus_out, ">%s\n", template_names[template]);
					if(ref_fsa) {
						for(i = 0; i < aln_len; i += 60) {
							for(j = i; j < i + 60 && aligned_assem->q[j]; j++) {
								if(aligned_assem->q[j]) {
									if(aligned_assem->q[j] != '-') {
										fprintf(consensus_out, "%c", aligned_assem->q[j]);
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
							for(j = i; j < i + 60 + bias && aligned_assem->q[j]; j++) {
								if(aligned_assem->q[j] && aligned_assem->q[j] != '-') {
									fprintf(consensus_out, "%c", aligned_assem->q[j]);
								} else {
									bias++;
								}
							}
							fprintf(consensus_out, "\n");
						}
					}
				}
				/* destroy this DB index */
				destroyPtr(template);
			} else {
				fseek(index_in, (template_lengths[template] << 1) * sizeof(int), SEEK_CUR);
				fseek(seq_in, ((template_lengths[template] >> 5) + 1) * sizeof(long unsigned), SEEK_CUR);
			}
		} else {
			fseek(index_in, (template_lengths[template] << 1) * sizeof(int), SEEK_CUR);
			fseek(seq_in, ((template_lengths[template] >> 5) + 1) * sizeof(long unsigned), SEEK_CUR);
		}	
	}
	free(qseq);
	free(header);
	if(aligned_assem->size != 0) {
		free(aligned_assem->t);
		free(aligned_assem->s);
		free(aligned_assem->q);
	}
	free(aligned_assem);
	
	/* Close files */
	fclose(index_in);
	fclose(seq_in);
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
		outputfilename[file_len] = 0;
	}
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# KMA-2.0 mapps raw reads to a template database, for optimal performance it is designed to use 3 threads.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-i\t\tInput/query file name\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t16\n");
	fprintf(helpOut, "#\t-e\t\tevalue\t\t\t\t0.05\n");
	fprintf(helpOut, "#\t-mem_mode\tUse kmers to choose best\n#\t\t\ttemplate, and save memory\tFalse\n");
	fprintf(helpOut, "#\t-ex_mode\tSearh kmers exhaustively\tFalse\n");
	fprintf(helpOut, "#\t-deCon\t\tRemove contamination\t\tFalse\n");
	fprintf(helpOut, "#\t-dense\t\tDo not allow insertions\n#\t\t\tin assembly\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ref_fsa\tConsensus sequnce will\n#\t\t\thave \"n\" instead of gaps\tFalse\n");
	fprintf(helpOut, "#\t-matrix\t\tPrint assembly matrix\t\tFalse\n");
	fprintf(helpOut, "#\t-a\t\tPrint all best mappings\t\tFalse\n");
	fprintf(helpOut, "#\t-mp\t\tMinimum phred score\t\t30\n");
	fprintf(helpOut, "#\t-5p\t\tCut a constant number of\n#\t\t\tnucleotides from the 5 prime.\t0\n");
	fprintf(helpOut, "#\t-Sparse\t\tRun KmerFinder\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ID\t\tMinimum ID\t\t\t1.0%%\n");
	fprintf(helpOut, "#\t-ss\t\tSparse sorting (q,c,d)\t\tq\n");
	fprintf(helpOut, "#\t-shm\t\tUse shared DB made by kma_shm\t0 (lvl)\n");
	fprintf(helpOut, "#\t-swap\t\tSwap DB to disk\t\t\t0 (lvl)\n");
	fprintf(helpOut, "#\t-1t1\t\tSkip HMM\t\t\tFalse\n");
	fprintf(helpOut, "#\t-boot\t\tBootstrap sequence\t\tFalse\n");
	fprintf(helpOut, "#\t-mrs\t\tMinimum alignment score,\n#\t\t\tnormalized to alignment length\t0.5\n");
	fprintf(helpOut, "#\t-reward\t\tScore for match\t\t\t1\n");
	fprintf(helpOut, "#\t-penalty\tPenalty for mismatch\t\t-2\n");
	fprintf(helpOut, "#\t-gapopen\tPenalty for gap opening\t\t-3\n");
	fprintf(helpOut, "#\t-gapextend\tPenalty for gap extension\t-1\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int i, j, args, stop, exe_len, minPhred, fiveClip, sparse_run, fileCounter;
	int mem_mode, one2one, step1, step2;
	char *exeBasic, **inputfiles, *outputfilename, *templatefilename, ss;
	FILE *templatefile;
	time_t t0, t1;
	
	if(sizeof(long unsigned) != 8) {
		fprintf(stderr, "Need a 64-bit system.\n");
		exit(-1);
	}
	
	/* SET DEFAULTS */
	assemblyPtr = &assemble_KMA;
	minPhred = 30;
	fiveClip = 0;
	sparse_run = 0;
	fileCounter = 1;
	outputfilename = 0;
	templatefilename = 0;
	print_matrix = 0;
	print_all = 0;
	ref_fsa = 0;
	contamination = -1;
	kmersize = 16;
	evalue = 0.05;
	delta = 2048;
	exhaustive = 0;
	shm = 0;
	diskDB = 0;
	scoreT = 0.5;
	step1 = 0;
	step2 = 0;
	ID_t = 1.0;
	one2one = 0;
	ss = 'q';
	deCon = 0;
	mem_mode = 0;
	M = 1;
	MM = -2;
	W1 = -3;
	U = -1;
	kmerScan = &save_kmers_HMM;
	printPtr = &print_ankers;
	deConPrintPtr = printPtr;
	ankerPtr = &ankerAndClean;
	alignLoadPtr = &alignLoad_fly;
	destroyPtr = &alignClean;
	printFsa_ptr = &printFsa;
	inputfiles = malloc(sizeof(char*));
	if(!inputfiles) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	inputfiles[0] = strdup("--");
	if(!inputfiles[0]) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			stop = 0;
			args++;
			if(args < argc && (strncmp(argv[args], "-", 1) != 0 || strncmp(argv[args], "--", 2) == 0)) {
				free(inputfiles[0]);
				inputfiles[0] = strdup(argv[args]);
				if(!inputfiles[0]) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				args++;
			}
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strncmp(argv[args], "--", 2) == 0) {
					fileCounter++;
					inputfiles = realloc(inputfiles, fileCounter * sizeof(char*));
					if(!inputfiles) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					inputfiles[fileCounter - 1] = strdup(argv[args]);
					if(!inputfiles[fileCounter - 1]) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
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
				outputfilename = malloc(strlen(argv[args]) + 64);
				if(!outputfilename) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				strcpy(outputfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-deCon") == 0) {
			deCon = 1;
			deConPrintPtr = &deConPrint;
		} else if(strcmp(argv[args], "-shm") == 0) {
			args++;
			if(args < argc && argv[args][0] != '-') {
				shm = atoi(argv[args]);
			} else {
				args--;
				shm = 1;
			}
		} else if(strcmp(argv[args], "-swap") == 0) {
			args++;
			if(args < argc && argv[args][0] != '-') {
				diskDB = atoi(argv[args]);
			} else {
				args--;
				diskDB = 1;
			}
		} else if(strcmp(argv[args], "-step1") == 0) {
			step1 = 1;
		} else if(strcmp(argv[args], "-step2") == 0) {
			step2 = 1;
		} else if(strcmp(argv[args], "-mem_mode") == 0) {
			mem_mode = 1;
			alignLoadPtr = &alignLoad_fly_mem;
			ankerPtr = &ankerAndClean_MEM;
		} else if(strcmp(argv[args], "-ex_mode") == 0) {
			exhaustive = 1;
		} else if(strcmp(argv[args], "-t_db") == 0) {
			args++;
			if(args < argc) {
				templatefilename = malloc(strlen(argv[args]) + 64);
				if(!templatefilename) {
					fprintf(stderr, "OOM\n");
					exit(-1);
				}
				strcpy(templatefilename, argv[args]);
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
		} else if(strcmp(argv[args], "-a") == 0) {
			print_all = 1;
		} else if(strcmp(argv[args], "-ref_fsa") == 0) {
			ref_fsa = 1;
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
		} else if(strcmp(argv[args], "-1t1") == 0) {
			kmerScan = &save_kmers;
			one2one = 1;
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
		} else if(strcmp(argv[args], "-boot") == 0) {
			printFsa_ptr = &bootFsa;
		} else if(strcmp(argv[args], "-NW") == 0 || strcmp(argv[args], "-SW") == 0) {
			fprintf(stderr, "As of version 2.0, recursive seeding has been disabled.\n");
			fprintf(stderr, "And Needleman-Wunsch is thus always enabled.\n");
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(-1);
		}
		args++;
	}
	if(outputfilename == 0 || templatefilename == 0) {
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
		/* set to2Bit conversion */
		to2Bit = malloc(128);
		if(!to2Bit) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		for(i = 0; i < 128; i++) {
			to2Bit[i] = 5;
		}
		to2Bit['A'] = 0;
		to2Bit['C'] = 1;
		to2Bit['G'] = 2;
		to2Bit['T'] = 3;
		to2Bit['N'] = 4;
		to2Bit['a'] = 0;
		to2Bit['c'] = 1;
		to2Bit['g'] = 2;
		to2Bit['t'] = 3;
		to2Bit['n'] = 4;
		
		/* IUPAC to N */
		to2Bit['Y'] = 4;
		to2Bit['R'] = 4;
		to2Bit['W'] = 4;
		to2Bit['S'] = 4;
		to2Bit['K'] = 4;
		to2Bit['M'] = 4;
		to2Bit['D'] = 4;
		to2Bit['V'] = 4;
		to2Bit['H'] = 4;
		to2Bit['B'] = 4;
		to2Bit['X'] = 4;
		to2Bit['y'] = 4;
		to2Bit['r'] = 4;
		to2Bit['w'] = 4;
		to2Bit['s'] = 4;
		to2Bit['k'] = 4;
		to2Bit['m'] = 4;
		to2Bit['d'] = 4;
		to2Bit['v'] = 4;
		to2Bit['h'] = 4;
		to2Bit['b'] = 4;
		to2Bit['x'] = 4;
		
		/* Future RNA encoding */
		to2Bit['U'] = 3;
		to2Bit['u'] = 3;
		
		if(sparse_run) {
			templates = malloc(sizeof(struct hashMapKMA));
			if(!templates) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
			exe_len = strlen(templatefilename);
			if(deCon) {
				strcat(templatefilename, ".decon.comp.b");
			} else {
				strcat(templatefilename, ".comp.b");
			}
			templatefile = fopen(templatefilename, "rb" );
			if(!templatefile) {
				fprintf(stderr, "Wrong format of DB, or DB does not exist. 2\n");
				exit(-1);
			}
			getPrefix(templates, templatefile);
			fclose(templatefile);
			templatefilename[exe_len] = 0;
			mask = 0;
			if(templates->prefix_len) {
				mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (templates->prefix_len << 1));
			} else {
				mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
			}
			run_input_sparse(inputfiles, fileCounter, minPhred, fiveClip);
		} else {
			run_input(inputfiles, fileCounter, minPhred, fiveClip);
		}
		fflush(stdout);
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for converting query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else if(step2) {
		exe_len = strlen("-step1");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step1");
		
		save_kmers_batch(templatefilename, exeBasic, one2one);
		fflush(stdout);
	} else if(sparse_run) {
		exe_len = strlen("-step1");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step1");
		
		save_kmers_sparse_batch(templatefilename, outputfilename, exeBasic, ss);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	} else if(mem_mode) {
		exe_len = strlen("-step2");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step2");
		
		runKMA_MEM(templatefilename, outputfilename, exeBasic);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	} else {
		exe_len = strlen("-step2");
		for(args = 0; args < argc; args++) {
			exe_len += strlen(argv[args]) + 1;
		}
		exeBasic = calloc((exe_len + 1), sizeof(char));
		if(!exeBasic) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		
		for(args = 0; args < argc; args++) {
			strcat(exeBasic, argv[args]);
			strcat(exeBasic, " ");
		}
		strcat(exeBasic, "-step2");
		runKMA(templatefilename, outputfilename, exeBasic);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	}
	
	return 0;
}
