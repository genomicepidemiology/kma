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
#define SetNuc(seq, nuc, pos)(seq[pos >> 5] |= (nuc << (62-((pos & 31) << 1))))

/*
	STRUCTURES
*/
struct compDNA {
	int seqlen;
	int size;
	int complen;
	//int nlen; // == first element
	//int nsize; // == seqsize
	long unsigned *seq;
	int *N;
};

struct KMA_SEQ {
	unsigned size;
	unsigned n;
	long unsigned *seq;
};

struct hashMapKMA {
	/* end product of script */
	unsigned kmersize;		// k
	long unsigned size;		// size of DB
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

struct hashTable {
	unsigned key;
	unsigned values;
	struct hashTable *next;
};

struct hashMap {
	/* open hash structure */
	unsigned kmersize;		// k
	long unsigned size;		// size of DB
	unsigned n;				// k-mers stored
	unsigned seq_size;		// size of seq
	unsigned seq_n;			// next seq index
	unsigned prefix_len;	// prefix length
	long unsigned prefix;	// prefix
	struct hashTable **table; // outer = exist, inner = key,value index
	long unsigned *seq;		// compressed sequence of k-mers
	int **values;		// uncompressed values
};

struct hashMap_index {
	unsigned len; // seqlen
	unsigned size; // size of index
	int *index; // k-mer posititions in seq
	long unsigned *seq; // 2-bit sequence
};

struct hashTable_kmers {
	long unsigned key;
	struct hashTable_kmers *next;
};

struct hashMap_kmers {
	unsigned size;
	unsigned n;
	struct hashTable_kmers **table;
};

struct valuesTable {
	unsigned values;
	unsigned value_index;
	struct valuesTable *next;
};

struct valuesHash {
	unsigned n;
	unsigned size;
	struct valuesTable **table;
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


/*
	GLOBAL VARIABLES
*/
struct hashMap *templates;
struct hashMap_kmers *foundKmers;
int kmersize, kmerindex, DB_size, prefix_len, MinLen, MinKlen, shifter;
unsigned shifterI, prefix_shifter, *Scores, *Scores_tot, *bestTemplates;
unsigned *template_lengths, *template_slengths, *template_ulengths;
long unsigned mask, INITIAL_SIZE, prefix;
char *to2Bit;
double homQ, homT;

/* 
	FUNCTION POINTERS
*/
int (*homcmp)(int, int);
int (*update_DB)(struct compDNA *);
void (*load_ptr)(struct hashMap*, char*);
int (*deConNode_ptr)(struct compDNA *, struct hashMapKMA *);
int (*addCont)(struct hashMapKMA *, long unsigned, int);
int (*QualCheck)(struct compDNA *);
void (*updateAnnotsPtr)(struct compDNA *, FILE *, FILE *);
int (*hashMap_add)(long unsigned, int, int);
int * (*hashMap_get)(long unsigned);
long unsigned (*getKmerP)(long unsigned *, unsigned);
unsigned (*addKmer)(long unsigned, int);

/*
	FUNCTIONS
*/

void OOM() {
	fprintf(stderr, "OOM\n");
	exit(errno);
}

int chomp(char *string) {
	/* remove trailing spaces and newlines */
	int k = strlen(string) - 1;
	/* isspace = ((string[k] >= 9  && string[k] <= 13) || string[k] == 32), in ASCII */
	while (isspace(string[k]))
		k--;
	k++;
	string[k] = 0;
	return k;
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
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		} else {
			return line;
		}
	}
	line[0] = 0;
	return line;
}

int uint_eq(const unsigned *s1, const unsigned *s2, int len) {
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

int int_eq(const int *s1, const int *s2, int len) {
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

void convertToNum(char *qseq, int q_len) {
	int i;
	for(i = 0; i < q_len; i++) {
		qseq[i] = to2Bit[qseq[i]];
	}
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

/*
	COMPRESSION FUNCTIONS
*/
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

void reallocComp(struct compDNA *compressor, int size) {
	
	if(size & 31) {
		size = (size >> 5) + 1;
		size <<= 5;
	}
	
	compressor->seq = realloc(compressor->seq, (size >> 5) * sizeof(long unsigned));
	compressor->N = realloc(compressor->N, (size + 1) * sizeof(int));
	if(!compressor->seq || !compressor->N) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	} else {
		memset(compressor->seq + (compressor->size >> 5), 0, ((size >> 5) - (compressor->size >> 5)) * sizeof(long unsigned));
	}
	
	compressor->size = size;
	
}

void compDNA(struct compDNA *compressor, char *seq, int seqlen) {
	
	int i, j, pos, end;
	char nuc;
	
	compressor->seqlen = seqlen;
	if(seqlen & 31) {
		compressor->complen = (seqlen >> 5) + 1;
	} else {
		compressor->complen = seqlen >> 5;
	}
	
	pos = 0;
	compressor->N[0] = 0;
	for(i = 0; i < seqlen; i += 32) {
		end = (i + 32 < seqlen) ? i + 32 : seqlen;
		pos = i >> 5;
		for(j = i; j < end; j++) {
			nuc = seq[j];
			if(nuc == 4) {
				compressor->seq[pos] <<= 2;
				compressor->N[0]++;
				compressor->N[compressor->N[0]] = j;
			} else {
				compressor->seq[pos] = (compressor->seq[pos] << 2) | nuc;
			}
		}
	}
	if(seqlen & 31) {
		compressor->seq[pos] <<= (64 - ((seqlen & 31) << 1));
	}
}

int compDNAref(struct compDNA *compressor, char *qseq, int seqlen) {
	
	int i, j, pos, end, bias;
	char nuc, *seq;
	
	/* trim leadin N's */
	seq = qseq;
	bias = 0;
	while(*seq == 4) {
		seq++;
		bias++;
	}
	seqlen -= bias;
	/* trim trailing N's */
	seqlen--;
	while(seq[seqlen] == 4) {
		seqlen--;
	}
	seqlen++;
	
	compressor->seqlen = seqlen;
	if(seqlen & 31) {
		compressor->complen = (seqlen >> 5) + 1;
	} else {
		compressor->complen = seqlen >> 5;
	}
	
	pos = 0;
	nuc = 0;
	compressor->N[0] = 0;
	for(i = 0; i < seqlen; i += 32) {
		end = (i + 32 < seqlen) ? i + 32 : seqlen;
		pos = i >> 5;
		for(j = i; j < end; j++) {
			if(seq[j] == 4) {
				compressor->seq[pos] = (compressor->seq[pos] << 2) | nuc;
				compressor->N[0]++;
				compressor->N[compressor->N[0]] = j;
			} else {
				nuc = seq[j];
				compressor->seq[pos] = (compressor->seq[pos] << 2) | nuc;
			}
		}
	}
	if(seqlen & 31) {
		compressor->seq[pos] <<= (64 - ((seqlen & 31) << 1));
	}
	
	return bias;
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

long unsigned getKmer(long unsigned *compressor, unsigned pos) {
	
	unsigned cPos, iPos;
	cPos = pos >> 5;
	iPos = (pos & 31) << 1;
	
	return (iPos <= shifter) ? ((compressor[cPos] << iPos) >> shifter) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifter);
}

long unsigned getK(long unsigned *compressor, unsigned pos) {
	return pos;
}

unsigned addKmerC(long unsigned key, int extend) {
	
	int i, pos;
	
	pos = templates->seq_n - kmersize + extend;
	key <<= ((kmersize - extend) << 1);
	for(i = kmersize - extend; i < kmersize; i++) {
		key <<= 2;
		SetNuc(templates->seq, ((key >> (kmersize << 1)) & 3), templates->seq_n);
		templates->seq_n++;
	}
	
	return pos;
}

unsigned addK(long unsigned key, int extend) {
	return key;
}

long unsigned getKmerIndex(long unsigned *compressor, unsigned pos) {
	
	unsigned cPos, iPos;
	cPos = pos >> 5;
	iPos = (pos & 31) << 1;
	
	return (iPos <= shifterI) ? ((compressor[cPos] << iPos) >> shifterI) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifterI);
}

long unsigned getPrefix(long unsigned *compressor, unsigned pos) {
	int cPos = pos >> 5, iPos = (pos & 31) << 1;
	
	if(prefix_len == 0) {
		return 0;
	} else if(iPos <= prefix_shifter) {
		return (compressor[cPos] << iPos) >> prefix_shifter;
	} else {
		return ((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> prefix_shifter;
	}
}

void setNuc(long unsigned *seq, long unsigned nuc, unsigned pos) {
	//int cPos = pos >> 5, iPos = (pos & 31) << 1;
	//int cPos = pos / 32, iPos = (pos % 32) * 2;
	
	//seq[cPos] |= (nuc << (62-iPos));
	seq[pos >> 5] |= (nuc << (62 - ((pos & 31) << 1)));
}

long unsigned binRev(long unsigned num) {
	
	unsigned int count = 62;
	long unsigned rev = num;
	
	num >>= 2;
	while(num) {
		rev <<= 2;
		rev |= num & 3;
		num >>=2;
		count -= 2;
	}
	rev <<= count;
	
	return rev;
}

long unsigned binRev2(long unsigned mer) {
	
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
		compressor_rc->seq[j] = binRev2(~compressor->seq[i]);
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
	compressor_rc->N[0] = compressor->N[0];
	r_shift = compressor->seqlen - 1;
	shift = compressor->N[0];
	for(i = 1, j = compressor->N[0]; i <= shift; i++, j--) {
		compressor_rc->N[i] = r_shift - compressor->N[j];
		//compressor_rc->N[i] = compressor->seqlen - compressor->N[i] - 1;
	}
	
}

void rcComp(struct compDNA *compressor) {
	
	int i, j, shift, r_shift, complen;
	long unsigned carry;
	
	/* reverse and complement*/
	complen = compressor->complen >> 1;
	for(i = 0, j = compressor->complen - 1; i < complen; i++, j--) {
		carry = binRev2(~compressor->seq[i]);
		compressor->seq[i] = binRev2(~compressor->seq[j]);
		compressor->seq[j] = carry;
	}
	if((compressor->complen & 1)) {
		compressor->seq[complen] = binRev2(~compressor->seq[complen]);
	}
	
	/* shift */
	if((compressor->seqlen & 31)) {
		shift = (((compressor->complen << 5) - compressor->seqlen) << 1);
		r_shift = 64 - shift;
		for(i = 0, j = 1; j < compressor->complen; i++, j++) {
			compressor->seq[i] = (compressor->seq[i] << shift) | (compressor->seq[j] >> r_shift);
		}
		compressor->seq[i] <<= shift;
	}
	
	/* Change N's */
	j = compressor->N[0] - 1;
	complen = compressor->N[0] >> 1;
	compressor->N++;
	r_shift = compressor->seqlen - 1;
	for(i = 0; i < complen; i++, j--) {
		shift = r_shift - compressor->N[i];
		compressor->N[i] = r_shift - compressor->N[j];
		compressor->N[j] = shift;
	}
	compressor->N--;
	if((compressor->N[0] & 1)) {
		complen++;
		compressor->N[complen] = r_shift - compressor->N[complen];
	}
	
}

void freeComp(struct compDNA *compressor) {
	
	compressor->seqlen = 0;
	compressor->complen = 0;
	compressor->size = 0;
	
	free(compressor->seq);
	free(compressor->N);
	
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
		if(seqsize + header->len > header->size) {
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
	if(seqsize + header->len > header->size) {
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
				dest->pos = destpos;
				seq->len = seqlen;
				return 1;
			}
			destbytes = dest->bytes;
		}
		
	}
	dest->pos = destpos;
	seq->len = seqlen;
	
	return 1;
}

/*
	HASHMAP FUNCTIONS
*/

struct hashMap * hashMap_initialize(long unsigned size) {
	
	struct hashMap *src;
	
	src = malloc(sizeof(struct hashMap));
	if(!src) {
		OOM();
	}
	
	src->kmersize = kmersize;
	src->size = size;
	src->n = 0;
	src->seq_size = size >> 5;
	src->seq_n = 0;
	src->prefix_len = 0;
	src->prefix = 0;
	
	if((size - 1) == mask) {
		src->table = 0;
		src->seq_size = 0;
		src->seq = 0;
		src->values = calloc(src->size, sizeof(unsigned *));
		if(!src->values) {
			OOM();
		}
	} else {
		src->table = calloc(src->size, sizeof(struct hashTable *));
		src->seq = calloc(src->seq_size, sizeof(long unsigned));
		src->values = malloc(src->size * sizeof(unsigned *));
		if(!src->table || !src->seq || !src->values) {
			OOM();
		}
	}
	
	/* masking */
	src->size--;
	
	return src;
}

struct hashMap * hashMap_load(char *filename) {
	
	unsigned i, index, pos[2];
	FILE *infile;
	struct hashMap *src;
	struct hashTable *node;
	
	infile = fopen(filename, "rb");
	if(!infile) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(errno);
	}
	
	src = malloc(sizeof(struct hashMap));
	if(!src) {
		OOM();
	}
	
	/* load content */
	fread(&src->kmersize, sizeof(unsigned), 1, infile);
	kmersize = src->kmersize;
	if(kmersize <= 16) {
		getKmerP = &getK;
		addKmer = &addK;
	}
	fread(&src->size, sizeof(long unsigned), 1, infile);
	fread(&src->n, sizeof(unsigned), 1, infile);
	fread(&src->seq_n, sizeof(unsigned), 1, infile);
	src->seq_size = src->seq_n >> 4;
	fread(&src->prefix_len, sizeof(unsigned), 1, infile);
	fread(&src->prefix, sizeof(long unsigned), 1, infile);
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	
	/* allocate */
	if((src->size - 1) == mask) {
		src->values = calloc(src->size, sizeof(unsigned *));
		if(!src->values) {
			OOM();
		}
		/* masking */
		src->size--;
		
		/* load */
		for(i = 0; i < src->n; i++) {
			fread(pos, sizeof(unsigned), 2, infile);
			src->values[pos[0]] = malloc((pos[1] + 1) * sizeof(unsigned));
			if(!src->values[pos[0]]) {
				OOM();
			}
			src->values[pos[0]][0] = pos[1];
			fread(src->values[pos[0]] + 1, sizeof(unsigned), pos[1], infile);
		}
	} else {
		src->table = calloc(src->size, sizeof(struct hashTable *));
		src->seq = calloc(src->seq_size, sizeof(long unsigned));
		src->values = malloc(src->size * sizeof(unsigned *));
		if(!src->table || !src->seq || !src->values) {
			OOM();
		}
		/* masking */
		src->size--;
		
		/* load */
		fread(src->seq, sizeof(long unsigned), (src->seq_n >> 5) + 1, infile);
		for(i = 0; i < src->n; i++) {
			fread(pos, sizeof(unsigned), 1, infile);
			src->values[i] = malloc((pos[0] + 1) * sizeof(unsigned));
			if(!src->values[i]) {
				OOM();
			}
			src->values[i][0] = pos[0];
			fread(src->values[i] + 1, sizeof(unsigned), pos[0], infile);
		}
		
		for(i = 0; i < src->n; i++) {
			//fread(pos, sizeof(unsigned), 2, infile);
			node = malloc(sizeof(struct hashTable));
			if(!node) {
				OOM();
			}
			/* push node */
			//fread(&index, sizeof(unsigned), 1, infile);
			fread(&node->key, sizeof(unsigned), 1, infile);
			fread(&node->values, sizeof(int), 1, infile);
			//node->key = pos[0];
			//node->values = pos[1];
			index = getKmerP(src->seq, node->key) & src->size;
			node->next = src->table[index];
			src->table[index] = node;
		}
	}
	
	fclose(infile);
	return src;
}

void hashMap_dump(struct hashMap *src, FILE *outfile) {
	
	long unsigned i;
	struct hashTable *node;
	
	/* dump content */
	src->size++;
	fwrite(&src->kmersize, sizeof(unsigned), 1, outfile);
	fwrite(&src->size, sizeof(long unsigned), 1, outfile);
	fwrite(&src->n, sizeof(unsigned), 1, outfile);
	fwrite(&src->seq_n, sizeof(unsigned), 1, outfile);
	fwrite(&src->prefix_len, sizeof(unsigned), 1, outfile);
	fwrite(&src->prefix, sizeof(long unsigned), 1, outfile);
	
	/* allocate */
	if((src->size - 1) == mask) {
		/* dump values */
		for(i = 0; i < src->size; i++) {
			if(src->values[i]) {
				fwrite(&(unsigned){i}, sizeof(unsigned), 1, outfile);
				fwrite(src->values[i], sizeof(unsigned), src->values[i][0] + 1, outfile);
			}
		}
	} else {
		/* dump */
		fwrite(src->seq, sizeof(long unsigned), (src->seq_n >> 5) + 1, outfile);
		for(i = 0; i < src->n; i++) {
			fwrite(src->values[i], sizeof(unsigned), src->values[i][0] + 1, outfile);
		}
		for(i = 0; i < src->size; i++) {
			for(node = src->table[i]; node != 0; node = node->next) {
				//fwrite(&(unsigned){i}, sizeof(unsigned), 1, outfile);
				fwrite(&node->key, sizeof(unsigned), 1, outfile);
				fwrite(&node->values, sizeof(int), 1, outfile);
			}
		}
	}
	
	/* masking */
	
	src->size--;
}

void deConMap_dump(struct hashMapKMA *dest, FILE *out) {
	
	long unsigned i;
	
	/* dump sizes */
	fwrite(&DB_size, sizeof(unsigned), 1, out);
	fwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	fwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	fwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	fwrite(&dest->size, sizeof(long unsigned), 1, out);
	fwrite(&dest->n, sizeof(unsigned), 1, out);
	fwrite(&dest->seqsize, sizeof(unsigned), 1, out); //seq size
	fwrite(&dest->v_index, sizeof(unsigned), 1, out);
	fwrite(&dest->null_index, sizeof(unsigned), 1, out);
	
	/* dump arrays */
	fwrite(dest->exist, sizeof(unsigned), dest->size, out);
	
	if(dest->seqsize != 0) {
		fwrite(dest->seq, sizeof(long unsigned), dest->seqsize, out);
		fwrite(dest->key_index, sizeof(unsigned), dest->n + 1, out);
		fwrite(dest->value_index, sizeof(unsigned), dest->n, out);
		for(i = 0; i < dest->n; i++) {
			fwrite(templates->values[i], sizeof(unsigned), templates->values[i][0] + 1, out);
		}
	} else {
		for(i = 0; i < dest->size; i++) {
			if(templates->values[i]) {
				fwrite(templates->values[i], sizeof(unsigned), templates->values[i][0] + 1, out);
			}
		}
	}
}

int CP(char *templatefilename, char *outputfilename) {
	
	int bytes, buffSize;
	char *buffer;
	FILE *file_in, *file_out;
	
	if(strcmp(templatefilename, outputfilename) == 0) {
		return 1;
	}
	
	file_in = fopen(templatefilename, "rb");
	file_out = fopen(outputfilename, "wb");
	if(!file_in || !file_out) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
	buffSize = 1024 * 1024;
	buffer = malloc(buffSize);
	if(!buffer) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	while((bytes = fread(buffer, 1, buffSize, file_in))) {
		fwrite(buffer, 1, bytes, file_out);
	}
	
	fclose(file_in);
	fclose(file_out);
	
	free(buffer);
	
	return 0;
}

int load_DBs(char *templatefilename, char *outputfilename) {
	
	int file_len, out_len, appender;
	FILE *infile;
	
	file_len = strlen(templatefilename);
	out_len = strlen(outputfilename);
	
	/* load hash */
	strcat(templatefilename, ".b");
	templates = hashMap_load(templatefilename);
	templatefilename[file_len] = 0;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	
	
	/* load lengths */
	strcat(templatefilename, ".length.b");
	infile = fopen(templatefilename, "rb");
	if(!infile) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	templatefilename[file_len] = 0;
	fread(&DB_size, sizeof(unsigned), 1, infile);
	
	if(prefix_len) {
		template_lengths = malloc((DB_size << 1) * sizeof(unsigned));
		template_slengths = malloc((DB_size << 1) * sizeof(unsigned));
		template_ulengths = malloc((DB_size << 1) * sizeof(unsigned));
		if(!template_lengths || !template_slengths || !template_ulengths) {
			OOM();
		}
		fread(template_slengths, sizeof(unsigned), DB_size, infile);
		fread(template_ulengths, sizeof(unsigned), DB_size, infile);
		fread(template_lengths, sizeof(unsigned), DB_size, infile);
		kmerindex = *template_lengths;
		template_ulengths[0] = DB_size << 1;
		template_slengths[0] = DB_size << 1;
	} else {
		template_lengths = malloc((DB_size << 1) * sizeof(unsigned));
		template_slengths = 0;
		template_ulengths = 0;
		if(!template_lengths) {
			OOM();
		}
		fread(template_lengths, sizeof(unsigned), DB_size, infile);
		kmerindex = *template_lengths;
		template_lengths[0] = DB_size << 1;
	}
	
	fclose(infile);
	
	/* cp name, seq and index */
	strcat(templatefilename, ".name");
	strcat(outputfilename, ".name");
	appender = CP(templatefilename, outputfilename);
	templatefilename[file_len] = 0;
	outputfilename[out_len] = 0;
	
	strcat(templatefilename, ".seq.b");
	strcat(outputfilename, ".seq.b");
	appender = CP(templatefilename, outputfilename);
	templatefilename[file_len] = 0;
	outputfilename[out_len] = 0;
	
	if(template_ulengths == 0) {
		strcat(templatefilename, ".index.b");
		strcat(outputfilename, ".index.b");
		appender = CP(templatefilename, outputfilename);
		templatefilename[file_len] = 0;
		outputfilename[out_len] = 0;
	}
	
	return 1;
	//return appender;
}

int megaMap_addKMA(long unsigned key, int value, int extend) {
	
	int *values;
	
	values = templates->values[key];
	
	if(values == 0) { 
		values = malloc(2 * sizeof(int));
		values[0] = 1;
		values[1] = value;
		templates->values[key] = values;
		templates->n++;
	} else if(values[values[0]] != value) {
		values[0]++;
		values = realloc(values, (values[0] + 1) * sizeof(int));
		if(!values) {
			OOM();
		}
		values[values[0]] = value;
		templates->values[key] = values;
	}
	
	return 0;
}

int megaMap_addKMA_sparse(long unsigned key, int value, int extend) {
	
	int *values;
	
	values = templates->values[key];
	
	if(values == 0) { 
		values = malloc(2 * sizeof(int));
		values[0] = 1;
		values[1] = value;
		templates->values[key] = values;
		templates->n++;
		template_ulengths[value]++;
	} else if(values[values[0]] != value) {
		values[0]++;
		values = realloc(values, (values[0] + 1) * sizeof(int));
		if(!values) {
			OOM();
		}
		values[values[0]] = value;
		templates->values[key] = values;
		template_ulengths[value]++;
	}
	
	return 0;
}

int * megaMap_getValue(long unsigned key) {
	
	return templates->values[key];
	
}

int hashMap_addCont(struct hashMapKMA *dest, long unsigned key, int value) {
	
	int *values;
	unsigned pos, kpos;
	long unsigned kmer;
	
	kpos = key & dest->size;
	pos = dest->exist[kpos];
	
	if(pos != dest->null_index) {
		kmer = getKmerP(dest->seq, dest->key_index[pos]);
		while(key != kmer) {
			pos++;
			if(kpos != (kmer & dest->size)) {
				return 0;
			}
			kmer = getKmerP(dest->seq, dest->key_index[pos]);
		}
		values = templates->values[dest->value_index[pos]];
		if(values[values[0]] != value) {
			values[0]++;
			values = realloc(values, (values[0] + 1) * sizeof(unsigned));
			if(!values) {
				OOM();
			}
			values[values[0]] = value;
			templates->values[dest->value_index[pos]] = values;
			return 1;
		}
	}
	
	return 0;
}

int megaMap_addCont(struct hashMapKMA *dest, long unsigned index, int value) {
	
	int *values;
	
	if(dest->exist[index] != dest->n) {
		values = templates->values[dest->exist[index]];
		if(values[values[0]] != value) {
			values[0]++;
			values = realloc(values, (values[0] + 1) * sizeof(unsigned));
			if(!values) {
				OOM();
			}
			values[values[0]] = value;
			templates->values[dest->exist[index]] = values;
			return 1;
		}
	}
	
	return 0;
}

void hashMap2megaMap(struct hashTable *table) {
	
	int **tmp;
	struct hashTable *node, *next;
	
	tmp = templates->values;
	templates->values = calloc(templates->size, sizeof(unsigned *));
	if(!templates->values) {
		OOM();
	}
	templates->size--;
	
	/* table sorted in semi descending order */
	for(node = table; node != 0; node = next) {
		next = node->next;
		
		/* move values */
		templates->values[getKmerP(templates->seq, node->key)] = tmp[node->values];
		
		/* clean */
		free(node);
	}
	
	/* clean */
	free(tmp);
	free(templates->seq);
	templates->seq = 0;
	templates->table = 0;
	
	/* set pointers */
	hashMap_add = &megaMap_addKMA;
	hashMap_get = &megaMap_getValue;
	addCont = &megaMap_addCont;
}

int hashMap_addKMA(long unsigned key, int value, int extend) {
	
	int *values;
	unsigned i, index;
	struct hashTable *node, *next, *table;
	
	index = key & templates->size;
	/* check if key exists */
	for(node = templates->table[index]; node != 0; node = node->next) {
		if(key == getKmerP(templates->seq, node->key)) {
			values = templates->values[node->values];
			if(values[*values] != value) {
				values[0]++;
				values = realloc(values, (values[0] + 1) * sizeof(unsigned));
				if(!values) {
					OOM();
				}
				values[*values] = value;
				templates->values[node->values] = values;
			}
			return 0;
		}
	}
	
	/* new value check if there is space */
	if(templates->n == templates->size) {
		templates->size++;
		/* link table */
		table = 0;
		for(i = 0; i < templates->size; i++) {
			for(node = templates->table[i]; node != 0; node = next) {
				next = node->next;
				node->next = table;
				table = node;
			}
		}
		free(templates->table);
		
		/* check for megamap */
		templates->size <<= 1;
		if((templates->size - 1) == mask) {
			hashMap2megaMap(table);
			return megaMap_addKMA(key, value, extend);
		}
		
		/* reallocate */
		templates->values = realloc(templates->values, templates->size * sizeof(unsigned *));
		templates->table = calloc(templates->size, sizeof(struct hashTable));
		if(!templates->values || !templates->table) {
			OOM();
		}
		templates->size--;
		
		for(node = table; node != 0; node = next) {
			next = node->next;
			i = getKmerP(templates->seq, node->key) & templates->size;
			node->next = templates->table[i];
			templates->table[i] = node;
		}
	}
	if(((templates->seq_n + kmersize) >> 5) >= (templates->seq_size - 1)) {
		/* resize hashMapSEQ */
		i = templates->seq_size;
		templates->seq_size <<= 1;
		templates->seq = realloc(templates->seq, templates->seq_size * sizeof(long unsigned));
		
		if(!templates->seq) {
			OOM();
		}
		/* nullify new chunk */
		for(; i < templates->seq_size; i++) {
			templates->seq[i] = 0;
		}
	}
	/* add new value */
	node = malloc(sizeof(struct hashTable));
	if(!node) {
		OOM();
	}
	/* key */
	node->key = addKmer(key, extend);
	/*if(extend) {
		node->key = templates->seq_n - kmersize + 1;
		SetNuc(templates->seq, (key & 3), templates->seq_n);
		templates->seq_n++;
	} else {
		node->key = templates->seq_n;
		for(i = 0; i < kmersize; i++) {
			key <<= 2;
			SetNuc(templates->seq, ((key >> (kmersize << 1)) & 3), templates->seq_n);
			templates->seq_n++;
		}
	}*/
	
	/* value */
	node->values = templates->n;
	templates->values[node->values] = malloc(2 * sizeof(unsigned));
	templates->values[node->values][0] = 1;
	templates->values[node->values][1] = value;
	
	/* push it */
	node->next = templates->table[index];
	templates->table[index] = node;
	
	templates->n++;
	
	return 1;
}

int hashMap_addKMASparse(long unsigned key, int value, int extend) {
	
	int *values;
	unsigned i, index;
	struct hashTable *node, *next, *table;
	
	index = key & templates->size;
	
	/* check if key exists */
	for(node = templates->table[index]; node != 0; node = node->next) {
		if(key == getKmerP(templates->seq, node->key)) {
			values = templates->values[node->values];
			if(values[*values] != value) {
				values[0]++;
				values = realloc(values, (values[0] + 1) * sizeof(unsigned));
				if(!values) {
					OOM();
				}
				values[*values] = value;
				templates->values[node->values] = values;
				template_ulengths[value]++;
			}
			return 0;
		}
	}
	
	/* new value check if there is space */
	if(templates->n == templates->size) {
		templates->size++;
		/* link table */
		table = 0;
		for(i = 0; i < templates->size; i++) {
			for(node = templates->table[i]; node != 0; node = next) {
				next = node->next;
				node->next = table;
				table = node;
			}
		}
		free(templates->table);
		
		/* check for megamap */
		templates->size <<= 1;
		if((templates->size - 1) == mask) {
			hashMap2megaMap(table);
			hashMap_add = &megaMap_addKMA_sparse;
			return megaMap_addKMA_sparse(key, value, extend);
		}
		
		/* reallocate */
		templates->values = realloc(templates->values, templates->size * sizeof(unsigned *));
		templates->table = calloc(templates->size, sizeof(struct hashTable));
		if(!templates->values || !templates->table) {
			OOM();
		}
		templates->size--;
		
		for(node = table; node != 0; node = next) {
			next = node->next;
			i = getKmerP(templates->seq, node->key) & templates->size;
			node->next = templates->table[i];
			templates->table[i] = node;
		}
	}
	if(((templates->seq_n + kmersize) >> 5) >= (templates->seq_size - 1)) {
		/* resize hashMapSEQ */
		i = templates->seq_size;
		templates->seq_size <<= 1;
		templates->seq = realloc(templates->seq, templates->seq_size * sizeof(long unsigned));
		
		if(!templates->seq) {
			OOM();
		}
		/* nullify new chunk */
		for(; i < templates->seq_size; i++) {
			templates->seq[i] = 0;
		}
	}
	
	/* add new value */
	node = malloc(sizeof(struct hashTable));
	if(!node) {
		OOM();
	}
	/* key */
	node->key = addKmer(key, extend);
	/*node->key = templates->seq_n - kmersize + extend;
	key <<= ((kmersize - extend) << 1);
	for(i = kmersize - extend; i < kmersize; i++) {
		key <<= 2;
		SetNuc(templates->seq, ((key >> (kmersize << 1)) & 3), templates->seq_n);
		templates->seq_n++;
	}*/
	
	/* value */
	node->values = templates->n;
	templates->values[node->values] = malloc(2 * sizeof(unsigned));
	templates->values[node->values][0] = 1;
	templates->values[node->values][1] = value;
	
	/* push it */
	node->next = templates->table[index];
	templates->table[index] = node;
	
	templates->n++;
	template_ulengths[value]++;
	
	return 1;
}

int * hashMap_getValue(long unsigned key) {
	
	struct hashTable *node;
	
	for(node = templates->table[key & templates->size]; node != 0; node = node->next) {
		if(key == getKmerP(templates->seq, node->key)) {
			return templates->values[node->values];
		}
	}
	
	return 0;
}

/*
	COMPRESSED HASHMAP FUNCTIONS
*/
void hashMapKMA_dump(struct hashMapKMA *dest, FILE *out) {
	
	/* dump sizes */
	fwrite(&DB_size, sizeof(unsigned), 1, out);
	fwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	fwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	fwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	fwrite(&dest->size, sizeof(long unsigned), 1, out);
	fwrite(&dest->n, sizeof(unsigned), 1, out);
	fwrite(&dest->seqsize, sizeof(unsigned), 1, out); //seq size
	fwrite(&dest->v_index, sizeof(unsigned), 1, out);
	fwrite(&dest->null_index, sizeof(unsigned), 1, out);
	
	/* dump arrays */
	fwrite(dest->exist, sizeof(unsigned), dest->size, out);
	fwrite(dest->seq, sizeof(long unsigned), dest->seqsize, out);
	fwrite(dest->values, sizeof(int), dest->v_index, out);
	fwrite(dest->key_index, sizeof(unsigned), dest->n + 1, out);
	fwrite(dest->value_index, sizeof(unsigned), dest->n, out);
	
}

void megaMapKMA_dump(struct hashMapKMA *dest, FILE *out) {
	
	/* dump sizes */
	fwrite(&DB_size, sizeof(unsigned), 1, out);
	fwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	fwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	fwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	fwrite(&dest->size, sizeof(long unsigned), 1, out);
	fwrite(&dest->n, sizeof(unsigned), 1, out);
	fwrite(&dest->seqsize, sizeof(unsigned), 1, out); //seq size
	fwrite(&dest->v_index, sizeof(unsigned), 1, out);
	fwrite(&dest->null_index, sizeof(unsigned), 1, out);
	
	/* dump arrays */
	fwrite(dest->exist, sizeof(unsigned), dest->size, out);
	fwrite(dest->values, sizeof(int), dest->v_index, out);
	
}

void hashMap_shm_detach(struct hashMapKMA *dest) {
	shmdt(dest->exist);
	shmdt(dest->seq);
	shmdt(dest->values);
	shmdt(dest->key_index);
	shmdt(dest->value_index);
}

void hashMapKMA_convertSHM(struct hashMapKMA *dest, const char *filename) {
	
	int shmid;
	void *tmp;
	key_t key;
	
	/* check shared memory, else load */
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap e\n");
	} else {
		/* found */
		tmp = shmat(shmid, NULL, 0);
		memcpy(tmp, dest->exist, dest->size * sizeof(unsigned));
		shmdt(tmp);
	}
	key = ftok(filename, 's');
	shmid = shmget(key, dest->seqsize * sizeof(long unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap s\n");
	} else {
		/* found */
		tmp = shmat(shmid, NULL, 0);
		memcpy(tmp, dest->seq, dest->seqsize * sizeof(long unsigned));
		shmdt(tmp);
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, dest->v_index * sizeof(int), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
	} else {
		/* found */
		tmp = shmat(shmid, NULL, 0);
		memcpy(tmp, dest->values, dest->v_index * sizeof(int));
		shmdt(tmp);
	}
	key = ftok(filename, 'k');
	shmid = shmget(key, (dest->n + 1) * sizeof(unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap k\n");
	} else {
		/* found */
		tmp = shmat(shmid, NULL, 0);
		memcpy(tmp, dest->key_index, (dest->n + 1) * sizeof(unsigned));
		shmdt(tmp);
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, dest->n * sizeof(unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap i\n");
	} else {
		/* found */
		tmp = shmat(shmid, NULL, 0);
		memcpy(tmp, dest->value_index, dest->n * sizeof(unsigned));
		shmdt(tmp);
	}
	
}

void hashMapKMA_setupSHM(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid;
	key_t key;
	
	/* load sizes */
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	/* check shared memory, else load */
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap e\n");
		fseek(file, dest->size * sizeof(unsigned), SEEK_CUR);
		dest->exist = 0;
	} else {
		dest->exist = shmat(shmid, NULL, 0);
		fread(dest->exist, sizeof(unsigned), dest->size, file);
	}
	key = ftok(filename, 's');
	shmid = shmget(key, dest->seqsize * sizeof(long unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap s\n");
		fseek(file, dest->seqsize * sizeof(long unsigned), SEEK_CUR);
		dest->seq = 0;
	} else {
		/* found */
		dest->seq = shmat(shmid, NULL, 0);
		fread(dest->seq, sizeof(long unsigned), dest->seqsize, file);
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, dest->v_index * sizeof(int), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
		fseek(file, dest->v_index * sizeof(int), SEEK_CUR);
		dest->values = 0;
	} else {
		/* found */
		dest->values = shmat(shmid, NULL, 0);
		fread(dest->values, sizeof(unsigned), dest->v_index, file);
	}
	key = ftok(filename, 'k');
	shmid = shmget(key, (dest->n + 1) * sizeof(unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap k\n");
		fseek(file, (dest->n + 1) * sizeof(unsigned), SEEK_CUR);
		dest->key_index = 0;
	} else {
		/* found */
		dest->key_index = shmat(shmid, NULL, 0);
		fread(dest->key_index, sizeof(unsigned), dest->n + 1, file);
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, dest->n * sizeof(unsigned), IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap i\n");
		fseek(file, dest->n * sizeof(unsigned), SEEK_CUR);
		dest->value_index = 0;
	} else {
		/* found */
		dest->value_index = shmat(shmid, NULL, 0);
		fread(dest->value_index, sizeof(unsigned), dest->n, file);
	}
}

void hashMapKMA_destroySHM(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid;
	key_t key;
	
	/* load sizes */
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
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
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	key = ftok(filename, 's');
	shmid = shmget(key, dest->seqsize * sizeof(long unsigned), 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, dest->v_index * sizeof(int), 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	key = ftok(filename, 'k');
	shmid = shmget(key, (dest->n + 1) * sizeof(unsigned), 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, dest->n * sizeof(unsigned), 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
}

void hashMapKMA_load(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid;
	long unsigned seekSize;
	key_t key;
	
	/* load sizes */
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	if(dest->kmersize <= 16) {
		getKmerP = &getK;
		addKmer = &addK;
	}
	
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	/* check shared memory, else load */
	seekSize = 0;
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), 0666);
	if(shmid < 0) {
		/* not shared, load */
		dest->exist = malloc(dest->size * sizeof(unsigned));
		if(!dest->exist) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		fread(dest->exist, sizeof(unsigned), dest->size, file);
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
		//fseek(file, dest->size * sizeof(unsigned), SEEK_CUR);
		seekSize += dest->size * sizeof(unsigned);
	}
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
		//fseek(file, dest->seqsize * sizeof(long unsigned), SEEK_CUR);
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
		//fseek(file, dest->v_index * sizeof(unsigned), SEEK_CUR);
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
		//fseek(file, dest->null_index * sizeof(unsigned), SEEK_CUR);
		seekSize += dest->n * sizeof(unsigned);
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

void hashMapKMA_load_old(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	/* load sizes */
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	/* allocate arrays */
	dest->exist = malloc(dest->size * sizeof(unsigned));
	dest->seq = malloc(dest->seqsize * sizeof(long unsigned));
	dest->values = malloc(dest->v_index * sizeof(unsigned));
	dest->key_index = malloc(dest->null_index * sizeof(unsigned));
	dest->value_index = malloc(dest->null_index * sizeof(unsigned));
	if(!dest->exist || !dest->seq || !dest->values || !dest->key_index || !dest->value_index) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	/* load arrays */
	fread(dest->exist, sizeof(unsigned), dest->size, file);
	fread(dest->seq, sizeof(long unsigned), dest->seqsize, file);
	fread(dest->values, sizeof(int), dest->v_index, file);
	fread(dest->key_index, sizeof(unsigned), dest->n + 1, file);
	fread(dest->value_index, sizeof(unsigned), dest->n, file);
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
	
	int index, pos;
	
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

int hashMap_index_getDub(struct hashMap_index *dest, long unsigned key, const char *qseq, int q_len) {
	
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
			if(getKmerIndex(dest->seq, pos - 1) == key) {
				dest->index[index] = 1 - dest->index[index];
				neg = -1;
			}
		} else {
			if(getKmerIndex(dest->seq, 1 - pos) == key) {
				neg = -1;
			}
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; index++) {
			if(pos > 0) {
				if(getKmerIndex(dest->seq, pos - 1) == key) {
					dest->index[index] = 1 - dest->index[index];
					neg = -1;
				}
			} else {
				if(getKmerIndex(dest->seq, 1 - pos) == key) {
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

void hashMap_index_dump(struct hashMap_index *src, FILE *seq, FILE *index) {
	
	fwrite(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
	fwrite(src->index, sizeof(int), src->size, index);
	
}

/* VALUES HASH */

struct valuesHash * initialize_hashValues(unsigned size) {
	
	struct valuesHash *dest;
	
	dest = malloc(sizeof(struct valuesHash));
	if(!dest) {
		OOM();
	}
	dest->n = 0;
	dest->size = size;
	
	dest->table = calloc(size, sizeof(struct valuesTable *));
	if(!dest->table) {
		OOM();
	}
	
	return dest;
}

void valuesHash_destroy(struct valuesHash *src) {
	
	unsigned i;
	struct valuesTable *node, *next;
	
	for(i = 0; i < src->size; i++) {
		for(node = src->table[i]; node != 0; node = next) {
			next = node->next;
			free(node);
		}
	}
	free(src->table);
	free(src);
}

unsigned valuesHash_add(struct valuesHash *dest, int *newValues, unsigned org_index, unsigned v_index) {
	
	long unsigned key;
	unsigned i, index;
	int *values;
	struct valuesTable *node;
	
	/* construct key */
	key = 0;
	for(i = 0; i <= newValues[0]; i++) {
		key = key * DB_size + newValues[i];
	}
	
	/* get index */
	index = key % dest->size;
	
	/* search for key */
	for(node = dest->table[index]; node != 0; node = node->next) {
		values = templates->values[node->values];
		if(values[0] == newValues[0] && int_eq(values + 1, newValues + 1, newValues[0])) { // Value exists
			return node->value_index;
		}
	}
	
	/* new values */
	dest->n++;
	node = malloc(sizeof(struct valuesTable));
	if(!node) {
		OOM();
	}
	node->value_index = v_index;
	node->values = org_index;
	node->next = dest->table[index];
	dest->table[index] = node;
	
	return v_index;
}

int hashMap_CountKmer(struct hashMap_kmers *dest, long unsigned key) {
	
	unsigned index;
	struct hashTable_kmers *node;
	
	index = key % dest->size;
	if(dest->table[index] == 0) {
		dest->table[index] = malloc(sizeof(struct hashTable_kmers));
		node = dest->table[index];
		node->key = key;
		node->next = 0;
		dest->n++;
		return 1;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(node->key == key) {
				return 0;
			} else if(node->next == 0) {
				node->next = malloc(sizeof(struct hashTable_kmers));
				node = node->next;
				node->key = key;
				node->next = 0;
				dest->n++;
				return 1;
			}
		}
	}
	return -1;
}

void emptyHash(struct hashMap_kmers *dest) {
	
	unsigned i;
	struct hashTable_kmers *node, *next;
	
	for(i = 0; i < dest->size; i++) {
		for(node = dest->table[i]; node != 0; node = next) {
			next = node->next;
			free(node);
		}
		dest->table[i] = 0;
	}
	dest->n = 0;
}

/*
	METHOD SPECIFIC
*/
int updateDBs(struct compDNA *qseq) {
	
	int i, j, end, extend;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	/* set parameters */
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	
	/* iterate sequence */
	for(i = 1, j = 0; i <= qseq->N[0]; i++) {
		extend = kmersize;
		end = qseq->N[i] - kmersize + 1;
		for(;j < end; j++) {
			/* update hashMap */
			extend = hashMap_add(getKmer(qseq->seq, j), DB_size, extend) ? 1 : kmersize;
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	
	return 1;
}

void makeIndexing(struct compDNA *compressor, FILE *seq_out, FILE *index_out) {
	
	int i, j, end;
	struct hashMap_index *template_index;
	
	/* allocate index */
	template_index = malloc(sizeof(struct hashMap_index));
	if(!template_index) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	template_index->len = compressor->seqlen;
	template_index->size = compressor->seqlen << 1;
	template_index->index = calloc(template_index->size, sizeof(int));
	if(!template_index->index) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	
	/* load index */
	template_index->seq = compressor->seq;
	compressor->N[0]++;
	compressor->N[compressor->N[0]] = compressor->seqlen + 1;
	j = 0;
	for(i = 1; i <= compressor->N[0]; i++) {
		end = compressor->N[i] - kmerindex;
		for(;j < end; j++) {
			hashMap_index_add(template_index, getKmerIndex(compressor->seq, j), j);
		}
		j = compressor->N[i] + 1;
	}
	compressor->N[0]--;
	
	/* dump index */
	hashMap_index_dump(template_index, seq_out, index_out);
	
	free(template_index->index);
	free(template_index);
}

int lengthCheck(struct compDNA *qseq) {
	
	int i, j, end, rc, thisKlen;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	thisKlen = MinKlen;
	
	for(rc = 0; rc < 2; rc++) {
		/* revers complement */
		if(rc) {
			rcComp(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0] && thisKlen != 0; i++) {
			end = qseq->N[i] - prefix_len - kmersize + 1;
			for(;j < end && thisKlen != 0; j++) {
				if(getPrefix(qseq->seq, j) == prefix) {
					thisKlen--;
				}
			}
			j = qseq->N[i] + 1;
		}
		qseq->N[0]--;
	}
	
	if(thisKlen) {
		return 0;
	} else {
		return 1;
	}
}

int queryCheck(struct compDNA *qseq) {
	
	int i, j, k, end, rc, thisKlen, *values;
	double bestQ, thisQ;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	thisKlen = 0;
	
	/* realloc */
	if(DB_size >= *Scores_tot) {
		free(Scores_tot);
		Scores_tot = calloc(2 * DB_size, sizeof(unsigned));
		free(bestTemplates);
		bestTemplates = malloc(2 * DB_size * sizeof(unsigned));
		if(!Scores_tot || !bestTemplates) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		Scores_tot[0] = 2 * DB_size;
	}
	
	/* get scores */
	bestTemplates[0] = 0;
	for(rc = 0; rc < 2; rc++) {
		/* revers complement */
		if(rc) {
			rcComp(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0]; i++) {
			end = qseq->N[i] - prefix_len - kmersize + 1;
			for(;j < end; j++) {
				if(getPrefix(qseq->seq, j) == prefix) {
					thisKlen++;
					if((values = hashMap_get(getKmer(qseq->seq, j + prefix_len)))) {
						for(k = 1; k <= *values; k++) {
							Scores_tot[values[k]]++;
							if(Scores_tot[values[k]] == 1) {
								bestTemplates[0]++;
								bestTemplates[bestTemplates[0]] = values[k];
							}
						}
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		qseq->N[0]--;
	}
	
	/* get query cov */
	bestQ = 0;
	for(i = 1; i <= *bestTemplates; i++) {
		thisQ = 1.0 * Scores_tot[bestTemplates[i]] / thisKlen;
		if(thisQ > bestQ) {
			bestQ = thisQ;
		}
		Scores_tot[bestTemplates[i]] = 0;
	}
	
	if(bestQ < homQ && thisKlen >= MinLen) {
		return 1;
	} else {
		return 0;
	}
}

int templateCheck(struct compDNA *qseq) {
	
	int i, j, k, end, rc, thisKlen, *values;
	double bestQ, thisQ, bestT, thisT;
	long unsigned key;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	thisKlen = 0;
	
	/* realloc */
	if(DB_size >= *Scores_tot) {
		free(Scores);
		Scores = calloc(2 * DB_size, sizeof(unsigned));
		free(Scores_tot);
		Scores_tot = calloc(2 * DB_size, sizeof(unsigned));
		free(bestTemplates);
		bestTemplates = malloc(2 * DB_size * sizeof(unsigned));
		if(!Scores || !Scores_tot || !bestTemplates) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		Scores_tot[0] = 2 * DB_size;
	}
	
	/* get scores */
	bestTemplates[0] = 0;
	for(rc = 0; rc < 2; rc++) {
		/* revers complement */
		if(rc) {
			rcComp(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0]; i++) {
			end = qseq->N[i] - prefix_len - kmersize + 1;
			for(;j < end; j++) {
				if(getPrefix(qseq->seq, j) == prefix) {
					thisKlen++;
					key = getKmer(qseq->seq, j + prefix_len);
					if((values = hashMap_get(key))) {
						for(k = 1; k <= *values; k++) {
							Scores_tot[values[k]]++;
							if(Scores_tot[values[k]] == 1) {
								bestTemplates[0]++;
								bestTemplates[bestTemplates[0]] = values[k];
							}
						}
						if(hashMap_CountKmer(foundKmers, key)) {
							for(k = 1; k <= *values; k++) {
								Scores[values[k]]++;
							}
						}
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		qseq->N[0]--;
	}
	
	/* get query cov */
	bestQ = 0;
	bestT = 0;
	for(i = 1; i <= *bestTemplates; i++) {
		thisQ = 1.0 * Scores_tot[bestTemplates[i]] / thisKlen;
		if(thisQ > bestQ) {
			bestQ = thisQ;
		}
		thisT = 1.0 * Scores[bestTemplates[i]] / template_ulengths[bestTemplates[i]];
		if(thisT > bestT) {
			bestT = thisT;
		}
		Scores_tot[bestTemplates[i]] = 0;
		Scores[bestTemplates[i]] = 0;
	}
	
	emptyHash(foundKmers);
	
	if(homcmp(bestT < homT, bestQ < homQ) && thisKlen >= MinLen) {
		/* realloc hash */
		if(thisKlen > foundKmers->size) {
			foundKmers->n = 0;
			foundKmers->size = 2 * thisKlen;
			free(foundKmers->table);
			foundKmers->table = calloc(foundKmers->size, sizeof(struct hashTable_kmers *));
			if(!foundKmers->table) {
				fprintf(stderr, "OOM\n");
				exit(-1);
			}
		}
		return 1;
	} else {
		return 0;
	}
}

int updateDBs_sparse(struct compDNA *qseq) {
	
	int i, j, end, extend, last, rc;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	/* test homology and length */
	if(QualCheck(qseq)) {
		template_slengths[DB_size] = 0;
		template_ulengths[DB_size] = 0;
		for(rc = 0; rc < 2; rc++) {
			/* revers complement */
			if(rc) {
				rcComp(qseq);
			}
			/* set last extender */
			last = -kmersize;
			extend = kmersize;
			
			/* iterate seq */
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			j = 0;
			for(i = 1; i <= qseq->N[0]; i++) {
				end = qseq->N[i] - prefix_len - kmersize + 1;
				for(;j < end; j++) {
					if(getPrefix(qseq->seq, j) == prefix) {
						/* add kmer */
						extend = kmersize < (j - last) ? kmersize : (j - last);
						last = hashMap_add(getKmer(qseq->seq, j + prefix_len), DB_size, extend) ? j : -kmersize;
						template_slengths[DB_size]++;
					}
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		return 1;
	}
	
	return 0;
}

struct hashMapKMA * compressKMA_DB_old(FILE *out) {
	
	long unsigned i, j;
	unsigned index, t_index, v_index, v_update, c_index, new_index, null_index;
	unsigned *tmp;
	int *values;
	struct hashMapKMA *finalDB;
	struct valuesHash *shmValues;
	struct hashTable *node, *next, *table;
	
	/* cut templates down */
	fprintf(stderr, "# Resizing DB.\n");
	templates->seq_size = (templates->seq_n / 32 + 1);
	templates->seq = realloc(templates->seq, templates->seq_size * sizeof(long unsigned));
	templates->values = realloc(templates->values, templates->n * sizeof(unsigned *));
	if(!templates->seq || !templates->values) {
		OOM();
	}
	table = 0;
	for(i = 0; i < templates->size; i++) {
		for(node = templates->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(templates->table);
	templates->table = 0;
	
	/* prepare final DB */
	fprintf(stderr, "# Preparing compressed DB.\n");
	finalDB = malloc(sizeof(struct hashMapKMA));
	if(!finalDB) {
		OOM();
	}
	/* Fill in known values */
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->seq = templates->seq;
	finalDB->seqsize = templates->seq_size;
	finalDB->prefix_len = prefix_len;
	finalDB->prefix = prefix;
	finalDB->kmersize = kmersize;
	
	/* allocate existence */
	finalDB->exist = malloc(finalDB->size * sizeof(unsigned));
	finalDB->key_index = malloc((finalDB->n + 1) * sizeof(unsigned));
	finalDB->value_index = malloc(finalDB->n * sizeof(unsigned));
	tmp = malloc(finalDB->n * sizeof(unsigned));
	if(!finalDB->exist || !finalDB->key_index || !finalDB->value_index || !tmp) {
		OOM();
	}
	null_index = finalDB->n;
	finalDB->null_index = null_index;
	
	/* mv table to finalDB */
	fprintf(stderr, "# Initialize cp of DB.\n");
	for(i = 0; i < finalDB->size; i++) {
		finalDB->exist[i] = null_index;
	}
	fprintf(stderr, "# Initial cp of DB.\n");
	node = table;
	finalDB->size--;
	t_index = 0;
	while(node != 0) {
		/* get index */
		index = getKmerP(finalDB->seq, node->key) & finalDB->size;
		finalDB->exist[index] = t_index;
		
		/* mv chain */
		while(node != 0 && (getKmerP(finalDB->seq, node->key) & finalDB->size) == index) {
			next = node->next;
			
			/* cp index */
			finalDB->key_index[t_index] = node->key;
			finalDB->value_index[t_index] = node->values;
			tmp[t_index] = node->values;
			t_index++;
			
			/* clean */
			free(node);
			node = next;
		}
	}
	
	/* get compressed indexes */
	fprintf(stderr, "# Compressing indexes.\n");
	v_index = 0;
	c_index = 0;
	shmValues = initialize_hashValues(null_index);
	for(i = 0; i < null_index; i++) {
		/* potential increase */
		v_update = templates->values[finalDB->value_index[i]][0] + 1;
		
		/* the actual index */
		new_index = valuesHash_add(shmValues, templates->values[finalDB->value_index[i]], finalDB->value_index[i], v_index);
		
		/* update to new index */
		finalDB->value_index[i] = new_index;
		
		/* update size of compression */
		if(new_index == v_index) {
			c_index += v_update;
			if(v_index < c_index) {
				v_index = c_index;
			} else {
				fprintf(stderr, "Compression overflow.\n");
				exit(-1);
			}
		}
	}
	/* destroy shmValues */
	valuesHash_destroy(shmValues);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->values = calloc(v_index, sizeof(unsigned));
	if(!finalDB->values) {
		OOM();
	}
	for(i = 0; i < finalDB->n; i++) {
		if(finalDB->values[finalDB->value_index[i]] == 0) {
			values = templates->values[tmp[i]];
			v_index = finalDB->value_index[i];
			for(j = 0; j <= values[0]; j++) {
				finalDB->values[v_index + j] = values[j];
			}
		}
	}
	
	/* add terminating key */
	i = 0;
	j = getKmerP(templates->seq, finalDB->key_index[finalDB->n - 1]) & finalDB->size;
	while(j == (getKmerP(templates->seq, i) & finalDB->size)) {
		i++;
	}
	finalDB->key_index[finalDB->n] = i;
	
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");	
	finalDB->size++;
	hashMapKMA_dump(finalDB, out);
	
	/* clean */
	free(finalDB->value_index);
	finalDB->value_index = tmp;
	
	return finalDB;
}

struct hashMapKMA * compressKMA_DB(FILE *out) {
	
	long unsigned i, j;
	unsigned index, t_index, v_index, v_update, c_index, new_index, null_index;
	unsigned *tmp;
	int *values;
	struct hashMapKMA *finalDB;
	struct valuesHash *shmValues;
	struct hashTable *node, *next, *table;
	
	/* cut templates down */
	fprintf(stderr, "# Resizing DB.\n");
	templates->seq_size = (templates->seq_n / 32 + 1);
	templates->seq = realloc(templates->seq, templates->seq_size * sizeof(long unsigned));
	templates->values = realloc(templates->values, templates->n * sizeof(unsigned *));
	if(!templates->seq || !templates->values) {
		OOM();
	}
	table = 0;
	for(i = 0; i < templates->size; i++) {
		for(node = templates->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(templates->table);
	templates->table = 0;
	
	/* prepare final DB */
	fprintf(stderr, "# Preparing compressed DB.\n");
	finalDB = malloc(sizeof(struct hashMapKMA));
	if(!finalDB) {
		OOM();
	}
	/* Fill in known values */
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->seq = templates->seq;
	finalDB->seqsize = templates->seq_size;
	finalDB->prefix_len = prefix_len;
	finalDB->prefix = prefix;
	finalDB->kmersize = kmersize;
	
	/* allocate existence */
	finalDB->exist = malloc(finalDB->size * sizeof(unsigned));
	finalDB->key_index = malloc((finalDB->n + 1) * sizeof(unsigned));
	finalDB->value_index = malloc(finalDB->n * sizeof(unsigned));
	tmp = malloc(finalDB->n * sizeof(unsigned));
	if(!finalDB->exist || !finalDB->key_index || !finalDB->value_index || !tmp) {
		OOM();
	}
	null_index = finalDB->n;
	finalDB->null_index = null_index;
	
	/* mv table to finalDB */
	fprintf(stderr, "# Initialize cp of DB.\n");
	for(i = 0; i < finalDB->size; i++) {
		finalDB->exist[i] = null_index;
	}
	fprintf(stderr, "# Initial cp of DB.\n");
	node = table;
	finalDB->size--;
	t_index = 0;
	while(node != 0) {
		/* get index */
		index = getKmerP(finalDB->seq, node->key) & finalDB->size;
		finalDB->exist[index] = t_index;
		
		/* mv chain */
		while(node != 0 && (getKmerP(finalDB->seq, node->key) & finalDB->size) == index) {
			next = node->next;
			
			/* cp index */
			finalDB->key_index[t_index] = node->key;
			finalDB->value_index[t_index] = node->values;
			tmp[t_index] = node->values;
			t_index++;
			
			/* clean */
			free(node);
			node = next;
		}
	}
	
	/* get compressed indexes */
	fprintf(stderr, "# Compressing indexes.\n");
	v_index = 0;
	c_index = 0;
	shmValues = initialize_hashValues(null_index);
	for(i = 0; i < null_index; i++) {
		/* potential increase */
		v_update = templates->values[finalDB->value_index[i]][0] + 1;
		
		/* the actual index */
		new_index = valuesHash_add(shmValues, templates->values[finalDB->value_index[i]], finalDB->value_index[i], v_index);
		
		/* update to new index */
		finalDB->value_index[i] = new_index;
		
		/* update size of compression */
		if(new_index == v_index) {
			c_index += v_update;
			if(v_index < c_index) {
				v_index = c_index;
			} else {
				fprintf(stderr, "Compression overflow.\n");
				exit(-1);
			}
		}
	}
	/* destroy shmValues */
	valuesHash_destroy(shmValues);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->values = calloc(v_index, sizeof(unsigned));
	if(!finalDB->values) {
		OOM();
	}
	for(i = 0; i < finalDB->n; i++) {
		if(finalDB->values[finalDB->value_index[i]] == 0) {
			values = templates->values[tmp[i]];
			v_index = finalDB->value_index[i];
			for(j = 0; j <= values[0]; j++) {
				finalDB->values[v_index + j] = values[j];
			}
		}
	}
	
	/* add terminating key */
	i = 0;
	j = getKmerP(templates->seq, finalDB->key_index[finalDB->n - 1]) & finalDB->size;
	while(j == (getKmerP(templates->seq, i) & finalDB->size)) {
		i++;
	}
	finalDB->key_index[finalDB->n] = i;
	
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");	
	finalDB->size++;
	hashMapKMA_dump(finalDB, out);
	
	/* clean */
	free(finalDB->value_index);
	finalDB->value_index = tmp;
	
	return finalDB;
}

struct hashMapKMA * compressKMA_megaDB(FILE *out) {
	
	long unsigned i, j;
	unsigned v_index, v_update, c_index, new_index, null_index, *tmp;
	int *values;
	struct hashMapKMA *finalDB;
	struct valuesHash *shmValues;
	
	/* cut templates down */
	fprintf(stderr, "# Resizing DB.\n");
	finalDB = malloc(sizeof(struct hashMapKMA));
	if(!finalDB) {
		OOM();
	}
	/* Fill in known values */
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->seq = 0;
	finalDB->seqsize = 0;
	finalDB->prefix_len = prefix_len;
	finalDB->prefix = prefix;
	finalDB->kmersize = kmersize;
	/* allocate existence */
	finalDB->exist = malloc(finalDB->size * sizeof(unsigned));
	finalDB->key_index = 0;
	finalDB->value_index = 0;
	if(!finalDB->exist) {
		OOM();
	}
	null_index = finalDB->n;
	finalDB->null_index = null_index;
	for(i = 0; i < finalDB->size; i++) {
		finalDB->exist[i] = null_index;
	}
	v_index = 0;
	while(templates->values[v_index] != 0) {
		v_index++;
	}
	for(i = v_index; i < finalDB->size; i++) {
		if(templates->values[i] != 0) {
			finalDB->exist[i] = v_index;
			templates->values[v_index] = templates->values[i];
			templates->values[i] = 0;
			v_index++;
		}
	}
	templates->values = realloc(templates->values, templates->n * sizeof(unsigned *));
	if(!templates->values) {
		OOM();
	}
	
	/* mv table to finalDB */
	fprintf(stderr, "# Initial cp of DB.\n");
	tmp = malloc(finalDB->size * sizeof(unsigned));
	if(!tmp) {
		OOM();
	}
	for(i = 0; i < finalDB->size; i++) {
		tmp[i] = finalDB->exist[i];
	}
	
	/* get compressed indexes */
	fprintf(stderr, "# Compressing indexes.\n");
	v_index = 0;
	c_index = 0;
	shmValues = initialize_hashValues(null_index);
	for(i = 0; i < finalDB->size; i++) {
		if(finalDB->exist[i] != null_index) {
			/* potential increase */
			v_update = templates->values[finalDB->exist[i]][0] + 1;
			
			/* the actual index */
			new_index = valuesHash_add(shmValues, templates->values[finalDB->exist[i]], finalDB->exist[i], v_index);
			
			/* update to new index */
			finalDB->exist[i] = new_index;
			
			/* update size of compression */
			if(new_index == v_index) {
				c_index += v_update;
				if(v_index < c_index) {
					v_index = c_index;
				} else {
					fprintf(stderr, "Compression overflow.\n");
					exit(-1);
				}
			}
		} else {
			finalDB->exist[i] = 1;
		}
	}
	/* destroy shmValues */
	valuesHash_destroy(shmValues);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->null_index = 1;
	finalDB->values = calloc(v_index, sizeof(unsigned));
	if(!finalDB->values) {
		OOM();
	}
	for(i = 0; i < finalDB->size; i++) {
		if(finalDB->exist[i] != finalDB->null_index && finalDB->values[finalDB->exist[i]] == 0) {
			values = templates->values[tmp[i]];
			v_index = finalDB->exist[i];
			for(j = 0; j <= values[0]; j++) {
				finalDB->values[v_index + j] = values[j];
			}
		}
	}
	
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");
	megaMapKMA_dump(finalDB, out);
	
	/* clean */
	free(finalDB->exist);
	finalDB->exist = tmp;
	
	return finalDB;
}

void compressKMA_deconDB(struct hashMapKMA *finalDB) {
	
	long unsigned i, j;
	unsigned v_index, c_index, v_update, new_index, *tmp;
	int *values;
	struct valuesHash *shmValues;
	
	fprintf(stderr, "# Compressing indexes.\n");
	tmp = malloc(finalDB->n * sizeof(unsigned));
	if(!tmp) {
		OOM();
	}
	v_index = 0;
	c_index = 0;
	shmValues = initialize_hashValues(finalDB->n);
	for(i = 0; i < finalDB->n; i++) {
		/* potential increase */
		v_update = templates->values[finalDB->value_index[i]][0] + 1;
		
		/* the actual index */
		new_index = valuesHash_add(shmValues, templates->values[finalDB->value_index[i]], finalDB->value_index[i], v_index);
		
		/* update to new index */
		tmp[i] = finalDB->value_index[i];
		finalDB->value_index[i] = new_index;
		
		/* update size of compression */
		if(new_index == v_index) {
			c_index += v_update;
			if(v_index < c_index) {
				v_index = c_index;
			} else {
				fprintf(stderr, "Compression overflow.\n");
				exit(-1);
			}
		}
	}
	/* destroy shmValues */
	valuesHash_destroy(shmValues);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->values = calloc(v_index, sizeof(unsigned));
	if(!finalDB->values) {
		OOM();
	}
	for(i = 0; i < finalDB->n; i++) {
		if(finalDB->values[finalDB->value_index[i]] == 0) {
			values = templates->values[tmp[i]];
			v_index = finalDB->value_index[i];
			for(j = 0; j <= values[0]; j++) {
				finalDB->values[v_index + j] = values[j];
			}
			free(templates->values[tmp[i]]);
		}
	}
	
	/* clean */
	free(tmp);
	free(templates->values);
	free(templates);
}

void compressKMA_deconMegaDB(struct hashMapKMA *finalDB) {
	
	long unsigned i, j;
	unsigned v_index, c_index, v_update, new_index, *tmp;
	int *values;
	struct valuesHash *shmValues;
	
	fprintf(stderr, "# Compressing indexes.\n");
	tmp = malloc(finalDB->size * sizeof(unsigned));
	if(!tmp) {
		OOM();
	}
	v_index = 0;
	c_index = 0;
	shmValues = initialize_hashValues(finalDB->n);
	for(i = 0; i < finalDB->size; i++) {
		if(finalDB->exist[i] != finalDB->n) {
			/* potential increase */
			v_update = templates->values[finalDB->exist[i]][0] + 1;
			
			/* the actual index */
			new_index = valuesHash_add(shmValues, templates->values[finalDB->exist[i]], finalDB->exist[i], v_index);
			
			/* update to new index */
			tmp[i] = finalDB->exist[i];
			finalDB->exist[i] = new_index;
			
			/* update size of compression */
			if(new_index == v_index) {
				c_index += v_update;
				if(v_index < c_index) {
					v_index = c_index;
				} else {
					fprintf(stderr, "Compression overflow.\n");
					exit(-1);
				}
			} else {
				free(templates->values[tmp[i]]);
			}
		} else {
			tmp[i] = finalDB->n;
		}
	}
	/* destroy shmValues */
	valuesHash_destroy(shmValues);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->values = calloc(v_index, sizeof(unsigned));
	if(!finalDB->values) {
		OOM();
	}
	for(i = 0; i < finalDB->size; i++) {
		if(finalDB->exist[i] != finalDB->n) {
			if(finalDB->values[finalDB->exist[i]] == 0) {
				values = templates->values[tmp[i]];
				v_index = finalDB->exist[i];
				for(j = 0; j <= values[0]; j++) {
					finalDB->values[v_index + j] = values[j];
				}
				free(templates->values[tmp[i]]);
			}
		} else {
			finalDB->exist[i] = finalDB->null_index;
		}
	}
	
	/* clean */
	free(tmp);
	free(templates->values);
	free(templates);
}

void hashMapKMA_destroy(struct hashMapKMA *dest) {
	free(dest->exist);
	free(dest->seq);
	free(dest->values);
	free(dest->key_index);
	free(dest->value_index);
	free(dest);
}

void hashMapKMA_NULL(struct hashMapKMA *dest) {
	dest->prefix = 0;
	dest->prefix_len = 0;
	free(dest->exist);
	free(dest->seq);
	free(dest->values);
	free(dest->key_index);
	free(dest->value_index);
	dest->size = 0;
	dest->n = 0;
}

void updateAnnots(struct compDNA *qseq, FILE *seq_out, FILE *index_out) {
	
	/* Dump annots */
	makeIndexing(qseq, seq_out, index_out);
	
	template_lengths[DB_size] = qseq->seqlen;
	DB_size++;
	if(DB_size >= template_lengths[0]) {
		template_lengths[0] *= 2;
		template_lengths = realloc(template_lengths, template_lengths[0] * sizeof(unsigned));
		if(!template_lengths) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
}

void updateAnnots_sparse(struct compDNA *qseq, FILE *seq_out, FILE *index_out) {
	
	/* Dump nibble seq */
	fwrite(qseq->seq, sizeof(long unsigned), (qseq->seqlen >> 5) + 1, seq_out);
	
	template_lengths[DB_size] = qseq->seqlen;
	DB_size++;
	if(DB_size >= template_ulengths[0]) {
		template_ulengths[0] *= 2;
		template_slengths = realloc(template_slengths, template_ulengths[0] * sizeof(unsigned));
		template_ulengths = realloc(template_ulengths, template_ulengths[0] * sizeof(unsigned));
		template_lengths = realloc(template_lengths, template_ulengths[0] * sizeof(unsigned));
		if(!template_lengths || !template_slengths || !template_ulengths) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
	}
	
}

int deConNode(struct compDNA *qseq, struct hashMapKMA *finalDB) {
	
	int i, j, end, mapped_cont;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	mapped_cont = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	j = 0;
	for(i = 1; i <= qseq->N[0]; i++) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end; j++) {
			mapped_cont += addCont(finalDB, getKmer(qseq->seq, j), DB_size);
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	return mapped_cont;
}

int deConNode_sparse(struct compDNA *qseq, struct hashMapKMA *finalDB) {
	
	int i, j, end, mapped_cont;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	mapped_cont = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	j = 0;
	for(i = 1; i <= qseq->N[0]; i++) {
		end = qseq->N[i] - prefix_len - kmersize + 1;
		for(;j < end; j++) {
			if(getPrefix(qseq->seq, j) == prefix) {
				mapped_cont += addCont(finalDB, getKmer(qseq->seq, j + prefix_len), DB_size);
			}
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	return mapped_cont;
}

unsigned deConDB(struct hashMapKMA *finalDB, char **inputfiles, int fileCount, char *outputfilename) {
	
	unsigned fileCounter, mapped_cont, file_len;
	char *filename, *zipped, *cmd;
	FILE *DB_update;
	struct qseqs *header, *qseq;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	
	/* allocate */
	compressor = malloc(sizeof(struct compDNA));
	cmd = malloc(1);
	zipped = strdup(".gz");
	if(!cmd || !zipped || !compressor) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* set variables */
	mapped_cont = 0;
	finalDB->size--;
	
	/* open files */
	file_len = strlen(outputfilename);
	strcat(outputfilename, 	".decon.b");
	DB_update = fopen(outputfilename, "wb");
	outputfilename[file_len] = 0;
	if(!DB_update) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	
	/* iterate inputfiles */
	for(fileCounter = 0; fileCounter < fileCount; fileCounter++) {
		/* open file */
		filename = inputfiles[fileCounter];
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
		
		/* Get first char and determine the format */
		buffFileBuff(inputfile);
		if(inputfile->buffer[0] != '>') { //FASTA
			fprintf(stderr, "Not Fasta format!!!\n");
			exit(-1);
		}
		
		/* parse the file */
		while(FileBuffgetFsa(inputfile, header, qseq)) {
			fprintf(stderr, "# Decon:\t%s\n", header->seq + 1);
			if(qseq->len > kmersize) {
				/* compress DNA */
				if(qseq->len >= compressor->size) {
					freeComp(compressor);
					allocComp(compressor, qseq->len);
				}
				compDNAref(compressor, qseq->seq, qseq->len);
				
				/* Add contamination */
				mapped_cont += deConNode_ptr(compressor, finalDB);
				/* rc */
				rcComp(compressor);
				mapped_cont += deConNode_ptr(compressor, finalDB);
			}
		}
		
		/* close file buffer */
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			pcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
	}
	finalDB->size++;
	
	/* dump DB */
	fprintf(stderr, "# Dumping DeCon DB.\n");
	//deConMap_dump(finalDB, DB_update);
	fclose(DB_update);
	
	/* clean */
	free(cmd);
	free(zipped);
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyFileBuff(inputfile);
	
	return mapped_cont;
}

int homcmp_or(int t, int q) {
	return (t || q);
}

int homcmp_and(int t, int q) {
	return (t && q);
}

void makeDB(char **inputfiles, int fileCount, char *outputfilename, int appender) {
	
	int fileCounter, file_len, bias;
	char *filename, *zipped, *cmd;
	FILE *index_out, *seq_out, *length_out, *name_out, *DB_update;
	struct qseqs *header, *qseq;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	
	/* allocate */
	compressor = malloc(sizeof(struct compDNA));
	cmd = malloc(1);
	zipped = strdup(".gz");
	if(!cmd || !zipped || !compressor) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* open files */
	file_len = strlen(outputfilename);
	strcat(outputfilename, 	".b");
	DB_update = fopen(outputfilename, "wb");
	outputfilename[file_len] = 0;
	strcat(outputfilename, 	".length.b");
	length_out = fopen(outputfilename, "wb");
	outputfilename[file_len] = 0;
	if(appender) {
		strcat(outputfilename, 	".name");
		name_out = fopen(outputfilename, "a");
		outputfilename[file_len] = 0;
	} else {
		strcat(outputfilename, 	".name");
		name_out = fopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	}
	if(!DB_update || !length_out || !name_out) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	}
	if(appender) {
		strcat(outputfilename, ".seq.b");
		seq_out = fopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
		if(!seq_out) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			exit(-1);
		}
	} else {
		strcat(outputfilename, ".seq.b");
		seq_out = fopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		if(!seq_out) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			exit(-1);
		}
	}
	
	if(template_ulengths == 0) {
		if(appender) {
			strcat(outputfilename, ".index.b");
			index_out = fopen(outputfilename, "ab");
			outputfilename[file_len] = 0;
			if(!index_out) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				exit(-1);
			}
		} else {
			strcat(outputfilename, ".index.b");
			index_out = fopen(outputfilename, "wb");
			outputfilename[file_len] = 0;
			if(!index_out) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
				exit(-1);
			}
			fwrite(&kmerindex, sizeof(int), 1, index_out);
		}
	} else {
		index_out = 0;
	}
	
	fprintf(stderr, "# Updating DBs\n");
	/* iterate inputfiles */
	for(fileCounter = 0; fileCounter < fileCount; fileCounter++) {
		/* open file */
		filename = inputfiles[fileCounter];
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
		
		/* Get first char and determine the format */
		buffFileBuff(inputfile);
		if(inputfile->buffer[0] != '>') { //FASTA
			fprintf(stderr, "Not Fasta format!!!\n");
			exit(-1);
		}
		
		/* parse the file */
		while(FileBuffgetFsa(inputfile, header, qseq)) {
			if(qseq->len >= compressor->size) {
				freeComp(compressor);
				allocComp(compressor, 2 * qseq->len);
			}
			bias = compDNAref(compressor, qseq->seq, qseq->len);
			if(qseq->len > MinLen && update_DB(compressor)) {
				/* Update annots */
				chomp(header->seq);
				if(bias > 0) {
					fprintf(name_out, "%s B%d\n", header->seq + 1, bias);
				} else {
					fprintf(name_out, "%s\n", header->seq + 1);
				}
				updateAnnotsPtr(compressor, seq_out, index_out);
				fprintf(stderr, "# Added:\t%s\n", header->seq + 1);
			} else {
				fprintf(stderr, "# Skipped:\t%s\n", header->seq + 1);
			}
		}
		
		/* close file buffer */
		if(strncmp(filename + (strlen(filename) - 3), zipped, 3) == 0) {
			pcloseFileBuff(inputfile);
		} else if(strncmp(filename, "--", 2) != 0) {
			closeFileBuff(inputfile);
		}
		
	}
	
	/* Dump annots */
	fwrite(&DB_size, sizeof(int), 1, length_out);
	if(template_ulengths != 0) {
		template_ulengths[0] = 0;
		template_slengths[0] = 0;
		fwrite(template_slengths, sizeof(unsigned), DB_size, length_out);
		fwrite(template_ulengths, sizeof(unsigned), DB_size, length_out);
		fwrite(template_lengths, sizeof(unsigned), DB_size, length_out);
		fclose(length_out);
	} else {
		template_lengths[0] = kmerindex;
		fwrite(template_lengths, sizeof(unsigned), DB_size, length_out);
		fclose(index_out);
		fclose(seq_out);
		fclose(length_out);
	}
	fclose(name_out);
	
	
	fprintf(stderr, "# Templates key-value pairs:\t%u.\n", templates->n);// / 1048576);
	hashMap_dump(templates, DB_update);
	fclose(DB_update);
	
	/* clean */
	free(cmd);
	free(zipped);
	freeComp(compressor);
	free(compressor);
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyFileBuff(inputfile);
	
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma_index creates the databases needed to run KMA, from a list of fasta files given.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\t\tDefault:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-i\t\tInput/query file name (STDIN: \"--\")\tNone\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\t\tInput/template file\n");
	fprintf(helpOut, "#\t-batch\t\tBatch input file\n");
	fprintf(helpOut, "#\t-deCon\t\tFile with contamination (STDIN: \"--\")\tNone/False\n");
	fprintf(helpOut, "#\t-batchD\t\tBatch decon file\n");
	fprintf(helpOut, "#\t-t_db\t\tAdd to existing DB\t\t\tNone/False\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t\t16\n");
	fprintf(helpOut, "#\t-k_t\t\tKmersize for template identification\t16\n");
	fprintf(helpOut, "#\t-k_i\t\tKmersize for indexing\t\t\t16\n");
	fprintf(helpOut, "#\t-ML\t\tMinimum length of templates\t\tkmersize (16)\n");
	fprintf(helpOut, "#\t-CS\t\tStart Chain size\t\t\t\t1 M\n");
	fprintf(helpOut, "#\t-ME\t\tMega DB\t\t\t\t\tFalse\n");
	fprintf(helpOut, "#\t-Sparse\t\tMake Sparse DB ('-' for no prefix)\tNone/False\n");
	fprintf(helpOut, "#\t-ht\t\tHomology template\t\t\t1.0\n");
	fprintf(helpOut, "#\t-hq\t\tHomology query\t\t\t\t1.0\n");
	fprintf(helpOut, "#\t-and\t\tBoth homolgy thresholds\n#\t\t\thas to be reached\t\t\tor\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int i, args, stop, filecount, deconcount, sparse_run;
	int line_size, l_len, mapped_cont, file_len, appender;
	unsigned megaDB;
	char **inputfiles, *outputfilename, *templatefilename, **deconfiles, *line;
	struct hashMapKMA *finalDB;
	FILE *inputfile, *out;
	
	if (argc == 1) {
		fprintf(stderr, "# Too few arguments handed.\n");
		helpMessage(-1);
	} else if(sizeof(long unsigned) != 8) {
		fprintf(stderr, "Need a 64-bit system.\n");
		exit(-1);
	}
	
	/* set defaults */
	INITIAL_SIZE = 1048576;
	kmersize = 16;
	kmerindex = 16;
	sparse_run = 0;
	appender = 0;
	MinLen = kmersize - 1;
	MinKlen = 1;
	prefix_len = 0;
	prefix = 0;
	homQ = 1;
	homT = 1;
	homcmp = &homcmp_or;
	getKmerP = &getKmer;
	addKmer = &addKmerC;
	template_ulengths = 0;
	template_slengths = 0;
	DB_size = 1;
	filecount = 0;
	deconcount = 0;
	outputfilename = 0;
	templatefilename = 0;
	megaDB = 0;
	inputfiles = malloc(sizeof(char*));
	deconfiles = malloc(sizeof(char*));
	line_size = 256;
	line = malloc(line_size);
	to2Bit = malloc(384);
	if(!inputfiles || !deconfiles || !line || !to2Bit) {
		fprintf(stderr, "OOM\n");
		exit(-1);
	}
	/* set to2Bit */
	for(i = 0; i < 384; i++) {
		to2Bit[i] = 5;
	}
	to2Bit += 128;
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
	/*
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
	*/
	to2Bit['R'] = 0;
	to2Bit['Y'] = 1;
	to2Bit['S'] = 2;
	to2Bit['W'] = 3;
	to2Bit['K'] = 2;
	to2Bit['M'] = 0;
	to2Bit['B'] = 1;
	to2Bit['D'] = 0;
	to2Bit['H'] = 3;
	to2Bit['V'] = 2;
	to2Bit['X'] = 4;
	to2Bit['r'] = 0;
	to2Bit['y'] = 1;
	to2Bit['s'] = 2;
	to2Bit['w'] = 3;
	to2Bit['k'] = 2;
	to2Bit['m'] = 0;
	to2Bit['b'] = 1;
	to2Bit['d'] = 0;
	to2Bit['h'] = 3;
	to2Bit['v'] = 2;
	to2Bit['x'] = 4;
	
	/* Future RNA encoding */
	to2Bit['U'] = 3;
	to2Bit['u'] = 3;
	
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			stop = 0;
			args++;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strcmp(argv[args], "--") == 0) {
					filecount++;
					inputfiles = realloc(inputfiles, filecount * sizeof(char*));
					if(inputfiles == NULL) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					inputfiles[filecount - 1] = strdup(argv[args]);
					if(!inputfiles[filecount - 1]) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					args++;
				} else {
					stop = 1;
				}
			}
			args--;
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
			stop = 0;
			args++;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strcmp(argv[args], "--") == 0) {
					deconfiles = realloc(deconfiles, (deconcount + 1) * sizeof(char*));
					if(deconfiles == NULL) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					deconfiles[deconcount] = strdup(argv[args]);
					if(deconfiles[deconcount] == NULL) {
						fprintf(stderr, "OOM\n");
						exit(-1);
					}
					deconcount++;
					args++;
				} else {
					stop = 1;
				}
			}
			if(deconcount == 0) {
				fprintf(stderr, "No deCon file specified.\n");
				exit(-1);
			}
			args--;
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
				} else if(kmersize > 32) {
					kmersize = 32;
				}
				kmerindex = kmersize;
			}
		} else if(strcmp(argv[args], "-k_t") == 0) {
			args++;
			if(args < argc) {
				kmersize = atoi(argv[args]);
				if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 32) {
					kmersize = 32;
				}
			}
		} else if(strcmp(argv[args], "-k_i") == 0) {
			args++;
			if(args < argc) {
				kmerindex = atoi(argv[args]);
				if(kmerindex == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmerindex = 16;
				} else if(kmerindex > 32) {
					kmerindex = 32;
				}
			}
		} else if(strcmp(argv[args], "-CS") == 0) {
			args++;
			if(args < argc) {
				
				INITIAL_SIZE = pow(2, ceil(log(atoi(argv[args]))/log(2))) + 0.5;
				INITIAL_SIZE *= 1048576;
				if(INITIAL_SIZE == 0) {
					fprintf(stderr, "# Invalid Chain Size parsed, using default\n");
					INITIAL_SIZE = 1048576;
				}
			}
		} else if(strcmp(argv[args], "-and") == 0) {
			homcmp = &homcmp_and;
		} else if(strcmp(argv[args], "-ML") == 0) {
			args++;
			if(args < argc) {
				MinLen = atoi(argv[args]);
				if(MinLen <= 0) {
					fprintf(stderr, "# Invalid minimum length parsed, using default\n");
				}
			}
		} else if(strcmp(argv[args], "-hq") == 0) {
			args++;
			if(args < argc) {
				homQ = atof(argv[args]);
				if(homQ < 0) {
					fprintf(stderr, "Invalid -hq\n");
					homQ = 1.0;
				}
			}
		} else if(strcmp(argv[args], "-ht") == 0) {
			args++;
			if(args < argc) {
				homT = atof(argv[args]);
				if(homT < 0) {
					fprintf(stderr, "Invalid -hq\n");
					homT = 1.0;
				}
			}
		} else if(strcmp(argv[args], "-batch") == 0) {
			args++;
			if(args < argc) {
				inputfile = fopen(argv[args], "r");
				if(!inputfile) {
					fprintf(stderr, "No such file:\t%s\n", argv[args]);
					exit(-1);
				}
				while(!feof(inputfile) && *(line = fget_line(line, &line_size, &l_len, inputfile))) {
					if(l_len != 0) {
						filecount++;
						inputfiles = realloc(inputfiles, filecount * sizeof(char*));
						if(inputfiles == NULL) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
						inputfiles[filecount - 1] = strdup(line);
						if(inputfiles[filecount - 1] == NULL) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
					}
				}
				fclose(inputfile);
			}
		} else if(strcmp(argv[args], "-batchD") == 0) {
			args++;
			if(args < argc) {
				inputfile = fopen(argv[args], "r");
				if(!inputfile) {
					fprintf(stderr, "No such file:\t%s\n", argv[args]);
					exit(-1);
				}
				while(!feof(inputfile) && *(line = fget_line(line, &line_size, &l_len, inputfile))) {
					if(l_len != 0) {
						deconcount++;
						deconfiles = realloc(deconfiles, deconcount * sizeof(char*));
						if(deconfiles == NULL) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
						deconfiles[deconcount - 1] = strdup(line);
						if(deconfiles[deconcount - 1] == NULL) {
							fprintf(stderr, "OOM\n");
							exit(-1);
						}
					}
				}
				fclose(inputfile);
			}
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
			args++;
			if(args < argc) {
				if(strcmp(argv[args], "-") == 0) {
					prefix_len = 0;
					prefix = 0;
				} else {
					prefix_len = strlen(argv[args]);
					prefix = 0;
					for(i = 0; i < prefix_len; i++) {
						prefix = (prefix << 2) | to2Bit[argv[args][i]];
						if(to2Bit[argv[args][i]] > 3) {
							fprintf(stderr, "Invalid prefix.\n");
							exit(-1);
						}
					}
					if(prefix_len == 0) {
						fprintf(stderr, "Invalid prefix.\n");
						exit(-1);
					}
				}
			}
		} else if(strcmp(argv[args], "-ME") == 0) {
			megaDB = 1;
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(-1);
		}
		args++;
	}
	
	/* check for sufficient input */
	if(filecount == 0 && deconcount == 0) {
		fprintf(stderr, "No inputfiles defined.\n");
		helpMessage(-1);
	} else if(outputfilename == 0 && templatefilename == 0) {
		fprintf(stderr, "Output destination not defined.\n");
		helpMessage(-1);
	} else if(outputfilename == 0 && templatefilename != 0) {
		outputfilename = malloc((strlen(templatefilename) + 64));
		if(!outputfilename) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		strcpy(outputfilename, templatefilename);
	}
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	if(megaDB) {
		INITIAL_SIZE = mask;
		INITIAL_SIZE++;
	}
	/* load or allocate DB */
	if(templatefilename != 0) {
		/* load */
		appender = load_DBs(templatefilename, outputfilename);
		
		/* determine params based on loaded DB */
		if(prefix_len == 0) {
			sparse_run = 0;
		} else {
			sparse_run = 1;
		}
		if(mask == templates->size) {
			megaDB = 1;
		} else {
			megaDB = 0;
		}
	} else {
		/* create */
		templates = hashMap_initialize(INITIAL_SIZE);
		if(templates->seq_size == 0) {
			megaDB = 1;
		}
		
		DB_size = 1;
		if(sparse_run) {
			templates->prefix = prefix;
			templates->prefix_len = prefix_len;
			
			template_lengths = malloc(1024 * sizeof(unsigned));;
			template_slengths = malloc(1024 * sizeof(unsigned));
			template_ulengths = malloc(1024 * sizeof(unsigned));
			if(!template_lengths || !template_slengths || !template_ulengths) {
				OOM();
			}
			*template_lengths = kmerindex;
			template_slengths[0] = 1024;
			template_ulengths[0] = 1024;
		} else {
			template_lengths = malloc(1024 * sizeof(unsigned));
			template_slengths = 0;
			template_ulengths = 0;
			if(!template_lengths) {
				OOM();
			}
			template_lengths[0] = 1024;
		}
	}
	
	/* set pointers and dependent global variables */
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	shifterI = sizeof(long unsigned) * sizeof(long unsigned) - (kmerindex << 1);
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	/* function pointers */
	if(sparse_run) {
		update_DB = &updateDBs_sparse;
		updateAnnotsPtr = &updateAnnots_sparse;
		hashMap_add = &hashMap_addKMASparse;
		hashMap_get = &hashMap_getValue;
		addCont = &hashMap_addCont;
	} else {
		update_DB = &updateDBs;
		updateAnnotsPtr = &updateAnnots;
		hashMap_add = &hashMap_addKMA;
		hashMap_get = &hashMap_getValue;
		addCont = &hashMap_addCont;
	}
	if(prefix_len != 0) {
		deConNode_ptr = &deConNode_sparse;
	} else {
		deConNode_ptr = &deConNode;
	}
	if(megaDB) {
		if(prefix_len) {
			hashMap_add = &megaMap_addKMA_sparse;
		} else {
			hashMap_add = &megaMap_addKMA;
		}
		hashMap_get = &megaMap_getValue;
		addCont = &megaMap_addCont;
	}
	if(kmersize <= 16) {
		getKmerP = &getK;
		addKmer = &addK;
	}
	
	/* set homology check */
	if(MinLen > (kmersize + prefix_len + 1)) {
		MinKlen = 2 * (MinLen - kmersize - prefix_len + 1);
		for(i = 0; i < prefix_len; i++) {
			MinKlen /= 4;
		}
	}
	if(homT < 1) {
		QualCheck = &templateCheck;
		foundKmers = malloc(sizeof(struct hashMap_kmers));
		if(!foundKmers) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		foundKmers->size = INITIAL_SIZE;
		foundKmers->table = calloc(foundKmers->size, sizeof(struct hashTable_kmers *));
		Scores = calloc(1024, sizeof(unsigned));
		Scores_tot = calloc(1024, sizeof(unsigned));
		bestTemplates = malloc(1024 * sizeof(unsigned));
		if(!foundKmers->table || !Scores || !Scores_tot || !bestTemplates) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		Scores_tot[0] = 1024;
	} else if(homQ < 1) {
		QualCheck = &queryCheck;
		Scores_tot = calloc(1024, sizeof(unsigned));
		bestTemplates = malloc(1024 * sizeof(unsigned));
		if(!Scores_tot || !bestTemplates) {
			fprintf(stderr, "OOM\n");
			exit(-1);
		}
		Scores_tot[0] = 1024;
	} else {
		QualCheck = &lengthCheck;
	}
	
	/* update DBs */
	if(filecount != 0) {
		makeDB(inputfiles, filecount, outputfilename, appender);	
	}
	
	/* compress db */
	fprintf(stderr, "# Compressing templates\n");
	file_len = strlen(outputfilename);
	strcat(outputfilename, ".comp.b");
	out = fopen(outputfilename, "wb");
	if(!out) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		exit(-1);
	} else if(templates->table != 0) {
		finalDB = compressKMA_DB(out);
	} else {
		finalDB = compressKMA_megaDB(out);
	}
	fclose(out);
	outputfilename[file_len] = 0;
	fprintf(stderr, "# Template database created.\n");
	
	/* decontaminate */
	if(deconcount != 0) {
		/* clear compression values,
		as only these are changed under decontamination */
		free(finalDB->values);
		
		/* get decontamination info */
		fprintf(stderr, "# Adding decontamination information\n");
		mapped_cont = deConDB(finalDB, deconfiles, deconcount, outputfilename);
		fprintf(stderr, "# Contamination information added.\n");
		fprintf(stderr, "# %d kmers mapped to the DB.\n", mapped_cont);
		fprintf(stderr, "# Contamination mapped to %f %% of the DB.\n", 100.0 * mapped_cont / templates->n);
		
		/* compress DB */
		fprintf(stderr, "# Compressing templates\n");
		if((finalDB->size - 1) != mask) {
			compressKMA_deconDB(finalDB);
		} else {
			compressKMA_deconMegaDB(finalDB); /* here */
		}
		
		/* dump DB */
		fprintf(stderr, "# Dumping DB.\n");
		strcat(outputfilename, ".decon.comp.b");
		out = fopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		if(!out) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			exit(-1);
		}
		if((finalDB->size - 1) != mask) {
			hashMapKMA_dump(finalDB, out);
		} else {
			megaMapKMA_dump(finalDB, out);
		}
		fclose(out);
		outputfilename[file_len] = 0;
	}
	
	return 0;
}
