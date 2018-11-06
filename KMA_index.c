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

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

#define HU_LIMIT 65535
#define U_LIMIT 4294967295
#define CHUNK 1048576
#define ENABLE_ZLIB_GZIP 32
#define getNuc(Comp,pos)((Comp[pos >> 5] << ((pos & 31) << 1)) >> 62)
#define SetNuc(seq, nuc, pos)(seq[pos >> 5] |= (nuc << (62-((pos & 31) << 1))))
#define MIN(X, Y) ((X < Y) ? X : Y)
#define MAX(X, Y) ((X < Y) ? Y : X)

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
	long unsigned size;				// size of DB
	long unsigned n;				// k-mers stored
	long unsigned null_index;		// null value
	long unsigned v_index;			// size of values
	unsigned kmersize;				// k
	unsigned prefix_len;			// prefix length
	long unsigned prefix;			// prefix
	unsigned *exist;				// size long
	long unsigned *exist_l;			// size long, big DBs
	unsigned *values;				// compressed values
	short unsigned *values_s;		// compressed values, few templates
	unsigned *key_index;			// Relative
	long unsigned *key_index_l;		// Relative, 16 < k
	unsigned *value_index;			// Relative
	long unsigned *value_index_l;	// Relative, big DBs
};

struct hashTable {
	long unsigned key;
	unsigned *values;
	struct hashTable *next;
};

struct hashMap {
	/* open hash structure */
	unsigned kmersize;			// k
	long unsigned size;			// size of DB
	long unsigned n;			// k-mers stored
	unsigned prefix_len;		// prefix length
	long unsigned prefix;		// prefix
	struct hashTable **table;	// org
	unsigned **values;			// ME
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
	long unsigned v_index;
	unsigned *values;
	struct valuesTable *next;
};

struct valuesHash {
	long unsigned n;
	long unsigned size;
	struct valuesTable **table;
};

struct qseqs {
	int size;
	int len;
	unsigned char *seq;
};

struct FileBuff {
	int bytes;
	int buffSize;
	unsigned char *buffer;
	unsigned char *inBuffer;
	unsigned char *next;
	FILE *file;
	z_stream *strm;
	int z_err;
};


/*
	GLOBAL VARIABLES
*/
int version[3] = {1, 0, 2};
struct hashMap *templates;
struct hashMap_kmers *foundKmers;
int kmersize, kmerindex, DB_size, prefix_len, MinLen, MinKlen, shifter;
unsigned shifterI, prefix_shifter, *Scores, *Scores_tot, *bestTemplates;
unsigned *template_lengths, *template_slengths, *template_ulengths;
long unsigned mask, INITIAL_SIZE, prefix;
double homQ, homT;

/* 
	FUNCTION POINTERS
*/
int (*homcmp)(int, int);
int (*update_DB)(struct hashMap *, struct compDNA *, unsigned);
void (*load_ptr)(struct hashMap*, char*);
int (*deConNode_ptr)(struct compDNA *, struct hashMapKMA *, unsigned **);
int (*addCont)(struct hashMapKMA *, long unsigned, int, unsigned **);
int (*QualCheck)(struct compDNA *);
void (*updateAnnotsPtr)(struct compDNA *, FILE *, FILE *);
int (*hashMap_add)(struct hashMap *, long unsigned, unsigned);
unsigned * (*hashMap_get)(struct hashMap *, long unsigned);
int (*buffFileBuff)(struct FileBuff *);
void (*dumpIndex)(struct compDNA *, FILE *, FILE *);
unsigned * (*updateValuePtr)(unsigned *, unsigned);
long unsigned (*valuesKeyPtr)(unsigned *);
int (*cmpValuesPtr)(unsigned *, unsigned *, unsigned);
unsigned (*valuesSize)(unsigned *);
void (*hashMapKMA_addKey_ptr)(struct hashMapKMA *, long unsigned, long unsigned);
void (*hashMapKMA_addValue_ptr)(struct hashMapKMA *, long unsigned, long unsigned);
void (*hashMapKMA_addExist_ptr)(struct hashMapKMA *, long unsigned, long unsigned);
void (*updateScoreAndTemplate_ptr)(unsigned *, unsigned *, unsigned *);
void (*addUscore_ptr)(unsigned *, unsigned *);
void (*addUniqueValues)(struct hashMap *, long unsigned, unsigned *);
long unsigned (*getExistPtr)(struct hashMapKMA *, long unsigned);
long unsigned (*getKeyPtr)(struct hashMapKMA *, long unsigned);
long unsigned (*getValueIndexPtr)(struct hashMapKMA *, long unsigned);
unsigned * (*getValuePtr)(struct hashMapKMA *, long unsigned);
int (*getSizePtr)(unsigned *, int);

/*
	FUNCTIONS
*/

void ERROR() {
	fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
	exit(errno);
}

void * smalloc(size_t size) {
	
	void *dest;
	
	dest = malloc(size);
	if(!dest) {
		ERROR();
	}
	
	return dest;
}

void * memdup(const void * src, size_t size) {
	
	void *dest;
	
	dest = smalloc(size);
	memcpy(dest, src, size);
	
	return dest;
}

FILE * sfopen(const char *filename, const char *mode) {
	
	FILE *file;
	
	file = fopen(filename, mode);
	if(!file) {
		fprintf(stderr, "File:\t%s\n", filename);
		ERROR();
	}
	
	return file;
}

void sfwrite(const void *src, size_t size, size_t nmemb, FILE *stream) {
	
	char *ptr;
	
	size *= nmemb;
	ptr = (char *)(src);
	do {
		nmemb = fwrite(ptr, 1, size, stream);
		if(nmemb == 0) {
			ERROR();
		}
		size -= nmemb;
		ptr += nmemb;
	} while(size);
	
}

int chomp(char *string) {
	/* remove trailing spaces and newlines */
	int k = strlen(string) - 1;
	/* isspace = ((string[k] >= 9  && string[k] <= 13) || string[k] == 32), in ASCII */
	while (isspace(string[k])) {
		--k;
	}
	++k;
	string[k] = 0;
	return k;
}

int uint_eq(const unsigned *s1, const unsigned *s2, int len) {
	
	if(len == 0) {
		return 1;
	}
	
	while(len--) {
		if(s1[len] != s2[len]) {
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
	for(i = 0; i < len; ++i) {
		if(s1[i] != s2[i]) {
			return 0;
		}
	}
	return 1;
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
		ERROR();
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
		ERROR();
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
		for(j = i; j < end; ++j) {
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

int compDNAref(struct compDNA *compressor, unsigned char *qseq, int seqlen) {
	
	int i, j, pos, end, bias;
	unsigned char *seq;
	
	/* trim leadin N's */
	seq = qseq;
	bias = 0;
	while(*seq == 4) {
		++seq;
		++bias;
	}
	seqlen -= bias;
	/* trim trailing N's */
	--seqlen;
	while(seq[seqlen] == 4) {
		--seqlen;
	}
	++seqlen;
	
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
		for(j = i; j < end; ++j) {
			if(seq[j] == 4) {
				compressor->seq[pos] = (compressor->seq[pos] << 2);
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
	
	return bias;
}

void unCompDNA(struct compDNA *compressor, char *seq) {
	
	int i;
	
	/* get nucs */
	for(i = 0; i < compressor->seqlen; ++i) {
		seq[i] = getNuc(compressor->seq, i);
	}
	
	/* get N's */
	for(i = 1; i <= compressor->N[0]; ++i) {
		seq[compressor->N[i]] = 4;
	}
	
}

long unsigned getKmer(long unsigned *compressor, unsigned cPos) {
	
	unsigned iPos = (cPos & 31) << 1;
	cPos >>= 5;
	
	return (iPos <= shifter) ? ((compressor[cPos] << iPos) >> shifter) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifter);
}

long unsigned getK(long unsigned *compressor, unsigned pos) {
	return pos;
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
	for(i = 0, j = compressor->complen - 1; i < compressor->complen; ++i, --j) {
		compressor_rc->seq[j] = binRev2(~compressor->seq[i]);
	}
	
	/* shift */
	if((compressor->seqlen & 31)) {
		shift = (((compressor->complen << 5) - compressor->seqlen) << 1);
		r_shift = 64 - shift;
		for(i = 0, j = 1; j < compressor->complen; i = j++) {
			compressor_rc->seq[i] = (compressor_rc->seq[i] << shift) | (compressor_rc->seq[j] >> r_shift);
		}
		compressor_rc->seq[i] <<= shift;
	}
	
	/* add N's */
	compressor_rc->N[0] = compressor->N[0];
	r_shift = compressor->seqlen - 1;
	shift = compressor->N[0];
	for(i = 1, j = compressor->N[0]; i <= shift; ++i, --j) {
		compressor_rc->N[i] = r_shift - compressor->N[j];
		//compressor_rc->N[i] = compressor->seqlen - compressor->N[i] - 1;
	}
	
}

void rcComp(struct compDNA *compressor) {
	
	int i, j, shift, r_shift, complen;
	long unsigned carry;
	
	/* reverse and complement*/
	complen = compressor->complen >> 1;
	for(i = 0, j = compressor->complen - 1; i < complen; ++i, --j) {
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
		for(i = 0, j = 1; j < compressor->complen; i = j++) {
			compressor->seq[i] = (compressor->seq[i] << shift) | (compressor->seq[j] >> r_shift);
		}
		compressor->seq[i] <<= shift;
	}
	
	/* Change N's */
	j = compressor->N[0] - 1;
	complen = compressor->N[0] >> 1;
	++compressor->N;
	r_shift = compressor->seqlen - 1;
	for(i = 0; i < complen; ++i, --j) {
		shift = r_shift - compressor->N[i];
		compressor->N[i] = r_shift - compressor->N[j];
		compressor->N[j] = shift;
	}
	--compressor->N;
	if((compressor->N[0] & 1)) {
		++complen;
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
		ERROR();
	}
	
	dest->len = 0;
	dest->size = size;
	dest->seq = malloc(size);
	if(!dest->seq) {
		ERROR();
	}
	
	return dest;
}

void destroyQseqs(struct qseqs *dest) {
	free(dest->seq);
	free(dest);
}

int BuffgzFileBuff(struct FileBuff *dest) {
	
	int status;
	z_stream *strm;
	
	/* check compressed buffer, and load it */
	strm = dest->strm;
	if(strm->avail_in == 0) {
		strm->avail_in = fread(dest->inBuffer, 1, dest->buffSize, dest->file);
		strm->next_in = (unsigned char*) dest->inBuffer;
		if(strm->avail_in == 0) {
			dest->bytes = 0;
			dest->next = dest->buffer;
			return 0;
		}
	}
	
	/* reset uncompressed buffer */
	strm->avail_out = dest->buffSize;
	strm->next_out = (unsigned char*) dest->buffer;
	
	/* uncompress buffer */
	status = inflate(strm, Z_NO_FLUSH);
	dest->z_err = status;
	
	/* concatenated file */
	if(status == Z_STREAM_END && strm->avail_out == dest->buffSize) {
		inflateReset(strm);
		return BuffgzFileBuff(dest);
	}
	
	if(status == Z_OK || status == Z_STREAM_END) {
		dest->bytes = dest->buffSize - strm->avail_out;
		dest->next = dest->buffer;
		if(status == Z_OK && dest->bytes == 0) {
			return BuffgzFileBuff(dest);
		}
	} else {
		dest->bytes = 0;
		dest->next = dest->buffer;
		fprintf(stderr, "Gzip error %d\n", status);
	}
	
	return dest->bytes;
}

void init_gzFile(struct FileBuff *inputfile) {
	
	int status;
	unsigned char *tmp;
	z_stream *strm;
	
	/* set inBuffer, for compressed format */
	if(inputfile->inBuffer) {
		tmp = inputfile->buffer;
		inputfile->buffer = inputfile->inBuffer;
		inputfile->inBuffer = tmp;
	} else {
		inputfile->inBuffer = inputfile->buffer;
		inputfile->buffer = malloc(CHUNK);
		if(!inputfile->buffer) {
			ERROR();
		}
	}
	inputfile->next = inputfile->buffer;
	
	/* set the compressed stream */
	strm = inputfile->strm;
	if(!strm && !(strm = malloc(sizeof(z_stream)))) {
		ERROR();
	}
	strm->zalloc = Z_NULL;
	strm->zfree  = Z_NULL;
	strm->opaque = Z_NULL;
	status = inflateInit2(strm, 15 | ENABLE_ZLIB_GZIP);
	if(status < 0) {
		fprintf(stderr, "Gzip error %d\n", status);
		exit(status);
	}
	strm->next_in = inputfile->inBuffer;
	strm->avail_in = inputfile->bytes;
	inputfile->strm = strm;
	inputfile->z_err = Z_OK;
	
	inputfile->bytes = BuffgzFileBuff(inputfile);
}

struct FileBuff * setFileBuff(int buffSize) {
	
	struct FileBuff *dest;
	
	dest = malloc(sizeof(struct FileBuff));
	if(!dest) {
		ERROR();
	}
	
	dest->file = 0;
	dest->inBuffer = 0;
	dest->strm = 0;
	dest->buffSize = buffSize;
	dest->buffer = malloc(buffSize);
	dest->next = dest->buffer;
	if(!dest->buffer) {
		ERROR();
	}
	
	return dest;
}

void closeFileBuff(struct FileBuff *dest) {
	fclose(dest->file);
	dest->file = 0;
}

void gzcloseFileBuff(struct FileBuff *dest) {
	
	int status;
	if((status = inflateEnd(dest->strm)) != Z_OK) {
		fprintf(stderr, "Gzip error %d\n", status);
	}
	if(dest->z_err != Z_STREAM_END) {
		fprintf(stderr, "Unexpected end of file\n");
	}
	dest->file = 0;
	dest->strm->avail_out = 0;
}

void openFileBuff(struct FileBuff *dest, char *filename, char *mode) {
	
	dest->file = fopen(filename, mode);
	if(!dest->file) {
		ERROR();
	}
}

void destroyFileBuff(struct FileBuff *dest) {
	free(dest->buffer);
	free(dest->inBuffer);
	free(dest->strm);
	free(dest);
}

int buff_FileBuff(struct FileBuff *dest) {
	dest->bytes = fread(dest->buffer, 1, dest->buffSize, dest->file);
	dest->next = dest->buffer;
	return dest->bytes;
}

int openAndDetermine(struct FileBuff *inputfile, char *filename) {
	
	int FASTQ;
	short unsigned *check;
	
	/* determine filetype and open it */
	FASTQ = 0;
	if(*filename == '-' && strcmp(filename + 1, "-") == 0) {
		inputfile->file = stdin;
	} else {
		openFileBuff(inputfile, filename, "rb");
	}
	if(buff_FileBuff(inputfile)) {
		check = (short unsigned *) inputfile->buffer;
		if(*check == 35615) {
			FASTQ = 4;
			init_gzFile(inputfile);
			buffFileBuff = &BuffgzFileBuff;
		} else {
			buffFileBuff = &buff_FileBuff;
		}
	} else {
		inputfile->buffer[0] = 0;
	}
	
	if(inputfile->buffer[0] == '>') { //FASTA
		FASTQ |= 2;
	} else {
		fprintf(stderr, "%s: is not in fasta format\n", filename);
		fclose(inputfile->file);
	}
	
	return FASTQ;
}

int chunkPos(char *seq, int start, int end) {
	
	while(start != end) {
		if(seq[start] == '\n') {
			return start;
		}
		++start;
	}
	return end;
}

int FileBuffgetFsa(struct FileBuff *src, struct qseqs *header, struct qseqs *qseq, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get header */
	seq = header->seq;
	size = header->size;
	while((*seq++ = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = header->size;
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				ERROR();
			}
			seq = header->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(isspace(*--seq)) {
		++size;
	}
	*++seq = 0;
	header->len = header->size - size + 1;
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while(*buff != '>') {
		*seq = trans[*buff++];
		if(((*seq) >> 3) == 0) {
			if(--size == 0) {
				size = qseq->size;
				qseq->size <<= 1;
				qseq->seq = realloc(qseq->seq, qseq->size);
				if(!qseq->seq) {
					ERROR();
				}
				seq = qseq->seq + size;
			} else {
				++seq;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				/* chomp seq */
				while(*--seq == 8) {
					++size;
				}
				*++seq = 0;
				qseq->len = qseq->size - size;
				
				src->bytes = 0;
				src->next = buff;
				return 1;
			}
			buff = src->buffer;
		}
	}
	
	/* chomp seq */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

/*
	HASHMAP FUNCTIONS
*/

struct hashMap * hashMap_initialize(long unsigned size) {
	
	struct hashMap *src;
	
	src = malloc(sizeof(struct hashMap));
	if(!src) {
		ERROR();
	}
	
	src->kmersize = kmersize;
	src->size = size;
	src->n = 0;
	src->prefix_len = 0;
	src->prefix = 0;
	
	if((size - 1) == mask) {
		src->table = 0;
		src->values = calloc(src->size, sizeof(unsigned *));
		if(!src->values) {
			ERROR();
		}
	} else {
		src->values = 0;
		src->table = calloc(src->size, sizeof(struct hashTable *));
		if(!src->table) {
			ERROR();
		}
	}
	
	/* masking */
	--src->size;
	
	return src;
}

struct hashMap * hashMap_load(char *filename) {
	
	return 0;
}

int CP(char *templatefilename, char *outputfilename) {
	
	int bytes, buffSize;
	char *buffer;
	FILE *file_in, *file_out;
	
	if(strcmp(templatefilename, outputfilename) == 0) {
		return 1;
	}
	
	file_in = fopen(templatefilename, "rb");
	if(!file_in) {
		return 2;
	}
	file_out = fopen(outputfilename, "wb");
	if(!file_out) {
		ERROR();
	}
	
	buffSize = 1024 * 1024;
	buffer = malloc(buffSize);
	if(!buffer) {
		ERROR();
	}
	
	while((bytes = fread(buffer, 1, buffSize, file_in))) {
		sfwrite(buffer, 1, bytes, file_out);
	}
	
	fclose(file_in);
	fclose(file_out);
	
	free(buffer);
	
	return 0;
}

int megaMap_addKMA(struct hashMap *templates, long unsigned key, unsigned value) {
	
	unsigned *values;
	
	values = templates->values[key];
	if(values == 0) {
		templates->values[key] = updateValuePtr(0, value);
		++templates->n;
		return 1;
	} else if((values = updateValuePtr(values, value))) {
		templates->values[key] = values;
		return 1;
	}
	
	return 0;
}

unsigned * megaMap_getValue(struct hashMap *templates, long unsigned key) {
	return templates->values[key];
}

void hashMapKMA_addKey(struct hashMapKMA *dest, long unsigned index, long unsigned key) {
	dest->key_index[index] = key;
}

void hashMapKMA_addKeyL(struct hashMapKMA *dest, long unsigned index, long unsigned key) {
	dest->key_index_l[index] = key;
}

void hashMapKMA_addValue(struct hashMapKMA *dest, long unsigned index, long unsigned v_index) {
	dest->value_index[index] = v_index;
}

void hashMapKMA_addValueL(struct hashMapKMA *dest, long unsigned index, long unsigned v_index) {
	dest->value_index_l[index] = v_index;
}

void hashMapKMA_addExist(struct hashMapKMA *dest, long unsigned index, long unsigned relative) {
	dest->exist[index] = relative;
}

void hashMapKMA_addExistL(struct hashMapKMA *dest, long unsigned index, long unsigned relative) {
	dest->exist_l[index] = relative;
}

int hashMap_addCont(struct hashMapKMA *dest, long unsigned key, int value, unsigned **Values) {
	
	unsigned pos, kpos, *values;
	long unsigned kmer;
	
	kpos = key & dest->size;
	pos = getExistPtr(dest, kpos);
	
	if(pos != dest->null_index) {
		kmer = getKeyPtr(dest, pos);
		while(key != kmer) {
			++pos;
			if(kpos != (kmer & dest->size)) {
				return 0;
			}
			kmer = getKeyPtr(dest, pos);
		}
		values = updateValuePtr(Values[getValueIndexPtr(dest, pos)], value);
		if(values) {
			Values[getValueIndexPtr(dest, pos)] = values;
			return 1;
		}
	}
	
	return 0;
}

int megaMap_addCont(struct hashMapKMA *dest, long unsigned index, int value, unsigned **Values) {
	
	long unsigned pos;
	unsigned *values;
	
	if((pos = getExistPtr(dest, index)) != dest->n) {
		values = updateValuePtr(Values[pos], value);
		if(values) {
			Values[pos] = values;
			return 1;
		}
	}
	
	return 0;
}

void hashMap2megaMap(struct hashMap *templates, struct hashTable *table) {
	
	struct hashTable *node, *next;
	
	templates->table = 0;
	templates->values = calloc(templates->size, sizeof(unsigned *));
	if(!templates->values) {
		ERROR();
	}
	--templates->size;
	
	/* convert table */
	for(node = table; node != 0; node = next) {
		next = node->next;
		
		/* move values */
		templates->values[node->key] = node->values;
		
		/* clean */
		free(node);
	}
	
	/* clean */
	templates->table = 0;
	
	/* set pointers */
	hashMap_add = &megaMap_addKMA;
	hashMap_get = &megaMap_getValue;
	addCont = &megaMap_addCont;
}

unsigned * updateValue(unsigned *values, unsigned value) {
	
	if(!values) {
		values = smalloc(2 * sizeof(unsigned));
		values[0] = 1;
		values[1] = value;
	} else if(values[*values] == value) {
		return 0;
	} else {
		values[0]++;
		values = realloc(values, (values[0] + 1) * sizeof(unsigned));
		if(!values) {
			ERROR();
		}
		values[*values] = value;
	}
	
	return values;
}

unsigned * updateShortValue(unsigned *valuesOrg, unsigned value) {
	
	short unsigned *values;
	
	values = (short unsigned *)(valuesOrg);
	if(!values) {
		values = malloc(2 * sizeof(short unsigned));
		values[0] = 1;
		values[1] = value;
	} else if(values[*values] == value) {
		return 0;
	} else {
		values[0]++;
		values = realloc(values, (values[0] + 1) * sizeof(short unsigned));
		if(!values) {
			ERROR();
		}
		values[*values] = value;
	}
	valuesOrg = (unsigned *)(values);
	
	return valuesOrg;
}

int hashMap_addKMA(struct hashMap *templates, long unsigned key, unsigned value) {
	
	unsigned *values;
	long unsigned index;
	struct hashTable *node, *next, *table;
	
	index = key & templates->size;
	/* check if key exists */
	for(node = templates->table[index]; node != 0; node = node->next) {
		if(key == node->key) {
			if((values = updateValuePtr(node->values, value))) {
				node->values = values;
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	/* new value check if there is space */
	if(templates->n == templates->size) {
		++templates->size;
		/* link table */
		table = 0;
		index = templates->size;
		while(index--) {
			for(node = templates->table[index]; node != 0; node = next) {
				next = node->next;
				node->next = table;
				table = node;
			}
		}
		free(templates->table);
		
		/* check for megamap */
		templates->size <<= 1;
		if((templates->size - 1) == mask) {
			hashMap2megaMap(templates, table);
			return megaMap_addKMA(templates, key, value);
		}
		
		/* reallocate */
		templates->table = calloc(templates->size, sizeof(struct hashTable));
		if(!templates->table) {
			ERROR();
		}
		--templates->size;
		
		for(node = table; node != 0; node = next) {
			next = node->next;
			index = node->key & templates->size;
			node->next = templates->table[index];
			templates->table[index] = node;
		}
		
		index = key & templates->size;
	}
	
	/* add new value */
	node = malloc(sizeof(struct hashTable));
	if(!node) {
		ERROR();
	}
	/* key */
	node->key = key;
	
	/* value */
	node->values = updateValuePtr(0, value);
	
	/* push it */
	node->next = templates->table[index];
	templates->table[index] = node;
	
	++templates->n;
	
	return 1;
}

unsigned * hashMap_getValue(struct hashMap *templates, long unsigned key) {
	
	struct hashTable *node;
	
	for(node = templates->table[key & templates->size]; node != 0; node = node->next) {
		if(key == node->key) {
			return node->values;
		}
	}
	
	return 0;
}

void hashMap_addUniqueValues(struct hashMap *dest, long unsigned key, unsigned *values) {
	
	long unsigned index;
	struct hashTable *node;
	
	node = smalloc(sizeof(struct hashTable));
	node->key = key;
	node->values = values;
	index = key & dest->size;
	node->next = dest->table[index];
	dest->table[index] = node;
	dest->n++;
}

void megaMap_addUniqueValues(struct hashMap *dest, long unsigned key, unsigned *values) {
	
	dest->values[key] = values;
	dest->n++;
}

/*
	COMPRESSED HASHMAP FUNCTIONS
*/
long unsigned getExist(struct hashMapKMA *dest, long unsigned pos) {
	return dest->exist[pos];
}

long unsigned getExistL(struct hashMapKMA *dest, long unsigned pos) {
	return dest->exist_l[pos];
}

long unsigned getKey(struct hashMapKMA *dest, long unsigned pos) {
	return dest->key_index[pos];
}

long unsigned getKeyL(struct hashMapKMA *dest, long unsigned pos) {
	return dest->key_index_l[pos];
}

long unsigned getValueIndex(struct hashMapKMA *dest, long unsigned pos) {
	return dest->value_index[pos];
}

long unsigned getValueIndexL(struct hashMapKMA *dest, long unsigned pos) {
	return dest->value_index[pos];
}

unsigned * getValue(struct hashMapKMA *dest, long unsigned pos) {
	return (dest->values + pos);
}

unsigned * getValueS(struct hashMapKMA *dest, long unsigned pos) {
	return (unsigned *)(dest->values_s + pos);
}

int getSize(unsigned *values, int pos) {
	return values[pos] * sizeof(unsigned) + sizeof(unsigned);
}

int getSizeS(unsigned *values, int pos) {
	return ((short unsigned *)(values))[pos] * sizeof(short unsigned) + sizeof(short unsigned);
}

struct hashMap * hashMapKMA_openChains(struct hashMapKMA *src) {
	
	long unsigned i, key;
	unsigned *values;
	struct hashMap *dest;
	
	if(mask != src->size) {
		free(src->exist);
	}
	dest = hashMap_initialize(src->size + 1);
	
	if(dest->size == mask) {
		addUniqueValues = &megaMap_addUniqueValues;
	} else {
		addUniqueValues = &hashMap_addUniqueValues;
	}
	
	if(mask != src->size) {
		/* norm */
		i = src->n;
		while(i--) {
			key = getKeyPtr(src, i);
			values = getValuePtr(src, getValueIndexPtr(src, i));
			values = memdup(values, getSizePtr(values, 0));
			addUniqueValues(dest, key, values);
		}
		free(src->key_index);
		free(src->value_index);
	} else {	
		/* mega */
		i = src->size + 1;
		while(i--) {
			if(getExistPtr(src, i) != 1) {
				values = getValuePtr(src, i);
				values = memdup(values, getSizePtr(values, 0));
				
				addUniqueValues(dest, i, values);
			}
		}
		free(src->exist);
	}
	free(src->values);
	free(src);
	
	return dest;
}

unsigned ** hashMapKMA_openValues(struct hashMapKMA *src) {
	
	long unsigned i, index, pos;
	unsigned *values, **Values;
	
	Values = smalloc(src->n * sizeof(unsigned *));
	index = src->n;
	if(mask != (src->size - 1)) {
		/* norm */
		i = src->n;
		while(i--) {
			/* get key and values */
			values = getValuePtr(src, getValueIndexPtr(src, i));
			Values[i] = memdup(values, getSizePtr(values, 0));
			hashMapKMA_addValue_ptr(src, i, i);
		}
	} else {
		/* mega */
		i = src->size;
		getSizePtr = &getSizeS;
		while(i--) {
			if((pos = getExistPtr(src, i)) != 1) {
				/* get key and values */
				values = getValuePtr(src, pos);
				Values[--index] = memdup(values, getSizePtr(values, 0));
				hashMapKMA_addExist_ptr(src, i, index);
			} else {
				hashMapKMA_addExist_ptr(src, i, src->n);
			}
		}
	}
	free(src->values);
	src->values = 0;
	
	return Values;
}

void hashMapKMA_dump(struct hashMapKMA *dest, FILE *out) {
	
	/* dump sizes */
	sfwrite(&DB_size, sizeof(unsigned), 1, out);
	sfwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	sfwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	sfwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	sfwrite(&dest->size, sizeof(long unsigned), 1, out);
	sfwrite(&dest->n, sizeof(long unsigned), 1, out);
	sfwrite(&dest->v_index, sizeof(long unsigned), 1, out);
	sfwrite(&dest->null_index, sizeof(long unsigned), 1, out);
	
	/* dump arrays */
	if(dest->n <= U_LIMIT) {
		sfwrite(dest->exist, sizeof(unsigned), dest->size, out);
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		sfwrite(dest->exist_l, sizeof(long unsigned), dest->size, out);
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	
	if(DB_size < HU_LIMIT) {
		sfwrite(dest->values_s, sizeof(short unsigned), dest->v_index, out);
		getValuePtr = &getValueS;
		getSizePtr = &getSizeS;
	} else {
		sfwrite(dest->values, sizeof(unsigned), dest->v_index, out);
		getValuePtr = &getValue;
		getSizePtr = &getSize;
	}
	
	if(dest->kmersize <= 16) {
		sfwrite(dest->key_index, sizeof(unsigned), dest->n + 1, out);
		getKeyPtr = &getKey;
	} else {
		sfwrite(dest->key_index_l, sizeof(long unsigned), dest->n + 1, out);
		getKeyPtr = &getKeyL;
	}
	
	if(dest->v_index < U_LIMIT) {
		sfwrite(dest->value_index, sizeof(unsigned), dest->n, out);
		hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
		getValueIndexPtr = &getValueIndex;
	} else {
		sfwrite(dest->value_index_l, sizeof(long unsigned), dest->n, out);
		hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
		getValueIndexPtr = &getValueIndexL;
	}
}

void megaMapKMA_dump(struct hashMapKMA *dest, FILE *out) {
	
	/* dump sizes */
	sfwrite(&DB_size, sizeof(unsigned), 1, out);
	sfwrite(&dest->kmersize, sizeof(unsigned), 1, out);
	sfwrite(&dest->prefix_len, sizeof(unsigned), 1, out);
	sfwrite(&dest->prefix, sizeof(long unsigned), 1, out);
	sfwrite(&dest->size, sizeof(long unsigned), 1, out);
	sfwrite(&dest->n, sizeof(long unsigned), 1, out);
	sfwrite(&dest->v_index, sizeof(long unsigned), 1, out);
	sfwrite(&dest->null_index, sizeof(long unsigned), 1, out);
	
	/* dump arrays */
	if(dest->v_index <= U_LIMIT) {
		sfwrite(dest->exist, sizeof(unsigned), dest->size, out);
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		sfwrite(dest->exist_l, sizeof(long unsigned), dest->size, out);
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	
	if(DB_size < HU_LIMIT) {
		sfwrite(dest->values_s, sizeof(short unsigned), dest->v_index, out);
		getValuePtr = &getValueS;
		getSizePtr = &getSizeS;
	} else {
		sfwrite(dest->values, sizeof(unsigned), dest->v_index, out);
		getValuePtr = &getValue;
		getSizePtr = &getSize;
	}
}

void * allocAndLoad(size_t size, size_t n, FILE *inputfile) {
	
	void *dest;
	
	dest = malloc(n * size);
	if(!dest) {
		ERROR();
	}
	
	fread(dest, size, n, inputfile);
	
	return dest;
}

int hashMapKMA_load(struct hashMapKMA *dest, FILE *file) {
	
	long unsigned check, size;
	
	/* load sizes */
	fread(&DB_size, sizeof(unsigned), 1, file);
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	fread(&dest->v_index, sizeof(long unsigned), 1, file);
	fread(&dest->null_index, sizeof(long unsigned), 1, file);
	kmersize = dest->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	
	/* exist */
	size = dest->size;
	if((dest->size - 1) == mask) {
		if(dest->v_index <= U_LIMIT) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
		}
	} else {
		if(dest->n <= U_LIMIT) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
			hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
		}
	}
	dest->exist = smalloc(size);
	check = fread(dest->exist, 1, size, file);
	if(check != size) {
		return 1;
	}
	dest->exist_l = (long unsigned *)(dest->exist);
	
	/* values */
	size = dest->v_index;
	if(DB_size < HU_LIMIT) {
		size *= sizeof(short unsigned);
		getValuePtr = &getValueS;
		getSizePtr = &getSizeS;
	} else {
		size *= sizeof(unsigned);
		getValuePtr = &getValue;
		getSizePtr = &getSize;
	}
	dest->values = smalloc(size);
	check = fread(dest->values, 1, size, file);
	if(check != size) {
		return 1;
	}
	dest->values_s = (short unsigned *)(dest->values);
	
	/* check for megaMap */
	if((dest->size - 1) == mask) {
		return 0;
	}
	
	/* kmers */
	size = dest->n + 1;
	if(dest->kmersize <= 16) {
		size *= sizeof(unsigned);
		getKeyPtr = &getKey;
	} else {
		size *= sizeof(long unsigned);
		getKeyPtr = &getKeyL;
	}
	dest->key_index = smalloc(size);
	check = fread(dest->key_index, 1, size, file);
	if(check != size) {
		return 1;
	}
	dest->key_index_l = (long unsigned *)(dest->key_index);
	
	/* value indexes */
	size = dest->n;
	if(dest->v_index < U_LIMIT) {
		size *= sizeof(unsigned);
		hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
		getValueIndexPtr = &getValueIndex;
	} else {
		size *= sizeof(long unsigned);
		hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
		getValueIndexPtr = &getValueIndexL;
	}
	dest->value_index = smalloc(size);
	check = fread(dest->value_index, 1, size, file);
	if(check != size) {
		return 1;
	}
	dest->value_index_l = (long unsigned *)(dest->value_index);
	
	return 0;
}

/* hashMap indexes */
void hashMap_index_initialize(struct hashMap_index *dest, int len) {
	
	dest->len = len;
	dest->size = len << 1;
	
	dest->index = malloc(dest->size * sizeof(int));
	dest->seq = malloc(((len >> 5) + 1) * sizeof(long unsigned));
	if(!dest->index || !dest->seq) {
		ERROR();
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
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(getKmer(dest->seq, abs(pos) - 1) == key) {
			return pos;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
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
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(pos < 0 && getKmer(dest->seq, (-1) - pos) == key) {
			pos = kmersize - pos + 1;
			score = 0;
			for(i = kmersize; i < q_len && pos < dest->len && getNuc(dest->seq, pos) == qseq[i]; ++i, ++pos) {
				++score;
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
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos < 0 && getKmer(dest->seq, (-1) - pos) == key) {
				pos = kmersize - pos + 1;
				score = 0;
				for(i = kmersize; i < q_len && pos < dest->len && getNuc(dest->seq, pos) == qseq[i]; ++i, ++pos) {
					++score;
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
	++newpos;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(pos > 0) {
			if(getKmerIndex(dest->seq, pos - 1) == key) {
				dest->index[index] = -pos;
				neg = -1;
			}
		} else {
			if(getKmerIndex(dest->seq, -pos - 1) == key) {
				neg = -1;
			}
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos > 0) {
				if(getKmerIndex(dest->seq, pos - 1) == key) {
					dest->index[index] = -pos;
					neg = -1;
				}
			} else {
				if(getKmerIndex(dest->seq, -pos - 1) == key) {
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
		ERROR();
	}
	hashMap_index_initialize(src, len);
	
	fread(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
	fread(src->index, sizeof(int), src->size, index);
	
	return src;
}

void hashMap_index_dump(struct hashMap_index *src, FILE *seq, FILE *index) {
	
	sfwrite(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
	sfwrite(src->index, sizeof(int), src->size, index);
	
}

/* VALUES HASH */

struct valuesHash * initialize_hashValues(long unsigned size) {
	
	struct valuesHash *dest;
	
	dest = malloc(sizeof(struct valuesHash));
	if(!dest) {
		ERROR();
	}
	dest->n = 0;
	dest->size = size;
	
	dest->table = calloc(size, sizeof(struct valuesTable *));
	if(!dest->table) {
		ERROR();
	}
	
	return dest;
}

void valuesHash_destroy(struct valuesHash *src) {
	
	long unsigned i;
	struct valuesTable *node, *next;
	
	i = src->size;
	while(i--) {
		for(node = src->table[i]; node != 0; node = next) {
			next = node->next;
			free(node);
		}
	}
	free(src->table);
	free(src);
}

long unsigned valuesKey(unsigned *values) {
	
	unsigned i;
	long unsigned key;
	
	key = 0;
	for(i = 0; i <= *values; ++i) {
		key = key * DB_size + values[i];
	}
	
	return key;
}

long unsigned huValuesKey(unsigned *valuesOrg) {
	
	unsigned i;
	long unsigned key;
	short unsigned *values;
	
	values = (short unsigned *)(valuesOrg);
	key = 0;
	for(i = 0; i <= *values; ++i) {
		key = key * DB_size + values[i];
	}
	
	return key;
}

unsigned uSize(unsigned *values) {
	return *values + 1;
}

unsigned huSize(unsigned *valuesOrg) {
	short unsigned *values;
	values = (short unsigned *)(valuesOrg);
	return *values + 1;
}

int cmpValues(unsigned *s1, unsigned *s2, unsigned len) {
	
	if(len == 0) {
		return 1;
	}
	
	while(len--) {
		if(s1[len] != s2[len]) {
			return 0;
		}
	}
	
	return 1;
}

int cmpHuValues(unsigned *s1_org, unsigned *s2_org, unsigned len_org) {
	
	short unsigned *s1, *s2, len;
	
	len = len_org;
	--len;
	
	if(len == 0) {
		return 1;
	}
	s1 = (short unsigned *)(s1_org);
	s2 = (short unsigned *)(s2_org);
	
	while(len--) {
		if(s1[len] != s2[len]) {
			return 0;
		}
	}
	
	return 1;
}

long unsigned valuesHash_add(struct valuesHash *src, unsigned *newValues, long unsigned v_index) {
	
	/* return 0 if values are new, 
	else return first index of seen value */
	
	unsigned *values;
	long unsigned index;
	struct valuesTable *node;
	
	/* get index */
	index = valuesKeyPtr(newValues) % src->size;
	
	/* search for values */
	for(node = src->table[index]; node != 0; node = node->next) {
		values = node->values;
		if(*values == *newValues && cmpValuesPtr(values + 1, newValues + 1, newValues[0])) { // Value exists
			return node->v_index;
		}
	}
	
	/* new values */
	++src->n;
	node = malloc(sizeof(struct valuesTable));
	if(!node) {
		ERROR();
	}
	node->v_index = v_index;
	node->values = newValues;
	node->next = src->table[index];
	src->table[index] = node;
	
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
		++dest->n;
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
				++dest->n;
				return 1;
			}
		}
	}
	return -1;
}

void emptyHash(struct hashMap_kmers *dest) {
	
	unsigned i;
	struct hashTable_kmers *node, *next;
	
	for(i = 0; i < dest->size; ++i) {
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
int updateDBs(struct hashMap *templates, struct compDNA *qseq, unsigned template) {
	
	int i, j, end;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	/* set parameters */
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	
	/* iterate sequence */
	for(i = 1, j = 0; i <= qseq->N[0]; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end; ++j) {
			/* update hashMap */
			hashMap_add(templates, getKmer(qseq->seq, j), template);
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	
	return 1;
}

int lengthCheck(struct compDNA *qseq) {
	
	int i, j, end, rc, thisKlen;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	} else if(prefix_len == 0) {
		if((qseq->seqlen - kmersize + 1) * 2 < MinKlen) {
			return 0;
		} else {
			return 1;
		}
	}
	
	thisKlen = MinKlen;
	
	for(rc = 0; rc < 2; ++rc) {
		/* revers complement */
		if(rc) {
			rcComp(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0] && thisKlen != 0; ++i) {
			end = qseq->N[i] - prefix_len - kmersize + 1;
			for(;j < end && thisKlen != 0; ++j) {
				if(getPrefix(qseq->seq, j) == prefix) {
					--thisKlen;
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

void updateScoreAndTemplate(unsigned *Scores_tot, unsigned *bestTemplates, unsigned *values) {
	
	unsigned i;
	
	i = *values + 1;
	while(--i) {
		if(Scores_tot[*++values] == 0) {
			bestTemplates[0]++;
			bestTemplates[*bestTemplates] = *values;
		}
		Scores_tot[*values]++;
	}
}

void updateScoreAndTemplateHU(unsigned *Scores_tot, unsigned *bestTemplates, unsigned *values_org) {
	
	unsigned i;
	short unsigned *values;
	
	values = (short unsigned *)(values_org);
	i = *values + 1;
	while(--i) {
		if(Scores_tot[*++values] == 0) {
			bestTemplates[0]++;
			bestTemplates[*bestTemplates] = *values;
		}
		Scores_tot[*values]++;
	}
}

void addUscore(unsigned *Scores, unsigned *values) {
	
	unsigned i;
	
	i = *values + 1;
	while(--i) {
		Scores[*++values]++;
	}
}

void addUscoreHU(unsigned *Scores, unsigned *values_org) {
	
	unsigned i;
	short unsigned *values;
	
	values = (short unsigned *)(values_org);
	i = *values + 1;
	while(--i) {
		Scores[*++values]++;
	}
}

int queryCheck(struct compDNA *qseq) {
	
	unsigned i, j, end, rc, thisKlen, *values;
	double bestQ, thisQ;
	long unsigned prefix;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	} else if(prefix_len == 0 && templates->prefix != 0) {
		prefix = 0;
	} else {
		prefix = templates->prefix;
	}
	
	thisKlen = 0;
	
	/* realloc */
	if(DB_size >= *Scores_tot) {
		free(Scores_tot);
		Scores_tot = calloc(2 * DB_size, sizeof(unsigned));
		free(bestTemplates);
		bestTemplates = malloc(2 * DB_size * sizeof(unsigned));
		if(!Scores_tot || !bestTemplates) {
			ERROR();
		}
		Scores_tot[0] = 2 * DB_size;
	}
	
	/* get scores */
	bestTemplates[0] = 0;
	for(rc = 0; rc < 2; ++rc) {
		/* revers complement */
		if(rc) {
			rcComp(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - prefix_len - kmersize + 1;
			for(;j < end; ++j) {
				if(getPrefix(qseq->seq, j) == prefix) {
					++thisKlen;
					if((values = hashMap_get(templates, getKmer(qseq->seq, j + prefix_len)))) {
						updateScoreAndTemplate_ptr(Scores_tot, bestTemplates, values);
					}
				}
			}
			j = qseq->N[i] + 1;
		}
		qseq->N[0]--;
	}
	
	/* get query cov */
	bestQ = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
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
	
	unsigned i, j, end, rc, thisKlen, *values;
	double bestQ, thisQ, bestT, thisT;
	long unsigned key, prefix;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	} else if(prefix_len == 0 && templates->prefix != 0) {
		prefix = 0;
	} else {
		prefix = templates->prefix;
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
			ERROR();
		}
		Scores_tot[0] = 2 * DB_size;
	}
	
	/* get scores */
	bestTemplates[0] = 0;
	for(rc = 0; rc < 2; ++rc) {
		/* revers complement */
		if(rc) {
			rcComp(qseq);
		}
		
		/* iterate seq */
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1, j = 0; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - prefix_len - kmersize + 1;
			for(;j < end; ++j) {
				if(getPrefix(qseq->seq, j) == prefix) {
					++thisKlen;
					key = getKmer(qseq->seq, j + prefix_len);
					if((values = hashMap_get(templates, key))) {
						updateScoreAndTemplate_ptr(Scores_tot, bestTemplates, values);
						if(hashMap_CountKmer(foundKmers, key)) {
							addUscore_ptr(Scores, values);
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
	for(i = 1; i <= *bestTemplates; ++i) {
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
				ERROR();
			}
		}
		return 1;
	} else {
		return 0;
	}
}

int updateDBs_sparse(struct hashMap *templates, struct compDNA *qseq, unsigned template) {
	
	int i, j, end, rc;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	/* test homology and length */
	if(QualCheck(qseq)) {
		template_slengths[template] = 0;
		template_ulengths[template] = 0;
		for(rc = 0; rc < 2; ++rc) {
			/* revers complement */
			if(rc) {
				rcComp(qseq);
			}
			
			/* iterate seq */
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			j = 0;
			if(prefix_len) {
				for(i = 1; i <= qseq->N[0]; ++i) {
					end = qseq->N[i] - prefix_len - kmersize + 1;
					for(;j < end; ++j) {
						if(getPrefix(qseq->seq, j) == prefix) {
							/* add kmer */
							if(hashMap_add(templates, getKmer(qseq->seq, j + prefix_len), template)) {
								template_ulengths[template]++;
							}
							template_slengths[template]++;
						}
					}
					j = qseq->N[i] + 1;
				}
				qseq->N[0]--;
			} else {
				for(i = 1; i <= qseq->N[0]; ++i) {
					end = qseq->N[i] - kmersize + 1;
					for(;j < end; ++j) {
						/* add kmer */
						if(hashMap_add(templates, getKmer(qseq->seq, j), template)) {
							template_ulengths[template]++;
						}
						template_slengths[template]++;
					}
					j = qseq->N[i] + 1;
				}
				qseq->N[0]--;
			}
		}
		return 1;
	}
	
	return 0;
}

struct hashMapKMA * compressKMA_DB(FILE *out) {
	
	long unsigned i, j, check;
	long unsigned index, t_index, v_index, new_index, null_index;
	unsigned *values;
	short unsigned *values_s;
	struct hashMapKMA *finalDB;
	struct valuesHash *shmValues;
	struct valuesTable *node, *next, *table;
	struct hashTable *node_t, *next_t, *table_t;
	
	/* convert templates to linked list */
	table_t = 0;
	i = templates->size + 1;
	while(i--) {
		for(node_t = templates->table[i]; node_t != 0; node_t = next_t) {
			next_t = node_t->next;
			node_t->next = table_t;
			table_t = node_t;
		}
	}
	free(templates->table);
	templates->table = 0;
	
	/* prepare final DB */
	check = 0;
	check = ~check;
	check >>= 32;
	fprintf(stderr, "# Preparing compressed DB.\n");
	finalDB = malloc(sizeof(struct hashMapKMA));
	if(!finalDB) {
		ERROR();
	}
	/* Fill in known values */
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->prefix_len = prefix_len;
	finalDB->prefix = prefix;
	finalDB->kmersize = kmersize;
	
	/* allocate existence */
	if(finalDB->n <= check) {
		finalDB->exist = smalloc(finalDB->size * sizeof(unsigned));
		finalDB->exist_l = 0;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		finalDB->exist = 0;
		finalDB->exist_l = smalloc(finalDB->size * sizeof(long unsigned));
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	
	if(kmersize <= 16) {
		finalDB->key_index = smalloc((finalDB->n + 1) * sizeof(unsigned));
		finalDB->key_index_l = 0;
		hashMapKMA_addKey_ptr = &hashMapKMA_addKey;
	} else {
		finalDB->key_index = 0;
		finalDB->key_index_l = smalloc((finalDB->n + 1) * sizeof(long unsigned));
		hashMapKMA_addKey_ptr = &hashMapKMA_addKeyL;
	}
	finalDB->value_index = smalloc(finalDB->n * sizeof(unsigned));
	
	null_index = finalDB->n;
	finalDB->null_index = null_index;
	/* fill with null_indexes */
	i = finalDB->size;
	while(i--) {
		hashMapKMA_addExist_ptr(finalDB, i, null_index);
	}
	
	/* get relative indexes */
	fprintf(stderr, "# Calculating relative indexes.\n");
	hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
	node_t = table_t;
	--finalDB->size;
	shmValues = initialize_hashValues(null_index);
	t_index = 0;
	v_index = 0;
	while(node_t != 0) {
		/* get index */
		index = (node_t->key & finalDB->size);
		hashMapKMA_addExist_ptr(finalDB, index, t_index);
		/* mv chain */
		while(node_t != 0 && (node_t->key & finalDB->size) == index) {
			next_t = node_t->next;
			
			/* add kmer */
			hashMapKMA_addKey_ptr(finalDB, t_index, node_t->key);
			
			/* the actual value index */
			new_index = valuesHash_add(shmValues, node_t->values, v_index);
			
			if(new_index == v_index) {
				v_index += valuesSize(node_t->values);
				if(check <= v_index) {
					fprintf(stderr, "# Compression overflow.\n");
					check = 0;
					check = ~check;
					hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
					getValueIndexPtr = &getValueIndexL;
					finalDB->value_index_l = realloc(finalDB->value_index, finalDB->n * sizeof(long unsigned));
					if(!finalDB->value_index_l) {
						ERROR();
					}
					finalDB->value_index = (unsigned *)(finalDB->value_index_l);
					j = finalDB->n;
					while(j--) {
						finalDB->value_index_l[j] = finalDB->value_index[j];
					}
					finalDB->value_index = 0;
				}
			} else {
				/* values were duplicated, clean up */
				free(node_t->values);
			}
			
			hashMapKMA_addValue_ptr(finalDB, t_index, new_index);
			++t_index;
			
			/* clean */
			free(node_t);
			node_t = next_t;
		}
	}
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	
	if(DB_size < HU_LIMIT) {
		finalDB->values = 0;
		finalDB->values_s = calloc(v_index, sizeof(short unsigned));
		if(!finalDB->values_s) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
		}
	} else {
		finalDB->values = calloc(v_index, sizeof(unsigned));
		finalDB->values_s = 0;
		if(!finalDB->values) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
		}
	}
	
	/* add terminating key */
	i = 0;
	if(finalDB->kmersize <= 16) { 
		j = finalDB->key_index[finalDB->n - 1] & finalDB->size;
		while(j == (finalDB->key_index[i] & finalDB->size)) {
			++i;
		}
		finalDB->key_index[finalDB->n] = finalDB->key_index[i];
	} else {
		j = finalDB->key_index_l[finalDB->n - 1] & finalDB->size;
		while(j == (finalDB->key_index_l[i] & finalDB->size)) {
			++i;
		}
		finalDB->key_index_l[finalDB->n] = finalDB->key_index_l[i];
	}
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");	
	++finalDB->size;
	hashMapKMA_dump(finalDB, out);
	
	return finalDB;
}

struct hashMapKMA * compressKMA_megaDB(FILE *out) {
	
	long unsigned i, j, v_index, new_index, null_index;
	unsigned check, *values;
	short unsigned *values_s;
	struct hashMapKMA *finalDB;
	struct valuesHash *shmValues;
	struct valuesTable *node, *next, *table;
	
	/* Fill in known values */
	hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	check = 0;
	check = ~check;
	finalDB = malloc(sizeof(struct hashMapKMA));
	if(!finalDB) {
		ERROR();
	}
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->prefix_len = prefix_len;
	finalDB->prefix = prefix;
	finalDB->kmersize = kmersize;
	
	/* allocate existence */
	finalDB->exist = malloc(finalDB->size * sizeof(unsigned));
	finalDB->exist_l = 0;
	finalDB->key_index = 0;
	finalDB->value_index = 0;
	if(!finalDB->exist) {
		ERROR();
	}
	
	/* get relative indexes */
	fprintf(stderr, "# Calculating relative indexes.\n");
	null_index = finalDB->n;
	v_index = 0;
	while(templates->values[v_index] != 0) {
		finalDB->exist[v_index] = v_index;
		++v_index;
	}
	for(i = v_index; i != finalDB->size; ++i) {
		if(templates->values[i]) {
			finalDB->exist[i] = v_index;
			templates->values[v_index] = templates->values[i];
			templates->values[i] = 0;
			++v_index;
		} else {
			finalDB->exist[i] = null_index;
		}
	}
	templates->values = realloc(templates->values, templates->n * sizeof(unsigned *));
	if(!templates->values) {
		ERROR();
	}
	
	/* get compressed indexes */
	fprintf(stderr, "# Compressing indexes.\n");
	v_index = 0;
	shmValues = initialize_hashValues(null_index);
	i = finalDB->size;
	j = 0;
	while(i--) {
		if(finalDB->exist[i] != null_index) {
			values = templates->values[finalDB->exist[i]];
			
			/* the actual index */
			new_index = valuesHash_add(shmValues, values, v_index);
			
			
			if(new_index == v_index) {
				v_index += valuesSize(values);
				if(check < v_index) {
					fprintf(stderr, "# Compression overflow.\n");
					j = 1;
					hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
					break;
				}
			} else {
				/* values were duplicated, clean up */
				free(values);
				templates->values[finalDB->exist[i]] = 0;
			}
			/* update to new index */
			finalDB->exist[i] = new_index;
		} else {
			finalDB->exist[i] = 1;
		}
	}
	if(j) {
		fprintf(stderr, "# Bypassing overflow.\n");
		finalDB->exist_l = realloc(finalDB->exist, finalDB->size * sizeof(long unsigned));
		if(!finalDB->exist_l) {
			ERROR();
		}
		finalDB->exist = (unsigned *)(finalDB->exist_l);
		j = finalDB->size;
		while(j--) {
			finalDB->exist_l[j] = finalDB->exist[j];
		}
		finalDB->exist = 0;
		finalDB->exist_l[i] = new_index;
		
		while(i--) {
			if(finalDB->exist_l[i] != null_index) {
				values = templates->values[finalDB->exist_l[i]];
				
				/* the actual index */
				new_index = valuesHash_add(shmValues, values, v_index);
				
				
				if(new_index == v_index) {
					v_index += valuesSize(values);
				} else {
					/* values were duplicated, clean up */
					free(values);
					templates->values[finalDB->exist_l[i]] = 0;
				}
				/* update to new index */
				finalDB->exist_l[i] = new_index;
				
			} else {
				finalDB->exist_l[i] = 1;
			}
		}
		fprintf(stderr, "# Overflow bypassed.\n");
	}
	free(templates->values);
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->null_index = 1;
	
	if(DB_size < HU_LIMIT) {
		finalDB->values = 0;
		finalDB->values_s = calloc(v_index, sizeof(short unsigned));
		if(!finalDB->values_s) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
		}
	} else {
		finalDB->values = calloc(v_index, sizeof(unsigned));
		finalDB->values_s = 0;
		if(!finalDB->values) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
		}
	}
	
	/* dump final DB */
	fprintf(stderr, "# Dumping compressed DB\n");
	megaMapKMA_dump(finalDB, out);
	
	return finalDB;
}

void compressKMA_deconDB(struct hashMapKMA *finalDB, unsigned **Values) {
	
	long unsigned i, j, v_index, new_index, check;
	unsigned *values;
	short unsigned *values_s;
	struct valuesHash *shmValues;
	struct valuesTable *node, *next, *table;
	
	/* prepare final DB */
	check = 0;
	check = ~check;
	if(finalDB->v_index < U_LIMIT) {
		check >>= 32;
		hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
		getValueIndexPtr = &getValueIndex;
	} else {
		hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
		getValueIndexPtr = &getValueIndexL;
	}
	i = finalDB->n;
	shmValues = initialize_hashValues(finalDB->n);
	v_index = 0;
	while(i--) {
		/* the actual value index */
		values = Values[i];
		new_index = valuesHash_add(shmValues, values, v_index);
		
		if(new_index == v_index) {
			v_index += valuesSize(values);
			if(check <= v_index) {
				fprintf(stderr, "# Compression overflow.\n");
				check = 0;
				check = ~check;
				hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
				getValueIndexPtr = &getValueIndexL;
				finalDB->value_index_l = realloc(finalDB->value_index, finalDB->n * sizeof(long unsigned));
				if(!finalDB->value_index_l) {
					ERROR();
				}
				finalDB->value_index = (unsigned *)(finalDB->value_index_l);
				j = finalDB->n;
				while(j--) {
					finalDB->value_index_l[j] = finalDB->value_index[j];
				}
				finalDB->value_index = 0;
			}
		} else {
			free(values);
		}
		
		hashMapKMA_addValue_ptr(finalDB, i, new_index);
	}
	free(Values);
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	if(DB_size < HU_LIMIT) {
		finalDB->values = 0;
		finalDB->values_s = calloc(v_index, sizeof(short unsigned));
		if(!finalDB->values_s) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
		}
	} else {
		finalDB->values = calloc(v_index, sizeof(unsigned));
		finalDB->values_s = 0;
		if(!finalDB->values) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
		}
	}
}

void compressKMA_deconMegaDB(struct hashMapKMA *finalDB, unsigned **Values) {
	
	long unsigned i, j, v_index, new_index, pos, check;
	unsigned *values;
	short unsigned *values_s;
	struct valuesHash *shmValues;
	struct valuesTable *node, *next, *table;
	
	fprintf(stderr, "# Compressing indexes.\n");
	check = 0;
	check = ~check;
	if(finalDB->v_index < U_LIMIT) {
		check >>= 32;
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	i = finalDB->size;
	shmValues = initialize_hashValues(finalDB->n);
	v_index = 0;
	while(i--) {
		if((pos = getExistPtr(finalDB, i)) != finalDB->n) {
			values = Values[pos];
			new_index = valuesHash_add(shmValues, values, v_index);
			
			if(new_index == v_index) {
				v_index += valuesSize(values);
				if(check <= v_index) {
					fprintf(stderr, "# Compression overflow.\n");
					check = 0;
					check = ~check;
					getExistPtr = &getExistL;
					hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
					finalDB->exist_l = realloc(finalDB->exist, finalDB->size * sizeof(long unsigned));
					if(!finalDB->value_index_l) {
						ERROR();
					}
					finalDB->exist = (unsigned *)(finalDB->exist_l);
					j = finalDB->size;
					while(j--) {
						finalDB->exist_l[j] = finalDB->exist[j];
					}
					finalDB->exist = 0;
				}
			} else {
				free(values);
			}
			hashMapKMA_addExist_ptr(finalDB, i, new_index);
		} else {
			hashMapKMA_addExist_ptr(finalDB, i, 1);
		}
	}
	free(Values);
	/* convert valuesHash to a linked list */
	table = 0;
	i = shmValues->size;
	while(i--) {
		for(node = shmValues->table[i]; node != 0; node = next) {
			next = node->next;
			node->next = table;
			table = node;
		}
	}
	free(shmValues->table);
	
	/* make compressed values */
	fprintf(stderr, "# Finalizing indexes.\n");
	finalDB->v_index = v_index;
	finalDB->null_index = 1;
	if(DB_size < HU_LIMIT) {
		finalDB->values = 0;
		finalDB->values_s = calloc(v_index, sizeof(short unsigned));
		if(!finalDB->values_s) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values_s = (short unsigned *)(node->values);
			for(i = node->v_index, j = 0; j <= *values_s; ++i, ++j) {
				finalDB->values_s[i] = values_s[j];
			}
			free(values_s);
			free(node);
		}
	} else {
		finalDB->values = calloc(v_index, sizeof(unsigned));
		finalDB->values_s = 0;
		if(!finalDB->values) {
			ERROR();
		}
		/* move values */
		for(node = table; node != 0; node = next) {
			next = node->next;
			values = node->values;
			for(i = node->v_index, j = 0; j <= *values; ++i, ++j) {
				finalDB->values[i] = values[j];
			}
			free(values);
			free(node);
		}
	}
}

void updateAnnots(struct compDNA *qseq, FILE *seq_out, FILE *index_out) {
	
	/* Dump annots */
	dumpIndex(qseq, seq_out, index_out);
	
	template_lengths[DB_size] = qseq->seqlen;
	++DB_size;
	if(DB_size >= template_lengths[0]) {
		template_lengths[0] *= 2;
		template_lengths = realloc(template_lengths, template_lengths[0] * sizeof(unsigned));
		if(!template_lengths) {
			ERROR();
		}
	}
}

void updateAnnots_sparse(struct compDNA *qseq, FILE *seq_out, FILE *index_out) {
	
	/* Dump annots */
	dumpIndex(qseq, seq_out, index_out);
	
	template_lengths[DB_size] = qseq->seqlen;
	++DB_size;
	if(DB_size >= template_ulengths[0]) {
		template_ulengths[0] *= 2;
		template_slengths = realloc(template_slengths, template_ulengths[0] * sizeof(unsigned));
		template_ulengths = realloc(template_ulengths, template_ulengths[0] * sizeof(unsigned));
		template_lengths = realloc(template_lengths, template_ulengths[0] * sizeof(unsigned));
		if(!template_lengths || !template_slengths || !template_ulengths) {
			ERROR();
		}
	}
	
}

int deConNode(struct compDNA *qseq, struct hashMapKMA *finalDB, unsigned **Values) {
	
	int i, j, end, mapped_cont;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	mapped_cont = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	j = 0;
	for(i = 1; i <= qseq->N[0]; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end; ++j) {
			mapped_cont += addCont(finalDB, getKmer(qseq->seq, j), DB_size, Values);
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	return mapped_cont;
}

int deConNode_sparse(struct compDNA *qseq, struct hashMapKMA *finalDB, unsigned **Values) {
	
	int i, j, end, mapped_cont;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	}
	
	mapped_cont = 0;
	qseq->N[0]++;
	qseq->N[qseq->N[0]] = qseq->seqlen;
	j = 0;
	for(i = 1; i <= qseq->N[0]; ++i) {
		end = qseq->N[i] - prefix_len - kmersize + 1;
		for(;j < end; ++j) {
			if(getPrefix(qseq->seq, j) == prefix) {
				mapped_cont += addCont(finalDB, getKmer(qseq->seq, j + prefix_len), DB_size, Values);
			}
		}
		j = qseq->N[i] + 1;
	}
	qseq->N[0]--;
	return mapped_cont;
}

unsigned deConDB(struct hashMapKMA *finalDB, char **inputfiles, int fileCount, char *trans, unsigned **Values) {
	
	int FASTQ;
	unsigned fileCounter, mapped_cont;
	char *filename;
	struct qseqs *header, *qseq;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	
	/* allocate */
	compressor = malloc(sizeof(struct compDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* set variables */
	mapped_cont = 0;
	--finalDB->size;
	
	/* iterate inputfiles */
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		/* open file */
		filename = inputfiles[fileCounter];
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
			/* parse the file */
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
				fprintf(stderr, "# Decon:\t%s\n", header->seq + 1);
				if(qseq->len > kmersize) {
					/* compress DNA */
					if(qseq->len >= compressor->size) {
						freeComp(compressor);
						allocComp(compressor, qseq->len);
					}
					compDNAref(compressor, qseq->seq, qseq->len);
					
					/* Add contamination */
					mapped_cont += deConNode_ptr(compressor, finalDB, Values);
					/* rc */
					rcComp(compressor);
					mapped_cont += deConNode_ptr(compressor, finalDB, Values);
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
	++finalDB->size;
	
	/* clean */
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


void makeIndexing(struct compDNA *compressor, FILE *seq_out, FILE *index_out) {
	
	int i, j, end;
	struct hashMap_index *template_index;
	
	/* allocate index */
	template_index = malloc(sizeof(struct hashMap_index));
	if(!template_index) {
		ERROR();
	}
	template_index->len = compressor->seqlen;
	template_index->size = compressor->seqlen << 1;
	template_index->index = calloc(template_index->size, sizeof(int));
	if(!template_index->index) {
		ERROR();
	}
	
	/* load index */
	template_index->seq = compressor->seq;
	compressor->N[0]++;
	compressor->N[compressor->N[0]] = compressor->seqlen + 1;
	j = 0;
	for(i = 1; i <= compressor->N[0]; ++i) {
		end = compressor->N[i] - kmerindex;
		for(;j < end; ++j) {
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

void dumpSeq(struct compDNA *qseq, FILE *seq_out, FILE *index_out) {
	sfwrite(qseq->seq, sizeof(long unsigned), (qseq->seqlen >> 5) + 1, seq_out);
}

struct hashMapKMA * load_DBs(char *templatefilename, char *outputfilename) {
	
	int file_len, out_len;
	FILE *infile;
	struct hashMapKMA *finalDB;
	
	file_len = strlen(templatefilename);
	out_len = strlen(outputfilename);
	
	/* load hash */
	finalDB = smalloc(sizeof(struct hashMapKMA));
	strcat(templatefilename, ".comp.b");
	infile = sfopen(templatefilename, "rb");
	if(hashMapKMA_load(finalDB, infile)) {
		fprintf(stderr, "Wrong format of DB\n");
		exit(1);
	} else {
		finalDB->size--;
	}
	templatefilename[file_len] = 0;
	fclose(infile);
	
	prefix = finalDB->prefix;
	prefix_len = finalDB->prefix_len;
	
	/* load lengths */
	strcat(templatefilename, ".length.b");
	infile = fopen(templatefilename, "rb");
	if(!infile) {
		ERROR();
	}
	templatefilename[file_len] = 0;
	fread(&DB_size, sizeof(unsigned), 1, infile);
	if(prefix_len) {
		template_lengths = smalloc((DB_size << 1) * sizeof(unsigned));
		template_slengths = smalloc((DB_size << 1) * sizeof(unsigned));
		template_ulengths = smalloc((DB_size << 1) * sizeof(unsigned));
		fread(template_slengths, sizeof(unsigned), DB_size, infile);
		fread(template_ulengths, sizeof(unsigned), DB_size, infile);
		fread(template_lengths, sizeof(unsigned), DB_size, infile);
		kmerindex = *template_lengths;
		template_ulengths[0] = DB_size << 1;
		template_slengths[0] = DB_size << 1;
	} else {
		template_lengths = smalloc((DB_size << 1) * sizeof(unsigned));
		template_slengths = 0;
		template_ulengths = 0;
		fread(template_lengths, sizeof(unsigned), DB_size, infile);
		kmerindex = *template_lengths;
		template_lengths[0] = DB_size << 1;
	}
	fclose(infile);
	
	/* cp name, seq and index */
	strcat(templatefilename, ".name");
	strcat(outputfilename, ".name");
	CP(templatefilename, outputfilename);
	templatefilename[file_len] = 0;
	outputfilename[out_len] = 0;
	
	strcat(templatefilename, ".seq.b");
	strcat(outputfilename, ".seq.b");
	CP(templatefilename, outputfilename);
	templatefilename[file_len] = 0;
	outputfilename[out_len] = 0;
	
	strcat(templatefilename, ".index.b");
	if(dumpIndex == &makeIndexing) {
		strcat(outputfilename, ".index.b");
		if(CP(templatefilename, outputfilename) == 2) {
			dumpIndex = &dumpSeq;
		}
		outputfilename[out_len] = 0;
	} else if((infile = fopen(templatefilename, "rb"))) {
		fclose(infile);
		remove(templatefilename);
	}
	templatefilename[file_len] = 0;
	
	return finalDB;
}

unsigned * HU2U(unsigned *values) {
	
	int i;
	short unsigned *hu_values;
	
	hu_values = (short unsigned *)(values);
	values = realloc(values, (hu_values[0] + 1) * sizeof(unsigned));
	if(!values) {
		ERROR();
	} else {
		hu_values = (short unsigned *)(values);
	}
	
	i = *hu_values + 1;
	while(i--) {
		values[i] = hu_values[i];
	}
	
	return values;
}

void convertToU(struct hashMap *templates) {
	
	long unsigned index;
	struct hashTable *node;
	
	/* convert values */
	index = templates->size + 1;
	if(templates->table) {
		while(index--) {
			for(node = templates->table[index]; node != 0; node = node->next) {
				node->values = HU2U(node->values);
			}
		}
	} else {
		while(index--) {
			if(templates->values[index]) {
				templates->values[index] = HU2U(templates->values[index]);
			}
		}	
	}
	
	/* set pointers */
	updateValuePtr = &updateValue;
	valuesKeyPtr = &valuesKey;
	cmpValuesPtr = &cmpValues;
	valuesSize = &uSize;
	updateScoreAndTemplate_ptr = &updateScoreAndTemplate;
	addUscore_ptr = &addUscore;
}

void makeDB(char **inputfiles, int fileCount, char *outputfilename, int appender, char *trans) {
	
	int fileCounter, file_len, bias, FASTQ;
	char *filename;
	unsigned char *seq;
	FILE *index_out, *seq_out, *length_out, *name_out;
	struct qseqs *header, *qseq;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	
	/* allocate */
	compressor = malloc(sizeof(struct compDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(1024);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* open files */
	file_len = strlen(outputfilename);
	strcat(outputfilename, 	".length.b");
	length_out = fopen(outputfilename, "wb");
	outputfilename[file_len] = 0;
	if(appender) {
		strcat(outputfilename, 	".name");
		name_out = fopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
	} else {
		strcat(outputfilename, 	".name");
		name_out = fopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
	}
	if(!length_out || !name_out) {
		ERROR();
	}
	if(appender) {
		strcat(outputfilename, ".seq.b");
		seq_out = fopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
		if(!seq_out) {
			ERROR();
		}
	} else {
		strcat(outputfilename, ".seq.b");
		seq_out = fopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		if(!seq_out) {
			ERROR();
		}
	}
	
	if(dumpIndex == &makeIndexing) {
		if(appender) {
			strcat(outputfilename, ".index.b");
			index_out = fopen(outputfilename, "ab");
			outputfilename[file_len] = 0;
			if(!index_out) {
				ERROR();
			}
		} else {
			strcat(outputfilename, ".index.b");
			index_out = fopen(outputfilename, "wb");
			outputfilename[file_len] = 0;
			if(!index_out) {
				ERROR();
			}
			sfwrite(&kmerindex, sizeof(int), 1, index_out);
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
				if(qseq->len > MinLen && update_DB(templates, compressor, DB_size)) {
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
					updateAnnotsPtr(compressor, seq_out, index_out);
					fprintf(stderr, "# Added:\t%s\n", header->seq + 1);
					
					if(DB_size == HU_LIMIT) {
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
	sfwrite(&DB_size, sizeof(int), 1, length_out);
	if(template_ulengths != 0) {
		template_ulengths[0] = 0;
		template_slengths[0] = 0;
		sfwrite(template_lengths, sizeof(unsigned), DB_size, length_out);
		sfwrite(template_slengths, sizeof(unsigned), DB_size, length_out);
		sfwrite(template_ulengths, sizeof(unsigned), DB_size, length_out);
	} else {
		template_lengths[0] = kmerindex;
		sfwrite(template_lengths, sizeof(unsigned), DB_size, length_out);
	}
	if(index_out) {
		fclose(index_out);
	}
	fclose(seq_out);
	fclose(length_out);
	fclose(name_out);
	
	fprintf(stderr, "# Templates key-value pairs:\t%lu.\n", templates->n);// / 1048576);
	
	/* clean */
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
	fprintf(helpOut, "#\t-CS\t\tStart Chain size\t\t\t1 M\n");
	fprintf(helpOut, "#\t-ME\t\tMega DB\t\t\t\t\tFalse\n");
	fprintf(helpOut, "#\t-NI\t\tDo not dump *.index.b\t\t\tFalse\n");
	fprintf(helpOut, "#\t-Sparse\t\tMake Sparse DB ('-' for no prefix)\tNone/False\n");
	fprintf(helpOut, "#\t-ht\t\tHomology template\t\t\t1.0\n");
	fprintf(helpOut, "#\t-hq\t\tHomology query\t\t\t\t1.0\n");
	fprintf(helpOut, "#\t-and\t\tBoth homolgy thresholds\n#\t\t\thas to be reached\t\t\tor\n");
	fprintf(helpOut, "#\t-v\t\tVersion\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int i, args, stop, filecount, deconcount, sparse_run;
	int size, mapped_cont, file_len, appender;
	unsigned megaDB, **Values;
	char **inputfiles, *outputfilename, *templatefilename, **deconfiles;
	char *to2Bit, *line, *exeBasic;
	unsigned char *update;
	struct hashMapKMA *finalDB;
	FILE *inputfile, *out;
	time_t t0, t1;
	
	if (argc == 1) {
		fprintf(stderr, "# Too few arguments handed.\n");
		helpMessage(-1);
	} else if(sizeof(long unsigned) != 8) {
		ERROR();
	}
	
	/* set defaults */
	INITIAL_SIZE = 1048576;
	templates = 0;
	kmersize = 16;
	kmerindex = 16;
	sparse_run = 0;
	appender = 0;
	MinLen = 0;
	MinKlen = 1;
	prefix_len = 0;
	prefix = 0;
	homQ = 1;
	homT = 1;
	homcmp = &homcmp_or;
	dumpIndex = &makeIndexing;
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
	to2Bit = malloc(384);
	if(!inputfiles || !deconfiles || !to2Bit) {
		ERROR();
	}
	/* set to2Bit */
	for(i = 0; i < 384; ++i) {
		to2Bit[i] = 8;
	}
	to2Bit += 128;
	to2Bit['\n'] = 16;
	
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
			++args;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strcmp(argv[args], "--") == 0) {
					++filecount;
					inputfiles = realloc(inputfiles, filecount * sizeof(char*));
					if(inputfiles == NULL) {
						ERROR();
					}
					inputfiles[filecount - 1] = strdup(argv[args]);
					if(!inputfiles[filecount - 1]) {
						ERROR();
					}
					++args;
				} else {
					stop = 1;
				}
			}
			--args;
		} else if(strcmp(argv[args], "-o") == 0) {
			++args;
			if(args < argc) {
				outputfilename = malloc(strlen(argv[args]) + 64);
				if(!outputfilename) {
					ERROR();
				}
				strcpy(outputfilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-deCon") == 0) {
			stop = 0;
			++args;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strcmp(argv[args], "--") == 0) {
					deconfiles = realloc(deconfiles, (deconcount + 1) * sizeof(char*));
					if(deconfiles == NULL) {
						ERROR();
					}
					deconfiles[deconcount] = strdup(argv[args]);
					if(deconfiles[deconcount] == NULL) {
						ERROR();
					}
					++deconcount;
					++args;
				} else {
					stop = 1;
				}
			}
			if(deconcount == 0) {
				fprintf(stderr, "No deCon file specified.\n");
				exit(1);
			}
			--args;
		} else if(strcmp(argv[args], "-t_db") == 0) {
			++args;
			if(args < argc) {
				templatefilename = malloc(strlen(argv[args]) + 64);
				if(!templatefilename) {
					ERROR();
				}
				strcpy(templatefilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-k") == 0) {
			++args;
			if(args < argc) {
				kmersize = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid kmersize parsed\n");
					exit(4);
				} else if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 32) {
					kmersize = 32;
				}
				kmerindex = kmersize;
			}
		} else if(strcmp(argv[args], "-k_t") == 0) {
			++args;
			if(args < argc) {
				kmersize = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid kmersize parsed\n");
					exit(4);
				} else if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 32) {
					kmersize = 32;
				}
			}
		} else if(strcmp(argv[args], "-k_i") == 0) {
			++args;
			if(args < argc) {
				kmerindex = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid kmersize parsed\n");
					exit(4);
				} else if(kmerindex == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmerindex = 16;
				} else if(kmerindex > 32) {
					kmerindex = 32;
				}
			}
		} else if(strcmp(argv[args], "-CS") == 0) {
			++args;
			if(args < argc) {
				
				size = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid start size parsed\n");
					exit(4);
				}
				INITIAL_SIZE = pow(2, ceil(log(size)/log(2))) + 0.5;
				INITIAL_SIZE *= 1048576;
				if(INITIAL_SIZE == 0) {
					fprintf(stderr, "# Invalid Chain Size parsed, using default\n");
					INITIAL_SIZE = 1048576;
				}
			}
		} else if(strcmp(argv[args], "-and") == 0) {
			homcmp = &homcmp_and;
		} else if(strcmp(argv[args], "-ML") == 0) {
			++args;
			if(args < argc) {
				MinLen = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum length parsed\n");
					exit(4);
				} else if(MinLen <= 0) {
					fprintf(stderr, "# Invalid minimum length parsed, using default\n");
					MinLen = 0;
				}
			}
		} else if(strcmp(argv[args], "-hq") == 0) {
			++args;
			if(args < argc) {
				homQ = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-hq\".\n");
					exit(4);
				} else if(homQ < 0) {
					fprintf(stderr, "Invalid -hq\n");
					homQ = 1.0;
				}
			}
		} else if(strcmp(argv[args], "-ht") == 0) {
			++args;
			if(args < argc) {
				homT = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-ht\".\n");
					exit(4);
				} else if(homT < 0) {
					fprintf(stderr, "Invalid -hq\n");
					homT = 1.0;
				}
			}
		} else if(strcmp(argv[args], "-batch") == 0) {
			++args;
			if(args < argc) {
				inputfile = fopen(argv[args], "rb");
				if(!inputfile) {
					ERROR();
				}
				fseek(inputfile, 0, SEEK_END);
				size = ftell(inputfile) + 1;
				rewind(inputfile);
				
				++filecount;
				inputfiles = realloc(inputfiles, filecount * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				inputfiles[filecount - 1] = malloc(size);
				if(!inputfiles[filecount - 1]) {
					ERROR();
				}
				
				/* get number of file */
				fread(inputfiles[filecount - 1], 1, size - 1, inputfile);
				fclose(inputfile);
				
				i = size;
				size = 0;
				line = inputfiles[filecount - 1];
				while(--i) {
					if(line[i] == '\n') {
						size++;
						while(isspace(line[i])) {
							line[i] = 0;
							--i;
						}
					}
				}
				--size;
				inputfiles = realloc(inputfiles, (filecount + size) * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				for(i = size; i; --i) {
					while(*line != 0) {
						++line;
					}
					while(*line == 0) {
						++line;
					}
					inputfiles[filecount] = line;
					++filecount;
				}
			}
		} else if(strcmp(argv[args], "-batchD") == 0) {
			++args;
			if(args < argc) {
				inputfile = fopen(argv[args], "rb");
				if(!inputfile) {
					ERROR();
				}
				fseek(inputfile, 0, SEEK_END);
				size = ftell(inputfile) + 1;
				rewind(inputfile);
				
				++deconcount;
				deconfiles = realloc(deconfiles, deconcount * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				deconfiles[deconcount - 1] = malloc(size);
				if(!deconfiles[deconcount - 1]) {
					ERROR();
				}
				
				/* get number of file */
				fread(deconfiles[deconcount - 1], 1, size - 1, inputfile);
				fclose(inputfile);
				
				i = size;
				size = 0;
				line = deconfiles[deconcount - 1];
				while(--i) {
					if(line[i] == '\n') {
						size++;
						while(isspace(line[i])) {
							line[i] = 0;
							--i;
						}
					}
				}
				--size;
				deconfiles = realloc(deconfiles, (deconcount + size) * sizeof(char*));
				if(!inputfiles) {
					ERROR();
				}
				for(i = size; i; --i) {
					while(*line != 0) {
						++line;
					}
					while(*line == 0) {
						++line;
					}
					deconfiles[deconcount] = line;
					++deconcount;
				}
			}
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
			++args;
			if(args < argc) {
				if(strcmp(argv[args], "-") == 0) {
					prefix_len = 0;
					prefix = 1;
				} else {
					prefix_len = strlen(argv[args]);
					prefix = 0;
					update = (unsigned char *) argv[args];
					for(i = 0; i < prefix_len; ++i) {
						prefix = (prefix << 2) | to2Bit[update[i]];
						if(to2Bit[update[i]] > 3) {
							fprintf(stderr, "Invalid prefix.\n");
							exit(1);
						}
					}
					if(prefix_len == 0) {
						fprintf(stderr, "Invalid prefix.\n");
						exit(1);
					}
				}
			}
		} else if(strcmp(argv[args], "-ME") == 0) {
			megaDB = 1;
		} else if(strcmp(argv[args], "-NI") == 0) {
			dumpIndex = &dumpSeq;
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA_index-%d.%d.%d\n", version[0], version[1], version[2]);
			exit(0);
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(1);
		}
		++args;
	}
	
	/* check for sufficient input */
	if(filecount == 0 && deconcount == 0) {
		fprintf(stderr, "No inputfiles defined.\n");
		helpMessage(-1);
	} else if(filecount == 0 && deconcount != 0 && templatefilename == 0) {
		fprintf(stderr, "Nothing to update.\n");
		exit(0);
	} else if(outputfilename == 0 && templatefilename == 0) {
		fprintf(stderr, "Output destination not defined.\n");
		helpMessage(-1);
	} else if(outputfilename == 0 && templatefilename != 0) {
		outputfilename = smalloc((strlen(templatefilename) + 64));
		strcpy(outputfilename, templatefilename);
	}
	file_len = strlen(outputfilename);
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	if(megaDB) {
		INITIAL_SIZE = mask;
		++INITIAL_SIZE;
	} else if(INITIAL_SIZE == (mask + 1)) {
		megaDB = 1;
	}
	
	/* load DB */
	if(templatefilename != 0) {
		/* load */
		fprintf(stderr, "# Loading database: %s\n", outputfilename);
		finalDB = load_DBs(templatefilename, outputfilename);
		
		/* determine params based on loaded DB */
		if(prefix_len == 0 && prefix == 0) {
			sparse_run = 0;
		} else {
			sparse_run = 1;
		}
		if(mask == finalDB->size) {
			megaDB = 1;
			INITIAL_SIZE = mask + 1;
		} else if(megaDB == 0) {
			INITIAL_SIZE = finalDB->size + 1;
		}
		appender = 1;
	} else {
		finalDB = 0;
		appender = 0;
	}
	
	/* set pointers and dependent global variables */
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	shifterI = sizeof(long unsigned) * sizeof(long unsigned) - (kmerindex << 1);
	prefix_shifter = sizeof(long unsigned) * sizeof(long unsigned) - (prefix_len << 1);
	
	/* function pointers */
	hashMap_add = &hashMap_addKMA;
	hashMap_get = &hashMap_getValue;
	addCont = &hashMap_addCont;
	
	if(DB_size < HU_LIMIT) {
		updateValuePtr = &updateShortValue;
		valuesKeyPtr = &huValuesKey;
		cmpValuesPtr = &cmpHuValues;
		valuesSize = &huSize;
		updateScoreAndTemplate_ptr = &updateScoreAndTemplateHU;
		addUscore_ptr = &addUscoreHU;
	} else {
		updateValuePtr = &updateValue;
		valuesKeyPtr = &valuesKey;
		cmpValuesPtr = &cmpValues;
		valuesSize = &uSize;
		updateScoreAndTemplate_ptr = &updateScoreAndTemplate;
		addUscore_ptr = &addUscore;
	}
	if(sparse_run) {
		update_DB = &updateDBs_sparse;
		updateAnnotsPtr = &updateAnnots_sparse;
		if(prefix_len == 0) {
			prefix = 1;
		}
	} else {
		update_DB = &updateDBs;
		updateAnnotsPtr = &updateAnnots;
	}
	
	if(prefix_len != 0) {
		deConNode_ptr = &deConNode_sparse;
	} else {
		deConNode_ptr = &deConNode;
	}
	if(megaDB) {
		hashMap_add = &megaMap_addKMA;
		hashMap_get = &megaMap_getValue;
		addCont = &megaMap_addCont;
	}
	
	/* set homology check */
	if(MinLen > (kmersize + prefix_len + 1)) {
		MinKlen = 2 * (MinLen - kmersize - prefix_len + 1);
		for(i = 0; i < prefix_len; ++i) {
			MinKlen /= 4;
		}
	} else {
		MinLen = MAX(kmersize, kmerindex);
	}
	if(homT < 1) {
		QualCheck = &templateCheck;
		foundKmers = malloc(sizeof(struct hashMap_kmers));
		if(!foundKmers) {
			ERROR();
		}
		foundKmers->size = INITIAL_SIZE;
		foundKmers->table = calloc(foundKmers->size, sizeof(struct hashTable_kmers *));
		Scores = calloc(1024, sizeof(unsigned));
		Scores_tot = calloc(1024, sizeof(unsigned));
		bestTemplates = malloc(1024 * sizeof(unsigned));
		if(!foundKmers->table || !Scores || !Scores_tot || !bestTemplates) {
			ERROR();
		}
		Scores_tot[0] = 1024;
	} else if(homQ < 1) {
		QualCheck = &queryCheck;
		Scores_tot = calloc(1024, sizeof(unsigned));
		bestTemplates = malloc(1024 * sizeof(unsigned));
		if(!Scores_tot || !bestTemplates) {
			ERROR();
		}
		Scores_tot[0] = 1024;
	} else {
		QualCheck = &lengthCheck;
	}
	
	/* update DBs */
	if(filecount != 0) {
		if(finalDB) {
			/* convert */
			templates = hashMapKMA_openChains(finalDB);
		} else {
			/* create */
			templates = hashMap_initialize(INITIAL_SIZE);
			
			DB_size = 1;
			if(sparse_run) {
				templates->prefix = prefix;
				templates->prefix_len = prefix_len;
				
				template_lengths = malloc(1024 * sizeof(unsigned));;
				template_slengths = malloc(1024 * sizeof(unsigned));
				template_ulengths = malloc(1024 * sizeof(unsigned));
				if(!template_lengths || !template_slengths || !template_ulengths) {
					ERROR();
				}
				*template_lengths = kmerindex;
				template_slengths[0] = 1024;
				template_ulengths[0] = 1024;
			} else {
				template_lengths = malloc(1024 * sizeof(unsigned));
				template_slengths = 0;
				template_ulengths = 0;
				if(!template_lengths) {
					ERROR();
				}
				template_lengths[0] = 1024;
			}
		}
		
		fprintf(stderr, "# Indexing databases.\n");
		t0 = clock();
		makeDB(inputfiles, filecount, outputfilename, appender, to2Bit);
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for DB indexing: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
		
		
		/* compress db */
		fprintf(stderr, "# Compressing templates\n");
		t0 = clock();
		strcat(outputfilename, ".comp.b");
		out = fopen(outputfilename, "wb");
		if(!out) {
			ERROR();
		} else if(templates->table != 0) {
			finalDB = compressKMA_DB(out);
		} else {
			finalDB = compressKMA_megaDB(out);
		}
		fclose(out);
		outputfilename[file_len] = 0;
		free(templates);
		fprintf(stderr, "# Template database created.\n");
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for DB compression: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else {
		++finalDB->size;
	}
	
	/* decontaminate */
	if(deconcount != 0) {
		/* open values */
		fprintf(stderr, "# Openning values\n");
		Values = hashMapKMA_openValues(finalDB);
		
		/* get decontamination info */
		fprintf(stderr, "# Adding decontamination information\n");
		t0 = clock();
		mapped_cont = deConDB(finalDB, deconfiles, deconcount, to2Bit, Values);
		fprintf(stderr, "# Contamination information added.\n");
		fprintf(stderr, "# %d kmers mapped to the DB.\n", mapped_cont);
		fprintf(stderr, "# Contamination mapped to %f %% of the DB.\n", 100.0 * mapped_cont / finalDB->n);
		
		/* compress DB */
		fprintf(stderr, "# Compressing templates\n");
		if((finalDB->size - 1) != mask) {
			compressKMA_deconDB(finalDB, Values);
		} else {
			compressKMA_deconMegaDB(finalDB, Values);
		}
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for DB decontamination: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
		
		/* dump DB */
		fprintf(stderr, "# Dumping DB.\n");
		strcat(outputfilename, ".decon.comp.b");
		out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
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
