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

#include <stdint.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <paths.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

#define HU_LIMIT 65535
#define I_LIMIT 2147483647
#define U_LIMIT 4294967295
#define CHUNK 1048576
#define windowBits 31
#define GZIP_ENCODING 16
#define ENABLE_ZLIB_GZIP 32
#define getNuc(Comp,pos) ((Comp[pos >> 5] << ((pos & 31) << 1)) >> 62)
#define lock(exclude) while(__sync_lock_test_and_set(exclude, 1)) {while(*exclude) {usleep(100);}}
#define lockTime(exclude, spin) while(__sync_lock_test_and_set(exclude, 1)) {while(*exclude) {usleep(spin);}}
#define unlock(exclude) (__sync_lock_release(exclude))
#define wait_atomic(src) while(src) {usleep(100);}
#define MIN(X, Y) ((X < Y) ? X : Y)
#define MAX(X, Y) ((X < Y) ? Y : X)
#define setEx(src, pos)(src[pos >> 3] |= (1 << (pos & 7)))
#define unsetEx(src, pos)(src[pos >> 3] ^= (1 << (pos & 7)))
#define getEx(src, pos)((src[pos >> 3] >> (pos & 7)) & 1)
#define ERROR() fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno)); exit(errno);

/*
 STRUCTURES
*/
struct aln {
	unsigned char *t;  /* template */
	char *s;  /* score */
	unsigned char *q;  /* query */
	unsigned pos; /* start of aln, relative to template */
	unsigned mapQ; /* mapping quality */
	int score; /* aln score */
};

struct assem {
	unsigned char *t;  /* template */
	char *s;  /* score */
	unsigned char *q;  /* query */
	unsigned cover;
	long unsigned depth;
	long unsigned depthVar;
	long unsigned score;
	unsigned len;
	unsigned aln_len;
	unsigned size;
};

struct frag {
	int buffer[6];
	unsigned char *qseq;
	unsigned char *header;
	struct frag *next;
};

struct hashTable {
	long unsigned key;
	unsigned *value;
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

static struct pid {
	FILE *fp;
	pid_t pid;
	struct pid *next;
} *pidlist;

struct compDNA {
	int seqlen;
	int size;
	int complen;
	long unsigned *seq;
	int *N;
};

struct compKmers {
	long unsigned n;
	long unsigned size;
	long unsigned *kmers;
};

struct hashMapKMA {
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

struct assembly {
	short unsigned counts[6];
	unsigned next;
};

struct assemInfo {
	int len;
	int size;
	struct assembly *assmb;
};

struct alnPoints {
	int size;
	int len;
	int *tStart;
	int *tEnd;
	int *qStart;
	int *qEnd;
	int *weight;
	int *score;
	int *next;
};

struct NWmat {
	long NW_s;
	long NW_q;
	unsigned char *E;
	int *D[2];
	int *P[2];
};

struct kmerScan_thread {
	pthread_t id;
	int num;
	int bestScore;
	int bestScore_r;
	int *bestTemplates;
	int *bestTemplates_r;
	FILE *inputfile;
	struct compDNA *qseq;
	struct compDNA *qseq_r;
	struct qseqs *header;
	struct kmerScan_thread *next;
};

struct assemble_thread {
	pthread_t id;
	int num;
	int template;
	int file_count;
	int spin;
	FILE **files;
	struct FileBuff *frag_out;
	struct assem *aligned_assem;
	struct aln *aligned, *gap_align;
	struct qseqs *qseq, *header;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	struct assemble_thread *next;
};

struct aln_thread {
	pthread_t id;
	int *matched_templates;
	int *bestTemplates;
	int *bestTemplates_r;
	int *best_start_pos;
	int *best_end_pos;
	long *index_indexes;
	long *seq_indexes;
	FILE *inputfile;
	FILE *frag_out_raw;
	FILE *index_in;
	FILE *seq_in;
	struct compDNA *qseq_comp;
	struct compDNA *qseq_r_comp;
	struct qseqs *qseq;
	struct qseqs *qseq_r;
	struct qseqs *header;
	struct qseqs *header_r;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	struct aln_thread *next;
};

struct spltDBbuff {
	int num;
	int size;
	unsigned char *buff;
	struct spltDBbuff *next;
	struct spltDBbuff *prev;
};

/*
 	GLOBAL VARIABLES
*/
int version[3] = {1, 1, 7};
struct hashMapKMA *templates;
struct hashMap_index **templates_index;
struct diskOffsets *templates_offsets;
char **template_names;
int *template_lengths, *template_ulengths;
int **RegionTemplates, *Score, *Score_r, **TmpNs, exhaustive, minLen, ref_fsa;
int **tScore, **tScore_r, **BestTemplates, **BestTemplates_r, one2one, diskDB;
int delta, deCon, contamination, print_matrix, print_all, DB_size, kmersize;
unsigned ***tVF_scores, ***tVR_scores, *valuesFile;
unsigned shifter, shifterS, r_shifter, shm, thread_num, countK, mq;
long unsigned *alignment_scores, *uniq_alignment_scores, *Mt1Score;
long unsigned mask;
double evalue, ID_t, scoreT, support, HMM_param[8];
int d[5][5], W1, U, M, MM, PE;
extern char **environ;
volatile int *excludeIn, *excludeOut;

/*
	Set function pointers 
*/
void (*printPtr)(int*, struct compDNA*, int, struct qseqs*);
void (*printPairPtr)(int*, struct compDNA*, int, struct qseqs*, struct compDNA*, int, struct qseqs*);
void (*deConPrintPtr)(int*, struct compDNA*, int, struct qseqs*);
void (*ankerPtr)(int*, int*, int*, unsigned**, unsigned**, int*, struct compDNA*, int, int, int, int, struct qseqs*);
void * (*assembly_KMA_Ptr)(void *);
void (*destroyPtr)(int);
struct hashMap_index * (*alignLoadPtr)(FILE*, FILE*, int, long unsigned, long unsigned);
unsigned * (*hashMap_get)(const struct hashMapKMA *, const long unsigned);
void (*kmerScan)(int*, int*, int*, int*, struct compDNA*, struct compDNA*, struct qseqs*, int*);
void (*save_kmers_pair)(int*, int*, int*, int*, int*, int*, struct compDNA*, struct compDNA*, struct qseqs*, struct qseqs*, int*);
void (*printFsa_ptr)(struct qseqs*, struct qseqs*, struct compDNA*);
void (*printFsa_pair_ptr)(struct qseqs*, struct qseqs*, struct qseqs*, struct qseqs*, struct compDNA*);
void (*alnFragsPE)(int*, struct compDNA*, struct compDNA*, unsigned char*, unsigned char*, struct qseqs*, struct qseqs*, int*, int*, int*, int*, FILE*, FILE*, long*, long*, FILE*, struct alnPoints *, struct NWmat *, volatile int *, volatile int *);
int (*cmp)(int, int);
int (*significantBase)(int, int);
unsigned char (*baseCall)(unsigned char, unsigned char, int, int, struct assembly*);
int (*buffFileBuff)(struct FileBuff *);
long unsigned (*getExistPtr)(const unsigned *, const long unsigned);
long unsigned (*getKeyPtr)(const unsigned *, const long unsigned);
long unsigned (*getValueIndexPtr)(const unsigned *, const long unsigned);
unsigned * (*getValuePtr)(const struct hashMapKMA *, const long unsigned);
int (*get_kmers_for_pair_ptr)(int *, int *, int *, int *, struct compDNA *, int *);
int (*intpos_bin_contaminationPtr)(const unsigned *, const int);
int (*chainSeedsPtr)(struct alnPoints *, int, int, unsigned *);

/*
 FUNCTIONS
*/

/* BASIC FUNCTIONS */
void * smalloc(const size_t size) {
	
	void *dest = malloc(size);
	if(!dest) {
		ERROR();
	}
	
	return dest;
}

FILE * sfopen(char *filename, char *mode) {
	
	FILE *file = fopen(filename, mode);
	if(!file) {
		fprintf(stderr, "Filename:\t%s\n", filename);
		ERROR();
	}
	
	return file;
}

#define sfwrite(ptr, size, nmemb, stream) if(fwrite(ptr, size, nmemb, stream) != nmemb) {ERROR();}
#define sfread(ptr, size, nmemb, stream) if(fread(ptr, size, nmemb, stream) != nmemb) {ERROR();}

FILE * kmaPopen(const char *cmd, const char *type) {
	
	/* kmaPopen does the same as popen, but allows for binary mode */
	
	int pdes[2];
	FILE *ioStream;
	pid_t pid;
	struct pid * volatile dest;
	char *argv[] = {"sh", "-c", NULL, NULL};
	
	/* check mode */
	if(*type != 'r' && *type != 'w') {
		errno = EINVAL;
		ERROR();
	}
	
	/* create pipe */
	dest = malloc(sizeof(struct pid));
	if(!dest) {
		ERROR();
	} else if(pipe(pdes) != 0) {
		ERROR();
	}
	
	/* spawn process */
	pid = vfork();
	if(pid < 0) {
		ERROR();
	} else if(pid == 0) {
		/* close filedescriptors */
		struct pid *pidPtr;
		for(pidPtr = pidlist; pidPtr; pidPtr = pidPtr->next) {
			close(fileno(pidPtr->fp));
		}
		
		/* error handling */
		if(*type == 'r') {
			close(pdes[0]);
			if(pdes[1] != STDOUT_FILENO) {
				dup2(pdes[1], STDOUT_FILENO);
				close(pdes[1]);
			}
		} else {
			close(pdes[1]);
			if(pdes[0] != STDIN_FILENO) {
				dup2(pdes[0], STDIN_FILENO);
				close(pdes[0]);
			}
		}
		
		/* start child work */
		argv[2] = (char *) cmd;
		execve(_PATH_BSHELL, argv, environ);
		
		/* kill child */
		_exit(127);
	}
	
	/* Parent work */
	if (*type == 'r') {
		ioStream = fdopen(pdes[0], type);
		close(pdes[1]);
	} else {
		ioStream = fdopen(pdes[1], type);
		close(pdes[0]);
	}
	if(!ioStream) {
		ERROR();
	}
	
	/* Link into list of file descriptors. */
	dest->fp = ioStream;
	dest->pid = pid;
	dest->next = pidlist;
	pidlist = dest;
	
	return ioStream;
}

int kmaPclose(FILE *ioStream) {
	
	struct pid *dest, *last;
	int status;
	pid_t pid;
	
	/* Get file pointer. */
	for (last = 0, dest = pidlist; dest->fp != ioStream; last = dest, dest = dest->next) {
		if(!dest) {
			return -1;
		}
	}
	
	/* close stream and get exit status */
	fclose(ioStream);
	status = -1;
	while ((pid = waitpid(dest->pid, &status, 0)) == -1 && errno == EINTR) {
		usleep(100);
	}
	
	/* Remove the entry from the linked list. */
	if (!last) {
		pidlist = dest->next;
	} else {
		last->next = dest->next;
	}
	free(dest);
	
	return status;
}

int cmp_or(int t, int q) {
	return (t || q);
}

int cmp_and(int t, int q) {
	return (t && q);
}

double fastp(long double q) {
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

double p_chisqr(long double q) {
	
	if(q < 0) {
		/* Handle negative overflow */
		return 1e-26;
	} else if(q > 49) {
		/* Get p-val from table, to avoid overflow */
		return fastp(q);
	}
	/* claculate p-value */
	return 1 - 1.772453850 * erf(sqrt(0.5 * q)) / tgamma(0.5);
}

double binP(int n, int k, double p) {
	
	int i, j, nk;
	double P, q, pq;
	
	/*
		P = n! / (k! (n-k)!) * p^k * q^(n - k)
	*/
	
	q = 1 - p;
	if(k == 0) {
		return pow(q, n);
	} else if(n == k) {
		return pow(p, n);
	} else if(p == 0 || q == 0) {
		return 0.0;
	}
	
	P = 1.0;
	nk = n - k;
	pq = p * q;
	i = n + 1;
	if(k < nk) {
		j = k + 1;
	} else {
		j = nk + 1;
	}
	
	while(--j) {
		P *= (--i * pq / j);
	}
	
	if(nk < k) {
		P *= pow(p, k - nk);
	} else if(k < nk) {
		P *= pow(q, nk - k);
	}
	
	return P != 0.0 ? P : 1.0e-308;
}

unsigned minimum(unsigned *src, unsigned n) {
	
	unsigned min;
	
	min = *src;
	while(--n) {
		if(*++src < min) {
			min = *src;
		}
	}
	
	return min;
}

long unsigned getKmer(long unsigned *compressor, unsigned cPos) {
	
	unsigned iPos = (cPos & 31) << 1;
	cPos >>= 5;
	
	return (iPos <= shifter) ? ((compressor[cPos] << iPos) >> shifter) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifter);
}

long unsigned getSmer(long unsigned *compressor, unsigned cPos) {
	
	unsigned iPos = (cPos & 31) << 1;
	cPos >>= 5;
	
	return (iPos <= shifterS) ? ((compressor[cPos] << iPos) >> shifterS) : (((compressor[cPos] << iPos) | (compressor[cPos + 1] >> (64-iPos))) >> shifterS);
}

long unsigned getK(long unsigned *compressor, unsigned pos) {
	return pos;
}

long unsigned makeKmer(const unsigned char *qseq, unsigned pos, unsigned size) {
	
	long unsigned key = qseq[pos];
	
	size += pos;
	for(++pos; pos < size; ++pos) {
		key = (key << 2) | qseq[pos];
	}
	
	return key;
}

int charpos(const unsigned char *src, unsigned char target, int start, int len) {
	
	int i;
	
	for(i = start; i < len; ++i) {
		if(src[i] == target) {
			return i;
		}
	}
	
	return -1;
}

int rcharpos(const unsigned char *src, unsigned char target, int start, int end) {
	
	int i;
	
	for(i = end; i >= 0; --i) {
		if(src[i] == target) {
			return i;
		}
	}
	
	return -1;
}

int strpos(const char* str1, const char* str2) {
	
	int i, len1, len2;
	char *strp;
	
	len1 = strlen(str1);
	len2 = strlen(str2);
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	strp = (char*)(str1);
	for(i = 0; i <= len1 - len2; ++i) {
		if(*strp == *str2) {
			if(strncmp(strp,str2,len2)==0)
				return i;
		}
		++strp;
	}
	return -1;
}

int strrpos(const char* str1, const char* str2) {
	
	int i, len1, len2;
	char *strp;
	
	len1 = strlen(str1);
	len2 = strlen(str2);
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	strp = (char*)(str1 + len1 - len2);
	for(i = len1 - len2; i >= 0; --i) {
		if(*strp == *str2) {
			if(strncmp(strp,str2,len2)==0)
				return i;
		}
		--strp;
	}
	return -1;
}

char * noFolder(const char *src) {
	
	int pos;
	
	pos = strlen(src) - 1;
	while(pos && src[pos] != '/') {
		--pos;
	}
	if(pos) {
		++pos;
	}
	
	return ((char *) src) + pos;
}

int intpos_bin(const unsigned *str1, const int str2) {
	
	int pos, upLim, downLim;
	
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

int intpos_bin_contamination(const unsigned *str1, const int str2) {
	
	int pos, upLim, downLim, template;
	
	upLim = *str1;
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) / 2;
	while(0 < (upLim - downLim)) {
		template = str1[pos];
		if(template == str2) {
			return pos;
		} else if(template < str2) {
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

int intpos_bin_contamination_s(const unsigned *Str1, const int str2) {
	
	int pos, upLim, downLim, template;
	short unsigned *str1 = (short unsigned *) Str1;
	
	upLim = *str1;
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) / 2;
	while(0 < (upLim - downLim)) {
		template = str1[pos];
		if(template == str2) {
			return pos;
		} else if(template < str2) {
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
	while(0 <= k && isspace(string[k])) {
		--k;
	}
	++k;
	string[k] = 0;
	return k;
}

void insert(char *dest, char src, int location, int dest_len) {
	int i;
	dest[dest_len + 1] = 0;
	for(i = dest_len; i > location; --i) {
		dest[i] = dest[i - 1];
	}
	dest[location] = src;
}

int replace_chars(char *dest, char src) {
	
	int i, bias, len;
	
	while(*dest == src) {
		++dest;
	}
	len = strlen(dest);
	if(len == 0)
		return len;
	
	bias = 0;
	for(i = 1; i < len && dest[i]; ++i) {
		if(dest[i] == src) {
			++bias;
		} else {
			dest[i - bias] = dest[i];
		}
	}
	len -= bias;
	return len;
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
	strm->avail_out = 0;
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

void openFileBuff(struct FileBuff *dest, char *filename, char *mode) {
	
	dest->file = fopen(filename, mode);
	if(!dest->file) {
		ERROR();
	}
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

z_stream * strm_init() {
	
	z_stream *strm;
	int status;
	
	strm = malloc(sizeof(z_stream));
	if(!strm) {
		ERROR();
	}
	
	strm->zalloc = Z_NULL;
	strm->zfree  = Z_NULL;
	strm->opaque = Z_NULL;
	
	status = deflateInit2(strm, 1, Z_DEFLATED, 31 | GZIP_ENCODING, 9, Z_FILTERED);
	if(status < 0) {
		fprintf(stderr, "Gzip error %d\n", status);
		exit(status);
	}
	
	return strm;
}

struct FileBuff * gzInitFileBuff(int size) {
	
	struct FileBuff *dest;
	
	dest = malloc(sizeof(struct FileBuff));
	if(!dest) {
		ERROR();
	}
	
	dest->file = 0;
	dest->bytes = size;
	dest->buffSize = size;
	dest->strm = strm_init();
	dest->buffer = malloc(size);
	dest->inBuffer = malloc(size);
	dest->next = dest->buffer;
	if(!dest->buffer || !dest->inBuffer) {
		ERROR();
	}
	
	return dest;
}

void resetGzFileBuff(struct FileBuff *dest, int size) {
	
	dest->bytes = size;
	dest->buffSize = size;
	free(dest->buffer);
	free(dest->inBuffer);
	dest->buffer = malloc(size);
	dest->inBuffer = malloc(size);
	dest->next = dest->buffer;
	if(!dest->buffer || !dest->inBuffer) {
		ERROR();
	}
}

void writeGzFileBuff(struct FileBuff *dest) {
	
	int check = Z_OK;
	z_stream *strm = dest->strm;
	strm->avail_in = dest->buffSize - dest->bytes;
	strm->next_in = dest->buffer;
	strm->avail_out = 0;
	
	while(strm->avail_out == 0 && check != Z_STREAM_END) {
		strm->avail_out = dest->buffSize;
		strm->next_out = dest->inBuffer;
		check = deflate(strm, Z_NO_FLUSH);
		sfwrite(dest->inBuffer, 1, dest->buffSize - strm->avail_out, dest->file);
	}
	dest->bytes = dest->buffSize;
	dest->next = dest->buffer;
}

void closeGzFileBuff(struct FileBuff *dest) {
	
	int check = Z_OK;
	z_stream *strm = dest->strm;
	strm->avail_in = dest->buffSize - dest->bytes;
	strm->next_in = dest->buffer;
	strm->avail_out = 0;
	
	while(strm->avail_out == 0 && check != Z_STREAM_END) {
		strm->avail_out = dest->buffSize;
		strm->next_out = dest->inBuffer;
		check = deflate(strm, Z_FINISH);
		sfwrite(dest->inBuffer, 1, dest->buffSize - strm->avail_out, dest->file);
	}
	deflateEnd(strm);
	fclose(dest->file);
}

void destroyGzFileBuff(struct FileBuff *dest) {
	
	closeGzFileBuff(dest);
	free(dest->strm);
	free(dest->buffer);
	free(dest->inBuffer);
	free(dest);
}

void updateMatrix(struct FileBuff *dest, char *template_name, long unsigned *template_seq, struct assemInfo *matrix, int t_len) {
	
	unsigned i, pos, check, avail, asm_len;
	char *update, bases[] = "ACGTN-";
	struct assembly *assembly;
	
	/* check buffer capacity */
	check = strlen(template_name) + 2;
	if(dest->bytes < check) {
		writeGzFileBuff(dest);
	}
	update = (char *) dest->next;
	avail = dest->bytes - check;
	
	/* fill in header */
	check -= 2;
	*update++ = '#';
	memcpy(update, template_name, check);
	update += check;
	*update++ = '\n';
	
	/* fill in rows */
	asm_len = matrix->len;
	assembly = matrix->assmb;
	i = 0;
	for(pos = 0; asm_len != 0; --asm_len, pos = assembly[pos].next) {
		/* check buffer capacity */
		if(avail < 38) {
			dest->bytes = avail;
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		
		/* update with row */
		if(pos < t_len) {
			check = sprintf(update, "%c\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", bases[getNuc(template_seq, i)], assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
			++i;
		} else {
			check = sprintf(update, "-\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
		}
		avail -= check;
		update += check;
	}
	
	/* update with last newline */
	if(avail == 0) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	*update++ = '\n';
	dest->next = (unsigned char *) update;
	dest->bytes = avail - 1;
}

void initialiseVcf(struct FileBuff *fileP, char *templateFilename) {
	
	unsigned check, avail;
	char *update;
	
	update = (char *) fileP->next;
	avail = fileP->bytes;
	
	/* header stuff */
	check = sprintf(update, "##fileformat=VCFv4.2\n");
	update += check; avail -= check;
	check = sprintf(update, "##kmaVersion=%d.%d.%d\n", version[0], version[1], version[2]);
	update += check; avail -= check;
	
	/* filter stuff */
	check = sprintf(update, "##FILTER=<ID=LowQual,Description=\"Low quality\">\n");
	update += check; avail -= check;
	
	/* INFO stuff */
	check = sprintf(update, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Revised Allele Fraction\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=DEL,Number=1,Type=Float,Description=\"Fraction of Reads Containing Spanning Deletions\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##INFO=<ID=AD6,Number=6,Type=Integer,Description=\"Count of all alternative alleles: A,C,G,T,N,-\">\n");
	update += check; avail -= check;
	
	/* FORMAT stuff */
	check = sprintf(update, "##FORMAT=<ID=Q,Number=1,Type=Float,Description=\"McNemar quantile\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##FORMAT=<ID=P,Number=1,Type=Float,Description=\"McNemar p-value\">\n");
	update += check; avail -= check;
	check = sprintf(update, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter\">\n");
	update += check; avail -= check;
	if(templateFilename) {
		check = sprintf(update, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", noFolder(templateFilename));
	} else {
		check = sprintf(update, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tspltDB\n");
	}
	update += check; avail -= check;
	
	fileP->next = (unsigned char *) update;
	fileP->bytes = avail;
}

void updateVcf(char *template_name, long unsigned *template_seq, int t_len, struct assemInfo *matrix, int filter, struct FileBuff *fileP) {
	
	static const char *PASS = "PASS", *FAIL = "FAIL", *LowQual = "LowQual", *UNKNOWN = ".";
	int i, j, pos, bestScore, depthUpdate, bestBaseScore, nucNum;
	int template_name_length, check, avail, DP, AD, DEL, QUAL;
	double AF, RAF, Q, P;
	const double lnConst = -10 / log(10);
	char *FILTER, **FILTER_ptr, *update;
	unsigned char nuc, bestNuc;
	const char bases[] = "ACGTN-";
	struct assembly *assembly;
	
	if(filter == 2) {
		FILTER_ptr = &FILTER;
	} else {
		FILTER_ptr = (char **) &UNKNOWN;
	}
	template_name_length = strlen(template_name);
	update = (char *) fileP->next;
	avail = fileP->bytes;
	assembly = matrix->assmb;
	pos = 0;
	i = 0;
	do {
		/* does not handle insertions yet */
		if(pos < t_len) {
			nuc = bases[getNuc(template_seq, pos)];
			++i;
		} else {
			nuc = '-';
		}
		/* call base */
		nucNum = 0;
		bestNuc = 5;
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 5; ++j) {
			if(bestScore < assembly[pos].counts[j]) {
				bestScore = assembly[pos].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[pos].counts[j];
		}
		if(bestScore < assembly[pos].counts[j]) {
			bestScore = assembly[pos].counts[j];
			bestNuc = j;
		}
		depthUpdate += assembly[pos].counts[j];
		nucNum = bestNuc;
		bestNuc = bases[bestNuc];
		
		/* Check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < assembly[pos].counts[j]) {
						bestBaseScore = assembly[pos].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[pos].counts[5];
		}
		
		if(bestScore) {
			/* determine base at current position */
			bestNuc = baseCall(bestNuc, nuc, bestScore, depthUpdate, &assembly[pos]);
			/* discard unimportant changes */
			if(nuc != toupper(bestNuc)) {
				/* INFO */
				DP = depthUpdate;
				AD = assembly[pos].counts[nucNum];
				AF = (double) AD / DP;
				RAF = (double) bestScore / DP;
				DEL = assembly[pos].counts[5];
				/* FORMAT */
				Q = pow(depthUpdate - (bestScore << 1), 2) / depthUpdate;
				P = p_chisqr(Q);
				/* QUAL */
				//QUAL = lnConst * log(P);
				QUAL = lnConst * log(binP(DP, AD, 0.25));
				
				/* FILTER */
				if(isupper(bestNuc)) {
					FILTER = (char *) PASS;
				} else if(P <= evalue || (9 * depthUpdate) <= (10 * bestScore)) {
					FILTER = (char *) LowQual;
				} else {
					FILTER = (char *) FAIL;
				}
				
				if(avail < template_name_length + 167) {
					fileP->bytes = avail;
					writeGzFileBuff(fileP);
					avail = fileP->bytes;
					update = (char *) fileP->next;
				}
				
				strcpy(update, template_name);
				update += template_name_length; avail -= template_name_length;
				*update++ = '\t'; --avail;
				
				if(nuc != '-') {
					check = sprintf(update, "%d", i);
					update += check; avail -= check;
				} else {
					*update++ = '0';
					--avail;
				}
				*update++ = '\t';
				*update++ = '.';
				*update++ = '\t';
				avail -= 3;
				
				if(nuc != '-') {
					*update++ = nuc;
					--avail;
				} else {
					*update++ = '<';
					*update++ = nuc;
					*update++ = '>';
					avail -= 3;
				}
				*update++ = '\t';
				--avail;
				if(nuc != '-' && bestNuc == '-') {
					*update++ = '<';
					*update++ = bestNuc;
					*update++ = '>';
					avail -= 3;
				} else {
					*update++ = bestNuc;
					--avail;
				}
				check = sprintf(update, "\t%d\t%s\tDP=%d;AD=%d;AF=%.2f;RAF=%.2f;DEL=%d;", QUAL, *FILTER_ptr, DP, AD, AF, RAF, DEL);
				update += check; avail -= check;
				check = sprintf(update, "AD6=%d,%d,%d,%d,%d,%d\t", assembly[pos].counts[0], assembly[pos].counts[1], assembly[pos].counts[2], assembly[pos].counts[3], assembly[pos].counts[4], assembly[pos].counts[5]);
				update += check; avail -= check;
				check = sprintf(update, "Q:P:FT\t%.2f:%4.1e:%s\n", Q, P, FILTER);
				update += check; avail -= check;
			}
		} else if(nuc != '-') {
			FILTER = (char *) FAIL;
			if(avail < template_name_length + 105) {
				fileP->bytes = avail;
				writeGzFileBuff(fileP);
				avail = fileP->bytes;
				update = (char *) fileP->next;
			}
			check = sprintf(update, "%s\t%d\t.\t%c\t%c\t%d\t%s\t", template_name, i, nuc, '.', 0, *FILTER_ptr);
			update += check; avail -= check;
			check = sprintf(update, "DP=%d;AD=%d;AF=%.2f;RAF=%.2f;DEL=%d;AD6=%d,%d,%d,%d,%d,%d\t", 0, 0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0);
			update += check; avail -= check;
			check = sprintf(update, "Q:P:FT\t%.2f:%4.1e:%s\n", 0.0, 1.0, FILTER);
			update += check; avail -= check;
		}
	} while((pos = assembly[pos].next) != 0);
	
	fileP->next = (unsigned char *) update;
	fileP->bytes = avail;
	
}

void updateFrags(struct FileBuff *dest, struct qseqs *qseq, struct qseqs *header, char *template_name, int *stats) {
	
	int check, avail;
	char *update;
	
	avail = dest->bytes;
	check = 47 + qseq->len + header->len + strlen(template_name);
	
	/* flush buffer */
	if(avail < check) {
		writeGzFileBuff(dest);
		
		/* seq is too big, reallocate buffer */
		if(dest->bytes < check) {
			resetGzFileBuff(dest, check << 1);
		}
		avail = dest->bytes;
	}
	
	/* update buffer with fragment */
	memcpy(dest->next, qseq->seq, qseq->len);
	dest->next += qseq->len;
	avail -= qseq->len;
	/* stats */
	update = (char *) dest->next;
	check = sprintf(update, "\t%d\t%d\t%d\t%d\t%s\t", stats[0], stats[1], stats[2], stats[3], template_name);
	dest->next += check;
	avail -= check;
	/* header */
	header->seq[header->len - 1] = '\n';
	memcpy(dest->next, header->seq, header->len);
	dest->next += header->len;
	avail -= header->len;
	
	dest->bytes = avail;
	
	/* equivalent with:
	fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
	*/
}

void updateAllFrag(unsigned char *qseq, int q_len, int bestHits, int best_read_score, int *best_start_pos, int *best_end_pos, int *bestTemplates, struct qseqs *header, struct FileBuff *dest) {
	
	int i, check, avail;
	char *update, bases[] = "ACGTN-";
	
	check = q_len;
	avail = dest->bytes;
	
	if(avail < check) {
		writeGzFileBuff(dest);
		
		/* seq is too big, reallocate buffer */
		if(dest->bytes < check) {
			resetGzFileBuff(dest, check << 1);
		}
		avail = dest->bytes;
	}
	
	/* copy seq */
	update = (char *) dest->next;
	++q_len;
	while(--q_len) {
		*update++ = bases[*qseq++];
	}
	avail -= check;
	
	check = 33;
	if(avail < check) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	check = sprintf(update, "\t%d\t%d\t%d", bestHits, best_read_score, *best_start_pos);
	avail -= check;
	update += check;
	
	for(i = 1; i < bestHits; ++i) {
		if(avail < 11) {
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		check = sprintf(update, ",%d", best_start_pos[i]);
		avail -= check;
		update += check;
	}
	
	if(avail < 11) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	check = sprintf(update, "\t%d", *best_end_pos);
	avail -= check;
	update += check;
	for(i = 1; i < bestHits; ++i) {
		if(avail < 11) {
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		check = sprintf(update, ",%d", best_end_pos[i]);
		avail -= check;
		update += check;
	}
	
	if(avail < 11) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	check = sprintf(update, "\t%d", *bestTemplates);
	avail -= check;
	update += check;
	for(i = 1; i < bestHits; ++i) {
		if(avail < 11) {
			writeGzFileBuff(dest);
			avail = dest->bytes;
			update = (char *) dest->next;
		}
		check = sprintf(update, ",%d", bestTemplates[i]);
		avail -= check;
		update += check;
	}
	
	check = header->len + 1;
	if(avail < check) {
		writeGzFileBuff(dest);
		avail = dest->bytes;
		update = (char *) dest->next;
	}
	*update++ = '\t';
	header->seq[header->len - 1] = '\n';
	memcpy(update, header->seq, header->len);
	update += header->len;
	
	dest->bytes = avail - check;
	dest->next = (unsigned char *) update;
}

int openAndDetermine(struct FileBuff *inputfile, char *filename) {
	
	unsigned FASTQ;
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
	
	if(inputfile->buffer[0] == '@') { //FASTQ
		FASTQ |= 1;
	} else if(inputfile->buffer[0] == '>') { //FASTA
		FASTQ |= 2;
	} else {
		fprintf(stderr, "Cannot determine format of file:\t%s\n", filename);
		fprintf(stderr, "%-.100s\n", inputfile->buffer);
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

unsigned char * ustrdup(unsigned char *src, size_t n) {
	
	unsigned char *dest;
	
	dest = malloc(n);
	if(!dest) {
		ERROR();
	}
	memcpy(dest, src, n);
	
	return dest;
}

int FileBuffgetFsa(struct FileBuff *src, struct qseqs *header, struct qseqs *qseq, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
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
				/* chomp header */
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
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFsaSeq(struct FileBuff *src, struct qseqs *qseq, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* skip header */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
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
				/* chomp header */
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
	
	/* chomp qseq */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFq(struct FileBuff *src, struct qseqs *header, struct qseqs *qseq, struct qseqs *qual, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
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
	while((*seq++ = trans[*buff++]) != 16) {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = qseq->size;
			qseq->size <<= 1;
			qseq->seq = realloc(qseq->seq, qseq->size);
			if(!qseq->seq) {
				ERROR();
			}
			seq = qseq->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	/* skip info */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get quality */
	if(qual->size != qseq->size) {
		qual->size = qseq->size;
		free(qual->seq);
		qual->seq = malloc(qual->size);
		if(!qual->seq) {
			ERROR();
		}
	}
	qual->len = qseq->len;
	if(qual->len < avail) {
		memcpy(qual->seq, buff, qual->len);
		avail -= qual->len;
		buff += qual->len;
	} else {
		seq = qual->seq;
		memcpy(seq, buff, avail);
		seq += avail;
		size = qual->len - avail;
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
		memcpy(seq, buff, size);
	}
	qual->seq[qual->len] = 0;
	
	/* skip newline */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 1;
		}
		src->bytes = avail;
		buff = src->buffer;
	}
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFqSeq(struct FileBuff *src, struct qseqs *qseq, struct qseqs *qual, char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* skip header */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while((*seq++ = trans[*buff++]) != 16) {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = qseq->size;
			qseq->size <<= 1;
			qseq->seq = realloc(qseq->seq, qseq->size);
			if(!qseq->seq) {
				ERROR();
			}
			seq = qseq->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	/* skip info */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get quality */
	if(qual->size != qseq->size) {
		qual->size = qseq->size;
		free(qual->seq);
		qual->seq = malloc(qual->size);
		if(!qual->seq) {
			ERROR();
		}
	}
	qual->len = qseq->len;
	if(qual->len < avail) {
		memcpy(qual->seq, buff, qual->len);
		avail -= qual->len;
		buff += qual->len;
	} else {
		seq = qual->seq;
		memcpy(seq, buff, avail);
		seq += avail;
		size = qual->len - avail;
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
		memcpy(seq, buff, size);
	}
	qual->seq[qual->len] = 0;
	
	/* skip newline */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 1;
		}
		src->bytes = avail;
		buff = src->buffer;
	}
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int getPhredFileBuff(struct FileBuff *dest) {
	
	int seek;
	unsigned char *buff;
	
	buff = dest->next;
	
	while(*buff != 0) {
		seek = 3;
		while(seek && *buff != 0) {
			if(*++buff == '\n') {
				--seek;
			}
		}
		
		seek = 1;
		while(seek) {
			if(*++buff == '\n') {
				seek = 0;
			} else if(*buff < 33) {
				return 0;
			} else if(53 < *buff && *buff < 59) {
				return 33;
			} else if(*buff > 84) {
				return 64;
			}
		}
	}
	
	return 0;
}

/* DNA SPECIFIC FUNCTIONS */
int find_contamination(int *out_Tem) {
	int i;
	for(i = *out_Tem; i != 0; --i) {
		if(out_Tem[i] == contamination) {
			return i;
		}
	}
	return -1;
}

int find_contamination2(int *out_Tem, int contamination_s) {
	
	int i;
	
	for(i = *out_Tem; i != 0; --i) {
		if(out_Tem[i] == contamination_s) {
			return i;
		}
	}
	return -1;
}

void strrc(unsigned char *qseq, int q_len) {
	
	int i, j, seqlen;
	unsigned char carry, comp[6] = {3, 2, 1, 0, 4, 5};
	
	seqlen = q_len >> 1;
	
	for(i = 0, j = q_len - 1; i < seqlen; ++i, --j) {
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
void allocCompKmers(struct compKmers *compressor, int size) {
	
	compressor->n = 0;
	compressor->size = size;
	compressor->kmers = malloc(size * sizeof(long unsigned));
	if(!compressor->kmers) {
		ERROR();
	}
	
}

void reallocCompKmers(struct compKmers *compressor, int size) {
	
	compressor->kmers = realloc(compressor->kmers, size * sizeof(long unsigned));
	if(!compressor->kmers) {
		ERROR();
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
		ERROR();
	}
	
	compressor->N[0] = 0;
	
}

int pushCompKmers(struct compKmers *compressor, char *qseq, int kmersize) {
	
	int i;
	long unsigned key = 0;
	
	for(i = 0; i < kmersize; ++i) {
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

void strtranslate(unsigned char *strp, char *trans) {
	
	while(*strp) {
		*strp = trans[*strp];
		++strp;
	}
}

int translateToKmersAndDump(long unsigned *Kmers, int n, int max, unsigned char *qseq, int seqlen, long unsigned prefix, int prefix_len) {
	
	int i, end, rc;
	long unsigned key;
	
	key = 0;
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			
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
					++i;
					if(key == prefix) {
						Kmers[n] = makeKmer(qseq, i, kmersize);
						++n;
						if(n == max) {
							sfwrite(Kmers, sizeof(long unsigned), n, stdout);
							n = 0;
						}
					}
				}
				i = end + kmersize + 1;
			}
		}
	} else {
		for(rc = 0; rc < 2; ++rc) {
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
					++i;
					Kmers[n] = key;
					++n;
					if(n == max) {
						sfwrite(Kmers, sizeof(long unsigned), n, stdout);
						n = 0;
					}
				}
				i = end + kmersize + 1;
			}
		}
	}
	return n;
}

void compDNA(struct compDNA *compressor, unsigned char *seq, int seqlen) {
	
	int i, j, pos, end;
	
	compressor->seqlen = seqlen;
	if(seqlen & 31) {
		compressor->complen = (seqlen >> 5) + 1;
	} else {
		compressor->complen = seqlen >> 5;
	}
	
	
	for(i = 0, pos = 0; i < seqlen; i += 32) {
		end = (i + 32 < seqlen) ? i + 32 : seqlen;
		pos = i >> 5;
		for(j = i; j < end; ++j) {
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

void unCompDNA(struct compDNA *compressor, unsigned char *seq) {
	
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
	for(i = 0, j = compressor->complen - 1; i < compressor->complen; ++i, --j) {
		compressor_rc->seq[j] = binRev(~compressor->seq[i]);
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
	shift = compressor->seqlen - 1;
	compressor_rc->N[0] = compressor->N[0];
	for(i = 1, j = compressor->N[0]; j != 0; ++i, --j) {
		compressor_rc->N[i] = shift - compressor->N[j];
	}
}

void comp_rc(struct compDNA *compressor) {
	
	int i, j, shift, r_shift;
	long unsigned carry;
	
	/* reverse and complement*/
	for(i = 0, j = compressor->complen - 1; i < j; ++i, --j) {
		carry = binRev(~compressor->seq[i]);
		compressor->seq[i] = binRev(~compressor->seq[j]);
		compressor->seq[j] = carry;
	}
	if(i == j) {
		compressor->seq[i] = binRev(~compressor->seq[i]);
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
	
	/* add N's */
	i = 1;
	j = compressor->N[0];
	shift = compressor->seqlen - 1;
	for(i = 1, j = compressor->N[0]; i < j; ++i, --j) {
		r_shift = shift - compressor->N[i];
		compressor->N[i] = shift - compressor->N[j];
		compressor->N[j] = r_shift;
	}
	if(i == j) {
		compressor->N[i] = shift - compressor->N[j];
	}
	
}

void dumpComp(struct compDNA *compressor, FILE* file) {
	
	sfwrite(&compressor->seqlen, sizeof(int), 1, file);
	sfwrite(&compressor->complen, sizeof(int), 1, file);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, file);
	sfwrite(compressor->N, sizeof(int), compressor->N[0] + 1, file);
	
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
			ERROR();
		}
		compressor->N[0] = 0;
	}
	
	fread(compressor->seq, sizeof(long unsigned), compressor->complen, file);
	fread(compressor->N, sizeof(int), 1, file);
	fread(compressor->N + 1, sizeof(int), compressor->N[0], file);
	
	return 1;
}

void dumpCompKmers(struct compKmers *compressor, FILE *file) {
	
	sfwrite(&compressor->n, sizeof(int), 1, file);
	sfwrite(compressor->kmers, sizeof(long unsigned), compressor->n, file);
	
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
long unsigned getExist(const unsigned *exist, const long unsigned pos) {
	return exist[pos];
}

long unsigned getExistL(const unsigned *exist, const long unsigned pos) {
	return ((long unsigned *) exist)[pos];
}

long unsigned getKey(const unsigned *key_index, const long unsigned pos) {
	return key_index[pos];
}

long unsigned getKeyL(const unsigned *key_index, const long unsigned pos) {
	return ((long unsigned *) key_index)[pos];
}

long unsigned getValueIndex(const unsigned *value_index, const long unsigned pos) {
	return value_index[pos];
}

long unsigned getValueIndexL(const unsigned *value_index, const long unsigned pos) {
	return ((long unsigned *) value_index)[pos];
}

unsigned * getValue(const struct hashMapKMA *dest, const long unsigned pos) {
	return (dest->values + pos);
}

unsigned * getValueS(const struct hashMapKMA *dest, const long unsigned pos) {
	return (unsigned *)(dest->values_s + pos);
}

unsigned * hashMap_getGlobal(const struct hashMapKMA *templates, const long unsigned key) {
	
	long unsigned pos, kpos, kmer;
	
	kpos = key & templates->size;
	pos = getExistPtr(templates->exist, kpos);
	
	if(pos != templates->null_index) {
		kmer = getKeyPtr(templates->key_index, pos);
		while(key != kmer) {
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			kmer = getKeyPtr(templates->key_index, ++pos);
		}
		return getValuePtr(templates, getValueIndexPtr(templates->value_index, pos));
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
		ERROR();
	}
	
	src->exist = 3 * sizeof(unsigned) + 5 * sizeof(long unsigned);
	src->values = src->exist;
	if((templates->size - 1) == mask) {
		if(templates->v_index <= U_LIMIT) {
			src->values += templates->size * sizeof(unsigned);
		} else {
			src->values += templates->size * sizeof(long unsigned);
		}
		src->key_index = src->values;
		src->value_index = src->key_index;
	} else {
		if(templates->n <= U_LIMIT) {
			src->values += templates->size * sizeof(unsigned);
		} else {
			src->values += templates->size * sizeof(long unsigned);
		}
		
		src->key_index = src->values;
		if(DB_size < HU_LIMIT) {
			src->key_index += (templates->n + 1) * sizeof(short unsigned);
		} else {
			src->key_index += (templates->n + 1) * sizeof(unsigned);
		}
		
		src->value_index = src->key_index;
		if(kmersize <= 16) {
			src->value_index += templates->n * sizeof(unsigned);
		} else {
			src->value_index += templates->n * sizeof(long unsigned);
		}
		
	}
	
	src->file = sfopen(filename, "rb");
	
	return src;
}

unsigned getUnsignedFromFile(long pos, FILE *file) {
	
	unsigned dest;
	
	fseek(file, pos, SEEK_SET);
	fread(&dest, sizeof(unsigned), 1, file);
	
	return dest;
}

long unsigned getKmerFromFile(long offset, unsigned pos, FILE *file) {
	
	unsigned cPos, iPos;
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

unsigned * hashMap_getGlobal_disk(const struct hashMapKMA *templates, const long unsigned key) {
	
	unsigned i, pos, kpos, key_index[4];
	long unsigned kmer;
	
	kpos = key & templates->size;
	pos = getUnsignedFromFile(templates_offsets->exist + kpos * sizeof(unsigned), templates_offsets->file);
	
	if(pos != templates->null_index) {
		/* get key indexes */
		fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
		fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
		i = 0;
		if(kmersize <= 16) {
			kmer = key_index[i];
		} else {
			kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
		}
		while(key != kmer) {
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			++i;
			++pos;
			if(i == 4) {
				fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
				fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
				i = 0;
			}
			if(kmersize <= 16) {
				kmer = key_index[i];
			} else {
				kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
			}
		}
		pos = getUnsignedFromFile(templates_offsets->value_index + pos * sizeof(unsigned), templates_offsets->file);
		fseek(templates_offsets->file, templates_offsets->values + pos * sizeof(unsigned), SEEK_SET);
		fread(valuesFile, sizeof(unsigned), 1, templates_offsets->file);
		fread(valuesFile + 1, sizeof(unsigned), *valuesFile, templates_offsets->file);
		
		return valuesFile;
	}
	
	return 0;
}

unsigned * megaMap_getGlobal_disk(const struct hashMapKMA *templates, const long unsigned key) {
	
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

unsigned * hashMap_getGlobal_semDisk(const struct hashMapKMA *templates, const long unsigned key) {
	
	unsigned i, pos, kpos, key_index[4];
	long unsigned kmer;
	
	kpos = key & templates->size;
	pos = getUnsignedFromFile(templates_offsets->exist + kpos * sizeof(unsigned), templates_offsets->file);
	
	if(pos != templates->null_index) {
		/* get key indices */
		fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
		fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
		i = 0;
		if(kmersize <= 16) {
			kmer = key_index[i];
		} else {
			kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
		}
		while(key != kmer) {
			if(kpos != (kmer & templates->size)) {
				return 0;
			}
			++i;
			++pos;
			if(i == 4) {
				fseek(templates_offsets->file, templates_offsets->key_index + pos * sizeof(unsigned), SEEK_SET);
				fread(key_index, sizeof(unsigned), 4, templates_offsets->file);
				i = 0;
			}
			if(kmersize <= 16) {
				kmer = key_index[i];
			} else {
				kmer = getKmerFromFile(templates_offsets->seq, key_index[i], templates_offsets->file);
			}
		}
		pos = getUnsignedFromFile(templates_offsets->value_index + pos * sizeof(unsigned), templates_offsets->file);
		return templates->values + pos;
	}
	
	return 0;
}

unsigned * megaMap_getGlobal_semDisk(const struct hashMapKMA *templates, const long unsigned key) {
	
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
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	fread(&dest->v_index, sizeof(long unsigned), 1, file);
	fread(&dest->null_index, sizeof(long unsigned), 1, file);
	kmersize = dest->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	
	dest->exist = 0;
	dest->exist_l = 0;
	dest->values = 0;
	dest->values_s = 0;
	dest->key_index = 0;
	dest->key_index_l = 0;
	dest->value_index = 0;
	dest->value_index_l = 0;
}

unsigned * megaMap_getGlobal(const struct hashMapKMA *templates, const long unsigned key) {
	
	long unsigned pos;
	
	if((pos = getExistPtr(templates->exist, key)) != 1) {
		return getValuePtr(templates, pos);
	}
	return 0;
}

int hashMapKMA_load(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	key_t key;
	int shmid;
	long unsigned check, size, seekSize;
	
	/* load sizes */
	sfread(&contamination, sizeof(unsigned), 1, file);
	sfread(&dest->kmersize, sizeof(unsigned), 1, file);
	sfread(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfread(&dest->prefix, sizeof(long unsigned), 1, file);
	sfread(&dest->size, sizeof(long unsigned), 1, file);
	sfread(&dest->n, sizeof(long unsigned), 1, file);
	sfread(&dest->v_index, sizeof(long unsigned), 1, file);
	sfread(&dest->null_index, sizeof(long unsigned), 1, file);
	DB_size = contamination;
	kmersize = dest->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	
	/* simple check for old indexing */
	if(dest->size < dest->n) {
		return 1;
	}
	
	/* check shared memory, else load */
	size = dest->size;
	if((dest->size - 1) == mask) {
		if(dest->v_index <= U_LIMIT) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
		}
	} else {
		if(dest->n <= U_LIMIT) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
		}
	}
	key = ftok(filename, 'e');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared, load */
		dest->exist = smalloc(size);
		check = fread(dest->exist, 1, size, file);
		if(check != size) {
			return 1;
		}
		seekSize = 0;
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
		seekSize = size;
	}
	dest->exist_l = (long unsigned *)(dest->exist);
	
	/* values */
	size = dest->v_index;
	if(DB_size < HU_LIMIT) {
		size *= sizeof(short unsigned);
		getValuePtr =&getValueS;
		intpos_bin_contaminationPtr = &intpos_bin_contamination_s;
	} else {
		size *= sizeof(unsigned);
		getValuePtr =&getValue;
		intpos_bin_contaminationPtr = &intpos_bin_contamination;
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared, load */
		dest->values = smalloc(size);
		if(seekSize) {
			fseek(file, seekSize, SEEK_CUR);
		}
		seekSize = 0;
		check = fread(dest->values, 1, size, file);
		if(check != size) {
			return 1;
		}
	} else {
		/* found */
		dest->values = shmat(shmid, NULL, 0);
		seekSize += size;
	}
	dest->values_s = (short unsigned *)(dest->values);
	
	/* check for megaMap */
	if((dest->size - 1) == mask) {
		hashMap_get = &megaMap_getGlobal;
		return 0;
	} else {
		hashMap_get = &hashMap_getGlobal;
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
	key = ftok(filename, 'k');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared, load */
		dest->key_index = smalloc(size);
		if(seekSize) {
			fseek(file, seekSize, SEEK_CUR);
		}
		seekSize = 0;
		check = fread(dest->key_index, 1, size, file);
		if(check != size) {
			return 1;
		}
	} else {
		/* found */
		dest->key_index = shmat(shmid, NULL, 0);
		seekSize += size;
	}
	dest->key_index_l = (long unsigned *)(dest->key_index);
	
	/* value indexes */
	size = dest->n;
	if(dest->v_index < U_LIMIT) {
		size *= sizeof(unsigned);
		getValueIndexPtr = &getValueIndex;
	} else {
		size *= sizeof(long unsigned);
		getValueIndexPtr = &getValueIndexL;
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared, load */
		dest->value_index = smalloc(size);
		if(seekSize) {
			fseek(file, seekSize, SEEK_CUR);
		}
		seekSize = 0;
		check = fread(dest->value_index, 1, size, file);
		if(check != size) {
			return 1;
		}
	} else {
		/* found */
		dest->value_index = shmat(shmid, NULL, 0);
		seekSize += size;
	}
	dest->value_index_l = (long unsigned *)(dest->value_index);
	
	return 0;
}

void hashMapKMA_load_shm(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	key_t key;
	int shmid;
	long unsigned size;
	
	/* load sizes */
	fread(&contamination, sizeof(unsigned), 1, file);
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	fread(&dest->v_index, sizeof(long unsigned), 1, file);
	fread(&dest->null_index, sizeof(long unsigned), 1, file);
	DB_size = contamination;
	kmersize = dest->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	
	/* check shared memory */
	size = dest->size;
	if((dest->size - 1) == mask) {
		if(dest->v_index <= U_LIMIT) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
		}
	} else {
		if(dest->n <= U_LIMIT) {
			size *= sizeof(unsigned);
			getExistPtr = &getExist;
		} else {
			size *= sizeof(long unsigned);
			getExistPtr = &getExistL;
		}
	}
	key = ftok(filename, 'e');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB not shared, see kma_shm\n");
		exit(2);
	} else {
		/* found */
		dest->exist = shmat(shmid, NULL, 0);
	}
	dest->exist_l = (long unsigned *)(dest->exist);
	
	/* values */
	size = dest->v_index;
	if(DB_size < HU_LIMIT) {
		size *= sizeof(short unsigned);
		getValuePtr =&getValueS;
		intpos_bin_contaminationPtr = &intpos_bin_contamination_s;
	} else {
		size *= sizeof(unsigned);
		getValuePtr =&getValue;
		intpos_bin_contaminationPtr = &intpos_bin_contamination;
	}
	key = ftok(filename, 'v');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB not shared, see kma_shm\n");
		exit(2);
	} else {
		/* found */
		dest->values = shmat(shmid, NULL, 0);
	}
	dest->values_s = (short unsigned *)(dest->values);
	
	/* check for megaMap */
	if((dest->size - 1) == mask) {
		hashMap_get = &megaMap_getGlobal;
		return;
	} else {
		hashMap_get = &hashMap_getGlobal;
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
	key = ftok(filename, 'k');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB not shared, see kma_shm\n");
		exit(2);
	} else {
		/* found */
		dest->key_index = shmat(shmid, NULL, 0);
	}
	dest->key_index_l = (long unsigned *)(dest->key_index);
	
	/* value indexes */
	size = dest->n;
	if(dest->v_index < U_LIMIT) {
		size *= sizeof(unsigned);
		getValueIndexPtr = &getValueIndex;
	} else {
		size *= sizeof(long unsigned);
		getValueIndexPtr = &getValueIndexL;
	}
	key = ftok(filename, 'i');
	shmid = shmget(key, size, 0666);
	if(shmid < 0) {
		/* not shared */
		fprintf(stderr, "DB not shared, see kma_shm\n");
		exit(2);
	} else {
		/* found */
		dest->value_index = shmat(shmid, NULL, 0);
	}
	dest->value_index_l = (long unsigned *)(dest->value_index);
}

/* hashMap indexes */
void hashMap_index_initialize(struct hashMap_index *dest, int len) {
	
	dest->len = len;
	dest->size = len << 1;
	
	dest->index = calloc(dest->size, sizeof(int));
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
	
	if(dest->index) {
		free(dest->index);
	}
	if(dest->seq) {
		free(dest->seq);
	}
	free(dest);
}

int hashMap_index_get(struct hashMap_index *dest, long unsigned key) {
	
	unsigned index;
	int pos;
	
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

int hashMap_index_getDub(struct hashMap_index *dest, long unsigned key, const char *qseq, int q_len, int *next) {
	
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
				*next = score;
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
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(min < abs(pos) && abs(pos) < max && getKmer(dest->seq, abs(pos) - 1) == key) {
			return pos;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(min < abs(pos) && abs(pos) < max && getKmer(dest->seq, abs(pos) - 1) == key) {
				return pos;
			}
		}
	}
	
	return 0;
}

int hashMap_index_getDubPos(struct hashMap_index *dest, long unsigned key, int value) {
	
	unsigned index;
	int pos;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(pos == value) {
			return index;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos == value) {
				return index;
			}
		}
	}
	
	return -1;
}

int hashMap_index_getNextDubPos(struct hashMap_index *dest, long unsigned key, int min, int max, unsigned index) {
	
	int pos;
	
	min = -min;
	max = -max;
	
	while(++index < dest->size && (pos = dest->index[index]) != 0) {
		if(max < pos && pos < min && getKmer(dest->seq, -pos - 1) == key) {
			return index;
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(max < pos && pos < min && getKmer(dest->seq, -pos - 1) == key) {
				return index;
			}
		}
	}
	
	return -1;
}

int hashMap_index_getDub_bound(struct hashMap_index *dest, long unsigned key, const char *qseq, int q_len, int *next, int min, int max) {
	
	int i, index, pos, maxScore, score, mPos;
	maxScore = 0;
	mPos = 0;
	
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(pos < 0 && max < pos && pos < min && getKmer(dest->seq, (-1) - pos) == key) {
			//pos = kmersize - pos + 1;
			pos = (-1) - pos + kmersize;
			score = 0;
			for(i = kmersize; i < q_len && pos < -max && getNuc(dest->seq, pos) == qseq[i]; ++i, ++pos) {
				++score;
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
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos < 0 && max < pos && pos < min && getKmer(dest->seq, (-1) - pos) == key) {
				//pos = kmersize - pos + 1;
				pos = (-1) - pos + kmersize;
				score = 0;
				for(i = kmersize; i < q_len && pos < -max && getNuc(dest->seq, pos) == qseq[i]; ++i, ++pos) {
					++score;
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
	
	int pos, neg;
	unsigned index;
	
	if(key == 0) {
		/* likely undefined region */
		return;
	}
	
	neg = 1;
	++newpos;
	for(index = key % dest->size; index < dest->size && (pos = dest->index[index]) != 0; ++index) {
		if(pos > 0) {
			if(getKmer(dest->seq, pos - 1) == key) {
				dest->index[index] = -pos;
				neg = -1;
			}
		} else {
			if(getKmer(dest->seq, -pos - 1) == key) {
				neg = -1;
			}
		}
	}
	
	if(index == dest->size) {
		for(index = 0; (pos = dest->index[index]) != 0; ++index) {
			if(pos > 0) {
				if(getKmer(dest->seq, pos - 1) == key) {
					dest->index[index] = -pos;
					neg = -1;
				}
			} else {
				if(getKmer(dest->seq, -pos - 1) == key) {
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

struct hashMap_index *hashMap_index_build(FILE *seq, int len) {
	
	int i, end;
	struct hashMap_index *src;
	
	src = malloc(sizeof(struct hashMap_index));
	if(!src) {
		ERROR();
	}
	hashMap_index_initialize(src, len);
	
	fread(src->seq, sizeof(long unsigned), (src->len >> 5) + 1, seq);
	
	end = len - kmersize + 1;
	for(i = 0; i < end; ++i) {
		hashMap_index_add(src, getKmer(src->seq, i), i);
	}
	
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
		ERROR();
	}
}

void hashMap_kmers_CountIndex(struct hashMap_kmers *dest, long unsigned key) {
	
	unsigned index;
	struct hashTable_kmers *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		++dest->n;
		dest->table[index] = malloc(sizeof(struct hashTable_kmers));
		if(!dest->table[index]) {
			ERROR();
		}					
		node = dest->table[index];
		node->value = 1;
		node->key = key;
		node->next = 0;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(key == node->key) { // Keys match change value
				++node->value;
				return;
			} else if(node->next == 0) { // This chain is filled, create next
				++dest->n;
				node->next = malloc(sizeof(struct hashTable_kmers));
				if(!node->next) {
					ERROR();
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

struct hashMap_index * alignLoad_fly_build(FILE *seq_in, FILE *index_in, int len, long unsigned seq_index, long unsigned index_index) {
	
	/* move file pointer */
	fseek(seq_in, seq_index, SEEK_SET);
	
	return hashMap_index_build(seq_in, len);
}

struct hashMap_index * alignLoad_fly_build_mem(FILE *seq_in, FILE *index_in, int len, long unsigned seq_index, long unsigned index_index) {
	
	return hashMap_index_build(seq_in, len);
}

struct hashMap_index * alignLoad_fly_shm(FILE *seq_in, FILE *index_in, int len, long unsigned seq_index, long unsigned index_index) {
	
	struct hashMap_index *dest;
	
	dest = malloc(sizeof(struct hashMap_index));
	if(!dest) {
		ERROR();
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
		ERROR();
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
		exit(2);
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
		exit(2);
	} else {
		dest->seq = shmat(shmid, NULL, 0);
	}
	rewind(seq_in);
	templatefilename[file_len] = 0;
	
	return dest;
}
void alignClean(int template) {
	
	if(templates_index[template]) {
		hashMap_index_destroy(templates_index[template]);
		templates_index[template] = 0;
	}
}

void alignClean_shm(int template) {
	free(templates_index[template]);
	templates_index[template] = 0;
}

void load_DBs_KMA(char *templatefilename) {
	
	/* load DBs needed for KMA */
	int file_len, shmid;
	long unsigned i, j, file_size;
	FILE *DB_file;
	key_t key;
	
	/* allocate DBs */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	DB_file = sfopen(templatefilename, "rb");
	sfread(&DB_size, sizeof(int), 1, DB_file);
	contamination = DB_size;
	if(shm & 4) {
		key = ftok(templatefilename, 'l');
		shmid = shmget(key, DB_size * sizeof(int), 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared length\n");
			exit(2);
		} else {
			template_lengths = shmat(shmid, NULL, 0);
		}
	} else {
		template_lengths = malloc(DB_size * sizeof(int));
		if(!template_lengths) {
			ERROR();
		}
		/* load lengths */
		sfread(template_lengths, sizeof(int), DB_size, DB_file);
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	/* allocate pointers */
	template_names = malloc(DB_size * sizeof(char*));
	templates_index = calloc(DB_size, sizeof(struct hashMap_index*));
	alignment_scores = calloc(DB_size, sizeof(long unsigned));
	uniq_alignment_scores = calloc(DB_size, sizeof(long unsigned));
	if(!templates_index || !template_names || !alignment_scores || !uniq_alignment_scores) {
		ERROR();
	}
	
	/* load names */
	strcat(templatefilename, ".name");
	DB_file = sfopen(templatefilename, "rb");
	
	/* get size of file */
	fseek(DB_file, 0, SEEK_END);
	file_size = ftell(DB_file);
	
	/* load file */
	if(shm & 16) {
		key = ftok(templatefilename, 'n');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared name\n");
			ERROR();
		} else {
			template_names[0] = shmat(shmid, NULL, 0);
		}
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; ++i) {
			if(template_names[0][i] == 0) {
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	} else {
		rewind(DB_file);
		template_names[0] = malloc(file_size);
		if(!template_names[0]) {
			ERROR();
		}
		sfread(template_names[0], 1, file_size, DB_file);
		template_names[0][file_size - 1] = 0;
		
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; ++i) {
			if(template_names[0][i] == '\n') {
				template_names[0][i] = 0;
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
}

void load_DBs_Sparse(char *templatefilename) {
	
	/* load DBs needed for KMA */
	int file_len, shmid;
	long unsigned i, j, file_size;
	FILE *DB_file;
	key_t key;
	
	
	/* Open DB */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	templatefilename[file_len + 9] = 0;
	DB_file = sfopen(templatefilename, "rb");
	
	/* allocate DBs */
	fread(&DB_size, sizeof(int), 1, DB_file);
	contamination = DB_size;
	if(shm & 4) {
		fseek(DB_file, 0, SEEK_END);
		file_size = ftell(DB_file) - sizeof(int);
		key = ftok(templatefilename, 'l');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared length\n");
			exit(errno);
		} else {
			template_lengths = shmat(shmid, NULL, 0);
		}
		template_lengths += DB_size;
		template_ulengths = template_lengths + DB_size;
	} else {
		template_lengths = malloc(DB_size * sizeof(int));
		template_ulengths = malloc(DB_size * sizeof(int));
		if(!template_lengths || !template_ulengths) {
			ERROR();
		}
		/* load lengths */
		fseek(DB_file, DB_size * sizeof(int), SEEK_CUR);
		fread(template_lengths, sizeof(int), DB_size, DB_file);
		fread(template_ulengths, sizeof(int), DB_size, DB_file);	
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	/* load names */
	strcat(templatefilename, ".name");
	templatefilename[file_len + 5] = 0;
	DB_file = sfopen(templatefilename, "rb");
	
	/* get size of file */
	fseek(DB_file, 0, SEEK_END);
	file_size = ftell(DB_file);
	template_names = malloc(DB_size * sizeof(char*));
	if(!template_names) {
		ERROR();
	}
	
	/* load file */
	if(shm & 16) {
		key = ftok(templatefilename, 'n');
		shmid = shmget(key, file_size, 0666);
		if(shmid < 0) {
			fprintf(stderr, "No shared name\n");
			ERROR();
		} else {
			template_names[0] = shmat(shmid, NULL, 0);
		}
		for(i = 0, j = 2; j < DB_size; ++i) {
			if(template_names[0][i] == 0) {
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	} else {
		rewind(DB_file);
		template_names[0] = malloc(file_size);
		if(!template_names[0]) {
			ERROR();
		}
		fread(template_names[0], 1, file_size, DB_file);
		template_names[0][file_size - 1] = 0;
		
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < DB_size; ++i) {
			if(template_names[0][i] == '\n') {
				template_names[0][i] = 0;
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
}

/* METHOD SPECIFIC METHODS */
void print_ankers(int *out_Tem, struct compDNA *qseq, int rc_flag, struct qseqs *header) {
	
	int infoSize[6];
	
	infoSize[0] = qseq->seqlen;
	infoSize[1] = qseq->complen;
	infoSize[2] = qseq->N[0];
	infoSize[3] = rc_flag;
	infoSize[4] = *out_Tem;
	infoSize[5] = header->len;
	sfwrite(infoSize, sizeof(int), 6, stdout);
	sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, stdout);
	if(qseq->N[0]) {
		sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], stdout);
	}
	sfwrite(out_Tem + 1, sizeof(int), *out_Tem, stdout);
	sfwrite(header->seq, 1, header->len, stdout);
}

void print_ankers_spltDB(int *out_Tem, struct compDNA *qseq, int rc_flag, struct qseqs *header) {
	
	static unsigned target = 1, allIn = 0;
	static struct spltDBbuff **buffers = 0;
	int num, size, tNum, infoSize[7];
	unsigned char *buff;
	struct spltDBbuff *node, *nodeN, *node0;
	
	if(!buffers) {
		buffers = calloc(thread_num, sizeof(struct spltDBbuff));
		if(!buffers) {
			ERROR();
		}
	} else if(qseq == 0) {
		/* catch remains after last call */
		while(allIn) {
			/* find and write next target */
			node = *buffers;
			target = node->next->num;
			num = thread_num;
			tNum = 0;
			while(--num) {
				if(buffers[num] && buffers[num]->next->num < target) {
					tNum = num;
					node = buffers[num]->next;
					target = node->next->num;
				}
			}
			sfwrite(node->buff, 1, node->size, stdout);
			
			/* update cylinder */
			nodeN = node->prev;
			node0 = node->next;
			if(nodeN == node) {
				--allIn;
				++target;
				nodeN = 0;
			} else {
				nodeN->next = node0;
				node0->prev = nodeN;
			}
			free(node->buff);
			free(node);
			buffers[tNum] = nodeN;
		}
		sfwrite(&(unsigned){U_LIMIT}, sizeof(unsigned), 1, stdout);
		sfwrite(&(int){out_Tem[1] - 1}, sizeof(int), 1, stdout);
		return;
	}
	
	num = out_Tem[-1];
	infoSize[0] = num;
	infoSize[1] = qseq->seqlen;
	infoSize[2] = qseq->complen;
	infoSize[3] = qseq->N[0];
	infoSize[4] = rc_flag;
	infoSize[5] = *out_Tem;
	infoSize[6] = header->len;
	
	if(num == target || thread_num == 1) {
		sfwrite(infoSize, sizeof(int), 7, stdout);
		sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, stdout);
		if(qseq->N[0]) {
			sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], stdout);
		}
		sfwrite(out_Tem + 1, sizeof(int), *out_Tem, stdout);
		sfwrite(header->seq, 1, header->len, stdout);
		++target;
	} else {
		
		node = smalloc(sizeof(struct spltDBbuff));
		size = header->len + sizeof(int) * (7 + qseq->N[0] + out_Tem[0]) + sizeof(long unsigned) * qseq->complen;
		buff = smalloc(size);
		node->num = num;
		node->size = size;
		node->buff = buff;
		memcpy(buff, infoSize, (size = sizeof(int) * 7));
		buff += size;
		memcpy(buff, qseq->seq, (size = sizeof(long unsigned) * qseq->complen));
		buff += size;
		memcpy(buff, qseq->N + 1, (size = sizeof(int) * qseq->N[0]));
		buff += size;
		memcpy(buff, out_Tem + 1, (size = sizeof(int) * (*out_Tem)));
		buff += size;
		memcpy(buff, header->seq, header->len);
		
		/* find cylinder piece */
		tNum = out_Tem[*out_Tem + 2];
		/* extend cylinder piece */
		if((nodeN = buffers[tNum])) {
			node0 = nodeN->next;
			node->prev = nodeN;
			nodeN->next = node;
			node->next = node0;
			node0->prev = node;
		} else {
			node->next = node;
			node->prev = node;
			++allIn;
		}
		buffers[tNum] = node;
	}
	
	/* collect remains */
	while(allIn == thread_num) {
		/* find and write next target */
		node = *buffers;
		target = node->next->num;
		num = thread_num;
		tNum = 0;
		while(--num) {
			if(buffers[num]->next->num < target) {
				tNum = num;
				node = buffers[num]->next;
				target = node->next->num;
			}
		}
		sfwrite(node->buff, 1, node->size, stdout);
		
		/* update cylinder */
		nodeN = node->prev;
		node0 = node->next;
		if(nodeN == node) {
			--allIn;
			++target;
			nodeN = 0;
		} else {
			nodeN->next = node0;
			node0->prev = nodeN;
		}
		free(node->buff);
		free(node);
		buffers[tNum] = nodeN;
	}
}

void print_ankers_Sparse(int *out_Tem, struct compDNA *qseq, int rc_flag, struct qseqs *header) {
	
	int infoSize[6];
	
	infoSize[0] = qseq->seqlen;
	infoSize[1] = qseq->complen;
	infoSize[2] = qseq->N[0];
	infoSize[3] = -(abs(rc_flag));
	infoSize[4] = *out_Tem;
	infoSize[5] = header->len;
	sfwrite(infoSize, sizeof(int), 6, stdout);
	
	sfwrite(qseq->seq, sizeof(long unsigned), qseq->complen, stdout);
	sfwrite(qseq->N + 1, sizeof(int), qseq->N[0], stdout);
	sfwrite(out_Tem + 1, sizeof(int), *out_Tem, stdout);
	sfwrite(header->seq, 1, header->len, stdout);
	
}

int get_ankers(int *out_Tem, struct compDNA *qseq, struct qseqs *header, FILE *inputfile) {
	
	int infoSize[6];
	
	if(fread(infoSize, sizeof(int), 6, inputfile)) {
		qseq->seqlen = infoSize[0];
		qseq->complen = infoSize[1];
		*out_Tem = infoSize[4];
		header->len = infoSize[5];
		
		/* reallocate */
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
		}
		
		qseq->N[0] = infoSize[2];
		if(header->size <= header->len) {
			free(header->seq);
			header->size = header->len << 1;
			header->seq = malloc(header->size);
			if(!header->seq) {
				ERROR();
			}
		}
		
		fread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		fread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		fread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		fread(header->seq, 1, header->len, inputfile);
	} else {
		infoSize[3] = 0;
	}
	
	/* return score */
	return infoSize[3];
}

int get_ankers_spltDB(int *infoSize, int *out_Tem, struct compDNA *qseq, struct qseqs *header, FILE *inputfile) {
	
	int num;
	
	if(qseq) {
		qseq->seqlen = infoSize[0];
		qseq->complen = infoSize[1];
		*out_Tem = infoSize[4];
		header->len = infoSize[5];
		
		/* reallocate */
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
		}
		
		qseq->N[0] = infoSize[2];
		if(header->size < header->len) {
			free(header->seq);
			header->size = header->len;
			header->seq = malloc(header->size);
			if(!header->seq) {
				ERROR();
			}
		}
		
		fread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		fread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		fread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		fread(header->seq, 1, header->len, inputfile);
	} else {
		/* in the case of equally well scoring DBs */
		fseek(inputfile, infoSize[1] * sizeof(long unsigned) + infoSize[2] * sizeof(int), SEEK_CUR);
		*out_Tem = infoSize[4];
		fread(out_Tem + 1, sizeof(int), *out_Tem, inputfile);
		fseek(inputfile, infoSize[5], SEEK_CUR);
	}
	
	/* get info for next read */
	fread(&num, sizeof(int), 1, inputfile);
	fread(infoSize, sizeof(int), 6, inputfile);
	return num;
}

void print_anker_pairs(int *out_Tem, struct compDNA *qseq, struct compDNA *qseq_r, int score, int score_r, struct qseqs *header, struct qseqs *header_r) {
	
	dumpComp(qseq, stdout);
	sfwrite(&(int){score}, sizeof(int), 1, stdout);
	sfwrite(out_Tem, sizeof(int), *out_Tem + 1, stdout);
	sfwrite(&(int){-header->len}, sizeof(int), 1, stdout);
	sfwrite(header->seq, 1, header->len, stdout);
	
	dumpComp(qseq_r, stdout);
	sfwrite(&(int){score_r}, sizeof(int), 1, stdout);
	sfwrite(&(int){-header_r->len}, sizeof(int), 1, stdout);
	sfwrite(header_r->seq, 1, header_r->len, stdout);
}

void ankerAndClean(int *regionTemplates, int *Score, int *Score_r, unsigned **VF_scores, unsigned **VR_scores, int *tmpNs, struct compDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, struct qseqs *header) {
	
	int k, l, bestHitsCov, template;
	unsigned *values, n, SU;
	short unsigned *values_s;
	double thisCov, bestCov;
	struct compDNA tmpQseq;
	
	if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	/* here */
	// make sure cuts isn't random seeds
	
	bestHitsCov = *regionTemplates;
	bestCov = 0;
	for(k = 1; k <= *regionTemplates; ++k) {
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
	
	for(k = start_cut; k <= end_cut; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				n = *values_s;
				for(l = 1; l <= n; ++l) {
					template = values_s[l];
					if(Score[template] != bestScore && template != contamination) {
						thisCov = 1.0 * Score[template] / template_lengths[template];
						if(thisCov > bestCov) {
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = template;
						} else if(thisCov == bestCov) {
							++bestHitsCov;
							regionTemplates[bestHitsCov] = template;
						}
					}
					Score[template]--;
				}
			} else {
				values = VF_scores[k];
				n = *values;
				for(l = 1; l <= n; ++l) {
					template = values[l];
					if(Score[template] != bestScore && template != contamination) {
						thisCov = 1.0 * Score[template] / template_lengths[template];
						if(thisCov > bestCov) {
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = template;
						} else if(thisCov == bestCov) {
							++bestHitsCov;
							regionTemplates[bestHitsCov] = template;
						}
					}
					Score[template]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				n = *values_s;
				for(l = 1; l <= n; ++l) {
					template = values_s[l];
					if(Score_r[template] != bestScore && template != contamination) {
						thisCov = 1.0 * Score_r[template] / template_lengths[template];
						if(thisCov > bestCov) {
							HIT = -1;
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = -template;
						} else if(thisCov == bestCov) {
							HIT = -1;
							++bestHitsCov;
							regionTemplates[bestHitsCov] = -template;
						}
					}
					Score_r[template]--;
				}
			} else {
				values = VR_scores[k];
				n = *values;
				for(l = 1; l <= n; ++l) {
					template = values[l];
					if(Score_r[template] != bestScore && template != contamination) {
						thisCov = 1.0 * Score_r[template] / template_lengths[template];
						if(thisCov > bestCov) {
							HIT = -1;
							bestCov = thisCov;
							bestHitsCov = *regionTemplates + 1;
							regionTemplates[bestHitsCov] = -template;
						} else if(thisCov == bestCov) {
							HIT = -1;
							++bestHitsCov;
							regionTemplates[bestHitsCov] = -template;
						}
					}
					Score_r[template]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	*regionTemplates = bestHitsCov;
	
	/* clear nearest templates on both sides of match */
	for(k = ((start_cut - 92) < 0) ? 0 : (start_cut - 92); k < start_cut; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
			}
			VR_scores[k] = 0;
		}
	}
	for(k = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92); k > end_cut; --k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
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
	
	for(k = 1, l = 0; k < qseq->N[0]; ++k) {
		if(start_cut <= qseq->N[k]) {
			++l;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				--l;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	--tmpQseq.seqlen;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		--tmpQseq.seqlen;
		--l;
	}
	++tmpQseq.seqlen;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	lock(excludeOut);
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header);
	unlock(excludeOut);
}

void ankerAndClean_MEM(int *regionTemplates, int *Score, int *Score_r, unsigned **VF_scores, unsigned **VR_scores, int *tmpNs, struct compDNA *qseq, int HIT, int bestScore, int start_cut, int end_cut, struct qseqs *header) {
	
	int k, l, SU;
	unsigned *values;
	short unsigned *values_s;
	struct compDNA tmpQseq;
	
	if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	/* clean up scores */
	start_cut = ((start_cut - 92) < 0) ? 0 : (start_cut - 92);
	end_cut = ((end_cut + 92) > qseq->seqlen) ? qseq->seqlen : (end_cut + 92);
	for(k = start_cut; k < end_cut; ++k) {
		if(VF_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VF_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]]--;
				}
			} else {
				values = VF_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score[values[l]]--;
				}
			}
			VF_scores[k] = 0;
		}
		if(VR_scores[k]) {
			if(SU) {
				values_s = (short unsigned *) VR_scores[k];
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]]--;
				}
			} else {
				values = VR_scores[k];
				l = (*values) + 1;
				while(--l) {
					Score_r[values[l]]--;
				}
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
	
	for(k = 1, l = 0; k < qseq->N[0]; ++k) {
		if(start_cut <= qseq->N[k]) {
			++l;
			tmpQseq.N[l] = qseq->N[k] - start_cut;
			if(tmpQseq.N[l] >= tmpQseq.seqlen) {
				--l;
				k = qseq->N[0];
			}
		}
	}
	
	/* trim trailing gaps */
	--tmpQseq.seqlen;
	while(tmpQseq.N[l] == tmpQseq.seqlen && l != 0) {
		--tmpQseq.seqlen;
		--l;
	}
	++tmpQseq.seqlen;
	tmpQseq.complen = (tmpQseq.seqlen >> 5) + 1;
	tmpQseq.N[0] = l;
	
	lock(excludeOut);
	deConPrintPtr(regionTemplates, &tmpQseq, HIT * bestScore, header);
	unlock(excludeOut);
}

void deConPrint(int *out_Tem, struct compDNA *qseq, int rc_flag, struct qseqs *header) {
	
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
		printPtr(out_Tem, qseq, rc_flag, header);
	}
}

void deConPrintPair(int *out_Tem, struct compDNA *qseq, int bestScore, struct qseqs *header, struct compDNA *qseq_r, int bestScore_r, struct qseqs *header_r) {
	
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
		contPos = *out_Tem;
		*out_Tem = 0;
		printPtr(out_Tem, qseq, bestScore, header);
		*out_Tem = contPos;
		printPtr(out_Tem, qseq_r, bestScore_r, header_r);
	}
}

void printPair(int *out_Tem, struct compDNA *qseq, int bestScore, struct qseqs *header, struct compDNA *qseq_r, int bestScore_r, struct qseqs *header_r) {
	
	int contPos;
	
	contPos = *out_Tem;
	*out_Tem = 0;
	printPtr(out_Tem, qseq, bestScore, header);
	*out_Tem = contPos;
	out_Tem[-1]++;
	printPtr(out_Tem, qseq_r, bestScore_r, header_r);
}

FILE * printFrags(struct frag **alignFrags) {
	
	int i;
	FILE *OUT;
	struct frag *alignFrag, *next;
	
	OUT = tmpfile();
	if(!OUT) {
		fprintf(stderr, "Could not create tmp files.\n");
		ERROR();
	}
	setvbuf(OUT, NULL, _IOFBF, CHUNK);
	for(i = 0; i < DB_size; ++i) {
		if(alignFrags[i] != 0) {
			for(alignFrag = alignFrags[i]; alignFrag != 0; alignFrag = next) {
				next = alignFrag->next;
				
				sfwrite(&i, sizeof(int), 1, OUT);
				sfwrite(alignFrag->buffer, sizeof(int), 6, OUT);
				sfwrite(alignFrag->qseq, 1, alignFrag->buffer[0], OUT);
				sfwrite(alignFrag->header, 1, alignFrag->buffer[5], OUT);
				
				free(alignFrag->qseq);
				free(alignFrag->header);
				free(alignFrag);
			}
			alignFrags[i] = 0;
		}
	}
	sfwrite(&(int){-1}, sizeof(int), 1, OUT);
	fflush(OUT);
	rewind(OUT);
	
	return OUT;
}

void bootFsa(struct qseqs *header, struct qseqs *qseq, struct compDNA *compressor) {
	
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
		
		sfwrite(buffer, sizeof(int), 4, stdout);
		sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, stdout);
		sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], stdout);
		sfwrite((header->seq + 1), 1, header->len, stdout);
		resetComp(compressor);
	}
	
	compDNA(compressor, qseq->seq + i, qseq->len - i);
	/* print */
	buffer[0] = compressor->seqlen;
	buffer[1] = compressor->complen;
	buffer[2] = compressor->N[0];
	sfwrite(buffer, sizeof(int), 4, stdout);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, stdout);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], stdout);
	sfwrite((header->seq + 1), 1, header->len, stdout);
	resetComp(compressor);
}

void printFsa(struct qseqs *header, struct qseqs *qseq, struct compDNA *compressor) {
	
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
	
	sfwrite(buffer, sizeof(int), 4, stdout);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, stdout);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], stdout);
	sfwrite((header->seq + 1), 1, header->len, stdout);
	
	
	resetComp(compressor);
}

void printFsaMt1(struct qseqs *header, struct qseqs *qseq, struct compDNA *compressor) {
	
	template_lengths[1] = qseq->len;
	template_lengths[6] = header->len;
	
	sfwrite(template_lengths, sizeof(int), 7, stdout);
	sfwrite(qseq->seq, 1, qseq->len, stdout);
	sfwrite(header->seq + 1, 1, header->len, stdout);
}

void printFsa_pairMt1(struct qseqs *header, struct qseqs *qseq, struct qseqs *header_r, struct qseqs *qseq_r, struct compDNA *compressor) {
	
	printFsaMt1(header, qseq, compressor);
	printFsaMt1(header_r, qseq_r, compressor);
	
}

void printFsa_pair(struct qseqs *header, struct qseqs *qseq, struct qseqs *header_r, struct qseqs *qseq_r, struct compDNA *compressor) {
	
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
	
	sfwrite(buffer, sizeof(int), 4, stdout);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, stdout);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], stdout);
	sfwrite((header->seq + 1), 1, header->len, stdout);
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
	
	sfwrite(buffer, sizeof(int), 4, stdout);
	sfwrite(compressor->seq, sizeof(long unsigned), compressor->complen, stdout);
	sfwrite(compressor->N + 1, sizeof(int), compressor->N[0], stdout);
	sfwrite((header_r->seq + 1), 1, header_r->len, stdout);
	resetComp(compressor);
}

int loadFsa(struct compDNA *qseq, struct qseqs *header, FILE *inputfile) {
	
	int buffer[4];
	
	if(fread(buffer, sizeof(int), 4, inputfile)) {
		qseq->seqlen = buffer[0];
		qseq->complen = buffer[1];
		/* if pair, header->len < 0 */
		header->len = abs(buffer[3]);
		
		if(qseq->size <= qseq->seqlen) {
			free(qseq->N);
			free(qseq->seq);
			if(qseq->seqlen & 31) {
				qseq->size = (qseq->seqlen >> 5) + 1;
				qseq->size <<= 6;
			} else {
				qseq->size = qseq->seqlen << 1;
			}
			
			qseq->seq = calloc(qseq->size >> 5, sizeof(long unsigned));
			qseq->N = malloc((qseq->size + 1) * sizeof(int));
			if(!qseq->seq || !qseq->N) {
				ERROR();
			}
		}
		qseq->N[0] = buffer[2];
		
		if(header->size <= header->len) {
			header->size = header->len << 1;
			free(header->seq);
			header->seq = malloc(header->size);
			if(!header->seq) {
				ERROR();
			}
		}
		fread(qseq->seq, sizeof(long unsigned), qseq->complen, inputfile);
		fread(qseq->N + 1, sizeof(int), qseq->N[0], inputfile);
		fread(header->seq, 1, header->len, inputfile);
	} else {
		qseq->seqlen = 0;
		return 0;
	}
	
	return buffer[3];
}

int get_kmers_for_pair(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, int *extendScore) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, l, rc, end, HIT, gaps, score, Ms, MMs, Us, W1s, template, SU;
	int hitCounter, bestSeqCount, *bests, *Scores;
	unsigned *values, *last, n;
	short unsigned *values_s;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	bests = bestTemplates;
	Scores = Score;
	bestSeqCount = 0;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			bests = bestTemplates_r;
			Scores = Score_r;
			comp_rc(qseq);
		}
		/* Make quick check of the qseq */
		HIT = exhaustive;
		hitCounter = 0;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		/* If deltamer qseq hits, then continue */
		if(HIT) {
			/* Scan the deltamer exhaustively, and collect scores in Score*/
			last = 0;
			gaps = 0;
			HIT = 0;
			Ms = 0;
			MMs = 0;
			Us = 0;
			W1s = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
						if(values == last) {
							if(kmersize < gaps) {
								Ms += kmersize;
								gaps -= kmersize;
								if(gaps) {
									/* go for best scenario */
									if(gaps == 1) {
										MMs += 2;
									} else {
										gaps -= 2;
										if((MM << 1) + gaps * M < 0) {
											Ms += gaps;
											MMs += 2;
										}
									}
								} else {
									++MMs;
								}
							} else if (gaps) {
								--gaps;
								++W1s;
								Us += gaps;
							} else {
								++Ms;
							}
							HIT = j;
							gaps = 0;
						} else {
							if(last) {
								if(HIT) {
									HIT += kmersize;
								} else {
									HIT = j + kmersize;
								}
								score = Ms * M + MMs * MM + Us * U + W1s * W1;
								if(SU) {
									values_s = (short unsigned *) last;
									l = (*values_s) + 1;
									while(--l) {
										Scores[(template = values_s[l])] += score;
										extendScore[template] = HIT;
									}
									
									score = kmersize * M;
									MMs = MM << 1;
									values_s = (short unsigned *) values;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										if(j < extendScore[(template = values_s[l])]) {
											if(extendScore[template] == HIT) {
												Scores[template] += M;
											} else {
												gaps = extendScore[template] - j - 1;
												Scores[template] += (W1 + gaps * U);
											}
										} else if(Scores[template] != 0) {
											Scores[template] += score;
											if((gaps = extendScore[template] - j)) {
												if(gaps == 1) {
													Scores[template] += MMs;
												} else {
													gaps -= 2;
													if((Ms = MMs + gaps * M) < 0) {
														Scores[template] += Ms;
													}
												}
											} else {
												Scores[template] += MM;
											}
										} else {
											Scores[template] = score;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								} else {
									l = (*last) + 1;
									while(--l) {
										Scores[(template = last[l])] += score;
										extendScore[template] = HIT;
									}
									
									score = kmersize * M;
									MMs = MM << 1;
									n = *values;
									for(l = 1; l <= n; ++l) {
										if(j < extendScore[(template = values[l])]) {
											if(extendScore[template] == HIT) {
												Scores[template] += M;
											} else {
												gaps = extendScore[template] - j - 1;
												Scores[template] += (W1 + gaps * U);
											}
										} else if(Scores[template] != 0) {
											Scores[template] += score;
											if((gaps = extendScore[template] - j)) {
												if(gaps == 1) {
													Scores[template] += MMs;
												} else {
													gaps -= 2;
													if((Ms = MMs + gaps * M) < 0) {
														Scores[template] += Ms;
													}
												}
											} else {
												Scores[template] += MM;
											}
										} else {
											Scores[template] = score;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								}
							} else if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								Ms = kmersize * M;
								for(l = 1; l <= n; ++l) {
									Scores[(template = values_s[l])] = Ms;
									bests[l] = template;
								}
								*bests = n;
							} else {
								n = *values;
								Ms = kmersize * M;
								for(l = 1; l <= n; ++l) {
									Scores[(template = values[l])] = Ms;
									bests[l] = template;
								}
								*bests = n;
							}
							
							HIT = 0;
							gaps = 0;
							Ms = 0;
							MMs = 0;
							Us = 0;
							W1s = 0;
							last = values;
						}
						++hitCounter;
					} else {
						++gaps;
					}
				}
				j = qseq->N[i] + 1;
			}
			
			if(last) {
				score = Ms * M + MMs * MM + Us * U + W1s * W1;
				if(SU) {
					values_s = (short unsigned *) last;
					l = (*values_s) + 1;
					while(--l) {
						Scores[values_s[l]] += score;
					}
				} else {
					l = (*last) + 1;
					while(--l) {
						Scores[last[l]] += score;
					}
				}
				for(l = *bests; l != 0; --l) {
					extendScore[(template = bests[l])] = 0;
					if(Scores[template] < 0) {
						Scores[template] = 0;
					}
				}
			}
			
			if(bestSeqCount < hitCounter) {
				bestSeqCount = hitCounter;
			}
		}
		qseq->N[0]--;
	}
	
	return bestSeqCount;
}

int get_kmers_for_pair_count(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, int *extendScore) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, l, rc, end, HIT, template, hitCounter, bestSeqCount, reps, SU;
	int *bests, *Scores;
	unsigned *values, *last, n;
	short unsigned *values_s;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	bests = bestTemplates;
	Scores = Score;
	bestSeqCount = 0;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			bests = bestTemplates_r;
			Scores = Score_r;
			comp_rc(qseq);
		}
		/* Make quick check of the qseq */
		HIT = exhaustive;
		hitCounter = 0;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		/* If deltamer qseq hits, then continue */
		if(HIT) {
			/* Scan the deltamer exhaustively, and collect scores in Score*/
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) last;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										if(Scores[(template = values_s[l])]) {
											Scores[template] += reps;
										} else {
											Scores[template] = reps;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								} else {
									n = *last;
									for(l = 1; l <= n; ++l) {
										if(Scores[(template = last[l])]) {
											Scores[template] += reps;
										} else {
											Scores[template] = reps;
											bests[0]++;
											bests[*bests] = template;
										}
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
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						if(Scores[(template = values_s[l])]) {
							Scores[template] += reps;
						} else {
							Scores[template] = reps;
							bests[0]++;
							bests[*bests] = template;
						}
					}
				} else {
					n = *last;
					for(l = 1; l <= n; ++l) {
						if(Scores[(template = last[l])]) {
							Scores[template] += reps;
						} else {
							Scores[template] = reps;
							bests[0]++;
							bests[*bests] = template;
						}
					}
				}
				hitCounter += reps;
			}
			
			reps = 0;
			if(bestSeqCount < hitCounter) {
				bestSeqCount = hitCounter;
			}
		}
		qseq->N[0]--;
	}
	
	return bestSeqCount;
}

int get_kmers_for_pair_Sparse(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, int *extendScore) {
	
	int i, j, k, l, n, end, rc, prefix_len, template, hitCounter, HIT, SU;
	int n_kmers, gaps, Ms, MMs, Us, W1s, score, reps, *bests, *Scores;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned prefix;
	
	if(qseq->seqlen < kmersize) {
		return 0;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	if(*extendScore) {
		bests = bestTemplates_r;
		Scores = Score_r;
	} else {
		bests = bestTemplates;
		Scores = Score;
	}
	
	*extendScore = 0;
	*bestTemplates = 0;
	*bestTemplates_r = 0;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	hitCounter = 0;
	n_kmers = 0;
	end = qseq->seqlen;
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			if(rc) {
				comp_rc(qseq);
			}
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize - prefix_len + 1;
				while(j < end) {
					if(getSmer(qseq->seq, j) == prefix) {
						if((values = hashMap_get(templates, getKmer(qseq->seq, j + prefix_len)))) {
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								for(k = 1; k <= n; ++k) {
									if(Scores[(template = values_s[k])] == 0) {
										(*bests)++;
										bests[*bests] = template;
									}
									Scores[template]++;
								}
							} else {
								n = *values;
								for(k = 1; k <= n; ++k) {
									if(Scores[(template = values[k])] == 0) {
										(*bests)++;
										bests[*bests] = template;
									}
									Scores[template]++;
								}
							}
							++hitCounter;
						}
						++n_kmers;
					}
					++j;
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		if(hitCounter) {
			hitCounter *= (((qseq->seqlen - kmersize + 1) << 1) / n_kmers);
		}
	} else if(countK) {
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			last = 0;
			reps = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
						if(values == last) {
							++reps;
						} else {
							if(last) {
								if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										if(Scores[(template = values_s[l])]) {
											Scores[template] += reps;
										} else {
											Scores[template] = reps;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								} else {
									n = *values;
									for(l = 1; l <= n; ++l) {
										if(Scores[(template = last[l])]) {
											Scores[template] += reps;
										} else {
											Scores[template] = reps;
											bests[0]++;
											bests[*bests] = template;
										}
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
				if(SU) {
					values_s = (short unsigned *) last;
					n = *values_s;
					for(l = 1; l <= n; ++l) {
						if(Scores[(template = values_s[l])]) {
							Scores[template] += reps;
						} else {
							Scores[template] = reps;
							bests[0]++;
							bests[*bests] = template;
						}
					}
				} else {
					n = *last;
					for(l = 1; l <= n; ++l) {
						if(Scores[(template = last[l])]) {
							Scores[template] += reps;
						} else {
							Scores[template] = reps;
							bests[0]++;
							bests[*bests] = template;
						}
					}
				}
				hitCounter += reps;
			}
		}
		qseq->N[0]--;
	} else {
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			last = 0;
			gaps = 0;
			HIT = 0;
			Ms = 0;
			MMs = 0;
			Us = 0;
			W1s = 0;
			j = 0;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize + 1;
				for(;j < end; ++j) {
					if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
						if(values == last) {
							if(kmersize < gaps) {
								Ms += kmersize;
								gaps -= kmersize;
								if(gaps) {
									/* go for best scenario */
									if(gaps == 1) {
										MMs += 2;
									} else {
										gaps -= 2;
										if((MM << 1) + gaps * M < 0) {
											Ms += gaps;
											MMs += 2;
										}
									}
								} else {
									++MMs;
								}
							} else if (gaps) {
								--gaps;
								++W1s;
								Us += gaps;
							} else {
								++Ms;
							}
							HIT = j;
							gaps = 0;
						} else {
							if(last) {
								if(HIT) {
									HIT += kmersize;
								} else {
									HIT = j + kmersize;
								}
								score = Ms * M + MMs * MM + Us * U + W1s * W1;
								if(SU) {
									values_s = (short unsigned *) last;
									l = (*values_s) + 1;
									while(--l) {
										Scores[(template = values_s[l])] += score;
										extendScore[template] = HIT;
									}
								} else {
									l = (*last) + 1;
									while(--l) {
										Scores[(template = last[l])] += score;
										extendScore[template] = HIT;
									}
								}
								score = kmersize * M;
								MMs = MM << 1;
								if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s;
									for(l = 1; l <= n; ++l) {
										template = values_s[l];
										if(j < extendScore[template]) {
											if(extendScore[template] == HIT) {
												Scores[template] += M;
											} else {
												gaps = extendScore[template] - j - 1;
												Scores[template] += (W1 + gaps * U);
											}
										} else if(Scores[template] != 0) {
											Scores[template] += score;
											if((gaps = extendScore[template] - j)) {
												if(gaps == 1) {
													Scores[template] += MMs;
												} else {
													gaps -= 2;
													if((Ms = MMs + gaps * M) < 0) {
														Scores[template] += Ms;
													}
												}
											} else {
												Scores[template] += MM;
											}
										} else {
											Scores[template] = score;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								} else {
									n = *values;
									for(l = 1; l <= n; ++l) {
										template = values[l];
										if(j < extendScore[template]) {
											if(extendScore[template] == HIT) {
												Scores[template] += M;
											} else {
												gaps = extendScore[template] - j - 1;
												Scores[template] += (W1 + gaps * U);
											}
										} else if(Scores[template] != 0) {
											Scores[template] += score;
											if((gaps = extendScore[template] - j)) {
												if(gaps == 1) {
													Scores[template] += MMs;
												} else {
													gaps -= 2;
													if((Ms = MMs + gaps * M) < 0) {
														Scores[template] += Ms;
													}
												}
											} else {
												Scores[template] += MM;
											}
										} else {
											Scores[template] = score;
											bests[0]++;
											bests[*bests] = template;
										}
									}
								}
							} else if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								Ms = kmersize * M;
								for(l = 1; l <= n; ++l) {
									Scores[(template = values_s[l])] = Ms;
									bests[l] = template;
								}
								*bests = n;
							} else {
								n = *values;
								Ms = kmersize * M;
								for(l = 1; l <= n; ++l) {
									Scores[(template = values[l])] = Ms;
									bests[l] = template;
								}
								*bests = n;
							}
							HIT = 0;
							gaps = 0;
							Ms = 0;
							MMs = 0;
							Us = 0;
							W1s = 0;
							last = values;
						}
						++hitCounter;
					} else {
						++gaps;
					}
				}
				j = qseq->N[i] + 1;
			}
			if(last) {
				score = Ms * M + MMs * MM + Us * U + W1s * W1;
				if(SU) {
					values_s = (short unsigned *) last;
					l = (*values_s) + 1;
					while(--l) {
						Scores[values_s[l]] += score;
					}
				} else {
					l = (*last) + 1;
					while(--l) {
						Scores[last[l]] += score;
					}
				}
				for(l = *bests; l != 0; --l) {
					extendScore[(template = bests[l])] = 0;
					if(Scores[template] < 0) {
						Scores[template] = 0;
					}
				}
			}
		}
		qseq->N[0]--;
	}
	
	return hitCounter;
}

void getFirstForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, bestHits;
	
	bestHits = 0;
	
	for(i = 1; i <= *bestTemplates; ++i) {
		++bestHits;
		regionTemplates[bestHits] = bestTemplates[i];
		regionScores[bestHits] = Score[bestTemplates[i]];
		Score[bestTemplates[i]] = 0;
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		++bestHits;
		regionTemplates[bestHits] = -bestTemplates_r[i];
		regionScores[bestHits] = Score_r[bestTemplates_r[i]];
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
}

int getSecondForce(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, score, bestHits, bestScore;
	
	bestHits = 0;
	bestScore = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			if((score = Score[regionTemplates[i]])) {
				score += regionScores[i];
				if(bestScore < score) {
					bestScore = score;
					bestHits = 1;
					regionTemplates[bestHits] = regionTemplates[i];
				} else if(bestScore == score) {
					++bestHits;
					regionTemplates[bestHits] = regionTemplates[i];
				}
			}
		} else {
			if((score = Score_r[-regionTemplates[i]])) {
				score += regionScores[i];
				if(bestScore < score) {
					bestScore = score;
					bestHits = 1;
					regionTemplates[bestHits] = regionTemplates[i];
				} else if(bestScore == score) {
					++bestHits;
					regionTemplates[bestHits] = regionTemplates[i];
				}
			}
		}
	}
	*regionTemplates = bestHits;
	
	for(i = *bestTemplates; i != 0; --i) {
		Score[bestTemplates[i]] = 0;
	}
	for(i = *bestTemplates_r; i != 0; --i) {
		Score_r[bestTemplates_r[i]] = 0;
	}
	
	return bestScore;
}

int getFirstPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores) {
	
	int i, score, bestScore, bestHits;
	
	bestScore = 0;
	bestHits = 0;
	/* save forward matches */
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore < (score = Score[bestTemplates[i]])) {
			bestScore = score;
		}
		++bestHits;
		regionTemplates[bestHits] = bestTemplates[i];
		regionScores[bestHits] = score;
		Score[bestTemplates[i]] = 0;
	}
	
	/* save reverse matches */
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore < (score = Score_r[bestTemplates_r[i]])) {
			bestScore = score;
		}
		++bestHits;
		regionTemplates[bestHits] = -bestTemplates_r[i];
		regionScores[bestHits] = score;
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getSecondPen(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, int bestScore) {
	
	int i, score, bestScore_r, compScore, bestHits;
	
	/* get best scoring tempates */
	bestScore_r = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < Score[bestTemplates[i]]) {
			bestScore_r = Score[bestTemplates[i]];
		}
	}
	bestHits = *bestTemplates;
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < Score_r[bestTemplates_r[i]]) {
			bestScore_r = Score_r[bestTemplates_r[i]];
		}
		++bestHits;
		bestTemplates[bestHits] = -bestTemplates_r[i];
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	if(bestScore_r) {
		compScore = bestScore + bestScore_r - PE;
		compScore = MAX(0, compScore);
		for(i = 1; i <= *regionTemplates; ++i) {
			if(0 < regionTemplates[i]) {
				/* we got one */
				if(0 < (score = Score_r[regionTemplates[i]])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
						bestHits = 1;
						regionTemplates[bestHits] = regionTemplates[i];
					} else if(compScore == score) {
						++bestHits;
						regionTemplates[bestHits] = regionTemplates[i];
					}
				}
			} else {
				/* we got one */
				if(0 < (score = Score[-regionTemplates[i]])) {
					score += regionScores[i];
					if(compScore < score) {
						compScore = score;
						bestHits = 1;
						regionTemplates[bestHits] = regionTemplates[i];
					} else if(compScore == score) {
						++bestHits;
						regionTemplates[bestHits] = regionTemplates[i];
					}
				}
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
		/* clear scores */
		for(i = *bestTemplates; i != 0; --i) {
			if(0 < bestTemplates[i]) {
				Score[bestTemplates[i]] = 0;
			} else {
				Score_r[-bestTemplates[i]] = 0;
			}
		}
	} else {
		/* get bestHits from each as SE */
		for(i = 1; i <= *regionTemplates; ++i) {
			if(bestScore == regionScores[i]) {
				++bestHits;
				regionTemplates[bestHits] = regionTemplates[i];
			}
		}
		*regionTemplates = bestHits;
		
		bestHits = 0;
		for(i = 1; i <= *bestTemplates; ++i) {
			if(0 < bestTemplates[i]) {
				if(bestScore_r == Score[bestTemplates[i]]) {
					++bestHits;
					bestTemplates[bestHits] = bestTemplates[i];
				}
				Score[bestTemplates[i]] = 0;
			} else {
				if(bestScore_r == Score_r[-bestTemplates[i]]) {
					++bestHits;
					bestTemplates[bestHits] = bestTemplates[i];
				}
				Score_r[-bestTemplates[i]] = 0;
			}
		}
		*bestTemplates = bestHits;
	}
	
	return bestScore_r;
}

int getF_Best(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	int i, score, bestScore, bestHits;
	
	bestScore = 0;
	bestHits = 0;
	
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore < (score = Score[bestTemplates[i]])) {
			bestScore = score;
			bestHits = 1;
			regionTemplates[bestHits] = bestTemplates[i];
		} else if(bestScore == score) {
			++bestHits;
			regionTemplates[bestHits] = bestTemplates[i];
		}
		Score[bestTemplates[i]] = 0;
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore < (score = Score_r[bestTemplates_r[i]])) {
			bestScore = score;
			bestHits = 1;
			regionTemplates[bestHits] = -bestTemplates_r[i];
		} else if(bestScore == score) {
			++bestHits;
			regionTemplates[bestHits] = -bestTemplates_r[i];
		}
		Score_r[bestTemplates_r[i]] = 0;
	}
	*regionTemplates = bestHits;
	
	return bestScore;
}

int getR_Best(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates) {
	
	int i, j, score, bestScore_r, bestHits;
	
	/* get best scoring tempates */
	bestScore_r = 0;
	bestHits = 0;
	for(i = 1; i <= *bestTemplates; ++i) {
		if(bestScore_r < (score = Score[bestTemplates[i]])) {
			for(j = bestHits; j != 0; --j) {
				Score[bestTemplates[j]] = 0;
			}
			bestScore_r = score;
			bestHits = 1;
			bestTemplates[bestHits] = bestTemplates[i];
		} else if(bestScore_r == score) {
			++bestHits;
			bestTemplates[bestHits] = bestTemplates[i];
		} else {
			Score[bestTemplates[i]] = 0;
		}
	}
	for(i = 1; i <= *bestTemplates_r; ++i) {
		if(bestScore_r < (score = Score_r[bestTemplates_r[i]])) {
			for(j = bestHits; j != 0; --j) {
				if(0 < bestTemplates[j]) {
					Score[bestTemplates[j]] = 0;
				} else {
					Score_r[-bestTemplates[j]] = 0;
				}
			}
			bestScore_r = score;
			bestHits = 1;
			bestTemplates[bestHits] = -bestTemplates_r[i];
		} else if(bestScore_r == score) {
			++bestHits;
			bestTemplates[bestHits] = -bestTemplates_r[i];
		} else {
			Score_r[bestTemplates_r[i]] = 0;
		}
	}
	*bestTemplates = bestHits;
	
	/* check union */
	bestHits = 0;
	for(i = 1; i <= *regionTemplates; ++i) {
		if(0 < regionTemplates[i]) {
			/* we got one */
			if(Score_r[regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		} else {
			/* we got one */
			if(Score[-regionTemplates[i]]) {
				++bestHits;
				score = regionTemplates[bestHits];
				regionTemplates[bestHits] = regionTemplates[i];
				regionTemplates[i] = score;
			}
		}
	}
	
	/* mark as PE */
	if(bestHits) {
		*regionTemplates = -bestHits;
	}
	
	/* clear scores */
	for(i = *bestTemplates; i != 0; --i) {
		if(0 < bestTemplates[i]) {
			Score[bestTemplates[i]] = 0;
		} else {
			Score_r[-bestTemplates[i]] = 0;
		}
	}
	
	return bestScore_r;
}

void save_kmers_Sparse(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs* header, int *extendScore) {
	
	int i, j, k, l, n, end, rc, prefix_len, template, hitCounter, HIT, SU;
	int n_kmers, gaps, Ms, MMs, Us, W1s, score, bestScore, bestHits, reps;
	unsigned *values, *last;
	short unsigned *values_s;
	long unsigned prefix;
	
	if(qseq->seqlen < kmersize) {
		return;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	*bestTemplates = 0;
	prefix = templates->prefix;
	prefix_len = templates->prefix_len;
	hitCounter = 0;
	n_kmers = 0;
	bestScore = 0;
	end = qseq->seqlen;
	if(prefix_len) {
		for(rc = 0; rc < 2; ++rc) {
			if(rc) {
				comp_rc(qseq);
			}
			j = 0;
			qseq->N[0]++;
			qseq->N[qseq->N[0]] = qseq->seqlen;
			for(i = 1; i <= qseq->N[0]; ++i) {
				end = qseq->N[i] - kmersize - prefix_len + 1;
				while(j < end) {
					if(getSmer(qseq->seq, j) == prefix) {
						if((values = hashMap_get(templates, getKmer(qseq->seq, j + prefix_len)))) {
							if(SU) {
								values_s = (short unsigned *) values;
								n = *values_s;
								for(k = 1; k <= n; ++k) {
									if(Score[(template = values_s[k])] == 0) {
										(*bestTemplates)++;
										bestTemplates[*bestTemplates] = template;
									}
									Score[template]++;
								}
							} else {
								n = *values;
								for(k = 1; k <= n; ++k) {
									if(Score[(template = values[k])] == 0) {
										(*bestTemplates)++;
										bestTemplates[*bestTemplates] = template;
									}
									Score[template]++;
								}
							}
							++hitCounter;
						}
						++n_kmers;
					}
					++j;
				}
				j = qseq->N[i] + 1;
			}
			qseq->N[0]--;
		}
		
		/* get best match(es) */
		if(hitCounter * kmersize > (n_kmers - hitCounter)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				score = Score[(template = bestTemplates[l])];
				score = score * kmersize * M + (n_kmers - score) * MM;
				if(score > bestScore) {
					bestScore = score;
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(score == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				Score[bestTemplates[l]] = 0;
			}
			*bestTemplates = 0;
		}
		end = n_kmers - hitCounter - bestScore;
	} else {
		HIT = exhaustive;
		j = 0;
		qseq->N[0]++;
		qseq->N[qseq->N[0]] = qseq->seqlen;
		for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end && !HIT; j += kmersize) {
				if(hashMap_get(templates, getKmer(qseq->seq, j))) {
					HIT = 1;
				}
			}
			j = qseq->N[i] + 1;
		}
		
		if(HIT) {
			
			if(countK) {
				last = 0;
				reps = 0;
				j = 0;
				for(i = 1; i <= qseq->N[0]; ++i) {
					end = qseq->N[i] - kmersize + 1;
					for(;j < end; ++j) {
						if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
							if(values == last) {
								++reps;
							} else {
								if(last) {
									if(SU) {
										values_s = (short unsigned *) values;
										n = *values_s;
										for(l = 1; l <= n; ++l) {
											if(Score[(template = values_s[l])]) {
												Score[template] += reps;
											} else {
												Score[template] = reps;
												bestTemplates[0]++;
												bestTemplates[*bestTemplates] = template;
											}
										}
									} else {
										n = *values;
										for(l = 1; l <= n; ++l) {
											if(Score[(template = last[l])]) {
												Score[template] += reps;
											} else {
												Score[template] = reps;
												bestTemplates[0]++;
												bestTemplates[*bestTemplates] = template;
											}
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
					if(SU) {
						values_s = (short unsigned *) last;
						n = *values_s;
						for(l = 1; l <= n; ++l) {
							if(Score[(template = values_s[l])]) {
								Score[template] += reps;
							} else {
								Score[template] = reps;
								bestTemplates[0]++;
								bestTemplates[*bestTemplates] = template;
							}
						}
					} else {
						n = *last;
						for(l = 1; l <= n; ++l) {
							if(Score[(template = last[l])]) {
								Score[template] += reps;
							} else {
								Score[template] = reps;
								bestTemplates[0]++;
								bestTemplates[*bestTemplates] = template;
							}
						}
					}
					hitCounter += reps;
				}
			} else {
				last = 0;
				gaps = 0;
				HIT = 0;
				Ms = 0;
				MMs = 0;
				Us = 0;
				W1s = 0;
				j = 0;
				for(i = 1; i <= qseq->N[0]; ++i) {
					end = qseq->N[i] - kmersize + 1;
					for(;j < end; ++j) {
						if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
							if(values == last) {
								if(kmersize < gaps) {
									Ms += kmersize;
									gaps -= kmersize;
									if(gaps) {
										/* go for best scenario */
										if(gaps == 1) {
											MMs += 2;
										} else {
											gaps -= 2;
											if((MM << 1) + gaps * M < 0) {
												Ms += gaps;
												MMs += 2;
											}
										}
									} else {
										++MMs;
									}
								} else if (gaps) {
									--gaps;
									++W1s;
									Us += gaps;
								} else {
									++Ms;
								}
								HIT = j;
								gaps = 0;
							} else {
								if(last) {
									if(HIT) {
										HIT += kmersize;
									} else {
										HIT = j + kmersize;
									}
									score = Ms * M + MMs * MM + Us * U + W1s * W1;
									if(SU) {
										values_s = (short unsigned *) last;
										l = (*values_s) + 1;
										while(--l) {
											Score[(template = values_s[l])] += score;
											extendScore[template] = HIT;
										}
									} else {
										l = (*last) + 1;
										while(--l) {
											Score[(template = last[l])] += score;
											extendScore[template] = HIT;
										}
									}
									
									score = kmersize * M;
									MMs = MM << 1;
									if(SU) {
										values_s = (short unsigned *) values;
										n = *values_s;
										for(l = 1; l <= n; ++l) {
											if(j < extendScore[(template = values_s[l])]) {
												if(extendScore[template] == HIT) {
													Score[template] += M;
												} else {
													gaps = extendScore[template] - j - 1;
													Score[template] += (W1 + gaps * U);
												}
											} else if(Score[template] != 0) {
												Score[template] += score;
												if((gaps = extendScore[template] - j)) {
													if(gaps == 1) {
														Score[template] += MMs;
													} else {
														gaps -= 2;
														if((Ms = MMs + gaps * M) < 0) {
															Score[template] += Ms;
														}
													}
												} else {
													Score[template] += MM;
												}
											} else {
												Score[template] = score;
												bestTemplates[0]++;
												bestTemplates[*bestTemplates] = template;
											}
										}
									} else {
										n = *values;
										for(l = 1; l <= n; ++l) {
											if(j < extendScore[(template = values[l])]) {
												if(extendScore[template] == HIT) {
													Score[template] += M;
												} else {
													gaps = extendScore[template] - j - 1;
													Score[template] += (W1 + gaps * U);
												}
											} else if(Score[template] != 0) {
												Score[template] += score;
												if((gaps = extendScore[template] - j)) {
													if(gaps == 1) {
														Score[template] += MMs;
													} else {
														gaps -= 2;
														if((Ms = MMs + gaps * M) < 0) {
															Score[template] += Ms;
														}
													}
												} else {
													Score[template] += MM;
												}
											} else {
												Score[template] = score;
												bestTemplates[0]++;
												bestTemplates[*bestTemplates] = template;
											}
										}
									}
								} else if(SU) {
									values_s = (short unsigned *) values;
									n = *values_s;
									Ms = kmersize * M;
									for(l = 1; l <= n; ++l) {
										Score[(template = values_s[l])] = Ms;
										bestTemplates[l] = template;
									}
									*bestTemplates = n;
								} else {
									n = *values;
									Ms = kmersize * M;
									for(l = 1; l <= n; ++l) {
										Score[(template = values[l])] = Ms;
										bestTemplates[l] = template;
									}
									*bestTemplates = n;
								}
								HIT = 0;
								gaps = 0;
								Ms = 0;
								MMs = 0;
								Us = 0;
								W1s = 0;
								last = values;
							}
							++hitCounter;
						} else {
							++gaps;
						}
					}
					j = qseq->N[i] + 1;
				}
				if(last) {
					score = Ms * M + MMs * MM + Us * U + W1s * W1;
					if(SU) {
						values_s = (short unsigned *) last;
						l = (*values_s) + 1;
						while(--l) {
							Score[values_s[l]] += score;
						}
					} else {
						l = (*last) + 1;
						while(--l) {
							Score[last[l]] += score;
						}
					}
					for(l = *bestTemplates; l != 0; --l) {
						extendScore[(template = bestTemplates[l])] = 0;
						if(Score[template] < 0) {
							Score[template] = 0;
						}
					}
				}
			}
		}
		qseq->N[0]--;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				Score[template] = 0;
			}
			*bestTemplates = 0;
		}
		end = qseq->seqlen + 1 - bestScore;
	}
	
	if(bestScore) {
		if(bestScore * kmersize > end) {
			lock(excludeOut);
			deConPrintPtr(bestTemplates, qseq, bestScore, header);
			unlock(excludeOut);
		}
	}
}

void save_kmers(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs* header, int *extendScore) {
	
	int i, j, l, end, HIT, gaps, score, Ms, MMs, Us, W1s, template;
	int bestHits, hitCounter, bestScore, bestScore_r;
	unsigned *values, *last, n, SU;
	short unsigned *values_s;
	
	if(qseq->seqlen < kmersize) {
		return;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
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
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq->seq, j))) {
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
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		end = qseq->seqlen;
		for(i = 1; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
					if(values == last) {
						if(kmersize < gaps) {
							Ms += kmersize;
							gaps -= kmersize;
							if(gaps) {
								/* go for best scenario */
								if(gaps == 1) {
									MMs += 2;
								} else {
									gaps -= 2;
									if((MM << 1) + gaps * M < 0) {
										Ms += gaps;
										MMs += 2;
									}
								}
							} else {
								++MMs;
							}
						} else if (gaps) {
							--gaps;
							++W1s;
							Us += gaps;
						} else {
							++Ms;
						}
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							if(HIT) {
								HIT += kmersize;
							} else {
								HIT = j + kmersize;
							}
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							
							if(SU) {
								values_s = (short unsigned *) last;
								l = (*values_s) + 1;
								while(--l) {
									Score[(template = values_s[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								values_s = (short unsigned *) values;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values_s[l])]) {
										if(extendScore[template] == HIT) {
											Score[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score[template] += (W1 + gaps * U);
										}
									} else if(Score[template] != 0) {
										Score[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score[template] += Ms;
												}
											}
										} else {
											Score[template] += MM;
										}
									} else {
										Score[template] = score;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							} else {
								l = (*last) + 1;
								while(--l) {
									Score[(template = last[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								n = *values;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values[l])]) {
										if(extendScore[template] == HIT) {
											Score[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score[template] += (W1 + gaps * U);
										}
									} else if(Score[template] != 0) {
										Score[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score[template] += Ms;
												}
											}
										} else {
											Score[template] += MM;
										}
									} else {
										Score[template] = score;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							}
						} else if(SU) {
							values_s = (short unsigned *) values;
							n = *values_s;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score[(template = values_s[l])] = Ms;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						} else {
							n = *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score[(template = values[l])] = Ms;
								bestTemplates[l] = template;
							}
							*bestTemplates = n;
						}
						HIT = 0;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
						last = values;
					}
					++hitCounter;
				} else {
					++gaps;
				}
			}
			j = qseq->N[i] + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score[last[l]] += score;
				}
			}
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				if(Score[template] < 0) {
					Score[template] = 0;
				}
			}
		}
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
				extendScore[(template = bestTemplates[l])] = 0;
				Score[template] = 0;
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
	for(i = 1; i <= qseq_r->N[0] && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq_r->seq, j))) {
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
		gaps = 0;
		HIT = 0;
		Ms = 0;
		MMs = 0;
		Us = 0;
		W1s = 0;
		j = 0;
		end = qseq_r->seqlen;
		for(i = 1; i <= qseq_r->N[0]; ++i) {
			end = qseq_r->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq_r->seq, j)))) {
					if(values == last) {
						if(kmersize < gaps) {
							Ms += kmersize;
							gaps -= kmersize;
							if(gaps) {
								/* go for best scenario */
								if(gaps == 1) {
									MMs += 2;
								} else {
									gaps -= 2;
									if((MM << 1) + gaps * M < 0) {
										Ms += gaps;
										MMs += 2;
									}
								}
							} else {
								++MMs;
							}
						} else if (gaps) {
							--gaps;
							++W1s;
							Us += gaps;
						} else {
							++Ms;
						}
						HIT = j;
						gaps = 0;
					} else {
						if(last) {
							if(HIT) {
								HIT += kmersize;
							} else {
								HIT = j + kmersize;
							}
							score = Ms * M + MMs * MM + Us * U + W1s * W1;
							if(SU) {
								values_s = (short unsigned *) last;
								l = (*values_s) + 1;
								while(--l) {
									Score_r[(template = values_s[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								values_s = (short unsigned *) values;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values_s[l])]) {
										if(extendScore[template] == HIT) {
											Score_r[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score_r[template] += (W1 + gaps * U);
										}
									} else if(Score_r[template] != 0) {
										Score_r[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score_r[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score_r[template] += Ms;
												}
											}
										} else {
											Score_r[template] += MM;
										}
									} else {
										Score_r[template] = score;
										bestTemplates_r[0]++;
										bestTemplates_r[*bestTemplates_r] = template;
									}
								}
							} else {
								l = (*last) + 1;
								while(--l) {
									Score_r[(template = last[l])] += score;
									extendScore[template] = HIT;
								}
								
								score = kmersize * M;
								MMs = MM << 1;
								n = *values;
								for(l = 1; l <= n; ++l) {
									if(j < extendScore[(template = values[l])]) {
										if(extendScore[template] == HIT) {
											Score_r[template] += M;
										} else {
											gaps = extendScore[template] - j - 1;
											Score_r[template] += (W1 + gaps * U);
										}
									} else if(Score_r[template] != 0) {
										Score_r[template] += score;
										if((gaps = extendScore[template] - j)) {
											if(gaps == 1) {
												Score_r[template] += MMs;
											} else {
												gaps -= 2;
												if((Ms = MMs + gaps * M) < 0) {
													Score_r[template] += Ms;
												}
											}
										} else {
											Score_r[template] += MM;
										}
									} else {
										Score_r[template] = score;
										bestTemplates_r[0]++;
										bestTemplates_r[*bestTemplates_r] = template;
									}
								}
							}
						} else if(SU) {
							values_s = (short unsigned *) values;
							n = *values_s;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score_r[(template = values_s[l])] = Ms;
								bestTemplates_r[l] = template;
							}
							*bestTemplates_r = n;
						} else {
							n = *values;
							Ms = kmersize * M;
							for(l = 1; l <= n; ++l) {
								Score_r[(template = values[l])] = Ms;
								bestTemplates_r[l] = template;
							}
							*bestTemplates_r = n;
						}
						HIT = 0;
						gaps = 0;
						Ms = 0;
						MMs = 0;
						Us = 0;
						W1s = 0;
						last = values;
					}
					++hitCounter;
				} else {
					++gaps;
				}
			}
			j = qseq_r->N[i] + 1;
		}
		if(last) {
			score = Ms * M + MMs * MM + Us * U + W1s * W1;
			if(SU) {
				values_s = (short unsigned *) last;
				l = (*values_s) + 1;
				while(--l) {
					Score_r[values_s[l]] += score;
				}
			} else {
				l = (*last) + 1;
				while(--l) {
					Score_r[last[l]] += score;
				}
			}
		}
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; ++l) {
				if(Score_r[(template = bestTemplates_r[l])] > bestScore_r) {
					bestScore_r = Score_r[template];
					bestHits = 1;
					bestTemplates_r[bestHits] = template;
				} else if(Score_r[template] == bestScore_r) {
					++bestHits;
					bestTemplates_r[bestHits] = template;
				}
				Score_r[template] = 0;
				extendScore[template] = 0;
			}
			*bestTemplates_r = bestHits;
		} else {
			for(l = *bestTemplates_r; l != 0; --l) {
				extendScore[(template = bestTemplates_r[l])] = 0;
				Score_r[template] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	if(bestScore > 0 || bestScore_r > 0) {
		end = qseq->seqlen + 1;
		if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r))) {
			if(bestScore > bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, bestScore, header);
				unlock(excludeOut);
			} else if(bestScore < bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header);
				unlock(excludeOut);
			} else {
				/* merge */
				for(i = 1; i <= *bestTemplates_r; ++i) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, -bestScore, header);
				unlock(excludeOut);
			}
		}
	}
}

void save_kmers_count(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs* header, int *extendScore) {
	
	int i, j, l, end, HIT, bestHits, hitCounter, bestScore, bestScore_r, reps;
	int n, template, SU;
	unsigned *values, *last;
	short unsigned *values_s;
	
	if(qseq->seqlen < kmersize) {
		return;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
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
	for(i = 1; i <= qseq->N[0] && !HIT; ++i) {
		end = qseq->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq->seq, j))) {
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
		end = qseq->seqlen;
		for(i = 1; i <= qseq->N[0]; ++i) {
			end = qseq->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq->seq, j)))) {
					if(values == last) {
						++reps;
					} else {
						if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(Score[(template = values_s[l])]) {
										Score[template] += reps;
									} else {
										Score[template] = reps;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
								}
							} else {
								n = *last;
								for(l = 1; l <= n; ++l) {
									if(Score[(template = last[l])]) {
										Score[template] += reps;
									} else {
										Score[template] = reps;
										bestTemplates[0]++;
										bestTemplates[*bestTemplates] = template;
									}
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
			if(SU) {
				values_s = (short unsigned *) last;
				n = *values_s;
				for(l = 1; l <= n; ++l) {
					if(Score[(template = values_s[l])]) {
						Score[template] += reps;
					} else {
						Score[template] = reps;
						bestTemplates[0]++;
						bestTemplates[*bestTemplates] = template;
					}
				}
			} else {
				n = *last;
				for(l = 1; l <= n; ++l) {
					if(Score[(template = last[l])]) {
						Score[template] += reps;
					} else {
						Score[template] = reps;
						bestTemplates[0]++;
						bestTemplates[*bestTemplates] = template;
					}
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates; ++l) {
				if(Score[(template = bestTemplates[l])] > bestScore) {
					bestScore = Score[template];
					bestHits = 1;
					bestTemplates[bestHits] = template;
				} else if(Score[template] == bestScore) {
					++bestHits;
					bestTemplates[bestHits] = template;
				}
				Score[template] = 0;
			}
			*bestTemplates = bestHits;
		} else {
			for(l = *bestTemplates; l != 0; --l) {
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
	for(i = 1; i <= qseq_r->N[0] && !HIT; ++i) {
		end = qseq_r->N[i] - kmersize + 1;
		for(;j < end && !HIT; j += kmersize) {
			if(hashMap_get(templates, getKmer(qseq_r->seq, j))) {
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
		end = qseq_r->seqlen;
		for(i = 1; i <= qseq_r->N[0]; ++i) {
			end = qseq_r->N[i] - kmersize + 1;
			for(;j < end; ++j) {
				if((values = hashMap_get(templates, getKmer(qseq_r->seq, j)))) {
					if(values == last) {
						++reps;
					} else {
						if(last) {
							if(SU) {
								values_s = (short unsigned *) last;
								n = *values_s;
								for(l = 1; l <= n; ++l) {
									if(Score_r[(template = values_s[l])]) {
										Score_r[template] += reps;
									} else {
										Score_r[template] = reps;
										bestTemplates_r[0]++;
										bestTemplates_r[*bestTemplates_r] = template;
									}	
								}
							} else {
								n = *last;
								for(l = 1; l <= n; ++l) {
									if(Score_r[(template = last[l])]) {
										Score_r[template] += reps;
									} else {
										Score_r[template] = reps;
										bestTemplates_r[0]++;
										bestTemplates_r[*bestTemplates_r] = template;
									}	
								}
							}
							hitCounter += reps;
						}
						last = values;
						reps = 1;
					}
				}
			}
			j = qseq_r->N[i] + 1;
		}
		if(last) {
			if(SU) {
				values_s = (short unsigned *) last;
				n = *values_s;
				for(l = 1; l <= n; ++l) {
					if(Score_r[(template = values_s[l])]) {
						Score_r[template] += reps;
					} else {
						Score_r[template] = reps;
						bestTemplates_r[0]++;
						bestTemplates_r[*bestTemplates_r] = template;
					}	
				}
			} else {
				n = *last;
				for(l = 1; l <= n; ++l) {
					if(Score_r[(template = last[l])]) {
						Score_r[template] += reps;
					} else {
						Score_r[template] = reps;
						bestTemplates_r[0]++;
						bestTemplates_r[*bestTemplates_r] = template;
					}	
				}
			}
			hitCounter += reps;
		}
		reps = 0;
		
		/* get best match(es) */
		if(hitCounter * kmersize > (end - hitCounter + kmersize)) {
			bestHits = 0;
			for(l = 1; l <= *bestTemplates_r; ++l) {
				if(Score_r[(template = bestTemplates_r[l])] > bestScore_r) {
					bestScore_r = Score_r[template];
					bestHits = 1;
					bestTemplates_r[bestHits] = template;
				} else if(Score_r[template] == bestScore_r) {
					++bestHits;
					bestTemplates_r[bestHits] = template;
				}
				Score_r[template] = 0;
			}
			*bestTemplates_r = bestHits;
		} else {
			for(l = *bestTemplates_r; l != 0; --l) {
				Score_r[bestTemplates_r[l]] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	qseq_r->N[0]--;
	
	/* Validate best match */
	if(bestScore > 0 || bestScore_r > 0) {
		end = qseq->seqlen + 1;
		if((bestScore >= bestScore_r && bestScore * kmersize > (end - bestScore)) || (bestScore < bestScore_r && bestScore_r * kmersize > (end - bestScore_r))) {
			if(bestScore > bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, bestScore, header);
				unlock(excludeOut);
			} else if(bestScore < bestScore_r) {
				lock(excludeOut);
				deConPrintPtr(bestTemplates_r, qseq_r, bestScore_r, header);
				unlock(excludeOut);
			} else {
				/* merge */
				for(i = 1; i <= *bestTemplates_r; ++i) {
					bestTemplates[0]++;
					bestTemplates[*bestTemplates] = -bestTemplates_r[i];
				}
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq, -bestScore, header);
				unlock(excludeOut);
			}
		}
	}
}

void save_kmers_unionPair(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs* header, struct qseqs* header_r, int *extendScore) {
	
	int i, bestScore, bestScore_r, hitCounter;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore)) && (qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
		/* got hits */
		bestScore = getF_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		
		if(bestScore * kmersize < (qseq->seqlen - bestScore)) {
			bestScore = 0;
		}
	} else {
		bestScore = 0;
		if(hitCounter) {
			i = *bestTemplates + 1;
			while(--i) {
				Score[bestTemplates[i]] = 0;
			}
			*bestTemplates = 0;
			i = *bestTemplates_r + 1;
			while(--i) {
				Score_r[bestTemplates_r[i]] = 0;
			}
			*bestTemplates_r = 0;
		}
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter = get_kmers_for_pair_ptr(bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore)) && (qseq_r->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
		if(bestScore) {
			bestScore_r = getR_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		} else {
			bestScore_r = getF_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		}
		if(bestScore_r * kmersize < (qseq_r->seqlen - bestScore_r)) {
			bestScore_r = 0;
			*regionTemplates = abs(*regionTemplates);
		}
	} else {
		bestScore_r = 0;
		i = *bestTemplates + 1;
		while(--i) {
			Score[bestTemplates[i]] = 0;
		}
		*bestTemplates = 0;
		i = *bestTemplates_r + 1;
		while(--i) {
			Score_r[bestTemplates_r[i]] = 0;
		}
		*bestTemplates_r = 0;
	}
	
	if(0 < bestScore && 0 < bestScore_r) {
		if(*regionTemplates < 0) {
			*regionTemplates = -(*regionTemplates);
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
					bestScore_r = -bestScore_r;
				}
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore_r, header_r);
				unlock(excludeOut);
			} else {
				comp_rc(qseq_r);
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq_r, bestScore_r, header_r, qseq, bestScore, header);
				unlock(excludeOut);
			}
		} else {
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
				}
			} else {
				for(i = 1; i <= *regionTemplates; ++i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(regionTemplates, qseq, bestScore, header);
			unlock(excludeOut);
			if(0 < bestTemplates[1]) {
				comp_rc(qseq_r);
				if(bestTemplates[*bestTemplates] < 0) {
					bestScore_r = -bestScore_r;
				}
			} else {
				for(i = 1; i <= *bestTemplates; ++i) {
					bestTemplates[i] = -bestTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(bestTemplates, qseq_r, bestScore_r, header_r);
			unlock(excludeOut);
		}
	} else if(bestScore) {
		if(0 < regionTemplates[1]) {
			comp_rc(qseq);
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore = -bestScore;
			}
		} else {
			for(i = 1; i <= *regionTemplates; ++i) {
				regionTemplates[i] = -regionTemplates[i];
			}
		}
		lock(excludeOut);
		deConPrintPtr(regionTemplates, qseq, bestScore, header);
		unlock(excludeOut);
	} else if(bestScore_r) {
		if(0 < regionTemplates[1]) {
			comp_rc(qseq_r);
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore_r = -bestScore_r;
			}
		} else {
			for(i = 1; i <= *regionTemplates; ++i) {
				regionTemplates[i] = -regionTemplates[i];
			}
		}
		lock(excludeOut);
		deConPrintPtr(regionTemplates, qseq_r, bestScore_r, header_r);
		unlock(excludeOut);
	}
}

void save_kmers_penaltyPair(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs* header, struct qseqs* header_r, int *extendScore) {
	
	int i, bestScore, bestScore_r, compScore, hitCounter, hitCounter_r;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore))) {
		/* got hits */
		bestScore = getFirstPen(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores);
	} else {
		bestScore = 0;
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter_r = get_kmers_for_pair_ptr(bestTemplates, bestTemplates_r, Score, Score_r, qseq_r, extendScore))) {
		if(0 < bestScore) {
			bestScore_r = getSecondPen(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores, bestScore);
		} else {
			bestScore_r = getF_Best(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates);
		}
	} else {
		bestScore_r = 0;
	}
	
	if(0 < bestScore && 0 < bestScore_r) {
		if(*regionTemplates < 0) {
			compScore = MIN((hitCounter + hitCounter_r), (bestScore + bestScore_r));
			if((qseq->seqlen + qseq_r->seqlen - compScore - (kmersize << 1)) < compScore * kmersize) {
				*regionTemplates = -(*regionTemplates);
				if(0 < regionTemplates[1]) {
					comp_rc(qseq);
					if(regionTemplates[*regionTemplates] < 0) {
						bestScore = -bestScore;
						bestScore_r = -bestScore_r;
					}
					lock(excludeOut);
					printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore_r, header_r);
					unlock(excludeOut);
				} else {
					comp_rc(qseq_r);
					for(i = *regionTemplates; i != 0; --i) {
						regionTemplates[i] = -regionTemplates[i];
					}
					lock(excludeOut);
					printPairPtr(regionTemplates, qseq_r, bestScore_r, header_r, qseq, bestScore, header);
					unlock(excludeOut);
				}
			}
		} else {
			hitCounter = MIN(hitCounter, bestScore);
			if((qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
				if(0 < regionTemplates[1]) {
					comp_rc(qseq);
					if(regionTemplates[*regionTemplates] < 0) {
						bestScore = -bestScore;
					}
				} else {
					for(i = *regionTemplates; i != 0; --i) {
						regionTemplates[i] = -regionTemplates[i];
					}
				}
				lock(excludeOut);
				deConPrintPtr(regionTemplates, qseq, bestScore, header);
				unlock(excludeOut);
			}
			hitCounter_r = MIN(hitCounter_r, bestScore_r);
			if((qseq_r->seqlen - hitCounter_r - kmersize) < hitCounter_r * kmersize) {
				if(0 < bestTemplates[1]) {
					comp_rc(qseq_r);
					if(bestTemplates[*bestTemplates] < 0) {
						bestScore_r = -bestScore_r;
					}
				} else {
					for(i = *bestTemplates; i != 0; --i) {
						bestTemplates[i] = -bestTemplates[i];
					}
				}
				lock(excludeOut);
				deConPrintPtr(bestTemplates, qseq_r, bestScore_r, header_r);
				unlock(excludeOut);
			}
		}
	} else if(0 < bestScore) {
		hitCounter = MIN(hitCounter, bestScore);
		if((qseq->seqlen - hitCounter - kmersize) < hitCounter * kmersize) {
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore = -bestScore;
				}
			} else {
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(regionTemplates, qseq, bestScore, header);
			unlock(excludeOut);
		}
	} else if(0 < bestScore_r) {
		hitCounter_r = MIN(hitCounter_r, bestScore_r);
		if((qseq_r->seqlen - hitCounter_r - kmersize) < hitCounter_r * kmersize) {
			if(0 < regionTemplates[1]) {
				comp_rc(qseq_r);
				if(regionTemplates[*regionTemplates] < 0) {
					bestScore_r = -bestScore_r;
				}
			} else {
				for(i = 1; i <= *regionTemplates; ++i) {
					regionTemplates[i] = -regionTemplates[i];
				}
			}
			lock(excludeOut);
			deConPrintPtr(regionTemplates, qseq_r, bestScore_r, header_r);
			unlock(excludeOut);
		}
	}
}

void save_kmers_forcePair(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, int *regionTemplates, int *regionScores, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs* header, struct qseqs* header_r, int *extendScore) {
	
	int i, bestScore, hitCounter, hitCounter_r;
	
	/* get_kmers_for_pair, returns a positive number if templates are found.
	zero otherwise */
	
	/* get forward */
	if((hitCounter = get_kmers_for_pair_ptr(bestTemplates, bestTemplates_r, Score, Score_r, qseq, extendScore))) {
		/* got hits */
		getFirstForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores);
	} else {
		return;
	}
	*extendScore = 1;
	
	/* get reverse */
	if((hitCounter_r = get_kmers_for_pair_ptr(bestTemplates_r, bestTemplates, Score_r, Score, qseq_r, extendScore)) && 
		(qseq->seqlen + qseq_r->seqlen - hitCounter - hitCounter_r - (kmersize << 1)) < (hitCounter + hitCounter_r) * kmersize && 
	(bestScore = getSecondForce(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores))) {
		
		if((qseq->seqlen + qseq_r->seqlen - bestScore) < bestScore * kmersize) {
			if(regionTemplates[*regionTemplates] < 0) {
				bestScore = -bestScore;
			}
			if(0 < regionTemplates[1]) {
				comp_rc(qseq);
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq, bestScore, header, qseq_r, bestScore, header_r);
				unlock(excludeOut);
			} else {
				comp_rc(qseq_r);
				for(i = *regionTemplates; i != 0; --i) {
					regionTemplates[i] = -regionTemplates[i];
				}
				lock(excludeOut);
				printPairPtr(regionTemplates, qseq_r, bestScore, header_r, qseq, bestScore, header);
				unlock(excludeOut);
			}
		}
	} else if(hitCounter || hitCounter_r) {
		i = *bestTemplates + 1;
		while(--i) {
			Score[bestTemplates[i]] = 0;
		}
		*bestTemplates = 0;
		i = *bestTemplates_r + 1;
		while(--i) {
			Score_r[bestTemplates_r[i]] = 0;
		}
		*bestTemplates_r = 0;
	}
}

void * save_kmers_threaded(void *arg) {
	
	static unsigned readNum = 0;
	struct kmerScan_thread *thread = arg;
	int *Score, *Score_r, *bestTemplates, *bestTemplates_r, *regionTemplates;
	int *regionScores, *extendScore, go, spltDB;
	FILE *inputfile;
	struct compDNA *qseq, *qseq_r;
	struct qseqs *header, *header_r;
	
	if(save_kmers_pair != &save_kmers_unionPair) {
		regionScores = calloc(DB_size << 1, sizeof(int));
		if(!regionScores) {
			ERROR();
		}
	} else {
		regionScores = 0;
	}
	extendScore = calloc((DB_size + 1), sizeof(int));
	header_r = malloc(sizeof(struct qseqs));
	if(!extendScore || !header_r) {
		ERROR();
	}
	header_r = setQseqs(256);
	
	qseq = thread->qseq;
	qseq_r = thread->qseq_r;
	header = thread->header;
	bestTemplates = thread->bestTemplates;
	bestTemplates_r = thread->bestTemplates_r;
	regionTemplates = RegionTemplates[thread->num];
	Score = tScore[thread->num];
	Score_r = tScore_r[thread->num];
	inputfile = thread->inputfile;
	*Score = thread->num;
	*bestTemplates = thread->num;
	*bestTemplates_r = thread->num;
	bestTemplates += 2;
	bestTemplates_r += 2;
	if(printPtr == &print_ankers_spltDB) {
		spltDB = 1;
	} else {
		spltDB = 0;
	}
	
	go = 1;
	while(go != 0) {
		
		/* load qseqs */
		lock(excludeIn);
		if((go = loadFsa(qseq, header, inputfile)) < 0) {
			/* PE */
			loadFsa(qseq_r, header_r, inputfile);
		}
		if(spltDB) {
			bestTemplates[-1] = ++readNum;
			bestTemplates_r[-1] = readNum;
		}
		unlock(excludeIn);
		
		/* allocate memory */
		if(qseq_r->size < qseq->size && 0 < go) {
			freeComp(qseq_r);
			allocComp(qseq_r, qseq->size);
			if(!one2one) {
				free(TmpNs[thread->num]);
				free(tVF_scores[thread->num]);
				free(tVR_scores[thread->num]);
				tVF_scores[thread->num] = calloc(qseq->size, sizeof(int*));
				tVR_scores[thread->num] = calloc(qseq->size, sizeof(int*));
				TmpNs[thread->num] = malloc(qseq->size * sizeof(int));
				if(!tVF_scores[thread->num] || !tVR_scores[thread->num] || !TmpNs[thread->num]) {
					ERROR();
				}
			}
		}
		
		/* find ankers */
		if(0 < go) {
			kmerScan(bestTemplates, bestTemplates_r, Score, Score_r, qseq, qseq_r, header, extendScore);
		} else if(go < 0) {
			save_kmers_pair(bestTemplates, bestTemplates_r, Score, Score_r, regionTemplates, regionScores, qseq, qseq_r, header, header_r, extendScore);
		}
	}
	
	return NULL;
}

void save_kmers_HMM(int *bestTemplates, int *bestTemplates_r, int *Score, int *Score_r, struct compDNA *qseq, struct compDNA *qseq_r, struct qseqs *header, int *extendScore) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	   and is the time determining step */
	int i, j, k, l, N, n, i_r, j_r, seqlen, seqend, end, HIT, Ncheck, SU;
	int hitCounter, template, bestScore, bestHits, start, stop;
	int start_cut, end_cut, *regionTemplates;
	unsigned *values, *last, *rlast, **VF_scores, **VR_scores, num;
	short unsigned *values_s;
	int *tmpNs, reps, rreps;
	double Ms, Ns, Ms_prev, Ns_prev;
	
	if(qseq->seqlen < kmersize) {
		return;
	} else if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	regionTemplates = RegionTemplates[*Score];
	VF_scores = tVF_scores[*Score];
	VR_scores = tVR_scores[*Score];
	tmpNs = TmpNs[*Score];
	
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
				if(hashMap_get(templates, getKmer(qseq->seq, i)) || hashMap_get(templates, getKmer(qseq_r->seq, i_r))) {
					HIT = 1;
				} else {
					++i;
					--i_r;
				}
			}
		} else {
			while(i < end && !HIT) {
				if(hashMap_get(templates, getKmer(qseq->seq, i)) || hashMap_get(templates, getKmer(qseq_r->seq, i_r))) {
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
			VF_scores[i] = hashMap_get(templates, getKmer(qseq->seq, i));
			VR_scores[i] = hashMap_get(templates, getKmer(qseq_r->seq, i_r));
			
			/* init HMM */
			Ms_prev = HMM_param[7] + HMM_param[2];
			Ns_prev = HMM_param[5] + HMM_param[0];
			Ms = 0;
			Ns = 0;
			
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
							--n;
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
						--k;
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
					VF_scores[j] = hashMap_get(templates, getKmer(qseq->seq, j));
					VR_scores[j] = hashMap_get(templates, getKmer(qseq_r->seq, j_r));
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						++hitCounter;
						
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
							--j;
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
							--j;
							break;
						}
					}
				}
				--j;
				++j_r;
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
							++N;
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
						++k;
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
					VF_scores[j] = hashMap_get(templates, getKmer(qseq->seq, j));
					VR_scores[j] = hashMap_get(templates, getKmer(qseq_r->seq, j_r));
					
					/* HMM */
					if(VF_scores[j] || VR_scores[j]) {
						++hitCounter;
						
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
							++j;
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
							++j;
							break;
						}
					}
				}
				++j;
				--j_r;
				Ns_prev = Ns;
				Ms_prev = Ms;
			}
			stop = j + kmersize - 1;
			
			/* evaluate hit */
			if(hitCounter > 0 && (hitCounter * kmersize > (stop - start - hitCounter + kmersize)) 
				&& ((stop - start) > minLen || start == 0 || stop == seqlen)) {
				if(deCon) {
					if(SU) {
						for(k = start; k < j; ++k) {
							if(((values_s = (short unsigned *) VF_scores[k]) && values_s[*values_s] == contamination) 
								|| ((values_s = (short unsigned *) VR_scores[k]) && values_s[*values_s] == contamination)) {
								--hitCounter;
							}
						}
					} else {
						for(k = start; k < j; ++k) {
							if(((values = VF_scores[k]) && values[*values] == contamination)
								|| ((values = VR_scores[k]) && values[*values] == contamination)) {
								--hitCounter;
							}
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
					for(k = start; k < j; ++k) {
						/* forward */
						if(VF_scores[k]) {
							if(VF_scores[k] == last) {
								++reps;
							} else {
								if(last) {
									if(SU) {
										values_s = (short unsigned *) last;
										num = *values_s;
										for(l = 1; l <= num; ++l) {
											if(Score[(template = values_s[l])]) {
												Score[template] += reps;
											} else {
												Score[template] = reps;
												bestTemplates[0]++;
												bestTemplates[*bestTemplates] = template;
											}
										}
									} else {
										num = *last;
										for(l = 1; l <= num; ++l) {
											if(Score[(template = last[l])]) {
												Score[template] += reps;
											} else {
												Score[template] = reps;
												bestTemplates[0]++;
												bestTemplates[*bestTemplates] = template;
											}
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
								++rreps;
							} else {
								if(rlast) {
									if(SU) {
										values_s = (short unsigned *) rlast;
										num = *values_s;
										for(l = 1; l <= num; ++l) {
											if(Score_r[(template = values_s[l])] > 0) {
												Score_r[template] += rreps;
											} else {
												Score_r[template] = rreps;
												bestTemplates_r[0]++;
												bestTemplates_r[*bestTemplates_r] = template;
											}
										}
									} else {
										num = *rlast;
										for(l = 1; l <= num; ++l) {
											if(Score_r[(template = rlast[l])] > 0) {
												Score_r[template] += rreps;
											} else {
												Score_r[template] = rreps;
												bestTemplates_r[0]++;
												bestTemplates_r[*bestTemplates_r] = template;
											}
										}
									}
								}
								rreps = 1;
								rlast = VR_scores[k];
							}
						}
						
					}
					if(last) {
						if(SU) {
							values_s = (short unsigned *) last;
							num = *values_s;
							for(l = 1; l <= num; ++l) {
								if(Score[(template = values_s[l])] > 0) {
									Score[template] += reps;
								} else {
									Score[template] = reps;
									bestTemplates[0]++;
									bestTemplates[*bestTemplates] = template;
								}
							}
						} else {
							num = *last;
							for(l = 1; l <= num; ++l) {
								if(Score[(template = last[l])] > 0) {
									Score[template] += reps;
								} else {
									Score[template] = reps;
									bestTemplates[0]++;
									bestTemplates[*bestTemplates] = template;
								}
							}
						}
					}
					if(rlast) {
						if(SU) {
							values_s = (short unsigned *) rlast;
							num = *values_s;
							for(l = 1; l <= num; ++l) {
								if(Score_r[(template = values_s[l])] > 0) {
									Score_r[template] += rreps;
								} else {
									Score_r[template] = rreps;
									bestTemplates_r[0]++;
									bestTemplates_r[*bestTemplates_r] = template;
								}
							}
						} else {
							num = *rlast;
							for(l = 1; l <= num; ++l) {
								if(Score_r[(template = rlast[l])] > 0) {
									Score_r[template] += rreps;
								} else {
									Score_r[template] = rreps;
									bestTemplates_r[0]++;
									bestTemplates_r[*bestTemplates_r] = template;
								}
							}
						}
					}
					
					/* cut out template hits */
					while(HIT != 0) {
						/* get best score */
						bestScore = 0;
						bestHits = 0;
						
						/* forward */
						for(k = 1; k <= *bestTemplates; ++k) {
							template = bestTemplates[k];
							if(Score[template] > bestScore) {
								bestScore = Score[template];
								bestHits = 1;
								regionTemplates[bestHits] = template;
							} else if(Score[template] == bestScore) {
								if(Score[template]) {
									++bestHits;
									regionTemplates[bestHits] = template;
								} else {
									bestTemplates[k] = bestTemplates[*bestTemplates];
									bestTemplates[0]--;
									--k;
								}
							}
						}
						
						/* rc */
						for(k = 1; k <= *bestTemplates_r; ++k) {
							template = bestTemplates_r[k];
							if(Score_r[template] > bestScore) {
								bestScore = Score_r[template];
								bestHits = 1;
								regionTemplates[bestHits] = -template;
							} else if(Score_r[template] == bestScore) {
								if(Score_r[template]) {
									++bestHits;
									regionTemplates[bestHits] = -template;
								} else {
									bestTemplates_r[k] = bestTemplates_r[*bestTemplates_r];
									bestTemplates_r[0]--;
									--k;
								}
							}
						}
						
						*regionTemplates = bestHits;
						if(bestScore > 0) {
							/* find limits of match */
							start_cut = j;
							for(k = 1; k <= bestHits; ++k) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = start; l < start_cut; ++l) {
									if(VR_scores[l] && intpos_bin_contaminationPtr(VR_scores[l], template) != -1) {
										start_cut = l;
									}
									if(VF_scores[l] && intpos_bin_contaminationPtr(VF_scores[l], template) != -1) {
										start_cut = l;
									}
								}
							}
							end_cut = start_cut;
							for(k = 1; k <= bestHits; ++k) {
								template = (regionTemplates[k] > 0) ? regionTemplates[k] : -regionTemplates[k];
								for(l = j; l > end_cut; --l) {
									if(VR_scores[l] && intpos_bin_contaminationPtr(VR_scores[l], template) != -1) {
										end_cut = l;
									}
									if(VF_scores[l] && intpos_bin_contaminationPtr(VF_scores[l], template) != -1) {
										end_cut = l;
									}
								}
							}
							
							/* evaluate best hit */
							if(bestScore * kmersize > (end_cut - start_cut - bestScore + kmersize)) {
								/* check for hits on rc */
								HIT = (regionTemplates[*regionTemplates] > 0) ? 1 : -1;
								/* print */
								if(start != 0 && j != qseq->seqlen) {
									ankerAndClean(regionTemplates, Score, Score_r, VF_scores, VR_scores, tmpNs, qseq, HIT, bestScore, start_cut, end_cut, header);
								} else {
									ankerPtr(regionTemplates, Score, Score_r, VF_scores, VR_scores, tmpNs, qseq, HIT, bestScore, start_cut, end_cut, header);
								}
							} else {
								/* clear scores */
								for(k = 1; k <= *bestTemplates; ++k) {
									Score[bestTemplates[k]] = 0;
								}
								for(k = 1; k <= *bestTemplates_r; ++k) {
									Score_r[bestTemplates_r[k]] = 0;
								}
								HIT = 0;
							}
						} else {
							/* clear scores */
							for(k = 1; k <= *bestTemplates; ++k) {
								Score[bestTemplates[k]] = 0;
							}
							for(k = 1; k <= *bestTemplates_r; ++k) {
								Score_r[bestTemplates_r[k]] = 0;
							}
							
							HIT = 0;
						}
					}
				}
			}
			
			/* clear scores */
			for(k = start; k < j; ++k) {
				VF_scores[k] = 0;
				VR_scores[k] = 0;
			}
			
			i = stop + 1;
			i_r = seqlen - kmersize - i;
		} else {
			++N;
		}
	}
}

void save_kmers_sparse(struct hashMap_kmers *foundKmers, struct compKmers *compressor) {
	
	/* save_kmers find ankering k-mers the in query sequence,
	and is the time determining step */
	
	int i;
	
	for(i = 0; i < compressor->n; ++i) {
		if(hashMap_get(templates, compressor->kmers[i])) {
			hashMap_kmers_CountIndex(foundKmers, compressor->kmers[i]);
		}
	}
	
}

struct hashTable * collect_Kmers(int *Scores, int *Scores_tot, struct hashMap_kmers *foundKmers, struct Hit *hits) {
	
	int template, SU;
	unsigned i, j, *value;
	short unsigned *values_s;
	struct hashTable_kmers *node, *node_next;
	struct hashTable *kmerNode, *kmerList;
	
	if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	for(i = 0; i < foundKmers->size; ++i) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			value = hashMap_get(templates, node->key);
			if(value) {
				++hits->n;
				hits->tot += node->value;
				
				kmerNode = malloc(sizeof(struct hashTable));
				if(!kmerNode) {
					ERROR();
				}
				if(SU) {
					values_s = (short unsigned *) value;
					kmerNode->value = smalloc(((*values_s) + 1) * sizeof(unsigned));
					kmerNode->value[0] = *values_s;
					j = kmerNode->value[0] + 1;
					while(--j) {
						template = values_s[j];
						kmerNode->value[j] = template;
						Scores[template]++;
						Scores_tot[template] += node->value;
					}
				} else {
					kmerNode->value = smalloc(((*value) + 1) * sizeof(unsigned));
					kmerNode->value[0] = *value;
					j = kmerNode->value[0] + 1;
					while(--j) {
						template = value[j];
						kmerNode->value[j] = template;
						Scores[template]++;
						Scores_tot[template] += node->value;
					}
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
	
	int template, SU;
	unsigned i, j, n, *value;
	short unsigned *value_s;
	struct hashTable_kmers *node, *node_next;
	struct hashTable *kmerNode, *kmerList;
	struct hashTable *decon_node, *deconList;
	struct hashTable **Returner;
	
	if(DB_size < HU_LIMIT) {
		SU = 1;
	} else {
		SU = 0;
	}
	value_s = 0;
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	deconList = 0;
	decon_node = 0;
	
	Returner = smalloc(sizeof(struct hashTable *) << 1);
	
	for(i = 0; i < foundKmers->size; ++i) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			if((value = hashMap_get(templates, node->key))) {
				/* check for contamination */
				++hits->n;
				hits->tot += node->value;
				
				if(SU) {
					value_s = (short unsigned *) value;
					n = *value_s;
					j = value_s[*value_s];
				} else {
					n = *value;
					j = value[*value];
				}
				
				if(j == contamination) {
					decon_node = malloc(sizeof(struct hashTable));
					if(!decon_node) {
						ERROR();
					}
					decon_node->value = malloc((n + 1) * sizeof(unsigned));
					if(!decon_node->value) {
						ERROR();
					}
					decon_node->value[0] = n;
					j = n + 1;
					if(SU) {
						while(--j) {
							decon_node->value[j] = value_s[j];
						}
					} else {
						while(--j) {
							decon_node->value[j] = value[j];
						}
					}
					decon_node->key = node->value;
					
					decon_node->next = deconList;
					deconList = decon_node;
				} else {
					kmerNode = malloc(sizeof(struct hashTable));
					if(!kmerNode) {
						ERROR();
					}
					kmerNode->value = malloc((n + 1) * sizeof(unsigned));
					if(!kmerNode->value) {
						ERROR();
					}
					kmerNode->value[0] = n;
					j = n + 1;
					if(SU) {
						while(--j) {
							template = value_s[j];
							kmerNode->value[j] = template;
							Scores[template]++;
							Scores_tot[template] += node->value;
						}
					} else {
						while(--j) {
							template = value[j];
							kmerNode->value[j] = template;
							Scores[template]++;
							Scores_tot[template] += node->value;
						}
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
	struct hashTable *node, *prev;
	prev = 0;
	
	if(kmerList == 0) {
		return 0;
	}
	
	node = kmerList;
	while(node != 0) {
		if(intpos_bin(node->value, template) != -1) {
			--hits->n;
			hits->tot -= node->key;
			for(i = 1; i <= node->value[0]; ++i) {
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
	
	unsigned i, belong;
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
			for(i = node->value[0]; i != 0 && !belong; --i) {
				if(node->value[i] == template) {
					belong = 1;
				}
			}
			if(belong) { //withdraw score
				--hits.n;
				hits.tot -= node->key;
				for(i = node->value[0]; i != 0; --i) {
					Scores[node->value[i]]--;
					Scores_tot[node->value[i]] -= node->key;
				}
				cont_node->value = node->value;
				cont_node->key = node->key;
				cont_node->next = malloc(sizeof(struct hashTable));
				if(!cont_node->next) {
					ERROR();
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

long unsigned run_input(char **inputfiles, int fileCount, int minPhred, int fiveClip, char *trans) {
	
	int fileCounter, phredCut, start, end;
	unsigned FASTQ;
	long unsigned count;
	char *filename;
	unsigned char *seq;
	struct qseqs *header, *qseq, *qual;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	
	compressor = malloc(sizeof(struct compDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(delta);
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
			phredCut = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			while(FileBuffgetFq(inputfile, header, qseq, qual, trans)) {
				/* trim */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1;
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
				if(qseq->len > kmersize) {
					/* dump seq */
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor);
					qseq->seq -= start;
					++count;
				}
			}
		} else if(FASTQ & 2) {
			while(FileBuffgetFsa(inputfile, header, qseq, trans)) {
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
				if(qseq->len > kmersize) {
					/* dump seq */
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor);
					qseq->seq -= start;
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
	destroyFileBuff(inputfile);
	
	return count;
}

long unsigned run_input_PE(char **inputfiles, int fileCount, int minPhred, int fiveClip, char *trans) {
	
	int fileCounter, phredCut, start, start2, end;
	unsigned FASTQ, FASTQ2;
	long unsigned count;
	char *filename;
	unsigned char *seq;
	struct qseqs *header, *qseq, *qual, *header2, *qseq2, *qual2;
	struct FileBuff *inputfile, *inputfile2;
	struct compDNA *compressor;
	
	
	compressor = malloc(sizeof(struct compDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(delta);
	qual = setQseqs(1024);
	header2 = setQseqs(256);
	qseq2 = setQseqs(delta);
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
			phredCut = getPhredFileBuff(inputfile);
			if(phredCut == 0) {
				phredCut = getPhredFileBuff(inputfile2);
			}
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			//while(FileBuffgetFq(inputfile, header, qseq, qual) && FileBuffgetFq(inputfile2, header2, qseq2, qual2)) {
			/* this ensures reading of truncated files */
			while((FileBuffgetFq(inputfile, header, qseq, qual, trans) | FileBuffgetFq(inputfile2, header2, qseq2, qual2, trans))) {
				/* trim forward */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1;
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
				
				/* print */
				if(qseq->len > kmersize && qseq2->len > kmersize) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > kmersize) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > kmersize) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor);
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
				if(qseq->len > kmersize && qseq2->len > kmersize) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > kmersize) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > kmersize) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor);
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

long unsigned run_input_INT(char **inputfiles, int fileCount, int minPhred, int fiveClip, char *trans) {
	
	int fileCounter, phredCut, start, start2, end;
	unsigned FASTQ;
	long unsigned count;
	char *filename;
	unsigned char *seq;
	struct qseqs *header, *qseq, *qual, *header2, *qseq2, *qual2;
	struct FileBuff *inputfile;
	struct compDNA *compressor;
	
	
	compressor = malloc(sizeof(struct compDNA));
	if(!compressor) {
		ERROR();
	}
	allocComp(compressor, 1024);
	header = setQseqs(256);
	qseq = setQseqs(delta);
	qual = setQseqs(1024);
	header2 = setQseqs(256);
	qseq2 = setQseqs(delta);
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
			phredCut = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			while((FileBuffgetFq(inputfile, header, qseq, qual, trans) | FileBuffgetFq(inputfile, header2, qseq2, qual2, trans))) {
				/* trim forward */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1;
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
				
				/* print */
				if(qseq->len > kmersize && qseq2->len > kmersize) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > kmersize) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > kmersize) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor);
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
				if(qseq->len > kmersize && qseq2->len > kmersize) {
					qseq->seq += start;
					qseq2->seq += start2;
					printFsa_pair_ptr(header, qseq, header2, qseq2, compressor);
					qseq->seq -= start;
					qseq2->seq -= start2;
					++count;
				} else if(qseq->len > kmersize) {
					qseq->seq += start;
					printFsa_ptr(header, qseq, compressor);
					qseq->seq -= start;
					++count;
				} else if(qseq2->len > kmersize) {
					qseq2->seq += start2;
					printFsa_ptr(header2, qseq2, compressor);
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

void run_input_sparse(char **inputfiles, int fileCount, int minPhred, int fiveClip, char *trans) {
	
	int FASTQ, fileCounter, phredCut, start, end;
	char *filename;
	unsigned char *seq;
	struct qseqs *qseq, *qual;
	struct FileBuff *inputfile;
	struct compKmers *Kmers;
	
	
	Kmers = malloc(sizeof(struct compKmers));
	if(!Kmers) {
		ERROR();
	}
	allocCompKmers(Kmers, delta);
	qseq = setQseqs(delta);
	qual = setQseqs(1024);
	inputfile = setFileBuff(CHUNK);
	
	for(fileCounter = 0; fileCounter < fileCount; ++fileCounter) {
		filename = (char*)(inputfiles[fileCounter]);
		
		/* determine filetype and open it */
		if((FASTQ = openAndDetermine(inputfile, filename)) & 3) {
			fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", filename);
		}
		
		/* parse the file */
		Kmers->n = 0;
		if(FASTQ & 1) {
			/* get phred scale */
			phredCut = getPhredFileBuff(inputfile);
			fprintf(stderr, "# Phred scale:\t%d\n", phredCut);
			phredCut += minPhred;
			
			/* parse reads */
			while(FileBuffgetFqSeq(inputfile, qseq, qual, trans)) {
				/* trim */
				seq = qual->seq;
				start = fiveClip;
				end = qseq->len - 1;
				while(end >= 0 && seq[end] < phredCut) {
					--end;
				}
				++end;
				while(start < end && seq[start] < phredCut) {
					++start;
				}
				qseq->len = end - start;
				
				/* print */
				if(qseq->len > kmersize) {
					/* translate to kmers */
					Kmers->n = translateToKmersAndDump(Kmers->kmers, Kmers->n, Kmers->size, qseq->seq + start, qseq->len, templates->prefix, templates->prefix_len);
				}
			}
			if(Kmers->n) {
				sfwrite(Kmers->kmers, sizeof(long unsigned), Kmers->n, stdout);
				Kmers->n = 0;
			}
		} else if(FASTQ & 2) {
			while(FileBuffgetFsaSeq(inputfile, qseq, trans)) {
				if(qseq->len > kmersize) {
					/* translate to kmers */
					Kmers->n = translateToKmersAndDump(Kmers->kmers, Kmers->n, Kmers->size, qseq->seq, qseq->len, templates->prefix, templates->prefix_len);
				}
			}
			if(Kmers->n) {
				sfwrite(Kmers->kmers, sizeof(long unsigned), Kmers->n, stdout);
				Kmers->n = 0;
			}
		}
		
		if(FASTQ & 4) {
			gzcloseFileBuff(inputfile);
		} else {
			closeFileBuff(inputfile);
		}
		
	}
	
	free(Kmers->kmers);
	free(Kmers);
	destroyQseqs(qseq);
	destroyQseqs(qual);
	destroyFileBuff(inputfile);
}

int save_kmers_batch(char *templatefilename, char *exePrev) {
	
	int i, file_len, shmid;
	FILE *inputfile, *templatefile;
	time_t t0, t1;
	key_t key;
	struct qseqs **Header;
	struct compDNA **Qseq, **Qseq_r;
	struct kmerScan_thread *threads, *thread;
	
	/* open pipe */
	//inputfile = popen(exePrev, "r");
	inputfile = kmaPopen(exePrev, "rb");
	if(!inputfile) {
		ERROR();
	}
	
	t0 = clock();
	
	/* initialize seqs */
	Qseq = malloc(thread_num * sizeof(struct compDNA *));
	Qseq_r = malloc(thread_num * sizeof(struct compDNA *));
	Header = malloc(thread_num * sizeof(struct qseqs *));
	if(!Qseq || !Qseq_r || !Header) {
		ERROR();
	}
	for(i = 0; i < thread_num; ++i) {
		Qseq[i] = malloc(sizeof(struct compDNA));
		Qseq_r[i] = malloc(sizeof(struct compDNA));
		Header[i] = setQseqs(256);
		if(!Qseq[i] || !Qseq_r[i]) {
			ERROR();
		}
		allocComp(Qseq[i], delta);
		allocComp(Qseq_r[i], delta);
	}
	
	/* load hashMap */
	file_len = strlen(templatefilename);
	if(deCon) {
		strcat(templatefilename, ".decon.comp.b");
	} else {
		strcat(templatefilename, ".comp.b");
	}
	templatefile = sfopen(templatefilename, "rb" );
	
	templates = malloc(sizeof(struct hashMapKMA));
	if(!templates) {
		ERROR();
	}
	
	if(diskDB) {
		getPrefix(templates, templatefile);
		fclose(templatefile);
		templates_offsets = getOffstets(templatefilename);
		if(1 || !one2one || diskDB & 2) {
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
					ERROR();
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
				ERROR();
			}
		}
	} else {
		hashMap_get = &hashMap_getGlobal;
		if((shm & 1) || (deCon && (shm & 2))) {
			hashMapKMA_load_shm(templates, templatefile, templatefilename);
		} else {
			if(hashMapKMA_load(templates, templatefile, templatefilename)) {
				fprintf(stderr, "Wrong format of DB.\n");
				exit(2);
			}
		}
		fclose(templatefile);
	}
	templatefilename[file_len] = 0;
	/* check if DB is sparse */
	if(templates->prefix_len != 0 || templates->prefix != 0) {
		/* set pointers to sparse detection */
		printPtr = &print_ankers_Sparse;
		if(deConPrintPtr != &deConPrint) {
			deConPrintPtr = printPtr;
		}
		one2one = 1;
		kmerScan = &save_kmers_Sparse;
		get_kmers_for_pair_ptr = &get_kmers_for_pair_Sparse;
	}
	
	/* make indexing a masking problem */
	--templates->size;
	
	kmersize = templates->kmersize;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	shifterS = sizeof(long unsigned) * sizeof(long unsigned) - (templates->prefix_len << 1);
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
	tScore = malloc(thread_num * sizeof(int *));
	tScore_r = malloc(thread_num * sizeof(int *));
	if(!tScore || !tScore_r) {
		ERROR();
	}
	for(i = 0; i < thread_num; ++i) {
		tScore[i] = calloc(DB_size, sizeof(int));
		tScore_r[i] = calloc(DB_size, sizeof(int));
		if(!tScore[i] || !tScore_r[i]) {
			ERROR();
		}
	}
	
	BestTemplates = malloc(thread_num * sizeof(int *));
	BestTemplates_r = malloc(thread_num * sizeof(int *));
	RegionTemplates = malloc(thread_num * sizeof(int *));
	if(!BestTemplates || !BestTemplates_r || !RegionTemplates) {
		ERROR();
	}
	for(i = 0; i < thread_num; ++i) {
		BestTemplates[i] = calloc((DB_size << 1) + 3, sizeof(int));
		BestTemplates_r[i] = calloc(DB_size + 3, sizeof(int));
		RegionTemplates[i] = calloc((DB_size << 1) + 1, sizeof(int));
		if(!BestTemplates[i] || !BestTemplates_r[i] || !RegionTemplates[i]) {
			ERROR();
		}
	}
	
	if(!one2one) {
		tVF_scores = malloc(thread_num * sizeof(int **));
		tVR_scores = malloc(thread_num * sizeof(int **));
		TmpNs = malloc(thread_num * sizeof(int *));
		if(!tVF_scores || !tVR_scores || !TmpNs) {
			ERROR();
		}
		for(i = 0; i < thread_num; ++i) {
			tVF_scores[i] = calloc(delta, sizeof(int*));
			tVR_scores[i] = calloc(delta, sizeof(int*));
			TmpNs[i] = malloc(delta * sizeof(int));
			if(!tVF_scores[i] || !tVR_scores[i] || !TmpNs[i]) {
				ERROR();
			}
		}
		
		/* load lengths */
		strcat(templatefilename, ".length.b");
		templatefile = sfopen(templatefilename, "rb");
		
		fread(&DB_size, sizeof(int), 1, templatefile);
		if(shm & 4) {
			key = ftok(templatefilename, 'l');
			shmid = shmget(key, DB_size * sizeof(int), 0666);
			if(shmid < 0) {
				fprintf(stderr, "No shared length\n");
				exit(2);
			} else {
				template_lengths = shmat(shmid, NULL, 0);
			}
		} else {
			template_lengths = malloc(DB_size * sizeof(int));
			if(!template_lengths) {
				ERROR();
			}
			fread(template_lengths, sizeof(int), DB_size, templatefile);
		}
		templatefilename[file_len] = 0;
		fclose(templatefile);
		minLen = template_lengths[1] / 2;
		for(i = 1; i < DB_size; ++i) {
			if(minLen > (template_lengths[i] / 2)) {
				minLen = template_lengths[i] / 2;
			}
		}
		minLen = (minLen < kmersize) ? kmersize : minLen;
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mer ankers\n");
	
	/* create exclusions */
	excludeIn = calloc(1, sizeof(int));
	excludeOut = calloc(1, sizeof(int));
	if(!excludeIn || !excludeOut) {
		ERROR();
	}
	
	/* initialize threads */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		
		thread = malloc(sizeof(struct kmerScan_thread));
		if(!thread) {
			ERROR();
		}
		thread->num = i;
		thread->bestScore = 0;
		thread->bestScore_r = 0;
		thread->bestTemplates = BestTemplates[i];
		thread->bestTemplates_r = BestTemplates_r[i];
		thread->qseq = Qseq[i];
		thread->qseq_r = Qseq_r[i];
		thread->header = Header[i];
		thread->inputfile = inputfile;
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, &save_kmers_threaded, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* start main thread */
	thread = malloc(sizeof(struct kmerScan_thread));
	if(!thread) {
		ERROR();
	}
	thread->num = 0;
	thread->bestScore = 0;
	thread->bestScore_r = 0;
	thread->bestTemplates = *BestTemplates;
	thread->bestTemplates_r = *BestTemplates_r;
	thread->qseq = *Qseq;
	thread->qseq_r = *Qseq_r;
	thread->header = *Header;
	thread->inputfile = inputfile;
	save_kmers_threaded(thread);
	
	/* print remaining buffer */
	if(printPtr == &print_ankers_spltDB) {
		print_ankers_spltDB(thread->bestTemplates, 0, 0, 0);
	}
	
	/* join threads */
	for(thread = threads; thread; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used ankering query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return kmaPclose(inputfile);
}

int save_kmers_sparse_batch(char *templatefilename, char *outputfilename, char *exePrev, char ss) {
	
	int i, file_len, stop, template, score, tmp_score;
	int score_add, score_tot_add, shmid, status;
	int *Scores, *Scores_tot, *w_Scores, *w_Scores_tot, *SearchList;
	unsigned Ntot;
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
	//inputfile = popen(exePrev, "r");
	inputfile = kmaPopen(exePrev, "rb");
	if(!inputfile) {
		ERROR();
	}
	
	/* Load hashMap */
	t0 = clock();
	file_len = strlen(templatefilename);
	if(deCon) {
		strcat(templatefilename, ".decon.comp.b");
	} else {
		strcat(templatefilename, ".comp.b");
	}
	templatefile = sfopen(templatefilename, "rb" );
	templates = malloc(sizeof(struct hashMapKMA));
	if(!templates) {
		ERROR();
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
					ERROR();
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
				ERROR();
			}
		}
	} else {
		hashMap_get = &hashMap_getGlobal;
		if((shm & 1) || (deCon && (shm & 2))) {
			hashMapKMA_load_shm(templates, templatefile, templatefilename);
		} else {
			if(hashMapKMA_load(templates, templatefile, templatefilename)) {
				fprintf(stderr, "Wrong format of DB.\n");
				exit(2);
			}
		}
		fclose(templatefile);
	}
	templatefilename[file_len] = 0;
	DB_size = contamination + 1;
	/* make indexing a masking problem */
	--templates->size;
	
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
		sparse_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
	}
	if(!sparse_out) {
		ERROR();
	}
	
	/* set hashMap for found kmers */
	foundKmers = malloc(sizeof(struct hashMap_kmers));
	if(!foundKmers) {
		ERROR();
	}
	
	foundKmers->size = templates->n;
	hashMap_kmers_initialize(foundKmers, foundKmers->size);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used for DB loading: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding k-mers\n");
	
	Kmers = malloc(sizeof(struct compKmers));
	if(!Kmers) {
		ERROR();
	}
	allocCompKmers(Kmers, delta);
	Ntot = 0;
	
	/* count kmers */
	while((Kmers->n = fread(Kmers->kmers, sizeof(long unsigned), Kmers->size, inputfile))) {
		Ntot += Kmers->n;
		save_kmers_sparse(foundKmers, Kmers);
	}
	status = kmaPclose(inputfile);
	
	t1 = clock();
	fprintf(stderr, "#\n# Total time used to identify k-mers in query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	t0 = clock();
	fprintf(stderr, "# Finding best matches and output results.\n");
	
	Scores = calloc(DB_size, sizeof(int));
	Scores_tot = calloc(DB_size, sizeof(int));
	SearchList = malloc(DB_size * sizeof(int));
	if(!Scores || !Scores_tot || !SearchList) {
		ERROR();
	}
	etta = 1.0e-6;
	expected = 0;
	q_value = 0;
	p_value = 0;
	
	fprintf(sparse_out, "#Template\tNum\tScore\tExpected\tTemplate_length\tQuery_Coverage\tTemplate_Coverage\tDepth\ttot_query_Coverage\ttot_template_Coverage\ttot_depth\tq_value\tp_value\n");
	stop = 0;
	if(deCon) {
		
		/* start by removing contamination and collect scores */
		Collecter = collect_Kmers_deCon(Scores, Scores_tot, foundKmers, &Nhits);
		kmerList = Collecter[0];
		deConTable = Collecter[1];
		free(Collecter);
		
		fprintf(stderr, "# total number of matches: %d of %d kmers\n", Nhits.tot, Ntot);
		/* copy scores */
		w_Scores = calloc(DB_size, sizeof(int));
		w_Scores_tot = calloc(DB_size, sizeof(int));
		if(!w_Scores || !w_Scores_tot) {
			ERROR();
		}
		
		for(i = 0; i < DB_size; ++i) {
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
			fprintf(stderr, " Failed at removing contamination\n");
			exit(1);
		}
		SearchList[contamination] = 0;
		
		/* get best matches */
		while(! stop) {
			/* get best match, depth first then coverage */
			depth = 0;
			cover = 0;
			score = 0;
			template = 0;
			if(ss == 'q') {
				for(i = 0; i < DB_size; ++i) {
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
				for(i = 0; i < DB_size; ++i) {
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
				for(i = 0; i < DB_size; ++i) {
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
			if(cover && cover >= ID_t) {
				/* with draw contamination k-mers matching this template */
				score_add = 0;
				score_tot_add = 0;
				if(deConTable != 0) {
					prev = 0;
					node = deConTable;
					while(node != 0) {
						if(intpos_bin(node->value, template) != -1) {
							++score_add;
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
			ERROR();
		}
		
		for(i = 0; i < DB_size; ++i) {
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
			template = 0;
			if(ss == 'q') {
				for(i = 0; i < DB_size; ++i) {
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
				for(i = 0; i < DB_size; ++i) {
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
				for(i = 0; i < DB_size; ++i) {
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
			if(cover && cover >= ID_t) {
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
	
	return status;
}

struct alnScore NW(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, struct aln *aligned, struct NWmat *matrices) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, template_length, pos[2];
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	struct alnScore Stat;
	
	q_len = q_e - q_s;
	t_len = t_e - t_s;
	template_length = aligned->pos;
	if(t_len < 0) {
		t_len += aligned->pos;
	}
	query = (unsigned char*)(queryOrg + q_s);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.gaps = 0;
			Stat.score = 0;
			aligned->s[0] = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.gaps = q_len;
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
			m = t_len;
			nuc_pos = t_e - 1;
			while(m--) {
				aligned->t[m] = getNuc(template, nuc_pos);
				if(--nuc_pos < 0) {
					nuc_pos = aligned->pos - 1;
				}
			}
		}
		return Stat;
	}
	
	/* check matrix size */
	if(matrices->NW_q <= q_len) {
		matrices->NW_q = q_len << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((q_len + 1) * (t_len + 1))) {
		matrices->NW_s = ((q_len + 2) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	Stat.pos = 0;
	if(0 < k) {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 0;
			E_ptr += (q_len + 1);
		}
		if(k == 1) {
			for(n = q_len - 1; n >= 0; --n) {
				D_prev[n] = W1 + (q_len - 1 - n) * U;
				P_prev[n] = thisScore;
				E_ptr[n] = 3;
			}
			E_ptr[q_len - 1] = 18;
			E_ptr[q_len] = 0;
			D_prev[q_len] = 0;
			P_prev[q_len] = 0;
		} else {
			for(n = q_len; n >= 0; --n) {
				D_prev[n] = 0;
				P_prev[n] = thisScore;
				E_ptr[n] = 0;
			}
		}
	} else {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 5;
			E_ptr += (q_len + 1);
		}
		E_ptr[-1] = 36;
		
		for(n = q_len - 1; n >= 0; --n) {
			D_prev[n] = W1 + (q_len - 1 - n) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[q_len - 1] = 18;
		E_ptr[q_len] = 0;
		D_prev[q_len] = 0;
		P_prev[q_len] = 0;
	}
	E_ptr -= (q_len + 1);
	
	/* Perform NW */
	pos[0] = 0;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1;
		}
		
		D_ptr[q_len] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
		Q_prev = (t_len + q_len) * (MM + U + W1);
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = q_len - 1; n >= 0; --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(e == 2) {
					D_ptr[n] = Q;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] < thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n + 1] + d[t_nuc][query[n]];
			if(D_ptr[n] < thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		E_ptr -= (q_len + 1);
		
		if(k < 0 && Stat.score <= *D_ptr) {
			Stat.score = *D_ptr;
			pos[0] = m;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	if(k < 0) {
		pos[1] = 0;
		if(k == -2) {
			for(n = 0; n < q_len; ++n) {
				if(D_prev[n] > Stat.score) {
					Stat.score = D_prev[n];
					pos[0] = 0;
					pos[1] = n;
				}
			}
		}
	} else {
		Stat.score = *D_prev;
		pos[0] = 0;
		pos[1] = 0;
	}
	
	/* make back tracking */
	m = pos[0];
	E_ptr = E + (m * (q_len + 1));
	n = pos[1];
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.gaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = (aligned->t[Stat.len] == aligned->q[Stat.len]) ? '|' : '_';
			++nuc_pos;
			E_ptr += (q_len + 1);
			++n;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				aligned->t[Stat.len] = getNuc(template, nuc_pos);
				aligned->q[Stat.len] = 5;
				aligned->s[Stat.len] = '_';	
				++nuc_pos;
				E_ptr += (q_len + 1);
				++Stat.len;
			}
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = 5;
			aligned->s[Stat.len] = '_';		
			++nuc_pos;
			E_ptr += (q_len + 1);
		} else {
			while(!(E_ptr[n] >> 3)) {
				aligned->t[Stat.len] = 5;
				aligned->q[Stat.len] = query[n];
				aligned->s[Stat.len] = '_';
				++Stat.gaps;
				++n;
				++Stat.len;
			}
			aligned->t[Stat.len] = 5;
			aligned->q[Stat.len] = query[n];
			aligned->s[Stat.len] = '_';
			++Stat.gaps;
			++n;
		}
		++Stat.len;
	}
	aligned->s[Stat.len] = 0;
	
	return Stat;
}

struct alnScore NW_band(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, struct aln *aligned, int band, struct NWmat *matrices) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, template_length, pos[2];
	int bq_len, halfBand, sn, en, sq, eq, q_pos, c_pos;
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	struct alnScore Stat;
	
	q_len = q_e - q_s;
	t_len = t_e - t_s;
	template_length = aligned->pos;
	if(t_len < 0) {
		t_len += aligned->pos;
	}
	query = (unsigned char*)(queryOrg + q_s);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.gaps = 0;
			Stat.score = 0;
			aligned->s[0] = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.gaps = q_len;
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
			m = t_len;
			nuc_pos = t_e - 1;
			while(m--) {
				aligned->t[m] = getNuc(template, nuc_pos);
				if(--nuc_pos < 0) {
					nuc_pos = aligned->pos - 1;
				}
			}
		}
		return Stat;
	}
	
	/* ensure that band is equal */
	if(band & 1) {
		++band;
	}
	halfBand = band >> 1;
	
	/* check matrix size */
	if(matrices->NW_q <= band) {
		matrices->NW_q = band << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((band + 2) * (t_len + 1))) {
		matrices->NW_s = ((band + 3) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	bq_len = band + 1; /* (band + 1) ~ q_len */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	Stat.pos = 0;
	E_ptr = E + (t_len * (bq_len + 1));
	c_pos = (t_len + q_len) >> 1;
	
	
	sn = q_len - 1 - (c_pos - halfBand);
	if(k != 2) {
		for(n = sn - 1; n >= 0; --n) {
			D_prev[n] = W1 + (sn - n - 1) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[sn - 1] = 18;
		E_ptr[sn] = 0;
		D_prev[sn] = 0;
		P_prev[sn] = 0;
	} else {
		for(n = sn; n >= 0; --n) {
			D_prev[n] = 0;
			P_prev[n] = thisScore;
			E_ptr[n] = 0;
		}
	}
	E_ptr -= (bq_len + 1);
	
	
	/* Perform banded NW */
	pos[0] = 0;
	en = 0;
	c_pos = (t_len + q_len) >> 1;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos, --c_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1; 
		}
		
		/* get banded boundaries, w.r.t. query */
		sq = c_pos + halfBand;
		eq = c_pos - halfBand;
		
		if(eq < 0) {
			eq = 0;
			++en;
		} else {
			en = 0;
		}
		
		/* get start penalties */
		Q_prev = (t_len + q_len) * (MM + U + W1);
		/* check boundaries */
		if(sq < (q_len - 1)) {
			sn = bq_len - 1;
			D_ptr[bq_len] = (t_len + q_len) * (MM + U + W1);
			E_ptr[bq_len] = 37;
		} else {
			sq = q_len - 1;
			sn = en + (q_len - eq);
			D_ptr[sn] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
			E_ptr[sn] = (0 < k) ? 0 : 37;
			--sn;
		}
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = sn, q_pos = sq; n > en; --q_pos, --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n - 1] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(e == 2) {
					D_ptr[n] = Q;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n - 1] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] < thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n] + d[t_nuc][query[q_pos]];
			if(D_ptr[n] < thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		/* handle banded boundary */
		E_ptr[n] = 0;
		
		/* update Q gap */
		Q = D_ptr[n + 1] + W1;
		thisScore = Q_prev + U;
		if(Q < thisScore) {
			Q = thisScore;
			e = 3;
		} else {
			e = 2;
			E_ptr[n] |= 16;
		}
		
		/* update unavailable P gap */
		P_ptr[n] = (t_len + q_len) * (MM + U + W1);
		
		/* Update D */
		D_ptr[n] = D_prev[n] + d[t_nuc][query[q_pos]];
		
		/* set D to max, and set E */
		if(Q < D_ptr[n]) {
			E_ptr[n] |= 1;
		} else {
			D_ptr[n] = Q;
			E_ptr[n] |= e;
		}
		
		/* continue as usual */
		E_ptr -= (bq_len + 1);
		
		if(k < 0 && Stat.score <= D_ptr[n]) {
			Stat.score = D_ptr[n];
			pos[0] = m;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	if(k < 0) {
		pos[1] = 0;
		q_pos = 0;
		if(k == -2) {
			for(n = en; n < bq_len; ++n) {
				if(D_prev[n] > Stat.score) {
					Stat.score = D_prev[n];
					pos[0] = 0;
					pos[1] = n;
				}
			}
		}
	} else {
		Stat.score = D_prev[en];
		q_pos = 0;
		pos[0] = 0;
		pos[1] = en;
	}
	
	
	/* back tracking */
	m = pos[0];
	n = pos[1];
	E_ptr = E + (m * (bq_len + 1));
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.gaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = query[q_pos];
			aligned->s[Stat.len] = (aligned->t[Stat.len] == aligned->q[Stat.len]) ? '|' : '_';
			++nuc_pos;
			E_ptr += (bq_len + 1);
			++q_pos;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				aligned->t[Stat.len] = getNuc(template, nuc_pos);
				aligned->q[Stat.len] = 5;
				aligned->s[Stat.len] = '_';	
				++nuc_pos;
				E_ptr += (bq_len + 1);
				--n;
				++Stat.len;
			}
			aligned->t[Stat.len] = getNuc(template, nuc_pos);
			aligned->q[Stat.len] = 5;
			aligned->s[Stat.len] = '_';		
			++nuc_pos;
			E_ptr += (bq_len + 1);
			--n;
		} else {
			while(!(E_ptr[n] >> 3)) {
				aligned->t[Stat.len] = 5;
				aligned->q[Stat.len] = query[q_pos];
				aligned->s[Stat.len] = '_';
				++Stat.gaps;
				++n;
				++q_pos;
				++Stat.len;
			}
			aligned->t[Stat.len] = 5;
			aligned->q[Stat.len] = query[q_pos];
			aligned->s[Stat.len] = '_';
			++Stat.gaps;
			++n;
			++q_pos;
		}
		++Stat.len;
	}
	aligned->s[Stat.len] = 0;
	
	return Stat;
}

struct alnScore NW_score(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, struct NWmat *matrices, int template_length) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, pos[2];
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	struct alnScore Stat;
	
	t_len = t_e - t_s;
	q_len = q_e - q_s;
	if(t_len < 0) {
		t_len += template_length;
	}
	query = (unsigned char*)(queryOrg + q_s);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.gaps = 0;
			Stat.score = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.gaps = q_len;
			Stat.score = W1 + (q_len - 1) * U;
		} else {
			Stat.len = t_len;
			Stat.gaps = 0;
			Stat.score = W1 + (t_len - 1) * U;
		}
		return Stat;
	}
	
	/* check matrix size */
	if(matrices->NW_q <= q_len) {
		matrices->NW_q = q_len << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((q_len + 1) * (t_len + 1))) {
		matrices->NW_s = ((q_len + 2) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	Stat.pos = 0;
	if(0 < k) {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 0;
			E_ptr += (q_len + 1);
		}
		if(k == 1) {
			for(n = q_len - 1; n >= 0; --n) {
				D_prev[n] = W1 + (q_len - 1 - n) * U;
				P_prev[n] = thisScore;
				E_ptr[n] = 3;
			}
			E_ptr[q_len - 1] = 18;
			E_ptr[q_len] = 0;
			D_prev[q_len] = 0;
			P_prev[q_len] = 0;
		} else {
			for(n = q_len; n >= 0; --n) {
				D_prev[n] = 0;
				P_prev[n] = thisScore;
				E_ptr[n] = 0;
			}
		}
	} else {
		E_ptr = E;
		for(m = 0; m < t_len; ++m) {
			E_ptr[q_len] = 5;
			E_ptr += (q_len + 1);
		}
		E_ptr[-1] = 36;
		
		for(n = q_len - 1; n >= 0; --n) {
			D_prev[n] = W1 + (q_len - 1 - n) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[q_len - 1] = 18;
		E_ptr[q_len] = 0;
		D_prev[q_len] = 0;
		P_prev[q_len] = 0;
	}
	E_ptr -= (q_len + 1);
	
	/* Perform NW */
	pos[0] = 0;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1;
		}
		
		D_ptr[q_len] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
		Q_prev = (t_len + q_len) * (MM + U + W1);
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = q_len - 1; n >= 0; --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(e == 2) {
					D_ptr[n] = Q;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] < thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n + 1] + d[t_nuc][query[n]];
			if(D_ptr[n] < thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		E_ptr -= (q_len + 1);
		
		if(k < 0 && Stat.score <= *D_ptr) {
			Stat.score = *D_ptr;
			pos[0] = m;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	if(k < 0) {
		pos[1] = 0;
		if(k == -2) {
			for(n = 0; n < q_len; ++n) {
				if(D_prev[n] > Stat.score) {
					Stat.score = D_prev[n];
					pos[0] = 0;
					pos[1] = n;
				}
			}
		}
	} else {
		Stat.score = *D_prev;
		pos[0] = 0;
		pos[1] = 0;
	}
	
	/* make back tracking */
	m = pos[0];
	E_ptr = E + (m * (q_len + 1));
	n = pos[1];
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.gaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			++nuc_pos;
			E_ptr += (q_len + 1);
			++n;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				++nuc_pos;
				E_ptr += (q_len + 1);
				++Stat.len;
			}
			++nuc_pos;
			E_ptr += (q_len + 1);
		} else {
			while(!(E_ptr[n] >> 3)) {
				++Stat.gaps;
				++n;
				++Stat.len;
			}
			++Stat.gaps;
			++n;
		}
		++Stat.len;
	}
	
	return Stat;
}

struct alnScore NW_band_score(const long unsigned *template, const unsigned char *queryOrg, int k, int t_s, int t_e, int q_s, int q_e, int band, struct NWmat *matrices, int template_length) {
	
	int m, n, t_len, q_len, thisScore, nuc_pos, pos[2];
	int bq_len, halfBand, sn, en, sq, eq, q_pos, c_pos;
	int *D_ptr, *D_prev, Q, Q_prev, *P_ptr, *P_prev, *tmp;
	unsigned char *query, t_nuc, *E, *E_ptr, e;
	struct alnScore Stat;
	
	t_len = t_e - t_s;
	q_len = q_e - q_s;
	if(t_len < 0) {
		t_len += template_length;
	}
	query = (unsigned char*)(queryOrg + q_s);
	
	if(t_len == 0 || q_len == 0) {
		if(t_len == q_len) {
			Stat.len = 0;
			Stat.gaps = 0;
			Stat.score = 0;
		} else if(t_len == 0) {
			Stat.len = q_len;
			Stat.gaps = q_len;
			Stat.score = W1 + (q_len - 1) * U;
		} else {
			Stat.len = t_len;
			Stat.gaps = 0;
			Stat.score = W1 + (t_len - 1) * U;
		}
		return Stat;
	}
	
	/* ensure that band is equal */
	if(band & 1) {
		++band;
	}
	halfBand = band >> 1;
	
	/* check matrix size */
	if(matrices->NW_q <= band) {
		matrices->NW_q = band << 1;
		free(matrices->D[0]);
		free(matrices->P[0]);
		matrices->D[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->P[0] = smalloc((matrices->NW_q << 1) * sizeof(int));
		matrices->D[1] = matrices->D[0] + matrices->NW_q;
		matrices->P[1] = matrices->P[0] + matrices->NW_q;
	}
	if(matrices->NW_s <= ((band + 2) * (t_len + 1))) {
		matrices->NW_s = ((band + 3) * (t_len + 2));
		free(matrices->E);
		matrices->E = smalloc(matrices->NW_s);
	}
	
	/* fill in start penalties */
	bq_len = band + 1; /* (band + 1) ~ q_len */
	D_ptr = matrices->D[0];
	D_prev = matrices->D[1];
	P_ptr = matrices->P[0];
	P_prev = matrices->P[1];
	E = matrices->E;
	thisScore = (t_len + q_len) * (MM + U + W1);
	Stat.score = thisScore;
	Stat.pos = 0;
	E_ptr = E + (t_len * (bq_len + 1));
	c_pos = (t_len + q_len) >> 1;
	
	
	sn = q_len - 1 - (c_pos - halfBand);
	if(k != 2) {
		for(n = sn - 1; n >= 0; --n) {
			D_prev[n] = W1 + (sn - n - 1) * U;
			P_prev[n] = thisScore;
			E_ptr[n] = 3;
		}
		E_ptr[sn - 1] = 18;
		E_ptr[sn] = 0;
		D_prev[sn] = 0;
		P_prev[sn] = 0;
	} else {
		for(n = sn; n >= 0; --n) {
			D_prev[n] = 0;
			P_prev[n] = thisScore;
			E_ptr[n] = 0;
		}
	}
	E_ptr -= (bq_len + 1);
	
	
	/* Perform banded NW */
	pos[0] = 0;
	en = 0;
	c_pos = (t_len + q_len) >> 1;
	for(m = t_len - 1, nuc_pos = t_e - 1; m >= 0; --m, --nuc_pos, --c_pos) {
		
		if(nuc_pos < 0) {
			nuc_pos = template_length - 1; 
		}
		
		/* get banded boundaries, w.r.t. query */
		sq = c_pos + halfBand;
		eq = c_pos - halfBand;
		
		if(eq < 0) {
			eq = 0;
			++en;
		} else {
			en = 0;
		}
		
		/* get start penalties */
		Q_prev = (t_len + q_len) * (MM + U + W1);
		/* check boundaries */
		if(sq < (q_len - 1)) {
			sn = bq_len - 1;
			D_ptr[bq_len] = (t_len + q_len) * (MM + U + W1);
			E_ptr[bq_len] = 37;
		} else {
			sq = q_len - 1;
			sn = en + (q_len - eq);
			D_ptr[sn] = (0 < k) ? 0 : (W1 + (t_len - 1 - m) * U);
			E_ptr[sn] = (0 < k) ? 0 : 37;
			--sn;
		}
		
		t_nuc = getNuc(template, nuc_pos);
		for(n = sn, q_pos = sq; n > en; --q_pos, --n) {
			E_ptr[n] = 0;
			
			/* update Q and P, gap openings */
			Q = D_ptr[n + 1] + W1;
			P_ptr[n] = D_prev[n - 1] + W1;
			if(Q < P_ptr[n]) {
				D_ptr[n] = P_ptr[n];
				e = 4;
			} else {
				D_ptr[n] = Q;
				e = 2;
			}
			
			/* update Q and P, gap extensions */
			/* mark bit 4 and 5 as possible gap-opennings, if necesarry */
			thisScore = Q_prev + U;
			if(Q < thisScore) {
				Q = thisScore;
				if(e == 2) {
					D_ptr[n] = Q;
					e = 3;
				}
			} else {
				E_ptr[n] |= 16;
			}
			thisScore = P_prev[n - 1] + U;
			if(P_ptr[n] < thisScore) {
				P_ptr[n] = thisScore;
				if(D_ptr[n] < thisScore) {
					D_ptr[n] = thisScore;
					e = 5;
				}
			} else {
				E_ptr[n] |= 32;
			}
			
			/* Update D, match */
			thisScore = D_prev[n] + d[t_nuc][query[q_pos]];
			if(D_ptr[n] < thisScore) {
				D_ptr[n] = thisScore;
				E_ptr[n] |= 1;
			} else {
				E_ptr[n] |= e;
			}
			
			Q_prev = Q;
		}
		
		/* handle banded boundary */
		E_ptr[n] = 0;
		
		/* update Q gap */
		Q = D_ptr[n + 1] + W1;
		thisScore = Q_prev + U;
		if(Q < thisScore) {
			Q = thisScore;
			e = 3;
		} else {
			e = 2;
			E_ptr[n] |= 16;
		}
		
		/* update unavailable P gap */
		P_ptr[n] = (t_len + q_len) * (MM + U + W1);
		
		/* Update D */
		D_ptr[n] = D_prev[n] + d[t_nuc][query[q_pos]];
		
		/* set D to max, and set E */
		if(Q < D_ptr[n]) {
			E_ptr[n] |= 1;
		} else {
			D_ptr[n] = Q;
			E_ptr[n] |= e;
		}
		
		/* continue as usual */
		E_ptr -= (bq_len + 1);
		
		if(k < 0 && Stat.score <= D_ptr[n]) {
			Stat.score = D_ptr[n];
			pos[0] = m;
		}
		
		tmp = D_ptr;
		D_ptr = D_prev;
		D_prev = tmp;
		
		tmp = P_ptr;
		P_ptr = P_prev;
		P_prev = tmp;
	}
	E_ptr = E;
	
	/* get start position of alignment */
	if(k < 0) {
		pos[1] = 0;
		q_pos = 0;
		if(k == -2) {
			for(n = en; n < bq_len; ++n) {
				if(D_prev[n] > Stat.score) {
					Stat.score = D_prev[n];
					pos[0] = 0;
					pos[1] = n;
				}
			}
		}
	} else {
		Stat.score = D_prev[en];
		q_pos = 0;
		pos[0] = 0;
		pos[1] = en;
	}
	
	
	/* back tracking */
	m = pos[0];
	n = pos[1];
	E_ptr = E + (m * (bq_len + 1));
	nuc_pos = m + t_s;
	Stat.len = 0;
	Stat.gaps = 0;
	while(E_ptr[n] != 0) {
		if(nuc_pos == template_length) {
			nuc_pos = 0;
		}
		if((E_ptr[n] & 7) == 1) {
			++nuc_pos;
			E_ptr += (bq_len + 1);
			++q_pos;
		} else if((E_ptr[n] & 7) >= 4) {
			while(!(E_ptr[n] >> 4)) {
				++nuc_pos;
				E_ptr += (bq_len + 1);
				--n;
				++Stat.len;
			}
			++nuc_pos;
			E_ptr += (bq_len + 1);
			--n;
		} else {
			while(!(E_ptr[n] >> 3)) {
				++Stat.gaps;
				++n;
				++q_pos;
				++Stat.len;
			}
			++Stat.gaps;
			++n;
			++q_pos;
		}
		++Stat.len;
	}
	
	return Stat;
}

struct alnPoints * seedPoint_init(int size) {
	
	struct alnPoints *dest;
	
	dest = malloc(sizeof(struct alnPoints));
	if(!dest) {
		ERROR();
	}
	dest->len = 0;
	dest->size = size;
	size *= sizeof(int);
	dest->tStart = malloc(size);
	dest->tEnd = malloc(size);
	dest->qStart = malloc(size);
	dest->qEnd = malloc(size);
	dest->weight = malloc(size);
	dest->score = malloc(size);
	dest->next = malloc(size);
	if(!dest->tStart || !dest->tEnd || !dest->qStart || !dest->qEnd || !dest->weight || !dest->score || !dest->next) {
		ERROR();
	}
	
	return dest;
}

void seedPoint_realloc(struct alnPoints *dest, int size) {
	
	dest->size = size;
	size *= sizeof(int);
	dest->tStart = realloc(dest->tStart, size);
	dest->tEnd = realloc(dest->tEnd, size);
	dest->qStart = realloc(dest->qStart, size);
	dest->qEnd = realloc(dest->qEnd, size);
	dest->weight = realloc(dest->weight, size);
	dest->score = realloc(dest->score, size);
	dest->next = realloc(dest->next, size);
	if(!dest->tStart || !dest->tEnd || !dest->qStart || !dest->qEnd || !dest->weight || !dest->score || !dest->next) {
		ERROR();
	}
}

void seedPoint_free(struct alnPoints *src) {
	
	free(src->tStart);
	free(src->tEnd);
	free(src->qStart);
	free(src->qEnd);
	free(src->weight);
	free(src->score);
	free(src->next);
	free(src);
}

int chainSeeds(struct alnPoints *points, int q_len, int t_len, unsigned *mapQ) {
	
	int i, j, nMems, weight, gap, score, bestScore, secondScore, bestPos;
	int tStart, tEnd, qEnd, tGap, qGap, nMin;
	
	//pM = MM + M;
	nMems = points->len;
	i = nMems - 1;
	bestPos = i;
	bestScore = points->weight[i] * M;
	secondScore = 0;
	points->score[i] = bestScore;
	points->next[i] = 0;
	points->weight[i] -= (kmersize - 1);
	
	while(i--) {
		weight = points->weight[i] * M;
		points->next[i] = 0;
		tEnd = points->tEnd[i];
		qEnd = points->qEnd[i];
		score = MIN(t_len - tEnd, q_len - qEnd);
		if(0 < --score) {
			score *= U;
			score += W1;
		} else if(score == 0) {
			score = W1;
		} else {
			score = 0;
		}
		score += weight;
		
		/* 64 is the bandwidth */
		nMin = MIN(nMems, i + 64);
		
		/* find best link */
		for(j = i + 1; j < nMin; ++j) {
			/* check compability */
			if(qEnd < points->qStart[j]) {
				if(tEnd < points->tStart[j]) { /* full compatability */
					tGap = points->tStart[j] - tEnd;
					qGap = points->qStart[j] - qEnd;
					
					/* calculate score for this chaining */
					if((gap = abs(tGap - qGap))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					gap += weight + points->score[j];
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} else if(tEnd < points->tEnd[j]) { /* semi compatability */
					/* calculate score for this chaining */
					if((gap = (points->qStart[j] - qEnd))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score */	
					gap += (weight + points->score[j] - (points->tStart[j] - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				}
			} else if(kmersize <= points->qEnd[j] - qEnd) {
				/* cut in query is reflected in template */
				tStart = points->tStart[j] + qEnd - points->qStart[j];
				if(tEnd < tStart) { /* full compatability */
					tGap = tStart - tEnd;
					
					/* calculate score for this chaining */
					if((gap = tGap)) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score */
					gap += (weight + points->score[j] - (tStart - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} /* if the modified tStart was invalid, then the entire mem is lost. */
			}
		}
		/* update seed */
		if(points->next[i]) {
			points->weight[i] += (points->weight[points->next[i]] - kmersize + 1);
		} else {
			points->weight[i] -= (kmersize - 1);
		}
		points->score[i] = score;
		
		/* update bestScore */
		if(bestScore < score) {
			if(points->next[i] != bestPos) {
				secondScore = bestScore;
			}
			bestScore = score;
			bestPos = i;
		} else if(bestScore == score && points->next[i] != bestPos) {
			secondScore = bestScore;
		}
	}
	/* calculate mapping quality */
	*mapQ = ceil(40 * (1 - 1.0 * secondScore / bestScore) * MIN(1, points->weight[bestPos] / 10.0) * log(bestScore));
	
	/* penalize start */
	/*
	if(bestPos != 0) {
		score = (MIN(points->qStart[bestPos], points->tStart[bestPos])) * pM;
		bestScore += MAX(W1, score);
		if(points->qStart[0] == 0 && bestScore < points->score[0]) {
			bestScore = points->score[0];
			bestPos = 0;
		}
	}
	*/
	return bestPos;
}

int chainSeeds_circular(struct alnPoints *points, int q_len, int t_len, unsigned *mapQ) {
	
	int i, j, nMems, weight, gap, score, bestScore, secondScore, bestPos;
	int tStart, tEnd, qEnd, tGap, qGap, nMin;
	
	//pM = MM + M;
	nMems = points->len;
	i = nMems - 1;
	bestPos = i;
	bestScore = points->weight[i] * M;
	secondScore = 0;
	points->score[i] = bestScore;
	points->next[i] = 0;
	points->weight[i] -= (kmersize - 1);
	
	while(i--) {
		weight = points->weight[i] * M;
		points->next[i] = 0;
		tEnd = points->tEnd[i];
		qEnd = points->qEnd[i];
		score = MIN(t_len - tEnd, q_len - qEnd);
		if(0 < --score) {
			score *= U;
			score += W1;
		} else if(score == 0) {
			score = W1;
		} else {
			score = 0;
		}
		score += weight;
		
		/* 64 is the bandwidth */
		nMin = MIN(nMems, i + 64);
		
		/* find best link */
		for(j = i + 1; j < nMin; ++j) {
			/* check compability */
			if(qEnd < points->qStart[j]) {
				if(tEnd < points->tStart[j]) { /* full compatability */
					tGap = points->tStart[j] - tEnd;
					qGap = points->qStart[j] - qEnd;
					
					/* calculate score for this chaining */
					if((gap = abs(tGap - qGap))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					gap += weight + points->score[j];
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} else if(tEnd < points->tEnd[j]) { /* semi compatability */
					/* calculate score for this chaining */
					if((gap = (points->qStart[j] - qEnd))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score*/	
					gap += (weight + points->score[j] - (points->tStart[j] - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} else { /* circular joinning */
					tGap = t_len - tEnd + points->tStart[j];
					qGap = points->qStart[j] - qEnd;
					
					/* calculate score for this chaining */
					if((gap = abs(tGap - qGap))) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					gap += weight + points->score[j];
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				}
			} else if(kmersize <= points->qEnd[j] - qEnd) {
				/* cut in query is reflected in template */
				tStart = points->tStart[j] + qEnd - points->qStart[j];
				if(t_len < tStart) {
					tStart -= t_len;
				}
				if(tEnd < tStart) { /* full compatability */
					tGap = tStart - tEnd;
					
					/* calculate score for this chaining */
					if((gap = tGap)) {
						--gap;
						gap *= U;
						gap += W1;
					}
					/* cut mem score */
					gap += (weight + points->score[j] - (tStart - tEnd) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				} else if(tEnd != tStart) { /* circular joining */
					tGap = t_len - tEnd + tStart;
					
					/* calculate score for this chaining */
					if((gap = tGap)) {
						--gap;
						gap *= U;
						gap += W1;
					}
					//gap += ((MIN(tGap, qGap)) * pM + weight + points->score[j]);
					gap += (weight + points->score[j] - (tEnd - tStart) * M);
					
					/* check if score is max */
					if(score < gap) {
						score = gap;
						points->next[i] = j;
					}
				}/* if the modified tStart was invalid, then the entire mem is lost. */
			}
		}
		/* update seed */
		if(points->next[i]) {
			points->weight[i] += (points->weight[points->next[i]] - kmersize + 1);
		} else {
			points->weight[i] -= (kmersize - 1);
		}
		points->score[i] = score;
		
		/* update bestScore */
		if(bestScore < score) {
			if(points->next[i] != bestPos) {
				secondScore = bestScore;
			}
			bestScore = score;
			bestPos = i;
		} else if(bestScore == score && points->next[i] != bestPos) {
			secondScore = bestScore;
		}
	}
	/* calculate mapping quality */
	*mapQ = ceil(40 * (1 - 1.0 * secondScore / bestScore) * MIN(1, points->weight[bestPos] / 10.0) * log(bestScore));
	
	/* penalize start */
	/*
	if(bestPos != 0) {
		score = (MIN(points->qStart[bestPos], points->tStart[bestPos])) * pM;
		bestScore += MAX(W1, score);
		if(points->qStart[0] == 0 && bestScore < points->score[0]) {
			bestScore = points->score[0];
			bestPos = 0;
		}
	}
	*/
	
	return bestPos;
}

struct alnScore KMA(const int template_name, const unsigned char *qseq, int q_len, struct aln *aligned, struct aln *Frag_align, int min, int max, struct alnPoints *points, struct NWmat *matrices) {
	
	int i, j, k, bias, prev, start, stop, t_len, value, end, band;
	int t_l, t_s, t_e, q_s, q_e, mem_count, score;
	long unsigned key;
	unsigned char nuc;
	struct alnScore Stat, NWstat;
	struct hashMap_index *template_index;
	
	/* Extract indexes and template sequence */
	template_index = templates_index[template_name];
	t_len = template_lengths[template_name];
	key = 0;
	
	/* find seeds */
	if(points->len) {
		mem_count = points->len;
	} else {
		mem_count = 0;
		i = 0;
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
			
			while(i < end) {
				key = ((key << 2) | qseq[i]) & mask;
				value = hashMap_index_get_bound(template_index, key, min, max);
				
				if(value == 0) {
					++i;
				} else if(0 < value) {
					i -= (kmersize - 1);
					
					/* backseed for overlapping seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
					}
					
					/* get start positions */
					points->qStart[mem_count] = j + 1;
					points->tStart[mem_count] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
					}
					
					/* get end positions */
					points->qEnd[mem_count] = i;
					points->tEnd[mem_count] = value + 1;
					
					/* calculate weight */
					points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
					++mem_count;
					
					/* realloc seeding points */
					if(mem_count == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				} else {
					/* move counter back */
					i -= (kmersize - 1);
					
					/* get position in hashmap */
					stop = hashMap_index_getDubPos(template_index, key, value);
					
					/* get all mems */
					bias = i;
					while(0 <= stop) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(template_index->index[stop]);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[mem_count] = j + 1;
						points->tStart[mem_count] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[mem_count] = k;
						points->tEnd[mem_count] = value + 1;
						
						/* calculate weight */
						points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
						++mem_count;
						
						/* realloc seeding points */
						if(mem_count == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						stop = hashMap_index_getNextDubPos(template_index, key, min, max, stop);
					}
					i = bias + 1;
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
					
				}
			}
			i = end + 1;
		}
	}
	
	if(mem_count) {
		points->len = mem_count;
	} else {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		aligned->s[0] = 0;
		points->len = 0;
		return Stat;
	}
	
	/* get best seed chain, returns best starting point */
	start = chainSeedsPtr(points, q_len, t_len, &aligned->mapQ);
	score = points->score[start];
	
	//if(score < kmersize || (score << 1) < scoreT * q_len) {
	if(aligned->mapQ < mq || score < kmersize || (score << 1) < scoreT * (points->qEnd[points->len - 1] - points->qStart[0])) {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		aligned->s[0] = 0;
		points->len = 0;
		
		return Stat;
	}
	
	/* initialize */
	Stat.len = 0;
	Stat.score = 0;
	Stat.gaps = 0;
	value = points->tStart[start] - 1;
	Stat.pos = value;
	i = points->qStart[start];
	/* align leading tail */
	if(i != 0) {
		/* get boundaries */
		t_s = 0;
		t_e = value;
		q_s = 0;
		q_e = i;
		if((q_e << 1) < t_e || (q_e + 64) < t_e) { // big leading template gap, cut down
			t_s = t_e - MIN(64, (q_e << 1));
		} else if((t_e << 1) < q_e || (t_e + 64) < q_e) { // big leading query gap, cut down
			q_s = q_e - MIN(64, (t_e << 1));
		}
		
		/* align */
		if(t_e - t_s > 0 && q_e - q_s > 0) {
			band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
			if(q_e - q_s <= band || t_e - t_s <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
				NWstat = NW(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, matrices);
			} else {
				NWstat = NW_band(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, Frag_align, band, matrices);
			}
			/* trim leading gaps */
			bias = 0;
			if(t_s == 0) {
				while(bias < NWstat.len && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
					if(Frag_align->t[bias] == 5) {
						--NWstat.gaps;
					}
					++bias;
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
		}
	}
	
	/* piece seeds together */
	stop = 1;
	while(stop) {
		/* MEM */
		q_s = points->qStart[start];
		end = points->qEnd[start] - q_s;
		memcpy(aligned->t + Stat.len, qseq + q_s, end);
		memset(aligned->s + Stat.len, '|', end);
		memcpy(aligned->q + Stat.len, qseq + q_s, end);
		Stat.len += end;
		end = points->qEnd[start];
		for(i = points->qStart[start]; i < end; ++i) {
			nuc = qseq[i];
			Stat.score += d[nuc][nuc];
		}
		
		/* join MEMs */
		if(points->next[start]) {
			/* get positions between seed-extends */
			q_s = points->qEnd[start];
			t_s = points->tEnd[start] - 1;
			start = points->next[start];
			
			/* check if next MEM is a semi match */
			if(points->qStart[start] < q_s) {
				points->tStart[start] += (q_s - points->qStart[start]);
				points->qStart[start] = q_s;
			}
			t_e = points->tStart[start] - 1;
			
			if(t_e < t_s) {
				if(t_s <= points->tEnd[start]) {
					points->qStart[start] += (t_s - t_e);
					t_e = t_s;
					t_l = t_e - t_s;
				} else {
					/* circular joining */
					Frag_align->pos = t_len;
					t_l = t_len - t_s + t_e;
				}
			} else {
				t_l = t_e - t_s;
			}
			q_e = points->qStart[start];
			
			/* piece seed-extends together */
			if(abs(t_l - q_e + q_s) * U > q_len * M || t_l > q_len || q_e - q_s > (q_len >> 1)) {
				/* gap is too big to give a positive score */
				Stat.score = 0;
				Stat.len = 1;
				Stat.gaps = 0;
				aligned->s[0] = 0;
				points->len = 0;
				return Stat;
			}
			if((t_l > 0 || q_e - q_s > 0)) {
				band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
				if(q_e - q_s <= band || t_l <= band) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
					NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, matrices);
				} else {
					NWstat = NW_band(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align, band, matrices);
					//NWstat = NW(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, Frag_align);
				}
				memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
				memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
				memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
				Stat.score += NWstat.score;
				Stat.len += NWstat.len;
				Stat.gaps += NWstat.gaps;
			}
		} else {
			stop = 0;
		}
	}
	
	/* align trailing tail */
	/* Get intervals in query and template to align */
	q_s = points->qEnd[start];
	t_s = points->tEnd[start] - 1;
	q_e = q_len;
	t_e = t_len;
	if(((q_len - q_s) << 1) < (t_len - t_s) || (q_len - q_s + 64) < (t_len - t_s)) { // big trailing template gap, cut down
		t_e = t_s + MIN(64, ((q_len - q_s) << 1));
	} else if(((t_len - t_s) << 1) < (q_len - q_s) || (t_len - t_s + 64) < (q_len - q_s)) { // big leading query gap, cut down
		q_e = q_s + MIN(64, ((t_len - t_s) << 1));
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
		if(q_e - q_s <= band || t_e - t_s <= band) {//|| abs(t_e - t_s - q_e - q_s) >= 32) {
			NWstat = NW(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, matrices);
		} else {
			NWstat = NW_band(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align, band, matrices);
			//NWstat = NW(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, Frag_align);
		}
		/* trim trailing gaps */
		if(t_e == t_len) {
			bias = NWstat.len - 1;
			while(bias && (Frag_align->t[bias] == 5 || Frag_align->q[bias] == 5)) {
				if(Frag_align->t[bias] == 5) {
					--NWstat.gaps;
				}
				--bias;
			}
			++bias;
			/*
			if(bias != NWstat.len) {
				NWstat.score -= (W1 + (NWstat.len - bias) * U);
				NWstat.len = bias;
			}
			*/
		}
		
		memcpy(aligned->t + Stat.len, Frag_align->t, NWstat.len);
		memcpy(aligned->s + Stat.len, Frag_align->s, NWstat.len);
		memcpy(aligned->q + Stat.len, Frag_align->q, NWstat.len);
		Stat.score += NWstat.score;
		Stat.len += NWstat.len;
		Stat.gaps += NWstat.gaps;
	}
	aligned->s[Stat.len] = 0;
	points->len = 0;
	
	return Stat;
}

struct alnScore KMA_score(const int template_name, const unsigned char *qseq, int q_len, const struct compDNA *qseq_comp, struct alnPoints *points, struct NWmat *matrices) {
	
	int i, j, k, l, bias, prev, start, stop, t_len, value, end, band;
	int t_l, t_s, t_e, q_s, q_e, mem_count, score;
	unsigned mapQ;
	long unsigned key;
	unsigned char nuc;
	struct alnScore Stat, NWstat;
	struct hashMap_index *template_index;
	
	/* Extract indexes and template sequence */
	template_index = templates_index[template_name];
	t_len = template_lengths[template_name];
	
	/* find seeds */
	mem_count = 0;
	for(i = 1, j = 0; i <= qseq_comp->N[0]; ++i) {
		end = qseq_comp->N[i] - kmersize + 1;
		while(j < end) {
			value = hashMap_index_get(template_index, getKmer(qseq_comp->seq, j));
			
			if(value == 0) {
				++j;
			} else if(0 < value) {
				/* backseed for ambiguos seeds */
				prev = value - 2;
				for(k = j - 1; 0 <= k && 0 <= prev && qseq[k] == getNuc(template_index->seq, prev); --k) {
					--prev;
				}
				
				/* get start positions */
				points->qStart[mem_count] = k + 1;
				points->tStart[mem_count] = prev + 2;
				
				
				/* skip k-mer bases */
				value += (kmersize - 1);
				j += kmersize;
				
				/* extend */
				end += (kmersize - 1);
				while(j < end && value < t_len && qseq[j] == getNuc(template_index->seq, value)) {
					++j;
					++value;
				}
				
				/* get end positions */
				points->qEnd[mem_count] = j;
				points->tEnd[mem_count] = value + 1;
				
				/* calculate weight */
				points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
				++mem_count;
				
				/* realloc seeding points */
				if(mem_count == points->size) {
					seedPoint_realloc(points, points->size << 1);
				}
			} else {
				/* get position in hashmap */
				key = getKmer(qseq_comp->seq, j);
				stop = hashMap_index_getDubPos(template_index, key, value);
				
				/* get all mems */
				bias = j;
				while(0 <= stop) {
					/* get mem info */
					l = j;
					/* backseed for overlapping seeds */
					value = abs(template_index->index[stop]);
					prev = value - 2;
					for(k = l - 1; 0 <= k && 0 <= prev && qseq[k] == getNuc(template_index->seq, prev); --k) {
						--prev;
					}
					
					/* get start positions */
					points->qStart[mem_count] = k + 1;
					points->tStart[mem_count] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					l += kmersize;
					
					/* extend */
					while(l < end && value < t_len && qseq[l] == getNuc(template_index->seq, value)) {
						++l;
						++value;
					}
					
					/* get end positions */
					points->qEnd[mem_count] = l;
					points->tEnd[mem_count] = value + 1;
					
					/* calculate weight */
					points->weight[mem_count] = (points->qEnd[mem_count] - points->qStart[mem_count]);
					++mem_count;
					
					/* realloc seeding points */
					if(mem_count == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					if(bias < l) {
						bias = l;
					}
					
					stop = hashMap_index_getNextDubPos(template_index, key, 0, t_len, stop);
				}
				j = bias + 1;
			}
		}
		j = qseq_comp->N[i] + 1;
	}
	
	if(mem_count) {
		if(points->qEnd[mem_count - 1] != q_len && points->tEnd[mem_count - 1] != t_len) {
			--mem_count;
			score = (MIN(q_len - points->qEnd[mem_count], t_len - points->tEnd[mem_count])) * (MM + M);
			/* W1 is deprecated, new pen for non-pairing is needed */
			points->weight[mem_count] += MAX(W1, score);
			++mem_count;
		}
		points->len = mem_count;
	} else {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		return Stat;
	}
	
	/* get best seed chain, returns best starting point */
	start = chainSeedsPtr(points, q_len, t_len, &mapQ);
	score = points->score[start];
	
	//if(score < (q_len >> 5)) {
	if(mapQ < mq || score < kmersize || (score << 1) < scoreT * (points->qEnd[points->len - 1] - points->qStart[0])) {
		Stat.score = 0;
		Stat.len = 1;
		Stat.gaps = 0;
		Stat.pos = 0;
		return Stat;
	}
	
	/* initialize */
	Stat.len = 0;
	Stat.score = 0;
	Stat.gaps = 0;
	value = points->tStart[start] - 1;
	Stat.pos = value;
	i = points->qStart[start];
	
	/* align leading tail */
	if(i != 0) {
		/* get boundaries */
		t_s = 0;
		t_e = value;
		q_s = 0;
		q_e = i;
		if(value > (i + 64) || value > (i << 1)) { // big leading template gap, cut down
			t_s = value - ((i < 64) ? (i << 1) : (i + 64));
		} else if (i > (value + 64) || i > (value << 1)) { // big leading query gap, cut down
			q_s = i - ((value < 64) ? (value << 1) : (value + 64));
		}
		
		/* align */
		if(t_e - t_s > 0 && q_e - q_s > 0) {
			band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
			if(q_e - q_s <= (band << 1) || t_e - t_s <= (band << 1)) {// || abs(t_e - t_s - q_e - q_s) >= 32) {
				NWstat = NW_score(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, matrices, t_len);
			} else {
				NWstat = NW_band_score(template_index->seq, qseq, -1 - (t_s == 0), t_s, t_e, q_s, q_e, band, matrices, t_len);
			}
			Stat.pos -= (NWstat.len - NWstat.gaps);
			Stat.score = NWstat.score;
			Stat.len = NWstat.len;
			Stat.gaps = NWstat.gaps;
		}
	}
	
	/* piece seeds together */
	stop = 1;
	while(stop) {
		/* MEM */
		q_s = points->qStart[start];
		end = points->qEnd[start] - q_s;
		Stat.len += end;
		end = points->qEnd[start];
		for(i = points->qStart[start]; i < end; ++i) {
			nuc = qseq[i];
			Stat.score += d[nuc][nuc];
		}
		
		/* join MEMs */
		if(points->next[start]) {
			/* get positions between seed-extends */
			q_s = points->qEnd[start];
			t_s = points->tEnd[start] - 1;
			start = points->next[start];
			
			/* check if next MEM is a semi match, or a circular joining */
			if(points->qStart[start] < q_s) {
				points->tStart[start] += (q_s - points->qStart[start]);
				points->qStart[start] = q_s;
			}
			t_e = points->tStart[start] - 1;
			
			if(t_e < t_s) {
				if(t_s <= points->tEnd[start]) {
					points->qStart[start] += (t_s - t_e);
					t_e = t_s;
					t_l = t_e - t_s;
				} else {
					/* circular joining */
					t_l = t_len - t_s + t_e;
				}
			} else {
				t_l = t_e - t_s;
			}
			q_e = points->qStart[start];
			
			/* piece seed-extends together */
			if(abs(t_l - q_e + q_s) * U > q_len * M || t_l > q_len || q_e - q_s > (q_len >> 1)) {
				/* gap is too big to give a positive score */
				Stat.score = 0;
				Stat.len = 1;
				Stat.gaps = 0;
				return Stat;
			}
			if((t_l > 0 || q_e - q_s > 0)) {
				band = 4 * abs(t_l - q_e + q_s) + 64;
				if(q_e - q_s <= (band << 1) || t_l <= (band << 1)) {
					NWstat = NW_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, matrices, t_len);
				} else {
					NWstat = NW_band_score(template_index->seq, qseq, 0, t_s, t_e, q_s, q_e, band, matrices, t_len);
				}
				Stat.score += NWstat.score;
				Stat.len += NWstat.len;
				Stat.gaps += NWstat.gaps;
			}
		} else {
			stop = 0;
		}
	}
	
	/* align trailing tail */
	/* Get intervals in query and template to align */
	q_s = points->qEnd[start];
	t_s = points->tEnd[start] - 1;
	q_e = q_len;
	t_e = t_len;
	if((t_len - t_s) > (q_len - q_s + 64) || (t_len - t_s) > ((q_len - q_s) << 1)) { // big trailing template gap, cut down
		t_e = t_s + ((((q_len - q_s)) < 64) ? ((q_len - q_s) << 1) : (q_len - q_s + 64));
	} else if ((q_len - q_s) > (t_len - t_s + 64) || (q_len - q_s) > ((t_len - t_s) << 1)) { // big leading query gap, cut down
		q_e = q_s + (((t_len - t_s) < 64) ? ((t_len - t_s) << 1) : ((t_len - t_s) + 64));
	}
	
	/* align trailing gap */
	if(t_e - t_s > 0 && q_e - q_s > 0) {
		band = 4 * abs(t_e - t_s - q_e + q_s) + 64;
		if(q_e - q_s <= (band << 1) || t_e - t_s <= (band << 1)) {//|| abs(t_e - t_s - q_e - q_s) >= 32) {
			NWstat = NW_score(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, matrices, t_len);
		} else {
			NWstat = NW_band_score(template_index->seq, qseq, 1 + (t_e == t_len), t_s, t_e, q_s, q_e, band, matrices, t_len);
		}
		Stat.score += NWstat.score;
		Stat.len += NWstat.len;
		Stat.gaps += NWstat.gaps;
	}
	
	return Stat;
}

int preseed(struct hashMap_index *template_index, unsigned char *qseq, int q_len) {
	
	int i;
	
	if(exhaustive) {
		return 0;
	}
	
	for(i = 0; i < q_len; i += kmersize) {
		if(hashMap_index_get_bound(template_index, makeKmer(qseq, i, kmersize), 0, template_index->len) != 0) {
			return 0;
		}
	}
	
	return i;
}

void intcpy(int *dest, int *src, int size) {
	
	while(size--) {
		*dest++ = *src++;
	}
}

int anker_rc(const int template_name, unsigned char *qseq, int q_len, struct alnPoints *points) {
	
	int i, j, k, rc, end, stop, score, score_r, value, t_len, prev, bias;
	int bestScore, mem_count, totMems;
	long unsigned key;
	struct hashMap_index *template_index;
	
	t_len = template_lengths[template_name];
	template_index = templates_index[template_name];
	key = 0;
	
	/* find seeds */
	bestScore = 0;
	score = 0;
	score_r = 0;
	mem_count = 0;
	totMems = 0;
	points->len = 0;
	for(rc = 0; rc < 2; ++rc) {
		if(rc) {
			strrc(qseq, q_len);
			score = score_r;
			points->len = mem_count;
		}
		score_r = 0;
		mem_count = 0;
		i = preseed(template_index, qseq, q_len);
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
			
			while(i < end) {
				key = ((key << 2) | qseq[i]) & mask;
				value = hashMap_index_get_bound(template_index, key, 0, t_len);
				if(value == 0) {
					++i;
				} else if(0 < value) {
					i -= (kmersize - 1);
					
					/* backseed for ambiguos seeds */
					prev = value - 2;
					for(j = i - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
						--prev;
						++score_r;
					}
					
					/* get start positions */
					points->qStart[totMems] = j + 1;
					points->tStart[totMems] = prev + 2;
					
					/* skip k-mer bases */
					value += (kmersize - 1);
					i += kmersize;
					score_r += kmersize;
					
					/* extend */
					while(i < end && value < t_len && qseq[i] == getNuc(template_index->seq, value)) {
						++i;
						++value;
						++score_r;
					}
					
					/* get end positions */
					points->qEnd[totMems] = i;
					points->tEnd[totMems] = value + 1;
					
					/* calculate weight */
					points->weight[totMems] = (points->tEnd[totMems] - points->tStart[totMems]);
					++mem_count;
					++totMems;
					
					/* realloc seeding points */
					if(totMems == points->size) {
						seedPoint_realloc(points, points->size << 1);
					}
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				} else {
					/* move counter back */
					i -= (kmersize - 1);
					score_r += kmersize;
					
					/* get position in hashmap */
					stop = hashMap_index_getDubPos(template_index, key, value);
					
					/* get all mems */
					bias = i;
					while(0 <= stop) {
						/* get mem info */
						k = i;
						/* backseed for overlapping seeds */
						value = abs(template_index->index[stop]);
						prev = value - 2;
						for(j = k - 1; 0 <= j && 0 <= prev && qseq[j] == getNuc(template_index->seq, prev); --j) {
							--prev;
						}
						
						/* get start positions */
						points->qStart[totMems] = j + 1;
						points->tStart[totMems] = prev + 2;
						
						/* skip k-mer bases */
						value += (kmersize - 1);
						k += kmersize;
						
						/* extend */
						while(k < end && value < t_len && qseq[k] == getNuc(template_index->seq, value)) {
							++k;
							++value;
						}
						
						/* get end positions */
						points->qEnd[totMems] = k;
						points->tEnd[totMems] = value + 1;
						
						/* calculate weight */
						points->weight[totMems] = (points->qEnd[totMems] - points->qStart[totMems]);
						++mem_count;
						++totMems;
						
						/* realloc seeding points */
						if(totMems == points->size) {
							seedPoint_realloc(points, points->size << 1);
						}
						
						if(bias < k) {
							bias = k;
						}
						
						stop = hashMap_index_getNextDubPos(template_index, key, 0, t_len, stop);
					}
					i = bias + 1;
					
					/* update position */
					if(i < end - kmersize) {
						key = makeKmer(qseq, i, kmersize - 1);
						i += (kmersize - 1);
					} else {
						i = end + 1;
					}
				}
			}
			i = end + 1;
		}
		if(bestScore < score_r) {
			bestScore = score_r;
		}
	}
	
	if(bestScore * kmersize < (q_len - kmersize - bestScore)) {
		bestScore = 0;
		points->len = 0;
	} else if(bestScore == score) {
		strrc(qseq, q_len);
		if(Mt1Score) {
			Mt1Score[template_name] += bestScore;
		}
	} else {
		/* move mems down */
		if(points->len) {
			intcpy(points->tStart, points->tStart + points->len, mem_count);
			intcpy(points->tEnd, points->tEnd + points->len, mem_count);
			intcpy(points->qStart, points->qStart + points->len, mem_count);
			intcpy(points->qEnd, points->qEnd + points->len, mem_count);
			intcpy(points->weight, points->weight + points->len, mem_count);
		}
		points->len = mem_count;
		if(Mt1Score) {
			Mt1Score[template_name] += bestScore;
		}
	}
	
	return bestScore;
}

int significantNuc(int X, int Y) {
	return (p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue && Y < X);
}

int significantAnd90Nuc(int X, int Y) {
	return (p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue && (9 * (X + Y) <= 10 * X));
}

int significantAndSupport(int X, int Y) {
	return (p_chisqr(pow(X - Y, 2) / (X + Y)) <= evalue && (support * (X + Y) <= X));
}

unsigned char baseCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, struct assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = '-';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore) == 0) {
			if(bestNuc == '-' && tNuc != '-') {
				bestNuc = 'n';
			} else {
				bestNuc = tolower(bestNuc);
			}
		}
	}
	
	return bestNuc;
}

unsigned char refCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, struct assembly *calls) {
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = 'n';
	} else if(bestNuc == '-' && tNuc != '-') {
		bestNuc = 'n';
	} else if(significantBase(bestScore, depthUpdate - bestScore) == 0) {
		bestNuc = tolower(bestNuc);
	}
	
	return bestNuc;
}

unsigned char nanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, struct assembly *calls) {
	
	int j, bestBaseScore;
	const char bases[] = "ACGTN-";
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = '-';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore) == 0) {
			if(bestNuc == '-' && tNuc != '-') {
				bestBaseScore = 0;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < calls->counts[j]) {
						bestBaseScore = calls->counts[j];
						bestNuc = j;
					}
				}
				if(bestBaseScore == 0) {
					bestNuc = '-';
				} else {
					bestNuc = tolower(bases[bestNuc]);
				}
			} else {
				bestNuc = tolower(bestNuc);
			}
		}
	}
	
	return bestNuc;
}

unsigned char refNanoCaller(unsigned char bestNuc, unsigned char tNuc, int bestScore, int depthUpdate, struct assembly *calls) {
	
	int j, bestBaseScore;
	const char bases[] = "ACGTN-";
	
	/* determine base at current position */
	if(depthUpdate == 0) {
		bestNuc = 'n';
	} else {
		/* Use MC Neymars test to test significance of the base call */
		if(significantBase(bestScore, depthUpdate - bestScore) == 0) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < calls->counts[j]) {
						bestBaseScore = calls->counts[j];
						bestNuc = j;
					}
				}
				if(bestBaseScore == 0) {
					bestNuc = 'n';
				} else {
					bestNuc = tolower(bases[bestNuc]);
				}
			} else {
				bestNuc = tolower(bestNuc);
			}
		} else if(bestNuc == '-') {
			bestNuc = 'n';
		}
	}
	
	return bestNuc;
}

void * assemble_KMA_threaded(void *arg) {
	
	static volatile int excludeIn = 0, excludeOut = 0, excludeMatrix = 0, mainTemplate = -2, thread_wait = 0;
	struct assemble_thread *thread = arg;
	int i, j, t_len, aln_len, start, end, bias, myBias, gaps, pos, asm_len;
	int read_score, depthUpdate, bestBaseScore, bestScore, template, spin;
	int nextTemplate, file_i, file_count, delta, stats[4], buffer[7];
	unsigned coverScore;
	long unsigned depth, depthVar;
	const char bases[] = "ACGTN-";
	double score;
	unsigned char bestNuc;
	FILE **files, *file;
	struct alnScore alnStat;
	struct assembly *assembly;
	struct FileBuff *frag_out;
	struct assem *aligned_assem;
	struct aln *aligned, *gap_align;
	struct qseqs *qseq, *header;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	spin = thread->spin;
	
	if(template != -2) {
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(&excludeMatrix);
			mainTemplate = template;
			unlock(&excludeMatrix);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(&excludeMatrix);
		t_len = template_lengths[template];
		matrix->len = t_len;
		if(matrix->size < (t_len << 1)) {
			matrix->size = (t_len << 1);
			free(matrix->assmb);
			matrix->assmb = malloc(matrix->size * sizeof(struct assembly));
			if(!matrix->assmb) {
				ERROR();
			}
		}
		
		/* cpy template seq */
		assembly = matrix->assmb;
		for(i = 0, j = 1; i < t_len; ++i, ++j) {
			assembly[i].counts[0] = 0;
			assembly[i].counts[1] = 0;
			assembly[i].counts[2] = 0;
			assembly[i].counts[3] = 0;
			assembly[i].counts[4] = 0;
			assembly[i].counts[5] = 0;
			assembly[i].next = j;
		}
		/* circularize */
		assembly[t_len - 1].next = 0;
		
		/* start threads */
		aligned_assem->score = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		unlock(&excludeMatrix);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(&excludeMatrix);
		template = mainTemplate;
		if(template != -1) {
			t_len = template_lengths[template];
		}
		unlock(&excludeMatrix);
		if(template == -1) {
			return NULL;
		}
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			lockTime(&excludeIn, spin);
			file = files[file_i];
			if(file != 0) {
				fread(buffer, sizeof(int), 7, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = malloc(qseq->size);
						if(!qseq->seq) {
							ERROR();
						}
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = malloc(header->size);
						if(!header->seq) {
							ERROR();
						}
					}
					fread(qseq->seq, 1, qseq->len, file);
					fread(header->seq, 1, header->len, file);
					unlock(&excludeIn);
					
					if(delta < qseq->len) {
						delta = qseq->len << 1;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = malloc((delta + 1) << 1);
						aligned->s = malloc((delta + 1) << 1);
						aligned->q = malloc((delta + 1) << 1);
						gap_align->t = malloc((delta + 1) << 1);
						gap_align->s = malloc((delta + 1) << 1);
						gap_align->q = malloc((delta + 1) << 1);
						if(!aligned->t || !aligned->s || !aligned->q || !gap_align->t || !gap_align->s || !gap_align->q) {
							ERROR();
						}
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template, qseq->seq, qseq->len, points)) {
						/* Start with alignment */
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						alnStat = KMA(template, qseq->seq, qseq->len, aligned, gap_align, stats[2], MIN(t_len, stats[3]), points, NWmatrices);
						
						/* get read score */
						aln_len = alnStat.len;
						start = alnStat.pos;
						end = start + aln_len - alnStat.gaps;
						
						/* Get normed score */
						read_score = alnStat.score;
						if(aln_len >= kmersize) {
							score = 1.0 * read_score / aln_len;
						} else {
							score = 0;
						}
						
						if(read_score > kmersize && score >= scoreT) {
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							//lock(&excludeMatrix);
							lockTime(&excludeMatrix, 10)
							aligned_assem->score += read_score;
							
							/* diff */
							i = 0;
							pos = start;
							assembly = matrix->assmb;
							while(i < aln_len) {
								if(aligned->t[i] == 5) { // Template gap, insertion
									if(t_len <= pos) {
										assembly[pos].counts[aligned->q[i]]++;
										++i;
										pos = assembly[pos].next;
									} else {
										/* get estimate for non insertions */
										myBias = 0;
										for(j = 0; j < 6; ++j) {
											myBias += assembly[pos].counts[j];
										}
										if(myBias > 0) {
											--myBias;
										}
										/* find position of insertion */
										gaps = pos;
										if(pos != 0) {
											--pos;
										} else {
											pos = t_len - 1;
										}
										while(assembly[pos].next != gaps) {
											pos = assembly[pos].next;
										}
										while(i < aln_len && aligned->t[i] == 5) {
											assembly[pos].next = matrix->len++;
											if(matrix->len == matrix->size) {
												matrix->size <<= 1;
												matrix->assmb = realloc(assembly, matrix->size * sizeof(struct assembly));
												if(!matrix->assmb) {
													matrix->size >>= 1;
													matrix->size += 1024;
													matrix->assmb = realloc(assembly, matrix->size * sizeof(struct assembly));
													if(!matrix->assmb) {
														ERROR();
													}
												}
												assembly = matrix->assmb;
											}
											pos = assembly[pos].next;
											assembly[pos].next = gaps;
											assembly[pos].counts[0] = 0;
											assembly[pos].counts[1] = 0;
											assembly[pos].counts[2] = 0;
											assembly[pos].counts[3] = 0;
											assembly[pos].counts[4] = 0;
											assembly[pos].counts[5] = myBias;
											assembly[pos].counts[aligned->q[i]]++;
											
											++i;
										}
										pos = assembly[pos].next;
									}
								} else if(t_len <= pos) { // Old template gap, not present in this read
									assembly[pos].counts[5]++;
									pos = assembly[pos].next;
								} else {
									assembly[pos].counts[aligned->q[i]]++;
									++i;
									pos = assembly[pos].next;
								}
								/* debug info */
								/*else if(aligned->t[i] == getNuc(templates_index[template]->seq, pos)) { // (Match, mismatch) and (Query gap, deletion)
									assembly[pos].counts[aligned->q[i]]++;
									++i;
									pos = assembly[pos].next;
								} else {
									fprintf(stderr, "KMA is out of sync\n");
									pos = start;
									
									for(i = 0; i < aln_len; ++i) {
										if(aligned->t[i] != 5) {
											fprintf(stderr, "%c\t%c\t%c\t%c\n", bases[getNuc(templates_index[template]->seq, pos)], bases[aligned->t[i]], aligned->s[i], bases[aligned->q[i]]);
											++pos;
										} else {
											fprintf(stderr, "%c\t%c\t%c\t%c\n", '-', bases[aligned->t[i]], aligned->s[i], bases[aligned->q[i]]);
										}
									}
									
									for(i = 0; i < aln_len; ++i) {
										fprintf(stderr, "%c", bases[aligned->t[i]]);
									}
									fprintf(stderr, "\n");
									for(i = 0; i < aln_len; ++i) {
										fprintf(stderr, "%c", aligned->s[i]);
									}
									fprintf(stderr, "\n");
									for(i = 0; i < aln_len; ++i) {
										fprintf(stderr, "%c", bases[aligned->q[i]]);
									}
									fprintf(stderr, "\n");
									exit(1);
									++i;
								}
								*/
							}
							
							unlock(&excludeMatrix);
							
							/* Convert fragment */
							for(i = 0; i < qseq->len; ++i) {
								 qseq->seq[i] = bases[qseq->seq[i]];
							}
							qseq->seq[qseq->len] = 0;
							
							/* Save fragment */
							//lock(&excludeOut);
							lockTime(&excludeOut, 10);
							updateFrags(frag_out, qseq, header, template_names[template], stats);
							unlock(&excludeOut);
							//fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
						}
					}
				} else if(nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						errno |= kmaPclose(file);
					}
					files[file_i] = 0;
					unlock(&excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					fseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(&excludeIn);
				} else {
					/* Move pointer back */
					fseek(file, (-7) * sizeof(int), SEEK_CUR);
					unlock(&excludeIn);
					++file_i;
				}
			} else {
				unlock(&excludeIn);
				++file_i;
			}
		}
		lock(&excludeIn);
		--thread_wait;
		unlock(&excludeIn);
	} while(thread->num != 0);
	
	wait_atomic(thread_wait);
	
	if(aligned_assem->score == 0) {
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		
		return NULL;
	}
	
	/* diff */
	/* pre on dense */
	/* Pepare and make alignment on consensus */
	asm_len = matrix->len;
	assembly = matrix->assmb;
	if(aligned_assem->size <= asm_len) {
		aligned_assem->size = (asm_len + 1) << 1;
		free(aligned_assem->t);
		free(aligned_assem->s);
		free(aligned_assem->q);
		aligned_assem->t = malloc(aligned_assem->size);
		aligned_assem->s = malloc(aligned_assem->size);
		aligned_assem->q = malloc(aligned_assem->size);
		if(!aligned_assem->t || !aligned_assem->s || !aligned_assem->q) {
			ERROR();
		}
	}
	
	/* Call nucleotides for the consensus */
	/* diff */
	i = 0;
	pos = 0;
	depth = 0;
	depthVar = 0;
	aln_len = 0;
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
		for(j = 0; j < 6; ++j) {
			if(bestScore < assembly[pos].counts[j]) {
				bestScore = assembly[pos].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[pos].counts[j];
		}
		bestNuc = bases[bestNuc];
		
		/* check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < assembly[pos].counts[j]) {
						bestBaseScore = assembly[pos].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[pos].counts[5];
		}
		
		/* determine base at current position */
		bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, &assembly[pos]);
		aligned_assem->q[i] = bestNuc;
		
		if(bestNuc != '-') {
			depth += depthUpdate;
			depthVar += (depthUpdate * depthUpdate);
			++aln_len;
		}
		
		++i;
		pos = assembly[pos].next;
	}
	
	/* Trim alignment on consensus */
	coverScore = 0;
	bias = 0;
	for(i = 0; i < asm_len; ++i) {
		if(aligned_assem->t[i] == '-' && aligned_assem->q[i] == '-') {
			++bias;
		} else {
			aligned_assem->t[i - bias] = aligned_assem->t[i];
			aligned_assem->q[i - bias] = aligned_assem->q[i];
			if(tolower(aligned_assem->t[i]) == tolower(aligned_assem->q[i])) {
				aligned_assem->s[i - bias] = '|';
				++coverScore;
			} else {
				aligned_assem->s[i - bias] = '_';
			}
		}
	}
	asm_len -= bias;
	aligned_assem->t[asm_len] = 0;
	aligned_assem->s[asm_len] = 0;
	aligned_assem->q[asm_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->depthVar = depthVar;
	aligned_assem->len = asm_len;
	aligned_assem->aln_len = aln_len;
	
	return NULL;
}

void * assemble_KMA_dense_threaded(void *arg) {
	
	static volatile int excludeIn = 0, excludeOut = 0, excludeMatrix = 0, mainTemplate = -2, thread_wait = 0;
	struct assemble_thread *thread = arg;
	int i, j, t_len, aln_len, start, end, file_i, file_count, template, spin;
	int pos, read_score, bestScore, depthUpdate, bestBaseScore, nextTemplate;
	int stats[4], buffer[7];
	unsigned coverScore, delta;
	long unsigned depth, depthVar;
	const char bases[] = "ACGTN-";
	double score;
	unsigned char bestNuc;
	FILE **files, *file;
	struct alnScore alnStat;
	struct assembly *assembly;
	struct FileBuff *frag_out;
	struct assem *aligned_assem;
	struct aln *aligned, *gap_align;
	struct qseqs *qseq, *header;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	
	/* get input */
	template = thread->template;
	file_count = thread->file_count;
	files = thread->files;
	frag_out = thread->frag_out;
	aligned_assem = thread->aligned_assem;
	aligned = thread->aligned;
	gap_align = thread->gap_align;
	qseq = thread->qseq;
	header = thread->header;
	matrix = thread->matrix;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	delta = qseq->size;
	spin = thread->spin;
	
	if(template != -2) {
		/* all assemblies done, 
		signal threads to return */
		if(template == -1) {
			lock(&excludeOut);
			mainTemplate = template;
			unlock(&excludeOut);
			return NULL;
		}
		
		/* Allocate assembly arrays */
		lock(&excludeOut);
		t_len = template_lengths[template];
		matrix->len = t_len;
		
		/* diff */
		if(aligned_assem->size <= t_len) {
			aligned_assem->size = t_len + 1;
			free(aligned_assem->t);
			free(aligned_assem->s);
			free(aligned_assem->q);
			aligned_assem->t = malloc(t_len + 1);
			aligned_assem->s = malloc(t_len + 1);
			aligned_assem->q = malloc(t_len + 1);
			if(!aligned_assem->t || !aligned_assem->s || !aligned_assem->q) {
				ERROR();
			}
		}
		if(matrix->size <= t_len) {
			matrix->size = t_len + 1;
			free(matrix->assmb);
			matrix->assmb = malloc(matrix->size * sizeof(struct assembly));
			if(!matrix->assmb) {
				ERROR();
			}
		}
		
		/* cpy template seq */
		assembly = matrix->assmb;
		for(i = 0, j = 1; i < t_len; ++i, ++j) {
			/* diff */
			aligned_assem->t[i] = getNuc(templates_index[template]->seq, i);
			assembly[i].counts[0] = 0;
			assembly[i].counts[1] = 0;
			assembly[i].counts[2] = 0;
			assembly[i].counts[3] = 0;
			assembly[i].counts[4] = 0;
			assembly[i].counts[5] = 0;
			assembly[i].next = j;
		}
		/* circularize */
		assembly[t_len - 1].next = 0;
		
		/* start threads */
		aligned_assem->score = 0;
		mainTemplate = template;
		thread_wait = thread_num;
		unlock(&excludeOut);
		template = -2;
	}
	
	do {
		while(template == mainTemplate) {
			usleep(100);
		}
		lock(&excludeOut);
		template = mainTemplate;
		if(template != -1) {
			t_len = template_lengths[template];
			assembly = matrix->assmb;
		}
		unlock(&excludeOut);
		if(template == -1) {
			return NULL;
		}
		
		/* load reads of this template */
		file_i = 0;
		while(file_i < file_count) {
			lockTime(&excludeIn, spin);
			file = files[file_i];
			if(file != 0) {
				fread(buffer, sizeof(int), 7, file);
				if((nextTemplate = buffer[0]) == template) {
					/* load frag */
					qseq->len = buffer[1];
					stats[0] = buffer[2];
					read_score = buffer[3];
					stats[2] = buffer[4];
					stats[3] = buffer[5];
					header->len = buffer[6];
					
					if(qseq->size < qseq->len) {
						free(qseq->seq);
						qseq->size = qseq->len << 1;
						qseq->seq = malloc(qseq->size);
						if(!qseq->seq) {
							ERROR();
						}
					}
					if(header->size < header->len) {
						header->size = header->len + 1;
						free(header->seq);
						header->seq = malloc(header->size);
						if(!header->seq) {
							ERROR();
						}
					}
					fread(qseq->seq, 1, qseq->len, file);
					fread(header->seq, 1, header->len, file);
					unlock(&excludeIn);
					
					if(delta < qseq->size) {
						delta = qseq->size;
						free(aligned->t);
						free(aligned->s);
						free(aligned->q);
						free(gap_align->t);
						free(gap_align->s);
						free(gap_align->q);
						aligned->t = malloc((delta + 1) << 1);
						aligned->s = malloc((delta + 1) << 1);
						aligned->q = malloc((delta + 1) << 1);
						gap_align->t = malloc((delta + 1) << 1);
						gap_align->s = malloc((delta + 1) << 1);
						gap_align->q = malloc((delta + 1) << 1);
						if(!aligned->t || !aligned->s || !aligned->q || !gap_align->t || !gap_align->s || !gap_align->q) {
							ERROR();
						}
					}
					
					/* Update assembly with read */
					if(read_score || anker_rc(template, qseq->seq, qseq->len, points)) {
						if(stats[3] <= stats[2]) {
							stats[2] = 0;
							stats[3] = t_len;
						}
						/* Start with alignment */
						alnStat = KMA(template, qseq->seq, qseq->len, aligned, gap_align, stats[2], MIN(t_len, stats[3]), points, NWmatrices);
						
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
							read_score = 0;
						}
						
						if(read_score > kmersize && score >= scoreT) {
							
							stats[1] = read_score;
							stats[2] = start;
							stats[3] = end;
							if(t_len < end) {
								stats[3] -= t_len;
							}
							/* Update backbone and counts */
							//lock(&excludeMatrix);
							lockTime(&excludeMatrix, 10)
							aligned_assem->score += read_score;
							
							/* diff */
							for(i = 0, pos = start; i < aln_len; ++i) {
								if(aligned->t[i] == aligned_assem->t[pos]) {
									assembly[pos].counts[aligned->q[i]]++;
									pos = assembly[pos].next;
								}
							}
							unlock(&excludeMatrix);
							
							/* Convert fragment */
							for(i = 0; i < qseq->len; ++i) {
								 qseq->seq[i] = bases[qseq->seq[i]];
							}
							qseq->seq[qseq->len] = 0;
							
							/* Save fragment */
							//lock(&excludeOut);
							lockTime(&excludeOut, 10);
							updateFrags(frag_out, qseq, header, template_names[template], stats);
							unlock(&excludeOut);
							//fprintf(frag_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", qseq->seq, stats[0], stats[1], stats[2], stats[3], template_names[template], header->seq);
						}
					}
				} else if (nextTemplate == -1) {
					if(template) {
						fclose(file);
					} else {
						errno |= kmaPclose(file);
					}
					files[file_i] = 0;
					unlock(&excludeIn);
					++file_i;
				} else if(nextTemplate < template) {
					/* Move pointer forward */
					fseek(file, buffer[1] + buffer[6], SEEK_CUR);
					unlock(&excludeIn);
				} else {
					/* Move pointer back */
					fseek(file, (-7) * sizeof(int), SEEK_CUR);
					unlock(&excludeIn);
					++file_i;
				}
			} else {
				unlock(&excludeIn);
				++file_i;
			}
		}
		
		lock(&excludeIn);
		--thread_wait;
		unlock(&excludeIn);
		
	} while(thread->num != 0);
	
	wait_atomic(thread_wait);
	
	if(aligned_assem->score == 0) {
		aligned_assem->cover = 0;
		aligned_assem->depth = 0;
		aligned_assem->depthVar = 0;
		aligned_assem->t[0] = 0;
		aligned_assem->s[0] = 0;
		aligned_assem->q[0] = 0;
		aligned_assem->len = 0;
		aligned_assem->aln_len = 0;
		
		return NULL;
	}
	
	/* Make consensus assembly by majority voting */
	/* diff */
	depth = 0;
	depthVar = 0;
	coverScore = 0;
	aln_len = 0;
	for(i = 0; i < t_len; ++i) {
		/* call template */
		aligned_assem->t[i] = bases[aligned_assem->t[i]];
		
		/* call query */
		bestNuc = 5;
		bestScore = 0;
		depthUpdate = 0;
		for(j = 0; j < 6; ++j) {
			if(bestScore < assembly[i].counts[j]) {
				bestScore = assembly[i].counts[j];
				bestNuc = j;
			}
			depthUpdate += assembly[i].counts[j];
		}
		bestNuc = bases[bestNuc];
		
		/* Check for minor base call */
		if((bestScore << 1) < depthUpdate) {
			if(bestNuc == '-') {
				bestBaseScore = 0;
				bestNuc = 4;
				for(j = 0; j < 5; ++j) {
					if(bestBaseScore < assembly[i].counts[j]) {
						bestBaseScore = assembly[i].counts[j];
						bestNuc = j;
					}
				}
				bestNuc = tolower(bases[bestNuc]);
			} else {
				bestNuc = tolower(bestNuc);
			}
			bestScore = depthUpdate - assembly[i].counts[5];
		}
		
		/* determine base at current position */
		bestNuc = baseCall(bestNuc, aligned_assem->t[i], bestScore, depthUpdate, &assembly[i]);
		aligned_assem->q[i] = bestNuc;
		
		if(bestNuc != '-') {
			depth += depthUpdate;
			depthVar += (depthUpdate * depthUpdate);
			++aln_len;
		}
		
		if(tolower(aligned_assem->q[i]) == tolower(aligned_assem->t[i])) {
			aligned_assem->s[i] = '|';
			++coverScore;
		} else {
			aligned_assem->s[i] = '_';
		}
	}
	aligned_assem->t[t_len] = 0;
	aligned_assem->s[t_len] = 0;
	aligned_assem->q[t_len] = 0;
	aligned_assem->cover = coverScore;
	aligned_assem->depth = depth;
	aligned_assem->depthVar = depthVar;
	aligned_assem->len = t_len;
	aligned_assem->aln_len = aln_len;
	
	return NULL;
}

void update_Scores(unsigned char *qseq, int q_len, int counter, int score, int *start, int *end, int *template, struct qseqs *header, FILE *frag_out_raw) {
	
	int i, buffer[4];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = score;
	buffer[3] = header->len;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 4, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { //Only one best match
		if(*template < 0) {
			template[0] = -template[0];
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter;
		while(i--) {
			alignment_scores[abs(template[i])] += score;
		}
	}
}

void update_Scores_pe(unsigned char *qseq, int q_len, unsigned char *qseq_r, int qr_len, int counter, int score, int *start, int *end, int *template, struct qseqs *header, struct qseqs *header_r, FILE *frag_out_raw) {
	
	int i, buffer[4];
	
	/* print frag */
	buffer[0] = q_len;
	buffer[1] = counter;
	buffer[2] = -score;
	buffer[3] = header->len;
	counter = abs(counter);
	sfwrite(buffer, sizeof(int), 4, frag_out_raw);
	sfwrite(qseq, 1, q_len, frag_out_raw);
	sfwrite(header->seq, 1, header->len, frag_out_raw);
	sfwrite(start, sizeof(int), counter, frag_out_raw);
	sfwrite(end, sizeof(int), counter, frag_out_raw);
	sfwrite(template, sizeof(int), counter, frag_out_raw);
	
	buffer[0] = qr_len;
	buffer[1] = header_r->len;
	sfwrite(buffer, sizeof(int), 2, frag_out_raw);
	sfwrite(qseq_r, 1, qr_len, frag_out_raw);
	sfwrite(header_r->seq, 1, header_r->len, frag_out_raw);
	
	/* update scores */
	if(counter == 1) { //Only one best match
		if(*template < 0) {
			template[0] = -template[0];
		}
		alignment_scores[*template] += score;
		uniq_alignment_scores[*template] += score;
	} else {
		i = counter;
		while(i--) {
			alignment_scores[abs(template[i])] += score;
		}
	}
}

void alnFragsSE(int *matched_templates, int rc_flag, struct compDNA *qseq_comp, struct compDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, int q_len, struct qseqs *header, int *bestTemplates, int *best_start_pos, int *best_end_pos, FILE *seq_in, FILE *index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, struct alnPoints *points, struct NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, bestHits, aln_len;
	int start, end;
	double score, bestScore;
	struct alnScore alnStat;
	
	/* reverse complement seq */
	if(rc_flag < 0) {
		if(qseq_r_comp->size < qseq_comp->size) {
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
	
	bestScore = 0;
	best_read_score = 0;
	bestHits = 0;
	
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		
		/* check if index DB is loaded */
		lock(excludeDB);
		if(template >= 0 && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(seq_in, index_in, template_lengths[-template], seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		/* align qseq */
		if(template < 0) {
			alnStat = KMA_score(-template, qseq_r, q_len, qseq_r_comp, points, NWmatrices);
		} else {
			alnStat = KMA_score(template, qseq, q_len, qseq_comp, points, NWmatrices);
		}
		
		/* get read score */
		aln_len = alnStat.len;
		start = alnStat.pos;
		end = start + aln_len - alnStat.gaps;
		if(template_lengths[abs(template)] < end) {
			end -= template_lengths[abs(template)];
		}
		
		/* penalty for non complete mapping */
		read_score = alnStat.score;
		/* full gene award */
		if((start == 0) && (end == template_lengths[abs(template)])) {
			read_score += abs(W1);
		}
		//read_score += (((start != 0) + (end != template_lengths[abs(template)])) * W1);
		
		/* Get normed score */
		if(aln_len > 0) {
			score = 1.0 * read_score / aln_len;
		} else {
			score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
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
				++bestHits;
			}
		}
	}
	if(best_read_score > kmersize) {
		lock(excludeOut);
		update_Scores(qseq, q_len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		unlock(excludeOut);
	}
}

void alnFragsUnionPE(int *matched_templates, struct compDNA *qseq_comp, struct compDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, struct qseqs *header, struct qseqs *header_r, int *bestTemplates, int *bestTemplates_r, int *best_start_pos, int *best_end_pos, FILE *seq_in, FILE *index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, struct alnPoints *points, struct NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, best_read_score_r;
	int compScore, bestHits, bestHits_r, aln_len, start, end, rc;
	double score;
	struct alnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	
	start = 0;
	end = 0;
	score = 0;
	best_read_score = 0;
	best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(seq_in, index_in, template_lengths[-template], seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		template = abs(template);
		
		/* align qseqs */
		alnStat = KMA_score(template, qseq, qseq_comp->seqlen, qseq_comp, points, NWmatrices);
		
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score <= read_score) {
			best_read_score = read_score;
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		alnStat = KMA_score(template, qseq_r, qseq_r_comp->seqlen, qseq_r_comp, points, NWmatrices);
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score_r <= read_score) {
			best_read_score_r = read_score;
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
		} else {
			bestTemplates_r[t_i] = 0;
			if(bestTemplates[t_i] != 0) {
				best_start_pos[t_i] = -1;
				best_end_pos[t_i] = -1;
			}
		}
		
		read_score += bestTemplates[t_i];
		if(best_start_pos[t_i] == 0 && best_end_pos[t_i] == template_lengths[abs(template)]) {
			read_score += abs(W1);
		}
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	if(best_read_score && best_read_score_r) {
		/* both matched */
		if(compScore && compScore == (best_read_score + best_read_score_r)) {
			/* proper pair */
			bestHits = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(compScore == (bestTemplates[t_i] + bestTemplates_r[t_i])) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header_r, header, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
				}
				lock(excludeOut);
				update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header, header_r, frag_out_raw);
				unlock(excludeOut);
			}
		} else {
			/* unmaided pair */
			bestHits = 0;
			bestHits_r = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] == best_read_score) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				} else if(bestTemplates_r[t_i] == best_read_score_r) {
					bestTemplates_r[bestHits_r] = matched_templates[t_i];
					best_start_pos[bestHits_r] = best_start_pos[t_i];
					best_end_pos[bestHits_r] = best_end_pos[t_i];
					++bestHits_r;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
			}
			lock(excludeOut);
			update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
			update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, frag_out_raw);
			unlock(excludeOut);
		}
	} else if(best_read_score) {
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates[t_i] == best_read_score) {
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = best_start_pos[t_i];
				best_end_pos[bestHits] = best_end_pos[t_i];
				++bestHits;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
		} else if(!rc) {
			strrc(qseq, qseq_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		unlock(excludeOut);
	} else if(best_read_score_r) {
		bestHits_r = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates_r[t_i] == best_read_score_r) {
				bestTemplates_r[bestHits_r] = matched_templates[t_i];
				best_start_pos[bestHits_r] = best_start_pos[t_i];
				best_end_pos[bestHits_r] = best_end_pos[t_i];
				++bestHits_r;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates_r < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates_r[t_i] = -bestTemplates_r[t_i];
			}
		} else if(!rc) {
			strrc(qseq_r, qseq_r_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, frag_out_raw);
		unlock(excludeOut);
	}
}

void alnFragsPenaltyPE(int *matched_templates, struct compDNA *qseq_comp, struct compDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, struct qseqs *header, struct qseqs *header_r, int *bestTemplates, int *bestTemplates_r, int *best_start_pos, int *best_end_pos, FILE *seq_in, FILE *index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, struct alnPoints *points, struct NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, best_read_score_r;
	int compScore, bestHits, bestHits_r, aln_len, start, end, rc;
	double score;
	struct alnScore alnStat;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	
	start = 0;
	end = 0;
	score = 0;
	best_read_score = 0;
	best_read_score_r = 0;
	compScore = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(seq_in, index_in, template_lengths[-template], seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		template = abs(template);
		
		/* align qseqs */
		alnStat = KMA_score(template, qseq, qseq_comp->seqlen, qseq_comp, points, NWmatrices);
		
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score <= read_score) {
			best_read_score = read_score;
			bestTemplates[t_i] = read_score;
			best_start_pos[t_i] = start;
			best_end_pos[t_i] = end;
		} else {
			bestTemplates[t_i] = 0;
			best_start_pos[t_i] = -1;
			best_end_pos[t_i] = -1;
		}
		
		alnStat = KMA_score(template, qseq_r, qseq_r_comp->seqlen, qseq_r_comp, points, NWmatrices);
		/* get read score */
		if(0 < alnStat.score) {
			aln_len = alnStat.len;
			
			start = alnStat.pos;
			end = alnStat.pos + alnStat.len - alnStat.gaps;
			
			read_score = alnStat.score;
			if(start == 0 && end == template_lengths[abs(template)]) {
				read_score += abs(W1);
			}
			score = 1.0 * read_score / aln_len;
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT && best_read_score_r <= read_score) {
			best_read_score_r = read_score;
			bestTemplates_r[t_i] = read_score;
			if(bestTemplates[t_i]) {
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(start < best_start_pos[t_i]) {
					best_start_pos[t_i] = start;
				} else {
					best_end_pos[t_i] = end;
				}
			} else {
				best_start_pos[t_i] = start;
				best_end_pos[t_i] = end;
			}
		} else {
			bestTemplates_r[t_i] = 0;
			if(bestTemplates[t_i] != 0) {
				best_start_pos[t_i] = -1;
				best_end_pos[t_i] = -1;
			}
		}
		
		read_score += (bestTemplates[t_i] + PE);
		if(compScore < read_score) {
			compScore = read_score;
		}
	}
	
	if(best_read_score && best_read_score_r) {
		/* both matched */
		if(compScore && (best_read_score + best_read_score_r) <= compScore) {
			/* proper pair */
			compScore -= PE;
			bestHits = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(compScore == (bestTemplates[t_i] + bestTemplates_r[t_i])) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
				lock(excludeOut);
				update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header_r, header, frag_out_raw);
				unlock(excludeOut);
			} else {
				if(!rc) {
					strrc(qseq, qseq_comp->seqlen);
					strrc(qseq_r, qseq_r_comp->seqlen);
				}
				lock(excludeOut);
				update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, compScore, best_start_pos, best_end_pos, bestTemplates, header, header_r, frag_out_raw);
				unlock(excludeOut);
			}
		} else {
			/* unmaided pair */
			bestHits = 0;
			bestHits_r = 0;
			for(t_i = 1; t_i <= *matched_templates; ++t_i) {
				if(bestTemplates[t_i] == best_read_score) {
					bestTemplates[bestHits] = matched_templates[t_i];
					best_start_pos[bestHits] = best_start_pos[t_i];
					best_end_pos[bestHits] = best_end_pos[t_i];
					++bestHits;
				} else if(bestTemplates_r[t_i] == best_read_score_r) {
					bestTemplates_r[bestHits_r] = matched_templates[t_i];
					best_start_pos[bestHits_r] = best_start_pos[t_i];
					best_end_pos[bestHits_r] = best_end_pos[t_i];
					++bestHits_r;
				}
			}
			/* check direction of qseqs */
			if(*bestTemplates < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates[t_i] = -bestTemplates[t_i];
				}
			} else if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
			}
			if(*bestTemplates_r < 0) {
				for(t_i = 0; t_i < bestHits; ++t_i) {
					bestTemplates_r[t_i] = -bestTemplates_r[t_i];
				}
			} else if(!rc) {
				strrc(qseq_r, qseq_r_comp->seqlen);
			}
			lock(excludeOut);
			update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
			update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, frag_out_raw);
			unlock(excludeOut);
		}
	} else if(best_read_score) {
		bestHits = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates[t_i] == best_read_score) {
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = best_start_pos[t_i];
				best_end_pos[bestHits] = best_end_pos[t_i];
				++bestHits;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
		} else if(!rc) {
			strrc(qseq, qseq_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
		unlock(excludeOut);
	} else if(best_read_score_r) {
		bestHits_r = 0;
		for(t_i = 1; t_i <= *matched_templates; ++t_i) {
			if(bestTemplates_r[t_i] == best_read_score_r) {
				bestTemplates_r[bestHits_r] = matched_templates[t_i];
				best_start_pos[bestHits_r] = best_start_pos[t_i];
				best_end_pos[bestHits_r] = best_end_pos[t_i];
				++bestHits_r;
			}
		}
		/* check direction of qseqs */
		if(*bestTemplates_r < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates_r[t_i] = -bestTemplates_r[t_i];
			}
		} else if(!rc) {
			strrc(qseq_r, qseq_r_comp->seqlen);
		}
		lock(excludeOut);
		update_Scores(qseq_r, qseq_r_comp->seqlen, bestHits_r, best_read_score_r, best_start_pos, best_end_pos, bestTemplates_r, header_r, frag_out_raw);
		unlock(excludeOut);
	}
}

void alnFragsForcePE(int *matched_templates, struct compDNA *qseq_comp, struct compDNA *qseq_r_comp, unsigned char *qseq, unsigned char *qseq_r, struct qseqs *header, struct qseqs *header_r, int *bestTemplates, int *bestTemplates_r, int *best_start_pos, int *best_end_pos, FILE *seq_in, FILE *index_in, long *seq_indexes, long *index_indexes, FILE *frag_out_raw, struct alnPoints *points, struct NWmat *NWmatrices, volatile int *excludeOut, volatile int *excludeDB) {
	
	int t_i, template, read_score, best_read_score, bestHits, aln_len;
	int start, end, rc;
	double score, bestScore;
	struct alnScore alnStat, alnStat_r;
	
	/* unpack qseqs */
	unCompDNA(qseq_comp, qseq);
	qseq_comp->N[0]++;
	qseq_comp->N[qseq_comp->N[0]] = qseq_comp->seqlen;
	unCompDNA(qseq_r_comp, qseq_r);
	qseq_r_comp->N[0]++;
	qseq_r_comp->N[qseq_r_comp->N[0]] = qseq_r_comp->seqlen;
	
	start = 0;
	end = 0;
	score = 0;
	bestScore = 0;
	best_read_score = 0;
	bestHits = 0;
	rc = 1;
	for(t_i = 1; t_i <= *matched_templates; ++t_i) {
		template = matched_templates[t_i];
		/* check if index DB is loaded, and rc */
		if(template < 0) {
			if(rc) {
				qseq_comp->N[0]--;
				comp_rc(qseq_comp);
				qseq_comp->N[0]++;
				qseq_r_comp->N[0]--;
				comp_rc(qseq_r_comp);
				qseq_r_comp->N[0]++;
				
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
				
				rc = 0;
			}
		}
		
		lock(excludeDB);
		if(0 <= template && templates_index[template] == 0) {
			templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], seq_indexes[template], index_indexes[template]);
		} else if(template < 0 && templates_index[-template] == 0) {
			templates_index[-template] = alignLoadPtr(seq_in, index_in, template_lengths[-template], seq_indexes[-template], index_indexes[-template]);
		}
		unlock(excludeDB);
		
		template = abs(template);
		
		/* align qseq */
		alnStat = KMA_score(template, qseq, qseq_comp->seqlen, qseq_comp, points, NWmatrices);
		if(0 < alnStat.score) {
			alnStat_r = KMA_score(template, qseq_r, qseq_r_comp->seqlen, qseq_r_comp, points, NWmatrices);
			
			/* get read score */
			if(0 < alnStat_r.score) {
				aln_len = alnStat.len + alnStat_r.len;
				
				/* Handle negative insertsizes caused by trimming,
				user stupidity or sample error. */
				if(alnStat.pos < alnStat_r.pos) {
					start = alnStat.pos;
					end = alnStat_r.pos + alnStat_r.len - alnStat_r.gaps;
				} else {
					start = alnStat_r.pos;
					end = alnStat.pos + alnStat.len - alnStat.gaps;
				}
				
				read_score = alnStat.score + alnStat_r.score;
				if(start == 0 && end == template_lengths[abs(template)]) {
					read_score += abs(W1);
				}
				score = 1.0 * read_score / aln_len;
			} else {
				read_score = 0;
			}
		} else {
			read_score = 0;
		}
		
		/* save best match(es) */
		if(read_score > kmersize && score >= scoreT) {
			if(score > bestScore) { // save as best match
				bestScore = score;
				best_read_score = read_score;
				*bestTemplates = matched_templates[t_i];
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score > best_read_score) { // save as best match
				bestScore = score;
				best_read_score = read_score;
				*bestTemplates = matched_templates[t_i];
				*best_start_pos = start;
				*best_end_pos = end;
				bestHits = 1;
			} else if(score == bestScore && read_score == best_read_score) { // update best match
				bestTemplates[bestHits] = matched_templates[t_i];
				best_start_pos[bestHits] = start;
				best_end_pos[bestHits] = end;
				++bestHits;
			}
		}
	}
	
	if(best_read_score) {
		/* check direction of qseqs */
		if(*bestTemplates < 0) {
			for(t_i = 0; t_i < bestHits; ++t_i) {
				bestTemplates[t_i] = -bestTemplates[t_i];
			}
			lock(excludeOut);
			update_Scores_pe(qseq_r, qseq_r_comp->seqlen, qseq, qseq_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header_r, header, frag_out_raw);
			unlock(excludeOut);
		} else {
			if(!rc) {
				strrc(qseq, qseq_comp->seqlen);
				strrc(qseq_r, qseq_r_comp->seqlen);
			}
			lock(excludeOut);
			update_Scores_pe(qseq, qseq_comp->seqlen, qseq_r, qseq_r_comp->seqlen, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, frag_out_raw);
			unlock(excludeOut);
		}
	}
}

void printConsensus(struct assem *aligned_assem, char *header, FILE *alignment_out, FILE *consensus_out) {
	
	int i, aln_len, bias;
	
	/* print alignment */
	aln_len = aligned_assem->len;
	fprintf(alignment_out, "# %s\n", header);
	for(i = 0; i < aln_len; i += 60) {
		fprintf(alignment_out, "%-10s\t%.60s\n", "template:", aligned_assem->t + i);
		fprintf(alignment_out, "%-10s\t%.60s\n", "", aligned_assem->s + i);
		fprintf(alignment_out, "%-10s\t%.60s\n\n", "query:", aligned_assem->q + i);
	}
	
	/* Prepare consensus */
	if(ref_fsa) {
		for(i = 0; i < aln_len; ++i) {
			if(aligned_assem->q[i] == '-') {
				aligned_assem->q[i] = 'n';
			}
		}
	} else {
		for(i = 0, bias = 0; i < aln_len; ++i, ++bias) {
			aligned_assem->q[bias] = aligned_assem->q[i];
			if(aligned_assem->q[i] == '-') {
				--bias;
			}
		}
		aln_len = bias;
		aligned_assem->q[aln_len] = 0;
	}
	/* Print consensus */
	fprintf(consensus_out, ">%s\n", header);
	for(i = 0; i < aln_len; i += 60) {
		fprintf(consensus_out, "%.60s\n", aligned_assem->q + i);
	}
}

void * alnFrags_threaded(void * arg) {
	
	static volatile int excludeIn = 0, excludeOut = 0, excludeDB = 0;
	struct aln_thread *thread = arg;
	int rc_flag, read_score, delta, *matched_templates;
	int *bestTemplates, *bestTemplates_r, *best_start_pos, *best_end_pos;
	long *index_indexes, *seq_indexes;
	FILE *inputfile, *frag_out_raw, *index_in, *seq_in;
	struct compDNA *qseq_comp, *qseq_r_comp;
	struct qseqs *qseq, *qseq_r, *header, *header_r;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	
	/* get input */
	matched_templates = thread->matched_templates;
	bestTemplates = thread->bestTemplates;
	bestTemplates_r = thread->bestTemplates_r;
	best_start_pos = thread->best_start_pos;
	best_end_pos = thread->best_end_pos;
	index_indexes = thread->index_indexes;
	seq_indexes = thread->seq_indexes;
	inputfile = thread->inputfile;
	frag_out_raw = thread->frag_out_raw;
	index_in = thread->index_in;
	seq_in = thread->seq_in;
	qseq_comp = thread->qseq_comp;
	qseq_r_comp = thread->qseq_r_comp;
	qseq = thread->qseq;
	qseq_r = thread->qseq_r;
	header = thread->header;
	header_r = thread->header_r;
	points = thread->points;
	NWmatrices = thread->NWmatrices;
	
	delta = qseq->size;
	read_score = 0;
	//lock(&excludeIn);
	lockTime(&excludeIn, 65536);
	while((rc_flag = get_ankers(matched_templates, qseq_comp, header, inputfile)) != 0) {
		if(*matched_templates) { // SE
			read_score = 0;
		} else { // PE
			read_score = get_ankers(matched_templates, qseq_r_comp, header_r, inputfile);
			read_score = labs(read_score);
			qseq_r->len = qseq_r_comp->seqlen;
		}
		unlock(&excludeIn);
		qseq->len = qseq_comp->seqlen;
		
		if(delta <= MAX(qseq->len, qseq_r->len)) {
			delta = MAX(qseq->len, qseq_r->len);
			delta <<= 1;
			qseq->size = delta;
			qseq_r->size = delta;
			free(qseq->seq);
			free(qseq_r->seq);
			qseq->seq = malloc(delta);
			qseq_r->seq = malloc(delta);
			if(!qseq->seq || !qseq_r->seq) {
				ERROR();
			}
		}
		
		if(kmersize <= qseq->len) {
			if(read_score && kmersize <= qseq_r->len) { // PE
				alnFragsPE(matched_templates, qseq_comp, qseq_r_comp, qseq->seq, qseq_r->seq, header, header_r, bestTemplates, bestTemplates_r, best_start_pos, best_end_pos, seq_in, index_in, seq_indexes, index_indexes, frag_out_raw, points, NWmatrices, &excludeOut, &excludeDB);
			} else { // SE
				alnFragsSE(matched_templates, rc_flag, qseq_comp, qseq_r_comp, qseq->seq, qseq_r->seq, qseq->len, header, bestTemplates, best_start_pos, best_end_pos, seq_in, index_in, seq_indexes, index_indexes, frag_out_raw, points, NWmatrices, &excludeOut, &excludeDB);
			}
		}
		lock(&excludeIn);
	}
	unlock(&excludeIn);
	
	return NULL;
}

void getExtendedFeatures(char *template_name, struct assemInfo *matrix, long unsigned *template_seq, int t_len, struct assem *aligned_assem, unsigned fragmentCount, unsigned readCount, FILE *outfile) {
	
	unsigned pos, depthUpdate, maxDepth, nucHighVarSum;
	long unsigned snpSum, insertSum, deletionSum;
	long double var, nucHighVar;
	struct assembly *assembly;
	
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
	
}

int runKMA(char *templatefilename, char *outputfilename, char *exePrev, int ConClave, int extendedFeatures, int vcf) {
	
	int i, j, tmp_template, tmp_tmp_template, file_len, bestTemplate, tot;
	int template, bestHits, t_len, start, end, aln_len, status, rand, sparse;
	int fragCount, fileCount, maxFrag, coverScore, tmp_start, tmp_end, score;
	int stats[4], *matched_templates, *bestTemplates, *bestTemplates_r;
	int *best_start_pos, *best_end_pos;
	unsigned randScore, *fragmentCounts, *readCounts;
	long read_score, best_read_score, *index_indexes, *seq_indexes;
	long unsigned Nhits, template_tot_ulen, bestNum, *w_scores;
	double tmp_score, bestScore, id, q_id, cover, q_cover, p_value;
	long double depth, expected, q_value;
	FILE *inputfile, *frag_in_raw, *index_in, *seq_in, *res_out;
	FILE *alignment_out, *consensus_out, *frag_out_raw, **template_fragments;
	FILE *extendedFeatures_out;
	time_t t0, t1;
	struct FileBuff *frag_out, *frag_out_all, *matrix_out, *vcf_out;
	struct aln *aligned, *gap_align;
	struct assem *aligned_assem;
	struct frag **alignFrags, *alignFrag;
	struct compDNA *qseq_comp, *qseq_r_comp;
	struct qseqs *qseq, *qseq_r, *header, *header_r;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	struct assemble_thread *threads, *thread;
	struct aln_thread *alnThreads, *alnThread;
	
	/* open pipe */
	//inputfile = popen(exePrev, "r");
	status = 0;
	inputfile = kmaPopen(exePrev, "rb");
	if(!inputfile) {
		ERROR();
	} else {
		setvbuf(inputfile, NULL, _IOFBF, CHUNK);
	}
	
	/* load databases */
	file_len = strlen(templatefilename);
	load_DBs_KMA(templatefilename);
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".seq.b");
	seq_in = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	
	strcat(templatefilename, ".index.b");
	index_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	if(!index_in) {
		alignLoadPtr = &alignLoad_fly_build;
		index_in = 0;
		if(kmersize < 4 || 32 < kmersize) {
			kmersize = 16;
		}
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
	index_indexes = smalloc((DB_size + 1) * sizeof(long));
	seq_indexes = smalloc((DB_size + 1) * sizeof(long));
	/* make file indexes of template indexing */
	*index_indexes = 0;
	*seq_indexes = 0;
	index_indexes[1] = sizeof(int);
	seq_indexes[1] = 0;
	for(i = 2; i < DB_size; ++i) {
		index_indexes[i] = index_indexes[i - 1] + (template_lengths[i - 1] << 1) * sizeof(int);
		seq_indexes[i] = seq_indexes[i - 1] + ((template_lengths[i - 1] >> 5) + 1) * sizeof(long unsigned);
	}
	
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag.gz");
		frag_out = gzInitFileBuff(CHUNK);
		openFileBuff(frag_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".aln");
		alignment_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".fsa");
		consensus_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		frag_out_raw = tmpfile();
		if(!frag_out_raw) {
			ERROR();
		}
		if(print_matrix) {
			matrix_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".mat.gz");
			openFileBuff(matrix_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(print_all) {
			strcat(outputfilename, ".frag_raw.gz");
			frag_out_all = gzInitFileBuff(CHUNK);
			openFileBuff(frag_out_all, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			frag_out_all = 0;
		}
		if(vcf) {
			vcf_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".vcf.gz");
			openFileBuff(vcf_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			vcf_out = 0;
		}
	} else {
		fprintf(stderr, " No output file specified!\n");
		exit(2);
	}
	
	fprintf(stderr, "# Running KMA.\n");
	t0 = clock();
	
	/* allocate stuff */
	i = 1;
	alnThreads = 0;
	while(i < thread_num) {
		/* allocate stuff */
		matched_templates = smalloc(((DB_size + 1) << 1) * sizeof(int));
		bestTemplates = smalloc(((DB_size + 1) << 1) * sizeof(int));
		bestTemplates_r = smalloc(((DB_size + 1) << 1) * sizeof(int));
		best_start_pos = smalloc((DB_size << 1) * sizeof(int));
		best_end_pos = smalloc((DB_size << 1) * sizeof(int));
		qseq_comp = smalloc(sizeof(struct compDNA));
		qseq_r_comp = smalloc(sizeof(struct compDNA));
		allocComp(qseq_comp, delta);
		allocComp(qseq_r_comp, delta);
		
		/* allocate matrcies for NW */
		NWmatrices = smalloc(sizeof(struct NWmat));
		NWmatrices->NW_s = delta * delta;
		NWmatrices->NW_q = delta;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		
		/* move it to the thread */
		alnThread = smalloc(sizeof(struct aln_thread));
		alnThread->matched_templates = matched_templates;
		alnThread->bestTemplates = bestTemplates;
		alnThread->bestTemplates_r = bestTemplates_r;
		alnThread->best_start_pos = best_start_pos;
		alnThread->best_end_pos = best_end_pos;
		alnThread->index_indexes = index_indexes;
		alnThread->seq_indexes = seq_indexes;
		alnThread->inputfile = inputfile;
		alnThread->frag_out_raw = frag_out_raw;
		alnThread->index_in = index_in;
		alnThread->seq_in = seq_in;
		alnThread->qseq_comp = qseq_comp;
		alnThread->qseq_r_comp = qseq_r_comp;
		alnThread->qseq = setQseqs(delta);
		alnThread->qseq_r = setQseqs(delta);
		alnThread->header = setQseqs(256);
		alnThread->header_r = setQseqs(256);
		alnThread->points = seedPoint_init(delta);
		alnThread->NWmatrices = NWmatrices;
		
		alnThread->next = alnThreads;
		alnThreads = alnThread;
		
		
		/* start thread */
		if((errno = pthread_create(&alnThread->id, NULL, &alnFrags_threaded, alnThread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			alnThreads = alnThread->next;
			free(alnThread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* allocate stuff */
	matched_templates = smalloc(((DB_size + 1) << 1) * sizeof(int));
	bestTemplates = smalloc(((DB_size + 1) << 1) * sizeof(int));
	bestTemplates_r = smalloc(((DB_size + 1) << 1) * sizeof(int));
	best_start_pos = smalloc((DB_size << 1) * sizeof(int));
	best_end_pos = smalloc((DB_size << 1) * sizeof(int));
	qseq = setQseqs(delta);
	qseq_r = setQseqs(delta);
	header = setQseqs(256);
	header_r = setQseqs(256);
	qseq_comp = smalloc(sizeof(struct compDNA));
	qseq_r_comp = smalloc(sizeof(struct compDNA));
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	points = seedPoint_init(delta);
	
	/* allocate matrcies for NW */
	NWmatrices = smalloc(sizeof(struct NWmat));
	NWmatrices->NW_s = delta * delta;
	NWmatrices->NW_q = delta;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	
	/* strat main thread */
	alnThread = smalloc(sizeof(struct aln_thread));
	alnThread->matched_templates = matched_templates;
	alnThread->bestTemplates = bestTemplates;
	alnThread->bestTemplates_r = bestTemplates_r;
	alnThread->best_start_pos = best_start_pos;
	alnThread->best_end_pos = best_end_pos;
	alnThread->index_indexes = index_indexes;
	alnThread->seq_indexes = seq_indexes;
	alnThread->inputfile = inputfile;
	alnThread->frag_out_raw = frag_out_raw;
	alnThread->index_in = index_in;
	alnThread->seq_in = seq_in;
	alnThread->qseq_comp = qseq_comp;
	alnThread->qseq_r_comp = qseq_r_comp;
	alnThread->qseq = qseq;
	alnThread->qseq_r = qseq_r;
	alnThread->header = header;
	alnThread->header_r = header_r;
	alnThread->points = points;
	alnThread->NWmatrices = NWmatrices;
	alnThread->next = 0;
	
	/* Get alignments */
	alnFrags_threaded(alnThread);
	free(alnThread);
	
	/* join threads */
	for(alnThread = alnThreads; alnThread != 0; alnThread = alnThread->next) {
		/* join thread */
		if((errno = pthread_join(alnThread->id, NULL))) {
			ERROR();
		}
	}
	status |= kmaPclose(inputfile);
	i = 0;
	sfwrite(&i, sizeof(int), 1, frag_out_raw);
	fflush(frag_out_raw);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(header->size < header_r->size) {
		free(header->seq);
		header->size = header_r->size;
		header->seq = smalloc(header->size);
	}
	if(qseq->size < qseq_r->size) {
		free(qseq->seq);
		qseq->size = qseq_r->size;
		qseq->seq = smalloc(qseq->size);
	}
	
	/* clean threads */
	i = 1;
	threads = 0;
	for(alnThread = alnThreads; alnThread != 0; alnThread = alnThreads) {
		alnThreads = alnThread->next;
		
		/* reuse what is possible */
		thread = smalloc(sizeof(struct assemble_thread));
		thread->num = i;
		thread->template = -2;
		thread->frag_out = frag_out;
		thread->NWmatrices = alnThread->NWmatrices;
		thread->qseq = alnThread->qseq;
		thread->header = alnThread->header;
		thread->points = alnThread->points;
		thread->points->len = 0;
		thread->next = threads;
		threads = thread;
		++i;
		
		/* check header and qseq */
		if(header->size < thread->header->size) {
			header->size = thread->header->size;
			free(header->seq);
			header->seq = smalloc(header->size);
		}
		if(qseq->size < thread->qseq->size) {
			qseq->size = thread->qseq->size;
			free(qseq->seq);
			qseq->seq = smalloc(qseq->size);
		}
		
		/* free the rest */
		free(alnThread->matched_templates);
		free(alnThread->bestTemplates);
		free(alnThread->bestTemplates_r);
		free(alnThread->best_start_pos);
		free(alnThread->best_end_pos);
		freeComp(alnThread->qseq_comp);
		free(alnThread->qseq_comp);
		freeComp(alnThread->qseq_r_comp);
		free(alnThread->qseq_r_comp);
		destroyQseqs(alnThread->qseq_r);
		destroyQseqs(alnThread->header_r);
		free(alnThread);
	}
	
	t1 = clock();
	fprintf(stderr, "#\n# KMA mapping time\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select KMA alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped read
	Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(long unsigned));
	if(!alignFrags || !w_scores) {
		ERROR();
	}
	outputfilename[file_len] = 0;
	frag_in_raw = frag_out_raw;
	rewind(frag_in_raw);
	outputfilename[file_len] = 0;
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		ERROR();
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	
	/* Patricks features */
	if(extendedFeatures) {
		fragmentCounts = calloc(DB_size, sizeof(unsigned));
		readCounts = calloc(DB_size, sizeof(unsigned));
		if(!fragmentCounts || !readCounts) {
			ERROR();
		}
	} else {
		fragmentCounts = 0;
		readCounts = 0;
	}
	
	/* Get expected values */
	sparse = 0;
	template_tot_ulen = 0;
	i = DB_size;
	while(--i) {
		template_tot_ulen += template_lengths[i];
	}
	
	/* ConClave */
	if(ConClave == 1) {
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			fread(qseq->seq, 1, qseq->len, frag_in_raw);
			fread(header->seq, 1, header->len, frag_in_raw);
			fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				start = 0;
				end = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					tmp_start = best_start_pos[i];
					tmp_end = best_end_pos[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					//if(tmp_score > bestScore) {
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
						start = tmp_start;
						end = tmp_end;
					} else if(alignment_scores[tmp_template] == best_read_score) {
					//} else if(tmp_score == bestScore) {
						if(tmp_score > bestScore) {
						//if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						//} else if(alignment_scores[tmp_template] == best_read_score) {
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
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(extendedFeatures) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = malloc(sizeof(struct frag));
			if(!alignFrag) {
				ERROR();
			}
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				fread(stats, sizeof(int), 2, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				fread(qseq->seq, 1, qseq->len, frag_in_raw);
				fread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					ERROR();
				}
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags);
		++fileCount;
	} else if(ConClave == 2) {
		/* find potential template candidates */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			/* best templates, skip rest */
			fseek(frag_in_raw, qseq->len + header->len + 2 * bestHits * sizeof(int), SEEK_CUR);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
					if(tmp_tmp_template < 0) {
						tmp_template = -tmp_tmp_template;
					} else {
						tmp_template = tmp_tmp_template;
					}
					tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
					//if(tmp_score > bestScore) {
					if(alignment_scores[tmp_template] > best_read_score) {
						bestTemplate = tmp_tmp_template;
						best_read_score = alignment_scores[tmp_template];
						bestScore = tmp_score;
						bestNum = uniq_alignment_scores[tmp_template];
					} else if(alignment_scores[tmp_template] == best_read_score) {
					//} else if(tmp_score == bestScore) {
						if(tmp_score > bestScore) {
						//if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						//} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
			}
			w_scores[abs(bestTemplate)] += read_score;
			
			if(stats[2] < 0) {
				fread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1], SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* discard insignifiacant templates */
		Nhits = 0;
		template = DB_size;
		while(--template) {
			Nhits += w_scores[template];
		}
		
		template = DB_size;
		while(--template) {
			if((read_score = w_scores[template])) {
				t_len = template_lengths[template];
				//expected = (Nhits - read_score) * (1.0 * t_len) / (template_tot_ulen - t_len + etta);
				expected = t_len;
				expected /= (template_tot_ulen - t_len);
				expected *= (Nhits - read_score);
				//q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= read_score - expected;
				p_value  = p_chisqr(q_value);
				if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
					w_scores[template] = 0;
				}
			}
		}
		
		/* identify sorting keys */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			/* best templates, skip rest */
			fseek(frag_in_raw, qseq->len + header->len + 2 * bestHits * sizeof(int), SEEK_CUR);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			if(bestHits != 1) {
				bestTemplate = 0;
				i = bestHits;
				while(i--) {
					template = abs(bestTemplates[i]);
					if(w_scores[template]) {
						if(bestTemplate) {
							bestTemplate = 0;
							break;
						} else {
							bestTemplate = template;
						}
					}
				}
				
				if(bestTemplate) {
					uniq_alignment_scores[bestTemplate] += read_score;
				}
			}
			
			if(stats[2] < 0) {
				fread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1], SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* choose the templates */
		memset(w_scores, 0, DB_size * sizeof(long unsigned));
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			fread(qseq->seq, 1, qseq->len, frag_in_raw);
			fread(header->seq, 1, header->len, frag_in_raw);
			fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best according to sorting keys */
			if(bestHits != 1) {
				
				bestTemplate = 0;
				bestScore = 0;
				start = 0;
				end = 0;
				
				tot = 0;
				i = bestHits;
				while(i--) {
					tot += uniq_alignment_scores[abs(bestTemplates[i])];
				}
				
				if(tot && 16 <= qseq->len) {
					/* get seed */
					rand = qseq->seq[0];
					i = -1;
					j = qseq->len;
					while(++i < 7) {
						rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
					}
					/* minimal standard */
					rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
					if (rand <= 0) {
						rand += 0x7fffffff;
					}
					
					tmp_score = rand;
					tmp_score /= INT_MAX;
					randScore = tmp_score * tot;
					
					score = 0;
					for(i = 0; i < bestHits; ++i) {
						score += uniq_alignment_scores[abs(bestTemplates[i])];
						if(randScore < score) {
							bestTemplate = bestTemplates[i];
							start = best_start_pos[i];
							end = best_end_pos[i];
							i = bestHits;
						}
					}
					
					if(bestTemplate == 0) {
						tot = 0;
					}
				} else {
					tot = 0;
				}
				
				if(tot == 0) {
					bestTemplate = -1;
					best_read_score = 0;
					bestNum = 0;
					
					/* iterate hits */
					for(i = 0; i < bestHits; ++i) {
						tmp_tmp_template = bestTemplates[i];
						tmp_start = best_start_pos[i];
						tmp_end = best_end_pos[i];
						if(tmp_tmp_template < 0) {
							tmp_template = -tmp_tmp_template;
						} else {
							tmp_template = tmp_tmp_template;
						}
						tmp_score = 1.0 * alignment_scores[tmp_template] / (template_lengths[tmp_template] - kmersize + 1);
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
							start = tmp_start;
							end = tmp_end;
						} else if(alignment_scores[tmp_template] == best_read_score) {
						//} else if(tmp_score == bestScore) {
							if(tmp_score > bestScore) {
							//if(alignment_scores[tmp_template] > best_read_score) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
								start = tmp_start;
								end = tmp_end;
							} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
							//} else if(alignment_scores[tmp_template] == best_read_score) {
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
				}
				
				
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(extendedFeatures) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = malloc(sizeof(struct frag));
			if(!alignFrag) {
				ERROR();
			}
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				fread(stats, sizeof(int), 2, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				fread(qseq->seq, 1, qseq->len, frag_in_raw);
				fread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					ERROR();
				}
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags);
		++fileCount;
	}
	
	fragCount = 0;
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(matched_templates);
	free(bestTemplates);
	destroyQseqs(qseq_r);
	fclose(frag_out_raw);
	if(frag_out_all) {
		destroyGzFileBuff(frag_out_all);
	}
	
	/* Get expected values */
	Nhits = 0;
	i = DB_size;
	while(--i) {
		Nhits += w_scores[i];
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(extendedFeatures) {
		strcat(outputfilename, ".mapstat");
		extendedFeatures_out = sfopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
		fprintf(extendedFeatures_out, "# refSequence\treadCount\tfragmentCount\tmapScoreSum\trefCoveredPositions\trefConsensusSum\tbpTotal\tdepthVariance\tnucHighDepthVariance\tdepthMax\tsnpSum\tinsertSum\tdeletionSum\n");
	} else {
		extendedFeatures_out = 0;
	}
	
	if(vcf) {
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* Get expected values */
	points->len = 0;
	
	/* preallocate assembly matrices */
	delta = MAX(delta, qseq->size);
	aligned_assem = smalloc(sizeof(struct assem));
	matrix = smalloc(sizeof(struct assemInfo));
	matrix->size = delta;
	for(i = 0; i < DB_size; ++i) {
		if(matrix->size < template_lengths[i]) {
			matrix->size = template_lengths[i];
		}
	}
	if(assembly_KMA_Ptr == &assemble_KMA_threaded) {
		matrix->size <<= 1;
	} else {
		matrix->size++;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(struct assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
	/* allocate matrcies for NW */
	for(thread = threads; thread != 0; thread = thread->next) {
		/* get remaining thread info */
		aligned = smalloc(sizeof(struct aln));
		gap_align = smalloc(sizeof(struct aln));
		aligned->t = smalloc((delta + 1) << 1);
		aligned->s = smalloc((delta + 1) << 1);
		aligned->q = smalloc((delta + 1) << 1);
		gap_align->t = smalloc((delta + 1) << 1);
		gap_align->s = smalloc((delta + 1) << 1);
		gap_align->q = smalloc((delta + 1) << 1);
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->file_count = fileCount;
		thread->files = template_fragments;
		thread->aligned_assem = aligned_assem;
		thread->matrix = matrix;
		thread->spin = (sparse < 0) ? 10 : 100;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
		}
	}
	
	/* start main thread */
	aligned = smalloc(sizeof(struct aln));
	gap_align = smalloc(sizeof(struct aln));
	aligned->t = smalloc((delta + 1) << 1);
	aligned->s = smalloc((delta + 1) << 1);
	aligned->q = smalloc((delta + 1) << 1);
	gap_align->t = smalloc((delta + 1) << 1);
	gap_align->s = smalloc((delta + 1) << 1);
	gap_align->q = smalloc((delta + 1) << 1);
	thread = smalloc(sizeof(struct assemble_thread));
	thread->num = 0;
	thread->template = 0;
	thread->file_count = fileCount;
	thread->files = template_fragments;
	thread->frag_out = frag_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->next = 0;
	thread->spin = (sparse < 0) ? 10 : 100;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	for(template = 1; template < DB_size; ++template) {
		if(w_scores[template] > 0) {
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = t_len;
			expected /= (template_tot_ulen - t_len);
			expected *= (Nhits - read_score);
			if(0 < expected) {
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= (read_score - expected);
			} else {
				q_value = read_score;
			}
			p_value  = p_chisqr(q_value);
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len))) {
				/* Do assembly */
				//status |= assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
				thread->template = template;
				assembly_KMA_Ptr(thread);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = aligned_assem->aln_len;
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
				}
				
				if(ID_t <= id && 0 < id) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_names[template], read_score, (unsigned) expected, t_len, id, cover, q_id, q_cover, (double) depth, (double) q_value, p_value);
					printConsensus(aligned_assem, template_names[template], alignment_out, consensus_out);
					/* print matrix */
					if(matrix_out) {
						updateMatrix(matrix_out, template_names[template], templates_index[template]->seq, matrix, t_len);
					}
					if(extendedFeatures) {
						getExtendedFeatures(template_names[template], matrix, templates_index[template]->seq, t_len, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
					}
					if(vcf) {
						updateVcf(template_names[template], templates_index[template]->seq, t_len, matrix, vcf, vcf_out);
					}
				}
			}
		}
	}
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	if(index_in) {
		fclose(index_in);
	}
	fclose(seq_in);
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	destroyGzFileBuff(frag_out);
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(extendedFeatures) {
		fclose(extendedFeatures_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}

int runKMA_MEM(char *templatefilename, char *outputfilename, char *exePrev, int ConClave, int extendedFeatures, int vcf) {
	
	/* runKMA_MEM is a memory saving version of runKMA,
	   at the cost it chooses best templates based on kmers
	   instead of alignment score. */
	
	int i, j, tmp_template, tmp_tmp_template, file_len, score, tot, sparse;
	int template, bestHits, t_len, start, end, aln_len, fragCount, maxFrag;
	int rc_flag, coverScore, tmp_start, tmp_end, bestTemplate, status, rand;
	int fileCount, progress, stats[4];
	int *matched_templates, *bestTemplates, *best_start_pos, *best_end_pos;
	unsigned randScore, *fragmentCounts, *readCounts;
	long best_read_score, read_score, seq_seeker, index_seeker;
	long unsigned Nhits, template_tot_ulen, bestNum, counter, *w_scores;
	double tmp_score, bestScore, id, cover, q_id, q_cover, p_value;
	long double depth, q_value, expected;
	FILE *inputfile, *frag_in_raw, *index_in, *seq_in, *res_out;
	FILE *alignment_out, *consensus_out, *frag_out_raw, **template_fragments;
	FILE *extendedFeatures_out;
	time_t t0, t1;
	struct FileBuff *frag_out, *frag_out_all, *matrix_out, *vcf_out;
	struct aln *aligned, *gap_align;
	struct assem *aligned_assem;
	struct frag **alignFrags, *alignFrag;
	struct compDNA *qseq_comp, *qseq_r_comp;
	struct qseqs *qseq, *qseq_r, *header, *header_r;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	struct assemble_thread *threads, *thread;
	
	/* open pipe */
	status = 0;
	//inputfile = popen(exePrev, "r");
	inputfile = kmaPopen(exePrev, "rb");
	if(!inputfile) {
		ERROR();
	} else {
		setvbuf(inputfile, NULL, _IOFBF, CHUNK);
	}
	
	/* load databases */
	file_len = strlen(templatefilename);
	load_DBs_KMA(templatefilename);
	templatefilename[file_len] = 0;
	strcat(templatefilename, ".index.b");
	index_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	if(index_in) {
		fread(&kmersize, sizeof(int), 1, index_in);
	} else {
		alignLoadPtr = &alignLoad_fly_build_mem;
		index_in = 0;
		if(kmersize < 4 || 32 < kmersize) {
			kmersize = 16;
		}
	}
	
	
	strcat(templatefilename, ".seq.b");
	seq_in = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	/* allocate stuff */
	file_len = strlen(outputfilename);
	qseq_comp = malloc(sizeof(struct compDNA));
	qseq_r_comp = malloc(sizeof(struct compDNA));
	if(!qseq_comp || !qseq_r_comp) {
		ERROR();
	}
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	qseq = setQseqs(delta);
	qseq_r = setQseqs(delta);
	header = setQseqs(256);
	header_r = setQseqs(256);
	points = seedPoint_init(delta);
	
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag.gz");
		frag_out = gzInitFileBuff(CHUNK);
		openFileBuff(frag_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".aln");
		alignment_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".fsa");
		consensus_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		frag_out_raw = tmpfile();
		if(!frag_out_raw) {
			ERROR();
		}
		if(print_matrix) {
			matrix_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".mat.gz");
			openFileBuff(matrix_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(print_all) {
			strcat(outputfilename, ".frag_raw.gz");
			frag_out_all = gzInitFileBuff(CHUNK);
			openFileBuff(frag_out_all, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			frag_out_all = 0;
		}
		if(vcf) {
			vcf_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".vcf.gz");
			openFileBuff(vcf_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			vcf_out = 0;
		}
	} else {
		fprintf(stderr, " No output file specified!\n");
		exit(2);
	}
	
	fprintf(stderr, "# Collecting k-mer scores.\n");
	t0 = clock();
	
	/* Get alignments */
	matched_templates = malloc(((DB_size + 1) << 1) * sizeof(int));
	best_start_pos = calloc((DB_size << 1), sizeof(int));
	best_end_pos = malloc((DB_size << 1) * sizeof(int));
	if(!matched_templates || !best_start_pos || !best_end_pos) {
		ERROR();
	}
	bestTemplates = (matched_templates + 1);
	
	/* consider printPair */
	t_len = 0;
	read_score = 0;
	while((rc_flag = get_ankers(matched_templates, qseq_comp, header, inputfile)) != 0) {
		if(*matched_templates) { // SE
			read_score = 0;
		} else { // PE
			read_score = get_ankers(matched_templates, qseq_r_comp, header_r, inputfile);
			read_score = labs(read_score);
			qseq_r->len = qseq_r_comp->seqlen;
		}
		qseq->len = qseq_comp->seqlen;
		
		if(kmersize <= qseq->len) {
			
			if(delta <= MAX(qseq->len, qseq_r->len)) {
				delta = MAX(qseq->len, qseq_r->len);
				delta <<= 1;
				qseq->size = delta;
				qseq_r->size = delta;
				free(qseq->seq);
				free(qseq_r->seq);
				qseq->seq = malloc(delta);
				qseq_r->seq = malloc(delta);
				if(!qseq->seq || !qseq_r->seq) {
					ERROR();
				}
			}
			unCompDNA(qseq_comp, qseq->seq);
			
			
			/* reverse complement seq */
			best_read_score = abs(rc_flag);
			
			for(i = 1, bestHits = 0; i <= *matched_templates; ++i, ++bestHits) {
				best_end_pos[bestHits] = template_lengths[abs(matched_templates[i])];
			}
			
			if(rc_flag < 0 && 0 < matched_templates[*matched_templates]) {
				bestHits = -bestHits;
			}
			
			if(read_score && kmersize <= qseq_r->len) {
				unCompDNA(qseq_r_comp, qseq_r->seq);
				update_Scores_pe(qseq->seq, qseq->len, qseq_r->seq, qseq_r->len, bestHits, best_read_score + read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, frag_out_raw);
			} else {
				update_Scores(qseq->seq, qseq->len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
			}
			
			/* dump seq to all */
			if(frag_out_all) {
				updateAllFrag(qseq->seq, qseq->len, abs(bestHits), best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_all);
				if(read_score) {
					updateAllFrag(qseq_r->seq, qseq_r->len, abs(bestHits), read_score, best_start_pos, best_end_pos, bestTemplates, header_r, frag_out_all);
				}
			}
		}
	}
	status |= kmaPclose(inputfile);
	i = 0;
	sfwrite(&i, sizeof(int), 1, frag_out_raw);
	fflush(frag_out_raw);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(header->size < header_r->size) {
		destroyQseqs(header);
		header = header_r;
	} else {
		destroyQseqs(header_r);
	}
	header_r = 0;
	if(qseq->size < qseq_r->size) {
		destroyQseqs(qseq);
		qseq = qseq_r;
	} else {
		destroyQseqs(qseq_r);
	}
	qseq_r = 0;
	
	if(frag_out_all) {
		destroyGzFileBuff(frag_out_all);
	}
	t1 = clock();
	fprintf(stderr, "#\n# Time for score collecting:\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select k-mer alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(long unsigned));
	Score = 0;
	if(!alignFrags || !w_scores) {
		ERROR();
	}
	outputfilename[file_len] = 0;
	frag_in_raw = frag_out_raw;
	rewind(frag_in_raw);
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!template_fragments) {
		ERROR();
	}
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	
	/* Patricks features */
	if(extendedFeatures) {
		fragmentCounts = calloc(DB_size, sizeof(unsigned));
		readCounts = calloc(DB_size, sizeof(unsigned));
		if(!fragmentCounts || !readCounts) {
			ERROR();
		}
	} else {
		fragmentCounts = 0;
		readCounts = 0;
	}
	
	/* Get expected values */
	sparse = 0;
	template_tot_ulen = 0;
	i = DB_size;
	while(--i) {
		template_tot_ulen += template_lengths[i];
	}
	
	/* ConClave */
	if(ConClave == 1) {
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			fread(qseq->seq, 1, qseq->len, frag_in_raw);
			fread(header->seq, 1, header->len, frag_in_raw);
			fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				start = 0;
				end = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; ++i) {
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
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(extendedFeatures) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = malloc(sizeof(struct frag));
			if(!alignFrag) {
				ERROR();
			}
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				fread(stats, sizeof(int), 2, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				fread(qseq->seq, 1, qseq->len, frag_in_raw);
				fread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					ERROR();
				}
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags);
		++fileCount;
	} else if(ConClave == 2) {
		/* find potential template candidates */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			/* best templates, skip rest */
			fseek(frag_in_raw, qseq->len + header->len + 2 * bestHits * sizeof(int), SEEK_CUR);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
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
					//} else if(alignment_scores[tmp_template] == best_read_score) {
					} else if(tmp_score == bestScore) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
			}
			w_scores[abs(bestTemplate)] += read_score;
			
			if(stats[2] < 0) {
				fread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1], SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* discard insignifiacant templates */
		Nhits = 0;
		template = DB_size;
		while(--template) {
			Nhits += w_scores[template];
		}
		
		template = DB_size;
		while(--template) {
			if((read_score = w_scores[template])) {
				t_len = template_lengths[template];
				//expected = (Nhits - read_score) * (t_len / (template_tot_ulen - t_len + etta));
				expected = t_len;
				expected /= (template_tot_ulen - t_len);
				expected *= (Nhits - read_score);
				//q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= read_score - expected;
				p_value  = p_chisqr(q_value);
				if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
					w_scores[template] = 0;
				}
			}
		}
		
		/* identify sorting keys */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			if(bestHits != 1) {
				/* best templates, skip rest */
				fseek(frag_in_raw, qseq->len + header->len + 2 * bestHits * sizeof(int), SEEK_CUR);
				fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
				bestTemplate = 0;
				i = bestHits;
				while(i--) {
					template = abs(bestTemplates[i]);
					if(w_scores[template]) {
						if(bestTemplate) {
							bestTemplate = 0;
							break;
						} else {
							bestTemplate = template;
						}
					}
				}
				
				if(bestTemplate) {
					uniq_alignment_scores[bestTemplate] += read_score;
				}
			} else {
				/* skip rest */
				fseek(frag_in_raw, qseq->len + header->len + 3 * sizeof(int), SEEK_CUR);
			}
			
			if(stats[2] < 0) {
				fread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1], SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* choose the templates */
		memset(w_scores, 0, DB_size * sizeof(long unsigned));
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			fread(qseq->seq, 1, qseq->len, frag_in_raw);
			fread(header->seq, 1, header->len, frag_in_raw);
			fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best according to sorting keys */
			if(bestHits != 1) {
				
				bestTemplate = 0;
				bestScore = 0;
				start = 0;
				end = 0;
				
				tot = 0;
				i = bestHits;
				while(i--) {
					tot += uniq_alignment_scores[abs(bestTemplates[i])];
				}
				
				if(tot && 16 <= qseq->len) {
					/* get seed */
					rand = qseq->seq[0];
					i = -1;
					j = qseq->len;
					while(++i < 7) {
						rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
					}
					/* minimal standard */
					rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
					if (rand <= 0) {
						rand += 0x7fffffff;
					}
					
					tmp_score = rand;
					tmp_score /= INT_MAX;
					randScore = tmp_score * tot;
					
					score = 0;
					for(i = 0; i < bestHits; ++i) {
						score += uniq_alignment_scores[abs(bestTemplates[i])];
						if(randScore < score) {
							bestTemplate = bestTemplates[i];
							start = best_start_pos[i];
							end = best_end_pos[i];
							i = bestHits;
						}
					}
					
					if(bestTemplate == 0) {
						tot = 0;
					}
				} else {
					tot = 0;
				}
				
				if(tot == 0) {
					bestTemplate = -1;
					best_read_score = 0;
					bestNum = 0;
					
					/* iterate hits */
					for(i = 0; i < bestHits; ++i) {
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
				}
				
				
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(extendedFeatures) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = malloc(sizeof(struct frag));
			if(!alignFrag) {
				ERROR();
			}
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				fread(stats, sizeof(int), 2, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				fread(qseq->seq, 1, qseq->len, frag_in_raw);
				fread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					ERROR();
				}
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
			
		}
		template_fragments[fileCount] = printFrags(alignFrags);
		++fileCount;
	}
	
	fragCount = 0;
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(matched_templates);
	fclose(frag_out_raw);
	
	/* Get expected values */
	Nhits = 0;
	i = DB_size;
	while(--i) {
		Nhits += w_scores[i];
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(extendedFeatures) {
		strcat(outputfilename, ".mapstat");
		extendedFeatures_out = sfopen(outputfilename, "ab");
		outputfilename[file_len] = 0;
		fprintf(extendedFeatures_out, "# refSequence\treadCount\tfragmentCount\tmapScoreSum\trefCoveredPositions\trefConsensusSum\tbpTotal\tdepthVariance\tnucHighDepthVariance\tdepthMax\tsnpSum\tinsertSum\tdeletionSum\n");
	} else {
		extendedFeatures_out = 0;
	}
	if(vcf) {
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* preallocate assembly matrices */
	matrix = smalloc(sizeof(struct assemInfo));
	aligned_assem = smalloc(sizeof(struct assem));
	matrix->size = delta;
	for(i = 0; i < DB_size; ++i) {
		if(matrix->size < template_lengths[i]) {
			matrix->size = template_lengths[i];
		}
	}
	if(assembly_KMA_Ptr == &assemble_KMA_threaded) {
		matrix->size <<= 1;
	} else {
		matrix->size++;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(struct assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
	/* allocate matrcies for NW */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		/* allocate matrices */
		NWmatrices = smalloc(sizeof(struct NWmat));
		NWmatrices->NW_s = 1024 * 1024;
		NWmatrices->NW_q = 1024;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		
		aligned = smalloc(sizeof(struct aln));
		gap_align = smalloc(sizeof(struct aln));
		aligned->t = smalloc((delta + 1) << 1);
		aligned->s = smalloc((delta + 1) << 1);
		aligned->q = smalloc((delta + 1) << 1);
		gap_align->t = smalloc((delta + 1) << 1);
		gap_align->s = smalloc((delta + 1) << 1);
		gap_align->q = smalloc((delta + 1) << 1);
		
		/* move it to the thread */
		thread = smalloc(sizeof(struct assemble_thread));
		thread->num = i;
		
		thread->template = -2;
		thread->file_count = fileCount;
		thread->files = template_fragments;
		thread->frag_out = frag_out;
		thread->aligned_assem = aligned_assem;
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->NWmatrices = NWmatrices;
		thread->matrix = matrix;
		thread->qseq = setQseqs(qseq->size);
		thread->header = setQseqs(header->size);
		thread->points = seedPoint_init(delta);
		thread->points->len = 0;
		thread->spin = (sparse < 0) ? 10 : 100;
		 
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* start main thread */
	NWmatrices = smalloc(sizeof(struct NWmat));
	NWmatrices->NW_s = 1024 * 1024;
	NWmatrices->NW_q = 1024;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	
	aligned = smalloc(sizeof(struct aln));
	gap_align = smalloc(sizeof(struct aln));
	aligned->t = smalloc((delta + 1) << 1);
	aligned->s = smalloc((delta + 1) << 1);
	aligned->q = smalloc((delta + 1) << 1);
	gap_align->t = smalloc((delta + 1) << 1);
	gap_align->s = smalloc((delta + 1) << 1);
	gap_align->q = smalloc((delta + 1) << 1);
	
	/* move it to the thread */
	thread = smalloc(sizeof(struct assemble_thread));
	thread->num = 0;
	thread->template = 0;
	thread->file_count = fileCount;
	thread->files = template_fragments;
	thread->frag_out = frag_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->next = 0;
	thread->spin = (sparse < 0) ? 10 : 100;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	seq_seeker = 0;
	index_seeker = 0;
	progress = 0;
	counter = 0;
	if(progress) {
		fprintf(stderr, "# Progress:\t%3d%%\r", 0);
		fflush(stderr);
	}
	
	for(template = 1; template < DB_size; ++template) {
		if(w_scores[template] > 0) {
			
			if(progress) {
				counter += w_scores[template];
				fprintf(stderr, "# Progress:\t%3lu%%\r", 100 * counter / Nhits);
				fflush(stderr);
			}
			
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = t_len;
			expected /= (template_tot_ulen - t_len);
			expected *= (Nhits - read_score);
			if(0 < expected) {
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= (read_score - expected);
			} else {
				q_value = read_score;
			}
			p_value  = p_chisqr(q_value);
			
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len))) {
				/* load DB */
				if(index_in) {
					index_seeker *= sizeof(int);
					fseek(index_in, index_seeker, SEEK_CUR);
					index_seeker = 0;
				}
				seq_seeker *= sizeof(long unsigned);
				fseek(seq_in, seq_seeker, SEEK_CUR);
				seq_seeker = 0;
				templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], 0, 0);
				
				/* Do assembly */
				//status |= assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
				thread->template = template;
				assembly_KMA_Ptr(thread);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = aligned_assem->aln_len;
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
				}
				
				if(ID_t <= id && 0 < id) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_names[template], read_score, (unsigned) expected, t_len, id, cover, q_id, q_cover, (double) depth, (double) q_value, p_value);
					printConsensus(aligned_assem, template_names[template], alignment_out, consensus_out);
					/* print matrix */
					if(matrix_out) {
						updateMatrix(matrix_out, template_names[template], templates_index[template]->seq, matrix, t_len);
					}
					if(extendedFeatures) {
						getExtendedFeatures(template_names[template], matrix, templates_index[template]->seq, t_len, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
					}
					if(vcf) {
						updateVcf(template_names[template], templates_index[template]->seq, t_len, matrix, vcf, vcf_out);
					}
				}
				/* destroy this DB index */
				destroyPtr(template);
			} else {
				if(index_in) {
					index_seeker += (template_lengths[template] << 1);
				}
				seq_seeker += ((template_lengths[template] >> 5) + 1);
			}
		} else {
			if(index_in) {
				index_seeker += (template_lengths[template] << 1);
			}
			seq_seeker += ((template_lengths[template] >> 5) + 1);
		}
	}
	
	if(progress) {
		fprintf(stderr, "\n");
	}
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	if(index_in) {
		fclose(index_in);
	}
	fclose(seq_in);
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	destroyGzFileBuff(frag_out);
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(extendedFeatures) {
		fclose(extendedFeatures_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}

void runKMA_Mt1(char *templatefilename, char *outputfilename, char *exePrev, int Mt1, int vcf) {
	
	int i, j, aln_len, t_len, coverScore, file_len;
	long unsigned read_score, seeker;
	double p_value, id, q_id, cover, q_cover;
	long double depth;
	FILE *res_out, *alignment_out, *consensus_out, *template_fragments;
	FILE *DB_file;
	time_t t0, t1;
	struct FileBuff *frag_out, *matrix_out, *vcf_out;
	struct aln *aligned, *gap_align;
	struct assem *aligned_assem;
	struct qseqs *qseq, *header;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	struct assemble_thread *threads, *thread;
	
	/* open pipe */
	//template_fragments = popen(exePrev, "r");
	template_fragments = kmaPopen(exePrev, "rb");
	if(!template_fragments) {
		ERROR();
	} else {
		setvbuf(template_fragments, NULL, _IOFBF, CHUNK);
	}
	
	file_len = strlen(outputfilename);
	header = setQseqs(256);
	qseq = setQseqs(delta);
	points = seedPoint_init(delta);
	
	/* open outputfiles */
	if(outputfilename) {
		strcat(outputfilename, ".res");
		res_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".frag.gz");
		frag_out = gzInitFileBuff(CHUNK);
		openFileBuff(frag_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".aln");
		alignment_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		strcat(outputfilename, ".fsa");
		consensus_out = sfopen(outputfilename, "w");
		outputfilename[file_len] = 0;
		if(print_matrix) {
			matrix_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".mat.gz");
			openFileBuff(matrix_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			matrix_out = 0;
		}
		if(vcf) {
			vcf_out = gzInitFileBuff(CHUNK);
			strcat(outputfilename, ".vcf.gz");
			openFileBuff(vcf_out, outputfilename, "wb");
			outputfilename[file_len] = 0;
		} else {
			vcf_out = 0;
		}
	} else {
		fprintf(stderr, " No output file specified!\n");
		exit(2);
	}
	
	/* load indexing */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".length.b");
	DB_file = sfopen(templatefilename, "rb");
	fread(&DB_size, sizeof(int), 1, DB_file);
	contamination = DB_size;
	/* load lengths */
	template_lengths = malloc(DB_size * sizeof(int));
	if(!template_lengths) {
		ERROR();
	}
	fread(template_lengths, sizeof(int), DB_size, DB_file);
	/*fseek(DB_file, (2 * DB_size) * sizeof(int), SEEK_CUR);
	if(fread(template_lengths, sizeof(int), DB_size, DB_file) == 0) {
		fseek(DB_file, sizeof(int), SEEK_SET);
		fread(template_lengths, sizeof(int), DB_size, DB_file);
	}*/
	templatefilename[file_len] = 0;
	fclose(DB_file);
	if(kmersize < 4) {
		kmersize = *template_lengths;
		if(32 < kmersize || kmersize < 4) {
			kmersize = 16;
		}
	}
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	/* get seq / index */
	seeker = 0;
	for(i = 2; i <= Mt1; ++i) {
		seeker += ((template_lengths[i - 1] >> 5) + 1) * sizeof(long unsigned);
	}
	templates_index = malloc(sizeof(struct hashMap_index*));
	if(!templates_index) {
		ERROR();
	}
	Mt1Score = calloc(1, sizeof(long unsigned));
	*templates_index = calloc(1, sizeof(struct hashMap_index));
	if(!templates_index[0] || !Mt1Score) {
		ERROR();
	}
	
	/* make index */
	*template_lengths = template_lengths[Mt1];
	template_lengths = realloc(template_lengths, sizeof(int));
	if(!template_lengths) {
		ERROR();
	}
	
	strcat(templatefilename, ".seq.b");
	DB_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	*templates_index = alignLoad_fly_build(DB_file, 0, *template_lengths, seeker, 0);
	fclose(DB_file);
	
	/* get name */
	strcat(templatefilename, ".name");
	DB_file = sfopen(templatefilename, "rb");
	/* get size of file */
	fseek(DB_file, 0, SEEK_END);
	seeker = ftell(DB_file);
	/* load file */
	rewind(DB_file);
	template_names = malloc(DB_size * sizeof(char*));
	if(!template_names) {
		ERROR();
	}
	template_names[0] = malloc(seeker);
	if(!template_names[0]) {
		ERROR();
	}
	fread(template_names[0], 1, seeker, DB_file);
	template_names[0][seeker - 1] = 0;
	
	template_names[1] = template_names[0];
	for(i = 0, j = 2; j < DB_size; ++i) {
		if(template_names[0][i] == '\n') {
			template_names[0][i] = 0;
			template_names[j] = template_names[0] + i + 1;
			++j;
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	*template_names = strdup(template_names[Mt1]);
	free(template_names[1]);
	template_names = realloc(template_names, sizeof(char*));
	if(!template_names || !template_names[0]) {
		ERROR();
	}
	
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(vcf) {
		initialiseVcf(vcf_out, templatefilename);
	}
	
	/* preallocate assembly matrices */
	matrix = smalloc(sizeof(struct assemInfo));
	aligned_assem = smalloc(sizeof(struct assem));
	if(assembly_KMA_Ptr == &assemble_KMA_threaded) {
		matrix->size = (*template_lengths) << 1;
	} else {
		matrix->size = (*template_lengths) + 1;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(struct assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
	/* allocate matrcies for NW */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		/* allocate matrices */
		NWmatrices = smalloc(sizeof(struct NWmat));
		NWmatrices->NW_s = 1024 * 1024;
		NWmatrices->NW_q = 1024;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		
		aligned = smalloc(sizeof(struct aln));
		gap_align = smalloc(sizeof(struct aln));
		aligned->t = smalloc((delta + 1) << 1);
		aligned->s = smalloc((delta + 1) << 1);
		aligned->q = smalloc((delta + 1) << 1);
		gap_align->t = smalloc((delta + 1) << 1);
		gap_align->s = smalloc((delta + 1) << 1);
		gap_align->q = smalloc((delta + 1) << 1);
		
		/* move it to the thread */
		thread = smalloc(sizeof(struct assemble_thread));
		thread->num = i;
		
		thread->template = -2;
		thread->file_count = 1;
		thread->files = &template_fragments;
		thread->frag_out = frag_out;
		thread->aligned_assem = aligned_assem;
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->NWmatrices = NWmatrices;
		thread->matrix = matrix;
		thread->qseq = setQseqs(qseq->size);
		thread->header = setQseqs(header->size);
		thread->points = seedPoint_init(delta);
		thread->points->len = 0;
		thread->spin = 10;
		
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* start main thread */
	NWmatrices = smalloc(sizeof(struct NWmat));
	NWmatrices->NW_s = 1024 * 1024;
	NWmatrices->NW_q = 1024;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	
	aligned = smalloc(sizeof(struct aln));
	gap_align = smalloc(sizeof(struct aln));
	aligned->t = smalloc((delta + 1) << 1);
	aligned->s = smalloc((delta + 1) << 1);
	aligned->q = smalloc((delta + 1) << 1);
	gap_align->t = smalloc((delta + 1) << 1);
	gap_align->s = smalloc((delta + 1) << 1);
	gap_align->q = smalloc((delta + 1) << 1);
	
	/* move it to the thread */
	thread = smalloc(sizeof(struct assemble_thread));
	thread->num = 0;
	thread->template = 0;
	thread->file_count = 1;
	thread->files = &template_fragments;
	thread->frag_out = frag_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->next = 0;
	thread->spin = 10;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	/* Do assembly */
	//assemblyPtr(aligned_assem, 0, &template_fragments, 1, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
	assembly_KMA_Ptr(thread);
	
	/* make p_value */
	read_score = *Mt1Score;
	t_len = *template_lengths;
	p_value  = p_chisqr(read_score);
	
	if(cmp((p_value <= evalue && read_score > 0), read_score >= scoreT * t_len)) {
		
		/* Depth, ID and coverage */
		if(aligned_assem->cover > 0) {
			coverScore = aligned_assem->cover;
			depth = aligned_assem->depth;
			depth /= t_len;
			id = 100.0 * coverScore / t_len;
			aln_len = aligned_assem->aln_len;
			q_id = 100.0 * coverScore / aln_len;
			cover = 100.0 * aln_len / t_len;
			q_cover = 100.0 * t_len / aln_len;
		} else {
			id = 0;
		}
		if(ID_t <= id && 0 < id) {
			/* Output result */
			fprintf(res_out, "%-12s\t%8lu\t%8d\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
				*template_names, read_score, 0, t_len, id, cover, q_id, q_cover, (double) depth, (double) read_score, p_value);
			printConsensus(aligned_assem, *template_names, alignment_out, consensus_out);
			/* print matrix */
			if(matrix_out) {
				updateMatrix(matrix_out, *template_names, templates_index[0]->seq, matrix, t_len);
			}
			if(vcf) {
				updateVcf(*template_names, templates_index[0]->seq, t_len, matrix, vcf, vcf_out);
			}
		}
		/* destroy this DB index */
		destroyPtr(0);
	}
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	destroyGzFileBuff(frag_out);
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
}

void getIndexName(char ***Template_names, struct hashMap_index ***Templates_index, int oldBias, int bias, int newBias, char *templatefilename) {
	
	int file_len, shmid, dbSize;
	long unsigned i, j, file_size;
	char **template_names;
	FILE *DB_file;
	key_t key;
	struct hashMap_index **templates_index;
	
	template_names = *Template_names;
	templates_index = *Templates_index;
	dbSize = newBias - bias;
	if(template_names && (shm & 8) == 0) {
		template_names += oldBias;
		free(*template_names);
		free(template_names);
	}
	if(templates_index && (shm & 16) == 0) {
		templates_index += oldBias;
		free(templates_index);
	}
	template_names = smalloc(dbSize * sizeof(char *));
	templates_index = calloc(dbSize, sizeof(struct hashMap_index*));
	if(!templates_index) {
		ERROR();
	}
	
	/* load names */
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".name");
	DB_file = sfopen(templatefilename, "rb");
	
	/* get size of file */
	fseek(DB_file, 0, SEEK_END);
	file_size = ftell(DB_file);
	
	/* load file */
	key = ftok(templatefilename, 'n');
	shmid = shmget(key, file_size, 0666);
	if(shm & 16 || 0 <= shmid) {
		key = ftok(templatefilename, 'n');
		shmid = shmget(key, file_size, 0666);
		template_names[0] = shmat(shmid, NULL, 0);
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < dbSize; ++i) {
			if(template_names[0][i] == 0) {
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	} else {
		rewind(DB_file);
		template_names[0] = smalloc(file_size);
		sfread(template_names[0], 1, file_size, DB_file);
		template_names[0][file_size - 1] = 0;
		template_names[1] = template_names[0];
		for(i = 0, j = 2; j < dbSize; ++i) {
			if(template_names[0][i] == '\n') {
				template_names[0][i] = 0;
				template_names[j] = template_names[0] + i + 1;
				++j;
			}
		}
	}
	templatefilename[file_len] = 0;
	fclose(DB_file);
	
	template_names -= bias;
	templates_index -= bias;
	*Template_names = template_names;
	*Templates_index = templates_index;
}

int runKMA_spltDB(char **templatefilenames, int targetNum, char *outputfilename, int argc, char **argv, int ConClave, int extendedFeatures, int vcf) {
	
	/* https://www.youtube.com/watch?v=LtXEMwSG5-8 */
	
	int i, j, k, tmp_template, tmp_tmp_template, t_len, file_len, score, tot;
	int template, bestHits, start, end, aln_len, fragCount, maxFrag, sparse;
	int rc_flag, coverScore, tmp_start, tmp_end, bestTemplate, status, rand;
	int fileCount, stats[4], *bestTargets, (*targetInfo)[6], (*ptrInfo)[6];
	int *matched_templates, *bestTemplates, *best_start_pos, *best_end_pos;
	int progress;
	unsigned randScore, num, target, targetScore, bias;
	unsigned *fragmentCounts, *readCounts, *nums, *uPtr, *dbBiases;
	long best_read_score, read_score, seq_seeker, index_seeker;
	long unsigned Nhits, template_tot_ulen, bestNum, counter, *w_scores;
	double tmp_score, bestScore, id, cover, q_id, q_cover, p_value;
	long double depth, q_value, expected;
	char *templatefilename, Date[11];
	FILE **inputfiles, *inputfile, *frag_in_raw, *index_in, *seq_in;
	FILE *res_out, *alignment_out, *consensus_out, *frag_out_raw;
	FILE *extendedFeatures_out, **template_fragments;
	time_t t0, t1;
	struct tm *tm;
	struct FileBuff *frag_out, *frag_out_all, *matrix_out, *vcf_out;
	struct aln *aligned, *gap_align;
	struct assem *aligned_assem;
	struct frag **alignFrags, *alignFrag;
	struct compDNA *qseq_comp, *qseq_r_comp;
	struct qseqs *qseq, *qseq_r, *header, *header_r;
	struct assemInfo *matrix;
	struct alnPoints *points;
	struct NWmat *NWmatrices;
	struct assemble_thread *threads, *thread;
	
	if(!outputfilename) {
		fprintf(stderr, " No output file specified!\n");
		exit(2);
	}
	
	/* load databases */
	status = 0;
	inputfiles = smalloc(targetNum * sizeof(FILE *));
	nums = smalloc(targetNum * sizeof(unsigned));
	dbBiases = smalloc((targetNum + 1) * sizeof(unsigned)); /* set these */
	bestTargets = smalloc((targetNum + 1) * sizeof(int));
	targetInfo = smalloc(targetNum * 6 * sizeof(int));
	DB_size = 0;
	template_lengths = NULL;
	template_names = NULL;
	templates_index = NULL;
	for(i = 0; i < targetNum; ++i) {
		/* get length, size and bias */
		templatefilename = templatefilenames[i];
		file_len = strlen(templatefilename);
		strcat(templatefilename, ".length.b");
		inputfile = sfopen(templatefilename, "rb");
		sfread(&bias, sizeof(int), 1, inputfile);
		dbBiases[i] = DB_size;
		DB_size += bias;
		template_lengths = realloc(template_lengths, DB_size * sizeof(int));
		if(!template_lengths) {
			ERROR();
		}
		fread(template_lengths + dbBiases[i], sizeof(int), DB_size, inputfile);
		templatefilename[file_len] = 0;
		fclose(inputfile);
	}
	dbBiases[i] = DB_size;
	if(kmersize < 4 || 32 < kmersize) {
		kmersize = 16;
	}
	
	/* get scoring arrays */
	matched_templates = smalloc(((DB_size + 1) << 1) * sizeof(int));
	bestTemplates = (matched_templates + 1);
	best_start_pos = calloc((DB_size << 1), sizeof(int));
	best_end_pos = smalloc((DB_size << 1) * sizeof(int));
	alignment_scores = calloc(DB_size, sizeof(long unsigned));
	uniq_alignment_scores = calloc(DB_size, sizeof(long unsigned));
	if(!best_start_pos || !alignment_scores || !uniq_alignment_scores) {
		ERROR();
	}
	
	/* allocate stuff */
	qseq_comp = malloc(sizeof(struct compDNA));
	qseq_r_comp = malloc(sizeof(struct compDNA));
	if(!qseq_comp || !qseq_r_comp) {
		ERROR();
	}
	allocComp(qseq_comp, delta);
	allocComp(qseq_r_comp, delta);
	qseq = setQseqs(delta);
	qseq_r = setQseqs(delta);
	header = setQseqs(256);
	header_r = setQseqs(256);
	points = seedPoint_init(delta);
	
	/* open outputfiles */
	file_len = strlen(outputfilename);
	strcat(outputfilename, ".res");
	res_out = sfopen(outputfilename, "w");
	outputfilename[file_len] = 0;
	strcat(outputfilename, ".frag.gz");
	frag_out = gzInitFileBuff(CHUNK);
	openFileBuff(frag_out, outputfilename, "wb");
	outputfilename[file_len] = 0;
	strcat(outputfilename, ".aln");
	alignment_out = sfopen(outputfilename, "w");
	outputfilename[file_len] = 0;
	strcat(outputfilename, ".fsa");
	consensus_out = sfopen(outputfilename, "w");
	outputfilename[file_len] = 0;
	frag_out_raw = tmpfile();
	if(!frag_out_raw) {
		ERROR();
	}
	if(print_matrix) {
		matrix_out = gzInitFileBuff(CHUNK);
		strcat(outputfilename, ".mat.gz");
		openFileBuff(matrix_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		matrix_out = 0;
	}
	if(print_all) {
		strcat(outputfilename, ".frag_raw.gz");
		frag_out_all = gzInitFileBuff(CHUNK);
		openFileBuff(frag_out_all, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		frag_out_all = 0;
	}
	if(vcf) {
		vcf_out = gzInitFileBuff(CHUNK);
		strcat(outputfilename, ".vcf.gz");
		openFileBuff(vcf_out, outputfilename, "wb");
		outputfilename[file_len] = 0;
	} else {
		vcf_out = 0;
	}
	if(extendedFeatures) {
		strcat(outputfilename, ".mapstat");
		extendedFeatures_out = sfopen(outputfilename, "wb");
		outputfilename[file_len] = 0;
		fprintf(extendedFeatures_out, "## method\tKMA\n");
		fprintf(extendedFeatures_out, "## version\t%d.%d.%d\n", version[0], version[1], version[2]);
		
		fprintf(extendedFeatures_out, "## databases %s", noFolder(*templatefilenames));
		for(i = 1; i < targetNum; ++i) {
			fprintf(extendedFeatures_out, ", %s", noFolder(templatefilenames[i]));
		}
		fprintf(extendedFeatures_out, "\n");
		
		time(&t1);
		tm = localtime(&t1);
		strftime(Date, sizeof(Date), "%Y-%m-%d", tm);
		fprintf(extendedFeatures_out, "## date\t%s\n", Date);
		fprintf(extendedFeatures_out, "## command\t%s", *argv);
		for(i = 1; i < argc; ++i) {
			fprintf(extendedFeatures_out, " %s", argv[i]);
		}
		fprintf(extendedFeatures_out, "\n");
	} else {
		extendedFeatures_out = 0;
	}
	
	/* open input streams */
	file_len = strlen(outputfilename);
	for(i = 0; i < targetNum; ++i) {
		sprintf(outputfilename + file_len, ".%d", i);
		while(!(inputfile = fopen(outputfilename, "rb"))) {
			usleep(100);
		}
		setvbuf(inputfile, NULL, _IOFBF, CHUNK);
		inputfiles[i] = inputfile;
		outputfilename[file_len] = 0;
	}
	
	fprintf(stderr, "# Collecting k-mer scores.\n");
	t0 = clock();
	
	/* collect from DB mappings */
	uPtr = nums + (i = targetNum);
	ptrInfo = targetInfo + i;
	while(i--) {
		sfread(--uPtr, sizeof(unsigned), 1, inputfiles[i]);
		if(*uPtr) {
			sfread(--ptrInfo, sizeof(int), 6, inputfiles[i]);
		} else {
			(*--ptrInfo)[3] = I_LIMIT;
		}
	}
	
	/* here */
	/* somthing is wrong in the next part,
	which cannot be debugged on macOs. */
	
	target = 0;
	targetScore = 0;
	rc_flag = 0;
	qseq->len = 0;
	qseq_comp->seqlen = 0;
	*bestTargets = 0;
	while(target != U_LIMIT) {
		/* join best templates */
		read_score = 0;
		*matched_templates = 0;
		qseq->len = 0;
		qseq_r->len = 0;
		for(i = 1; i <= *bestTargets; ++i) {
			num = bestTargets[i];
			if(i == 1) {
				nums[num] = get_ankers_spltDB(targetInfo[num], matched_templates + (*matched_templates + 1), qseq_comp, header, inputfiles[num]);
			} else {
				nums[num] = get_ankers_spltDB(targetInfo[num], matched_templates + (*matched_templates + 1), 0, 0, inputfiles[num]);
			}
			qseq->len = qseq_comp->seqlen;
			if(matched_templates[*matched_templates + 1]) {
				read_score = 0;
			} else { // PE
				if(qseq_r->len == 0) {
					nums[num] = get_ankers_spltDB(targetInfo[num], matched_templates + (*matched_templates + 1), qseq_r_comp, header_r, inputfiles[num]);
					qseq_r->len = qseq_r_comp->seqlen;
				} else {
					nums[num] = get_ankers_spltDB(targetInfo[num], matched_templates + (*matched_templates + 1), 0, 0, inputfiles[num]);
				}
				read_score = 1;
			}
			
			/* bias the templates */
			bias = dbBiases[num];
			j = *matched_templates + 1;
			*matched_templates += matched_templates[*matched_templates + 1];
			while((k = j++) <= *matched_templates) {
				matched_templates[k] = matched_templates[j] + bias;
			}
		}
		
		if(kmersize <= qseq->len) {
			if(delta <= MAX(qseq->len, qseq_r->len)) {
				delta = MAX(qseq->len, qseq_r->len);
				delta <<= 1;
				qseq->size = delta;
				qseq_r->size = delta;
				free(qseq->seq);
				free(qseq_r->seq);
				qseq->seq = malloc(delta);
				qseq_r->seq = malloc(delta);
				if(!qseq->seq || !qseq_r->seq) {
					ERROR();
				}
			}
			unCompDNA(qseq_comp, qseq->seq);
			
			best_read_score = targetScore;
			for(i = 1, bestHits = 0; i <= *matched_templates; ++i, ++bestHits) {
				best_end_pos[bestHits] = template_lengths[abs(matched_templates[i])];
			}
			
			if(rc_flag < 0 && 0 < matched_templates[*matched_templates]) {
				bestHits = -bestHits;
			}
			
			if(read_score && kmersize <= qseq_r->len) {
				unCompDNA(qseq_r_comp, qseq_r->seq);
				update_Scores_pe(qseq->seq, qseq->len, qseq_r->seq, qseq_r->len, bestHits, best_read_score + read_score, best_start_pos, best_end_pos, bestTemplates, header, header_r, frag_out_raw);
			} else {
				update_Scores(qseq->seq, qseq->len, bestHits, best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_raw);
			}
			
			/* dump seq to all */
			if(frag_out_all) {
				updateAllFrag(qseq->seq, qseq->len, abs(bestHits), best_read_score, best_start_pos, best_end_pos, bestTemplates, header, frag_out_all);
				if(read_score) {
					updateAllFrag(qseq_r->seq, qseq_r->len, abs(bestHits), read_score, best_start_pos, best_end_pos, bestTemplates, header_r, frag_out_all);
				}
			}
		}
		
		if(*matched_templates) {
			++target;
			i = targetNum;
			uPtr = nums + (i - 1);
			while(i--) {
				if(target == *uPtr) {
					fseek(inputfiles[i], targetInfo[i][1] * sizeof(long unsigned) + (targetInfo[i][2] + targetInfo[i][4]) * sizeof(int) + targetInfo[i][5], SEEK_CUR);
					sfread(uPtr, sizeof(unsigned), 1, inputfiles[i]);
					if(*uPtr != U_LIMIT) {
						sfread(targetInfo[i], sizeof(int), 6, inputfiles[i]);
					} else {
						targetInfo[i][3] = I_LIMIT;
					}
				}
				--uPtr;
			}
		}
		
		/* get best templates for next read */
		uPtr = nums;
		target = U_LIMIT;
		targetScore = U_LIMIT;
		rc_flag = 0;
		*bestTargets = 0;
		for(i = 0; i < targetNum; ++i) {
			if(*uPtr < target) {
				target = *uPtr;
				targetScore = abs(targetInfo[i][3]);
				*bestTargets = 1;
				bestTargets[1] = i;
			} else if(*uPtr == target) {
				if(targetScore < abs(targetInfo[i][3])) {
					target = *uPtr;
					rc_flag = targetInfo[i][3];
					targetScore = abs(rc_flag);
					*bestTargets = 1;
					bestTargets[1] = i;
				} else if(targetScore == abs(targetInfo[i][3])) {
					bestTargets[++*bestTargets] = i;
					rc_flag = rc_flag < 0 ? rc_flag : targetInfo[i][3];
				} else if(*uPtr == U_LIMIT) {
					targetInfo[i][3] = I_LIMIT;
				} else {
					fseek(inputfiles[i], targetInfo[i][1] * sizeof(long unsigned) + (targetInfo[i][2] + targetInfo[i][4]) * sizeof(int) + targetInfo[i][5], SEEK_CUR);
					if(matched_templates[*matched_templates + 1]) { // PE
						nums[i] = get_ankers_spltDB(targetInfo[i], matched_templates + (*matched_templates + 1), qseq_comp, header, inputfiles[i]);
						qseq_r->len = qseq_r_comp->seqlen;
						read_score = 1;
					}
					sfread(uPtr, sizeof(unsigned), 1, inputfiles[i]);
					if(*uPtr != U_LIMIT) {
						sfread(targetInfo[i], sizeof(int), 6, inputfiles[i]);
					} else {
						targetInfo[i][3] = I_LIMIT;
					}
				}
			}
			++uPtr;
		}
	}
	
	/* get fragmentCount */
	if(extendedFeatures) {
		fprintf(extendedFeatures_out, "## fragmentCount\t%u\n", **targetInfo);
		fprintf(extendedFeatures_out, "# refSequence\treadCount\tfragmentCount\tmapScoreSum\trefCoveredPositions\trefConsensusSum\tbpTotal\tdepthVariance\tnucHighDepthVariance\tdepthMax\tsnpSum\tinsertSum\tdeletionSum\n");
	}
	
	/* close files */
	i = targetNum;
	while(i--) {
		fclose(inputfiles[i]);
	}
	
	i = 0;
	sfwrite(&i, sizeof(int), 1, frag_out_raw);
	fflush(frag_out_raw);
	freeComp(qseq_comp);
	free(qseq_comp);
	freeComp(qseq_r_comp);
	free(qseq_r_comp);
	if(header->size < header_r->size) {
		destroyQseqs(header);
		header = header_r;
	} else {
		destroyQseqs(header_r);
	}
	header_r = 0;
	if(qseq->size < qseq_r->size) {
		destroyQseqs(qseq);
		qseq = qseq_r;
	} else {
		destroyQseqs(qseq_r);
	}
	qseq_r = 0;
	if(frag_out_all) {
		destroyGzFileBuff(frag_out_all);
	}
	t1 = clock();
	fprintf(stderr, "#\n# Time for score collecting:\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Sort, output and select k-mer alignments.\n");
	t0 = clock();
	
	/* Get best template for each mapped deltamer/read */
	/* Best hit chosen as: highest mapping score then higest # unique maps */
	Score = 0;
	alignFrags = calloc(DB_size, sizeof(struct frag*));
	w_scores = calloc(DB_size, sizeof(long unsigned));
	template_fragments = calloc(DB_size, sizeof(FILE*));
	if(!alignFrags || !w_scores || !template_fragments) {
		ERROR();
	}
	frag_in_raw = frag_out_raw;
	rewind(frag_in_raw);
	fragCount = 0;
	fileCount = 0;
	maxFrag = 1000000;
	
	/* Patricks features */
	if(extendedFeatures) {
		fragmentCounts = calloc(DB_size, sizeof(unsigned));
		readCounts = calloc(DB_size, sizeof(unsigned));
		if(!fragmentCounts || !readCounts) {
			ERROR();
		}
	} else {
		fragmentCounts = 0;
		readCounts = 0;
	}
	
	/* Get expected values */
	sparse = 0;
	template_tot_ulen = 0;
	i = DB_size;
	while(--i) {
		template_tot_ulen += template_lengths[i];
	}
	
	/* ConClave */
	if(ConClave == 1) {
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			fread(qseq->seq, 1, qseq->len, frag_in_raw);
			fread(header->seq, 1, header->len, frag_in_raw);
			fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				start = 0;
				end = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; ++i) {
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
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(extendedFeatures) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = malloc(sizeof(struct frag));
			if(!alignFrag) {
				ERROR();
			}
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				fread(stats, sizeof(int), 2, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				fread(qseq->seq, 1, qseq->len, frag_in_raw);
				fread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					ERROR();
				}
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
		}
		template_fragments[fileCount] = printFrags(alignFrags);
		++fileCount;
	} else if(ConClave == 2) {
		/* find potential template candidates */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			/* best templates, skip rest */
			fseek(frag_in_raw, qseq->len + header->len + 2 * bestHits * sizeof(int), SEEK_CUR);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best */
			if(bestHits > 1) {
				bestTemplate = -1;
				bestScore = 0;
				best_read_score = 0;
				bestNum = 0;
				/* iterate hits */
				for(i = 0; i < bestHits; ++i) {
					tmp_tmp_template = bestTemplates[i];
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
					//} else if(alignment_scores[tmp_template] == best_read_score) {
					} else if(tmp_score == bestScore) {
						//if(tmp_score > bestScore) {
						if(alignment_scores[tmp_template] > best_read_score) {
							bestTemplate = tmp_tmp_template;
							best_read_score = alignment_scores[tmp_template];
							bestScore = tmp_score;
							bestNum = uniq_alignment_scores[tmp_template];
						//} else if(tmp_score == bestScore && alignment_scores[tmp_template] > bestNum) {
						} else if(alignment_scores[tmp_template] == best_read_score) {
							if(uniq_alignment_scores[tmp_template] > bestNum) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							} else if(uniq_alignment_scores[tmp_template] == bestNum && tmp_template < abs(bestTemplate)) {
								bestTemplate = tmp_tmp_template;
								best_read_score = alignment_scores[tmp_template];
								bestScore = tmp_score;
								bestNum = uniq_alignment_scores[tmp_template];
							}
						}
					}
				}
			} else {
				bestTemplate = *bestTemplates;
			}
			w_scores[abs(bestTemplate)] += read_score;
			
			if(stats[2] < 0) {
				fread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1], SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* discard insignifiacant templates */
		Nhits = 0;
		template = DB_size;
		while(--template) {
			Nhits += w_scores[template];
		}
		
		template = DB_size;
		while(--template) {
			if((read_score = w_scores[template])) {
				t_len = template_lengths[template];
				//expected = (Nhits - read_score) * (t_len / (template_tot_ulen - t_len + etta));
				expected = t_len;
				expected /= (template_tot_ulen - t_len);
				expected *= (Nhits - read_score);
				//q_value = pow(read_score - expected, 2) / (expected + read_score + etta);
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= read_score - expected;
				p_value  = p_chisqr(q_value);
				if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len)) == 0) {
					w_scores[template] = 0;
				}
			}
		}
		
		/* identify sorting keys */
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			if(bestHits != 1) {
				/* best templates, skip rest */
				fseek(frag_in_raw, qseq->len + header->len + 2 * bestHits * sizeof(int), SEEK_CUR);
				fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
				bestTemplate = 0;
				i = bestHits;
				while(i--) {
					template = abs(bestTemplates[i]);
					if(w_scores[template]) {
						if(bestTemplate) {
							bestTemplate = 0;
							break;
						} else {
							bestTemplate = template;
						}
					}
				}
				
				if(bestTemplate) {
					uniq_alignment_scores[bestTemplate] += read_score;
				}
			} else {
				/* skip rest */
				fseek(frag_in_raw, qseq->len + header->len + 3 * sizeof(int), SEEK_CUR);
			}
			
			if(stats[2] < 0) {
				fread(stats, sizeof(int), 2, frag_in_raw);
				fseek(frag_in_raw, stats[0] + stats[1], SEEK_CUR);
			}
		}
		rewind(frag_in_raw);
		
		/* choose the templates */
		memset(w_scores, 0, DB_size * sizeof(long unsigned));
		while(fread(stats, sizeof(int), 4, frag_in_raw) && stats[0] != 0) {
			qseq->len = stats[0];
			sparse = stats[1];
			bestHits = abs(sparse);
			read_score = abs(stats[2]);
			header->len = stats[3];
			
			fread(qseq->seq, 1, qseq->len, frag_in_raw);
			fread(header->seq, 1, header->len, frag_in_raw);
			fread(best_start_pos, sizeof(int), bestHits, frag_in_raw);
			fread(best_end_pos, sizeof(int), bestHits, frag_in_raw);
			fread(bestTemplates, sizeof(int), bestHits, frag_in_raw);
			
			/* Several mapped templates, choose best according to sorting keys */
			if(bestHits != 1) {
				
				bestTemplate = 0;
				bestScore = 0;
				start = 0;
				end = 0;
				
				tot = 0;
				i = bestHits;
				while(i--) {
					tot += uniq_alignment_scores[abs(bestTemplates[i])];
				}
				
				if(tot && 16 <= qseq->len) {
					/* get seed */
					rand = qseq->seq[0];
					i = -1;
					j = qseq->len;
					while(++i < 7) {
						rand = (((rand << 2) | qseq->seq[i]) << 2) | qseq->seq[--j];
					}
					/* minimal standard */
					rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
					if (rand <= 0) {
						rand += 0x7fffffff;
					}
					
					tmp_score = rand;
					tmp_score /= INT_MAX;
					randScore = tmp_score * tot;
					
					score = 0;
					for(i = 0; i < bestHits; ++i) {
						score += uniq_alignment_scores[abs(bestTemplates[i])];
						if(randScore < score) {
							bestTemplate = bestTemplates[i];
							start = best_start_pos[i];
							end = best_end_pos[i];
							i = bestHits;
						}
					}
					
					if(bestTemplate == 0) {
						tot = 0;
					}
				} else {
					tot = 0;
				}
				
				if(tot == 0) {
					bestTemplate = -1;
					best_read_score = 0;
					bestNum = 0;
					
					/* iterate hits */
					for(i = 0; i < bestHits; ++i) {
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
				}
				
				
			} else {
				bestTemplate = *bestTemplates;
				start = *best_start_pos;
				end = *best_end_pos;
			}
			
			/* reverse complement seq */
			if(bestTemplate < 0) {
				bestTemplate = -bestTemplate;
				strrc(qseq->seq, qseq->len);
			}
			w_scores[bestTemplate] += read_score;
			if(extendedFeatures) {
				fragmentCounts[bestTemplate]++;
				readCounts[bestTemplate]++;
			}
			
			/* dump frag info */
			alignFrag = malloc(sizeof(struct frag));
			if(!alignFrag) {
				ERROR();
			}
			alignFrag->buffer[0] = qseq->len;
			alignFrag->buffer[1] = bestHits;
			alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
			alignFrag->buffer[3] = start;
			alignFrag->buffer[4] = end;
			alignFrag->buffer[5] = header->len;
			alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
			alignFrag->header = ustrdup(header->seq, header->len);
			alignFrag->next = alignFrags[bestTemplate];
			alignFrags[bestTemplate] = alignFrag;
			
			++fragCount;
			
			if(stats[2] < 0) {
				if(extendedFeatures) {
					readCounts[bestTemplate]++;
				}
				fread(stats, sizeof(int), 2, frag_in_raw);
				qseq->len = stats[0];
				header->len = stats[1];
				fread(qseq->seq, 1, qseq->len, frag_in_raw);
				fread(header->seq, 1, header->len, frag_in_raw);
				/* dump frag info */
				alignFrag = malloc(sizeof(struct frag));
				if(!alignFrag) {
					ERROR();
				}
				alignFrag->buffer[0] = qseq->len;
				alignFrag->buffer[1] = bestHits;
				alignFrag->buffer[2] = (sparse < 0) ? 0 : read_score;
				alignFrag->buffer[3] = start;
				alignFrag->buffer[4] = end;
				alignFrag->buffer[5] = header->len;
				alignFrag->qseq = ustrdup(qseq->seq, qseq->len);
				alignFrag->header = ustrdup(header->seq, header->len);
				alignFrag->next = alignFrags[bestTemplate];
				alignFrags[bestTemplate] = alignFrag;
				
				++fragCount;
			}
			
			if(fragCount >= maxFrag) {
				template_fragments[fileCount] = printFrags(alignFrags);
				++fileCount;
				fragCount = 0;
				/* control fileamount */
				if(fileCount >= DB_size) {
					template_fragments = realloc(template_fragments, (fileCount + 1) * sizeof(FILE*));
					if(!template_fragments) {
						ERROR();
					}
				}
			}
			
		}
		template_fragments[fileCount] = printFrags(alignFrags);
		++fileCount;
	}
	
	fragCount = 0;
	free(alignFrags);
	free(best_start_pos);
	free(best_end_pos);
	free(matched_templates);
	fclose(frag_out_raw);
	
	/* Get expected values */
	Nhits = 0;
	i = DB_size;
	while(--i) {
		Nhits += w_scores[i];
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time for sorting and outputting KMA alignment\t%.2f s.\n", difftime(t1, t0) / 1000000);
	fprintf(stderr, "#\n# Doing local assemblies of found templates, and output results\n");
	t0 = clock();
	
	/* print heading of resistance file: */
	fprintf(res_out, "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n");
	if(vcf) {
		/* here */
		/* might be done more elegant */
		templatefilename = 0;
		initialiseVcf(vcf_out, templatefilename);
	}
	seq_in = 0;
	index_in = 0;
	
	/* preallocate assembly matrices */
	matrix = smalloc(sizeof(struct assemInfo));
	aligned_assem = smalloc(sizeof(struct assem));
	matrix->size = delta;
	for(i = 0; i < DB_size; ++i) {
		if(matrix->size < template_lengths[i]) {
			matrix->size = template_lengths[i];
		}
	}
	if(assembly_KMA_Ptr == &assemble_KMA_threaded) {
		matrix->size <<= 1;
	} else {
		matrix->size++;
	}
	matrix->assmb = smalloc(matrix->size * sizeof(struct assembly));
	aligned_assem->size = matrix->size;
	aligned_assem->t = smalloc(aligned_assem->size);
	aligned_assem->s = smalloc(aligned_assem->size);
	aligned_assem->q = smalloc(aligned_assem->size);
	
	/* allocate matrcies for NW */
	i = 1;
	threads = 0;
	while(i < thread_num) {
		/* allocate matrices */
		NWmatrices = smalloc(sizeof(struct NWmat));
		NWmatrices->NW_s = 1024 * 1024;
		NWmatrices->NW_q = 1024;
		NWmatrices->E = smalloc(NWmatrices->NW_s);
		NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
		NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
		NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
		
		aligned = smalloc(sizeof(struct aln));
		gap_align = smalloc(sizeof(struct aln));
		aligned->t = smalloc((delta + 1) << 1);
		aligned->s = smalloc((delta + 1) << 1);
		aligned->q = smalloc((delta + 1) << 1);
		gap_align->t = smalloc((delta + 1) << 1);
		gap_align->s = smalloc((delta + 1) << 1);
		gap_align->q = smalloc((delta + 1) << 1);
		
		/* move it to the thread */
		thread = smalloc(sizeof(struct assemble_thread));
		thread->num = i;
		
		thread->template = -2;
		thread->file_count = fileCount;
		thread->files = template_fragments;
		thread->frag_out = frag_out;
		thread->aligned_assem = aligned_assem;
		thread->aligned = aligned;
		thread->gap_align = gap_align;
		thread->NWmatrices = NWmatrices;
		thread->matrix = matrix;
		thread->qseq = setQseqs(qseq->size);
		thread->header = setQseqs(header->size);
		thread->points = seedPoint_init(delta);
		thread->points->len = 0;
		thread->spin = (sparse < 0) ? 10 : 100;
		
		thread->next = threads;
		threads = thread;
		
		/* start thread */
		if((errno = pthread_create(&thread->id, NULL, assembly_KMA_Ptr, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", i);
			threads = thread->next;
			free(thread);
			i = thread_num;
		} else {
			++i;
		}
	}
	
	/* start main thread */
	NWmatrices = smalloc(sizeof(struct NWmat));
	NWmatrices->NW_s = 1024 * 1024;
	NWmatrices->NW_q = 1024;
	NWmatrices->E = smalloc(NWmatrices->NW_s);
	NWmatrices->D[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->P[0] = smalloc((NWmatrices->NW_q << 1) * sizeof(int));
	NWmatrices->D[1] = NWmatrices->D[0] + NWmatrices->NW_q;
	NWmatrices->P[1] = NWmatrices->P[0] + NWmatrices->NW_q;
	
	aligned = smalloc(sizeof(struct aln));
	gap_align = smalloc(sizeof(struct aln));
	aligned->t = smalloc((delta + 1) << 1);
	aligned->s = smalloc((delta + 1) << 1);
	aligned->q = smalloc((delta + 1) << 1);
	gap_align->t = smalloc((delta + 1) << 1);
	gap_align->s = smalloc((delta + 1) << 1);
	gap_align->q = smalloc((delta + 1) << 1);
	
	/* move it to the thread */
	thread = smalloc(sizeof(struct assemble_thread));
	thread->num = 0;
	thread->template = 0;
	thread->file_count = fileCount;
	thread->files = template_fragments;
	thread->frag_out = frag_out;
	thread->aligned_assem = aligned_assem;
	thread->aligned = aligned;
	thread->gap_align = gap_align;
	thread->NWmatrices = NWmatrices;
	thread->matrix = matrix;
	thread->qseq = qseq;
	thread->header = header;
	thread->points = points;
	thread->points->len = 0;
	thread->next = 0;
	thread->spin = (sparse < 0) ? 10 : 100;
	
	/* Do local assemblies of fragments mapping to the same template */
	depth = 0;
	q_id = 0;
	cover = 0;
	q_cover = 0;
	seq_seeker = 0;
	index_seeker = 0;
	progress = 0;
	counter = 0;
	if(progress) {
		fprintf(stderr, "# Progress:\t%3d%%\r", 0);
		fflush(stderr);
	}
	
	/* get index, seq and names on the fly */
	bias = *dbBiases++;
	templatefilename = *templatefilenames++;
	getIndexName(&template_names, &templates_index, 0, bias, *dbBiases, templatefilename);
	file_len = strlen(templatefilename);
	strcat(templatefilename, ".seq.b");
	seq_in = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	
	strcat(templatefilename, ".index.b");
	index_in = fopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	if(!index_in) {
		alignLoadPtr = &alignLoad_fly_build;
		destroyPtr = &alignClean;
		index_in = 0;
		if(kmersize < 4 || 32 < kmersize) {
			kmersize = 16;
		}
	} else if(shm & 8) {
		templates_index[0] = alignLoad_shm_initial(templatefilename, file_len, seq_in, index_in);
		alignLoadPtr = &alignLoad_fly_shm;
		destroyPtr = &alignClean_shm;
	} else {
		alignLoadPtr = &alignLoad_fly_mem;
		destroyPtr = &alignClean;
		fread(&kmersize, sizeof(int), 1, index_in);
	}
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
	shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
	r_shifter = (kmersize << 1) - 2;
	
	for(template = 1; template < DB_size; ++template) {
		if(template == bias) {
			/* swap indexes */
			templatefilename = *templatefilenames++;
			getIndexName(&template_names, &templates_index, bias - dbBiases[-1], *dbBiases, dbBiases[1], templatefilename);
			bias = *dbBiases++;
			file_len = strlen(templatefilename);
			strcat(templatefilename, ".seq.b");
			seq_in = sfopen(templatefilename, "rb");
			templatefilename[file_len] = 0;
			
			strcat(templatefilename, ".index.b");
			index_in = fopen(templatefilename, "rb");
			templatefilename[file_len] = 0;
			if(!index_in) {
				alignLoadPtr = &alignLoad_fly_build;
				destroyPtr = &alignClean;
				index_in = 0;
				if(kmersize < 4 || 32 < kmersize) {
					kmersize = 16;
				}
			} else if(shm & 8) {
				templates_index[0] = alignLoad_shm_initial(templatefilename, file_len, seq_in, index_in);
				alignLoadPtr = &alignLoad_fly_shm;
				destroyPtr = &alignClean_shm;
			} else {
				alignLoadPtr = &alignLoad_fly_mem;
				destroyPtr = &alignClean;
				fread(&kmersize, sizeof(int), 1, index_in);
			}
			mask = 0;
			mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
			shifter = sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1);
			r_shifter = (kmersize << 1) - 2;
			seq_seeker = 0;
			index_seeker = 0;
		} else if(w_scores[template] > 0) {
			
			if(progress) {
				counter += w_scores[template];
				fprintf(stderr, "# Progress:\t%3lu%%\r", 100 * counter / Nhits);
				fflush(stderr);
			}
			
			/* make p_value to see whether assembly is feasable */
			read_score = w_scores[template];
			t_len = template_lengths[template];
			expected = t_len;
			expected /= (template_tot_ulen - t_len);
			expected *= (Nhits - read_score);
			if(0 < expected) {
				q_value = read_score - expected;
				q_value /= (expected + read_score);
				q_value *= (read_score - expected);
			} else {
				q_value = read_score;
			}
			p_value  = p_chisqr(q_value);
			
			if(cmp((p_value <= evalue && read_score > expected), (read_score >= scoreT * t_len))) {
				/* load DB */
				if(index_in) {
					index_seeker *= sizeof(int);
					fseek(index_in, index_seeker, SEEK_CUR);
					index_seeker = 0;
				}
				seq_seeker *= sizeof(long unsigned);
				fseek(seq_in, seq_seeker, SEEK_CUR);
				seq_seeker = 0;
				templates_index[template] = alignLoadPtr(seq_in, index_in, template_lengths[template], 0, 0);
				
				/* Do assembly */
				//status |= assemblyPtr(aligned_assem, template, template_fragments, fileCount, frag_out, aligned, gap_align, qseq, header, matrix, points, NWmatrices);
				thread->template = template;
				assembly_KMA_Ptr(thread);
				
				/* Depth, ID and coverage */
				if(aligned_assem->cover > 0) {
					coverScore = aligned_assem->cover;
					depth = aligned_assem->depth;
					depth /= t_len;
					id = 100.0 * coverScore / t_len;
					aln_len = aligned_assem->aln_len;
					q_id = 100.0 * coverScore / aln_len;
					cover = 100.0 * aln_len / t_len;
					q_cover = 100.0 * t_len / aln_len;
				} else {
					id = 0;
				}
				
				if(ID_t <= id && 0 < id) {
					/* Output result */
					fprintf(res_out, "%-12s\t%8ld\t%8u\t%8d\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%4.1e\n",
						template_names[template], read_score, (unsigned) expected, t_len, id, cover, q_id, q_cover, (double) depth, (double) q_value, p_value);
					printConsensus(aligned_assem, template_names[template], alignment_out, consensus_out);
					/* print matrix */
					if(matrix_out) {
						updateMatrix(matrix_out, template_names[template], templates_index[template]->seq, matrix, t_len);
					}
					if(extendedFeatures) {
						getExtendedFeatures(template_names[template], matrix, templates_index[template]->seq, t_len, aligned_assem, fragmentCounts[template], readCounts[template], extendedFeatures_out);
					}
					if(vcf) {
						updateVcf(template_names[template], templates_index[template]->seq, t_len, matrix, vcf, vcf_out);
					}
				}
				/* destroy this DB index */
				destroyPtr(template);
			} else {
				if(index_in) {
					index_seeker += (template_lengths[template] << 1);
				}
				seq_seeker += ((template_lengths[template] >> 5) + 1);
			}
		} else {
			if(index_in) {
				index_seeker += (template_lengths[template] << 1);
			}
			seq_seeker += ((template_lengths[template] >> 5) + 1);
		}
	}
	template_names += bias;
	templates_index += bias;
	
	if(progress) {
		fprintf(stderr, "\n");
	}
	
	/* join threads */
	thread->template = -1;
	assembly_KMA_Ptr(thread);
	for(thread = threads; thread != 0; thread = thread->next) {
		/* join thread */
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
	}
	
	/* Close files */
	if(index_in) {
		fclose(index_in);
	}
	fclose(seq_in);
	fclose(res_out);
	fclose(alignment_out);
	fclose(consensus_out);
	destroyGzFileBuff(frag_out);
	if(matrix_out) {
		destroyGzFileBuff(matrix_out);
	}
	if(extendedFeatures) {
		fclose(extendedFeatures_out);
	}
	if(vcf) {
		destroyGzFileBuff(vcf_out);
	}
	
	t1 = clock();
	fprintf(stderr, "# Total time used for local assembly: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	
	return status;
}

char * strjoin(char **strings, int len) {
	
	int i, new_len, escape;
	char *newStr, *stringPtr;
	
	new_len = len + 16;
	escape = 0;
	for(i = 0; i < len; ++i) {
		if(*strings[i] == '-') {
			escape = 0;
		} else if(escape) {
			new_len += 2;
		}
		new_len += strlen(strings[i]);
		if(strncmp(strings[i], "-i", 2) == 0) {
			escape = 1;
		}
	}
	
	newStr = malloc(new_len);
	if(!newStr) {
		ERROR();
	}
	
	*newStr = 0;
	escape = 0;
	stringPtr = newStr;
	for(i = 0; i < len; ++i) {
		if(*strings[i] == '-') {
			escape = 0;
		}
		
		if(escape) {
			*stringPtr = '\"';
			++stringPtr;
		}
		new_len = strlen(strings[i]);
		strcpy(stringPtr, strings[i]);
		stringPtr += new_len;
		if(escape) {
			*stringPtr = '\"';
			++stringPtr;
		}
		*stringPtr = ' ';
		++stringPtr;
		
		if(*strings[i] == '-' && (strings[i][1] == 'i' || strings[i][1] == 'o')) {
			escape = 1;
		}
	}
	
	return newStr;
}

int sparseCheck(char *filename) {
	
	unsigned prefix_len, file_len;
	long unsigned prefix;
	FILE *file;
	
	file_len = strlen(filename);
	strcat(filename, ".comp.b");
	file = sfopen(filename, "rb");
	filename[file_len] = 0;
	
	/* get prefix */
	fseek(file, 2 * sizeof(unsigned), SEEK_SET);
	fread(&prefix_len, sizeof(unsigned), 1, file);
	fread(&prefix, sizeof(long unsigned), 1, file);
	fclose(file);
	
	if(prefix_len != 0 || prefix != 0) {
		alignLoadPtr = &alignLoad_fly_mem;
		ankerPtr = &ankerAndClean_MEM;
		return 1;
	}
	
	return 0;
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# KMA-%d.%d.%d mapps raw reads to a template database.\n", version[0], version[1], version[2]);
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-i\t\tInput file name(s)\t\tSTDIN\n");
	fprintf(helpOut, "#\t-ipe\t\tInput paired end file name(s)\n");
	fprintf(helpOut, "#\t-int\t\tInput interleaved file name(s)\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t%s\n", "DB defined");
	fprintf(helpOut, "#\t-e\t\tevalue\t\t\t\t%1.2f\n", evalue);
	fprintf(helpOut, "#\t-ConClave\tConClave version\t\t%d\n", 1);
	fprintf(helpOut, "#\t-mem_mode\tUse kmers to choose best\n#\t\t\ttemplate, and save memory\tFalse\n");
	fprintf(helpOut, "#\t-ex_mode\tSearh kmers exhaustively\tFalse\n");
	fprintf(helpOut, "#\t-ef\t\tPrint additional features\tFalse\n");
	fprintf(helpOut, "#\t-vcf\t\tMake vcf file, 2 to apply FT\tFalse/0\n");
	fprintf(helpOut, "#\t-deCon\t\tRemove contamination\t\tFalse\n");
	fprintf(helpOut, "#\t-dense\t\tDo not allow insertions\n#\t\t\tin assembly\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ref_fsa\tConsensus sequnce will\n#\t\t\thave \"n\" instead of gaps\tFalse\n");
	fprintf(helpOut, "#\t-matrix\t\tPrint assembly matrix\t\tFalse\n");
	fprintf(helpOut, "#\t-a\t\tPrint all best mappings\t\tFalse\n");
	fprintf(helpOut, "#\t-mp\t\tMinimum phred score\t\t30\n");
	fprintf(helpOut, "#\t-5p\t\tCut a constant number of\n#\t\t\tnucleotides from the 5 prime.\t0\n");
	fprintf(helpOut, "#\t-Sparse\t\tOnly count kmers\t\tFalse\n");
	fprintf(helpOut, "#\t-Mt1\t\tMap only to \"num\" template.\t0 / False\n");
	fprintf(helpOut, "#\t-ID\t\tMinimum ID\t\t\t%3.1f%%\n", ID_t);
	fprintf(helpOut, "#\t-ss\t\tSparse sorting (q,c,d)\t\tq\n");
	fprintf(helpOut, "#\t-pm\t\tPairing method (p,u,f)\t\tu\n");
	fprintf(helpOut, "#\t-fpm\t\tFine Pairing method (p,u,f)\tu\n");
	fprintf(helpOut, "#\t-apm\t\tSets both pm and fpm\t\tu\n");
	fprintf(helpOut, "#\t-shm\t\tUse shared DB made by kma_shm\t%d (lvl)\n", shm);
	fprintf(helpOut, "#\t-swap\t\tSwap DB to disk\t\t\t%d (lvl)\n", diskDB);
	fprintf(helpOut, "#\t-1t1\t\tSkip HMM\t\t\tFalse\n");
	fprintf(helpOut, "#\t-ck\t\tCount kmers instead of\n#\t\t\tpseudo alignment\t\tFalse\n");
	fprintf(helpOut, "#\t-ca\t\tMake circular alignments\tFalse\n");
	fprintf(helpOut, "#\t-boot\t\tBootstrap sequence\t\tFalse\n");
	fprintf(helpOut, "#\t-bc\t\tBase calls should be\n#\t\t\tsignificantly overrepresented.\t[True]\n");
	fprintf(helpOut, "#\t-bc90\t\tBase calls should be both\n#\t\t\tsignificantly overrepresented,\n#\t\t\tand have 90%% agreement.\t\tFalse\n");
	fprintf(helpOut, "#\t-bcNano\t\tCall bases at suspicious\n#\t\t\tdeletions, made for nanopore.\tFalse\n");
	fprintf(helpOut, "#\t-and\t\tBoth mrs and p_value thresholds\n#\t\t\thas to reached to in order to\n#\t\t\treport a template hit.\t\tor\n");
	fprintf(helpOut, "#\t-mq\t\tMinimum mapping quality\t\t%u\n", mq);
	fprintf(helpOut, "#\t-mrs\t\tMinimum alignment score,\n#\t\t\tnormalized to alignment length\t%1.2f\n", scoreT);
	fprintf(helpOut, "#\t-reward\t\tScore for match\t\t\t%d\n", M);
	fprintf(helpOut, "#\t-penalty\tPenalty for mismatch\t\t%d\n", MM);
	fprintf(helpOut, "#\t-gapopen\tPenalty for gap opening\t\t%d\n", W1);
	fprintf(helpOut, "#\t-gapextend\tPenalty for gap extension\t%d\n", U);
	fprintf(helpOut, "#\t-per\t\tReward for pairing reads\t%d\n", PE);
	fprintf(helpOut, "#\t-cge\t\tSet CGE penalties and rewards\tFalse\n");
	fprintf(helpOut, "#\t-t\t\tNumber of threads\t\t%d\n", thread_num);
	fprintf(helpOut, "#\t-v\t\tVersion\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int i, j, args, exe_len, minPhred, fiveClip, sparse_run, mem_mode, Mt1;
	int step1, step2, fileCounter, fileCounter_PE, fileCounter_INT, status;
	int ConClave, extendedFeatures, vcf, targetNum, size, escape, spltDB;
	long unsigned totFrags;
	char *exeBasic, *outputfilename, *templatefilename, **templatefilenames;
	char **inputfiles, **inputfiles_PE, **inputfiles_INT, *to2Bit, Date[11];
	char ss;
	FILE *templatefile;
	time_t t0, t1;
	struct tm *tm;
	
	if(sizeof(long unsigned) != 8) {
		fprintf(stderr, "Need a 64-bit system.\n");
		exit(3);
	}
	
	/* SET DEFAULTS */
	ConClave = 1;
	totFrags = 0;
	vcf = 0;
	targetNum = 0;
	spltDB = 0;
	extendedFeatures = 0;
	status = 0;
	countK = 0;
	assembly_KMA_Ptr = &assemble_KMA_threaded;
	cmp = &cmp_or;
	minPhred = 30;
	fiveClip = 0;
	sparse_run = 0;
	fileCounter = 0;
	fileCounter_PE = 0;
	fileCounter_INT = 0;
	outputfilename = 0;
	templatefilename = 0;
	print_matrix = 0;
	print_all = 0;
	ref_fsa = 0;
	contamination = -1;
	kmersize = 0;
	evalue = 0.05;
	delta = 2048;
	exhaustive = 0;
	shm = 0;
	diskDB = 0;
	mq = 0;
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
	PE = 7;
	thread_num = 1;
	kmerScan = &save_kmers_HMM;
	get_kmers_for_pair_ptr = &get_kmers_for_pair;
	save_kmers_pair = &save_kmers_unionPair;
	alnFragsPE = &alnFragsUnionPE;
	printPairPtr = &printPair;
	printPtr = &print_ankers;
	printFsa_pair_ptr = &printFsa_pair;
	deConPrintPtr = printPtr;
	ankerPtr = &ankerAndClean;
	alignLoadPtr = &alignLoad_fly;
	destroyPtr = &alignClean;
	printFsa_ptr = &printFsa;
	inputfiles_PE = 0;
	inputfiles_INT = 0;
	inputfiles = 0;
	templatefilenames = 0;
	Mt1 = 0;
	significantBase = &significantNuc; //-bc
	baseCall = &baseCaller;
	chainSeedsPtr = &chainSeeds;
	inputfiles = 0;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			if(++args < argc) {
				templatefilename = malloc(strlen(argv[args]) + 64);
				if(!templatefilename) {
					ERROR();
				}
				strcpy(templatefilename, argv[args]);
				++targetNum;
				templatefilenames = realloc(templatefilenames, targetNum * sizeof(char *));
				templatefilenames[targetNum - 1] = templatefilename;
			}
			while(++args < argc && *argv[args] != '-') {
				++targetNum;
				templatefilenames = realloc(templatefilenames, targetNum * sizeof(char *));
				templatefilenames[targetNum - 1] = argv[args];
			}
			--args;
		} else if(strcmp(argv[args], "-i") == 0) {
			++args;
			status = fileCounter;
			for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
				++fileCounter;
			}
			if(fileCounter == 0) {
				fprintf(stderr, "No files were specified.\n");
				exit(3);
			} else {
				inputfiles = realloc(inputfiles, fileCounter * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
			}
			
			for(i = status; i < fileCounter; ++i, ++args) {
				inputfiles[i] = strdup(argv[args]);
				if(!inputfiles[i]) {
					ERROR();
				}
			}
			--args;
		} else if(strcmp(argv[args], "-ipe") == 0) {
			++args;
			status = fileCounter_PE;
			for(i = args; i < argc && strncmp(argv[i], "-", 1) != 0; ++i) {
				++fileCounter_PE;
			}
			if(fileCounter_PE % 2) {
				fprintf(stderr, "Uneven number of paired end files.\n");
				exit(3);
			} else if(fileCounter_PE == 0) {
				fprintf(stderr, "No paired end files were specified.\n");
				exit(3);
			} else {
				inputfiles_PE = realloc(inputfiles_PE, fileCounter_PE * sizeof(char *));
				if(!inputfiles_PE) {
					ERROR();
				}
			}
			
			for(i = status; i < fileCounter_PE; ++i, ++args) {
				inputfiles_PE[i] = strdup(argv[args]);
				if(!inputfiles_PE[i]) {
					ERROR();
				}
			}
			--args;
		} else if(strcmp(argv[args], "-int") == 0) {
			++args;
			status = fileCounter_INT;
			for(i = args; i < argc && (strncmp(argv[i], "-", 1) != 0 || strcmp(argv[i], "--") == 0); ++i) {
				++fileCounter_INT;
			}
			if(fileCounter_INT == 0) {
				fprintf(stderr, "No interleaved files were specified.\n");
				exit(3);
			}
			inputfiles_INT = realloc(inputfiles_INT, fileCounter_INT * sizeof(char *));
			if(!inputfiles_INT) {
				ERROR();
			}
			for(i = status; i < fileCounter_INT; ++i, ++args) {
				inputfiles_INT[i] = strdup(argv[args]);
				if(!inputfiles_INT[i]) {
					ERROR();
				}
			}
			--args;
		} else if(strcmp(argv[args], "-pm") == 0) {
			++args;
			if(args < argc) {
				if(*(argv[args]) == 'p') {
					save_kmers_pair = &save_kmers_penaltyPair;
				} else if(*(argv[args]) == 'u') {
					save_kmers_pair = &save_kmers_unionPair;
				} else if(*(argv[args]) == 'f') {
					save_kmers_pair = &save_kmers_forcePair;
				} else {
					fprintf(stderr, "Invalid argument at pairing method: \"-pm\"\n");
					fprintf(stderr, "Options are:\n");
					fprintf(stderr, "p:\tReward for pairing.\n");
					fprintf(stderr, "u:\tUnion of best hits.\n");
					fprintf(stderr, "f:\tForce paring.\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-fpm") == 0) {
			++args;
			if(args < argc) {
				if(*(argv[args]) == 'p') {
					alnFragsPE = &alnFragsPenaltyPE;
				} else if(*(argv[args]) == 'u') {
					alnFragsPE = &alnFragsUnionPE;
				} else if(*(argv[args]) == 'f') {
					alnFragsPE = &alnFragsForcePE;
				} else {
					fprintf(stderr, "Invalid argument at fine pairing method: \"-fpm\"\n");
					fprintf(stderr, "Options are:\n");
					fprintf(stderr, "p:\tReward for pairing.\n");
					fprintf(stderr, "u:\tUnion of best hits.\n");
					fprintf(stderr, "f:\tForce paring.\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-apm") == 0) {
			++args;
			if(args < argc) {
				if(*(argv[args]) == 'p') {
					alnFragsPE = &alnFragsPenaltyPE;
					save_kmers_pair = &save_kmers_penaltyPair;
				} else if(*(argv[args]) == 'u') {
					alnFragsPE = &alnFragsUnionPE;
					save_kmers_pair = &save_kmers_unionPair;
				} else if(*(argv[args]) == 'f') {
					alnFragsPE = &alnFragsForcePE;
					save_kmers_pair = &save_kmers_forcePair;
				} else {
					fprintf(stderr, "Invalid argument at fine pairing method: \"-fpm\"\n");
					fprintf(stderr, "Options are:\n");
					fprintf(stderr, "p:\tReward for pairing.\n");
					fprintf(stderr, "u:\tUnion of best hits.\n");
					fprintf(stderr, "f:\tForce paring.\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-delta") == 0) {
			++args;
			if(args < argc) {
				delta = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, " Invalid delta specified.\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-ConClave") == 0) {
			++args;
			if(args < argc) {
				ConClave = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0 || ConClave < 0 || 2 < ConClave) {
					fprintf(stderr, " Invalid ConClave version specified.\n");
					exit(4);
				}
			}
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
			deCon = 1;
			deConPrintPtr = &deConPrint;
			printPairPtr = &deConPrintPair;
		} else if(strcmp(argv[args], "-shm") == 0) {
			++args;
			if(args < argc && argv[args][0] != '-') {
				shm = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid shm-lvl specified.\n");
					exit(4);
				}
			} else {
				--args;
				shm = 3;
			}
		} else if(strcmp(argv[args], "-t") == 0) {
			++args;
			if(args < argc && argv[args][0] != '-') {
				thread_num = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid number of threads specified.\n");
					exit(4);
				}
			} else {
				--args;
			}
			if(thread_num < 1) {
				thread_num = 1;
			}
		} else if(strcmp(argv[args], "-swap") == 0) {
			++args;
			if(args < argc && argv[args][0] != '-') {
				/*
				diskDB = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-swap\".\n");
					exit(4);
				}
				*/
				diskDB = 0;
			} else {
				--args;
				/*diskDB = 1;*/
				diskDB = 0;
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
			}
		} else if(strcmp(argv[args], "-mp") == 0) {
			++args;
			if(args < argc) {
				minPhred = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum phred score parsed\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-mq") == 0) {
			++args;
			if(args < argc) {
				mq = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "# Invalid minimum mapping quality parsed\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-5p") == 0) {
			++args;
			if(args < argc) {
				fiveClip = strtoul(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-5p\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-dense") == 0) {
			assembly_KMA_Ptr = &assemble_KMA_dense_threaded;
		} else if(strcmp(argv[args], "-matrix") == 0) {
			print_matrix = 1;
		} else if(strcmp(argv[args], "-a") == 0) {
			print_all = 1;
			mem_mode = 1;
			alignLoadPtr = &alignLoad_fly_mem;
			ankerPtr = &ankerAndClean_MEM;
		} else if(strcmp(argv[args], "-ref_fsa") == 0) {
			ref_fsa = 1;
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
		} else if(strcmp(argv[args], "-1t1") == 0) {
			kmerScan = &save_kmers;
			one2one = 1;
		} else if(strcmp(argv[args], "-ck") == 0) {
			countK = 1;
		} else if(strcmp(argv[args], "-ca") == 0) {
			chainSeedsPtr = &chainSeeds_circular;
		} else if(strcmp(argv[args], "-ss") == 0) {
			++args;
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
			++args;
			if(args < argc) {
				evalue = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-e\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-bc") == 0) {
			if(++args < argc && argv[args][0] != '-') {
				significantBase = &significantAndSupport;
				support = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0 || 1 < support) {
					fprintf(stderr, "Invalid argument at \"-bc\".\n");
					exit(4);
				}
			} else {
				--args;
				significantBase = &significantNuc;
			}
		} else if(strcmp(argv[args], "-bc90") == 0) {
			significantBase = &significantAnd90Nuc;
		} else if(strcmp(argv[args], "-bcNano") == 0) {
			if(significantBase == &significantNuc) {
				significantBase = &significantAnd90Nuc;
			}
			baseCall = &nanoCaller;
		} else if(strcmp(argv[args], "-ID") == 0) {
			++args;
			if(args < argc) {
				ID_t = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-ID\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-mrs") == 0) {
			++args;
			if(args < argc) {
				scoreT = strtod(argv[args], &exeBasic);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-mrs\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-reward") == 0) {
			++args;
			if(args < argc) {
				M = strtol(argv[args], &exeBasic, 10);
				M = abs(M);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-reward\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-penalty") == 0) {
			++args;
			if(args < argc) {
				MM = strtol(argv[args], &exeBasic, 10);
				MM = MIN(-MM, MM);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-penalty\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-gapopen") == 0) {
			++args;
			if(args < argc) {
				W1 = strtol(argv[args], &exeBasic, 10);
				W1 = MIN(-W1, W1);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-gapopen\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-gapextend") == 0) {
			++args;
			if(args < argc) {
				U = strtol(argv[args], &exeBasic, 10);
				U = MIN(-U, U);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-gapextend\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-per") == 0) {
			++args;
			if(args < argc) {
				PE = strtol(argv[args], &exeBasic, 10);
				PE = abs(PE);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-per\".\n");
					exit(4);
				}
			}
		} else if(strcmp(argv[args], "-and") == 0) {
			cmp = &cmp_and;
		} else if(strcmp(argv[args], "-boot") == 0) {
			printFsa_ptr = &bootFsa;
		} else if(strcmp(argv[args], "-Mt1") == 0) {
			++args;
			if(args < argc) {
				Mt1 = strtol(argv[args], &exeBasic, 10);
				if(*exeBasic != 0) {
					fprintf(stderr, "Invalid argument at \"-Mt1\".\n");
					exit(4);
				}
			}
			if(Mt1 < 1) {
				fprintf(stderr, "Invalid template specified at \"-Mt1\"\n");
				exit(3);
			}
			printFsa_ptr = &printFsaMt1;
			printFsa_pair_ptr = &printFsa_pairMt1;
		} else if(strcmp(argv[args], "-ef") == 0) {
			extendedFeatures = 1;
		} else if(strcmp(argv[args], "-vcf") == 0) {
			vcf = 1;
			if(++args < argc) {
				if(argv[args][0] != '-') {
					vcf = strtol(argv[args], &exeBasic, 10);
					if(*exeBasic != 0) {
						fprintf(stderr, "Invalid argument at \"-vcf\".\n");
						exit(4);
					}
				} else {
					--args;
				}
			}
		} else if(strcmp(argv[args], "-cge") == 0) {
			scoreT = 0.75;
			M = 1;
			MM = -3;
			W1 = -5;
			U = -1;
			PE = 17;
		} else if(strcmp(argv[args], "-spltDB") == 0) {
			spltDB = 1;
		} else if(strcmp(argv[args], "-v") == 0) {
			fprintf(stdout, "KMA-%d.%d.%d\n", version[0], version[1], version[2]);
			exit(0);
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, " Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, " Printing help message:\n");
			helpMessage(1);
		}
		++args;
	}
	
	if(spltDB || targetNum != 1) {
		printPtr = &print_ankers_spltDB;
		if(deConPrintPtr != &deConPrint) {
			deConPrintPtr = printPtr;
		}
		kmerScan = &save_kmers;
		one2one = 1;
	}
	
	if(one2one && countK) {
		kmerScan = &save_kmers_count;
		get_kmers_for_pair_ptr = &get_kmers_for_pair_count;
	}
	
	if(ref_fsa) {
		if(baseCall == nanoCaller) {
			baseCall = &refNanoCaller;
		} else {
			baseCall = &refCaller;
		}
	}
	
	if(outputfilename == 0 || templatefilename == 0) {
		fprintf(stderr, " Too few arguments handed\n");
		fprintf(stderr, " Printing help message:\n");
		helpMessage(1);
	}
	
	if(fileCounter == 0 && fileCounter_PE == 0 && fileCounter_INT == 0) {
		inputfiles = malloc(sizeof(char*));
		if(!inputfiles) {
			ERROR();
		}
		inputfiles[0] = strdup("--");
		if(!inputfiles[0]) {
			ERROR();
		}
		fileCounter = 1;
	}
	status = 0;
	
	/* set scoring matrix */
	for(i = 0; i < 4; ++i) {
		for(j = 0; j < 4; ++j) {
			d[i][j] = MM;
		}
		d[i][i] = M;
	}
	for(i = 0; i < 5; ++i) {
		d[4][i] = U;
		d[i][4] = U;
	}
	d[4][4] = 0;
	
	if(spltDB && targetNum != 1) {
		/* allocate space for commands */
		escape = 0;
		size = argc + strlen(outputfilename) + 32;
		for(args = 0; args < argc; ++args) {
			if(*argv[args] == '-') {
				escape = 0;
			} else if(escape) {
				size += 2;
			}
			size += strlen(argv[i]);
			if(strncmp(argv[i], "-i", 2) == 0) {
				escape = 1;
			}
		}
		exeBasic = smalloc(size);
		
		fprintf(stderr, "# Map\n");
		for(i = 0; i < targetNum; ++i) {
			to2Bit = exeBasic;
			*to2Bit = 0;
			args = -1;
			while(++args < argc) {
				if(strcmp(argv[args], "-t_db") == 0) {
					escape = 1;
					while(escape && ++args < argc) {
						if(*argv[args] == '-') {
							escape = 0;
						}
					}
					--args;
				} else {
					if(*argv[args] == '-') {
						escape = 0;
					}
					
					if(escape) {
						*to2Bit = '\"';
						++to2Bit;
					}
					
					exe_len = strlen(argv[args]);
					strcpy(to2Bit, argv[args]);
					to2Bit += exe_len;
					
					if(escape) {
						*to2Bit = '\"';
						++to2Bit;
					}
					*to2Bit = ' ';
					++to2Bit;
					
					if(strncmp(argv[args], "-i", 2) == 0) {
						escape = 1;
					}
				}
			}
			fprintf(stdout, "%s-t_db %s -step2 > %s.%d\n", exeBasic, templatefilenames[i], outputfilename, i);
		}
		
		fprintf(stderr, "# Reduce:\n");
		to2Bit = exeBasic;
		*to2Bit = 0;
		args = -1;
		while(++args < argc) {
			if(strcmp(argv[args], "-spltDB") != 0) {
				if(*argv[args] == '-') {
					escape = 0;
				}
				
				if(escape) {
					*to2Bit = '\"';
					++to2Bit;
				}
				
				exe_len = strlen(argv[args]);
				strcpy(to2Bit, argv[args]);
				to2Bit += exe_len;
				
				if(escape) {
					*to2Bit = '\"';
					++to2Bit;
				}
				*to2Bit = ' ';
				++to2Bit;
				
				if(strncmp(argv[args], "-i", 2) == 0) {
					escape = 1;
				}
			}
		}
		fprintf(stdout, "%s\n", exeBasic);
		
		return 0;
	} else {
		freopen(NULL, "wb", stdout);
		setvbuf(stdout, NULL, _IOFBF, CHUNK);
	}
	
	if(step1) {
		t0 = clock();
		/* set to2Bit conversion */
		to2Bit = malloc(384); /* 128 * 3 = 384 -> OS independent */
		if(!to2Bit) {
			ERROR();
		}
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
		
		if(sparse_run) {
			templates = malloc(sizeof(struct hashMapKMA));
			if(!templates) {
				ERROR();
			}
			exe_len = strlen(templatefilename);
			if(deCon) {
				strcat(templatefilename, ".decon.comp.b");
			} else {
				strcat(templatefilename, ".comp.b");
			}
			templatefile = sfopen(templatefilename, "rb" );
			getPrefix(templates, templatefile);
			fclose(templatefile);
			templatefilename[exe_len] = 0;
			mask = 0;
			if(templates->prefix_len) {
				mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (templates->prefix_len << 1));
			} else {
				mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (kmersize << 1));
			}
			
			/* merge reads */
			if(fileCounter_PE > 0) {
				inputfiles = realloc(inputfiles, (fileCounter + fileCounter_PE) * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
				for(i = 0; i < fileCounter_PE; ++i, ++fileCounter) {
					inputfiles[fileCounter] = inputfiles_PE[i];
				}
				free(inputfiles_PE);
				fprintf(stderr, "Paired end information is not considered in Sparse mode.\n");
			}
			if(fileCounter_INT > 0) {
				inputfiles = realloc(inputfiles, (fileCounter + fileCounter_INT) * sizeof(char *));
				if(!inputfiles) {
					ERROR();
				}
				for(i = 0; i < fileCounter_INT; ++i, ++fileCounter) {
					inputfiles[fileCounter] = inputfiles_INT[i];
				}
				free(inputfiles_INT);
				fprintf(stderr, "Interleaved information is not considered in Sparse mode.\n");
			}
			
			run_input_sparse(inputfiles, fileCounter, minPhred, fiveClip, to2Bit);
		} else {
			
			if(Mt1) {
				template_lengths = malloc(7 * sizeof(int));
				if(!template_lengths) {
					ERROR();
				}
				strcat(templatefilename, ".length.b");
				templatefile = sfopen(templatefilename, "rb");
				/* DB_size s_lengths u_length lengths */
				template_lengths[0] = 0;
				template_lengths[2] = 1;
				template_lengths[3] = 0;
				template_lengths[4] = 0;
				fread(&DB_size, sizeof(int), 1, templatefile);
				fseek(templatefile, Mt1 * sizeof(int), SEEK_CUR);
				fread(template_lengths + 5, sizeof(int), 1, templatefile);
				/*
				fseek(templatefile, (2 * DB_size + Mt1) * sizeof(int), SEEK_CUR);
				if(fread(template_lengths + 5, sizeof(int), 1, templatefile) == 0) {
					fseek(templatefile, (Mt1 + 1) * sizeof(int), SEEK_SET);
				}
				*/
				fclose(templatefile);
				DB_size = Mt1;
			}
			kmersize = 16;
			
			/* SE */
			if(fileCounter > 0) {
				totFrags += run_input(inputfiles, fileCounter, minPhred, fiveClip, to2Bit);
			}
			
			/* PE */
			if(fileCounter_PE > 0) {
				totFrags += run_input_PE(inputfiles_PE, fileCounter_PE, minPhred, fiveClip, to2Bit);
			}
			
			/* INT */
			if(fileCounter_INT > 0) {
				totFrags += run_input_INT(inputfiles_INT, fileCounter_INT, minPhred, fiveClip, to2Bit);
			}
			
			if(Mt1) {
				Mt1 = -1;
				sfwrite(&Mt1, sizeof(int), 1, stdout);
			}
			
			if(extendedFeatures && targetNum == 1) {
				strcat(outputfilename, ".mapstat");
				templatefile = sfopen(outputfilename, "wb");
				fprintf(templatefile, "## method\tKMA\n");
				fprintf(templatefile, "## version\t%d.%d.%d\n", version[0], version[1], version[2]);
				fprintf(templatefile, "## database %s\n", noFolder(templatefilename));
				fprintf(templatefile, "## fragmentCount\t%lu\n", totFrags);
				time(&t1);
				tm = localtime(&t1);
				strftime(Date, sizeof(Date), "%Y-%m-%d", tm);
				fprintf(templatefile, "## date\t%s\n", Date);
				//fprintf(templatefile, "## date\t%s", ctime(&t1));
				fprintf(templatefile, "## command\t%s", *argv);
				argc -= step1;
				argc -= step2;
				for(args = 1; args < argc; ++args) {
					fprintf(templatefile, " %s", argv[args]);
				}
				fprintf(templatefile, "\n");
				fclose(templatefile);
			}
		}
		fflush(stdout);
		t1 = clock();
		fprintf(stderr, "#\n# Total time used for converting query: %.2f s.\n#\n", difftime(t1, t0) / 1000000);
	} else if(Mt1) {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-step1");
		
		runKMA_Mt1(templatefilename, outputfilename, exeBasic, Mt1, vcf);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	} else if(step2) {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-step1");
		status = save_kmers_batch(templatefilename, exeBasic);
		fflush(stdout);
	} else if(sparse_run) {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-step1");
		
		status = save_kmers_sparse_batch(templatefilename, outputfilename, exeBasic, ss);
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	} else {
		exeBasic = strjoin(argv, argc);
		strcat(exeBasic, "-step2");
		
		if(spltDB == 0 && targetNum != 1) {
			status = runKMA_spltDB(templatefilenames, targetNum, outputfilename, argc, argv, ConClave, extendedFeatures, vcf);
		} else if(mem_mode || sparseCheck(templatefilename)) {
			status = runKMA_MEM(templatefilename, outputfilename, exeBasic, ConClave, extendedFeatures, vcf);
		} else {
			status = runKMA(templatefilename, outputfilename, exeBasic, ConClave, extendedFeatures, vcf);
		}
		fprintf(stderr, "# Closing files\n");
		fflush(stdout);
	}
	
	return status;
}
