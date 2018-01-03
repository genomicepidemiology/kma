/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

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
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/*
 STRUCTURES
*/
union DNA_tree {
	int value;
	union DNA_tree *next;
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

/*
  GLOBAL VARIABLES
*/
//union DNA_tree templates; // hash for linking templates to k-mers
//int **templates_values;
struct hashMap templates;
struct hashMap_index *foundKmers;
union DNA_tree *templates_align; // hash-list for linking k-mers to posistions
char **template_seqs; //Whole sequences
char **template_names; //Convertion from int to string name
char **templates_descriptions; //Description of templates
int *template_lengths, *template_ulengths; // lengths of sequences
double **alignment_scores; // global scores for each tmeplate
int *Scores, *Scores_tot, *bestTemplates;
unsigned kmersize, DB_size, DB_mem;
int deCon, contamination;
char bases[] = "ATCGN-";
int com_bases[] = {1, 0, 3, 2, 4};
long unsigned *sizeArray;
int mapped_cont = 0;
long unsigned *convertNum;
unsigned INITIAL_SIZE = 1048576;
long unsigned MAX_SIZE;
char *input_convertion;
int (*homcmp)(int, int);
//unsigned INITIAL_SIZE = 83752;
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
	if(len1 == 0 || len2 == 0 || len1 - len2 <= 0) {
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
	if(len1 == 0 || len2 == 0 || len1 - len2 <= 0) {
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

int chomp(char *string) {
	/* remove trailing spaces and newlines */
	int k = strlen(string) - 1;
	/* isspace = ((string[k] >= 9  && string[k] <= 13) || string[k] == 32), in ASCII */
	while (isspace(string[k]))
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

int hasN(int *s1, int len) {
	int i;
	for(i = 0; i < len; i++) {
		if(s1[i] == 4) {
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

/* DNA SPECIFIC FUNCTIONS */
void decode_DNA(char *qseq, int seqlen) {
	int i;
	for(i = 0; i < seqlen; i++) {
		qseq[i] = input_convertion[qseq[i]];
	}
}

void convert_DNA(char *qseq, int *qseq_int, int seqlen) {
	int i;
	for(i = 0; i < seqlen; i++) {
		qseq_int[i] = chrpos(bases, qseq[i]);
	}
}

void rc_int(int *qseq_int, int seqlen) {
	int i, carry;
	for(i = 0; i < seqlen / 2; i++) {
		carry = com_bases[qseq_int[i]];
		qseq_int[i] = com_bases[qseq_int[seqlen - i - 1]];
		qseq_int[seqlen - i - 1] = carry;
	}
	if(seqlen % 2) {
		i = seqlen / 2 + 1;
		qseq_int[i] = com_bases[qseq_int[i]];
	}
	
}

int nextNoN(int *qseq_int, int len) {
	
	int i;
	for(i = 0; i < len; i++) {
		if(qseq_int[i] != 4) {
			return i;
		}
	}
	return len;
}

int nextN(int *qseq_int, int len) {
	
	int i;
	for(i = 0; i < len; i++) {
		if(qseq_int[i] == 4) {
			return i;
		}
	}
	return -1;
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
			if(node->next == NULL) {
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
			if(node->next == NULL) {
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
			if(node->next == NULL) {
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
			return -1;
		}
		node = (union DNA_tree*)(&(node->next[key[i]]));
	}
	
	if(node->next == NULL) {
		return -1;
	} else {
		return node->value;
	}
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
		if(node->next == NULL) {
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
			exit(1);
		}
	}
}

void load_values(int **values, FILE *file) {
	int i, size;
	fread(&size, sizeof(int), 1, file);
	values = malloc((size + 1) * sizeof(int*));
	if(values == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	values[0] = malloc(sizeof(int));
	if(values[0] == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	values[0][0] = size;
	for(i = 1; i < values[0][0]; i++) {
		fread(&size, sizeof(int), 1, file);
		values[i] = malloc((size + 1) * sizeof(int));
		if(values[i] == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		fread(values[i], (size + 1) * sizeof(int), 1, file);
	}
}

/* HASHMAP FUNCTIONS */
void initialize_hashMap(struct hashMap *dest, unsigned newSize) {
	/* set hashMap */
	dest->size = newSize;
	dest->n = 0;
	/* set hashTable */
	dest->table = calloc(newSize, sizeof(struct hashTable*));
	if(dest->table == NULL) {
		fprintf(stderr, "Out of memory\n");
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
		if(dest->table[index] == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node = dest->table[index];
		node->value = calloc(2, sizeof(unsigned));
		if(node->value == NULL) {
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

void hashMap_addIndex(struct hashMap *dest, long unsigned key, unsigned value) {
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
		if(dest->table[index] == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node = dest->table[index];
		node->value = calloc(2, sizeof(unsigned));
		if(node->value == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node->value[0] = 1;
		node->value[1] = value;
		node->key = key;
		node->next = 0;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(key == node->key) { // Keys match change value
				if(node->value[node->value[0]] != value) {
					node->value[0]++;
					node->value = realloc(node->value, (node->value[0] + 1) * sizeof(unsigned));
					if(node->value == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					node->value[node->value[0]] = value;
				}
				return;
			} else if(node->next == 0) { // This chain is filled, create next
				node->next = malloc(sizeof(struct hashTable));
				if(node->next == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				node = node->next;
				node->next = 0;
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
	dest->table = calloc(dest->size, sizeof(struct hashTable*));
	if(dest->table == NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	
	/* load hashTable */
	for(i = 0; i < dest->n; i++) {
		fread(&key, sizeof(long unsigned), 1, file);
		index = key % dest->size;
		if(dest->table[index] == 0) {
			dest->table[index] = malloc(sizeof(struct hashTable));
			if(dest->table[index] == NULL) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			node = dest->table[index];
			node->key = key;
			fread(&tmp_size, sizeof(unsigned), 1, file);
			node->value = malloc(tmp_size * sizeof(unsigned));
			if(node->value == NULL) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			fread(node->value, tmp_size * sizeof(unsigned), 1, file);
			node->next = 0;
		} else {
			for(node = dest->table[index]; node != 0; node = node->next) {
				if(node->value == 0) {
					node->key = key;
					fread(&tmp_size, sizeof(unsigned), 1, file);
					node->value = malloc(tmp_size * sizeof(unsigned));
					if(node->value == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					fread(node->value, tmp_size * sizeof(unsigned), 1, file);
				} else if(key == node->key) {
					break;
				} else if(node->next == 0) {
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
}

void initialize_hashMap_index(struct hashMap_index *dest, unsigned newSize) {
	/* set hashMap */
	dest->size = newSize;
	dest->n = 0;
	/* set hashTable */
	dest->table = calloc(newSize, sizeof(struct hashTable_index*));
	if(dest->table == NULL) {
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

void hashMap_CountIndex(struct hashMap_index *dest, long unsigned key) {
	unsigned i, j, index;
	struct hashTable_index *node;
	
	/* get index */
	index = key % dest->size;
	
	/* find pos */
	if(dest->table[index] == 0) { // New value, no collision
		dest->n++;
		dest->table[index] = malloc(sizeof(struct hashTable_index));
		if(dest->table[index] == NULL) {
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
				if(node->next == NULL) {
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

/* METHOD SPECIFIC FUNCTIONS */
int find_contamination(int *bestTemplates) {
	int i;
	for(i = bestTemplates[0]; i > 0; i--) {
		if(bestTemplates[i] == contamination) {
			return i;
		}
	}
	return -1;
}

int homcmp_or(int t, int q) {
	return (t || q);
}

int homcmp_and(int t, int q) {
	return (t && q);
}

void update_DBs(char *header, char *qseq, int q_len) {
	
	int i, j, end, Ncheck, checkNcheck, dupCheck;
	unsigned *value;
	long unsigned key;
	int *qseq_int;
	qseq_int = malloc(q_len * sizeof(int));
	if(qseq_int == NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	decode_DNA(qseq, q_len);
	
	checkNcheck = 0;
	Ncheck = 1;
	for(i = 0; i < q_len; i++) {
		qseq_int[i] = qseq[i] - '0';
	}
	
	/* fill hashMap */
	end = q_len - kmersize + 1;
	
	i = nextNoN(qseq_int, end);
	
	while((j = nextN(qseq_int + i, kmersize)) != -1) {
		i += j + 1;
	}
	if(i < end) {
		key = quinaryToDecimal(qseq_int, i, kmersize);
		hashMap_addIndex(&templates, key, DB_size);
	}
	for(i = 1; i < end; i++) {
		if(qseq_int[i + kmersize - 1] != 4) {
			key = (key - qseq_int[i - 1] * convertNum[0]) * 5 + qseq_int[i + kmersize - 1];
			hashMap_addIndex(&templates, key, DB_size);
		} else {
			i += kmersize;
			i += nextNoN(qseq_int + i, end);
			while((j = nextN(qseq_int + i, kmersize)) != -1) {
				i += j + 1;
			}
			if(i < end) {
				key = quinaryToDecimal(qseq_int, i, kmersize);
				hashMap_addIndex(&templates, key, DB_size);
			}
		}
	}
	
	free(qseq_int);
}

int update_DBs_sparse(char *qseq, int q_len, int *prefix, int prefix_len, int MinKlen, double homQ, double homT) {
	int i, j, k, end, hits, len, template;
	int *qseq_int;
	long unsigned key;
	unsigned *value;
	double score, bestT, bestQ;
	/* convert seq to numbers */
	qseq_int = calloc(2 * q_len, sizeof(int));
	if(qseq_int == NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	decode_DNA(qseq, q_len);
	for(i = 0; i < q_len; i++) {
		qseq_int[i] = qseq[i] - '0';
	}
	/* get rc */
	for(i = q_len; i > 0; i--) {
		qseq_int[2 * q_len - i] = com_bases[qseq[i - 1] - '0'];
	}
	end = q_len - kmersize + 1 - prefix_len;
	
	/* check if seq meets minimum length */
	if(prefix_len != 0) {
		for(i = 0; i < 2 * q_len && MinKlen > 0; i += q_len) {
			for(j = i; j < i + end && MinKlen > 0; j++) {
				if(int_eq((qseq_int + j), prefix, prefix_len)) {
					MinKlen--;
				}
			}
		}
		if(MinKlen > 0) {
			free(qseq_int);
			return 0;
		}
	} else if(2 * end < MinKlen){
		free(qseq_int);
		return 0;
	}
	
	/* check for homology */
	if(homQ < 1 || homT < 1) {
		struct hashTable_index *node, *node_next;
		len = 0;
		bestTemplates[0] = 0;
		for(i = 0; i < 2 * q_len; i += q_len) {
			for(j = i; j < i + end; j++) {
				/* check if prefix matches */
				if(int_eq((qseq_int + j), prefix, prefix_len)) {
					/* check if k-mer exists */
					key = quinaryToDecimal(qseq_int, j + prefix_len, kmersize);
					if((value = hashMap_getValue(&templates, key))) {
						for(k = 1; k <= value[0]; k++) {
							Scores_tot[value[k]]++;
							if(Scores_tot[value[k]] == 1) {
								bestTemplates[0]++;
								bestTemplates[bestTemplates[0]] = value[k];
							}
						}
						hashMap_CountIndex(foundKmers, key);
					}
					len++;
				}
			}
		}
		
		/* get template cov */
		for(i = 0; i < foundKmers->size; i++) {
			for(node = foundKmers->table[i]; node != 0; node = node_next) {
				node_next = node->next;
				if((value = hashMap_getValue(&templates, node->key))) {
					for(j = 1; j <= value[0]; j++) {
						Scores[value[j]]++;
					}
				}
				free(node);
			}
			foundKmers->table[i] = 0;
		}
		foundKmers->n = 0;
		
		/* get best hit */
		bestT = 0.0;
		bestQ = 0.0;
		for(i = 1; i <= bestTemplates[0]; i++) {
			template = bestTemplates[i];
			score = 1.0 * Scores[template] / template_ulengths[template];
			if(score > bestT) {
				bestT = score;
			}
			if(Scores_tot[template] > bestQ) {
				bestQ = Scores_tot[template];
			}
			Scores[template] = 0;
			Scores_tot[template] = 0;
		}
		
		/* discard */
		if((*homcmp)(bestT > homT, (1.0 * bestQ / len) > homQ)) {
			//fprintf(stderr, "%f\t%f\n", bestT, (1.0 * bestQ / len));
			free(qseq_int);
			return 0;
		}
	}
	
	/* go through seq and rc seq */
	for(i = 0; i < 2 * q_len; i += q_len) {
		for(j = i; j < i + end; j++) {
			/* check if prefix matches */
			if(int_eq((qseq_int + j), prefix, prefix_len) && hasN((qseq_int + j + prefix_len), kmersize)) {
				/* check if k-mer exists */
				key = quinaryToDecimal(qseq_int, j + prefix_len, kmersize);
				value = hashMap_getValue(&templates, key);
				if(value) { // k-mer exists
					/* check if k-mer is a duplicate, remember that the values are ordered!! */
					if(value[value[0]] == DB_size) {
						template_lengths[DB_size]++;
					} else {
						template_lengths[DB_size]++;
						template_ulengths[DB_size]++;
						hashMap_addKey(&templates, key, DB_size);
					}
				} else { // New k-mer
					template_lengths[DB_size]++;
					template_ulengths[DB_size]++;
					hashMap_addKey(&templates, key, DB_size);
				}
			}
		}
	}
	free(qseq_int);
	return 1;
}

void update_DB(int *qseq_int, int q_len) {
	
	int i;
	unsigned *value;
	long unsigned key;
	/* check and update forward strand */
	for(i = 0; i < q_len; i++) {
		key = quinaryToDecimal(qseq_int, i, kmersize);
		value = hashMap_getValue(&templates, key);
		if(value != NULL && value[value[0]] != contamination) {
			hashMap_addKey(&templates, quinaryToDecimal(qseq_int, i, kmersize), contamination);
			mapped_cont++;
		}
	}
	/* check and update reverse strand */
	rc_int(qseq_int, q_len);
	for(i = 0; i < q_len; i++) {
		key = quinaryToDecimal(qseq_int, i, kmersize);
		value = hashMap_getValue(&templates, key);
		if(value != NULL && value[value[0]] != contamination) {
			hashMap_addKey(&templates, quinaryToDecimal(qseq_int, i, kmersize), contamination);
			mapped_cont++;
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

void alignDump(char *qseq, int q_len, union DNA_tree *align, int template, const char *align_path) {
	
	int out_len;
	char *out_name;
	FILE *align_out;
	/* open outputfile */
	out_len = strlen(align_path);
	out_name = calloc(out_len + 256, sizeof(char));
	if(out_name == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	sprintf(out_name, "%s%d.b", align_path, template);
	align_out = fopen(out_name, "wb");
	if(!align_out) {
		fprintf(stderr, "File corrution: %s\n", out_name);
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

void alignDump_fly(char *qseq, int q_len, union DNA_tree *align, int template, const char *align_path) {
	
	int out_len;
	char *out_name;
	FILE *align_out;
	/* open outputfile */
	out_len = strlen(align_path);
	out_name = calloc(out_len + 256, sizeof(char));
	if(out_name == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	sprintf(out_name, "%s%d.b", align_path, template);
	align_out = fopen(out_name, "wb");
	if(!align_out) {
		fprintf(stderr, "File corrution: %s\n", out_name);
		exit(1);
	}
	/* dump seq */
	fwrite(qseq, q_len * sizeof(char), 1, align_out);
	/* dump index */
	
	/* close up */
	fclose(align_out);
	destroy_tree(align, 0);
	free(out_name);
}

/* DB LOADING */
void load_DBs_KMA(char *filename) {
	/* load DBs needed for KMA */
	int i, file_len;
	char *templatefilename;
	FILE *DB_file;
	
	/* reallocate filename, to get all DBs */
	templatefilename = strdup(filename);
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
	if(!DB_file) {
		fprintf(stderr, "File corrution: %s\n", templatefilename);
		exit(1);
	}
	templatefilename[file_len] = '\0';
	fread(&DB_size, sizeof(int), 1, DB_file);
	template_lengths = malloc(DB_size * sizeof(int));
	template_seqs = calloc(DB_size, sizeof(char*));
	template_names = malloc(DB_size * sizeof(char*));
	templates_descriptions = malloc(DB_size * sizeof(char*));
	templates_align = malloc(DB_size * sizeof(union DNA_tree));
	alignment_scores = malloc(DB_size * sizeof(double*));
	if(!template_lengths || !template_seqs || !template_names || !templates_descriptions || !templates_align || !alignment_scores) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	for(i = 0; i < DB_size; i++) {
		alignment_scores[i] = calloc(2, sizeof(double));
		if(alignment_scores[i] == NULL) {
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
	if(!DB_file) {
		fprintf(stderr, "File corrution: %s\n", templatefilename);
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
		fprintf(stderr, "File corrution: %s\n", templatefilename);
		exit(1);
	}
	templatefilename[file_len] = '\0';
	
	/* allocate DBs */
	fread(&DB_size, sizeof(int), 1, DB_file);
	DB_mem = 2 * DB_size;
	template_lengths = calloc(DB_mem, sizeof(int));
	template_ulengths = calloc(DB_mem, sizeof(int));
	template_names = malloc(DB_mem * sizeof(char*));
	if(!template_lengths || !template_ulengths || !template_names) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	
	/* load lengths */
	fread(template_lengths, sizeof(int), DB_size, DB_file);
	fread(template_ulengths, sizeof(int), DB_size, DB_file);
	fclose(DB_file);
	
	/* load names */
	strcat(templatefilename, ".name.b");
	templatefilename[file_len + 7] = '\0';
	DB_file = fopen(templatefilename, "rb");
	if(!DB_file) {
		fprintf(stderr, "File corrution: %s\n", templatefilename);
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

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# makeKMA_DB creates the databases needed to run KMA, from a list of fasta files given.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\t\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-i\t\tInput/query file name\t\tSTDIN\n");
	fprintf(helpOut, "#\t-o\t\tOutput file\t\t\t\t\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-batch\t\tBatch input file\n");
	fprintf(helpOut, "#\t-deCon\t\tFile name of contamination\tNone/False\n");
	fprintf(helpOut, "#\t-t_db\t\tAdd to existing DB\t\tNone/False\n");
	fprintf(helpOut, "#\t-k\t\tKmersize\t\t\t16\n");
	fprintf(helpOut, "#\t-ML\t\tMinimum length of templates\t0\n");
	fprintf(helpOut, "#\t-CS\t\tChain size\t\t\t8 MB\n");
	fprintf(helpOut, "#\t-MS\t\tMax chain size\t\t\t14 GB\n");
	fprintf(helpOut, "#\t-SM\t\tStart at max chain size\t\tFalse\n");
	fprintf(helpOut, "#\t-QS\t\tStart query size\t\t1M\n");
	fprintf(helpOut, "#\t-Sparse\t\tMake Sparse DB\n#\t\t\t('-' for no prefix)\t\tNone/False\n");
	fprintf(helpOut, "#\t-ht\t\tHomology template\t\t1.0\n");
	fprintf(helpOut, "#\t-hq\t\tHomology query\t\t\t1.0\n");
	fprintf(helpOut, "#\t-and\t\tBoth homolgy thresholds\n#\t\t\thas to be reached\t\tor\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int args, stop, entryCounter, filecount, filecounter;
	int i, j, q_size, q_len, out_len, l_len, zipped_len, *qseq_int, file_len, path_len;
	double homQ, homT;
	char *file, *line, *qseq, *header, *cmd, *align_path;
	int line_size = 256 * sizeof(char);
	char *zipped = ".gz";
	line = malloc(line_size);
	if(line == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	FILE *inputfile, *outputfile;
	time_t t0, t1;
	/* SET DEFAULTS */
	char **inputfiles = 0;
	filecount = 0;
	char *outputfilename = 0;
	char *templatefilename = 0;
	kmersize = 16;
	deCon = 0;
	DB_size = 1;
	MAX_SIZE = 14 * (long unsigned)(pow(2, 30) + 0.5) / 8; //14 = #GB, 2^30 = GB, 8 = sizeof pointer
	int SM = 0;
	int sparse_run = 0;
	int MinLen = 0;
	int MinKlen = 1;
	int *prefix, prefix_len;
	char *deCon_file;
	q_size = 1000000;
	contamination = -1;
	FILE *templatefile;
	homQ = 1.0;
	homT = 1.0;
	if (argc == 1) {
		fprintf(stderr, "# Too few arguments handed.\n");
		helpMessage(1);
	}
	inputfiles = malloc(sizeof(char*));
	/* set input conversion */
	input_convertion = malloc(128 * sizeof(char));
	if(input_convertion == NULL || inputfiles == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	//memset(input_convertion, 4, 128);
	for(i = 0; i < 128; i++) {
		input_convertion[i] = '4';
	}
	input_convertion['A'] = '0';
	input_convertion['T'] = '1';
	input_convertion['C'] = '2';
	input_convertion['G'] = '3';
	input_convertion['a'] = '0';
	input_convertion['t'] = '1';
	input_convertion['c'] = '2';
	input_convertion['g'] = '3';
	
	homcmp = &homcmp_or;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-i") == 0) {
			stop = 0;
			args++;
			while(stop == 0 && args < argc) {
				if(strncmp(argv[args], "-", 1) != 0 || strncmp(argv[args], "--", 2) == 0) {
					filecount++;
					inputfiles = realloc(inputfiles, filecount * sizeof(char*));
					if(inputfiles == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					inputfiles[filecount - 1] = malloc((strlen(argv[args]) + 1) * sizeof(char));
					if(inputfiles[filecount - 1] == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					strcpy(inputfiles[filecount - 1], argv[args]);
					inputfiles[filecount - 1][strlen(argv[args])] = '\0';
					args++;
				} else {
					stop = 1;
				}
			}
			args--;
		} else if(strcmp(argv[args], "-o") == 0) {
			args++;
			if(args < argc) {
				outputfilename = malloc((strlen(argv[args]) + 1) * sizeof(char));
				if(outputfilename == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(outputfilename, argv[args]);
				outputfilename[strlen(outputfilename)] = '\0';
			}
		} else if(strcmp(argv[args], "-deCon") == 0) {
			args++;
			if(args < argc) {
				deCon_file = strdup(argv[args]);
				if(deCon_file == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				deCon = 1;
			}
		} else if(strcmp(argv[args], "-t_db") == 0) {
			args++;
			if(args < argc) {
				templatefilename = strdup(argv[args]);
				if(templatefilename == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
			}
		} else if(strcmp(argv[args], "-k") == 0) {
			args++;
			if(args < argc) {
				kmersize = atoi(argv[args]);
				if(kmersize == 0) {
					fprintf(stderr, "# Invalid kmersize parsed, using default\n");
					kmersize = 16;
				} else if(kmersize > 27) {
					kmersize = 27;
				}
			}
		} else if(strcmp(argv[args], "-QS") == 0) {
			args++;
			if(args < argc) {
				q_size = atoi(argv[args]);
				if(kmersize == 0) {
					fprintf(stderr, "# Invalid QS parsed, using default\n");
					q_size = 1000000;
				}
			}
		} else if(strcmp(argv[args], "-CS") == 0) {
			args++;
			if(args < argc) {
				INITIAL_SIZE = atoi(argv[args]) * 1048576 / 8;
				if(INITIAL_SIZE == 0) {
					fprintf(stderr, "# Invalid Chain Size parsed, using default\n");
					INITIAL_SIZE = 1048576;
				}
				//INITIAL_SIZE = strtoul(argv[args], NULL, 10);
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
		} else if(strcmp(argv[args], "-MS") == 0) {
			args++;
			if(args < argc) {
				MAX_SIZE = atoi(argv[args]) * (long unsigned)(pow(2, 30) + 0.5) / 8;
				if(MAX_SIZE == 0) {
					fprintf(stderr, "# Invalid Max Size parsed, using default\n");
					MAX_SIZE = 14 * (long unsigned)(pow(2, 30) + 0.5) / 8;
				}
				//MAX_SIZE = strtoul(argv[args], NULL, 10);	
			}
		} else if(strcmp(argv[args], "-SM") == 0) {
			SM = 1;
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
					exit(1);
				}
				while(!feof(inputfile)) {
					line = fget_line(line, &line_size, &l_len, inputfile);
					if(l_len != 0) {
						filecount++;
						inputfiles = realloc(inputfiles, filecount * sizeof(char*));
						if(inputfiles == NULL) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
						inputfiles[filecount - 1] = strdup(line);
						if(inputfiles[filecount - 1] == NULL) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					}
				}
				fclose(inputfile);
			}
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			args++;
			if(args < argc) {
				sparse_run = 1;
				if(strcmp(argv[args], "-") == 0) {
					prefix_len = 0;
					prefix = malloc(sizeof(int));
					if(prefix == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					prefix[0] = 4;
				} else {
					prefix_len = strlen(argv[args]);
					prefix = malloc(prefix_len * sizeof(int));
					if(prefix == NULL) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					for(i = 0; i < prefix_len; i++) {
						prefix[i] = input_convertion[argv[args][i]] - '0';
						if(prefix[i] == 4) {
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
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(1);
		}
		args++;
	}
	if(filecount == 0) {
		inputfiles[0] = malloc(3 * sizeof(char));
		if(inputfiles[0] == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		strncpy(inputfiles[0], "--", 2);
		inputfiles[0][2] = '\0';
		filecount = 1;
	}
	
	/* validate that the needed arguments was passed */
	if(outputfilename == 0) {
		fprintf(stderr, "# Too few arguments handed.\n");
		helpMessage(1);
	}
	
	/* Set variables needed */
	zipped_len = strlen(zipped);
	cmd = malloc(sizeof(char));
	q_len = 0;
	qseq = malloc(q_size * sizeof(char));
	header = malloc(256 * sizeof(char));
	if(!cmd || !qseq || !header) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	DB_size = 1;
	if(templatefilename != 0) {
		file_len = strlen(templatefilename);
		templatefilename = realloc(templatefilename, (file_len + 10) * sizeof(char));
		if(templatefilename == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		if(sparse_run) {
			strcat(templatefilename, ".sparse.b");
			templatefilename[file_len + 9] = '\0';
			templatefile = fopen(templatefilename, "rb");
			if(!templatefile) {
				fprintf(stderr, "File corrution: %s\n", templatefilename);
				exit(1);
			}
			fread(&DB_size, sizeof(int), 1, templatefile);
			fread(&prefix_len, sizeof(int), 1, templatefile);
			prefix = malloc(prefix_len * sizeof(int));
			if(prefix == NULL) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			fread(prefix, prefix_len * sizeof(int), 1, templatefile);
		} else {
			strcat(templatefilename, ".b");
			templatefilename[file_len + 2] = '\0';
			templatefile = fopen(templatefilename, "rb");
			if(!templatefile) {
				fprintf(stderr, "File corrution: %s\n", templatefilename);
				exit(1);
			}
			fread(&DB_size, sizeof(int), 1, templatefile);
		}
		hashMap_load(&templates, templatefile);
		fclose(templatefile);
		templatefilename[file_len] = '\0';
	} else if(SM) {
		initialize_hashMap(&templates, MAX_SIZE);
	} else {
		initialize_hashMap(&templates, INITIAL_SIZE);
	}
	/* initialize convertNum */
	convertNum = malloc(kmersize * sizeof(long unsigned));
	if(convertNum == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	convertNum[kmersize - 1] = 1;
	for(i = kmersize - 2; i >= 0; i--) {
		convertNum[i] = convertNum[i+1] * 5;
	}
	
	DB_mem = 1024;
	if(sparse_run) {
		/* calculate minimum kmer length */
		if(MinLen > 0) {
			MinKlen = 2 * (MinLen - kmersize - prefix_len + 1);
			for(i = 0; i < prefix_len; i++) {
				MinKlen /= 4;
			}
		}
		/* set outputfiles */
		FILE *length_out, *name_out;
		out_len = strlen(outputfilename);
		outputfilename = realloc(outputfilename, (out_len + 10) * sizeof(char));
		if(outputfilename == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		strcat(outputfilename, ".spaLen.b");
		outputfilename[out_len + 9] = '\0';
		length_out = fopen(outputfilename, "wb");
		if(!length_out) {
			fprintf(stderr, "File corrution: %s\n", outputfilename);
			exit(1);
		}
		outputfilename[out_len] = '\0';
		
		strcat(outputfilename, ".name.b");
		outputfilename[out_len + 7] = '\0';
		name_out = fopen(outputfilename, "wb");
		if(!name_out) {
			fprintf(stderr, "File corrution: %s\n", outputfilename);
			exit(1);
		}
		outputfilename[out_len] = '\0';
		
		if(templatefilename != 0) {
			load_DBs_Sparse(templatefilename);
			DB_mem = 2 * DB_size;
			for(i = 0; i < DB_size; i++) {
				fwrite(template_names[i], 256 * sizeof(char), 1, name_out);
				free(template_names[i]);
				template_names[i] = NULL;
			}
			free(template_names);
			if(homQ < 1 || homT < 1) {
				Scores = calloc(DB_mem, sizeof(int));
				Scores_tot = calloc(DB_mem, sizeof(int));
				bestTemplates = malloc((DB_mem + 1) * sizeof(int));
				foundKmers = malloc(sizeof(struct hashMap_index));
				if(!Scores || !Scores_tot || !bestTemplates || !foundKmers) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				initialize_hashMap_index(foundKmers, 1024 * 1024);
			}
		} else {
			DB_size = 0;
			template_lengths = calloc(DB_mem, sizeof(int));
			template_ulengths = calloc(DB_mem, sizeof(int));
			if(!template_lengths || !template_ulengths) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			if(homQ < 1 || homT < 1) {
				Scores = calloc(DB_mem, sizeof(int));
				Scores_tot = calloc(DB_mem, sizeof(int));
				bestTemplates = malloc((DB_mem + 1) * sizeof(int));
				foundKmers = malloc(sizeof(struct hashMap_index));
				if(!Scores || !Scores_tot || !bestTemplates || !foundKmers) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				initialize_hashMap_index(foundKmers, 1024 * 1024);
			}
			template_lengths[0] = 0;
			template_ulengths[0] = 0;
		}
		/* Parse inputfiles */
		fprintf(stderr, "# Updating DBs\n");
		for(filecounter = 0; filecounter < filecount; filecounter++) {
			file = inputfiles[filecounter];
			if(strcmp(file, "--") == 0) {
				inputfile = stdin;
			} else if(strncmp(file + (strlen(file) - zipped_len), zipped, zipped_len) == 0) {
				cmd = realloc(cmd, (strlen(file) + strlen("gunzip -c ") + 1) * sizeof(char));
				if(!cmd) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				sprintf(cmd, "gunzip -c %s", file);
				cmd[strlen(file) + strlen("gunzip -c ") + 1] = '\0';
				inputfile = popen(cmd, "r");
				if(!inputfile) {
					fprintf(stderr, "File corrution: %s\n", file);
					exit(1);
				}
			} else {
				inputfile = fopen(file, "r");
				if(!inputfile) {
					fprintf(stderr, "File corrution: %s\n", file);
					exit(1);
				}
			}
			if(inputfile == NULL) {
				fprintf(stderr, "Cannot acces file:\t%s\n", file);
				exit(1);
			}
			/* inputfile */
			q_len = 0;
			while(!feof(inputfile)) {
				if(inputfile && fgets(line, line_size, inputfile) != NULL) {
					if(strncmp(line, ">", 1) == 0) { //Header
						if(q_len > MinLen && update_DBs_sparse(qseq, q_len, prefix, prefix_len, MinKlen, homQ, homT)) {
							fprintf(stderr, "# Adding:\t%s\n", header);
							fwrite(header, 256 * sizeof(char), 1, name_out);
							DB_size++;
							if(DB_size >= DB_mem) {
								template_lengths = realloc(template_lengths, 2 * DB_mem * sizeof(int));
								template_ulengths = realloc(template_ulengths, 2 * DB_mem * sizeof(int));
								if(!template_lengths || !template_ulengths) {
									fprintf(stderr, "OOM\n");
									exit(1);
								}
								for(i = DB_mem; i < DB_mem * 2; i++) {
									template_lengths[i] = 0;
									template_ulengths[i] = 0;
								}
								if(homQ < 1 || homT < 1) {
									free(Scores);
									Scores = calloc(2 * DB_mem, sizeof(int));
									free(Scores_tot);
									Scores_tot = calloc(2 * DB_mem, sizeof(int));
									free(bestTemplates);
									bestTemplates = malloc((2 * DB_mem + 1) * sizeof(int));
									if(!Scores || !Scores_tot || !bestTemplates) {
										fprintf(stderr, "OOM\n");
										exit(1);
									}
								}
								DB_mem *= 2;
							}
							if((homQ < 1 || homT < 1) && q_len > foundKmers->size) {
								free(foundKmers->table);
								initialize_hashMap_index(foundKmers, 2 * q_len);
							}
						} else if(q_len != 0) {
							fprintf(stderr, "# Not added:\t%s\n", header);
						}
						/* start over */
						q_len = 0;
						qseq[0] = '\0';
						l_len = chomp(line);
						strcpy(header, (line + 1));
						header[l_len - 1] = '\0';
					} else { //Seq
						l_len = chomp(line);
						if((q_len + l_len) >= q_size) { //realloc qseq if needed
							q_size *= 2;
							qseq = realloc(qseq, q_size * sizeof(char));
							if(qseq == NULL) {
								fprintf(stderr, "OOM\n");
								exit(1);
							}
						}
						strncpy(qseq + q_len, line, l_len);
						q_len += l_len;
					}
				}
			}
			if(q_len > MinLen && update_DBs_sparse(qseq, q_len, prefix, prefix_len, MinKlen, homQ, homT)) {
				fprintf(stderr, "# Adding:\t%s\n", header);
				fwrite(header, 256 * sizeof(char), 1, name_out);
				DB_size++;
				if(DB_size >= DB_mem) {
					template_lengths = realloc(template_lengths, 2 * DB_mem * sizeof(int));
					template_ulengths = realloc(template_ulengths, 2 * DB_mem * sizeof(int));
					if(!template_lengths || !template_ulengths) {
						fprintf(stderr, "OOM\n");
						exit(1);
					}
					for(i = DB_mem; i < DB_mem * 2; i++) {
						template_lengths[i] = 0;
						template_ulengths[i] = 0;
					}
					if(homQ < 1 || homT < 1) {
						free(Scores);
						Scores = calloc(2 * DB_mem, sizeof(int));
						free(Scores_tot);
						Scores_tot = calloc(2 * DB_mem, sizeof(int));
						free(bestTemplates);
						bestTemplates = malloc((2 * DB_mem + 1) * sizeof(int));
						if(!Scores || !Scores_tot || !bestTemplates) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					}
					DB_mem *= 2;
				}
				if((homQ < 1 || homT < 1) && q_len > foundKmers->size) {
					free(foundKmers->table);
					initialize_hashMap_index(foundKmers, 2 * q_len);
				}
				q_len = 0;
				qseq[0] = '\0';
			} else if(q_len != 0) {
				fprintf(stderr, "# Not added:\t%s\n", header);
			}
			q_len = 0;
			qseq[0] = '\0';
			if(strncmp(file + (strlen(file) - zipped_len), zipped, zipped_len) == 0) {
				pclose(inputfile);
			} else {
				fclose(inputfile);
			}
		}
		
		/* Dump remaining DBs */
		fprintf(stderr, "# Dumping DB\n");
		fwrite(&(int){DB_size}, sizeof(int), 1, length_out);
		fwrite(template_lengths, DB_size * sizeof(int), 1, length_out);
		fwrite(template_ulengths, DB_size * sizeof(int), 1, length_out);
		
		outputfilename[out_len] = '\0';
		strcat(outputfilename, ".sparse.b");
		outputfilename[out_len + 9] = '\0';
		outputfile = fopen(outputfilename, "wb");
		if(!outputfile) {
			fprintf(stderr, "File corrution: %s\n", outputfilename);
			exit(1);
		}
		fwrite(&(int){DB_size}, sizeof(int), 1, outputfile);
		fwrite(&(int){prefix_len}, sizeof(int), 1, outputfile);
		fwrite(prefix, prefix_len * sizeof(int), 1, outputfile);
		//dump_tree(&templates, outputfile, 0);
		hashMap_dump(&templates, outputfile);
		fclose(outputfile);
		outputfilename[out_len] = '\0';
		fclose(length_out);
		fclose(name_out);
		
		fprintf(stderr, "# templates key-value pairs:\t%lu.\n", templates.n);// / 1048576);
		if(deCon == 0) {
			fprintf(stderr, "Done\n");
			exit(0);
		}
		fprintf(stderr, "Done\n");
		exit(0);
	}
	
	/* set outputfiles */
	FILE *length_out, *name_out, *desc_out, *align_in, *align_out;
	char *in_name, *align_path_in;
	int in_path_len;
	out_len = strlen(outputfilename);
	outputfilename = realloc(outputfilename, (out_len + 17) * sizeof(char));
	if(outputfilename == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	strcat(outputfilename, ".align.b");
	outputfilename[out_len + 8] = '\0';
	align_out = fopen(outputfilename, "wb");
	if(!align_out) {
		fprintf(stderr, "File corrution: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[out_len] = '\0';
	
	strcat(outputfilename, ".length.b");
	outputfilename[out_len + 9] = '\0';
	length_out = fopen(outputfilename, "wb");
	if(!length_out) {
		fprintf(stderr, "File corrution: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[out_len] = '\0';
	
	strcat(outputfilename, ".name.b");
	outputfilename[out_len + 7] = '\0';
	name_out = fopen(outputfilename, "wb");
	if(!name_out) {
		fprintf(stderr, "File corrution: %s\n", outputfilename);
		exit(1);
	}
	outputfilename[out_len] = '\0';
	
	if(templatefilename != 0) {
		load_DBs_KMA(templatefilename);
		DB_mem = DB_size;
		/* cp existing db to new DB */
		for(i = 0; i < DB_size; i++) {
			fwrite(template_names[i], 256 * sizeof(char), 1, name_out);
			free(template_names[i]);
			template_names[i] = NULL;
		}
		free(template_names);
		
		strcat(templatefilename, ".align.b");
		templatefilename[file_len + 8] = '\0';
		align_in = fopen(templatefilename, "rb");
		if(!align_in) {
			fprintf(stderr, "File corrution: %s\n", templatefilename);
			exit(1);
		}
		for(i = 0; i < DB_size; i++) {
			/* load seq */
			template_seqs[i] = malloc((template_lengths[i] + 1) * sizeof(char));
			if(template_seqs[i] == NULL) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			fread(template_seqs[i], template_lengths[i] * sizeof(char), 1, align_in);
			template_seqs[i][template_lengths[i]] = '\0';
			
			/* dump seq */
			fwrite(template_seqs[i], template_lengths[i] * sizeof(char), 1, align_out);
			free(template_seqs[i]);
			template_seqs[i] = NULL;
		}
		free(template_seqs);
		/* alloc for future DB */
		template_lengths = realloc(template_lengths, 2 * DB_mem * sizeof(int));
		if(template_lengths == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		DB_mem *= 2;
	} else {
		DB_size = 1;
		template_lengths = malloc(DB_mem * sizeof(int));
		if(template_lengths == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		template_lengths[0] = 0;
		fwrite(header, 256 * sizeof(char), 1, name_out);
	}
	
	/* Parse inputfiles */
	fprintf(stderr, "# Updating DBs\n");
	for(filecounter = 0; filecounter < filecount; filecounter++) {
		file = inputfiles[filecounter];
		if(strcmp(file, "--") == 0) {
			inputfile = stdin;
		} else if(strncmp(file + (strlen(file) - zipped_len), zipped, zipped_len) == 0) {
			cmd = realloc(cmd, (strlen(file) + strlen("gunzip -c ") + 1) * sizeof(char));
			if(cmd == NULL) {
				fprintf(stderr, "OOM\n");
				exit(1);
			}
			sprintf(cmd, "gunzip -c %s", file);
			cmd[strlen(file) + strlen("gunzip -c ") + 1] = '\0';
			inputfile = popen(cmd, "r");
			if(!inputfile) {
				fprintf(stderr, "File corrution: %s\n", file);
				exit(1);
			}
		} else {
			inputfile = fopen(file, "r");
			if(!inputfile) {
				fprintf(stderr, "File corrution: %s\n", file);
				exit(1);
			}
		}
		if(inputfile == NULL) {
			fprintf(stderr, "Cannot acces file:\t%s\n", file);
			exit(1);
		}
		/* inputfile */
		while(!feof(inputfile)) {
			if(inputfile && fgets(line, line_size, inputfile) != NULL) {
				if(strncmp(line, ">", 1) == 0) { //Header
					if(q_len > MinLen) {
						fprintf(stderr, "# Adding:\t%s\n", header);
						update_DBs(header, qseq, q_len); //Add seq to DB hashes
						/* Dump annots & tree */
						fwrite(qseq, q_len * sizeof(char), 1, align_out);
						fwrite(header, 256 * sizeof(char), 1, name_out);
						template_lengths[DB_size] = q_len;
						DB_size++;
						if(DB_size >= DB_mem) {
							template_lengths = realloc(template_lengths, 2 * DB_mem * sizeof(int));
							if(template_lengths == NULL) {
								fprintf(stderr, "OOM\n");
								exit(1);
							}
							DB_mem *= 2;
						}
						
					}
					
					q_len = 0;
					qseq[0] = '\0';
					l_len = chomp(line);
					strcpy(header, (line + 1));
					header[l_len] = '\0';
				} else { //Seq
					l_len = chomp(line);
					if((q_len + l_len) >= q_size) { //realloc qseq if needed
						q_size *= 2;
						qseq = realloc(qseq, q_size * sizeof(char));
						if(qseq == NULL) {
							fprintf(stderr, "OOM\n");
							exit(1);
						}
					}
					strncpy(qseq + q_len, line, l_len);
					q_len += l_len;
					qseq[q_len] = '\0';
				}
			}
		}
		if(q_len > MinLen) {
			fprintf(stderr, "# Adding:\t%s\n", header);
			update_DBs(header, qseq, q_len); //Add seq to DB hashes
			/* Dump annots & tree */
			fwrite(qseq, q_len * sizeof(char), 1, align_out);
			fwrite(header, 256 * sizeof(char), 1, name_out);
			template_lengths[DB_size] = q_len;
			DB_size++;
			if(DB_size >= DB_mem) {
				template_lengths = realloc(template_lengths, 2 * DB_mem * sizeof(int));
				if(template_lengths == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				DB_mem *= 2;
			}
		}
		q_len = 0;
		qseq[0] = '\0';
		fclose(inputfile);
	}
	
	/* Dump remaining DBs */
	fprintf(stderr, "# Dumping DB\n");
	fwrite(&(int){DB_size}, sizeof(int), 1, length_out);
	fwrite(template_lengths, DB_size * sizeof(int), 1, length_out);
	
	outputfilename[out_len] = '\0';
	strcat(outputfilename, ".b");
	outputfilename[out_len + 2] = '\0';
	outputfile = fopen(outputfilename, "wb");
	if(!outputfile) {
		fprintf(stderr, "File corrution: %s\n", outputfilename);
		exit(1);
	}
	fwrite(&(int){DB_size}, sizeof(int), 1, outputfile);
	
	hashMap_dump(&templates, outputfile);
	fclose(outputfile);
	outputfilename[out_len] = '\0';
	fclose(length_out);
	fclose(name_out);
	
	fprintf(stderr, "# templates key-value pairs:\t%lu.\n", templates.n);// / 1048576);
	if(deCon == 0) {
		fprintf(stderr, "Done\n");
		return 0;
	}
	contamination = DB_size;
	q_len = 0;
	entryCounter = 0;
	//q_size = 250000000;
	qseq_int = malloc(q_size * sizeof(int));
	if(qseq_int == NULL) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	/* open file */
	file = deCon_file;
	if(strncmp(file + (strlen(file) - zipped_len), zipped, zipped_len) == 0) {
		cmd = realloc(cmd, (strlen(file) + strlen("gunzip -c ") + 1) * sizeof(char));
		if(cmd == NULL) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		sprintf(cmd, "gunzip -c %s", file);
		cmd[strlen(file) + strlen("gunzip -c ") + 1] = '\0';
		inputfile = popen(cmd, "r");
		if(!inputfile) {
			fprintf(stderr, "File corrution: %s\n", file);
			exit(1);
		}
	} else if(strcmp(file, "--") == 0) {
		inputfile = stdin;
	} else {
		inputfile = fopen(file, "r");
		if(!inputfile) {
			fprintf(stderr, "File corrution: %s\n", file);
			exit(1);
		}
	}
	if(inputfile == NULL) {
		fprintf(stderr, "Cannot read file:\t#%s#\n", file);
		exit( 1 );
	} else {
		fprintf(stderr, "%s\t%s\n", "# Reading inputfile: ", file);
	}
	
	/* parse the file */
	while(!feof(inputfile) && fgets(line, line_size, inputfile) != NULL) {
		if(line[0] == '>') { // new seq
			fprintf(stderr, "# Adding entry number: %d\n", entryCounter);
			entryCounter++;
			update_DB(qseq_int, q_len);
			q_len = 0;
		} else { // load seq
			l_len = chomp(line);
			/* realloc qseq if needed */
			if(l_len + q_len > q_size) {
				q_size *= 2;
				qseq_int = realloc(qseq_int, q_size * sizeof(int));
				if(qseq_int == NULL) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
			}
			/* Update qseq */
			for(j = 0; j < l_len; j++) {
				qseq_int[q_len + j] = input_convertion[line[j]] - '0';
				if(qseq_int[q_len + j] == -1) {
					qseq_int[q_len + j] = 4;
				}
			}
			q_len += l_len;
		}
	}
	if(q_len > 0) {
		fprintf(stderr, "# Adding entry number: %d\n", entryCounter);
		entryCounter++;
		update_DB(qseq_int, q_len);
		q_len = 0;
	}
	/* close file */
	if(strncmp(file + (strlen(file) - zipped_len), zipped, zipped_len) == 0) {
		pclose(inputfile);
	} else {
		fclose(inputfile);
	}
	/* dump DB */
	fprintf(stderr, "# Contamination information added.\n");
	fprintf(stderr, "# %d kmers mapped to the DB.\n", mapped_cont);
	fprintf(stderr, "# Contamination mapped to %f %% of the DB.\n", 100.0 * mapped_cont / templates.n);
	fprintf(stderr, "# Dumping DB.\n");
	outputfilename[out_len] = '\0';
	strcat(outputfilename, ".decontaminated.b");
	outputfilename[out_len + strlen(".decontaminated.b")] = '\0';
	outputfile = fopen(outputfilename, "wb");
	if(!outputfile) {
		fprintf(stderr, "File corrution: %s\n", outputfilename);
		exit(1);
	}
	fwrite(&(int){DB_size}, sizeof(int), 1, outputfile);
	hashMap_dump(&templates, outputfile);
	fclose(outputfile);
	
	fprintf(stderr, "Done\n");
	fflush(0);
	return 0;
}
