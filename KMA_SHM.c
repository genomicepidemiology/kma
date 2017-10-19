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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

/*
 STRUCTURES
*/
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

struct dump {
	long unsigned size;
	long unsigned n;
	unsigned valueSize;
	key_t e_id; //size long
	key_t t_id; //n long
	key_t v_id; // valueSize long
};

struct valuesTable {
	unsigned *values;
	unsigned value_index;
	struct valuesTable *next;
};

struct valuesHash {
	unsigned n;
	unsigned size;
	struct valuesTable **table;
};

unsigned kmersize, contamination, sparse_run, magic_number;

/*
 FUNCTIONS
*/

/* BASIC FUNCTIONS */
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
	/* remove trailing spaces and newlines */
	int k = strlen(string) - 1;
	/* isspace = ((string[k] >= 9  && string[k] <= 13) || string[k] == 32), in ASCII */
	//while ((string[k] >= 9  && string[k] <= 13) || string[k] == 32)
	while(string[k] == '\n' || isspace(string[k]))
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

int intArray_pos(const int *s1, const int *s2, int len1, int len2) {
	
	int *ptr, i, end;
	
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	ptr = (int*)(s1);
	end = len1 - len2;
	for(i = 0; i <= end; i++) {
		if(*ptr == *s2 && int_eq(ptr + 1, s2 + 1, s2[0])) {
			return i;
		}
		ptr++;
	}
	return -1;
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

int uintArray_pos(const unsigned *s1, const unsigned *s2, int len1, int len2) {
	
	int i, end;
	unsigned *ptr;
	
	if(len1 == 0 || len2 == 0 || len1 - len2 < 0) {
		return -1;
	}
	
	ptr = (unsigned*)(s1);
	end = len1 - len2;
	for(i = 0; i <= end; i++) {
		if(*ptr == *s2 && uint_eq(ptr + 1, s2 + 1, s2[0])) {
			return i;
		}
		ptr++;
	}
	return -1;
}

long unsigned fget_seq(char *line, long unsigned line_size, FILE *file) {
	
	int i = 0;
	long unsigned grow = line_size;
	
	while(!feof(file) && fgets((line + i), grow, file) != NULL) {
		if(line[strlen(line) - 1] != '\n') {
			i = strlen(line) - 1;
			line_size += grow;
			line = realloc(line, line_size);
			if(line == NULL) {
				fprintf(stderr, "# Memory error.\n");
				exit(-1);
			}
		} else {
			return line_size;
		}
	}
	return line_size;
}

long unsigned fget_to(char *line, long unsigned line_size, FILE *file, char stopChar) {
	int c, i;
	
	i = 0;
	while((c = fgetc(file)) != EOF && (line[i] = c) != stopChar) {
		if((i++) == line_size) {
			line_size *= 2;
			line = realloc(line, line_size);
			if(line == NULL) {
				fprintf(stderr, "Memory error\n");
				exit(-1);
			}
		}
	}
	line[i] = '\0';
	
	return line_size;
}

/* VALUES HASH */

struct valuesHash * initialize_hashValues(unsigned size) {
	struct valuesHash *dest;
	
	dest = malloc(sizeof(struct valuesHash));
	if(!dest) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	dest->n = 0;
	dest->size = size;
	
	dest->table = calloc(size, sizeof(struct valuesTable *));
	if(dest->table == NULL) {
		fprintf(stderr, "calloc ERROR\n");
		exit(-1);
	}
	
	return dest;
}

unsigned valuesHash_add(struct valuesHash *dest, unsigned *newValues, unsigned v_index) {
	
	long unsigned key;
	unsigned i, index;
	struct valuesTable *node;
	
	/* construct key */
	key = 0;
	for(i = 0; i <= newValues[0]; i++) {
		key = key * magic_number + newValues[i];
	}
	
	/* get index */
	index = key % dest->size;
	
	if(dest->table[index] == 0) { //New value
		dest->n++;
		dest->table[index] = malloc(sizeof(struct valuesTable));
		if(!dest->table[index]) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		node = dest->table[index];
		node->values = malloc((newValues[0] + 1) * sizeof(unsigned));
		if(!node->values) {
			fprintf(stderr, "OOM\n");
			exit(1);
		}
		for(i = 0; i <= newValues[0]; i++) {
			node->values[i] = newValues[i];
		}
		node->value_index = v_index;
		node->next = 0;
		return v_index;
	} else {
		for(node = dest->table[index]; node != 0; node = node->next) {
			if(node->values[0] == newValues[0] && uint_eq(node->values + 1, newValues + 1, newValues[0])) { // Value exists
				return node->value_index;
			} else if(node->next == 0) { // Make room for new value
				dest->n++;
				node->next = malloc(sizeof(struct valuesTable));
				if(!node->next) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				node = node->next;
				node->values = malloc((newValues[0] + 1) * sizeof(unsigned));
				if(!node->values) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				for(i = 0; i <= newValues[0]; i++) {
					node->values[i] = newValues[i];
				}
				node->next = 0;
				node->value_index = v_index;
				return v_index;
			}
		}
	}
	fprintf(stderr, "Something is wrong...\nCall Philip.\n");
	exit(-1);
}


/* HASHMAP FUNCTIONS */

void hashMap_shm_detach(struct hashMap_shm *dest) {
	shmdt(dest->exist);
	shmdt(dest->table);
	shmdt(dest->values);
}

void hashMap_shm_detach_some(struct hashMap_shm *dest) {
	shmdt(dest->exist);
	shmdt(dest->table);
}

struct dump initialize_hashMap_shm(struct hashMap_shm *dest, const char *filename) {
	
	key_t key;
	int shmid;
	struct dump info;
	
	/* Attach sharings */
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap e\n");
		exit(-1);
	} else {
		info.e_id = key;
	}
	dest->exist = (unsigned *) shmat(shmid, NULL, 0);
	
	key = ftok(filename, 't');
	shmid = shmget(key, dest->n * sizeof(struct hashTable_shm), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap t\n");
		shmid = shmget(info.e_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		exit(-1);
	} else {
		info.t_id = key;
	}
	dest->table = (struct hashTable_shm *) shmat(shmid, NULL, 0);
	
	/*key = ftok(filename, 'v');
	shmid = shmget(key, valueSize * sizeof(unsigned), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
		shmid = shmget(info.e_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		key = ftok(filename, 't');
		shmid = shmget(info.t_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		exit(-1);
	} else {
		info.v_id = key;
	}
	dest->values = (unsigned *) shmat(shmid, NULL, 0);
	*/
	info.size = dest->size;
	info.n = dest->n;
	//info.valueSize = valueSize;
	
	return info;
}

struct dump initialize_values_shm(struct hashMap_shm *dest, const char *filename, struct dump info) {
	
	key_t key;
	int shmid;
	
	key = ftok(filename, 'v');
	shmid = shmget(key, info.valueSize * sizeof(unsigned), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
		shmid = shmget(info.e_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		key = ftok(filename, 't');
		shmid = shmget(info.t_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		exit(-1);
	} else {
		info.v_id = key;
	}
	dest->values = (unsigned *) shmat(shmid, NULL, 0);
	
	return info;
}

void hashMap_shm_dump(struct dump info, FILE* file) {
	
	fwrite(&contamination, sizeof(unsigned), 1, file);
	fwrite(&kmersize, sizeof(unsigned), 1, file);
	fwrite(&info.size, sizeof(long unsigned), 1, file);
	fwrite(&info.n, sizeof(long unsigned), 1, file);
	fwrite(&info.valueSize, sizeof(unsigned), 1, file);
	fwrite(&info.e_id, sizeof(key_t), 1, file);
	fwrite(&info.t_id, sizeof(key_t), 1, file);
	fwrite(&info.v_id, sizeof(key_t), 1, file);
	
}

void hashMap_shm_load(struct hashMap_shm *dest, FILE* file) {
	
	key_t key;
	int id, i;
	unsigned valueSize;
	
	fread(&contamination, sizeof(unsigned), 1, file);
	fread(&kmersize, sizeof(unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	fread(&valueSize, sizeof(unsigned), 1, file);
	
	/* get tables */
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, dest->size * sizeof(unsigned), 0777);
	dest->exist = (unsigned *) shmat(id, NULL, 0);
	
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, dest->n * sizeof(struct hashTable_shm), 0777);
	dest->table = (struct hashTable_shm *) shmat(id, NULL, 0);
	
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, valueSize * sizeof(unsigned), 0777);
	dest->values = (unsigned *) shmat(id, NULL, 0);
	
}

unsigned getValueSize(FILE* file) {
	
	unsigned v_update, valueSize, n, i;
	
	/* jump basic info */
	if(sparse_run) {
		fseek(file, sizeof(unsigned), SEEK_SET);
		fread(&i, sizeof(unsigned), 1, file);
		if(i != 0) {
			fseek(file, i * sizeof(int), SEEK_CUR);
		} else {
			fseek(file, sizeof(int), SEEK_CUR);
		}
	} else {
		fseek(file, 2 * sizeof(unsigned) + sizeof(long unsigned), SEEK_SET);
	}
	
	/* get number of key-value pairs */
	fread(&n, sizeof(long unsigned), 1, file);
	
	/* get valueSize */
	valueSize = 0;
	for(i = 0; i < n; i++) {
		fseek(file, sizeof(long unsigned), SEEK_CUR);
		fread(&v_update, sizeof(unsigned), 1, file);
		fseek(file, v_update * sizeof(unsigned), SEEK_CUR);
		valueSize += v_update;
	}
	
	/* rewind */
	rewind(file);
	
	return valueSize;
}

struct dump hashMap_shm_fill(struct hashMap_shm *dest, const char *filename) {
	
	int i, j;
	unsigned t_index, v_index, index, new_index, node_index, valueSize_in, v_update;
	long unsigned key;
	FILE* file;
	struct dump dump_info;
	unsigned *exist; //size long
	struct hashTable_shm *table; //n long
	unsigned *values; // undef long
	struct valuesHash *shmValues;
	struct valuesTable *node, *prev;
	
	
	/* get hash */
	file = fopen(filename, "rb");
	if(!file) {
		fprintf(stderr, "File corruption: %s\n", filename);
		exit(1);
	}
	
	if(sparse_run) {
		fread(&contamination, sizeof(unsigned), 1, file);
		fread(&i, sizeof(int), 1, file);
		if(i != 0) {
			fseek(file, i * sizeof(int), SEEK_CUR);
		} else {
			fseek(file, sizeof(int), SEEK_CUR);
		}
	} else {
		fread(&contamination, sizeof(unsigned), 1, file);
	}
	fread(&kmersize, sizeof(unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(long unsigned), 1, file);
	magic_number = contamination;
	shmValues = initialize_hashValues(dest->size);
	/* allocate hash */
	dump_info = initialize_hashMap_shm(dest, filename);
	/* fill arrays */
	exist = dest->exist;
	table = dest->table;
	valueSize_in = 512;
	values = malloc(valueSize_in * sizeof(unsigned));
	if(values == NULL) {
		fprintf(stderr, "Cannot allocate enough memory.\n");
		hashMap_shm_detach_some(dest);
		exit(-1);
	}
	for(i = 0; i < dest->size; i++) {
		exist[i] = dest->n + 1;
	}
	t_index = 0;
	v_index = 0;
	
	/* load key value airs */
	for(i = 0; i < dest->n; i++) {
		fread(&key, sizeof(long unsigned), 1, file);
		
		/* get index */
		index = key % dest->size;
		
		/* chenck index */
		if(exist[index] == dest->n + 1) { // index is unused
			exist[index] = t_index;
		} else { // index in use
			node_index = exist[index];
			/* link through the list */
			while(table[node_index].next != 0) {
				node_index = table[node_index].next;
			}
			/* change index */
			table[node_index].next = t_index;
		}
		
		/* update table */
		table[t_index].key = key;
		
		/* update value */
		fread(&v_update, sizeof(unsigned), 1, file);
		
		/* realloc if neccessary */
		if(v_update > valueSize_in) {
			valueSize_in = 2 * v_update;
			values = realloc(values, valueSize_in * sizeof(unsigned));
			if(values == NULL) {
				fprintf(stderr, "Cannot allocate enough memory.\n");
				hashMap_shm_detach_some(dest);
				exit(-1);
			}
		}
		
		fread(values, v_update * sizeof(unsigned), 1, file);
		new_index = valuesHash_add(shmValues, values, v_index);
		table[t_index].value = new_index;
		
		if(new_index == v_index) {
			v_index += v_update;
		}
		
		t_index++;
		table[t_index].next = 0;
		
	}
	fclose(file);
	
	/* cpy values to shared memory */
	dump_info.valueSize = v_index;
	dump_info = initialize_values_shm(dest, filename, dump_info);
	values = dest->values;
	for(i = 0; i < shmValues->size; i++) {
		node = shmValues->table[i];
		while(node != 0) {
			for(j = 0; j <= node->values[0]; j++) {
				values[node->value_index + j] = node->values[j];
			}
			prev = node;
			node = node->next;
			free(prev->values);
			free(prev);
		}
	}
	
	fprintf(stderr, "# Reduced values from: %lu to %d\n", dest->n, shmValues->n);
	free(shmValues->table);
	free(shmValues);
	
	return dump_info;
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

void hashMap_shm_destroy(FILE* file, char *filename) {
	
	key_t key;
	int id;
	long unsigned size, n;
	unsigned valueSize;
	
	fseek(file, 2 * sizeof(unsigned), SEEK_SET);
	fread(&size, sizeof(long unsigned), 1, file);
	fread(&n, sizeof(long unsigned), 1, file);
	fread(&valueSize, sizeof(unsigned), 1, file);
	
	/* get tables */
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, size * sizeof(unsigned), 0777);
	if(id >= 0) {
		shmctl(id, IPC_RMID, NULL);
	} else {
		key = ftok(filename, 'e');
		id = shmget(key, size * sizeof(unsigned), 0777);
		if(id >= 0)
			shmctl(id, IPC_RMID, NULL);
	}
	
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, n * sizeof(struct hashTable_shm), 0777);
	shmctl(id, IPC_RMID, NULL);
	if(id >= 0) {
		shmctl(id, IPC_RMID, NULL);
	} else {
		key = ftok(filename, 't');
		id = shmget(key, n * sizeof(struct hashTable_shm), 0777);
		if(id >= 0)
			shmctl(id, IPC_RMID, NULL);
	}
	
	fread(&key, sizeof(key_t), 1, file);
	id = shmget(key, valueSize * sizeof(unsigned), 0777);
	shmctl(id, IPC_RMID, NULL);
	if(id >= 0) {
		shmctl(id, IPC_RMID, NULL);
	} else {
		key = ftok(filename, 'v');
		id = shmget(key, valueSize * sizeof(unsigned), 0777);
		if(id >= 0)
			shmctl(id, IPC_RMID, NULL);
	}
	
}

void DumpConversion(struct hashMap_shm *dest, unsigned valueSize, FILE *file) {
	
	fwrite(&(int){dest->size}, sizeof(unsigned), 1, file);
	fwrite(dest->exist, dest->size * sizeof(unsigned), 1, file);
	fwrite(&(int){dest->n}, sizeof(unsigned), 1, file);
	fwrite(dest->table, dest->n * sizeof(struct hashTable_shm), 1, file);
	fwrite(&(int){valueSize}, sizeof(unsigned), 1, file);
	fwrite(dest->values, valueSize * sizeof(unsigned), 1, file);
	
}

struct dump LoadConversion(struct hashMap_shm *dest, const char *filename) {
	
	int valueSize, shmid;
	key_t key;
	FILE *file;
	struct dump dump_info;
	char *filename_mb;
	
	filename_mb = malloc(strlen(filename) + 2);
	if(!filename_mb) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	strcpy(filename_mb, filename);
	filename_mb[strlen(filename) - 1] = '\0';
	strcat(filename_mb, "mb");
	file = fopen(filename_mb, "rb");
	if(!file) {
		fprintf(stderr, "File corruption: %s\n", filename_mb);
		exit(1);
	}
	/* get exist */
	fread(&dest->size, sizeof(unsigned), 1, file);
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap e\n");
		exit(-1);
	} else {
		dump_info.e_id = key;
		dump_info.size = dest->size;
	}
	dest->exist = (unsigned *) shmat(shmid, NULL, 0);
	fread(dest->exist, dest->size * sizeof(unsigned), 1, file);
	
	/* get table */
	fread(&dest->n, sizeof(unsigned), 1, file);
	key = ftok(filename, 't');
	shmid = shmget(key, dest->n * sizeof(struct hashTable_shm), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
		shmid = shmget(dump_info.e_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		exit(-1);
	} else {
		dump_info.t_id = key;
		dump_info.n = dest->n;
	}
	dest->table = (struct hashTable_shm *) shmat(shmid, NULL, 0);
	fread(dest->table, dest->n * sizeof(struct hashTable_shm), 1, file);
	
	/* get values */
	fread(&valueSize, sizeof(unsigned), 1, file);
	key = ftok(filename, 'v');
	shmid = shmget(key, valueSize * sizeof(unsigned), IPC_CREAT | 0777);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared hashMap v\n");
		shmid = shmget(dump_info.e_id, dest->size * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		shmid = shmget(dump_info.t_id, dest->n * sizeof(unsigned), IPC_CREAT | 0777);
		shmctl(shmid, IPC_RMID, NULL);
		exit(-1);
	} else {
		dump_info.v_id = key;
		dump_info.valueSize = valueSize;
	}
	dest->values = (unsigned *) shmat(shmid, NULL, 0);
	fread(dest->values, valueSize * sizeof(unsigned), 1, file);
	
	fclose(file);
	
	return dump_info;
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma_shm sets up a shared database for mapping with KMA.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-deCon\t\tSetup contamination DB\t\tFalse\n");
	fprintf(helpOut, "#\t-destroy\tDestroy shared DB\t\tFalse\n");
	fprintf(helpOut, "#\t-Sparse\t\tSparse DB\t\t\tFalse\n");
	fprintf(helpOut, "#\t-IOD\t\tDump converted DB\t\tFalse\n");
	fprintf(helpOut, "#\t-IOL\t\tLoad previous converted DB\tFalse\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int args, file_len, deCon, destroy, IOL, IOD;
	time_t t0, t1;
	char *templatefilename;
	FILE *file;
	struct dump dump_info;
	struct hashMap_shm *templates;
	/* SET DEFAULTS */
	templatefilename = 0;
	deCon = 0;
	destroy = 0;
	sparse_run = 0;
	IOD = 0;
	IOL = 0;
	
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			args++;
			if(args < argc) {
				templatefilename = malloc((strlen(argv[args]) + 1) * sizeof(char));
				if(!templatefilename) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(templatefilename, argv[args]);
				templatefilename[strlen(argv[args])] = '\0';
			}
		} else if(strcmp(argv[args], "-deCon") == 0) {
			deCon = 1;
		} else if(strcmp(argv[args], "-destroy") == 0) {
			destroy = 1;
		} else if(strcmp(argv[args], "-Sparse") == 0) {
			sparse_run = 1;
		} else if(strcmp(argv[args], "-IOD") == 0) {
			IOD = 1;
		} else if(strcmp(argv[args], "-IOL") == 0) {
			IOL = 1;
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else {
			fprintf(stderr, "# Invalid option:\t%s\n", argv[args]);
			fprintf(stderr, "# Printing help message:\n");
			helpMessage(-1);
		}
		args++;
	}
	if(templatefilename == 0) {
		fprintf(stderr, "# Too few arguments handed\n");
		fprintf(stderr, "# Printing help message:\n");
		helpMessage(-1);
	}
	
	file_len = strlen(templatefilename);
	
	templatefilename = realloc(templatefilename, (file_len + strlen(".sparse.decontaminated.b") + 3) * sizeof(char));
	templates = malloc(sizeof(struct hashMap_shm));
	if(!templatefilename || !templates) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	if(sparse_run) {
		strcpy((templatefilename + file_len), ".sparse");
		file_len += strlen(".sparse");
		templatefilename[file_len] = '\0';
	}
	if(deCon) {
		strcpy((templatefilename + file_len), ".decontaminated");
		file_len += strlen(".decontaminated");
		templatefilename[file_len] = '\0';
	}
	
	if(destroy) {
		/* load sharing info */
		templatefilename[file_len] = '\0';
		strcat(templatefilename, ".m");
		file = fopen(templatefilename, "rb");
		if(file) {
			fprintf(stderr, "# Destroying shared DB: %s\n", templatefilename);
			t0 = clock();
			hashMap_shm_destroy(file, templatefilename);
			fclose(file);
			remove(templatefilename);
			t1 = clock();
			fprintf(stderr, "# DB destroyed in: %d s.\n", (int) difftime(t1, t0) / 1000000);
		} else {
			fprintf(stderr, "# Shared DB: %s, is already destroyed.\n", templatefilename);
		}
	} else {
		strcat(templatefilename, ".b");
		fprintf(stderr, "# Creating shared DB for: %s\n", templatefilename);
		t0 = clock();
		if(IOL) {
			dump_info = LoadConversion(templates, templatefilename);
		} else {
			dump_info = hashMap_shm_fill(templates, templatefilename);
		}
		/* dump sharing info */
		templatefilename[file_len] = '\0';
		strcat(templatefilename, ".m");
		file = fopen(templatefilename, "wb");
		if(!file) {
			fprintf(stderr, "File corruption: %s\n", templatefilename);
			exit(1);
		}
		hashMap_shm_dump(dump_info, file);
		fclose(file);
		templatefilename[file_len] = '\0';
		t1 = clock();
		fprintf(stderr, "# Shared DB created in: %d s.\n", (int) difftime(t1, t0) / 1000000);
		
		/* dump conversion */
		if(IOD) {
			strcat(templatefilename, ".mb");
			file = fopen(templatefilename, "wb");
			if(!file) {
				fprintf(stderr, "File corruption: %s\n", templatefilename);
				exit(1);
			}
			DumpConversion(templates, dump_info.valueSize, file);
			fclose(file);
			templatefilename[file_len] = '\0';
		}
		
		/* detach the DB */
		hashMap_shm_detach(templates);
	}
	
	/* Set sys shm to 1 GB
	* sudo sysctl -w kern.sysv.shmmax=1073741824
	* sudo sysctl -w kern.sysv.shmall=1073741824
	*
	* Set sys shm to 2 GB
	* sudo sysctl -w kern.sysv.shmmax=2147483648
	* sudo sysctl -w kern.sysv.shmall=2147483648
	*
	* Set sys shm to 3 GB
	* sudo sysctl -w kern.sysv.shmmax=3221225472
	* sudo sysctl -w kern.sysv.shmall=3221225472
	*
	* check status:
	* ipcs -a
	*/
	
	return 0;
}
