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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "compress.h"
#include "hashmap.h"
#include "hashmapkma.h"
#include "pherror.h"
#include "valueshash.h"

HashMapKMA * compressKMA_DB(HashMap *templates, FILE *out) {
	
	long unsigned i, j, check;
	long unsigned index, t_index, v_index, new_index, null_index;
	unsigned *values;
	short unsigned *values_s;
	HashMapKMA *finalDB;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	HashTable *node_t, *next_t, *table_t;
	
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
	finalDB = smalloc(sizeof(HashMapKMA));
	/* Fill in known values */
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->mask = templates->mask;
	finalDB->prefix_len = templates->prefix_len;
	finalDB->prefix = templates->prefix;
	finalDB->kmersize = templates->kmersize;
	finalDB->DB_size = templates->DB_size;
	
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
	
	if(finalDB->kmersize <= 16) {
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
	shmValues = initialize_hashValues(null_index, finalDB->DB_size);
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
			new_index = valuesHash_add(shmValues, node_t->value, v_index);
			
			if(new_index == v_index) {
				v_index += valuesSize(node_t->value);
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
				free(node_t->value);
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
	
	if(finalDB->DB_size < USHRT_MAX) {
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

HashMapKMA * compressKMA_megaDB(HashMap *templates, FILE *out) {
	
	long unsigned i, j, v_index, new_index, null_index;
	unsigned check, *values;
	short unsigned *values_s;
	HashMapKMA *finalDB;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	
	/* Fill in known values */
	hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	check = 0;
	check = ~check;
	finalDB = smalloc(sizeof(HashMapKMA));
	finalDB->size = templates->size + 1;
	finalDB->n = templates->n;
	finalDB->mask = templates->mask;
	finalDB->prefix_len = templates->prefix_len;
	finalDB->prefix = templates->prefix;
	finalDB->kmersize = templates->kmersize;
	finalDB->DB_size = templates->DB_size;
	
	/* allocate existence */
	finalDB->exist = smalloc(finalDB->size * sizeof(unsigned));
	finalDB->exist_l = 0;
	finalDB->key_index = 0;
	finalDB->value_index = 0;
	
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
	shmValues = initialize_hashValues(null_index, finalDB->DB_size);
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
	
	if(finalDB->DB_size < USHRT_MAX) {
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

void compressKMA_deconDB(HashMapKMA *finalDB, unsigned **Values) {
	
	long unsigned i, j, v_index, new_index, check;
	unsigned *values;
	short unsigned *values_s;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	
	/* prepare final DB */
	check = 0;
	check = ~check;
	if(finalDB->v_index < UINT_MAX) {
		check >>= 32;
		hashMapKMA_addValue_ptr = &hashMapKMA_addValue;
		getValueIndexPtr = &getValueIndex;
	} else {
		hashMapKMA_addValue_ptr = &hashMapKMA_addValueL;
		getValueIndexPtr = &getValueIndexL;
	}
	i = finalDB->n;
	shmValues = initialize_hashValues(finalDB->n, finalDB->DB_size);
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
	if(finalDB->DB_size < USHRT_MAX) {
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

void compressKMA_deconMegaDB(HashMapKMA *finalDB, unsigned **Values) {
	
	long unsigned i, j, v_index, new_index, pos, check;
	unsigned *values;
	short unsigned *values_s;
	ValuesHash *shmValues;
	ValuesTable *node, *next, *table;
	
	fprintf(stderr, "# Compressing indexes.\n");
	check = 0;
	check = ~check;
	if(finalDB->v_index < UINT_MAX) {
		check >>= 32;
		getExistPtr = &getExist;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExist;
	} else {
		getExistPtr = &getExistL;
		hashMapKMA_addExist_ptr = &hashMapKMA_addExistL;
	}
	i = finalDB->size;
	shmValues = initialize_hashValues(finalDB->n, finalDB->DB_size);
	v_index = 0;
	while(i--) {
		if((pos = getExistPtr(finalDB->exist, i)) != finalDB->n) {
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
	if(finalDB->DB_size < USHRT_MAX) {
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
