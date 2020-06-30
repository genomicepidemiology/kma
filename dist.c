/* Philip T.L.C. Clausen Jan 2020 plan@dtu.dk */

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

#define _XOPEN_SOURCE 600
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashmapkma.h"
#include "matrix.h"
#include "pherror.h"
#include "runkma.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

HashMapKMA * loadValues(const char *filename) {
	
	long unsigned check, size;
	FILE *file;
	HashMapKMA *dest;
	
	/* init */
	file = sfopen(filename, "rb");
	dest = smalloc(sizeof(HashMapKMA));
	
	/* load sizes */
	sfread(&dest->DB_size, sizeof(unsigned), 1, file);
	sfread(&dest->kmersize, sizeof(unsigned), 1, file);
	sfread(&dest->prefix_len, sizeof(unsigned), 1, file);
	sfread(&dest->prefix, sizeof(long unsigned), 1, file);
	sfread(&dest->size, sizeof(long unsigned), 1, file);
	sfread(&dest->n, sizeof(long unsigned), 1, file);
	sfread(&dest->v_index, sizeof(long unsigned), 1, file);
	sfread(&dest->null_index, sizeof(long unsigned), 1, file);
	
	dest->mask = 0;
	dest->mask = (~dest->mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (dest->kmersize << 1));
	dest->shmFlag = 0;
	
	/* simple check for old indexing */
	if(dest->size < dest->n) {
		free(dest);
		fclose(file);
		return 0;
	}
	
	/* exist */
	size = dest->size;
	if((dest->size - 1) == dest->mask) {
		/* mega */
		if(dest->v_index <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		
		/* load */
		dest->exist = smalloc(size);
		check = fread(dest->exist, 1, size, file);
		if(check != size) {
			free(dest);
			free(dest->exist);
			fclose(file);
			return 0;
		}
		dest->exist_l = (long unsigned *)(dest->exist);
		dest->shmFlag |= 1;
	} else {
		if(dest->n <= UINT_MAX) {
			size *= sizeof(unsigned);
		} else {
			size *= sizeof(long unsigned);
		}
		
		/* skip */
		dest->exist = 0;
		dest->exist_l = 0;
		fseek(file, size, SEEK_CUR);
	}
	
	/* values */
	size = dest->v_index;
	if(dest->DB_size < USHRT_MAX) {
		size *= sizeof(short unsigned);
	} else {
		size *= sizeof(unsigned);
	}
	dest->values = smalloc(size);
	check = fread(dest->values, 1, size, file);
	if(check != size) {
		free(dest);
		free(dest->exist);
		free(dest->values);
		fclose(file);
		return 0;
	}
	dest->values_s = (short unsigned *)(dest->values);
	dest->shmFlag |= 2;
	
	/* check for megaMap */
	if(dest->exist) {
		fclose(file);
		return dest;
	}
	
	/* skip kmers */
	size = dest->n + 1;
	if(dest->kmersize <= 16) {
		size *= sizeof(unsigned);
	} else {
		size *= sizeof(long unsigned);
	}
	dest->key_index = 0;
	dest->key_index_l = 0;
	fseek(file, size, SEEK_CUR);
	
	/* value indexes */
	size = dest->n;
	if(dest->v_index < UINT_MAX) {
		size *= sizeof(unsigned);
	} else {
		size *= sizeof(long unsigned);
	}
	dest->exist = smalloc(size);
	check = fread(dest->exist, 1, size, file);
	if(check != size) {
		free(dest);
		free(dest->exist);
		free(dest->values);
		fclose(file);
		return 0;
	}
	dest->exist_l = (long unsigned *)(dest->value_index);
	dest->shmFlag |= 1;
	
	fclose(file);
	return dest;
}

void destroyValues(HashMapKMA *src) {
	
	free(src->exist);
	free(src->values);
	free(src);
}

void kmerSimilarity(HashMapKMA *DB, Matrix *Dist, int *N) {
	
	int i, j, el, vs, **D;
	unsigned *exist, *values_i, *values_j;
	long unsigned n, pos, *exist_l;
	short unsigned *values_si, *values_sj;
	
	/* init */
	el = DB->v_index < UINT_MAX;
	vs = DB->DB_size < USHRT_MAX;
	exist = DB->exist - 1;
	exist_l = DB->exist_l - 1;
	D = Dist->mat;
	
	/* get values */
	n = DB->n;
	while(n) {
		pos = el ? *++exist : *++exist_l;
		if(pos != 1) {
			if(vs) {
				values_si = DB->values_s + pos;
				i = *values_si + 1;
				values_si += i;
				while(--i) {
					j = i;
					values_sj = --values_si;
					while(--j) {
						++D[*values_si - 1][*--values_sj - 1];
					}
					++N[*values_si - 1];
				}
			} else {
				values_i = DB->values + pos;
				i = *values_i + 1;
				values_i += i;
				while(--i) {
					j = i;
					values_j = --values_i;
					while(--j) {
						++D[*values_i - 1][*--values_j - 1];
					}
					++N[*values_i - 1];
				}
			}
			--n;
		}
	}
	
	Dist->n = DB->DB_size - 1;
}

void kmerDist(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format) {
	
	int i, j, *D;
	char *name;
	
	/* init */
	D = *(Dist->mat) - 1;
	
	if(format & 4) {
		fprintf(outfile, "#%s\n", "k-mer distance");
	}
	fprintf(outfile, "%10d\n", Dist->n);
	for(i = 0; i < Dist->n; ++i) {
		name = nameLoad(template_name, name_file);
		if(format & 1) {
			fprintf(outfile, "%s", name);
		} else {
			fprintf(outfile, "%-10.10s", name);
		}
		
		for(j = 0; j < i; ++j) {
			fprintf(outfile, "\t%d", N[i] + N[j] - 2 * *++D);
		}
		fprintf(outfile, "\n");
	}
	fseek(name_file, 0, SEEK_SET);
}

void kmerShared(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format) {
	
	int i, j, *D;
	char *name;
	
	/* init */
	D = *(Dist->mat) - 1;
	
	if(format & 4) {
		fprintf(outfile, "#%s\n", "shared k-mers");
	}
	fprintf(outfile, "%10d\n", Dist->n);
	for(i = 0; i < Dist->n; ++i) {
		name = nameLoad(template_name, name_file);
		if(format & 1) {
			fprintf(outfile, "%s", name);
		} else {
			fprintf(outfile, "%-10.10s", name);
		}
		
		for(j = 0; j < i; ++j) {
			fprintf(outfile, "\t%d", *++D);
		}
		fprintf(outfile, "\n");
	}
	fseek(name_file, 0, SEEK_SET);
}

void kmerQuery(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format) {
	
	int i, j, **D;
	double d;
	char *name;
	
	/* init */
	D = Dist->mat;
	
	if(format & 4) {
		fprintf(outfile, "#%s\n", "Avg. k-mer coverage");
	}
	fprintf(outfile, "%10d\n", Dist->n);
	for(i = 0; i < Dist->n; ++i) {
		name = nameLoad(template_name, name_file);
		if(format & 1) {
			fprintf(outfile, "%s", name);
		} else {
			fprintf(outfile, "%-10.10s", name);
		}
		
		for(j = 0; j < Dist->n - 1; ++j) {
			d = i < j ? D[j][i] : D[i][j];
			fprintf(outfile, "\t%.2f", 100.0 * d / N[i]);
		}
		fprintf(outfile, "\n");
	}
	fseek(name_file, 0, SEEK_SET);
}

void kmerTemplate(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format) {
	
	int i, j, **D;
	double d;
	char *name;
	
	/* init */
	D = Dist->mat;
	
	if(format & 4) {
		fprintf(outfile, "#%s\n", "Avg. k-mer coverage");
	}
	fprintf(outfile, "%10d\n", Dist->n);
	for(i = 0; i < Dist->n; ++i) {
		name = nameLoad(template_name, name_file);
		if(format & 1) {
			fprintf(outfile, "%s", name);
		} else {
			fprintf(outfile, "%-10.10s", name);
		}
		
		for(j = 0; j < Dist->n - 1; ++j) {
			d = i < j ? D[j][i] : D[i][j];
			fprintf(outfile, "\t%.2f", 100.0 * d / N[j]);
		}
		fprintf(outfile, "\n");
	}
	fseek(name_file, 0, SEEK_SET);
}

void kmerAvg(FILE *outfile, Matrix *Dist, int *N, FILE *name_file, Qseqs *template_name, unsigned format) {
	
	int i, j, *D;
	char *name;
	
	/* init */
	D = *(Dist->mat) - 1;
	
	if(format & 4) {
		fprintf(outfile, "#%s\n", "Avg. k-mer coverage");
	}
	fprintf(outfile, "%10d\n", Dist->n);
	for(i = 0; i < Dist->n; ++i) {
		name = nameLoad(template_name, name_file);
		if(format & 1) {
			fprintf(outfile, "%s", name);
		} else {
			fprintf(outfile, "%-10.10s", name);
		}
		
		for(j = 0; j < i; ++j) {
			fprintf(outfile, "\t%.2f", 200.0 * *++D / (N[i] + N[j]));
		}
		fprintf(outfile, "\n");
	}
	fseek(name_file, 0, SEEK_SET);
}

void runDist(char *templatefilename, char *outputfilename, int flag, int format) {
	
	int file_len, *N;
	unsigned DB_size;
	FILE *outfile, *name_file;
	HashMapKMA *DB;
	Matrix *Dist;
	Qseqs *template_name;
	
	/* init */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	file_len = strlen(templatefilename);
	strcpy(templatefilename + file_len, ".comp.b");
	if(!(DB = loadValues(templatefilename))) {
		fprintf(stderr, "Wrong format of DB.\n");
		exit(1);
	}
	templatefilename[file_len] = 0;
	DB_size = DB->DB_size;
	N = calloc(DB_size, sizeof(int));
	if(!N) {
		ERROR();
	}
	Dist = ltdMatrix_init(DB_size);
	
	/* get kmer similarities and lengths */
	kmerSimilarity(DB, Dist, N);
	destroyValues(DB);
	
	/* load names */
	strcpy(templatefilename + file_len, ".name");
	name_file = sfopen(templatefilename, "rb");
	templatefilename[file_len] = 0;
	template_name = setQseqs(1024);
	
	/* k-mer dist, lt */
	if(flag & 1) {
		kmerDist(outfile, Dist, N, name_file, template_name, format);
	}
	
	/* k-mer shared, lt */
	if(flag & 2) {
		kmerShared(outfile, Dist, N, name_file, template_name, format);
	}
	
	/* query cov, asym */
	if(flag & 4) {
		kmerQuery(outfile, Dist, N, name_file, template_name, format);
	}
	
	/* template cov, asym */
	if(flag & 8) {
		kmerTemplate(outfile, Dist, N, name_file, template_name, format);
	}
	
	/* avg. cov, lt */
	if(flag & 16) {
		kmerAvg(outfile, Dist, N, name_file, template_name, format);
	}
	
	/* clean */
	fclose(outfile);
	fclose(name_file);
	free(N);
	Matrix_destroy(Dist);
	destroyQseqs(template_name);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#kma dist calculates distances between templates from a kma index\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t_db", "Template DB", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

/* main */
int dist_main(int argc, char *argv[]) {
	
	int args, flag, format, file_len;
	char *arg, *errorMsg, *templatefilename, *outputfilename;
	
	/* init */
	flag = 1;
	format = 1;
	templatefilename = 0;
	outputfilename = 0;
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			 if(strcmp(arg, "t_db") == 0) {
				if(++args < argc) {
					file_len = strlen(argv[args]);
					templatefilename = smalloc(file_len + 64);
					strcpy(templatefilename, argv[args]);
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					format = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-f\"");
					}
				} else {
					missArg("\"-f\"");
				}
			} else if(strcmp(arg, "fh") == 0) {
				fprintf(stdout, "# Format flags output, add them to combine them.\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#   1:\tRelaxed Phylip\n");
				fprintf(stdout, "#   4:\tInclude template name in phylip file\n");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					flag = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-f\"");
					}
				} else {
					missArg("\"-d\"");
				}
			} else if(strcmp(arg, "dh") == 0) {
				fprintf(stdout, "# Distance calculation methods, add them to combine them:\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#   1:\tk-mer distance\n");
				fprintf(stdout, "#   2:\tShared k-mers\n");
				fprintf(stdout, "#   4:\tk-mer query coverage\n");
				fprintf(stdout, "#   8:\tk-mer template coverage\n");
				fprintf(stdout, "#   16:\tk-mer avg. coverage\n");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(arg, "h") == 0) {
				return helpMessage(stdout);
			} else {
				fprintf(stderr, "Unknown option:%s\n", arg - 1);
				return helpMessage(stderr);
			}
		} else {
			fprintf(stderr, "Unknown argument:%s\n", arg - 1);
			return helpMessage(stderr);
		}
		++args;
	}
	
	if(!templatefilename) {
		fprintf(stderr, "Too few arguments handed.\n");
		exit(1);
	}
	if(!outputfilename) {
		outputfilename = "--";
	}
	
	runDist(templatefilename, outputfilename, flag, format);
	
	return 0;
}
