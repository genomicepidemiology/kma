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
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <unistd.h>
#include <errno.h>

/*
 STRUCTURES
*/
struct hashMapKMA {
	/* end product of script */
	unsigned kmersize;		// k
	unsigned size;			// size of DB
	unsigned n;				// k-mers stored
	unsigned null_index;	// null value
	unsigned seqsize;		// size of seq
	unsigned v_index;		// size of values
	unsigned prefix_len;	// prefix length
	long unsigned *prefix;	// prefix
	unsigned *exist;		// size long
	long unsigned *seq;		// compressed sequence of k-mers
	unsigned *values;		// compressed values
	unsigned *key_index	;	// Relative
	unsigned *value_index;	// Relative
};

/*
 FUNCTIONS
*/
void hashMap_shm_detach(struct hashMapKMA *dest) {
	shmdt(dest->exist);
	shmdt(dest->seq);
	shmdt(dest->values);
	shmdt(dest->key_index);
	shmdt(dest->value_index);
}

void hashMapKMA_setupSHM(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid;
	long unsigned mask;
	key_t key;
	
	/* load sizes */
	fseek(file, sizeof(int), SEEK_CUR);
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
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (dest->kmersize << 1));
	
	if((dest->size - 1) == mask) {
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(unsigned), IPC_CREAT | 0666);
		if(shmid < 0) {
			fprintf(stderr, "Could not setup the shared hashMap v\n");
			fseek(file, dest->v_index * sizeof(unsigned), SEEK_CUR);
			dest->values = 0;
		} else {
			/* found */
			dest->values = shmat(shmid, NULL, 0);
			fread(dest->values, sizeof(unsigned), dest->v_index, file);
		}
	} else {
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
		shmid = shmget(key, dest->v_index * sizeof(unsigned), IPC_CREAT | 0666);
		if(shmid < 0) {
			fprintf(stderr, "Could not setup the shared hashMap v\n");
			fseek(file, dest->v_index * sizeof(unsigned), SEEK_CUR);
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
}

void hashMapKMA_destroySHM(struct hashMapKMA *dest, FILE *file, const char *filename) {
	
	int shmid;
	long unsigned mask;
	key_t key;
	
	/* load sizes */
	fseek(file, sizeof(int), SEEK_CUR);
	fread(&dest->kmersize, sizeof(unsigned), 1, file);
	fread(&dest->prefix_len, sizeof(unsigned), 1, file);
	fread(&dest->prefix, sizeof(long unsigned), 1, file);
	fread(&dest->size, sizeof(long unsigned), 1, file);
	fread(&dest->n, sizeof(unsigned), 1, file);
	fread(&dest->seqsize, sizeof(unsigned), 1, file);
	fread(&dest->v_index, sizeof(unsigned), 1, file);
	fread(&dest->null_index, sizeof(unsigned), 1, file);
	
	/* check shared memory, and destroy */
	key = ftok(filename, 'e');
	shmid = shmget(key, dest->size * sizeof(unsigned), 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * sizeof(long unsigned) - (dest->kmersize << 1));
	
	if((dest->size - 1) == mask) {
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(unsigned), 0666);
		if(shmid >= 0) {
			shmctl(shmid, IPC_RMID, NULL);
		}
	} else {
		key = ftok(filename, 's');
		shmid = shmget(key, dest->seqsize * sizeof(long unsigned), 0666);
		if(shmid >= 0) {
			shmctl(shmid, IPC_RMID, NULL);
		}
		key = ftok(filename, 'v');
		shmid = shmget(key, dest->v_index * sizeof(unsigned), 0666);
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
}

int * length_setupSHM(int *template_lengths, FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	int DB_size;
	fseek(file, 0, SEEK_END);
	size = ftell(file) - sizeof(int);
	fseek(file, sizeof(int), SEEK_SET);
	
	key = ftok(filename, 'l');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		template_lengths = 0;
	} else {
		template_lengths = shmat(shmid, NULL, 0);
		fread(template_lengths, sizeof(unsigned), size / sizeof(unsigned), file);
	}
	
	return template_lengths;
}

void length_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file) - sizeof(int);
	
	key = ftok(filename, 'l');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
	
}

int * index_setupSHM(int *index, FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file) - sizeof(int);
	fseek(file, sizeof(int), SEEK_SET);
	
	key = ftok(filename, 'i');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		index = 0;
	} else {
		index = shmat(shmid, NULL, 0);
		fread(index, sizeof(int), size / sizeof(int), file);
	}
	
	return index;
}

void index_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file) - sizeof(int);
	
	key = ftok(filename, 'i');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
}

long unsigned * seq_setupSHM(long unsigned *seq, FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file);
	rewind(file);
	
	key = ftok(filename, 's');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		seq = 0;
	} else {
		seq = shmat(shmid, NULL, 0);
		fread(seq, sizeof(long unsigned), size / sizeof(long unsigned), file);
	}
	
	return seq;
}

void seq_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file);
	
	key = ftok(filename, 's');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
}

char * name_setupSHM(char *template_names, FILE *file, const char *filename) {
	
	int i, shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file);
	rewind(file);
	
	key = ftok(filename, 'n');
	shmid = shmget(key, size, IPC_CREAT | 0666);
	if(shmid < 0) {
		fprintf(stderr, "Could not setup the shared length\n");
		template_names = 0;
	} else {
		template_names = shmat(shmid, NULL, 0);
		fread(template_names, 1, size, file);
		for(i = 0; i < size; i++) {
			if(template_names[i] == '\n') {
				template_names[i] = 0;
			}
		}
	}
	
	return template_names;
}

void name_destroySHM(FILE *file, const char *filename) {
	
	int shmid;
	long unsigned size;
	key_t key;
	
	/* load size */
	fseek(file, 0, SEEK_END);
	size = ftell(file);
	
	key = ftok(filename, 'n');
	shmid = shmget(key, size, 0666);
	if(shmid >= 0) {
		shmctl(shmid, IPC_RMID, NULL);
	}
}

void helpMessage(int exeStatus) {
	FILE *helpOut;
	if(exeStatus == 0) {
		helpOut = stdout;
	} else {
		helpOut = stderr;
	}
	fprintf(helpOut, "# kma_shm sets up a shared database (sysV) for mapping with KMA.\n");
	fprintf(helpOut, "# Options are:\t\tDesc:\t\t\t\tDefault:\tRequirements:\n");
	fprintf(helpOut, "#\n");
	fprintf(helpOut, "#\t-t_db\t\tTemplate DB\t\t\tNone\t\tREQUIRED\n");
	fprintf(helpOut, "#\t-destroy\tDestroy shared DB\t\tFalse\n");
	fprintf(helpOut, "#\t-shmLvl\t\tLevel of shared memory\t\t1\n");
	fprintf(helpOut, "#\t-shm-h\t\tExplain shm levels\n");
	fprintf(helpOut, "#\t-h\t\tShows this help message\n");
	fprintf(helpOut, "#\n");
	exit(exeStatus);
}

int main(int argc, char *argv[]) {
	
	int args, file_len, destroy, *template_lengths, *index;
	unsigned shmLvl;
	long unsigned *seq;
	char *templatefilename, *template_names;
	struct hashMapKMA *templates;
	time_t t0, t1;
	FILE *file;
	
	/* SET DEFAULTS */
	templatefilename = 0;
	destroy = 0;
	shmLvl = 1;
	
	/* PARSE COMMAND LINE OPTIONS */
	args = 1;
	while(args < argc) {
		if(strcmp(argv[args], "-t_db") == 0) {
			args++;
			if(args < argc) {
				templatefilename = malloc(strlen(argv[args]) + 64);
				if(!templatefilename) {
					fprintf(stderr, "OOM\n");
					exit(1);
				}
				strcpy(templatefilename, argv[args]);
			}
		} else if(strcmp(argv[args], "-destroy") == 0) {
			destroy = 1;
		} else if(strcmp(argv[args], "-shmLvl") == 0) {
			args++;
			if(args < argc) {
				shmLvl = atoi(argv[args]);
				if(!shmLvl) {
					fprintf(stderr, "Invalid shmLvl\n");
					exit(0);
				}
			}
		} else if(strcmp(argv[args], "-h") == 0) {
			helpMessage(0);
		} else if(strcmp(argv[args], "-shm-h") == 0) {
			fprintf(stderr, "# Flags for shared memory, add them to combine them\n");
			fprintf(stderr, "# After shm is setup, the DB-files should not be changed.\n");
			fprintf(stderr, "#\n");
			fprintf(stderr, "#\tDB piece\t\tFlag\n");
			fprintf(stderr, "#\t*.comp.b\t\t1\n");
			fprintf(stderr, "#\t*.decon.comp.b\t\t2\n");
			fprintf(stderr, "#\t*.length.b\t\t4\n");
			fprintf(stderr, "#\t*.seq.b *.index.b\t8\n");
			fprintf(stderr, "#\t*.name\t\t\t16\n");
			//fprintf(stderr, "#\tall\t\t\t31\n");
			exit(0);
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
	templates = malloc(sizeof(struct hashMapKMA));
	if(!templates) {
		fprintf(stderr, "OOM\n");
		exit(1);
	}
	
	/* setup or destroy shm */
	if(destroy) {
		/* *comp.b */
		if(shmLvl & 1) {
			strcat(templatefilename, ".comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				hashMapKMA_destroySHM(templates, file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *decon.comp.b */
		if(shmLvl & 2) {
			strcat(templatefilename, ".decon.comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				hashMapKMA_destroySHM(templates, file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.length.b */
		if(shmLvl & 4) {
			strcat(templatefilename, ".length.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				length_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.seq.b *.index.b */
		if(shmLvl & 8) {
			/* *.seq.b */
			strcat(templatefilename, ".seq.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				seq_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
			
			/* *.index.b */
			strcat(templatefilename, ".index.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				index_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.name */
		if(shmLvl & 16) {
			strcat(templatefilename, ".name");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				name_destroySHM(file, templatefilename);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
	} else {
		/* *.comp.b */
		if(shmLvl & 1) {
			strcat(templatefilename, ".comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				hashMapKMA_setupSHM(templates, file, templatefilename);
				hashMap_shm_detach(templates);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.decon.comp.b */
		if(shmLvl & 2) {
			strcat(templatefilename, ".decon.comp.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				hashMapKMA_setupSHM(templates, file, templatefilename);
				hashMap_shm_detach(templates);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.length.b */
		if(shmLvl & 4) {
			strcat(templatefilename, ".length.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				template_lengths = length_setupSHM(template_lengths, file, templatefilename);
				if(template_lengths)
					shmdt(template_lengths);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.seq.b *.index.b */
		if(shmLvl & 8) {
			/* *.seq.b */
			strcat(templatefilename, ".seq.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				seq = seq_setupSHM(seq, file, templatefilename);
				if(seq)
					shmdt(seq);
				fclose(file);
			}
			templatefilename[file_len] = 0;
			
			/* *.index.b */
			strcat(templatefilename, ".index.b");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				index = index_setupSHM(index, file, templatefilename);
				if(index)
					shmdt(index);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
		
		/* *.name */
		if(shmLvl & 16) {
			strcat(templatefilename, ".name");
			file = fopen(templatefilename, "rb");
			if(!file) {
				fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			} else {
				template_names = name_setupSHM(template_names, file, templatefilename);
				if(template_names)
					shmdt(template_names);
				fclose(file);
			}
			templatefilename[file_len] = 0;
		}
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
