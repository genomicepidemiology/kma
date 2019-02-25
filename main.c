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
#define _XOPEN_SOURCE 600
#include <string.h>
#include "kma.h"
#include "index.h"
#include "shm.h"
#include "seq2fasta.h"
#include "update.h"

static int helpmessage() {
	
	fprintf(stderr, "# For help please use:\n");
	fprintf(stderr, "#\tkma -h\n");
	fprintf(stderr, "#\tkma index -h\n");
	fprintf(stderr, "#\tkma shm -h\n");
	fprintf(stderr, "#\tkma seq2fasta -h\n");
	fprintf(stderr, "#\tkma update -h\n");
	return 1;
}

int main(int argc, char *argv[]) {
	
	int status;
	
	if(argc != 1) {
		if(**++argv == '-') {
			status = kma_main(argc, --argv);
		} else if(--argc && strcmp(*argv, "index") == 0) {
			status = index_main(argc, argv);
		} else if(strcmp(*argv, "shm") == 0) {
			status = shm_main(argc, argv);
		} else if(strcmp(*argv, "seq2fasta") == 0) {
			status = seq2fasta_main(argc, argv);
		} else if(strcmp(*argv, "update") == 0) {
			status = update_main(argc, argv);
		} else {
			fprintf(stderr, "Invalid option:\t%s\n", *argv);
			status = helpmessage();
		}
	} else {
		fprintf(stderr, "Too few arguments handed\n");
		status = helpmessage();
	}
	
	return status;
}
