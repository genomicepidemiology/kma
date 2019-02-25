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
#include <stdlib.h>
#include <stdio.h>
#include "kmapipe.h"
#include "pherror.h"
#ifdef _WIN32
#define _PATH_BSHELL "/bin/sh" 
#include <windows.h>
#else
#include <paths.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

FILE * kmaPipe(const char *cmd, const char *type, FILE *ioStream, int *status) {
	
	/* kmaPipe is a combination of popen and pclose, but allows for binary mode */
	static Pid *pidlist;
	int pdes[2];
	pid_t pid;
	Pid *src, *last;
	Pid * volatile dest;
	char *argv[] = {"sh", "-c", NULL, NULL};
	
	if(cmd && type) {
		/* check mode */
		if(*type != 'r' && *type != 'w') {
			errno = EINVAL;
			ERROR();
		}
		
		/* create pipe */
		dest = malloc(sizeof(Pid));
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
			Pid *pidPtr;
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
	} else {
		*status = 0;
		/* Get file pointer. */
		for (last = 0, src = pidlist; src->fp != ioStream; last = src, src = src->next) {
			if(!src) {
				*status = 1;
				return 0;
			}
		}
		
		/* close stream and get exit status */
		fclose(ioStream);
		*status = 1;
		#ifndef _WIN32
		while ((pid = waitpid(src->pid, status, 0)) == -1 && errno == EINTR) {
			usleep(100);
		}
		*status = 0;
		#else
		WaitForSingleObject(src->pid, INFINITE);
		#endif
		
		/* Remove the entry from the linked list. */
		if (!last) {
			pidlist = src->next;
		} else {
			last->next = src->next;
		}
		free(src);
		
		return 0;
	}
}
