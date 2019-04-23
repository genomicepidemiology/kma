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

#include <stdlib.h>
#include <string.h>
#include "blabla.h"
#include "pherror.h"
#include "stddef.h"
#include "threader.h"

Bla * blablastart(size_t size) {
	
	Bla *dest;
	
	dest = smalloc(sizeof(Bla));
	dest->read = 0;
	dest->write = 0;
	dest->size = size;
	dest->lock = 0;
	dest->buff = smalloc(size);
	
	return dest;
}

void blablabreak(Bla *stream) {
	
	stream->write = -1;
}

void blablastop(Bla *stream) {
	
	free(stream->buff);
	stream->buff = 0;
}

int blablawrite(void *msg, size_t size, Bla *stream) {
	
	volatile unsigned *lock, wait;
	
	if(stream->buff == 0) {
		return 0;
	} else if(size == 0) {
		return 1;
	}
	
	lock = &stream->lock;
	if(stream->write == stream->size) {
		lock(lock);
		wait = stream->read;
		while(wait == 0) {
			unlock(lock);
			wait = 0;
			while(wait == stream->read) {
				usleep(100);
			}
			lock(lock);
			wait = stream->read;
		}
		stream->write = 0;
		unlock(lock);
	}
	
	if(stream->write + size <= stream->size) { /* one cpy */
		/* check if reader is in the way */
		lock(lock);
		wait = stream->read;
		while(stream->write < wait && wait <= stream->write + size) {
			unlock(lock);
			while(wait == stream->read) {
				usleep(100);
			}
			lock(lock);
			wait = stream->read;
		}
		unlock(lock);
		
		memcpy(stream->buff + stream->write, msg, size);
		lock(lock);
		stream->write += size;
		unlock(lock);
	} else { /* splt cpy */
		/* check if reader is in the way */
		lock(lock);
		wait = stream->read;
		while(stream->write < wait) {
			unlock(lock);
			while(wait == stream->read) {
				usleep(100);
			}
			lock(lock);
			wait = stream->read;
		}
		unlock(lock);
		
		/* cpy first piece */
		memcpy(stream->buff + stream->write, msg, stream->size - stream->write);
		size = stream->write + size - stream->size;
		
		/* check if piece fits */
		if(size <= stream->size) {
			lock(lock);
			wait = 0;
			while(wait < stream->read && stream->read <= wait + size) {
				wait = stream->read;
				unlock(lock);
				while(wait == stream->read) {
					usleep(100);
				}
				lock(lock);
				wait = 0;
			}
			unlock(lock);
			
			memcpy(stream->buff, msg, size);
		} else {
			/* realloc buff */
			lock(lock);
			stream->write = 0;
			wait = 0;
			while(wait != stream->read) {
				unlock(lock);
				while(wait != stream->read) {
					usleep(100);
				}
				lock(lock);
			}
			unlock(lock);
			
			free(stream->buff);
			stream->buff = smalloc(size);
			stream->size = size;
			memcpy(stream->buff, msg, size);
		}
		lock(lock);
		stream->write = size;
		unlock(lock);
	}
	
	return 1;
}

int blablaread(void *msg, size_t size, Bla *stream) {
	
	volatile unsigned *lock, wait;
	
	lock = &stream->lock;
	if(stream->write < 0) {
		return 0;
	} else if(size == 0) {
		return 1;
	} else if(stream->read + size <= stream->size) { /* one cpy */
		/* check if writer is in the way */
		lock(lock);
		wait = stream->write;
		while(stream->read <= wait && wait <= stream->read + size) {
			unlock(lock);
			while(wait == stream->write) {
				usleep(100);
			}
			lock(lock);
			wait = stream->write;
		}
		unlock(lock);
		
		memcpy(msg, stream->buff + stream->read, size);
		lock(lock);
		if((stream->read += size) == stream->size) {
			stream->read = 0;
		}
		unlock(lock);
	} else { /* splt cpy */
		/* check if reader is in the way */
		lock(lock);
		wait = stream->write;
		while(stream->read <= wait && wait != stream->size) {
			unlock(lock);
			while(wait == stream->write) {
				usleep(100);
			}
			lock(lock);
			wait = stream->write;
		}
		unlock(lock);
		
		/* cpy first piece */
		memcpy(msg, stream->buff + stream->read, stream->size - stream->read);
		size = stream->read + size - stream->size;
		
		/* check if piece fits */
		if(size < stream->size) {
			lock(lock);
			stream->read = 0;
			wait = stream->write;
			while(0 <= wait && wait <= size) {
				unlock(lock);
				while(wait == stream->write) {
					usleep(100);
				}
				lock(lock);
				wait = stream->write;
			}
			unlock(lock);
			
			memcpy(stream->buff, msg, size);
			lock(lock);
			stream->read = size;
			unlock(lock);
		} else {
			lock(lock);
			stream->read = 0;
			unlock(lock);
			blablaread(msg, size, stream);
		}
	}
		return 1;
}

unsigned blablaseek(Bla *stream, size_t size) {
	
	if(stream->size < (stream->read += size)) {
		stream->read -=  stream->size;
	}
	
	return stream->read;
}
