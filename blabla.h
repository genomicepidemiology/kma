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

#include "stddef.h"

#ifndef BLABLA
typedef struct bla Bla;
struct bla {
	unsigned read; /* reader */
	unsigned write; /* writer */
	unsigned size;
	volatile unsigned lock;
	void *buff;
};
#define BLABLA 1
#endif

Bla * blablastart(size_t size);
void blablabreak(Bla *stream);
void blablastop(Bla *stream);
int blablawrite(void *msg, size_t size, Bla *stream);
int blablaread(void *msg, size_t size, Bla *stream);
unsigned blablaseek(Bla *stream, size_t size);
