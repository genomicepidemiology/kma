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

#define getNuc(Comp,pos) ((Comp[pos >> 5] << ((pos & 31) << 1)) >> 62)
#define setEx(src, pos)(src[pos >> 3] |= (1 << (pos & 7)))
#define unsetEx(src, pos)(src[pos >> 3] ^= (1 << (pos & 7)))
#define getEx(src, pos)((src[pos >> 3] >> (pos & 7)) & 1)
#define getKmer_macro(kmer, Comp, pos, cPos, iPos, shifter) \
		iPos = (pos & 31) << 1;\
		cPos = pos >> 5;\
		kmer = (iPos <= shifter) ? ((Comp[cPos] << iPos) >> shifter) : (((Comp[cPos] << iPos) | (Comp[cPos + 1] >> (64-iPos))) >> shifter);

long unsigned getKmer(long unsigned *compressor, unsigned cPos, const unsigned shifter);
long unsigned makeKmer(const unsigned char *qseq, unsigned pos, unsigned size);
int charpos(const unsigned char *src, unsigned char target, int start, int len);
void strrc(unsigned char *qseq, int q_len);
void strtranslate(unsigned char *strp, char *trans);
void nibble2base(unsigned char *seq, int len);
