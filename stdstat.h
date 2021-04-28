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

#define MIN(X, Y) ((X < Y) ? X : Y)
#define MAX(X, Y) ((X < Y) ? Y : X)

extern int (*cmp)(int, int);
int cmp_or(int t, int q);
int cmp_and(int t, int q);
int cmp_true(int t, int q);
double fastp(long double q);
double p_chisqr(long double q);
double power(double x, unsigned n);
double binP(int n, int k, double p);
unsigned minimum(unsigned *src, unsigned n);
double eQual(unsigned char *qual, const int len, const int minQ, const double *prob);
