#pragma once

/* vector<Moo> goo;
 * foreach(g, goo) {
 *   if (g->pritzical)
 *     g->frobnicate();
 * }
 */
#define foreach(I,VECTOR) \
  for (decltype(&*VECTOR.begin()) _i=0, I=&(VECTOR)[0]; \
       ((unsigned long int)_i) < (VECTOR).size(); \
       (*((unsigned long int*)&_i))++,I=&(VECTOR)[(unsigned long int)_i])

