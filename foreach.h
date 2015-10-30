#pragma once

/*  vector<Moo> goo;
 * or:
 *  list<Moo> goo;
 * foreach(g, goo) {
 *   if (g->pritzical)
 *     g->frobnicate();
 * }
 */
#define foreach(I,LIST) \
  for (decltype(LIST.begin()) I=LIST.begin(); \
       I != (LIST).end(); \
       ++I)

