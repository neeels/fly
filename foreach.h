#pragma once

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(*X))

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

