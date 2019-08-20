#ifndef READ_SPLIT
#define READ_SPLIT

#include "nary_tree.h"
#include "read.h"

// functions
int sequenceDetectSplitGeneome(trie *T, char* genomefname, char* pam, FILE *outputFile, faread_struct *fas, int printGRNAsOnly, int splitGenome);


#endif