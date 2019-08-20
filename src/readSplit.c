#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <ctype.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <ctype.h>
#include <assert.h>

#include "read.h"
#include "nary_tree.h"
#include "vcf.h"
#include "readSplit.h"


int sequenceDetectSplitGeneome(trie *T, char* genomefname, char* pam, FILE *outputFile, faread_struct *fas, int printGRNAsOnly, int splitGenome) {
    /* Identifies all gRNAs in reference sequence and loads them to Trie. */
    grna_list *g = NULL;
    int chr_number;
    int guidelen = PROTOSPACER_LENGTH + strlen(pam);
    FILE *outfh = outputFile;

    int fastaReaderOutput = 0;
    // First check if number of headers read is >= splitGenome to prevent fastaReader(fas) from
    // running in that case
    while (splitGenome && (fastaReaderOutput = fastaReader(fas)))
    {
        fprintf(stdout,"[crisflash] Processing %s (length %lld) in %s ...\n", fas->header, fas->slen, genomefname);
        fastaReaderImproveSequence(fas);
        // faSequenceToTrie(T, fas->header, fas->s, fas->sr, fas->slen, pam, outputfname, uppercaseOnly);
        g = fastaSequenceToGRNAs(fas->header, fas->s, fas->sr, fas->slen, pam);
        if (printGRNAsOnly)
        {
        printGRNAs(g, guidelen, outfh);
        }
        else
        {
        chr_number = addChr(T, g->chr, strlen(g->chr)); // we get the chr number in order of appearence in fasta file
        GRNAsToTrie(T, g, guidelen, chr_number);
        }
        freeGRNAs(g);

        splitGenome--;
    }

    return fastaReaderOutput;
}
