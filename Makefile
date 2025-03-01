CC=gcc
RM=rm
MAKE=make
CFLAGS=-std=c99 -O3 -pthread
all:
	$(MAKE) --no-print-directory ./bin/crisflash
	$(MAKE) --no-print-directory crisflashVcf

./bin/crisflash: ./src/main_crisflash.c ./src/read.c ./src/nary_tree.c ./src/vcf.c ./src/readSplit.c
	mkdir -p bin
	$(CC) $(CFLAGS) ./src/main_crisflash.c ./src/read.c ./src/nary_tree.c ./src/vcf.c ./src/readSplit.c -o ./bin/crisflash

crisflashVcf: ./src/main_vcf.c ./src/vcf.c ./src/nary_tree.c ./src/read.c ./src/readSplit.c
	$(CC) $(CFLAGS) ./src/main_vcf.c ./src/vcf.c ./src/read.c ./src/nary_tree.c ./src/readSplit.c -o ./bin/crisflashVcf

clean: ./bin/crisflash
	$(RM) ./bin/crisflash
	$(RM) ./bin/crisflashVcf
