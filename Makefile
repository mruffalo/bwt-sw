DEFINE = 
CC = gcc
CFLAGS = -O3 -funroll-loops -maccumulate-outgoing-args -Wall

all:	BWTFormatdb BWTBlast BWTSW BlastCompare

BWTBlast.o:	BWTBlast.c MiscUtilities.h TypeNLimit.h MemManager.h BWT.h TextConverter.h HSP.h iniparser.h dictionary.h r250.h HSPstatistic.h Timing.h blast_dust.h Socket.h Makefile

BWT.o:	BWT.c BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h Makefile

BWTConstruct.o:	BWTConstruct.c BWTConstruct.h TypeNLimit.h BWT.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h QSufSort.h r250.h Makefile

BWTFormatdb.o:	BWTFormatdb.c BWTConstruct.h TypeNLimit.h BWT.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h iniparser.h dictionary.h BWTFormatdb.h Timing.h Makefile

BWTSW.o:     BWTSW.c MiscUtilities.h TypeNLimit.h MemManager.h BWT.h TextConverter.h HSP.h iniparser.h dictionary.h r250.h HSPstatistic.h Timing.h Socket.h Makefile

BlastCompare.o:	BlastCompare.c MiscUtilities.h TypeNLimit.h iniparser.h dictionary.h Makefile

dictionary.o:	dictionary.c dictionary.h inistrlib.h dictionary.h Makefile

blast_dust.o:	blast_dust.c blast_dust.h Makefile

DNACount.o:	DNACount.c DNACount.h TypeNLimit.h MiscUtilities.h Makefile

HSP.o:	HSP.c TextConverter.h TypeNLimit.h MemManager.h MiscUtilities.h r250.h HSP.h HSPstatistic.h Makefile

HSPstatistic.o:	HSPstatistic.c karlin.h HSPstatistic.h Makefile

iniparser.o:	iniparser.c iniparser.h dictionary.h inistrlib.h Makefile

inistrlib.o:	inistrlib.c inistrlib.h Makefile

karlin.o:	karlin.c karlin.h Makefile

MemManager.o:	MemManager.c MiscUtilities.h TypeNLimit.h MemManager.h Makefile

MiscUtilities.o:	MiscUtilities.c MiscUtilities.h TypeNLimit.h Makefile

QSufSort.o:	QSufSort.c QSufSort.h TypeNLimit.h MiscUtilities.h Makefile

r250.o:	r250.c r250.h Makefile

Socket.o:	Socket.c Socket.h TypeNLimit.h MemManager.h MiscUtilities.h Makefile

TextConverter.o:	TextConverter.c TextConverter.h TypeNLimit.h MemManager.h MiscUtilities.h r250.h Makefile

Timing.o:	Timing.c Timing.h Makefile

BWTBlast:	BWTBlast.o BWT.o MiscUtilities.o MemManager.o TextConverter.o r250.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o HSPstatistic.o karlin.o Timing.o blast_dust.o Socket.o Makefile
	$(CC) $(CFLAGS) BWTBlast.o BWT.o MiscUtilities.o MemManager.o TextConverter.o r250.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o HSPstatistic.o karlin.o Timing.o blast_dust.o Socket.o -o bwtblast -lm

BWTSW:	BWTSW.o BWT.o MiscUtilities.o MemManager.o TextConverter.o r250.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o HSPstatistic.o karlin.o blast_dust.o Timing.o Socket.o Makefile
	$(CC) $(CFLAGS) BWTSW.o BWT.o MiscUtilities.o MemManager.o TextConverter.o r250.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o HSPstatistic.o karlin.o blast_dust.o Timing.o Socket.o -o bwtsw -lm

BWTFormatdb:	BWTFormatdb.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o Timing.o Socket.o HSP.o HSPstatistic.o karlin.o Makefile
	$(CC) $(CFLAGS) BWTFormatdb.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o HSPstatistic.o karlin.o Timing.o Socket.o -o bwtformatdb -lm

BlastCompare:	BlastCompare.o MiscUtilities.o iniparser.o inistrlib.o dictionary.o Makefile
	$(CC) $(CFLAGS) BlastCompare.o MiscUtilities.o iniparser.o inistrlib.o dictionary.o -o blastcompare -lm

