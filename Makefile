CFLAGS = -Wall -O3 -std=c99
LIBS = align.o alnfrags.o ankers.o assembly.o chain.o compdna.o compkmers.o compress.o decon.o ef.o filebuff.o frags.o hashmap.o hashmapindex.o hashmapkma.o hashmapkmers.o hashtable.o index.o kma.o kmapipe.o kmers.o loadupdate.o makeindex.o mt1.o nw.o pherror.o printconsensus.o qseqs.o qualcheck.o runinput.o runkma.o savekmers.o seq2fasta.o seqparse.o shm.o sparse.o spltdb.o stdnuc.o stdstat.o update.o updateindex.o updatescores.o valueshash.o vcf.o
PROGS = kma kma_index kma_shm kma_update

.c .o:
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(PROGS)

kma: main.c libkma.a
	$(CC) $(CFLAGS) -o $@ main.c libkma.a -lm -lpthread -lz

kma_index: kma_index.c libkma.a
	$(CC) $(CFLAGS) -o $@ kma_index.c libkma.a -lm -lz

kma_shm: kma_shm.c libkma.a
	$(CC) $(CFLAGS) -o $@ kma_shm.c libkma.a

kma_update: kma_update.c libkma.a
	$(CC) $(CFLAGS) -o $@ kma_update.c libkma.a

libkma.a: $(LIBS)
	$(AR) -csru $@ $(LIBS)

clean:
	$(RM) $(LIBS) $(PROGS) libkma.a

align.o: align.h chain.h compdna.h hashmapindex.h nw.h stdnuc.h stdstat.h
alnfrags.o: alnfrags.h align.h ankers.h compdna.h hashmapindex.h qseqs.h threader.h updatescores.h
ankers.o: ankers.h compdna.h pherror.h qseqs.h
assembly.o: assembly.h align.h filebuff.h kmapipe.h pherror.h stdnuc.h stdstat.h threader.h
chain.o: chain.h penalties.h pherror.h stdstat.h
compdna.o: compdna.h pherror.h stdnuc.h
compkmers.o: compkmers.h pherror.h
compress.o: compress.h hashmap.h hashmapkma.h pherror.h valueshash.h
decon.o: decon.h compdna.h filebuff.h hashmapkma.h seqparse.h stdnuc.h qseqs.h updateindex.h
ef.o: ef.h assembly.h stdnuc.h vcf.h version.h
filebuff.o: filebuff.h pherror.h qseqs.h
frags.o: frags.h filebuff.h pherror.h qseqs.h
hashmap.o: hashmap.h hashtable.h pherror.h
hashmapindex.o: hashmapindex.h pherror.h stdnuc.h
hashmapkma.o: hashmapkma.h pherror.h
hashmapkmers.o: hashmapkmers.h pherror.h
hashtable.o: hashtable.h hashmapkma.h hashmapkmers.h pherror.h
index.o: index.h compress.h decon.h hashmap.h hashmapkma.h loadupdate.h makeindex.h pherror.h stdstat.h version.h
kma.o: kma.h ankers.h assembly.h chain.h hashmapkma.h kmers.h mt1.h penalties.h pherror.h qseqs.h runinput.h runkma.h savekmers.h sparse.h spltdb.h version.h
kmapipe.o: kmapipe.h pherror.h
kmers.o: kmers.h ankers.h compdna.h hashmapkma.h kmapipe.h pherror.h qseqs.h savekmers.h spltdb.h
loadupdate.o: loadupdate.h pherror.h hashmap.h hashmapkma.h updateindex.h
makeindex.o: makeindex.h compdna.h filebuff.h hashmap.h pherror.h qseqs.h seqparse.h updateindex.h
mt1.o: mt1.h assembly.h chain.h filebuff.h hashmapindex.h kmapipe.h nw.h penalties.h pherror.h printconsensus.h qseqs.h runkma.h stdstat.h vcf.h
nw.o: nw.h pherror.h stdnuc.h penalties.h
pherror.o: pherror.h
printconsensus.o: printconsensus.h assembly.h
qseqs.o: qseqs.h pherror.h
qualcheck.o: qualcheck.h compdna.h hashmap.h pherror.h stdnuc.h stdstat.h
runinput.o: runinput.h compdna.h filebuff.h pherror.h qseqs.h seqparse.h
runkma.o: runkma.h align.h alnfrags.h assembly.h chain.h compdna.h ef.h filebuff.h frags.h hashmapindex.h kmapipe.h nw.h pherror.h printconsensus.h qseqs.h stdnuc.h stdstat.h vcf.h
savekmers.o: savekmers.h ankers.h compdna.h hashmapkma.h penalties.h pherror.h qseqs.h stdnuc.h stdstat.h threader.h
seq2fasta.o: seq2fasta.h pherror.h qseqs.h runkma.h stdnuc.h
seqparse.o: seqparse.h filebuff.h qseqs.h
shm.o: shm.h pherror.h hashmapkma.h version.h
sparse.o: sparse.h compkmers.h hashtable.h kmapipe.h pherror.h runinput.h savekmers.h stdnuc.h stdstat.h
spltdb.o: spltdb.h align.h alnfrags.h assembly.h chain.h compdna.h ef.h filebuff.h frags.h hashmapindex.h kmapipe.h nw.h pherror.h printconsensus.h qseqs.h runkma.h stdnuc.h stdstat.h vcf.h
stdnuc.o: stdnuc.h
stdstat.o: stdstat.h
update.o: update.h hashmapkma.h pherror.h stdnuc.h
updateindex.o: updateindex.h compdna.h hashmap.h hashmapindex.h pherror.h qualcheck.h stdnuc.h pherror.h
updatescores.o: updatescores.h qseqs.h
valueshash.o: valueshash.h pherror.h
vcf.o: vcf.h assembly.h filebuff.h stdnuc.h stdstat.h version.h
