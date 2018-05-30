CFLAGS = -std=c99 -w -O3
BINS = kma kma_index kma_shm

all: $(BINS)

kma: KMA.c
	$(CC) $(CFLAGS) -o $@ $< -lm -lpthread

kma_index: KMA_index.c
	$(CC) $(CFLAGS) -o $@ $< -lm

kma_shm: KMA_SHM.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(BINS)
