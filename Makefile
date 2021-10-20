CC=gcc
SRC=cli.c decoder.c decoder_bp.c error_floor.c qcmdpc_decoder.c sparse_cyclic.c threshold.c weak.c xoroshiro128plus.c
OBJ=$(SRC:%.c=%.o)
DEP=$(SRC:%.c=%.d)
PROF=$(SRC:%.c=%.gcda)
LFLAGS=-lm
CFLAGS=-Wall -std=gnu11 $(OPT) $(EXTRA)
ifdef AVX
    CFLAGS+=-DAVX
endif
ifdef PROFGEN
    CFLAGS+=-fprofile-generate
endif
ifdef PROFUSE
    CFLAGS+=-fprofile-use -fprofile-correction
endif

default: avx2

noavx:
	make OPT="-Ofast -march=native -flto" qcmdpc_decoder

avx2:
	make OPT="-Ofast -march=native -flto" AVX=1 qcmdpc_decoder_avx2

format:
	clang-format -i -style=file $(SRC) *.h

qcmdpc_decoder: $(OBJ)
	$(CC) $(CFLAGS) -pthread $^ -o $@ $(LFLAGS)

qcmdpc_decoder_avx2: $(OBJ)
	$(CC) $(CFLAGS) -pthread $^ -o $@ $(LFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -MMD -c -o $@ $<

-include $(DEP)

clean:
	- /bin/rm qcmdpc_decoder qcmdpc_decoder_avx2 $(OBJ) $(DEP) $(PROF)
