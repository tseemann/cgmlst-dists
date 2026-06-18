
CC = gcc
CFLAGS = -Wall -Wextra -Ofast -fopenmp
LIBSDIR =
LIBS = -lm

EXE = cgmlst-dists
PREFIX = $(HOME)/.local
TESTDIR = test

.PHONY: all check clean
.DEFAULT: all

all: $(EXE)

$(EXE): main.c
	$(CC) $(CFLAGS) -o $(EXE) $< $(LIBSDIR) $(LIBS)

install: $(EXE)
	install -v -t $(PREFIX)/bin $(EXE)

clean:
	$(RM) *~ *.o $(EXE)

check: $(EXE)
	./$(EXE) -v
	./$(EXE) /dev/null || true
	./$(EXE) -q $(TESTDIR)/chewie.tab
	./$(EXE) -c $(TESTDIR)/boring.tab
	./$(EXE) -m 1 $(TESTDIR)/100.tab > /dev/null
	! ./$(EXE) -q $(TESTDIR)/ragged.tab > /dev/null 2>&1

	
