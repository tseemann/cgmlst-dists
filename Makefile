
CC = gcc
CFLAGS = -Wall -Wextra -Ofast
LIBSDIR =
LIBS = -lm

EXE = cgmlst-dists
PREFIX = /usr/local
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
	
	
