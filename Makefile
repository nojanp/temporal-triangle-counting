INCDIR =./ETTC
CC = g++
CFLAGS=-I$(INCDIR) -Wall -std=c++11 -g -O3

ODIR=./
# LDIR =../lib

_DEPS = Graph.h t_triangle_counting.h
DEPS = $(patsubst %,$(INCDIR)/%,$(_DEPS))

_OBJ = Graph.o t_triangle_counting.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ettc: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

