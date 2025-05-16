INC = -I$(HOME)/local/include
LIBINC = -L$(HOME)/local/lib
#CC = g++ -Ofast
CC = g++ -std=c++17 -pthread
OBJhist = histograms.o hist.o vec3D.o string_pad.o files.o
OBJfil = main.o filaments_main_functions.o filament.o fitting.o interpolation.o files.o vec3D.o string_pad.o analyze.o argparser.o options.o
#TARGET = main

all: main hist

%.o : %.cpp
	$(CC) $< -c $(INC)

# -lhdf5 -lhdf5_cpp
main: $(OBJfil)
	$(CC) -o flmnt $(OBJfil) $(LIBINC) -lhdf5 -lhdf5_cpp -lcfitsio

hist: $(OBJhist)
	$(CC) -o hist $(OBJhist) $(LIBINC) -lhdf5 -lhdf5_cpp -lcfitsio

clean:
	rm -f $(OBJhist) $(OBJfil) flmnt hist