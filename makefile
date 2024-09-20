INC = -I$(HOME)/local/include
LIBINC = -L$(HOME)/local/lib
#CC = g++ -Ofast
CC = g++ -std=c++17 -pthread
OBJhist = hist.o histograms.o vec3D.o string_pad.o files.o
OBJfil = main.o filament.o fitting.o interpolation.o files.o vec3D.o string_pad.o analyze.o
#TARGET = main

all: main hist

%.o : %.cpp
	$(CC) $< -c $(INC)

main: $(OBJfil)
	$(CC) -o main $(OBJfil) $(LIBINC) -lhdf5 -lhdf5_cpp -lcfitsio

hist: $(OBJhist)
	$(CC) -o hist $(OBJhist) $(LIBINC) -lhdf5 -lhdf5_cpp -lcfitsio

clean:
	rm -f $(OBJhist) $(OBJfil) main hist
