INC = -I$(HOME)/local/include
LIBINC = -L$(HOME)/local/lib
#CC = g++ -Ofast
CC = g++ -std=c++17 -pthread
OBJFILES = main.o filament.o fitting.o interpolation.o files.o vec3D.o
TARGET = main

all: $(TARGET)

%.o : %.cpp
	$(CC) $< -c $(INC)

$(TARGET): $(OBJFILES)
	$(CC) -o $(TARGET) $(OBJFILES) $(LIBINC) -lhdf5 -lhdf5_cpp

clean:
	rm -f $(OBJFILES) $(TARGET)
