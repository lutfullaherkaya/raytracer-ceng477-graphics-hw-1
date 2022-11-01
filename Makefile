OBJS	= parser.o ppm.o raytracer.o tinyxml2.o
SOURCE	= parser.cpp ppm.cpp raytracer.cpp tinyxml2.cpp
HEADER	= parser.h ppm.h tinyxml2.h bvh.h
OUT	= raytracer
CC	 = g++
FLAGS	 = -g -c -Wall  -O3 # todo: o3 ekle
LFLAGS	 = 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

parser.o: parser.cpp
	$(CC) $(FLAGS) parser.cpp -std=c++11

ppm.o: ppm.cpp
	$(CC) $(FLAGS) ppm.cpp -std=c++11

raytracer.o: raytracer.cpp
	$(CC) $(FLAGS) raytracer.cpp -std=c++11

tinyxml2.o: tinyxml2.cpp
	$(CC) $(FLAGS) tinyxml2.cpp -std=c++11


clean:
	rm -f $(OBJS) $(OUT)

run: $(OUT)
	./$(OUT)
