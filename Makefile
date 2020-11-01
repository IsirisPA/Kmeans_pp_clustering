OBJS	= main_lsh.o utilities.o lsh.o 
SOURCE	= main_lsh.cpp utilities.cpp lsh.cpp
HEADER	= utilities.hpp lsh.hpp
OUT		= lsh
OBJS1	= main_cube.o utilities.o cube.o 
SOURCE1	= main_cube.cpp utilities.cpp cube.cpp
HEADER1	= utilities.hpp cube.hpp
OUT1    = cube
OBJS2	= clustermain.o clusterkmeans.o lsh.o cube.o utilities.o 
SOURCE2 = clustermain.cpp clusterkmeans.cpp lsh.cpp cube.cpp utilities.cpp
HEADER2 = clusterkmeans.hpp lsh.hpp cube.hpp utilities.hpp
OUT2	= cluster
FLAGS 	= -g -c
CC 		= g++ -std=c++11



cluster: $(OBJS2)
	$(CC) -g   $(OBJS2) -o $(OUT2)

lsh: $(OBJS)
	$(CC) -g   $(OBJS) -o $(OUT)

cube: $(OBJS1)
	$(CC) -g   $(OBJS1) -o $(OUT1)

clustermain.o: clustermain.cpp
	$(CC)	$(FLAGS) clustermain.cpp

clusterkmeans.o: clusterkmeans.cpp
	$(CC)   $(FLAGS) clusterkmeans.cpp

main_lsh.o: main_lsh.cpp
	$(CC)	$(FLAGS) main_lsh.cpp

main_cube.o: main_cube.cpp
	$(CC)	$(FLAGS) main_cube.cpp

cube.o: cube.cpp
	$(CC)	$(FLAGS) cube.cpp

lsh.o: lsh.cpp
	$(CC)	$(FLAGS) lsh.cpp

utilities.o: utilities.cpp
	$(CC)	$(FLAGS) utilities.cpp 

clean:
	rm -f $(OBJS) $(OUT)
	rm -f $(OBJS1) $(OUT1)
	rm -f $(OBJS2) $(OUT2)

count:
	wc $(SOURCE) $(HEADER)