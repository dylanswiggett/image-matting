FLAGS = -O2

main: main.o GridSampler.o ImageManager.o FGBGMatte.o
	g++ -Wall -g -std=c++11 -o main main.o GridSampler.o ImageManager.o FGBGMatte.o -lSDL

main.o: main.cpp GridSampler.hpp ImageManager.hpp FGBGMatte.hpp 
	g++ -Wall -g -std=c++11 -c -o main.o main.cpp $(FLAGS)

GridSampler.o: GridSampler.cpp GridSampler.hpp
	g++ -Wall -g -std=c++11 -c -o GridSampler.o GridSampler.cpp $(FLAGS)

ImageManager.o: ImageManager.cpp ImageManager.hpp
	g++ -Wall -g -std=c++11 -c -o ImageManager.o ImageManager.cpp $(FLAGS)

FGBGMatte.o: FGBGMatte.cpp FGBGMatte.hpp
	g++ -Wall -g -std=c++11 -c -o FGBGMatte.o FGBGMatte.cpp $(FLAGS)

clean:
	rm -f *.o