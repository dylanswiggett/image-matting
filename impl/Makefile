main: main.o GridSolver.o GridSampler.o ImageManager.o FGBGMatte.o
	g++ -Wall -g -std=c++11 -o main main.o GridSolver.o GridSampler.o ImageManager.o FGBGMatte.o -lSDL

main.o: main.cpp GridSolver.hpp GridSampler.hpp ImageManager.hpp FGBGMatte.hpp
	g++ -Wall -g -std=c++11 -c -o main.o main.cpp

GridSolver.o: GridSolver.cpp GridSolver.hpp GridSampler.hpp
	g++ -Wall -g -std=c++11 -c -o GridSolver.o GridSolver.cpp

GridSampler.o: GridSampler.cpp GridSampler.hpp
	g++ -Wall -g -std=c++11 -c -o GridSampler.o GridSampler.cpp

ImageManager.o: ImageManager.cpp ImageManager.hpp
	g++ -Wall -g -std=c++11 -c -o ImageManager.o ImageManager.cpp

FGBGMatte.o: FGBGMatte.cpp FGBGMatte.hpp
	g++ -Wall -g -std=c++11 -c -o FGBGMatte.o FGBGMatte.cpp