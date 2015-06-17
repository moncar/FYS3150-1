LDFLAGS = -larmadillo
CXXFLAGS = -g -std=c++11

main: randomNumbers.cpp walker.cpp methods.cpp main.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm main

remall:
	rm average_1D.txt test_1D.txt average_2D.txt test_2D.txt mc2D.txt expEuler2D.txt impJacobi2D.txt
