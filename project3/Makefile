LDFLAGS = -larmadillo
CXXFLAGS = -g -O3 -std=c++11

main: body.cpp solarsystem.cpp methods.cpp main.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm main
