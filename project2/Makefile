LDFLAGS = -larmadillo
CXXFLAGS = -O3 -std=c++11

new_prog: new_prog.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f new_prog
