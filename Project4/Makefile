LDFLAGS = -larmadillo
CXXFLAGS = -g -O3 -std=c++11

main: triDiagonalSolver.cpp methods.cpp main.cpp
	g++ $(CXXFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm main

remall:
	rm t1.txt t2.txt t3.txt t4.txt t5.txt

allEE:
	./main 0.002 EE >> t1.txt
	./main 0.05 EE >> t2.txt
	./main 0.1 EE >> t3.txt
	./main 0.2 EE >> t4.txt
	./main 0.5 EE >> t5.txt
	
allIE:
	./main 0.002 IE >> t1.txt
	./main 0.05 IE >> t2.txt
	./main 0.1 IE >> t3.txt
	./main 0.2 IE >> t4.txt
	./main 0.5 IE >> t5.txt

allCN:
	./main 0.002 CN >> t1.txt
	./main 0.05 CN >> t2.txt
	./main 0.1 CN >> t3.txt
	./main 0.2 CN >> t4.txt
	./main 0.5 CN >> t5.txt
