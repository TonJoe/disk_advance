all: Disk_Math.h Polynomial.h test.cpp
	g++ Disk_Math.h Polynomial.h test.cpp -I .
clean:
	