all: Disk_Math.h Polynomial.h test.cpp
	g++ -o test.exe Disk_Math.h Polynomial.h test.cpp -I .
clean:
	
