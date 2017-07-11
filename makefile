all: Disk_Math.h Polynomial.h test.cpp Monca.h
	g++ Disk_Math.h Polynomial.h Monca.h test.cpp -I.
clean:
	
