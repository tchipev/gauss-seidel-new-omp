CXXFLAGS =	-fopenmp -std=c++11 -O3 -march=native -Wall -fmessage-length=0
#CXXFLAGS =	-g -fopenmp -std=c++11 -O0 -march=native -Wall -fmessage-length=0

OBJS =		another-gauss-seidel.o GaussSeidel2D.o

LIBS = -fopenmp

TARGET =	another-gauss-seidel

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
