CXXFLAGS =	-std=c++11 -O3 -march=native -Wall -fmessage-length=0

OBJS =		another-gauss-seidel.o GaussSeidel2D.o

LIBS =

TARGET =	another-gauss-seidel

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
