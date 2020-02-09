CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		another-gauss-seidel.o GaussSeidel2D.o

LIBS =

TARGET =	another-gauss-seidel

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
