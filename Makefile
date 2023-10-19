CXX = g++
CXXFLAGS += -std=c++17 -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG

INCLUDES =
LIBS =

# atlas / blas
CXXFLAGS += -DUSE_BLAS
LIBS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# likwid
CXXFLAGS += -DUSE_LIKWID -pthread
INCLUDES += -I/apps/likwid/5.2.1/include
LDFLAGS += -L/apps/likwid/5.2.1/lib
LIBS += -llikwid

TARGET = ex01_example
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

$(TARGET).o: $(TARGET).cpp Timer.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp

clean:
	@$(RM) -rf *.o $(TARGET)
