.PHONY: clean all
LINK.o = $(LINK.cc)
CXXFLAGS += -O3 -std=c++11 -Wall -Impfrc++-3.6.2
LDLIBS += -lmpfr -lgmp

LORENZ_OBJS = lorenz.o
OBJS = $(LORENZ_OBJS)

all: lorenz

lorenz: $(LORENZ_OBJS)

clean:
	rm -f $(OBJS)
