.PHONY: clean all
LINK.o = $(LINK.cc)
CXXFLAGS += -O3 -std=c++11 -Wall -Ieigen-3.3.4
LDLIBS += -lmpfr -lgmp

LORENZ_OBJS = lorenz.o
TEST_OBJS = test.o
THURSTON_OBJS = thurston.o
OBJS = $(LORENZ_OBJS) $(TEST_OBJS) $(THURSTON_OBJS)

all: lorenz test thurston

lorenz: $(LORENZ_OBJS)

test: $(TEST_OBJS)

thurston: $(THURSTON_OBJS)

clean:
	rm -f $(OBJS)
