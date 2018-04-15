.PHONY: clean all
LINK.o = $(LINK.cc)
CXXFLAGS += -O3 -std=c++11 -Wall -Ieigen-3.3.4
LDLIBS += -lmpfr -lgmp

FIXED_OBJS = fixed.o
LORENZ_OBJS = lorenz.o
RENORM_OBJS = renorm.o
TEST_OBJS = test.o
THURSTON_OBJS = thurston.o
OBJS = $(FIXED_OBJS) $(LORENZ_OBJS) $(RENORM_OBJS) $(TEST_OBJS) $(THURSTON_OBJS)

all: fixed lorenz renorm test thurston

fixed: $(FIXED_OBJS)

lorenz: $(LORENZ_OBJS)

renorm: $(RENORM_OBJS)

test: $(TEST_OBJS)

thurston: $(THURSTON_OBJS)

clean:
	rm -f $(OBJS)
