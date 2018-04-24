.PHONY: clean all
LINK.o = $(LINK.cc)
CXXFLAGS += -O3 -std=c++11 -Wall -Ieigen-3.3.4
LDLIBS += -lmpfr -lgmp

FIXEDPT_OBJS = fixedpt.o
OBJS = $(FIXEDPT_OBJS)

all: fixedpt

fixedpt: $(FIXEDPT_OBJS)

clean:
	rm -f $(OBJS)
