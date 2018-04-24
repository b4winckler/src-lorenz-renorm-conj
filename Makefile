.PHONY: clean all
LINK.o = $(LINK.cc)
CXXFLAGS += -O3 -std=c++11 -Wall -Ieigen-3.3.4
LDLIBS += -lmpfr -lgmp

FIXEDPT_OBJS = fixedpt.o
RENORM_OBJS = renorm.o
OBJS = $(FIXEDPT_OBJS) $(RENORM_OBJS)

all: fixedpt renorm

fixedpt: $(FIXEDPT_OBJS)

renorm: $(RENORM_OBJS)

clean:
	rm -f $(OBJS)
