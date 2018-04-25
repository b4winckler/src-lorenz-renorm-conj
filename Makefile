.PHONY: clean all
LINK.o = $(LINK.cc)
CXXFLAGS += -O3 -std=c++11 -Wall -Ieigen-3.3.4
LDLIBS += -lmpfr -lgmp

DERIV_OBJS = deriv.o
FIXEDPT_OBJS = fixedpt.o
RENORM_OBJS = renorm.o
OBJS = $(DERIV_OBJS) $(FIXEDPT_OBJS) $(RENORM_OBJS)

all: deriv fixedpt renorm

deriv: $(DERIV_OBJS)

fixedpt: $(FIXEDPT_OBJS)

renorm: $(RENORM_OBJS)

clean:
	rm -f $(OBJS)

%.o: %.cpp lorenz.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<
