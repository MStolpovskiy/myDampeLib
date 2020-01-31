# Makefile for the myDampeLib
# Author : Mikhail Stolpovskiy

# ROOT specific instructions ------------------------------------------

RC     := root-config
ifeq ($(shell which $(RC) 2>&1 | sed -ne "s@.*/$(RC)@$(RC)@p"),$(RC))
MKARCH := $(wildcard $(shell $(RC) --etcdir)/Makefile.arch)
RCONFIG := $(wildcard $(shell $(RC) --incdir)/RConfigure.h)
endif
ifneq ($(MKARCH),)
include $(MKARCH)
else
ifeq ($(ROOTSYS),)
ROOTSYS = ..
endif
include $(ROOTSYS)/etc/Makefile.arch
endif

#------------------------------------------------------------------------------

CXXFLAGS += -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
CXXFLAGS += -I $(DMPSWSYS)/include
CXXFLAGS += -I include

LIBS += -L $(DMPSWSYS)/lib -lDmpEvent -lDmpService
LIBS += -lTMVA

#------------------------------------------------------------------------------

SrcSuf = cpp

CLUMONIO       = $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(wildcard *.$(SrcSuf)))
CLUMONIS       = $(wildcard *.$(SrcSuf))
CLUMONI        = vacorrvalidation$(ExeSuf)

MYLIBDS       = source
MYLIBDO       = obj
MYLIBS        = $(wildcard $(MYLIBDS)/*.$(SrcSuf))
MYLIBO        = $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(notdir $(MYLIBS)))
MYLIBOBJS     = $(addprefix $(MYLIBDO)/, $(MYLIBO))

OBJS          = $(MYLIBOBJS)
SRCS          = $(MYLIBS)

$(info OBJS is $(OBJS))
$(info SRCS is $(SRCS))

OUTPUTFILE    = bin/libmydampe.a

ifeq ($(ARCH),aix5)
MAKESHARED    = /usr/vacpp/bin/makeC++SharedLib
endif

#------------------------------------------------------------------------------


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

.PHONY: all

all:            $(OUTPUTFILE)

# $(OBJS):        $(MYLIBS)
# 				$(CXX) $(LDFLAGS) $(CXXFLAGS) -c $<

$(OUTPUTFILE):  $(OBJS)
				$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
				@rm -f $(OBJS) $(TRACKMATHSRC) $(OUTPUTFILE)


# .SUFFIXES: .$(SrcSuf)

# .$(SrcSuf).$(ObjSuf):
#         $(CXX)  $(CXXFLAGS) -c $<
