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
ObjSuf = o
IncSuf = hpp

MYLIBDS       = source
MYLIBDO       = obj
MYLIBDI       = include

MYLIBS        = $(wildcard $(MYLIBDS)/*.$(SrcSuf))
MYLIBO_       = $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(notdir $(MYLIBS)))
MYLIBO        = $(addprefix $(MYLIBDO)/, $(MYLIBO_))
MYLIBI_       = $(patsubst %.$(SrcSuf),%.$(IncSuf),$(notdir $(MYLIBS)))
MYLIBI        = $(addprefix $(MYLIBDI)/, $(MYLIBI_))

OBJS          = $(MYLIBO)
SRCS          = $(MYLIBS)
INCS          = $(MYLIBI)

$(info OBJS is $(OBJS))
$(info SRCS is $(SRCS))
$(info INCS is $(INCS))

OUTPUTFILE    = bin/libmydampe.a

ifeq ($(ARCH),aix5)
MAKESHARED    = /usr/vacpp/bin/makeC++SharedLib
endif

#------------------------------------------------------------------------------

.PHONY: all

all:            $(OUTPUTFILE)

$(OBJS):        $(SRCS) $(INCS)
				$(CXX) $(CXXFLAGS) -c $(SRCS)

$(OUTPUTFILE):  $(OBJS)
				$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
				@rm -f $(OBJS) $(TRACKMATHSRC) $(OUTPUTFILE)
