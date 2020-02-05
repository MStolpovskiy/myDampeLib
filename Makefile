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

CXXFLAGS += -std=c++14 -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
CXXFLAGS += -I $(DMPSWSYS)/include
CXXFLAGS += -I include

LIBS += -L $(DMPSWSYS)/lib -lDmpEvent -lDmpService
LIBS += -lTMVA

#------------------------------------------------------------------------------

SrcSuf := cpp
ObjSuf := o
IncSuf := hpp

MYLIBDS := source
MYLIBDO := obj
MYLIBDI := include
MYLIBDB := bin

MYLIBS  := $(wildcard $(MYLIBDS)/*.$(SrcSuf))
MYLIBO_ := $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(notdir $(MYLIBS)))
MYLIBO  := $(addprefix $(MYLIBDO)/, $(MYLIBO_))
MYLIBI_ := $(patsubst %.$(SrcSuf),%.$(IncSuf),$(notdir $(MYLIBS)))
MYLIBI  := $(addprefix $(MYLIBDI)/, $(MYLIBI_))

OBJS    := $(MYLIBO)
SRCS    := $(MYLIBS)
INCS    := $(MYLIBI)

$(info OBJS is $(OBJS))
$(info SRCS is $(SRCS))
$(info INCS is $(INCS))

OUTPUTFILE    = $(MYLIBDB)/libmydampe.a

ifeq ($(ARCH),aix5)
MAKESHARED    = /usr/vacpp/bin/makeC++SharedLib
endif

#------------------------------------------------------------------------------

.PHONY: all clean

all:            $(OUTPUTFILE)

# $(OBJS):        $(SRCS) $(INCS)
# 				$(CXX) $(CXXFLAGS) -c $(SRCS)

$(MYLIBDO)/%.$(ObjSuf): $(MYLIBDS)/%.$(SrcSuf) | $(MYLIBDO)
				$(CXX) $(CXXFLAGS) -c $< -o $@
				@echo "Compilation complete!"

$(MYLIBDO):
				mkdir $@

$(OUTPUTFILE):  $(OBJS) | $(MYLIBDB)
				ar rcs $@ $^
				@echo "Linking complete!"
# 				$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

$(MYLIBDB):
				mkdir $@

clean:
				@rm -f $(OBJS) $(TRACKMATHSRC) $(OUTPUTFILE)
