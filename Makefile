src = $(wildcard *.C)
obj = $(src:.C=.o)

# include the library options determined by configure
# include $(LIBMESH_DIR)/Make.common
libmesh_CXX      := $(shell libmesh-config --cxx)
libmesh_CC       := $(shell libmesh-config --cc)
libmesh_F77      := $(shell libmesh-config --fc)
libmesh_F90      := $(shell libmesh-config --fc)
libmesh_INCLUDE  := $(shell libmesh-config --include)
libmesh_CPPFLAGS := $(shell libmesh-config --cppflags)
libmesh_CXXFLAGS := $(shell libmesh-config --cxxflags)
libmesh_CFLAGS   := $(shell libmesh-config --cflags)
# libmesh_FFLAGS   := $(shell libmesh-config --fflags)
libmesh_LIBS     := $(shell libmesh-config --libs)
# libmesh_shared   := $(shell $(libmesh_LIBTOOL) --config | grep build_libtool_libs | cut -d'=' -f2)
# libmesh_static   := $(shell $(libmesh_LIBTOOL) --config | grep build_old_libs | cut -d'=' -f2)


CXXFLAGS=${libmesh_CXXFLAGS}
# CXXFLAGS+=${libmesh_INCLUDE}
CXXFLAGS+=-I$(OPENMPI_4_INC)
CXXFLAGS+=-I$(LIBMESH_INC)
CXXFLAGS+=-I$(PETSC3p13_INC)
CXXFLAGS+=-I.
LDFLAGS=$(libmesh_LIBS)
TARGET=fc
$(TARGET): $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
# 	echo ciao
# 
%.o: %.c
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# ex1:$(obj)
# 	@echo "Linking "$@"..."
# 	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link \
# 	$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS)

.PHONY: clean
clean:
	rm -f $(obj) $(TARGET)
