#-----------------------------------------------------------------------#
#- GLOBAL DEFS ---------------------------------------------------------#
#-----------------------------------------------------------------------#

# Keep this as generic as possible.
NAME=riss
VERSION=7


SIG=$(NAME)_$(VERSION)
DIR=$(NAME)_$(VERSION)
TARGET=libipasir$(SIG).a
SHARED=libipasir$(SIG).so


#-----------------------------------------------------------------------#
#- REQUIRED TOP RULES --------------------------------------------------#
#-----------------------------------------------------------------------#

all:	$(TARGET)
shared: $(SHARED)

clean:
	rm -f  *.a *.so*
	rm -rf $(DIR)


#-----------------------------------------------------------------------#
#- INVISIBLE INTERNAL SUB RULES ----------------------------------------#
#-----------------------------------------------------------------------#

libipasir$(SIG).a: prepare
	cd $(DIR)/build; make riss-coprocessor-lib-static
	mv $(DIR)/build/lib/libriss-coprocessor.a .
	ln -sf libriss-coprocessor.a $@

libipasir$(SIG).so: prepare
	cd $(DIR)/build; make coprocessor-lib-shared riss-lib-shared
	mv $(DIR)/build/lib/lib*.so* .
	ln -sf libriss.so $@

prepare: .FORCE
	@#
	@# extract library
	@#
	rm -rf $(DIR)
	tar xvf $(DIR).tar.gz
	@#
	@# configure and build library
	@#
	mkdir -p $(DIR)/build
	cd $(DIR)/build; cmake -DQUIET=ON -DCMAKE_BUILD_TYPE=Release -DSTATIC_BINARIES=ON ..

.FORCE:
.PHONY: all shared prepare clean
