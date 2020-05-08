# 1. Put this file in the same folder as your 'driver' code 
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it)
#LIBRARY_DIR=$(HOME)/workspace/c_c++/ITensor/itensor
LIBRARY_DIR=external/itensor

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=mps_mylib

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=

# 5. For any additional .cc (source) files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cc

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

## google test ##

TEST_DIR=test

TEST_SRCS=$(TEST_DIR)/test.cc
TEST_OBJS=$(TEST_DIR)/test.o
TEST_TARGET=$(TEST_DIR)/test

GTEST_DIR=extsrc/googletest/googletest

GTEST_HEADERS=$(GTEST_DIR)/include/gtest/*.h\
				  $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS=$(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

$(TEST_DIR)/gtest-all.o : $(GTEST_SRCS)
	$(CCCOM) $(CCFLAGS) -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-o $@ $(GTEST_DIR)/src/gtest-all.cc

$(TEST_DIR)/gtest_main.o : $(GTEST_SRCS)
	$(CCCOM) $(CCFLAGS) -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-o $@ $(GTEST_DIR)/src/gtest_main.cc

$(TEST_DIR)/gtest.a : $(TEST_DIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

$(TEST_DIR)/gtest_main.a : $(TEST_DIR)/gtest-all.o $(TEST_DIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^


.PHONY: test
test: $(TEST_TARGET) 
	test/test

$(TEST_TARGET) : $(TEST_OBJS) $(ITENSOR_LIBS) $(TEST_DIR)/gtest_main.a
	$(CCCOM) $(CCFLAGS) $(TEST_OBJS) -o $(TEST_TARGET) $(TEST_DIR)/gtest_main.a $(LIBFLAGS) -lpthread

$(TEST_OBJS) : $(TEST_SRCS) $(GTEST_HEADERS) $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -I$(GTEST_DIR)/include -c $(CCFLAGS) -o $@ $<

##             ##

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g
	rm -fr test/*.o test/test

mkdebugdir:
	mkdir -p .debug_objs

