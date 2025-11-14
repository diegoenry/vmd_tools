# Makefile for testing GROMACS topology plugin

# VMD plugin include path
VMDPLUGIN_INC = /Users/deb0054/gitlab/vmd2prototype/plugins/include

# Compiler settings
CC = gcc
CXX = g++
CFLAGS = -Wall -g -I$(VMDPLUGIN_INC)
CXXFLAGS = -Wall -g -I$(VMDPLUGIN_INC)

# Target executables
TARGET = test_grotop
TARGET2 = test_grotop_to_psf
TARGET3 = test_grotop_to_js

# Source files
SRCS = test_grotop.c
SRCS2 = test_grotop_to_psf.c
SRCS3 = test_grotop_to_js.c
OBJS = $(SRCS:.c=.o)
OBJS2 = $(SRCS2:.c=.o)
OBJS3 = $(SRCS3:.c=.o)
GROMACS_WRAPPER_OBJ = gromacs_wrapper.o

# Default target
all: $(TARGET) $(TARGET2) $(TARGET3)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

$(TARGET2): $(OBJS2)
	$(CC) $(CFLAGS) -o $(TARGET2) $(OBJS2)

$(TARGET3): $(OBJS3) $(GROMACS_WRAPPER_OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(OBJS3) $(GROMACS_WRAPPER_OBJ)

%.o: %.c grotopplugin.c
	$(CC) $(CFLAGS) -c $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Test targets
test: $(TARGET)
	@echo "=== Testing with example_topol.top ==="
	./$(TARGET) example_topol.top
	@echo ""
	@echo "=== Testing with insane.top ==="
	./$(TARGET) insane.top
	@echo ""
	@echo "=== Testing with BIGtopol.top ==="
	./$(TARGET) BIGtopol.top

# Individual test targets
test-example: $(TARGET)
	./$(TARGET) example_topol.top

test-insane: $(TARGET)
	./$(TARGET) insane.top

test-big: $(TARGET)
	./$(TARGET) BIGtopol.top

# Test PSF conversion
test-psf: $(TARGET2)
	@echo "=== Converting example_topol.top to PSF ==="
	./$(TARGET2) example_topol.top example_output.psf
	@echo ""
	@echo "=== Comparing with Python script output ==="
	python3 gromacs_to_psf.py example_topol.top python_output.psf
	@echo ""
	@echo "=== Checking file sizes ==="
	ls -lh example_output.psf python_output.psf

# Test JS conversion
test-js: $(TARGET3)
	@echo "=== Converting topol.top + bilayer.gro to JS ==="
	./$(TARGET3) topol.top bilayer.gro output.js
	@echo ""
	@echo "=== Checking output file ==="
	ls -lh output.js
	@echo ""
	@echo "=== File header (first 100 bytes) ==="
	head -c 100 output.js | od -c | head -10

# Clean
clean:
	rm -f $(TARGET) $(TARGET2) $(TARGET3) $(OBJS) $(OBJS2) $(OBJS3) $(GROMACS_WRAPPER_OBJ) *.psf *.js

.PHONY: all test test-example test-insane test-big test-psf test-js clean
