PROJECT_BASE_DIR = $(shell dirname `pwd`)

.PHONY: all documentation clean

all: documentation

documentation:
	PROJECT_BASE_DIR=$(PROJECT_BASE_DIR) doxygen Doxyfile_linux

clean:
	rm -rf html/
