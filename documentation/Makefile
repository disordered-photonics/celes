PROJECT_BASE_DIR = $(shell dirname `pwd`)

.PHONY: all documentation clean push

all: documentation

documentation:
	PROJECT_BASE_DIR=$(PROJECT_BASE_DIR) doxygen Doxyfile_linux

clean:
	rm -rf html/

push: documentation
	npm run deploy
