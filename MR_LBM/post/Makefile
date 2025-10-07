# Compiler and flags
CC_COMP = g++
FLAGS = -O3

# Files and output
SRC = post_stats.cpp
OUT = post

# Build and run
all:
	$(CC_COMP) $(FLAGS) $(SRC) -o $(OUT)
	./$(OUT)

# Clean rule
.PHONY: clean
clean:
	rm -f $(OUT) *.o *.mod