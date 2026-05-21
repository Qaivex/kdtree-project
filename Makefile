
CC      ?= cc
STD     := -std=c11
WARN    := -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Wno-unused-parameter
OPT     := -O2
DBG     := -O0 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer
INC     := -Iinclude
LIBS    := -lm

SRC_DIR := src
BLD_DIR := build
DATA    := data

SRCS    := $(wildcard $(SRC_DIR)/*.c)
OBJS    := $(patsubst $(SRC_DIR)/%.c,$(BLD_DIR)/%.o,$(SRCS))
DOBJS   := $(patsubst $(SRC_DIR)/%.c,$(BLD_DIR)/%.dbg.o,$(SRCS))

BIN     := $(BLD_DIR)/kdtree
DBIN    := $(BLD_DIR)/kdtree_dbg

all: $(BIN)

$(BIN): $(OBJS) | $(BLD_DIR)
	$(CC) $(STD) $(WARN) $(OPT) $(INC) $^ -o $@ $(LIBS)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(STD) $(WARN) $(OPT) $(INC) -c $< -o $@

debug: $(DBIN)

$(DBIN): $(DOBJS) | $(BLD_DIR)
	$(CC) $(STD) $(WARN) $(DBG) $(INC) $^ -o $@ $(LIBS)

$(BLD_DIR)/%.dbg.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(STD) $(WARN) $(DBG) $(INC) -c $< -o $@

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

demo: $(BIN)
	@echo "=== KD nearest ==="
	./$(BIN) $(DATA)/lidar.csv -kd_nearest 1.0,2.0
	@echo
	@echo "=== Fuzzy C-means, c=3 ==="
	./$(BIN) $(DATA)/blobs.csv -cmeans 3
	@echo
	@echo "=== Elbow ==="
	./$(BIN) $(DATA)/blobs.csv -cmeans_elbow 8 --no-ansi
	@echo
	@echo "=== DBSCAN ==="
	./$(BIN) $(DATA)/moons.csv -dbscan 0.25,5

clean:
	rm -rf $(BLD_DIR)

.PHONY: all debug demo clean
