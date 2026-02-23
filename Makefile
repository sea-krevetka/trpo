CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -g
LDFLAGS = -lm

OPENBLAS_DIR = ./OpenBLAS
OPENBLAS_INC = -I$(OPENBLAS_DIR)
OPENBLAS_LIB = -L$(OPENBLAS_DIR) -lopenblas

GREEN = \033[0;32m
RED = \033[0;31m
NC = \033[0m

all: test_openblas test_fake

test_openblas: test_cblas_level3.c
	@echo "$(GREEN)Building tests with OpenBLAS...$(NC)"
	$(CC) $(CFLAGS) -DUSE_OPENBLAS $(OPENBLAS_INC) \
		-o $@ $< $(OPENBLAS_LIB) $(LDFLAGS)
	@echo "$(GREEN)Done: test_openblas$(NC)"

fake_cblas.o: fake_cblas.c cblas.h
	$(CC) $(CFLAGS) -fPIC -c -o $@ fake_cblas.c

libfakeblas.a: fake_cblas.o
	ar rcs $@ $^

test_fake: test_cblas_level3.c libfakeblas.a
	@echo "$(RED)Building tests with FAKE library...$(NC)"
	$(CC) $(CFLAGS) -I. -o $@ test_cblas_level3.c -L. -lfakeblas $(LDFLAGS)
	@echo "$(RED)Done: test_fake$(NC)"

run_openblas: test_openblas
	@echo "\n$(GREEN)=== Running tests with OpenBLAS ===$(NC)\n"
	LD_LIBRARY_PATH=$(OPENBLAS_DIR):$$LD_LIBRARY_PATH ./test_openblas

run_fake: test_fake
	@echo "\n$(RED)=== Running tests with FAKE library ===$(NC)\n"
	LD_LIBRARY_PATH=.$$LD_LIBRARY_PATH ./test_fake

run_all: run_openblas run_fake

clean:
	rm -f *.o *.a *.so test_openblas test_fake
	@echo "Cleaned"

check_openblas:
	@if [ ! -f $(OPENBLAS_DIR)/libopenblas.a ] && [ ! -f $(OPENBLAS_DIR)/libopenblas.so ]; then \
		echo "$(RED)Error: OpenBLAS not found in $(OPENBLAS_DIR)$(NC)"; \
		echo "Please build OpenBLAS first:"; \
		echo "  cd OpenBLAS"; \
		echo "  make"; \
		exit 1; \
	else \
		echo "$(GREEN)OpenBLAS found in $(OPENBLAS_DIR)$(NC)"; \
	fi

.PHONY: all run_openblas run_fake run_all clean check_openblas