CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic -O2 -Iinclude
DEBUG_FLAGS := -g -O0 -fsanitize=address,undefined
TEST_FLAGS := -std=c++20 -Wall -Wextra -Wpedantic -Iinclude

TEST_DIR := tests
TEST_BUILD := build/tests
EXAMPLE_DIR := examples
EXAMPLE_BUILD := build/examples

TEST_SOURCES := $(shell find $(TEST_DIR) -name '*.cpp' 2>/dev/null)
TEST_BINS := $(patsubst $(TEST_DIR)/%.cpp,$(TEST_BUILD)/%,$(TEST_SOURCES))

EXAMPLE_SOURCES := $(wildcard $(EXAMPLE_DIR)/*.cpp)
EXAMPLE_BINS := $(patsubst $(EXAMPLE_DIR)/%.cpp,$(EXAMPLE_BUILD)/%,$(EXAMPLE_SOURCES))

.PHONY: all clean test examples debug help

all: test

$(TEST_BUILD)/%: $(TEST_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(TEST_FLAGS) $< -o $@

$(EXAMPLE_BUILD)/%: $(EXAMPLE_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $< -o $@

define run_tests
	@echo "Running tests"
	@failed=0; \
	for test in $(1); do \
		echo "$$test"; \
		if $$test; then \
			echo "  PASS"; \
		else \
			echo "  FAIL"; \
			failed=$$((failed + 1)); \
		fi; \
	done; \
	if [ $$failed -eq 0 ]; then \
		echo "All tests passed"; \
	else \
		echo "$$failed test(s) failed"; \
		exit 1; \
	fi
endef

test: $(TEST_BINS)
	$(call run_tests,$(TEST_BINS))

test-linalg: $(filter $(TEST_BUILD)/test_linalg%,$(TEST_BINS))
	$(call run_tests,$(filter $(TEST_BUILD)/test_linalg%,$(TEST_BINS)))

debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: TEST_FLAGS += $(DEBUG_FLAGS)
debug: test

examples: $(EXAMPLE_BINS)

clean:
	rm -rf build

run-%: $(EXAMPLE_BUILD)/%
	./$

help:
	@echo "Available targets:"
	@echo "  make test          - Run all tests"
	@echo "  make test-[topic]  - Run all tests within selected topic/concept yk selective testing"
	@echo "  make debug         - Run tests with sanitizers"
	@echo "  make examples      - Build all examples"
	@echo "  make clean         - Remove build artifacts"
