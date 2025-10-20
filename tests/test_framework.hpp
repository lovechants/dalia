#ifndef MATH_TEST_FRAMEWORK_HPP
#define MATH_TEST_FRAMEWORK_HPP

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <source_location>

namespace math::test {

struct TestFailure {
    std::string message;
    std::source_location location;
};

struct TestRegistry {
    std::vector<std::pair<std::string, void(*)()>> tests;
    
    static TestRegistry& instance() {
        static TestRegistry registry;
        return registry;
    }
    
    void add(std::string name, void(*fn)()) {
        tests.push_back({std::move(name), fn});
    }
    
    int run() {
        int passed = 0;
        int failed = 0;
        
        for (const auto& [name, fn] : tests) {
            try {
                fn();
                std::cout << "  PASS: " << name << "\n";
                passed++;
            } catch (const TestFailure& f) {
                std::cout << "  FAIL: " << name << "\n";
                std::cout << "    " << f.message << " at " 
                          << f.location.file_name() << ":" << f.location.line() << "\n";
                failed++;
            }
        }
        
        std::cout << passed << " passed, " << failed << " failed\n";
        return failed > 0 ? 1 : 0;
    }
};

struct TestAdder {
    TestAdder(const char* name, void(*fn)()) {
        TestRegistry::instance().add(name, fn);
    }
};

template<typename T>
void assert_eq(const T& expected, const T& actual, 
               std::source_location loc = std::source_location::current()) {
    if (!(expected == actual)) {
        throw TestFailure{"Values not equal", loc};
    }
}

template<typename T>
void assert_near(const T& expected, const T& actual, const T& epsilon,
                 std::source_location loc = std::source_location::current()) {
    if (std::abs(expected - actual) > epsilon) {
        throw TestFailure{"Values not within epsilon", loc};
    }
}

void assert_true(bool condition, 
                 std::source_location loc = std::source_location::current()) {
    if (!condition) {
        throw TestFailure{"Condition was false", loc};
    }
}

}

#define TEST(name) \
    void test_##name(); \
    static math::test::TestAdder adder_##name(#name, test_##name); \
    void test_##name()

#define RUN_ALL_TESTS() \
    int main() { \
        return math::test::TestRegistry::instance().run(); \
    }

#endif
