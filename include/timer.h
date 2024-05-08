#pragma once

#include <chrono>

// timer for host functions
#define TIMER(func)                                                                                                 \
  do {                                                                                                              \
    auto start = std::chrono::high_resolution_clock::now();                                                         \
    func;                                                                                                           \
    auto end = std::chrono::high_resolution_clock::now();                                                           \
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();                         \
    printf("FUNCTION CALL takes %d microseconds, '%s'in file %s line %d: \n", duration, #func, __FILE__, __LINE__); \
  } while (0)

#define TIMER_SPECIFIC_ITER(func, realIter, timerIter) \
  do {                                                 \
    if (realIter == timerIter) {                       \
      TIMER(func);                                     \
    } else {                                           \
      func;                                            \
    }                                                  \
  } while (0)
