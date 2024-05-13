#pragma once

#include <chrono>

#define OUTPUT_TIME

#ifdef OUTPUT_TIME
// timer for host functions

// TIMER 对某个程序段/函数调用进行计时
// TIMER_SPECIFIC_ITER 对某个程序段/函数调用进行计时，但是只在特定的迭代次数时进行计时(如CG算法只在第10轮进行计算)
#define TIMER(func)                                                                                                 \
  do {                                                                                                              \
    auto start = std::chrono::high_resolution_clock::now();                                                         \
    func;                                                                                                           \
    auto end = std::chrono::high_resolution_clock::now();                                                           \
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();                         \
    printf("FUNCTION CALL takes %e seconds, '%s'in file %s line %d: \n", double(duration) / 1e9, #func, __FILE__, __LINE__); \
  } while (0)

#define TIMER_SPECIFIC_ITER(func, realIter, timerIter) \
  do {                                                 \
    if (realIter == timerIter) {                       \
      TIMER(func);                                     \
    } else {                                           \
      func;                                            \
    }                                                  \
  } while (0)
#else

#define TIMER(func) func

#define TIMER_SPECIFITIMER_SPECIFIC_ITER(func, realIter, timerIter) func

#endif