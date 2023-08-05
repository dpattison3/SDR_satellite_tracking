#include <iostream>


template <typename T>
void print_internal(T t) {
  std::cout << t << " ";
}

template<typename T, typename... Args>
void print_internal(T t, Args... args) {
  std::cout << t << " ";
  print_internal(args...);
}

// A print function that allows calling print without manually typing std::cout and << between
// every argument. This is marginally easier to use than printf which requires specifically
// typing each argument (ex: "%d").
template<typename T, typename... Args>
void print(T t, Args... args) {
  print_internal(t, args...);
  std::cout << std::endl;
}

// Same as above but red text.
// NOTE: This isn't guaranteed to work in all terminal environments (mac?).
template<typename T, typename... Args>
void print_error(T t, Args... args) {
  // 31 is for red
  std::cout << "\033[" << 31 << "m";
  print_internal(t, args...);
  // 49 sets back to default
  std::cout << "\033[" << 39 << "m" << std::endl;
}