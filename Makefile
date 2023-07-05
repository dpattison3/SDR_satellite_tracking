all: c_utils io_utils

c_utils:
	gcc -shared -Wl,-soname,phase_locked_loop -o build/phase_locked_loop.so -fPIC util/phase_locked_loop.c -O3 -ffast-math

io_utils:
	g++ -shared -Wl,-soname,data_io -o build/data_io.so -fPIC util/data_io.cc -O3 -ffast-math -std=c++17
