all: pll_util soapy_util

pll_util:
	gcc -shared -Wl,-soname,phase_locked_loop -o build/phase_locked_loop.so -fPIC util/phase_locked_loop.c -O3 -ffast-math

soapy_util:
	g++ util/soapy_test.cc -o build/soapy_test -lSoapySDR
