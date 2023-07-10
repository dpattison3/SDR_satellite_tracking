all: pll_util sdr_record

pll_util:
	gcc -shared -Wl,-soname,phase_locked_loop -o build/phase_locked_loop.so -fPIC util/phase_locked_loop.c -O3 -ffast-math

sdr_record:
	g++ util/sdr_record.cc -o build/sdr_record -lSoapySDR
