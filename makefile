all: magisome_sim

hpalib:
	make -C hpalib
	make install -C hpalib

magisome_sim: magisome_sim.cpp hpalib
	g++ -g magisome_sim.cpp -L/usr/local/lib -lhpa -lm -o magisome_sim.exe

clean:
	rm magisome_sim.exe
