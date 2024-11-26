#/usr/bin/env bash
mkdir build
cd build
../configure --disable-mcs51-port  --disable-z80-port  --disable-z180-port  --disable-r2k-port  --disable-r2ka-port  --disable-r3ka-port  --disable-sm83-port  --disable-tlcs90-port  --disable-ez80_z80-port  --disable-z80n-port  --disable-ds390-port  --disable-ds400-port  --disable-pic14-port  --disable-pic16-port  --disable-s08-port  --disable-stm8-port  --disable-pdk13-port  --disable-pdk14-port  --disable-pdk15-port  --disable-mos6502-port  --disable-ucsim  --disable-packihx
IN_BUILD=1 make -j$(nproc)
cd support/regression
make test-hc08 | tee output.txt
cp /home/spl/spl_g.py ./
python spl_g.py

