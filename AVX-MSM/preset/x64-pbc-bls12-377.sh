# #!/bin/sh
cmake -DWSIZE=64 -DRAND=UDEV -DSHLIB=OFF -DSTBIN=ON -DTIMER=HREAL -DCHECK=off -DVERBS=off -DARITH=x64-asm-6l -DFP_PRIME=377 -DFP_METHD="INTEG;INTEG;INTEG;MONTY;JMPDS;JMPDS;SLIDE" -DCFLAGS="-O3 -funroll-loops -fomit-frame-pointer -finline-small-functions -march=native -mtune=native" -DFP_PMERS=off -DFP_QNRES=off -DFPX_METHD="INTEG;INTEG;LAZYR" -DEP_PLAIN=off -DEP_SUPER=off -DPP_METHD="LAZYR;OATEP" $1

#!/bin/sh
# cmake -DWSIZE=64 -DRAND=UDEV -DSHLIB=OFF -DSTBIN=ON -DTIMER=CYCLE -DCHECK=off -DVERBS=off  -DFP_PRIME=377 -DFP_METHD="INTEG;INTEG;INTEG;MONTY;JMPDS;JMPDS;SLIDE" -DCFLAGS="-O3 -funroll-loops -fomit-frame-pointer -finline-small-functions -march=native -mtune=native" -DFP_PMERS=off -DFP_QNRES=off -DFPX_METHD="INTEG;INTEG;LAZYR" -DEP_PLAIN=off -DEP_SUPER=off -DPP_METHD="LAZYR;OATEP" $1
