# wcsTABproj
WCSlib spatial reprojection from TAB coordinate description.

From Calabetta WCSlib utils: it takes a FITS image which spatial coordinates are described using TAB algorithm and resamples it in one of the other standard FITS projections (CAR, MER, etc.,)

It needs version 6.2 of Calabretta libwcs.

Run
```
$ gcc wcsTABproj.c -o wcsTABproj -I/usr/local/include/wcslib -g -O2 -I/usr/include/cfitsio -lcfitsio -L/usr/lib/gcc/x86_64-redhat-linux/8 -L/usr/local/lib64 -L/usr/local/lib/ -lgfortran -lquadmath -lm -lwcs
```
for compilation on Fedora29 with libwcs-6.2 installed in /usr/local
