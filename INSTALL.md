

- Requirements
  - a modern C++ compiler.  Known to work:
    - GCC 6.3
  
  - htslib
    - from https://github.com/samtools/htslib or package
    - found with pkg-config by default
    - can explicitly specify htslib install directory --with-htslib=<htsdir>

   - running tests require
    - samtools program
    - bedToBigBed from UCSC browser

- Building
  ./configure
  make
  make install

- Developer
  The configure script and built documentation files are checked into the source tree,
  as they are only modified by developers, with many users wanting to build from a
  github checkout.
  
  To rebuild the configure script, GNU autoconf is required, then run:
     `autoreconf -fi`
     
  To rebuild the documentation, pandoc is require.  The generated
  files are check into the tree.   Rebuild documentation
  creates built-in help files (src/*.man.h) and man pages (man/*.1)
  and is run with:
      `make doc`
