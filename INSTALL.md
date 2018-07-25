

- Requirements
  - a modern C++ compiler.  Known to work:
    - GCC 6.3
  
  - htslib
    - from https://github.com/samtools/htslib or package
    - can compile against either built source tree or installed package

- Building
  ./configure --with-htslib=<htsdir>
  make
  make test
  
  
- Developer
  The configure script and built documentation files are checked into the source tree,
  as they are only modified by developers, with many users wanting to build from a
  github checkout.  
  
  To rebuild the configure script, GNU autoconf is required, then run:
     `autoconf`
     
  To rebuild the documentation, pandoc is require.  Rebuild documentation
  with:
      `make doc`
  
  


  
  
