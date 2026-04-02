# Development

## Making a release:
- update CHANGELOG.md
- update version in configure.ac
- autoreconf -fi
- make doc 
- test:
  - make distclean
  - ./configure --prefix=${MED_OPT}
  - make
  - make install
  - make test
- git commit -am 'v1.5.0 release'
- tag in the form v1.5.0
