# Development

## Making a release:
- update CHANGELOG.md
- update version in configure.ac
- autoreconf -fi
- ./configure
- make doc
  this builds help and manual pages from markdown
- test:
  - make distclean
  - ./configure --prefix=/some/test/dir
  - make -j 32
  - make -j 32 test
  - make install
- git commit -am 'v1.5.0 release'
- git tag v1.5.0
- git push
- git push --tags
- github draft a release with CHANGELOG text as description
