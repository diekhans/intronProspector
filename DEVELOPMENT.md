# Development

Notes:
* The configure script is kept current it git, which allows builds from github
  generated release download archives.

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
- in fork of https://github.com/bioconda/bioconda-recipes
  update the version and sha256 in
     bioconda-recipes/recipes/intron-prospector/meta.yaml
  make a pull request to bioconda/bioconda-recipes

