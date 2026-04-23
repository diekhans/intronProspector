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
- git commit -am '1.x.x release'
- git tag v1.x.x
- git push
- git push --tags
- github draft a release with CHANGELOG text as description
- BioConda
  - first time or dependencies changes
    in fork of https://github.com/bioconda/bioconda-recipes
    update the version and sha256 in
       bioconda-recipes/recipes/intron-prospector/meta.yaml
    make a pull request to bioconda/bioconda-recipes
  - wait for and  automatic pull request to be generate
    Once the tests have passed, add this comment to the pull request:
	    @BiocondaBot please add label
    This should take care of the red 'Review Required' and 'Merging is Blocked' notifications

