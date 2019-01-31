#!/bin/bash

addToDrat(){


  ## Set up Repo parameters
  git init
  git config user.name "Richard White"
  git config user.email "r@rwhite.no"
  git config --global push.default simple

  ## Get drat repo
  git remote add upstream "https://$GITHUB_PAT@github.com/folkehelseinstituttet/drat.git"
  git fetch upstream 2>err.txt
  git checkout gh-pages

  Rscript -e "drat::insertPackage('$PKG_REPO/$PKG_TARBALL', \
    repodir = '.', \
    commit='Travis update $PKG_REPO: build $TRAVIS_BUILD_NUMBER')"
  Rscript -e "saveRDS(read.dcf('src/contrib/PACKAGES'),'src/contrib/PACKAGES.rds')"
  git commit -a -m "Travis update $PKG_REPO: build $TRAVIS_BUILD_NUMBER"
  git push 2>err.txt

}


set -o errexit -o nounset
PKG_REPO=$PWD
cd ..
mkdir drat; cd drat

addToDrat

rm $PKG_REPO/$PKG_TARBALL
cd $PKG_REPO
