#!/bin/bash

set -e

# create a temporary directory
tmpdir1=$(mktemp -d)
tmpdir2=$(mktemp -d)

function finish()
{
    rm -rf "${tmpdir1}"
    rm -rf "${tmpdir2}"
}

trap finish EXIT

# create HTML files
cd "${tmpdir1}"
git clone https://github.com/cschwan/hep-mc.git
cd hep-mc
meson build -Ddoxygen=true
cd build
ninja

version=$(git describe)

cd "${tmpdir2}"
git clone --depth 1 --branch gh-pages https://github.com/cschwan/hep-mc.git
cd hep-mc
git rm -rf *.{css,html,js,png} search
mv "${tmpdir1}"/hep-mc/build/doc/html/* .
git add .

if [[ $(git commit -m "Update to ${version}") != 0 ]]; then
    echo "No new commits"
fi
