#!/bin/bash -e
set -beEu -o pipefail

# wrapper to start intronProspector programs with distributed shared libraries

distexec=$(dirname $(dirname $(realpath $0)))/libexec
distprog=${distexec}/bin/$(basename $0)

export LD_LIBRARY_PATH=${distexec}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

exec ${distprog} "$@"
