#!/usr/bin/env bash
for f in $1/*.fasta; do (cat "${f}"; echo) >> $1/multi.fasta; done
