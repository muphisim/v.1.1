#!/bin/sh
parentdir=$PWD
echo $parentdir

VAR1="export PYTHONPATH=$parentdir:\$PYTHONPATH"
echo "$VAR1" >> $HOME/.bashrc
VAR2="export PATH=$parentdir:\$PATH"
echo "$VAR2" >> $HOME/.bashrc
. $HOME/.bashrc
