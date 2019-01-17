#/bin/bash
[[ ! -d eps256i ]] && mkdir eps256i
for exe in dkd kdk kdkdk2 kdkdk4 dkdc kdkc kdkdk2c kdkdk4c
do
  rm -f eps256i/$exe.dat
  for tick in 512 1024 2048 4096 8192 16384
  do
    ./$exe pl1k.dat $tick | grep err >> eps256i/$exe.dat
  done
done
