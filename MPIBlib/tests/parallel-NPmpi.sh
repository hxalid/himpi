#!/bin/bash

for i in capricorne-37.lyon.grid5000.fr#capricorne-39.lyon.grid5000.fr capricorne-4.lyon.grid5000.fr#capricorne-40.lyon.grid5000.fr capricorne-41.lyon.grid5000.fr#capricorne-42.lyon.grid5000.fr capricorne-43.lyon.grid5000.fr#capricorne-45.lyon.grid5000.fr ; do
a=`echo $i | sed  "s/#[a-z0-9.-]*$//"`
b=`echo $i | sed  "s/^[a-z0-9.-]*#//"`
echo $a
echo $b
mpirun -n 1 -host $a  /home/alastovetsky/KIRIL/NetPIPE-3.7.1/NPmpi -o NPmpi.$a.out : -n 1  -host $b /home/alastovetsky/KIRIL/NetPIPE-3.7.1/NPmpi  -o NPmpi.$a.out &

done
