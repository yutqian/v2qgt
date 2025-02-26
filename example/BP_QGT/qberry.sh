#bin/bash
main=`pwd`
#
for file in *
do
if [ -d $file ]
then
cd $file
echo $file
#Quantum_metric
v2qgt=/path/yourcode
$v2qgt  -ii 1 -if 10 -kubo 1 -qm 1 -kp 1 -s 1 > berryout
sleep 3   #waiting for 5 sec to submit next job 
cd $main
fi
done
