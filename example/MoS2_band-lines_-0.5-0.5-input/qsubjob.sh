#bin/bash
main=`pwd`
#
for file in *
do
if [ -d $file ]
then
cd $file
echo $file
cp ../vasp.sh .
qsub vasp.sh
sleep 5   #waiting for 5 sec to submit next job 
cd $main
fi
done
