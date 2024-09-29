#bin/bash
main=`pwd`
#
for file in *
do
if [ -d $file ]
then
cd $file
echo $file
#sleep 5   #waiting for 5 sec to submit next job 
tail vasp.out -n 3
cd $main
fi
done
