#bin/bash
main=`pwd`
#
for file in *
do
if [ -d "$file" ] && [ "$file" != "BC_result" ]; then
cd $file
echo $file
#rm BERRYCURV_KUBO* BerryCurv/
#mkdir BerryCurv
#mv BERRYCURV_KUBO*  berryout BerryCurv/
mkdir QM
mv BERRYCURV_KUBO*  berryout QM/
v2qgt=/ioss/ytqian/soft/v2qgt/v2qgt_v2/v2qgt.x
$v2qgt -kx 11 -ky 11 -ii 1 -if 18 -kubo 1 -kp 1 > berryout
#$v2qgt -kx 11 -ky 11 -ii 1 -if 18 -kubo 1 -qm 1 -kp 1 > berryout
sleep 3   #waiting for 5 sec to submit next job 
cd $main
fi
done
