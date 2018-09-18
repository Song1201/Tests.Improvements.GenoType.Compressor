totalRuntime=0
for i in {1..10}
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null  
	allTime=`(time taskset 0x40 gtc view -b -c 0 archive18 > /dev/null) 2>&1 | grep -E "user|sys" | sed s/[a-z]//g` 
	Runtime=0 
	for j in $allTime; do Runtime=`python -c "print($Runtime+$j)"`; done 
	totalRuntime=`python -c "print($totalRuntime + $Runtime)"`	
done 
printf 'Decompressing Time: ' 
printf `python -c "print('{0:.3f}'.format($totalRuntime/10))"`'\n\n'


totalRuntime=0
for i in {1..10}
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null  
	allTime=`(time taskset 0x40 gtc view -b -c 0 -r 22:16050075-17050075 archive18 > /dev/null) 2>&1 | grep -E "user|sys" | sed s/[a-z]//g` 
	Runtime=0 
	for j in $allTime; do Runtime=`python -c "print($Runtime+$j)"`; done 
	totalRuntime=`python -c "print($totalRuntime + $Runtime)"`	
done 
printf 'Query for 10^6 Bases: ' 
printf `python -c "print('{0:.3f}'.format($totalRuntime/10))"`'\n\n'


