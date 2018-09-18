suffix='1_50'
numberOfExp=10

allTime=0
for i in $(seq 1 $numberOfExp)
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null
	timeOutput=`(time taskset 0x40 gtc view -b -c 0 archive$suffix > /dev/null) 2>&1 | grep real` 
	minutes=`printf "$timeOutput" | sed -e 's/.*real\(.*\)m.*/\1/'`
	seconds=`printf "$timeOutput" | sed -e 's/.*m\(.*\)s.*/\1/'` 
	time=`python -c "print($minutes*60 + $seconds)"`
	allTime=`python -c "print($allTime + $time)"`
	printf "$time "
done
printf '\nDecompressing Time: ' 
python -c "print($allTime/$numberOfExp)"
printf "\n"


allTime=0
for i in $(seq 1 $numberOfExp)
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null
	timeOutput=`(time taskset 0x40 gtc view -b -c 0 -r 22:16976761-17000000 archive$suffix > /dev/null) 2>&1 | grep real` 
	minutes=`printf "$timeOutput" | sed -e 's/.*real\(.*\)m.*/\1/'`
	seconds=`printf "$timeOutput" | sed -e 's/.*m\(.*\)s.*/\1/'` 
	time=`python -c "print($minutes*60 + $seconds)"`
	allTime=`python -c "print($allTime + $time)"`
	printf "$time "
done
printf '\nQuery for 23,239 Bases: ' 
python -c "print($allTime/$numberOfExp)"
printf "\n"


allTime=0
for i in $(seq 1 $numberOfExp)
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null
	timeOutput=`(time taskset 0x40 gtc view -b -c 0 -r 22:16050075-17050075 archive$suffix > /dev/null) 2>&1 | grep real` 
	minutes=`printf "$timeOutput" | sed -e 's/.*real\(.*\)m.*/\1/'`
	seconds=`printf "$timeOutput" | sed -e 's/.*m\(.*\)s.*/\1/'` 
	time=`python -c "print($minutes*60 + $seconds)"`
	allTime=`python -c "print($allTime + $time)"`
	printf "$time "
done
printf '\nQuery for 10^6 Bases: ' 
python -c "print($allTime/$numberOfExp)"	
printf "\n"


allTime=0
for i in $(seq 1 $numberOfExp)
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null
	timeOutput=`(time taskset 0x40 gtc view -b -c 0 -r 22:17135480-17135480 archive$suffix > /dev/null) 2>&1 | grep real` 
	minutes=`printf "$timeOutput" | sed -e 's/.*real\(.*\)m.*/\1/'`
	seconds=`printf "$timeOutput" | sed -e 's/.*m\(.*\)s.*/\1/'` 
	time=`python -c "print($minutes*60 + $seconds)"`
	allTime=`python -c "print($allTime + $time)"`
	printf "$time "
done
printf '\nQuery for A Single Variant: ' 
python -c "print($allTime/$numberOfExp)"	
printf "\n"


00allTime=0
for i in $(seq 1 $numberOfExp)
do
	sudo bash -c 'free && sync && echo 3 > /proc/sys/vm/drop_caches && free' > /dev/null
	timeOutput=`(time taskset 0x40 gtc view -b -c 0 -s S_175000 archive$suffix > /dev/null) 2>&1 | grep real` 
	minutes=`printf "$timeOutput" | sed -e 's/.*real\(.*\)m.*/\1/'`
	seconds=`printf "$timeOutput" | sed -e 's/.*m\(.*\)s.*/\1/'` 
	time=`python -c "print($minutes*60 + $seconds)"`
	allTime=`python -c "print($allTime + $time)"`
	printf "$time "
done
printf '\nQuery for A Single Sample: ' 
python -c "print($allTime/$numberOfExp)"
printf "\n"
