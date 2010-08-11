exec=texhcod
exec=trandom

ref=0

nbTest=1e4
#3060 782 5051 5477 7019 4666 4172 6728 1994 3922 3573 5447 7116 7689

for k in `seq 0 $nbTest`
do
for j in `seq 0 $nbTest`
do
for i in `seq 0 $nbTest`
do
    ./unitTesting/$exec &> tmp.txt
    res=$?;
    seed=`cat tmp.txt | grep seed`
    rank=`printf %04d $k``printf %04d $j``printf %04d $i`
    if [ $res != "0" ]
    then
	echo $rank": "$seed " ---> MISS."
	mv tmp.txt bugrepport_$ref.txt
	ref=`expr $ref + 1 `
    else
	echo $rank": "$seed " OK."
	rm tmp.txt
    fi
done
done
done

