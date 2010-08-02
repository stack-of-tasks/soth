exec=texhcod
exec=trandom

ref=0

nbTest=1e3
for i in `seq 1 $nbTest`
do
    ./unitTesting/$exec &> tmp.txt
    res=$?;
    seed=`cat tmp.txt | grep seed`
    if [ $res != "0" ]
    then
	echo $i": "$seed " ---> MISS."
	mv tmp.txt $ref.txt
	ref=`expr $ref + 1 `
    else
	echo $i": "$seed " OK."
	rm tmp.txt
    fi
done

