#!/bin/bash
#NAME=/baby/one/more/time
#echo $NAME>name.log
#FULLNAME=$(cat name.log)
#echo $FULLNAME
#POS=$(dirname $FULLNAME)
#echo $POS
#NAME=$(basename $POS)
#echo $NAME

echo "++Note:"
echo "		./lzw.sh \$YourFilePath \$Yourchannel"
path=$1
pathname=$(basename $path)
inp=$2
def=0
startid=${inp:-$def}
echo "the start ID is $startid"
pre=''
flag=0

subdir_num=$(ls -l $path | grep "^d" | wc -l)
if [ $subdir_num -gt $flag ]
then

	for name in `ls $path`
	do

		echo $name
		if [ -d $path/$name ]
		then
			Num=$(ls -l $path/$name | grep "^d" | wc -l)
			nfile=$(ls -l $path/$name/wave* | grep "^-" | wc -l)
			#nfile=$(expr $(ls -l $path/$name/${channel}* | grep "^-" | wc -l) - 1)
			echo $Num
			echo $nfile

			if [ $Num -eq $flag ]
			then
				echo $pre
				echo "Nfile =" $((nfile))
				./analyze $path/$name $((nfile)) $path/$pathname$name.root $startid
				echo "./analyze $path/$name $((nfile)) $path/$pathname$name.root $startid"
			else
				for subname in `ls $path/$name`
				do
					nsubfile=$(ls -l $path/$name/$subname/wave* | grep "^-" | wc -l)
					echo "Nfile =" $((nsubfile))
					./analyze $path/$name/$subname $((nsubfile)) $name_$subname.root $startid
					echo "./analyze $path/$name/$subname $((nsubfile)) $name_$subname.root $startid"
				done
			fi
		fi

	done
else

	nfile=$(ls -l $path/wave* | grep "^-" | wc -l)
	#nfile=$(expr $(ls -l $path/${channel}* | grep "^-" | wc -l) - 1)
	echo "Nfile =" $((nfile))
	echo $pre
	./analyze $path/ $((nfile)) $path/../${pathname}.root $startid
	echo "./analyze $path/ $((nfile)) $path/../${pathname}.root $startid"
fi

