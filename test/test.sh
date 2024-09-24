#Code to run all tests in test/bin/
#N.B. test should end by .x

set -e

if [ ! -d "bin" ]
then
    echo "\e[31m ERROR \e[0m"
    echo " There is no *bin* directory"
    echo " Try  'make all' before testing"
    return 1
fi

cd bin/


HERE=`pwd`
echo $HERE

while read DIR; do
    if [ -d $DIR ]; then
	if ls $DIR/*.x  1> /dev/null 2>&1; then
	    echo "TESTING $DIR"
	    cd $DIR
	    pwd
	    for exe in *.x
	    do
		echo "Running $exe:"
		time ./$exe
		echo ""
	    done
	    cd $HERE
	fi
    fi
done<list_dir
    

