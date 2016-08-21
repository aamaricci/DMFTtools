#!/bin/bash
NAME=DMFTT
UNAME=`echo $NAME |tr [:lower:] [:upper:]`
LNAME=`echo $NAME |tr [:upper:] [:lower:]`
LIBNAME=lib$LNAME.a
LOG=install.log
>$LOG
exec >  >(tee -a $LOG)
exec 2> >(tee -a $LOG >&2)

echo ""
echo "Invoked as:"
echo "$0 $*"
echo ""

USER=$(id -un)
GUSER=$(id -gn)


#>>> GET THE ACTUAL DIRECTORY
HERE=$(pwd)


#>>> USAGE FUNCTION
usage(){
    echo ""
    echo "usage:"
    echo ""
    echo "$0  --plat=FC_PLAT --prefix=PREFIX_DIR  --sfprefix=SciFor_root ($SFROOT) [  -d,--debug   -h,--help   -m,--no-mpi   -q,--quiet ]"
    echo ""
    echo ""
    echo "mandatory arguments:" 
    echo "    --plat       : specifies the actual platform/compiler to use [intel,gnu]"
    echo "    --prefix     : specifies the target directory "
    echo "    --sfprefix   : specifies the SciFor root directory "
    echo ""    
    echo "optional arguments:"
    echo "    -d,--debug  : debug flag"
    echo "    -h,--help   : this help"
    echo "    -m,--no-mpi : do not use MPI compiler: produce a serial library"
    echo "    -q,--quiet  : assume Y to all questions."
    echo ""
    exit
}



#>>> GET Nth TIMES PARENT DIRECTORY
nparent_dir(){
    local DIR=$1
    local N=$2
    for i in `seq 1 $N`;
    do 
	DIR=$(dirname $DIR)
    done
    echo $DIR
}

#>>> GET THE ENTIRE LIST OF ARGUMENTS PASSED TO STDIN
LIST_ARGS=$*

#>>> GET LONG & SHORT OPTIONS
params="$(getopt -n "$0" --options dhmq --longoptions plat:,prefix:,sfprefix:,debug,help,no-mpi,quiet -- "$@")"
if [ $? -ne 0 ];then
    usage
fi
eval set -- "$params"
unset params

#>>> CHECK THE NUMBER OF ARGUMENTS. IF NONE ARE PASSED, PRINT HELP AND EXIT.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
    usage
fi


#>>> SET SOME DEFAULTS VARIABLES AND OTHER ONES. 0=True; 1=False
WPLAT=1
DEBUG=1
QUIET=1
WMPI=0				# assume with MPI (WMPI=True)
VERSION=$(git describe --tags 2>/dev/null)
WRK_INSTALL=$(pwd)
BIN_INSTALL=$WRK_INSTALL/bin
ETC_INSTALL=$WRK_INSTALL/etc
OPT_INSTALL=$WRK_INSTALL/opt
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/src


#>>> THE LISTS OF ALLOWED PLAT
LIST_FC="gnu intel"


#>>> GO THROUGH THE INPUT ARGUMENTS. FOR EACH ONE IF REQUIRED TAKE ACTION BY SETTING VARIABLES.
while true
do
    case $1 in
	--plat)
	    WPLAT=0
	    PLAT=$2
	    shift 2
	    [[ ! $LIST_FC =~ (^|[[:space:]])"$PLAT"($|[[:space:]]) ]] && {
		echo "Incorrect Fortran PLAT: $PLAT";
		echo " available values are: $LIST_FC"
		exit 1
	    }
	    ;;
	--prefix)
	    PREFIX=$2;
	    shift 2
	    ;;
	--sfprefix)
	    SFROOT=$2
	    shift 2
	    ;;
	-d|--debug) DEBUG=0;shift ;;
        -h|--help) usage ;;
	-m|--no-mpi) WMPI=1;shift ;;
	-q|--quiet) QUIET=0;shift ;;
        --) shift; break ;;
        *) usage ;;
    esac
done

#>>> CHECK THAT THE MANDATORY OPTION -p,-plat IS PRESENT:
[[ $WPLAT == 0 ]] && [[ ! -z $PREFIX ]] || usage
[[ ! -z $SFROOT ]]  || usage


case $PLAT in
    gnu) _FC_=gfortran ;;
    intel) _FC_=ifort ;;
esac


#RENAME WITH DEBUG IF NECESSARY 
[[ $DEBUG == 0 ]] && PLAT=${PLAT}_debug
[[ $WMPI == 1 ]] && PLAT=${PLAT}_nompi

#>>> SET STANDARD NAMES FOR THE TARGET DIRECTORY
DIR_TARGET=$PREFIX/$PLAT
BIN_TARGET=$DIR_TARGET/bin
ETC_TARGET=$DIR_TARGET/etc
LIB_TARGET=$DIR_TARGET/lib
INC_TARGET=$DIR_TARGET/include
DIR_TARGET_W=0
if [ -d $PREFIX ];then
    TEST_W_DIR=$PREFIX
else
    TEST_W_DIR=$(nparent_dir $PREFIX 1)
fi
if [ ! -w $TEST_W_DIR ];then
    DIR_TARGET_W=1
    echo "Can not create $DIR_TARGET: $TEST_W_DIR has no write access"
    sleep 1
    echo "Try to grant root privileges to create $DIR_TARGET for $USER:$GROUP"
    sudo -v
fi

#TEST SCIFOR DIRECTORY EXISTS
if [ ! -d $SFROOT ];then echo "$0: can not find SciFor root directory at $SCIFOR";exit;fi
if [ ! -d $SFROOT/$PLAT ];then echo "$0: can not find $SFROOT/$PLAT directory";exit;fi


#>>> TEST MPI
if [ $WMPI == 0 ];then
    MPIFC=mpif90
    which $MPIFC >/dev/null 2>&1
    CHECK_MPIFC=$?
    if [ $CHECK_MPIFC == 0 ];then
	echo "Using the following MPI compiler:"
	mpif90 -show
	sleep 0.5
	echo "On platform: $PLAT"
	sleep 0.5
	_MPIFC_=$(mpif90 -show | awk '{print $1}') # get the MPI f90 compiler
	if [ "$_MPIFC_" != "$_FC_" ];then
	    echo "Possible mismatch between MPI:+$_MPIFC_ and Platform $PLAT:+$_FC_"
	    echo "Check your MPI installation... exiting"
	    sleep 2
	    exit 1
	fi
    else
	echo "Can not find the $MPIFC compiler in the system."
	sleep 0.5
	echo "Trying installing with the --no-mpi option"
	sleep 2
	exit 1 
    fi
    CPP_FLAG='MPI'
else
    CPP_FLAG=''
fi

create_makeinc(){
    local PLAT=$1
    cd $WRK_INSTALL
    case $PLAT in
	intel)
	    local FC=ifort
	    local FCOPT="-O2 -ftz -static-intel -fpp -D_$CPP_FLAG "
	    local MOPT="-module "
	    ;;
	gnu)
	    local FC=gfortran
	    local FCOPT="-O2 -funroll-all-loops -static -cpp -D_$CPP_FLAG "
	    local MOPT=-J
	    ;;
	intel_debug)
	    local FC=ifort
	    local FCOPT="-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel -fpp -D_$CPP_FLAG "
	    local MOPT="-module "
	    ;;
	gnu_debug)
	    local FC=gfortran
	    local FCOPT="-O0 -p -g -Wall -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -static -cpp -D_$CPP_FLAG "
	    local MOPT=-J
	    ;;
	*)
	    usage
	    ;;
    esac

    if [ $WMPI == 0 ];then
	FC=$MPIFC		# set Fortran Compiler to the MPI Compiler (we checked it exists)
    fi
    
    cat << EOF > make.inc
FC=$FC
MPIFC=$MPIFC
FCOPT=$FCOPT
MOPT=$MOPT
PLAT=$PLAT
LIB_DMFTT=$LIB_TARGET/libdmftt.a
INC_SCIFOR=$SFROOT/$PLAT/include
INC_TARGET=$INC_TARGET
EOF
}


#>>> GET THE ACTUAL DIRECTORY
HERE=$(pwd)


create_makeinc $PLAT
sleep 2



if [ $QUIET == 1 ];then		# if NO quiet
    _DIR=Y
    echo -n "Installing in $DIR_TARGET. Continue [Y/n]: "
    read _DIR;
    _DIR=`echo $_DIR |tr [:lower:] [:upper:]`
    [[ $_DIR == Y ]] || exit 1
else				# if quiet 
    echo "Installing DMFT_Tools in $DIR_TARGET (quiet mode): "
    sleep 2
fi


# >>> CREATE THE DIRECTORY HIERARCHY:
echo "Creating directories:"
sleep 1
if [ $DIR_TARGET_W -eq 0 ];then
    echo "mkdir -pv $DIR_TARGET"
    mkdir -pv $DIR_TARGET
else
    echo "sudo mkdir -pv $DIR_TARGET && sudo chown $USER:$GROUP $DIR_TARGET"
    sudo mkdir -pv $DIR_TARGET && sudo chown $USER:$GROUP $DIR_TARGET
fi
mkdir -pv $BIN_TARGET
mkdir -pv $ETC_TARGET/modules/$LNAME
mkdir -pv $LIB_TARGET
mkdir -pv $INC_TARGET
sleep 2


echo "Copying init script for $UNAME" 
cp -fv $BIN_INSTALL/configvars.sh $BIN_TARGET/configvars.sh
cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${PREFIX}/${PLAT}
EOF
echo "" 
sleep 1

echo "Generating environment module file for $UNAME" 
cat <<EOF > $ETC_TARGET/modules/$LNAME/$PLAT
#%Modules
set	root	$PREFIX
set	plat	$PLAT
set	version	"($PLAT)"
EOF
cat $ENVMOD_INSTALL/module >> $ETC_TARGET/modules/$LNAME/$PLAT
echo "" 
sleep 1

echo "Compiling $UNAME library on platform $PLAT:"
echo "" 
sleep 1






#rm -fv $LIB_TARGET/$LIBNAME
make all
sleep 1
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.dmfttols
else
    echo "Error from Makefile. STOP here."
    exit 1
fi


#LAST TOUCH COPY THE CONFIGVARS AND CREATE THE USER MODULES FILE. PRINT USAGE DETAILS.
CONFIGFILE=$PREFIX/$PLAT/bin/configvars.sh
MODULEFILE=$PREFIX/$PLAT/etc/modules/$LNAME/$PLAT
mkdir -pv $HOME/.modules.d/other/$LNAME
cp -vf $MODULEFILE $HOME/.modules.d/other/$LNAME/$PLAT


echo "" 
echo "USAGE:" 
echo "" 
echo " To add DMFT_Tools to your system:" 
echo "   $ source $CONFIGFILE" 
echo " (or copy this line into your bash profile [e.g. .bashrc])" 
echo ""
module avail >/dev/null 2>&1 
if [ $? == 0 ];then
    echo " or load the $UNAME modules:" 
echo "   $ module use $HOME/.modules.d/applications" 
echo "   $ module load $LNAME/$PLAT" 
echo "(or copy these lines into your bash profile [e.g. .bashrc])" 
echo ""
fi
echo ""
echo "Enjoy... (for info: adriano.amaricciATgmail.com)"
echo ""


exit 0
