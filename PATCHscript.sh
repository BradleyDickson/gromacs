#!/bin/bash
nowname=`date +%s`
LOG=`echo "PATCH.LOG."$nowname`
flag=`echo $@ |grep WTmetaD |wc -l`
cflag=`echo $@ |grep -i cylinder |wc -l`
alaflag=`echo $@ |grep alanine |wc -l`

if [ $alaflag -gt 0 ] ; then
    cd src/programs/mdrun
    cp ALANINEMDdotC md.c
    cp ALANINE.f hello.f
else

if [ -z "$1" ] ; then
    echo "Please specify inputs. Run this scripts as:"
    echo "./PATCHscript.sh 480 8 4 4 6 5.5"
    echo "where:"
    echo "BMAX=480 is number of bins"
    echo "NPARTS=8 is particle number"
    echo "NCV1=4 is number of particles in CV1"
    echo "NCV2=4 is number of particles in CV2"
    echo "CVMAX=12 is largest allowable CV value"
    echo "CVREST=12 is Restraint-edge CV value"
    echo "ALL UNITS ARE nanometers as determined by GROMACS."
else
    cd src/programs/mdrun
    if [ $flag -gt 0 ] ; then 
	if [ $cflag -gt 0 ] ; then
	    sed -e 's/cWTmetaD//g' Cylinder.f > tempfile.f
	else
	    sed -e 's/cWTmetaD//g' PatchingRMSD.f > tempfile.f
	fi
    else
	if [ $cflag -gt 0 ] ; then
	    sed -e 's/cmABP//g' Cylinder.f > tempfile.f
	else
	    sed -e 's/cmABP//g' PatchingRMSD.f > tempfile.f
	fi
    fi
    cmd=`echo "sed -i 's/BMAX/"$1"/g' tempfile.f"`
    eval $cmd
    cmd=`echo "sed -i 's/NPARTS/"$2"/g' tempfile.f"`
    eval $cmd 
    cmd=`echo "sed -i 's/NCV1/"$3"/g' tempfile.f"`
    eval $cmd 
    cmd=`echo "sed -i 's/NCV2/"$4"/g' tempfile.f"`
    eval $cmd
    cmd=`echo "sed -i 's/CVMAX/"$5"/g' tempfile.f"`
    eval $cmd
    cmd=`echo "sed -i 's/CVREST/"$6"/g' tempfile.f"`
    eval $cmd
    cmd=`echo "sed -e 's/NPARTS/"$2"/g' PatchingMDdotc > tempMD.c"`
    eval $cmd 
    NPTHREE=`echo 3*$2 |bc`
    cmd=`echo "sed -i 's/NPTHREE/"$NPTHREE"/g' tempMD.c"`
    eval $cmd 
    cd ../../..

    file="./list"
    if [ -e $file ] ; then
	cd src/programs/mdrun
	atoms=`cat ../../../list`
	MINE=`echo $atoms |sed -e 's/ /, /g'`
	cmd=`echo "sed -i 's/MINE/"$MINE"/g' tempMD.c"`
	eval $cmd 
	#commit the patches here
	mv tempMD.c md.c
	mv tempfile.f hello.f
	cd ../../.. 
	echo "You are finished patching." | tee -a $LOG
	echo "Your params.in file should specify:" | tee -a $LOG
	echo "system temp in KELVIN" | tee -a $LOG
	echo "bias parameter b" | tee -a $LOG
	echo "bias parameter c" | tee -a $LOG
	echo "hill width \"a\" AS NUMBER OF BINS" | tee -a $LOG
	echo "YOUR-BOX-EDGE in nanometers" | tee -a $LOG
	echo "Cylinder or Sphere radius in nanometers" | tee -a $LOG
	echo "value of p" | tee -a $LOG
	echo "restart status: 0=not a restart, 1=restart" | tee -a $LOG
	echo '-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\' | tee -a $LOG
	echo "An working example of params.in is:" | tee -a $LOG
	echo "300.0 #" | tee -a $LOG
	echo "0.9" | tee -a $LOG
	echo "0.1" | tee -a $LOG
	echo "10.0" | tee -a $LOG
	echo "7.0265" | tee -a $LOG
	echo "1.3" | tee -a $LOG
	echo "20.0" | tee -a $LOG
	echo "0" | tee -a $LOG
	echo "You patched with options: " $@ | tee -a $LOG
	echo "All of this output was recorded in the file " $LOG
	echo "...the numeric part is the date in unix time."
        echo "...the date in human time is recovered by \"date --date='@NUMERIC-PART'\""
	echo "...to produce: "
        cmd=`echo "date --date='@"$nowname"'"`
	eval $cmd
    else
	echo "You need to make an atom list file."
	echo "Just write your biased atoms (the indexes) into a file called list\n either as a row or column."
    fi
fi
fi
