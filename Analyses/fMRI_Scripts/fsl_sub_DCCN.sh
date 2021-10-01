#!/bin/sh

# Copyright (C) 2007 University of Oxford
# Authors: Dave Flitney & Stephen Smith

#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 4.0 (c) 2007, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/1112.






###########################################################################
# Edit this file in order to setup FSL to use your local compute
# cluster.
###########################################################################


###########################################################################
# The following section determines what to do when fsl_sub is called
# by an FSL program. If it finds a local cluster if will pass the
# commands onto the cluster. Otherwise it will run the commands
# itself. There are three values for the METHOD variable, "SGE", "TORQUE" and
# "NONE". You should setup the tests to look for whether the calling
# computer can see your cluster setup scripts, and run them (if that's
# what you want, i.e. if you haven't already run them in the user's
# login scripts). Note that these tests look for the environment
# variable SGE_ROOT, which a user can unset if they don't want the
# cluster to be used.
###########################################################################

METHOD=SGE
if [ "x$SGE_ROOT" = "x" ] ; then
    if [ -f /usr/local/share/sge/default/common/settings.sh ] ; then
	. /usr/local/share/sge/default/common/settings.sh
    elif [ -f /usr/local/sge/default/common/settings.sh ] ; then
	. /usr/local/sge/default/common/settings.sh
    elif [ -f /opt/cluster/old/torque.sh ] ; then
	. /opt/cluster/old/torque.sh
	METHOD=TORQUE
        WALLTIME="86400"
        MailOpts="ab"
    else
	METHOD=NONE
    fi
fi

###########################################################################
# The following auto-decides what cluster queue to use. The calling
# FSL program will probably use the -T option when calling fsl_sub,
# which tells fsl_sub how long (in minutes) the process is expected to
# take (in the case of the -t option, how long each line in the
# supplied file is expected to take). You need to setup the following
# list to map ranges of timings into your cluster queues - it doesn't
# matter how many you setup, that's up to you.
###########################################################################

map_qname ()
{
    if [ $1 -le 20 ] ; then
	queue=veryshort.q
    elif [ $1 -le 120 ] ; then
	queue=short.q
    elif [ $1 -le 1440 ] ; then
	queue=long.q
    else
	queue=verylong.q
    fi
    #echo "Estimated time was $1 mins: queue name is $queue"

}


###########################################################################
# Don't change the following (but keep scrolling down!)
###########################################################################

POSIXLY_CORRECT=1
export POSIXLY_CORRECT
command=`basename $0`

usage ()
{
  cat <<EOF

$command V1.0beta - wrapper for job control system such as SGE or TORQUE

Usage: $command [options] <command>

$command gzip *.img *.hdr
$command -q short.q gzip *.img *.hdr
$command -a darwin regscript rawdata outputdir ...

  -T <minutes>          Estimated job length in minutes, used to auto-set queue name (default 20 minutes)
  -q <queuename>        Possible values for <queuename> are "batch", "verylong.q", "long.q" 
                        and "short.q". See below for details
                        Default is "batch".
  -M <email-address>    Who to email, default = none 
  -j <jid>              Place a hold on this task until job jid has completed
  -t <filename>         Specify a task file of commands to execute in parallel
  -N <jobname>          Specify jobname as it will appear on queue
  -n <nCPUs>            Number of CPUs that job will use
  -l <logdirname>       Where to output logfiles
  -m <mailoptions>      Change the notification mail options, see qsub for details
  -F                    Verbose mode.
  -u			Memory setting for qsub usage (No default set!)

Queues:

There are several batch queues configured on the cluster, each with defined CPU
time limits. In almost every case, if the job isn't a matlab job it should be run in the batch queue.

EOF

  exit 1
}

nargs=$#
if [ $nargs -eq 0 ] ; then
  usage
fi

set -- `getopt T:q:M:j:t:N:n:Fvm:l:ru: $*`
result=$?
if [ $result != 0 ] ; then
  echo "What? Your arguments make no sense!"
fi

if [ $nargs -eq 0 ] || [ $result != 0 ] ; then
  usage
fi


###########################################################################
# The following sets up the default queue name, which you may want to
# change. It also sets up the basic emailing control.
###########################################################################

queue=batch
WALLTIME=24000 #1200
mailto=`whoami`
MailOpts="n"
TORQUEDEPENDANCYMODE="w" #is one of "w" "l"
# memoryLimit="10gb"

###########################################################################
# In the following, you might want to change the behaviour of some
# flags so that they prepare the right arguments for the actual
# cluster queue submission program, in our case "qsub".
#
# -a sets is the cluster submission flag for controlling the required
# hardware architecture (normally not set by the calling program)
#
# -p set the priority of the job - ignore this if your cluster
# environment doesn't have priority control in this way.
#
# -j tells the cluster not to start this job until cluster job ID $jid
# has completed. You will need this feature.
#
# -t will pass on to the cluster software the name of a text file
# containing a set of commands to run in parallel; one command per
# line.
#
# -N option determines what the command will be called when you list
# running processes.
#
# -l tells the cluster what to call the standard output and standard
# -error logfiles for the submitted program.
###########################################################################

if [ -z $FSLSUBVERBOSE ] ; then
    verbose=0
else
    verbose=$FSLSUBVERBOSE;
    echo "METHOD=$METHOD : args=$@" >&2
fi

scriptmode=0

while [ $1 != -- ] ; do
  case $1 in
    -T)
	  time=$2
      WALLTIME=`echo "60 * $time" | bc -l`
	  shift;;
    -q)
      queue=$2
 	  shift;;
    -M)
      mailto=$2
      if [ $METHOD = "TORQUE" ] ; then
        mailto=`whoami`@wusm.wustl.edu
      fi
      shift;;
    -j)
      jid="$2"
      if [ $TORQUEDEPENDANCYMODE = "l" ] ; then
		length1=${#jid}
		temp=${jid%[]}
		length2=${#temp}
		if [ $length1 -eq $length2 ] ; then
			torque_hold=",depend=afterok:$jid"
		else
			torque_hold=",depend=afterokarray:$jid"
		fi
      elif [ $TORQUEDEPENDANCYMODE = "w" ] ; then
        length1=${#jid}
		temp=${jid%[]}
		length2=${#temp}
		if [ $length1 -eq $length2 ] ; then
			torque_hold=" -W depend=afterok:$jid"
		else
			torque_hold=" -W depend=afterokarray:$jid"
		fi
	  fi
      sge_hold="-hold_jid $jid"
      shift;;
    -t)
      taskfile=$2
      tasks=`wc -l $taskfile | awk '{print $1}'`
      sge_tasks="-t 1-$tasks"
      shift;;
    -N)
      JobName=$2;
      shift;;
    -n)
      if [ $METHOD = "SGE" ] ; then
        NumCPUs="-pe make $2"
      elif [ $METHOD = "TORQUE" ] ; then
        NumCPUs=",nodes=1:ppn=$2"
      fi
      shift;;
    -m)
      MailOpts=$2;
      if [[ $METHOD = "TORQUE" && $MailOpts = "as" ]] ; then
        MailOpts=ab
      fi 
      shift;;
    -l)
      LogOpts="-o $2 -e $2";
      LogDir="${2}/";
      mkdir -p $2;
      shift;;
    -u)
	  memoryLimit=$2;
	  shift;;
    -v)
      verbose=1
      ;;
	
  esac
  shift  # next flag
 
done
shift

###########################################################################
# Make sure a job runs locally when memory is not set for DCCN torque compatibility
###########################################################################

if [ "x$memoryLimit" = x ] ; then
	METHOD=NONE
	echo "running locally"
fi


###########################################################################
# Don't change the following (but keep scrolling down!)
###########################################################################

if [ "x$JobName" = x ] ; then 
    if [ "x$taskfile" != x ] ; then
	JobName=`basename $taskfile`
    else
	JobName=`basename $1`
    fi
fi

if [ "x$tasks" != x ] && [ ! -f "$taskfile" ] ; then
    echo $taskfile: invalid input!
    echo Should be a text file listing all the commands to run!
    exit -1
fi

if [ "x$tasks" != "x" ] && [ "x$@" != "x" ] ; then
    echo $@
    echo Spurious input after parsing command line!
    exit -1
fi

## REMOVE after testing
# verbose=1

case $METHOD in

###########################################################################
# The following is the main call to the cluster, using the "qsub" SGE
# program. If $tasks has not been set then qsub is running a single
# command, otherwise qsub is processing a text file of parallel
# commands.
###########################################################################

    SGE)
	if [ "x$tasks" = "x" ] ; then
	    if [ $scriptmode -ne 1 ] ; then
		sge_command="qsub -V -cwd -shell n -b y -r y -q $queue -M $mailto -N $JobName -m $MailOpts $LogOpts $sge_arch $sge_hold $NumCPUs"
	    else
		sge_command="qsub $LogOpts $sge_arch $sge_hold $NumCPUs"
	    fi
	    if [ $verbose -eq 1 ] ; then 
		echo sge_command: $sge_command >&2
		echo executing: $@ >&2
	    fi
	    exec $sge_command $@ | awk '{print $3}'
	else
	    sge_command="qsub -V -cwd -q $queue -M $mailto -N $JobName -m $MailOpts $LogOpts $sge_arch $sge_hold $sge_tasks $NumCPUs"
	    if [ $verbose -eq 1 ] ; then 
		echo sge_command: $sge_command >&2
		echo control file: $taskfile >&2
	    fi
	    exec $sge_command <<EOF | awk '{print $3}' | awk -F. '{print $1}'
            

#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskfile\`

exec /bin/sh -c "\$command"
EOF
	fi
	;;

###########################################################################
# The following is the main call to the cluster, using the "qsub" TORQUE
# program. If $tasks has not been set then qsub is running a single
# command, otherwise qsub is processing a text file of parallel
# commands. This script is compatible with MOAB 5.3.7.s15113 and higher.
###########################################################################

    TORQUE)
    # SGE takes args after script, TORQUE does not. Tempscript stores the command and arguements to be used
    
    scriptbase=$( basename $0 )
    tempscript=$( mktemp ${scriptbase}.XXXXXXXX ) || exit 1
    queue=$( basename $queue .q )

 	if [ "x$tasks" = "x" ] ; then
            echo $@ > $tempscript
	    if [ $scriptmode -ne 1 ] ; then
                torque_command="qsub -V -d . -b y -r y -q $queue -M $mailto -N $JobName -m $MailOpts -l mem=$memoryLimit $LogOpts $sge_arch -l walltime=$WALLTIME$NumCPUs$torque_hold"	
            else
		torque_command="qsub -l mem=$memoryLimit $LogOpts $sge_arch -l walltime=$WALLTIME$NumCPUs$torque_hold"
	    fi
	    if [ $verbose -eq 1 ] ; then 
		echo torque_command: $torque_command >&2
                Tempscript=`cat $tempscript`
                echo tempscript: $Tempscript >&2
	    fi
	    exec $torque_command $tempscript | awk '{print $1}' | awk -F. '{print $1}'
            rm $tempscript 
            # sleep 2
	else
            # echo "command=\`cat "$taskfile" | head -\$PBS_ARRAYID | tail -1\` ; exec \$command" > $tempscript
			echo "#!/bin/bash
			commandList=\`cat "$taskfile" | head -n \$PBS_ARRAYID | tail -n 1\` "$'\n'"
				
				commands=\$(grep -o ';' <<< \$commandList | wc -l)
				if [ \$commands -ne 0 ] ; then "$'\n'"
					command_exec=\`echo \$commandList | awk -F ';' '{print \$1}'\` "$'\n'"
					\$command_exec "$'\n'"
					commandList=\${commandList#*;} "$'\n'"
					temp=\"\${commandList//[^;]/}\" "$'\n'"
					length=\`echo \${#temp}\` "$'\n'"
					while [ \$length -gt 0 ] "$'\n'"
						do "$'\n'"
						command_exec=\`echo \$commandList | awk -F ';' '{print \$1}'\` "$'\n'"
						\$command_exec "$'\n'"
						commandList=\${commandList#*;} "$'\n'"
						temp=\"\${commandList//[^;]/}\" "$'\n'"
						length=\`echo \${#temp}\` "$'\n'"
					done "$'\n'"
					\$commandList "$'\n'"
				else "$'\n'"
					\$commandList "$'\n'"
				fi" > $tempscript # create script to excecute each line of the task file
            torque_command="qsub -V -d . -q $queue -M $mailto -N $JobName -m $MailOpts -l mem=$memoryLimit $LogOpts $sge_arch -l walltime=$WALLTIME$NumCPUs$torque_hold $sge_tasks"
			chmod 755 $tempscript
	    if [ $verbose -eq 1 ] ; then 
		echo torque_command: $torque_command >&2
		echo control file: $taskfile >&2
                Tempscript=`cat $tempscript`
                echo tempscript: $Tempscript >&2
	    fi
			exec $torque_command $tempscript | awk '{print $1}' | awk -F. '{print $1}' # formats the output to return the correct jobID
			rm $tempscript
            # sleep 2
	fi
	;;


###########################################################################
# Don't change the following - this runs the commands directly if a
# cluster is not being used.
###########################################################################

    NONE)
	if [ "x$tasks" = "x" ] ; then
	    if [ $verbose -eq 1 ] ; then 
		echo executing: $@ >&2
	    fi

	    /bin/sh <<EOF1 > ${LogDir}${JobName}.o$$ 2> ${LogDir}${JobName}.e$$
$@
EOF1
	else
	    if [ $verbose -eq 1 ] ; then 
		echo "Running commands in: $taskfile" >&2
	    fi

	    n=1
	    while [ $n -le $tasks ] ; do
		line=`sed -n -e ''${n}'p' $taskfile`
		if [ $verbose -eq 1 ] ; then 
		    echo executing: $line >&2
		fi
		/bin/sh <<EOF2 > ${LogDir}${JobName}.o$$.$n 2> ${LogDir}${JobName}.e$$.$n
$line
EOF2
		n=`expr $n + 1`
	    done
	fi	
	echo $$
	;;

esac

###########################################################################
# Done.
###########################################################################

