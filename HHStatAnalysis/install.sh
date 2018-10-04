#!/bin/bash
# Install cms-hh framework.
# This file is part of https://github.com/cms-hh/HHStatAnalysis.

INSTALL_MODES=(full plotting)
DEFAULT_MODE=full
DEFAULT_N_JOBS=8
DEFAULT_RELEASE="CMSSW_7_4_7"

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

if [ $# -gt 3 -o "x$1" = "x-h" -o "x$1" = "x--help" ] ; then
    echo "Usage: [mode] [n_jobs] [cmssw_release]"
    printf "\n\tmode\t\t\tinstallation mode. Supported modes: "
    join_by ", " "${INSTALL_MODES[@]}"
    printf ". Default: $DEFAULT_MODE."
    printf "\n\tn_jobs\t\t\tthe number of jobs to run simultaneous during the compilation. Default: $DEFAULT_N_JOBS.\n"
    printf "\tcmssw_release\t\tCMSSW release."
    printf " Default: $DEFAULT_RELEASE.\n"
    exit 1
fi

MODE=$1
if [ "x$MODE" = "x" ] ; then
    MODE=$DEFAULT_MODE
elif [[ ! ${INSTALL_MODES[@]} =~ $MODE ]] ; then
    echo "ERROR: unsupported installation mode '$MODE'."
    echo "Supported installation modes: ${INSTALL_MODES[@]}"
    exit 1
fi

N_JOBS=$2
if [ "x$N_JOBS" = "x" ] ; then N_JOBS=$DEFAULT_N_JOBS ; fi
if ! [ $N_JOBS -eq $N_JOBS ] 2>/dev/null ; then
    echo "ERROR: invalid number of jobs '$N_JOBS'."
    exit 1
fi

RELEASE=$3
if [ "x$RELEASE" = "x" ] ; then
    RELEASE=$DEFAULT_RELEASE
fi
if [ -e $RELEASE ] ; then
    echo "ERROR: Working area for $RELEASE already exists."
    exit 1
fi

export SCRAM_ARCH=slc6_amd64_gcc491
scramv1 project CMSSW $RELEASE
RESULT=$?
if [ $RESULT -ne 0 ] ; then
    echo "ERROR: unable to create working area for CMSSW release '$RELEASE'."
    exit 2
fi

cd $RELEASE/src
eval `scramv1 runtime -sh`
RESULT=$?
if [ $RESULT -ne 0 ] ; then
    echo "ERROR: unable to setup the environment for CMSSW release '$RELEASE'."
    exit 2
fi

if [ $MODE = "full" ] ; then
    # Combine tool
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git checkout v6.3.0
    cd ../..

    # CombineHarvester package
    mkdir CombineHarvester
    cd CombineHarvester
    git init
    git remote add origin git@github.com:cms-analysis/CombineHarvester.git
    git config core.sparsecheckout true
    echo -e "CombineTools/\nCombinePdfs/" >> .git/info/sparse-checkout
    git pull origin master
    cd ..

    # HHStatAnalysis package
    git clone -o cms-hh git@github.com:cms-hh/HHStatAnalysis.git
else
    # HHStatAnalysis package
    mkdir HHStatAnalysis
    cd HHStatAnalysis
    git init
    git remote add cms-hh git@github.com:cms-hh/HHStatAnalysis.git
    git config core.sparsecheckout true
    echo -e "Core/" >> .git/info/sparse-checkout
    git pull cms-hh master
    cd ..
fi

git clone -o cms-hh git@github.com:cms-hh/Plotting.git HHStatAnalysis/Plotting
git clone -o cms-hh git@github.com:cms-hh/Resources.git HHStatAnalysis/Resources

GITHUB_USER=$(git config user.github)
if ! [ "x$GITHUB_USER" = "x" ] ; then
    cd HHStatAnalysis
    git ls-remote git@github.com:$GITHUB_USER/HHStatAnalysis.git &> /dev/null
    RESULT=$?
    if [ $RESULT -eq 0 ] ; then
        git remote add $GITHUB_USER git@github.com:$GITHUB_USER/HHStatAnalysis.git
        git fetch $GITHUB_USER
    fi
    cd ..

    cd HHStatAnalysis/Plotting
    git ls-remote git@github.com:$GITHUB_USER/Plotting.git &> /dev/null
    RESULT=$?
    if [ $RESULT -eq 0 ] ; then
        git remote add $GITHUB_USER git@github.com:$GITHUB_USER/Plotting.git
        git fetch $GITHUB_USER
    fi
    cd ../..

    cd HHStatAnalysis/Resources
    git ls-remote git@github.com:$GITHUB_USER/Resources.git &> /dev/null
    RESULT=$?
    if [ $RESULT -eq 0 ] ; then
        git remote add $GITHUB_USER git@github.com:$GITHUB_USER/Resources.git
        git fetch $GITHUB_USER
    fi
    cd ../..
fi

scram b -j$N_JOBS
