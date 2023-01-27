# Setup the PATH and LD_LIBRARY_PATH environment variables 
# and run all setup scripts
HERE=`pwd` 
export ND280_SYSTEM=`nd280-system` 
unset -f path_remove 
unset -f path_append 
unset -f ld_remove 
unset -f ld_append 
path_remove ()  { export PATH=`echo -n $PATH | awk -v RS=: -v ORS=: '$0 != "'$1'"' | sed 's/:$//'`; } 
path_append ()  { path_remove $1; export PATH="$PATH:$1" ;} 
ld_remove ()  { export LD_LIBRARY_PATH=`echo -n $LD_LIBRARY_PATH | awk -v RS=: -v ORS=: '$0 != "'$1'"' | sed 's/:$//'`; } 
ld_append ()  { ld_remove $1; export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1" ;} 
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/LDMAnalysis/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/LDMAnalysis/${ND280_SYSTEM}/lib 
export ND280SOFTWAREPOLICYROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/nd280SoftwarePolicy_master/modules/.. 
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/MYSQL_5.6.20.01/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/MYSQL_5.6.20.01/${ND280_SYSTEM}/lib 
export MYSQLROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/MYSQL_5.6.20.01 
export MYSQLCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/GSL_1.15.0.00/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/GSL_1.15.0.00/${ND280_SYSTEM}/lib 
export GSLROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/GSL_1.15.0.00 
export GSLCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/ROOT_5.34.34.00/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/ROOT_5.34.34.00/${ND280_SYSTEM}/lib 
export ROOTROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/ROOT_5.34.34.00 
export ROOTCONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/ROOT_5.34.34.00/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheCore_3.45/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheCore_3.45/${ND280_SYSTEM}/lib 
export PSYCHECOREROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheCore_3.45 
export PSYCHECORECONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheCore_3.45/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheEventModel_3.41/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheEventModel_3.41/${ND280_SYSTEM}/lib 
export PSYCHEEVENTMODELROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheEventModel_3.41 
export PSYCHEEVENTMODELCONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheEventModel_3.41/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheUtils_3.34/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheUtils_3.34/${ND280_SYSTEM}/lib 
export PSYCHEUTILSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheUtils_3.34 
export PSYCHEUTILSCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheND280Utils_3.63.1/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheND280Utils_3.63.1/${ND280_SYSTEM}/lib 
export PSYCHEND280UTILSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheND280Utils_3.63.1 
export PSYCHEND280UTILSCONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheND280Utils_3.63.1/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheIO_3.34/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheIO_3.34/${ND280_SYSTEM}/lib 
export PSYCHEIOROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheIO_3.34 
export PSYCHEIOCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCore_2.40/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCore_2.40/${ND280_SYSTEM}/lib 
export HIGHLANDCOREROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCore_2.40 
export HIGHLANDCORECONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCore_2.40/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/oaAnalysisReader_2.23/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/oaAnalysisReader_2.23/${ND280_SYSTEM}/lib 
export OAANALYSISREADERROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/oaAnalysisReader_2.23 
export OAANALYSISREADERCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandEventModel_2.34/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandEventModel_2.34/${ND280_SYSTEM}/lib 
export HIGHLANDEVENTMODELROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandEventModel_2.34 
export HIGHLANDEVENTMODELCONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandEventModel_2.34/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandUtils_2.39/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandUtils_2.39/${ND280_SYSTEM}/lib 
export HIGHLANDUTILSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandUtils_2.39 
export HIGHLANDUTILSCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheSelections_3.51/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheSelections_3.51/${ND280_SYSTEM}/lib 
export PSYCHESELECTIONSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheSelections_3.51 
export PSYCHESELECTIONSCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheSystematics_3.54/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheSystematics_3.54/${ND280_SYSTEM}/lib 
export PSYCHESYSTEMATICSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/psycheSystematics_3.54 
export PSYCHESYSTEMATICSCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCorrections_2.25/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCorrections_2.25/${ND280_SYSTEM}/lib 
export HIGHLANDCORRECTIONSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandCorrections_2.25 
export HIGHLANDCORRECTIONSCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandIO_2.45/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandIO_2.45/${ND280_SYSTEM}/lib 
export HIGHLANDIOROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandIO_2.45 
export HIGHLANDIOCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandTools_2.29/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandTools_2.29/${ND280_SYSTEM}/lib 
export HIGHLANDTOOLSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandTools_2.29 
export HIGHLANDTOOLSCONFIG=${ND280_SYSTEM}  
for file in /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandTools_2.29/setup_script/*.sh ; do [ -f $file ] && . $file ; done
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseAnalysis_2.33/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseAnalysis_2.33/${ND280_SYSTEM}/lib 
export BASEANALYSISROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseAnalysis_2.33 
export BASEANALYSISCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseTrackerAnalysis_2.25/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseTrackerAnalysis_2.25/${ND280_SYSTEM}/lib 
export BASETRACKERANALYSISROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseTrackerAnalysis_2.25 
export BASETRACKERANALYSISCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandSystematics_2.10/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandSystematics_2.10/${ND280_SYSTEM}/lib 
export HIGHLANDSYSTEMATICSROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/highlandSystematics_2.10 
export HIGHLANDSYSTEMATICSCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/nueCCAnalysis_2.29/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/nueCCAnalysis_2.29/${ND280_SYSTEM}/lib 
export NUECCANALYSISROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/nueCCAnalysis_2.29 
export NUECCANALYSISCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseP0DAnalysis_2.7/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseP0DAnalysis_2.7/${ND280_SYSTEM}/lib 
export BASEP0DANALYSISROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/baseP0DAnalysis_2.7 
export BASEP0DANALYSISCONFIG=${ND280_SYSTEM}  
path_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/LDMAnalysis/${ND280_SYSTEM}/bin 
ld_append /meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/LDMAnalysis/${ND280_SYSTEM}/lib 
export LDMANALYSISROOT=/meg/scratch/suslov/work_t2k/T2K/work/highland2Software_ldm0/LDMAnalysis 
export LDMANALYSISCONFIG=${ND280_SYSTEM}  
export LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | sed 's/^://g'`
cd $HERE 
