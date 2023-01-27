"""
 **********************************************************************************
 * Project: HistFitter - A ROOT-based package for statistical data analysis       *
 * Package: HistFitter                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      Simple example configuration with input trees                             *
 *                                                                                *
 * Authors:                                                                       *
 *      HistFitter group, CERN, Geneva                                            *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in the file          *
 * LICENSE.                                                                       *
 **********************************************************************************
"""

#Example: HistFitter.py  -t -w -D "before" -f -V -i ./MyConfigExample.py -L DEBUG
#Example: HistFitter.py -t  -w -D "before" -f -F excl -i ./MyConfigExample.py -L DEBUG
#         HistFitter.py -t  -w -D "before,after" -f -F excl -l -i ./MyConfigExample.py -L DEBUG

#Example: YieldsTable.py -s nuEleElastic,ccnue,otherBG -c SR1,VR1 -w results/MyLdmConfigExample/BkgOnly_combined_NormalMeasurement_model_afterFit.root -o test.tex
#          YieldsTable.py -s nuEleElastic,ccnue,otherBG -c SR1,VR1,VR2 -w results/MyLdmConfigExample/BkgOnly_combined_NormalMeasurement_model_afterFit.root -o test.tex
#          YieldsTable.py -s LDM020,nuEleElastic,ccnue,otherBG -c SR1_FGD1_FHC,SR1_FGD2_FHC,SR1_FGD1_RHC,SR1_FGD2_RHC -w results/MyLdmConfigExample/SPlusB_combined_NormalMeasurement_model_afterFit.root -o test.tex
#          YieldsTable.py -s LDM020,nuEleElastic,ccnue,otherBG -c SR1_FGD1_FHC,SR1_FGD2_FHC,SR1_FGD1_RHC,SR1_FGD2_RHC -w results/MyLdmConfigExample/SPlusB_combined_NormalMeasurement_model_afterFit.root -b -S -o test.tex
#
#     SysTable.py -s nuEleElastic,otherBG -c SR1_FGD1_FHC -w results/MyLdmConfigExample/SPlusB_combined_NormalMeasurement_model_afterFit.root -o tmp.tex
#     SysTable.py -s nuEleElastic,otherBG -c SR1_FGD1_FHC -w results/MyLdmConfigExample/SPlusB_combined_NormalMeasurement_model_afterFit.root -b -o tmp.tex
#     SysTable.py  -c SR1_FGD1_FHC,SR1_FGD2_FHC,SR1_FGD1_RHC,SR1_FGD2_RHC -w results/MyLdmConfigExample/SPlusB_combined_NormalMeasurement_model_afterFit.root -b -o tmp.tex

################################################################
## In principle all you have to setup is defined in this file ##
################################################################

## This configuration performs a simplified version of the "soft lepton" fits documented in ATLAS-CONF-2012-041.
## Only two systematics are considered:
##   -JES (Tree-based) conservatively treated like an MC stat error
##   -Alpgen Kt scale (weight-based)
##
## For the real complete implementation, see: HistFitterUser/MET_jets_leptons/python/MyOneLeptonKtScaleFit_mergerSoftLep.py

from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange,kDashed,kSolid,kDotted
#from ROOT import *
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt

from ROOT import gROOT, TLegend, TLegendEntry, TCanvas
#gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
#ROOT.SetAtlasStyle()

# function to calculate POTs for given sample
def nPOT(list_of_files):
    pot = 0.
    for ffile in list_of_files:
        pot += configMgr.cppMgr.GetPOT(ffile,True)
    return pot



#---------------------------------------------------------------------------------------------
# Some flags for overridding normal execution and telling ROOT to shut up... use with caution!
#---------------------------------------------------------------------------------------------
#gROOT.ProcessLine("gErrorIgnoreLevel=10001;")
#configMgr.plotHistos = True

#---------------------------------------
# Flags to control which fit is executed
#---------------------------------------
useStat=True #True
doValidation=False #use or use not validation regions to check exptrapolation to signal regions

#-------------------------------
# Parameters for hypothesis test
#-------------------------------
#configMgr.doHypoTest=False
#configMgr.nTOYs=1000
configMgr.calculatorType=2 #2   # 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=3     # 3=one-sided profile likelihood test statistic (LHC default)
configMgr.nPoints=20
#configMgr.ATLASLabelText = "Work in progress"

#--------------------------------
# Now we start to build the model
#--------------------------------
# First define HistFactory attributes
configMgr.analysisName = "MyLdmConfigExample"

# Scaling calculated by outputLumi / inputLumi
#configMgr.inputLumi = 1. # Luminosity of input TTree after weighting
#configMgr.outputLumi = 1. # Luminosity required for output histograms
#configMgr.setLumiUnits("fb-1")
configMgr.treeName = "default"

configMgr.histCacheFile = "data/"+configMgr.analysisName+".root"

configMgr.outputFileName = "results/"+configMgr.analysisName+"_Output.root"

#----------------------------------------------
# Set the files to read from
fhc_fgd1_sig_files = []
fhc_fgd1_mc_files = []
fhc_fgd1_data_files = []
fhc_fgd2_sig_files = []
fhc_fgd2_mc_files = []
fhc_fgd2_data_files = []
rhc_fgd1_sig_files = []
rhc_fgd1_mc_files = []
rhc_fgd1_data_files = []
rhc_fgd2_sig_files = []
rhc_fgd2_mc_files = []
rhc_fgd2_data_files = []
if configMgr.readFromTree:
    #FHC
    fhc_fgd1_mc_files.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run2_anal_fgd1.root")
    fhc_fgd1_mc_files.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run4_anal_fgd1.root")

    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run2_air_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run2_water_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run3b_air_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run3c_air_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run4_air_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run4_water_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run8_air_fhc_fgd1.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run8_water_fhc_fgd1.root")
    
    fhc_fgd1_sig_files.append("samples/ldm/fgd1_pi0_decay_prerun7geom.root")
    #
    fhc_fgd2_mc_files.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run2_anal_fgd2.root")
    fhc_fgd2_mc_files.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run4_anal_fgd2.root")

    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run2_air_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run2_water_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run3b_air_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run3c_air_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run4_air_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run4_water_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run8_air_fhc_fgd2.root")
    fhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run8_water_fhc_fgd2.root")
    
    fhc_fgd2_sig_files.append("samples/ldm/fgd2_pi0_decay_prerun7geom.root")
    #RHC
    rhc_fgd1_mc_files.append("samples/ldm/production006B_mcp_anti-genie_2010-11-air_magnet_run6_anal_fgd1.root")
    rhc_fgd1_mc_files.append("samples/ldm/production006B_mcp_anti-genie_2010-11-water_magnet_run5_anal_fgd1.root")

    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run5c_rhc_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6b_air_rhc_fgd1.roo")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6c_air_rhc_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6d_air_rhc_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6e_air_rhc_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run7b_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run9b_water_rhc_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run9c_water_rhc_fgd1.root")
    rhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run9d_water_rhc_fgd1.root")
    
    rhc_fgd1_sig_files.append("samples/ldm/fgd1_pi0_decay_prerun7geom.root")
    #
    rhc_fgd2_mc_files.append("samples/ldm/production006B_mcp_anti-genie_2010-11-air_magnet_run6_anal_fgd2.root")
    rhc_fgd2_mc_files.append("samples/ldm/production006B_mcp_anti-genie_2010-11-water_magnet_run5_anal_fgd2.root")

    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run5c_rhc_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6b_air_rhc_fgd2.roo")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6c_air_rhc_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6d_air_rhc_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run6e_air_rhc_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run7b_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run9b_water_rhc_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run9c_water_rhc_fgd2.root")
    rhc_fgd2_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run9d_water_rhc_fgd2.root")
    
    rhc_fgd2_sig_files.append("samples/ldm/fgd2_pi0_decay_prerun7geom.root")

#----------------------------------------------
#Systematics
ucb = Systematic("ucb", configMgr.weights, 1.3,0.7, "user","userOverallSys")
ucs = Systematic("ucs", configMgr.weights, 1.2,0.8, "user","userOverallSys")

#----------------------------------------------
# List of samples and their plotting colours
topSample = Sample("nuEleElastic",kGreen-9)
#topSample.setNormFactor("mu_Top",1.,0.,5.)
topSample.setStatConfig(useStat)
#topSample.setNormRegions([("SR1","selelec_mom"),("VR1","selelec_mom")])
topSample.addSystematic(ucb)
#
wzSample = Sample("ccnue",kAzure+1)
#wzSample.setNormFactor("mu_WZ",1.,0.,5.)
wzSample.setStatConfig(useStat)
#wzSample.setNormRegions([("SR1","selelec_mom"),("VR1","selelec_mom")])
wzSample.addSystematic(ucb)
#
bgSample = Sample("otherBG",kYellow-3)
#bgSample.setNormFactor("mu_BG",1.,0.,5.)
bgSample.setStatConfig(useStat)
#bgSample.setNormRegions([("SR1","selelec_mom"),("VR1","selelec_mom")])
bgSample.addSystematic(ucb)
#
dataSample = Sample("Data",kBlack)
dataSample.setData()
#
if myFitType==FitType.Exclusion:
    sigSamples=["LDM020"]
    for sig in sigSamples:
         sigSample = Sample(sig,kPink)
         #sigSample.setNormByTheory()
         sigSample.setStatConfig(useStat)
         sigSample.setNormFactor("mu_SIG",0.,0.,5.)                    
         sigSample.addSystematic(ucs)




# Set the files to read from
if configMgr.readFromTree:
    dataSample.setFileList_fgd1_fhc(fhc_fgd1_data_files)
    dataSample.setFileList_fgd2_fhc(fhc_fgd2_data_files)
    dataSample.setFileList_fgd1_rhc(rhc_fgd1_data_files)
    dataSample.setFileList_fgd2_rhc(rhc_fgd2_data_files)
    pot_data_fgd1_fhc = nPOT(dataSample.files_fgd1_fhc)
    pot_data_fgd2_fhc = nPOT(dataSample.files_fgd2_fhc)
    pot_data_fgd1_rhc = nPOT(dataSample.files_fgd1_rhc)
    pot_data_fgd2_rhc = nPOT(dataSample.files_fgd2_rhc)
    # set the file from which the samples should be taken
    for sam in [topSample, wzSample, bgSample]:
            sam.setFileList_fgd1_fhc(fhc_fgd1_mc_files)
            sam.setFileList_fgd2_fhc(fhc_fgd2_mc_files)
            sam.setFileList_fgd1_rhc(rhc_fgd1_mc_files)
            sam.setFileList_fgd2_rhc(rhc_fgd2_mc_files)
            pot_fgd1_fhc = nPOT(sam.files_fgd1_fhc)
            pot_fgd2_fhc = nPOT(sam.files_fgd2_fhc)
            pot_fgd1_rhc = nPOT(sam.files_fgd1_rhc)
            pot_fgd2_rhc = nPOT(sam.files_fgd2_rhc)
            if sam==topSample:
                sam.weights_fgd1_fhc = ("1","weight_corr_total",str(pot_data_fgd1_fhc/pot_fgd1_fhc),"(nuElectronElasticTopo==10)")
                sam.weights_fgd2_fhc = ("1","weight_corr_total",str(pot_data_fgd2_fhc/pot_fgd2_fhc),"(fgd2nuElectronElasticTopo==10)")
                sam.weights_fgd1_rhc = ("1","weight_corr_total",str(pot_data_fgd1_rhc/pot_fgd1_rhc),"(nuElectronElasticTopo==10)")
                sam.weights_fgd2_rhc = ("1","weight_corr_total",str(pot_data_fgd2_rhc/pot_fgd2_rhc),"(fgd2nuElectronElasticTopo==10)")
            elif sam==wzSample:
                sam.weights_fgd1_fhc = ("1","weight_corr_total",str(pot_data_fgd1_fhc/pot_fgd1_fhc),"(nuElectronElasticTopo==1)")            
                sam.weights_fgd2_fhc = ("1","weight_corr_total",str(pot_data_fgd2_fhc/pot_fgd2_fhc),"(fgd2nuElectronElasticTopo==1)")            
                sam.weights_fgd1_rhc = ("1","weight_corr_total",str(pot_data_fgd1_rhc/pot_fgd1_rhc),"(nuElectronElasticTopo==1)")            
                sam.weights_fgd2_rhc = ("1","weight_corr_total",str(pot_data_fgd2_rhc/pot_fgd2_rhc),"(fgd2nuElectronElasticTopo==1)")            
            else:
                sam.weights_fgd1_fhc = ("1","weight_corr_total",str(pot_data_fgd1_fhc/pot_fgd1_fhc),"!(nuElectronElasticTopo==10 || nuElectronElasticTopo==1)")            
                sam.weights_fgd2_fhc = ("1","weight_corr_total",str(pot_data_fgd2_fhc/pot_fgd2_fhc),"!(fgd2nuElectronElasticTopo==10 || fgd2nuElectronElasticTopo==1)")            
                sam.weights_fgd1_rhc = ("1","weight_corr_total",str(pot_data_fgd1_rhc/pot_fgd1_rhc),"!(nuElectronElasticTopo==10 || nuElectronElasticTopo==1)")            
                sam.weights_fgd2_rhc = ("1","weight_corr_total",str(pot_data_fgd2_rhc/pot_fgd2_rhc),"!(fgd2nuElectronElasticTopo==10 || fgd2nuElectronElasticTopo==1)")            
            print(dataSample.name,pot_data_fgd1_fhc,pot_data_fgd2_fhc,pot_data_fgd1_rhc,pot_data_fgd2_rhc)
            print(sam.name,pot_fgd1_fhc,pot_fgd2_fhc,sam.name,pot_fgd1_rhc,pot_fgd2_rhc)
            raw_input("Press Enter to continue.")
            
            if myFitType==FitType.Exclusion:
                for sig in sigSamples:
                    sigSample.setFileList_fgd1_fhc(fhc_fgd1_sig_files)
                    sigSample.setFileList_fgd2_fhc(fhc_fgd2_sig_files)
                    sigSample.setFileList_fgd1_rhc(rhc_fgd1_sig_files)
                    sigSample.setFileList_fgd2_rhc(rhc_fgd2_sig_files)
                    sigSample.weights_fgd1_fhc = ("1./3000.","weight_corr_total","ldm_expectedEvents","(ldm_Ma>1.99e-2 && ldm_Ma<2.01e-2)",str(pot_data_fgd1_fhc/1.8e+21))
                    sigSample.weights_fgd2_fhc = ("1./3000.","weight_corr_total","ldm_expectedEvents","(ldm_Ma>1.99e-2 && ldm_Ma<2.01e-2)",str(pot_data_fgd2_fhc/1.8e+21))
                    sigSample.weights_fgd1_rhc = ("1./3000.","weight_corr_total","ldm_expectedEvents","(ldm_Ma>1.99e-2 && ldm_Ma<2.01e-2)",str(pot_data_fgd1_rhc/1.8e+21))
                    sigSample.weights_fgd2_rhc = ("1./3000.","weight_corr_total","ldm_expectedEvents","(ldm_Ma>1.99e-2 && ldm_Ma<2.01e-2)",str(pot_data_fgd2_rhc/1.8e+21))

else:
    workFile = ["data/"+configMgr.analysisName+".root"]
    if myFitType==FitType.Exclusion:
        for sam in [topSample, wzSample, bgSample, dataSample,sigSample]:
            sam.setFileList(workFile)
    else:
        for sam in [topSample, wzSample, bgSample, dataSample]:
            sam.setFileList(workFile)


#--------------------
# Dictionnary of cuts for Tree->hist
#--------------------
configMgr.cutsDict["SR1_FGD1_FHC"] = "accum_level>14 && weight_corr_total>=0"
configMgr.cutsDict["VR1_FGD1_FHC"] = "accum_level==14 && weight_corr_total>=0" 
configMgr.cutsDict["VR2_FGD1_FHC"] = "accum_level==14 && selelec_costheta>0.9 && weight_corr_total>=0" 
#
configMgr.cutsDict["SR1_FGD2_FHC"] = "accum_level>14 && weight_corr_total>=0"
configMgr.cutsDict["VR1_FGD2_FHC"] = "accum_level==14 && weight_corr_total>=0" 
configMgr.cutsDict["VR2_FGD2_FHC"] = "accum_level==14 && selelec_costheta>0.9 && weight_corr_total>=0" 
#RHC
configMgr.cutsDict["SR1_FGD1_RHC"] = "accum_level>14 && weight_corr_total>=0"
configMgr.cutsDict["VR1_FGD1_RHC"] = "accum_level==14 && weight_corr_total>=0" 
configMgr.cutsDict["VR2_FGD1_RHC"] = "accum_level==14 && selelec_costheta>0.9 && weight_corr_total>=0" 
#
configMgr.cutsDict["SR1_FGD2_RHC"] = "accum_level>14 && weight_corr_total>=0"
configMgr.cutsDict["VR1_FGD2_RHC"] = "accum_level==14 && weight_corr_total>=0" 
configMgr.cutsDict["VR2_FGD2_RHC"] = "accum_level==14 && selelec_costheta>0.9 && weight_corr_total>=0" 

#--------------------
# List of systematics
#--------------------
# name of nominal histogram for systematics
configMgr.nomName = ""
#---------------------------------------------




#############################################################
#************
#Bkg only fit
#************
if myFitType==FitType.Background:
    bkt = configMgr.addFitConfig("BkgOnly")
    if useStat:
        bkt.statErrThreshold=0.005 
    else:
        bkt.statErrThreshold=None
    bkt.addSamples([topSample,wzSample,bgSample,dataSample])
    
#************
#Exclusion fit
#************
elif myFitType==FitType.Exclusion:
    bkt = configMgr.addFitConfig("SPlusB")
    if useStat:
        bkt.statErrThreshold=0.005 
    else:
        bkt.statErrThreshold=None
    bkt.addSamples([topSample,wzSample,bgSample,dataSample,sigSample])
    bkt.setSignalSample(sigSample)

    # Systematics to be applied globally within this topLevel
    #bkt.getSample("Top").addSystematic(topKtScale)
    #bkt.getSample("WZ").addSystematic(wzKtScale)


#-----------------------------
#Measurement
#-----------------------------
meas=bkt.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=0.0001)
meas.addPOI("mu_SIG")
#meas.addParamSetting("mu_Top",True,1)
#meas.addParamSetting("mu_BG",True,1)
#meas.addParamSetting("Lumi",True,1)

    
#-------------------------------------------------
# Constraining regions - statistically independent
#-------------------------------------------------
Nbins = 1
mm2J_fgd1_fhc = bkt.addChannel("selelec_mom",["SR1_FGD1_FHC"],Nbins,0.,1200.)
mm2J_fgd1_fhc.useOverflowBin=True

mm2J_fgd2_fhc = bkt.addChannel("selelec_mom",["SR1_FGD2_FHC"],Nbins,0.,1200.)
mm2J_fgd2_fhc.useOverflowBin=True

mm2J_fgd1_rhc = bkt.addChannel("selelec_mom",["SR1_FGD1_RHC"],Nbins,0.,1200.)
mm2J_fgd1_rhc.useOverflowBin=True

mm2J_fgd2_rhc = bkt.addChannel("selelec_mom",["SR1_FGD2_RHC"],Nbins,0.,1200.)
mm2J_fgd2_rhc.useOverflowBin=True

#--------------------------------------------------------------
# Validation regions - not necessarily statistically independent
#--------------------------------------------------------------

if doValidation:
    #emom1 = bkt.addChannel("selelec_mom",["SR1"],2,0.,1200.)
    #emom1.useOverflowBin=True
    #emom1.addSystematic(jes)
    
    costheta1 = bkt.addChannel("selelec_costheta",["VR1_FGD1_FHC"],10,0.,1.)
    costheta1.useOverflowBin=True
    #costheta1.addSystematic(jes)
    
    emom2 = bkt.addChannel("selelec_mom",["VR2_FGD1_FHC"],20,0.,2000.)
    emom2.useOverflowBin=True

    costheta2 = bkt.addChannel("selelec_costheta",["VR2_FGD1_FHC"],10,0.9,1.)
    costheta2.useOverflowBin=True

    #bkt.setValidationChannels([costheta1,emom2,costheta2])
   

if myFitType==FitType.Background:
     bkt.setBkgConstrainChannels([mm2J_fgd1_fhc,mm2J_fgd2_fhc,mm2J_fgd1_rhc,mm2J_fgd2_rhc])
     if doValidation:
         bkt.setValidationChannels([costheta1,emom2,costheta2])
elif myFitType==FitType.Exclusion:
     bkt.setSignalChannels([mm2J_fgd1_fhc,mm2J_fgd2_fhc,mm2J_fgd1_rhc,mm2J_fgd2_rhc])
     #bkt.setValidationChannels([mm2J_fgd2_fhc])
     if doValidation:
         bkt.setValidationChannels([costheta1,emom2,costheta2])

###################
#                                               #
#    Example new cosmetics     #
#                                               #
###################

# Set global plotting colors/styles
bkt.dataColor = dataSample.color
bkt.totalPdfColor = kBlue
bkt.errorFillColor = kBlue-5
bkt.errorFillStyle = 3004
bkt.errorLineStyle = kDashed
bkt.errorLineColor = kBlue-5

# Set Channel titleX, titleY, minY, maxY, logY
#emom.minY = 0.5
#emom.maxY = 5000
#emom.titleX = "n jets"
#emom.titleY = "Entries"
#emom.logY = True
#emom.ATLASLabelX = 0.25
#emom.ATLASLabelY = 0.85
#emom.ATLASLabelText = "Work in progress"


   


#**************
# Discovery fit
#**************

if myFitType==FitType.Discovery:
    discovery = configMgr.addFitConfigClone(bkt,"Discovery")
    
    # s1l2jT = signal region/channel
    #ssChannel = discovery.addChannel("cuts",["SS"],srNBins,srBinLow,srBinHigh)
    #ssChannel.addSystematic(jes)
    #ssChannel.addDiscoverySamples(["SS"],[1.],[0.],[100.],[kMagenta])
    #discovery.setSignalChannels([ssChannel])


	
	
	
# Create TLegend (AK: TCanvas is needed for that, but it gets deleted afterwards)
c = TCanvas()
compFillStyle = 1001 # see ROOT for Fill styles
leg = TLegend(0.6,0.475,0.9,0.925,"")
leg.SetFillStyle(0)
leg.SetFillColor(0)
leg.SetBorderSize(0)
#
entry = TLegendEntry()
entry = leg.AddEntry("","DATA","pl") 
entry.SetMarkerColor(bkt.dataColor)
entry.SetMarkerStyle(20)
#
entry = leg.AddEntry("","Total pdf","lf") 
entry.SetLineColor(bkt.totalPdfColor)
entry.SetLineWidth(2)
entry.SetFillColor(bkt.errorFillColor)
entry.SetFillStyle(bkt.errorFillStyle)
#
entry = leg.AddEntry("",topSample.name,"lf") 
entry.SetLineColor(topSample.color)
entry.SetFillColor(topSample.color)
entry.SetFillStyle(compFillStyle)
#
entry = leg.AddEntry("",wzSample.name,"lf") 
entry.SetLineColor(wzSample.color)
entry.SetFillColor(wzSample.color)
entry.SetFillStyle(compFillStyle)
#
#entry = leg.AddEntry("","multijets (data estimate)","lf") 
#entry.SetLineColor(qcdSample.color)
#entry.SetFillColor(qcdSample.color)
#entry.SetFillStyle(compFillStyle)
#
entry = leg.AddEntry("",bgSample.name,"lf") 
entry.SetLineColor(bgSample.color)
entry.SetFillColor(bgSample.color)
entry.SetFillStyle(compFillStyle)
#
if myFitType==FitType.Exclusion:
    entry = leg.AddEntry("","signal","lf") 
    entry.SetLineColor(kPink)
    entry.SetFillColor(kPink)
    entry.SetFillStyle(compFillStyle)


# Set legend for fitConfig
bkt.tLegend = leg
if myFitType==FitType.Exclusion:
    bkt.tLegend = leg
c.Close()
