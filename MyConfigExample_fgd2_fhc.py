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
#Example: YieldsTable.py -s nuEleElastic,ccnue,otherBG -c SR1,VR1 -w results/MyLdmConfigExample/BkgOnly_combined_NormalMeasurement_model_afterFit.root -o test.tex
#            YieldsTable.py -s nuEleElastic,ccnue,otherBG -c SR1,VR1,VR2 -w results/MyLdmConfigExample/BkgOnly_combined_NormalMeasurement_model_afterFit.root -o test.tex
 
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
def nPOT(sam):
    pot = 0.
    for ffile in sam.files:
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
useStat=True
doValidation=True #use or use not validation regions to check exptrapolation to signal regions

#-------------------------------
# Parameters for hypothesis test
#-------------------------------
#configMgr.doHypoTest=False
#configMgr.nTOYs=1000
configMgr.calculatorType=2
configMgr.testStatType=3
configMgr.nPoints=20
#configMgr.ATLASLabelText = "Work in progress"

#--------------------------------
# Now we start to build the model
#--------------------------------
# First define HistFactory attributes
configMgr.analysisName = "MyLdm_fgd2_fhc"

# Scaling calculated by outputLumi / inputLumi
#configMgr.inputLumi = 1. # Luminosity of input TTree after weighting
#configMgr.outputLumi = 1. # Luminosity required for output histograms
#configMgr.setLumiUnits("fb-1")
configMgr.treeName = "default"

configMgr.histCacheFile = "data/"+configMgr.analysisName+".root"

configMgr.outputFileName = "results/"+configMgr.analysisName+"_Output.root"

#----------------------------------------------
# Set the files to read from
fhc_fgd1_mc_files = []
fhc_fgd1_data_files = []
if configMgr.readFromTree:
    fhc_fgd1_mc_files.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run2_anal_fgd2.root")
    fhc_fgd1_mc_files.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run4_anal_fgd2.root")

    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run2_air_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run2_water_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run3b_air_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run3c_air_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run4_air_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run4_water_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run8_air_fhc_fgd2.root")
    fhc_fgd1_data_files.append("samples/ldm/nd280Highland2_v2r45_1_data_run8_water_fhc_fgd2.root")

#----------------------------------------------
# List of samples and their plotting colours
topSample = Sample("nuEleElastic",kGreen-9)
topSample.setNormFactor("mu_Top",1.,0.,5.)
#topSample.setStatConfig(useStat)
#topSample.setNormRegions([("SR1","selelec_mom"),("VR1","selelec_mom")])
wzSample = Sample("ccnue",kAzure+1)
#wzSample.setNormFactor("mu_WZ",1.,0.,5.)
#wzSample.setStatConfig(useStat)
#wzSample.setNormRegions([("SR1","selelec_mom"),("VR1","selelec_mom")])
bgSample = Sample("otherBG",kYellow-3)
#bgSample.setNormFactor("mu_BG",1.,0.,5.)
#bgSample.setStatConfig(useStat)
#bgSample.setNormRegions([("SR1","selelec_mom"),("VR1","selelec_mom")])
#
dataSample = Sample("Data",kBlack)
#dataSample.cutsDict["SR1"] = "accum_level==14 && weight_corr_total>=0 "
dataSample.setData()


# Set the files to read from
bgdFiles = []
sigFiles = []
if configMgr.readFromTree:
    #bgdFiles.append("samples/ldm/production006B_mcp_genie_2010-11-water_magnet_run2_anal_fgd1.root")
    dataSample.setFileList(fhc_fgd1_data_files)
    pot_data = nPOT(dataSample)
    print(dataSample.name,pot_data,dataSample.weights)
    raw_input("Press Enter to continue.")
    # set the file from which the samples should be taken
    for sam in [topSample, wzSample, bgSample]:
            sam.setFileList(fhc_fgd1_mc_files)
            pot = nPOT(sam)
            #sam.addWeight(str(pot_data/pot))
            if sam==topSample:
                sam.weights = ("1","weight_corr_total",str(pot_data/pot),"(fgd2nuElectronElasticTopo==10)")
            elif sam==wzSample:
                sam.weights = ("1","weight_corr_total",str(pot_data/pot),"(fgd2nuElectronElasticTopo==1)")            
            else:
                sam.weights = ("1","weight_corr_total",str(pot_data/pot),"!(fgd2nuElectronElasticTopo==10 || fgd2nuElectronElasticTopo==1)")            
            print(sam.name,pot,sam.weights)
            raw_input("Press Enter to continue.")
else:
    workFile = ["data/"+configMgr.analysisName+".root"]
    for sam in [topSample, wzSample, bgSample, dataSample]:
            sam.setFileList(workFile)












# Dictionnary of cuts for Tree->hist
configMgr.cutsDict["SR1"] = "accum_level>14 && weight_corr_total>=0"
configMgr.cutsDict["VR1"] = "accum_level==14 && weight_corr_total>=0" 
configMgr.cutsDict["VR2"] = "accum_level==14 && selelec_costheta>0.9 && weight_corr_total>=0" 


# Tuples of nominal weights without and with b-jet selection
#configMgr.weights = ("1","0.333","weight_corr_total")


#--------------------
# List of systematics
#--------------------
# name of nominal histogram for systematics
configMgr.nomName = ""
#---------------------------------------------




#####################################
#************
#Bkg only fit
#************

bkt = configMgr.addFitConfig("BkgOnly")
if useStat:
    bkt.statErrThreshold=0.05 
else:
    bkt.statErrThreshold=None
bkt.addSamples([topSample,wzSample,bgSample,dataSample])

# Systematics to be applied globally within this topLevel
#bkt.getSample("Top").addSystematic(topKtScale)
#bkt.getSample("WZ").addSystematic(wzKtScale)

meas=bkt.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=0.039)
meas.addPOI("mu_SIG")
meas.addParamSetting("mu_BG",True,1)
meas.addParamSetting("Lumi",True,1)

#-------------------------------------------------
# Constraining regions - statistically independent
#-------------------------------------------------

emom = bkt.addChannel("selelec_mom",["SR1"],2,0.,1200.)
emom.useOverflowBin = True
#emom.addSystematic(jes)

bkt.setBkgConstrainChannels([emom])


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
emom.ATLASLabelText = "Work in progress"


#--------------------------------------------------------------
# Validation regions - not necessarily statistically independent
#--------------------------------------------------------------

if doValidation:
    emom1 = bkt.addChannel("selelec_mom",["VR1"],20,0.,2000.)
    emom1.useOverflowBin=True
    #emom1.addSystematic(jes)
    
    costheta1 = bkt.addChannel("selelec_costheta",["VR1"],10,0.,1.)
    costheta1.useOverflowBin=True
    #costheta1.addSystematic(jes)
    
    emom2 = bkt.addChannel("selelec_mom",["VR2"],20,0.,2000.)
    emom2.useOverflowBin=True

    costheta2 = bkt.addChannel("selelec_costheta",["VR2"],10,0.9,1.)
    costheta2.useOverflowBin=True

    bkt.setValidationChannels([emom1,costheta1,emom2,costheta2])
   
   


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


#-----------------------------
# Exclusion fits (1-step simplified model in this case)
#-----------------------------
if myFitType==FitType.Exclusion:
    sigSamples=["SM_GG_onestepCC_425_385_345"]
                        
    for sig in sigSamples:
        myTopLvl = configMgr.addFitConfigClone(bkt,"Sig_%s"%sig)

        sigSample = Sample(sig,kPink)
        sigSample.setFileList(sigFiles)
        sigSample.setNormByTheory()
        sigSample.setStatConfig(useStat)
        sigSample.setNormFactor("mu_SIG",1.,0.,5.)                    
        myTopLvl.addSamples(sigSample)
        myTopLvl.setSignalSample(sigSample)
    

        # s1l2j using met/meff
        if doValidation:
            mm2J = myTopLvl.getChannel("met/meff2Jet",["SS"])
            iPop=myTopLvl.validationChannels.index("SS_metmeff2Jet")
            myTopLvl.validationChannels.pop(iPop)
        else:
            mm2J = myTopLvl.addChannel("met/meff2Jet",["SS"],5,0.2,0.7)
            mm2J.useOverflowBin=True
            mm2J.addSystematic(jes)
            pass
        myTopLvl.setSignalChannels([mm2J])
	
	
	
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
    myTopLvl.tLegend = leg
c.Close()
