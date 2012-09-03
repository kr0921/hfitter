#!/usr/bin/python
# -*- coding: utf-8 -*-
from ROOT import TFile, TMath, RooRandom, TH1, TH1F
from ROOT import kBlack, kWhite, kGray, kRed, kPink, kMagenta, kViolet, \
    kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kOrange, \
    kDashed, kSolid, kDotted
from os import system
from math import fabs

import generateToys

TH1.SetDefaultSumw2(True)

from copy import deepcopy, copy
from configManager import configMgr


class Sample(object):

    """
    Defines a Sample in a channel XML file
    """

    def __init__(self, name, color=1.):
        """
        Store configuration, set sample name, and if to normalize by theory

        Scales histograms to luminosity set in configuration
        """

        self.name = name
        self.color = color
        self.isData = False
        self.isQCD = False
        self.isDiscovery = False
        self.write = True
        self.normByTheory = False
        self.statConfig = False
        self.histoSystList = []
        self.shapeSystList = []
        self.overallSystList = []
        self.shapeFactorList = []
        self.systList = []
        self.weights = []
        self.systDict = {}
        self.normFactor = []
        self.qcdSyst = None
        self.unit = 'GeV'
        self.cutsDict = {}
        self.files = []
        self.treeName = ''
        self.xsecWeight = None
        self.xsecUp = None
        self.xsecDown = None

    def buildHisto(
        self,
        binValues,
        region,
        var,
        ):
        """
        Allow user to give bin values eg. for checking stats in papers
        """

        try:
            self.binValues[(region, var)] = binValues
        except AttributeError:
            self.binValues = {}
            self.binValues[(region, var)] = binValues
        if not self.isData:
            self.histoName = 'h' + self.name + 'Nom_' + region + '_obs_'\
                 + var
        else:
            self.histoName = 'h' + self.name + '_' + region + '_obs_'\
                 + var
        configMgr.hists[self.histoName] = TH1F(self.histoName,
                self.histoName, len(self.binValues[(region, var)]), 0.,
                float(len(self.binValues[(region, var)])))
        for (iBin, val) in enumerate(self.binValues[(region, var)]):
            configMgr.hists[self.histoName].SetBinContent(iBin + 1.,
                    val)

    def buildStatErrors(
        self,
        binStatErrors,
        region,
        var,
        ):
        """
        Allow user to give bin stat errors eg. for checking stats in papers
        """

        try:
            self.binStatErrors[(region, var)] = binStatErrors
        except AttributeError:
            self.binStatErrors = {}
            self.binStatErrors[(region, var)] = binStatErrors
        if not len(self.binStatErrors[(region, var)])\
             == len(self.binValues[(region, var)]):
            raise Exception('Length of errors list in region %s and variable %s does not match the nominal histogram!'
                             % (region, var))
        if not self.isData:
            self.histoName = 'h' + self.name + 'Nom_' + region + '_obs_'\
                 + var
        else:
            self.histoName = 'h' + self.name + '_' + region + '_obs_'\
                 + var
        for (iBin, err) in enumerate(self.binStatErrors[(region, var)]):
            try:
                configMgr.hists[self.histoName].SetBinError(iBin + 1.,
                        err)
            except:
                raise Exception('Errors specified without building histogram!'
                                )

    def Clone(self):
        newInst = deepcopy(self)

        # for (key,val) in self.systDict.items():
        #    newInst.systDict[key] = val

        return newInst

    def setUnit(self, unit):
        self.unit = unit
        return

    def setCutsDict(self, cutsDict):
        self.cutsDict = cutsDict
        return

    def setData(self, isData=True):
        self.isData = isData
        return

    def setWeights(self, weights):
        """
        Set the weights for this sample - overrides
        """

        self.weights = deepcopy(weights)
        return

    def addWeight(self, weight):
        """
        Add a weight to this sample and propagate
        """

        if not weight in self.weights:
            self.weights.append(weight)
        else:
            raise RuntimeError('Weight %s already defined in sample %s'
                                % (weight, self.name))
        for syst in self.systDict.values():
            if syst.type == 'weight':
                if not weight in syst.high:
                    syst.high.append(weight)
                if not weight in syst.low:
                    syst.low.append(weight)
        return

    def removeWeight(self, weight):
        if weight in self.weights:
            self.weights.remove(weight)
        for syst in self.systDict.values():
            if syst.type == 'weight':
                if weight in syst.high:
                    syst.high.remove(weight)
                if weight in syst.low:
                    syst.low.remove(weight)
        return

    def setQCD(self, isQCD=True, qcdSyst='uncorr'):
        self.isQCD = isQCD
        self.qcdSyst = qcdSyst
        return

    def setDiscovery(self, isDiscovery=True):
        self.isDiscovery = isDiscovery
        return

    def setNormByTheory(self, normByTheory=True):
        self.normByTheory = normByTheory
        return

    def setStatConfig(self, statConfig):
        self.statConfig = statConfig
        return

    def setWrite(self, write=True):
        self.write = write
        return

    def setHistoName(self, histoName):
        """
        Set the name of the nominal histogram for this sample
        """

        self.histoName = histoName
        return

    def setTreeName(self, treeName):
        self.treeName = treeName
        return

    def propagateTreeName(self, treeName):
        if self.treeName == '':
            self.treeName = treeName

        # ## MAB: Propagate treeName down to systematics of sample
        # for (systName,systList) in self.systDict.items():
        #    for syst in systList:
        #        syst.propagateTreeName(self.treeName)
        #        pass

        return

    def addHistoSys(
        self,
        systName,
        nomName,
        highName,
        lowName,
        includeOverallSys,
        normalizeSys,
        symmetrize=True,
        oneSide=False,
        samName='',
        normString='',
        ):
        """
        Add a HistoSys entry using the nominal, high and low histograms, set if to include OverallSys

        If includeOverallSys then extract scale factors

        If normalizeSys then normalize shapes to nominal
        """

        if normalizeSys:
            highIntegral = configMgr.hists['h' + samName + systName
                     + 'High_' + normString + 'Norm'].Integral()
            lowIntegral = configMgr.hists['h' + samName + systName
                     + 'Low_' + normString + 'Norm'].Integral()
            nomIntegral = configMgr.hists['h' + samName + 'Nom_'
                     + normString + 'Norm'].Integral()

            try:
                high = highIntegral / nomIntegral
                low = lowIntegral / nomIntegral
            except ZeroDivisionError:
                print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g remove from fit.'\
                     % (nomName, systName, nomIntegral, highIntegral,
                        lowIntegral)
                return

            configMgr.hists[highName + 'Norm'] = \
                configMgr.hists[highName].Clone(highName + 'Norm')
            configMgr.hists[lowName + 'Norm'] = \
                configMgr.hists[lowName].Clone(lowName + 'Norm')

            try:
                configMgr.hists[highName + 'Norm'].Scale(1. / high)
                configMgr.hists[lowName + 'Norm'].Scale(1. / low)
            except ZeroDivisionError:

                print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g remove from fit.'\
                     % (nomName, systName, nomIntegral, highIntegral,
                        lowIntegral)
                return
            if oneSide:
                if configMgr.hists[highName].Integral()\
                     > configMgr.hists[nomName].Integral():
                    self.histoSystList.append((
                        systName,
                        highName + 'Norm',
                        nomName,
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
                else:
                    self.histoSystList.append((
                        systName,
                        nomName,
                        lowName + 'Norm',
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
            else:
                self.histoSystList.append((
                    systName,
                    highName + 'Norm',
                    lowName + 'Norm',
                    configMgr.histCacheFile,
                    '',
                    '',
                    '',
                    '',
                    ))
        elif includeOverallSys:
            nomIntegral = configMgr.hists[nomName].Integral()
            lowIntegral = configMgr.hists[lowName].Integral()
            highIntegral = configMgr.hists[highName].Integral()
            try:
                high = highIntegral / nomIntegral
                low = lowIntegral / nomIntegral
            except ZeroDivisionError:
                print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g remove from fit.'\
                     % (nomName, systName, nomIntegral, highIntegral,
                        lowIntegral)
                return

            configMgr.hists[highName + 'Norm'] = \
                configMgr.hists[highName].Clone(highName + 'Norm')
            configMgr.hists[lowName + 'Norm'] = \
                configMgr.hists[lowName].Clone(lowName + 'Norm')
            try:
                configMgr.hists[highName + 'Norm'].Scale(1. / high)
                configMgr.hists[lowName + 'Norm'].Scale(1. / low)
            except ZeroDivisionError:
                print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g keeping in fit (offending histogram should be empty).'\
                     % (nomName, systName, nomIntegral, highIntegral,
                        lowIntegral)
            if high < 1. and low < 1.:
                print '    WARNING addHistoSys: high=%f is < 1.0 in %s. Taking symmetric value from low %f %f'\
                     % (high, systName, low, 2. - low)
                high = 2. - low
            if low > 1. and high > 1.:
                print '    WARNING addHistoSys: low=%f is > 1.0 in %s. Taking symmetric value from high %f %f'\
                     % (low, systName, high, 2. - high)
                low = 2. - high
            if low < 0.:
                print '    WARNING addHistoSys: low=%f < 0.0 in %s. Setting low=0.0.'\
                     % (low, systName)
                low = 0.

            self.histoSystList.append((
                systName,
                highName + 'Norm',
                lowName + 'Norm',
                configMgr.histCacheFile,
                '',
                '',
                '',
                '',
                ))
            self.addOverallSys(systName, high, low)
        else:
            if symmetrize:
                nomIntegral = configMgr.hists[nomName].Integral()
                lowIntegral = configMgr.hists[lowName].Integral()
                highIntegral = configMgr.hists[highName].Integral()
                try:
                    high = highIntegral / nomIntegral
                    low = lowIntegral / nomIntegral
                except ZeroDivisionError:
                    print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g remove from fit.'\
                         % (nomName, systName, nomIntegral,
                            highIntegral, lowIntegral)
                    return
                if high < 1. and low < 1.:
                    print '    WARNING addHistoSys: high=%f is < 1.0 in %s. Taking symmetric value from low %f %f'\
                         % (high, systName, low, 2. - low)
                    configMgr.hists[highName + 'Norm'] = \
                        configMgr.hists[highName].Clone(highName
                             + 'Norm')
                    try:
                        configMgr.hists[highName + 'Norm'].Scale((2.
                                 - low) / high)
                    except ZeroDivisionError:
                        print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g remove from fit.'\
                             % (nomName, systName, nomIntegral,
                                highIntegral, lowIntegral)
                        return
                    self.histoSystList.append((
                        systName,
                        highName + 'Norm',
                        lowName,
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
                    if not systName in configMgr.systDict.keys():
                        self.systList.append(systName)
                    return
                if low > 1. and high > 1.:
                    print '    WARNING addHistoSys: low=%f is > 1.0 in %s. Taking symmetric value from high %f %f'\
                         % (low, systName, high, 2. - high)
                    configMgr.hists[lowName + 'Norm'] = \
                        configMgr.hists[lowName].Clone(lowName + 'Norm')
                    try:
                        configMgr.hists[lowName + 'Norm'].Scale((2.
                                 - high) / low)
                    except ZeroDivisionError:
                        print '    ERROR: generating HistoSys for %s syst=%s nom=%g high=%g low=%g remove from fit.'\
                             % (nomName, systName, nomIntegral,
                                highIntegral, lowIntegral)
                        return
                    self.histoSystList.append((
                        systName,
                        highName,
                        lowName + 'Norm',
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
                    if not systName in configMgr.systDict.keys():
                        self.systList.append(systName)
                    return
                if low < 0.:
                    print '    WARNING addHistoSys: low=%f is < 0.0 in %s. Setting negative bins to 0.0.'\
                         % (low, systName)
                    configMgr.hists[lowName + 'Norm'] = \
                        configMgr.hists[lowName].Clone(lowName + 'Norm')
                    for iBin in xrange(1., configMgr.hists[lowName
                             + 'Norm'].GetNbinsX() + 1.):
                        if configMgr.hists[lowName + 'Norm'
                                ].GetBinContent(iBin) < 0.:
                            configMgr.hists[lowName + 'Norm'
                                    ].SetBinContent(iBin, 0.)
                    self.histoSystList.append((
                        systName,
                        highName,
                        lowName + 'Norm',
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
                    if not systName in configMgr.systDict.keys():
                        self.systList.append(systName)
                        return
                self.histoSystList.append((
                    systName,
                    highName,
                    lowName,
                    configMgr.histCacheFile,
                    '',
                    '',
                    '',
                    '',
                    ))
                return
            elif oneSide:
                if configMgr.hists[highName].Integral()\
                     > configMgr.hists[nomName].Integral():
                    self.histoSystList.append((
                        systName,
                        highName,
                        nomName,
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
                else:

                # elif configMgr.hists[lowName].Integral() < configMgr.hists[nomName].Integral():

                    self.histoSystList.append((
                        systName,
                        nomName,
                        lowName,
                        configMgr.histCacheFile,
                        '',
                        '',
                        '',
                        '',
                        ))
                return
            else:
                self.histoSystList.append((
                    systName,
                    highName,
                    lowName,
                    configMgr.histCacheFile,
                    '',
                    '',
                    '',
                    '',
                    ))
                return
        if not systName in configMgr.systDict.keys():
            self.systList.append(systName)

    def addShapeSys(
        self,
        systName,
        nomName,
        highName,
        lowName,
        constraintType='Gaussian',
        ):
        """
        Add a ShapeSys entry using the nominal, high and low histograms
        """

        overallSystHigh = 1.
        overallSystLow = 1.

        configMgr.hists[highName + 'Norm'] = \
            configMgr.hists[highName].Clone(highName + 'Norm')
        configMgr.hists[lowName + 'Norm'] = \
            configMgr.hists[lowName].Clone(lowName + 'Norm')
        configMgr.hists[nomName + 'Norm'] = \
            configMgr.hists[nomName].Clone(nomName + 'Norm')

        highIntegral = configMgr.hists[highName].Integral()
        lowIntegral = configMgr.hists[lowName].Integral()
        nomIntegral = configMgr.hists[nomName].Integral()

        for iBin in xrange(configMgr.hists[highName + 'Norm'
                           ].GetNbinsX() + 1.):
            try:
                configMgr.hists[highName + 'Norm'].SetBinContent(iBin,
                        fabs(configMgr.hists[highName + 'Norm'
                        ].GetBinContent(iBin)
                         / configMgr.hists[nomName].GetBinContent(iBin)
                         - 1.))
                configMgr.hists[highName + 'Norm'].SetBinError(iBin, 0.)
            except ZeroDivisionError:
                configMgr.hists[highName + 'Norm'].SetBinContent(iBin,
                        0.)
                configMgr.hists[highName + 'Norm'].SetBinError(iBin, 0.)
        for iBin in xrange(configMgr.hists[lowName + 'Norm'].GetNbinsX()
                            + 1.):
            try:
                configMgr.hists[lowName + 'Norm'].SetBinContent(iBin,
                        fabs(configMgr.hists[lowName + 'Norm'
                        ].GetBinContent(iBin)
                         / configMgr.hists[nomName].GetBinContent(iBin)
                         - 1.))
                configMgr.hists[lowName + 'Norm'].SetBinError(iBin, 0.)
            except ZeroDivisionError:
                configMgr.hists[lowName + 'Norm'].SetBinContent(iBin,
                        0.)
                configMgr.hists[lowName + 'Norm'].SetBinError(iBin, 0.)

        for iBin in xrange(configMgr.hists[nomName + 'Norm'].GetNbinsX()
                            + 1.):
            try:
                configMgr.hists[nomName + 'Norm'].SetBinContent(iBin,
                        max(configMgr.hists[highName + 'Norm'
                        ].GetBinContent(iBin), configMgr.hists[lowName
                         + 'Norm'].GetBinContent(iBin)))
                if configMgr.verbose > 1.:
                    print '!!!!!! shapeSys %s bin %g value %g'\
                         % (systName, iBin, configMgr.hists[nomName
                             + 'Norm'].GetBinContent(iBin))
                configMgr.hists[nomName + 'Norm'].SetBinError(iBin, 0.)
            except ZeroDivisionError:
                configMgr.hists[nomName + 'Norm'].SetBinContent(iBin,
                        0.)
                configMgr.hists[nomName + 'Norm'].SetBinError(iBin, 0.)

        if not systName in configMgr.systDict.keys():
            self.systList.append(systName)
        return

    def addOverallSys(
        self,
        systName,
        high,
        low,
        ):
        """
        Add an OverallSys entry using the high and low values
        """

        if high == 1. and low == 1.:
            print '    WARNING: addOverallSys for %s high==1.0 and low==1.0 remove from fit'\
                 % systName
            return

        if high == 0. and low == 0.:
            print '    WARNING: addOverallSys for %s high=%g low=%g remove from fit.'\
                 % (systName, systName, high, low)
            return

# #         if high<1.0 and low<1.0:
# #             highOld=high
# #             high = 2.0 - low
# #             print "WARNING addOverallSys: high=%f is < 1.0 in %s. Taking symmetric value from low %f %f"%(highOld,systName,low,high)

# #         if low>1.0 and high>1.0:
# #             lowOld=low
# #             low = 2.0 - high
# #             print "WARNING addOverallSys: low=%f is > 1.0 in %s. Taking symmetric value from high %f %f"%(lowOld,systName,low,high)

        if high == 1. and low > 0. and low != 1.:
            highOld = high
            high = 2. - low
            print '    WARNING addOverallSys: high=%f in %s. Taking symmetric value from low %f %f'\
                 % (highOld, systName, low, high)

        if low == 1. and high > 0. and high != 1.:
            lowOld = low
            low = 2. - high
            print '    WARNING addOverallSys: low=%f in %s. Taking symmetric value from high %f %f'\
                 % (lowOld, systName, low, high)

        if low < 0.:
            print '    WARNING addOverallSys: low=%f is < 0.0 in %s. Setting to low=0.0. High=%f.'\
                 % (low, systName, high)
            low = 0.

        if high < 0.:
            print '    WARNING addOverallSys: high=%f is < 0.0 in %s. Setting to high=0.0. Low=%f.'\
                 % (high, systName, low)
            high = 0.

        self.overallSystList.append((systName, high, low))
        if not systName in configMgr.systDict.keys():
            self.systList.append(systName)
        return

    def addNormFactor(
        self,
        name,
        val,
        high,
        low,
        const=False,
        ):
        """
        Add a normlization factor
        """

        self.normFactor.append((name, val, high, low, str(const)))
        if not name in configMgr.normList:
            configMgr.normList.append(name)
        return

    def setNormFactor(
        self,
        name,
        val,
        low,
        high,
        const=False,
        ):
        """
        Set normalization factor
        """

        self.normFactor = []
        self.normFactor.append((name, val, high, low, str(const)))
        if not name in configMgr.normList:
            configMgr.normList.append(name)
        return

    def setFileList(self, filelist):
        """
        Set file list for this Sample directly
        """

        self.files = filelist

    def setFile(self, file):
        """
        Set file for this Sample directly
        """

        self.files = [file]

    def propagateFileList(self, fileList):
        """
        Propagate the file list downwards.
        """

        # if we don't have our own file list, use the one given to us

        if self.files == []:
            self.files = fileList

        # we are the leaves of the configmgr->toplevelxml->channel->sample tree,
        # so no propagation necessary

    def addShapeFactor(self, name):
        """
        Bin-by-bin factors to build histogram eg. for data-driven estimates
        """

        self.shapeFactorList.append(name)

    def addSystematic(self, syst):
        """
        Add a systematic to this Sample directly
        """

        if syst.name in self.systDict.keys():
            raise Exception('Attempt to overwrite systematic %s in Sample %s'
                             % (syst.name, self.name))
        else:
            self.systDict[syst.name] = syst.Clone()
            return

    def getSystematic(self, systName):
        try:
            return self.systDict[systName]
        except KeyError:
            raise KeyError('Could not find systematic %s in topLevel %s'
                            % (systName, self.name))

    def removeSystematic(self, name):
        """
        Remove a systematic
        """

        del self.systDict[name]

    def clearSystematics(self):
        """
        Remove a systematic
        """

        self.systDict.clear()

    def __str__(self):
        """
        Convert instance to XML string
        """

        self.sampleString = \
            '  <Sample Name="%s" HistoName="%s" InputFile="%s" NormalizeByTheory="%s">\n'\
             % (self.name, self.histoName, configMgr.histCacheFile,
                self.normByTheory)
        if self.statConfig:
            self.sampleString += '    <StatError Activate="%s"/>\n'\
                 % self.statConfig
        for histoSyst in self.histoSystList:
            self.sampleString += \
                '    <HistoSys Name="%s" HistoNameHigh="%s" HistoNameLow="%s" />\n'\
                 % (histoSyst[0.], histoSyst[1.], histoSyst[2.])
        for shapeSyst in self.shapeSystList:
            self.sampleString += \
                '    <ShapeSys Name="%s" HistoName="%s" ConstraintType="%s"/>\n'\
                 % (shapeSyst[0.], shapeSyst[1.], shapeSyst[2.])
        for overallSyst in self.overallSystList:
            self.sampleString += \
                '    <OverallSys Name="%s" High="%g" Low="%g" />\n'\
                 % (overallSyst[0.], overallSyst[1.], overallSyst[2.])
        for shapeFact in self.shapeFactorList:
            self.sampleString += '    <ShapeFactor Name="%s" />\n'\
                 % shapeFact
        if len(self.normFactor) > 0.:
            for normFactor in self.normFactor:
                self.sampleString += \
                    '    <NormFactor Name="%s" Val="%g" High="%g" Low="%g" Const="%s" />\n'\
                     % (normFactor[0.], normFactor[1.], normFactor[2.],
                        normFactor[3], normFactor[4])
                pass
        self.sampleString += '''  </Sample>

'''
        return self.sampleString

