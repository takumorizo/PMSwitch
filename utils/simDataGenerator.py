#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import re
import copy
import math
import array
import distutils.util
import random
import contextlib
import ConfigParser

import numpy as np
import numpy.random

"""
I : 固定
Js[i]
for i in range(I):
    Js(i) ~ poisson(JLambda)

T  : 固定
N  : 固定
L  : 固定
Ml : 固定

alpha : 固定
beta  : 固定
gamma ~ P_dir([gammaBase]*N)
gamma *= gammaFactor

eta_l ~ P_dir([etaBase]*Ml[l])
eta_l *= etaFactor

f_k_l ~ P_dir([eta_l])

# hyper parameter for g_n_l
phi_l ~ P_dir([phiBase]*Ml[l])
phi_l *= phiFactor

これらのハイパーパラメタを元にして,
x[i,j,l]
g_n_l
を生成する

残りの作業: デバッグ

python ~/Desktop/All/work/sftp_scripts/github/packages/PMSwitch/utils/simDataGenerator.py \
~/Desktop/All/work/sftp_scripts/github/packages/PMSwitch/utils/sampleConfig.ini \
~/Desktop/All/work/sftp_scripts/github/packages/PMSwitch/tests/data/simDB.txt \
~/Desktop/All/work/sftp_scripts/github/packages/PMSwitch/tests/data/simFV.txt \
~/Desktop/All/work/sftp_scripts/github/packages/PMSwitch/tests/data/param.txt

"""

class SignatureSimulator(object):
    def __init__(self, generalConfPath):
        super(SignatureSimulator, self).__init__()
        inifile = ConfigParser.SafeConfigParser()
        inifile.optionxform = str
        inifile.read(generalConfPath)

        self.__prepared = False

        self.__I = inifile.getint('settings', 'I')
        self.__Js = []
        self.__JLambda = inifile.getfloat('settings', 'JLambda')
        for i in range(self.__I):
            self.__Js.append( numpy.random.poisson(self.__JLambda) )


        self.__T = inifile.getint('settings', 'T')
        self.__N = inifile.getint('settings', 'N')
        self.__L = inifile.getint('settings', 'L')
        self.__Ml = []
        for l in range(self.__L):
            m = inifile.getint('settings', 'ml_'+ str(l) )
            self.__Ml.append(m)

        self.__alpha = inifile.getfloat('settings', 'alpha')
        self.__beta0 = inifile.getfloat('settings', 'beta0')
        self.__beta1 = inifile.getfloat('settings', 'beta1')

        self.__etaFactor = [inifile.getfloat('settings', 'etaFactor' )]
        self.__etaBases = []
        self.__eta = []
        for l in range(self.__L):
            self.__etaBases.append( [inifile.getfloat('settings', 'eta' )] * self.__Ml[l] )
            self.__eta.append( numpy.random.dirichlet(self.__etaBases[l]) * self.__etaFactor )


        self.__f = []
        for k in range(self.__T):
            self.__f.append([])

        for k in range(self.__T):
            for l in range(self.__L):
                self.__f[k].append(numpy.random.dirichlet(self.__eta[l]) )

        self.__gammaFactor = inifile.getfloat('settings', 'gammaFactor' )
        self.__gammaBases  = [inifile.getfloat('settings', 'gamma' )] * self.__N
        self.__gamma       = numpy.random.dirichlet(self.__gammaBases) * self.__gammaFactor

        self.__phiFactor = [inifile.getfloat('settings', 'phiFactor' )]
        self.__phiBases  = []
        for l in range(self.__L):
            self.__phiBases.append([inifile.getfloat('settings', 'phi' )] * self.__Ml[l] )

        self.__phi = []
        for l in range(self.__L):
            self.__phi.append(numpy.random.dirichlet(self.__phiBases[l]) * self.__phiFactor)

        self.__g = []
        for n in range(self.__N):
            self.__g.append([])

        for n in range(self.__N):
            for l in range(self.__L):
                self.__g[n].append( numpy.random.dirichlet(self.__phi[l]) )


    def genDBData(self, dbPath):
        with open(dbPath, 'w') as dbFile:
            """
            .db ファイルフォーマット
            N L m1 m2, ... mL
            n 1 g_{n,1,1}, ..., g_{n,1,m1}
            n . .
            n . .
            n L g_{n,L,1}, ..., g_{n,L,mL}
            """
            headerList = [self.__N, self.__L]
            for l in range(self.__L):
                headerList.append(self.__Ml[l])
            dbFile.writelines( self.__makeOutputStr(headerList) )

            for n in range(self.__N):
                for l in range(self.__L):
                    outputList = [n, l]
                    for m in range(self.__Ml[l]):
                        outputList.append( self.__g[n][l][m] )
                    dbFile.writelines( self.__makeOutputStr(outputList) )

    def genFVData(self, fvPath):
        with open(fvPath, 'w') as fvFile:
            """
            .fv ファイルフォーマット
            L m1 ... mL //header
            sample, x1, ..., xL //eachSample
            """
            headerList = [self.__L]
            for l in range(self.__L):
                headerList.append(self.__Ml[l])
            fvFile.writelines( self.__makeOutputStr(headerList) )

            Xpattern = []
            for l in range(self.__L):
                Xpattern.append( [  m for m in range(self.__Ml[l]) ] )

            Ypattern = [ n for n in range(self.__N) ]
            Zpattern = [ k for k in range(self.__T) ]
            Spattern = [0,1]
            for i in range(self.__I):
                theta = numpy.random.dirichlet([self.__beta0, self.__beta1])
                pi    = numpy.random.dirichlet(self.__gamma)
                vList = []

                for k in range(self.__T - 1):
                    vList.append( numpy.random.dirichlet( [1.0, self.__alpha] )[0] )

                sticks = self.__makeStickBreaks(vList)
                for j in range(self.__Js[i]):
                    Z = numpy.random.choice(Zpattern, p=np.asarray(sticks))
                    S = numpy.random.choice(Spattern, p=np.asarray(theta))
                    Y = numpy.random.choice(Ypattern, p=np.asarray(pi))

                    Xlist = ['sample:'+str(i)]
                    for l in range(self.__L):
                        if   S == 0:
                            """
                                X is generated from f
                            """
                            X = numpy.random.choice(Xpattern[l], p=self.__f[Z][l])
                        elif S == 1:
                            """
                                X is generated from g
                            """
                            X = numpy.random.choice(Xpattern[l], p=self.__g[Y][l])
                        Xlist.append(X)

                    line = self.__makeOutputStr(Xlist)
                    fvFile.writelines(line)

    def __makeStickBreaks(self, vList):
        sticks = [1.0] * (len(vList)+1)
        for i in range(len(vList)):
            sticks[i] *= vList[i]
            for ii in range(i):
                sticks[i] *= 1.0 - vList[ii]
        for i in range(len(vList)):
            sticks[len(sticks)-1] *= 1.0 - vList[i]
        return sticks

    def __makeOutputStr(self, outputList, isReturn = True):
        outputString = '\t'.join(map(str, outputList))
        outputString = outputString.replace('\n', '')
        if isReturn:
            outputString = outputString + '\n'
        return outputString

    """
    .db ファイルフォーマット
    T L m1 m2, ... mL
    k 1 f_{k,1,1}, ..., f_{k,1,m1}
    k . .
    k . .
    k L f_{k,L,1}, ..., f_{k,L,mL}
    """
    def genParamData(self, tempParameterPath):
        with open(tempParameterPath, 'w') as paramFile:
            headerList = [self.__T, self.__L]
            for l in range(self.__L):
                headerList.append(self.__Ml[l])
            paramFile.writelines( self.__makeOutputStr(headerList) )

            for k in range(self.__T):
                for l in range(self.__L):
                    outputList = [k, l]
                    for m in range(self.__Ml[l]):
                        outputList.append( self.__f[k][l][m] )
                    paramFile.writelines( self.__makeOutputStr(outputList) )


def main():
    argvs = sys.argv
    config = argvs[1]
    dbPath = argvs[2]
    fvPath = argvs[3]
    tempParameterPath = argvs[4]
    simulator = SignatureSimulator(config)
    simulator.genDBData(dbPath)
    simulator.genFVData(fvPath)
    simulator.genParamData(tempParameterPath)


if __name__ == '__main__':
    main()
