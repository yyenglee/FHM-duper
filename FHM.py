"""
Created on Fri Feb 22 17:53:14 2019
@author: leeyy
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def usage():
    print("python FHM.py")

def readVCFtoSNP(fnSNP):
    source = open(fnSNP,"r")
    snpList = {'hom':[], 'het':[]}
    for line in source:
        if line[0] == "#":
            continue
        lineInfo = line.strip().split("\t")
        if not lineInfo[9] == "./." and not lineInfo[9] == "0/0": 
            SNPInfo = lineInfo[9].split(":")
            cov = int(SNPInfo[2])
            if cov >= 5:
                chrom = lineInfo[0]
                pos = lineInfo[1]
                snpID = chrom + "_" + pos
                if SNPInfo[0] == "1/1":
                    snpList['hom'].append(snpID)
                elif SNPInfo[0] == "0/1":
                    snpList['het'].append(snpID)
    return snpList

## Return a dict showing 0: common snv, 1: unique snv in lakeview, 2:unique snv in bio14.6
def defineSNVstrain(snvp1,snvp2):
    snvDict = {}
    commSnvHom = list(set(snvp1['hom'])&set(snvp2['hom']))
    combineSnvHet = list(set(snvp1['het'] + snvp2['het']))
    
    for s in commSnvHom:
        snvDict[s] = 2
    for s in combineSnvHet:
        snvDict[s] = 0.5
    for s in list(set(snvp1['hom']) - set(snvp2['hom']) - set(snvp2['het'])):   # Label snp unique in snvp1 as 1
        snvDict[s] = 1
    for s in list(set(snvp2['hom']) - set(snvp1['hom']) - set(snvp1['het'])):   # Label snp unique in snvp2 as 0
        snvDict[s] = 0
    return snvDict

def readContigSize(contigfn):
    contigDict = {}
    source = open(contigfn,"r")
    for line in source:
        info = line.strip().split("\t")
        cSize = int(info[2])
        contigDict[info[0]] = [cSize, math.ceil(cSize/WINSIZE)]
    return contigDict
    
def windowHomScore(iStart,iEnd,contigScoreDict):
    return np.mean(contigScoreDict[iStart:iEnd])

## Global variables
WINFRAME = [100]    ## list: window size (int) to caclculate homozygosity score
highcutoff = 0.9    ## float: homozygous score cutoff for Lakeview linkage
lowcutoff = 0.1     ## float: homozygous score cutoff for Bio14.6 linkage
plotting = False    ## boolean: if plotting the homozygosity score figure as output

## working path and path to vcf file
pdir = "~/duper/finalize_analysis_20210712/"
snvfn1 = pdir + "SNP_withsmallercontig/F0LakeviewF.pass.vcf" # lakeview
snvfn2 = pdir + "SNP_withsmallercontig/F0Bio14_6M.pass.vcf" # bio14.6
mutSNV = pdir + "SNP_withsmallercontig/Pool45.pass.vcf" # duper

## penalizedCommonHomSNPs - boolean:True/False, whether to include common homozygous SNPs from parental strains
# True to exclude - pro: regions with common homozygous will be included if duper pool match with the few Lakeview unique SNPs
# False to include - pro: selected regions are with higher confident, region with densed common homozygous snp (likely mapping error or reference error) will be excluded
penalizedCommonHomSNPs = False

## read snp from parental strains
snvp1 = readVCFtoSNP(snvfn1)
snvp2 = readVCFtoSNP(snvfn2)
snvDict = defineSNVstrain(snvp1, snvp2)

## Information
SNPtype = {"1/1":"hom", "0/1":"het", "1/2":"het_notRef"}
source = open(mutSNV, "r")

## assign score to each snp in duper pool with coverage >=5
snvMappedDict = {}
contigScoreDict = {}
duperSnpDict = {}
homScoreDict = {}
for line in source:
    if line[0] == "#":
        continue
    lineInfo = line.strip().split("\t")
    if not lineInfo[9] == "./." and not lineInfo[9] == "0/0": 
        SNPInfo = lineInfo[9].split(":")
        cov = int(SNPInfo[2])
        if cov >= 5:
            chrom = lineInfo[0]
            pos = lineInfo[1]
            snpID = chrom + "_" + pos
            snvType = SNPtype[SNPInfo[0]]
            try:
                ## Linked to parental if snp is found(inherited) in the parental Snps pool
                SNPsource = snvDict[snpID]
                
                if penalizedCommonHomSNPs or SNPsource < 2:        
                    ##  homozygous snp in both parental strands are skipped if penalizedCommonHomSNPs=FALSE, otherwise, assigned a score of 0.5 and increase a snp in the window
                    if snvType == "hom":
                        if SNPsource == 2:
                            score = 0.5
                        else:
                            score = SNPsource
                    elif snvType == "het":
                        ## reduce score but keep the trend for het SNPs in duper pool
                        if SNPsource == 0:
                            score = 0.2
                        elif SNPsource == 1:
                            score = 0.8
                        else:
                            score = 0.5
                
                    try:
                        snvMappedDict[chrom].append(int(pos))
                        contigScoreDict[chrom].append(score)     
                    except KeyError:
                        snvMappedDict[chrom] = [int(pos)]
                        contigScoreDict[chrom] = [score]

            except KeyError:
                SNPsource = -1      ## -1 for unique snv in duper progenitors
                try:
                    duperSnpDict[chrom].append(int(pos))
                except KeyError:
                    duperSnpDict[chrom] = [int(pos)]

source.close()

## calculate average score in each window
for WINSIZE in WINFRAME:
    outFN = pdir + "homoScore_win" + str(WINSIZE) + "_out.txt"
    wf = open(outFN,"w")
    wf.write("#chr\tstart\tend\thom.score\n")
    j = 1
    n = 0
    fig = []
    fig.append(plt.figure())
    fig[n].subplots_adjust(hspace=0.4, wspace=0.4)
    for chrom in list(contigScoreDict.keys()):
        nSNV = len(contigScoreDict[chrom])
        if nSNV <= WINSIZE:           
            continue
        
        homScoreDict[chrom] = []
        for i in range(nSNV-WINSIZE):
            homScoreDict[chrom].append(windowHomScore(i,i+WINSIZE,contigScoreDict[chrom]))
            wf.write("%s\t%s\t%s\t%s\n"%(chrom,snvMappedDict[chrom][i],snvMappedDict[chrom][i+WINSIZE],homScoreDict[chrom][i]))
        ax = fig[n].add_subplot(2,2,j)
        ax.plot([snvMappedDict[chrom][0]]+snvMappedDict[chrom][WINSIZE:],[homScoreDict[chrom][0]]+homScoreDict[chrom],'b-')
        ax.text(snvMappedDict[chrom][0],0.95,chrom)
        plt.ylim(0,1.1)
        j += 1
        if j == 5:
            fig.append(plt.figure())
            n += 1
            j = 1
    wf.close()

    if plotting == True:
        for i in range(len(fig)):
            fig[i].savefig(pdir + "/homoScoreFig/homoScore_win"+str(WINSIZE)+"_"+str(i)+".pdf")  ## please manually create a folder homoScoreFig at the working directory

## merge windows if homozygous score fit cutoff
for varN in WINFRAME:
    inFile = pdir + "homoScore_win" + str(varN) + "_out.txt"
    outFile = pdir + "homoScore_win" + str(varN) + "_merged_out.txt"
    curContig = -1
    curStart = -1
    curEnd = -1
    curGroup = -1

    source = open(inFile, "r")
    wf = open(outFile, "w")
    wf.write("#chr\tstart\tend\tsize\thom.score\n")

    for line in source:
        if line[0] == "#":
            continue
        info = line.strip().split("\t")
        contig = info[0]
        start = int(info[1])
        end = int(info[2])
        if float(info[3]) <= lowcutoff:
            group = 0
        elif float(info[3]) >= highcutoff:
            group = 1
        else:
            group = 0.5
        if curContig == -1:
            curContig = contig
            curStart = start
            curEnd = end
            curGroup = group

        elif curContig == contig and curGroup == group:
            curEnd = end
        else:
            wf.write("%s\t%s\t%s\t%s\t%s\n" % (curContig, curStart, curEnd, str(curEnd - curStart + 1), curGroup))
            curContig = contig
            curStart = start
            curEnd = end
            curGroup = group
    wf.close()
    source.close()
