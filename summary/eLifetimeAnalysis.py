# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 19:31:31 2023

@author: Richie-LHEP
"""

from datetime import datetime, timedelta
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
tmax=2.25
textstring="DUNE:ProtoDUNE-SP"

def combineErrr(lifetime,lifetimeAlt,statErr):
    sceErr=[]
    i=0
    while i<len(statErr):
        sceErr.append(abs(lifetime[i]-lifetimeAlt[i]))

        i=i+1
    return np.array(sceErr)

def tauUncertainty(lifetimes,lifetimesAlt,statErr,recombUp, recombDown,cathode=1, statOnly=0):
    qaqcHigh=[]
    qaqcLow=[]
    qaqcMid=[]
    diffConst=0.024
    if (cathode==0): diffConst=0.036
    if (statOnly==1): diffConst=0
    i=0  
    lifetimesHigh=[]
    lifetimesLow=[]
    while i<len(lifetimes):
        sigmaRecombUp=abs(lifetimes[i]-recombUp[i])
        sigmaRecombDown=abs(lifetimes[i]-recombDown[i])
        sigmaStat=statErr[i]
        
        sigmaSCE=abs(lifetimes[i]-lifetimesAlt[i])
        diffLifetime=-tmax/np.log(np.exp(-tmax/lifetimes[i])-diffConst)
        sigmaDiff=abs(lifetimes[i]-diffLifetime)
        if (statOnly==1):
            sigmaRecombDown=0
            sigmaSCE=0 
            sigmaRecombUp=0 
            sigmaDiff=0 
        errHigh=np.sqrt(np.power(sigmaRecombDown,2)+np.power(sigmaSCE,2)+np.power(sigmaStat,2)+np.power(sigmaDiff,2))
        errLow=np.sqrt(np.power(sigmaRecombUp,2)+np.power(sigmaSCE,2)+np.power(sigmaStat,2)+np.power(sigmaDiff,2))
        if (errLow>lifetimes[i]): errLow=lifetimes[i]
        lifetimesHigh.append(lifetimes[i]+errHigh)
        lifetimesLow.append(lifetimes[i]-errLow)
        #print(sigmaStat,sigmaSCE,sigmaRecombUp,sigmaRecombDown)
        i=i+1
    return np.array(lifetimes), np.array(lifetimesHigh), np.array(lifetimesLow)
def eLifetimeToQaQc(lifetimes,lifetimesAlt,statErr,recombUp, recombDown,cathode=1, statOnly=0,printTable=False, runNumbers=[]):
    qaqcHigh=[]
    qaqcLow=[]
    qaqcMid=[]
    i=0
    diffConst=0.024
    if (cathode==0): diffConst=0.036
    if (statOnly==1): diffConst=0
    while i<len(lifetimes):
        qaqcMidTemp=np.exp(-tmax/lifetimes[i])
        qaqcMidAltTemp=np.exp(-tmax/lifetimesAlt[i])
        qaqcMidRecombUp=np.exp(-tmax/recombUp[i])
        qaqcMidRecombDown=np.exp(-tmax/recombDown[i])

        diffQaQcAlt=abs(qaqcMidTemp-qaqcMidAltTemp)
        diffQaQcRecombUp=abs(qaqcMidTemp-qaqcMidRecombUp)
        diffQaQcRecombDown=abs(qaqcMidTemp-qaqcMidRecombDown)
        sigmaRecombUp=abs(lifetimes[i]-recombUp[i])
        sigmaRecombDown=abs(lifetimes[i]-recombDown[i])
        sigmaStat=statErr[i]
        sigmaSCE=abs(lifetimes[i]-lifetimesAlt[i])
        diffLifetime=-tmax/np.log(np.exp(-tmax/lifetimes[i])-diffConst)
        sigmaDiff=abs(lifetimes[i]-diffLifetime)
        if (statOnly==1):
            sigmaRecombDown=0
            sigmaSCE=0 
            sigmaRecombUp=0 
            sigmaDiff=0 
        errHigh=np.sqrt(np.power(sigmaRecombDown,2)+np.power(sigmaSCE,2)+np.power(sigmaStat,2)+np.power(sigmaDiff,2))
        errLow=np.sqrt(np.power(sigmaRecombUp,2)+np.power(sigmaSCE,2)+np.power(sigmaStat,2)+np.power(sigmaDiff,2))
        if (printTable!=True and printTable!=False):
            print(runNumbers[i], " & ", "{:.2f}".format(lifetimes[i])," & ", "{:.2f}".format(errHigh), " & ", "{:.2f}".format(errLow), " & ", "{:.2f}".format(sigmaStat) ," & ", "{:.2f}".format(sigmaSCE)," & ", "{:.2f}".format(sigmaDiff)," & ", "{:.2f}".format(sigmaRecombUp)," & ","{:.2f}".format(sigmaRecombDown), r" \\ \hline " )
        if (qaqcMidTemp>1): qaqcMidTemp=1 
        qaqcMid.append(qaqcMidTemp)
        qaqcHighTemp=qaqcMidTemp+(tmax*errHigh/(lifetimes[i]*lifetimes[i]))
        qaqcLowTemp=qaqcMidTemp-(tmax*errLow/(lifetimes[i]*lifetimes[i]))
        #print(tmax*errHigh/(lifetimes[i]+lifetimes[i]))
        if (qaqcHighTemp>1): qaqcHighTemp=1
        qaqcHigh.append(qaqcHighTemp)
        if (qaqcLowTemp<0): qaqcLowTemp=0 
        qaqcLow.append(qaqcLowTemp)
        
        i=i+1    
    return np.array(qaqcMid),np.array(qaqcHigh),np.array(qaqcLow)


def chi2_calcuator(x,y,xerrHigh,xerrLow,yerrHigh,yerrLow,deltaT):
    i=0
    chi2=0
    count=0
    print(len(x),len(y))
    while i<len(x):
        if (deltaT[i]>timedelta(hours=6)): 
            i=i+1
            continue
        count=count+1
        xerr=xerrHigh[i]
        yerr=yerrLow[i]
        if (x[i]>y[i]): 
            yerr=yerrHigh[i]
            xerr=xerrLow[i]
        num=(x[i]-y[i])**2
        denom=(xerr**2+yerr**2)
        chi2=chi2+(num/denom)
        
        i=i+1
    if (chi2>1000):
        i=0

        while i<len(x):
            if (deltaT[i]>timedelta(hours=6)): 
                i=i+1
                continue
            count=count+1
            xerr=xerrHigh[i]
            yerr=yerrLow[i]
            if (x[i]>y[i]): 
                yerr=yerrHigh[i]
                xerr=xerrLow[i]
            num=(x[i]-y[i])**2
            denom=(xerr**2+yerr**2)
            i=i+1
            print(num/denom)
    return(chi2,count)
def closest_dates(date_TPC, date_Pur):
    matchedLifetime=[]
    matchedLifetimeTPC=[]
    matchedLifetimeErr=[]
    lowest_deltas=[]
    i=0
    for d in date_TPC:
        i=1+i
        lowest_delta=timedelta(hours=6)
        lowest_deltaArray=0
        j=0
        for c in date_Pur:
            j=j+1
            if (abs(d-c)<lowest_delta):
                lowest_delta=abs(d-c)
                if (j-1 in matchedLifetime): continue
                lowest_deltaArray=j-1
        if (lowest_deltaArray!=0):
            
            lowest_deltas.append(lowest_delta)
            matchedLifetime.append(lowest_deltaArray)
            matchedLifetimeTPC.append(i-1)
    return(matchedLifetimeTPC,matchedLifetime, lowest_deltas)
def getPurMonDataLifetime(fileName,date_CRT,invertDate=False):
    lifetime_pur=[]
    lifetime_purLow=[]
    lifetime_purHigh=[]
    lifetime_purErr=[]
    date_pur=[]
    qaqc_pur=[]
    qaqc_purErr=[]
    qaqc_purTop=[]
    qaqc_purLow=[]
    with open(fileName, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if (row[1][0]=="T"):
                continue
            splitsTime=row[1].split(":")
            splitsDate=row[0].split("-")
            if (row[0][2]=="/"):
                splitsDate=row[0].split("/")
            if (invertDate==False): date_pur.append(datetime(int(splitsDate[0]),int(splitsDate[1]),int(splitsDate[2]),int(splitsTime[0]),int(splitsTime[1])))
            else: date_pur.append(datetime(int(splitsDate[2]),int(splitsDate[1]),int(splitsDate[0]),int(splitsTime[0]),int(splitsTime[1])))

            lifetime_pur.append(float(row[2]))
            
            lifetime_purLow.append(float(row[3]))
            lifetime_purHigh.append(float(row[4]))

            qaqc_pur.append(float(row[5]))
            qaqc_purTop.append(float(row[7]))
            qaqc_purLow.append(float(row[6]))

            lifetime_purErr.append(abs(lifetime_purLow[-1]-lifetime_purHigh[-1])/2)
    qaqc_purAdjLow=[np.exp(-tmax/i) for i in lifetime_purLow]
    qaqc_purAdjTop=[np.exp(-tmax/i) for i in lifetime_purHigh]
    qaqc_purAdj=[np.exp(-tmax/i) for i in lifetime_pur]
    return date_pur,lifetime_pur,lifetime_purLow,lifetime_purHigh

def getPurMonData(fileName,date_CRT,invertDate=False):
    lifetime_pur=[]
    lifetime_purLow=[]
    lifetime_purHigh=[]
    lifetime_purErr=[]
    date_pur=[]
    qaqc_pur=[]
    qaqc_purErr=[]
    qaqc_purTop=[]
    qaqc_purLow=[]
    with open(fileName, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if (row[1][0]=="T"):
                continue
            splitsTime=row[1].split(":")
            splitsDate=row[0].split("-")
            if (row[0][2]=="/"):
                splitsDate=row[0].split("/")
            if (invertDate==False): date_pur.append(datetime(int(splitsDate[0]),int(splitsDate[1]),int(splitsDate[2]),int(splitsTime[0]),int(splitsTime[1])))
            else: date_pur.append(datetime(int(splitsDate[2]),int(splitsDate[1]),int(splitsDate[0]),int(splitsTime[0]),int(splitsTime[1])))
            """
            if (float(row[2])>200):
                lifetime_pur.append(200)
            else:
                lifetime_pur.append(float(row[2]))
            """    
            lifetime_pur.append(float(row[2]))
            
            lifetime_purLow.append(float(row[3]))
            lifetime_purHigh.append(float(row[4]))

            qaqc_pur.append(float(row[5]))
            qaqc_purTop.append(float(row[7]))
            qaqc_purLow.append(float(row[6]))

            lifetime_purErr.append(abs(lifetime_purLow[-1]-lifetime_purHigh[-1])/2)
    qaqc_purAdjLow=[np.exp(-tmax/i) for i in lifetime_purLow]
    qaqc_purAdjTop=[np.exp(-tmax/i) for i in lifetime_purHigh]
    qaqc_purAdj=[np.exp(-tmax/i) for i in lifetime_pur]
    return date_pur,qaqc_purAdj,qaqc_purAdjLow,qaqc_purAdjTop 
def getCRTData(fileName):
    lifetime_CRT=[]
    lifetime_CRTErr=[]
    date_CRT=[]
    runNumber_CRT=[]
    with open(fileName, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:

            lifetime_CRT.append(float(row[5]))
            lifetime_CRTErr.append(float(row[6]))
            minute=int(row[1][2:3])
            hr=int(row[1][:2])
            if len(row[1])==3:
                hr=int(row[1][0])
                minute=int(row[1][1])
            if len(row[1])==4 and row[1][0]!=1:
                hr=int(row[1][0])
            if len(row[1])==5 and row[1][0]!=1:
                hr=int(row[1][0])
            runNumber_CRT.append(int(row[0]))
            date_CRT.append(datetime(int(row[2][:4]),int(row[2][4:6]),int(row[2][6:]),hr,minute))
    return date_CRT,lifetime_CRT,lifetime_CRTErr,runNumber_CRT
    

def getDataCathodeCrosser(fileName):
    lifetime_cathodeBL=[]
    lifetime_cathodeErrBL=[]
    lifetime_cathodeBR=[]
    lifetime_cathodeErrBR=[]
    lifetime_cathodeAverage=[]
    lifetime_cathodeErrAvg=[]
    date_cathode=[]
    runNumber_cathode=[]
    with open(fileName, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if (int(row[0])==5768): continue
        
            if (float(row[-8])<0 or float(row[-8])>1000):
                lifetime_cathodeBR.append(200)
                lifetime_cathodeErrBR.append(70)
                a=1
            else:
                lifetime_cathodeBR.append(float(row[-8]))
                lifetime_cathodeErrBR.append(float(row[-7]))
            if (float(row[3])<0 or float(row[3])>1000):
                lifetime_cathodeAverage.append(200)
                lifetime_cathodeErrAvg.append(70)
                a=1
            else:
                lifetime_cathodeAverage.append(float(row[3]))
                lifetime_cathodeErrAvg.append(float(row[4]))
            if (float(row[-4])<0 or float(row[-4])>1000):
                lifetime_cathodeBL.append(200)
                lifetime_cathodeErrBL.append(70)
                a=1
            else:
                lifetime_cathodeBL.append(float(row[-4]))
                lifetime_cathodeErrBL.append(float(row[-3]))
            minute=int(row[2][2:3])
            hr=int(row[2][:2])
            if len(row[2])==3:
                hr=int(row[2][0])
                minute=int(row[2][1])
            if len(row[2])==4 and row[2][0]!=1:
                hr=int(row[2][0])
            if len(row[2])==5 and row[2][0]!=1:
                hr=int(row[2][0])
            runNumber_cathode.append(int(row[0]))
            date_cathode.append(datetime(int(row[1][:4]),int(row[1][4:6]),int(row[1][6:]),hr,minute))  
    return date_cathode, lifetime_cathodeBL,lifetime_cathodeErrBL, lifetime_cathodeBR, lifetime_cathodeErrBR, lifetime_cathodeAverage,lifetime_cathodeErrAvg, runNumber_cathode
def plotMidBL(novDates=True):
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2.txt")          
    date_cathode=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_alt.txt")       
    date_cathodeAlt=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return
    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp,lifetime_cathodeBLRecombDown,cathode=1)
    dates_cathode = matplotlib.dates.date2num(date_cathode)


    crtLifetime_data=getCRTData("crt_lifetimeresults_Nov2018.txt")          
    date_crt=crtLifetime_data[0]
    crtLifetime=crtLifetime_data[1]
    crtLifetimeErr=crtLifetime_data[2]
    crt_runNumber=crtLifetime_data[3]
    
    crtLifetime_dataAltSCE=getCRTData("crt_lifetimeresults_Nov2018_alt.txt")          
    date_crtAlt=crtLifetime_dataAltSCE[0]
    crtLifetimeAlt=crtLifetime_dataAltSCE[1]
    crtLifetimeErrAlt=crtLifetime_dataAltSCE[2]
    crt_runNumberAlt=crtLifetime_dataAltSCE[3]
    crtLifetime_dataRecombUp=getCRTData("crt_lifetimeresults_Nov2018_recombUp.txt")          
    date_crtRecombUp=crtLifetime_dataRecombUp[0]
    crtLifetimeRecombUp=crtLifetime_dataRecombUp[1]
    crtLifetimeErrRecombUp=crtLifetime_dataRecombUp[2]
    crt_runNumberRecombUp=crtLifetime_dataRecombUp[3]

    crtLifetime_dataRecombDown=getCRTData("crt_lifetimeresults_Nov2018_recombDown.txt")          
    date_crtRecombDown=crtLifetime_dataRecombDown[0]
    crtLifetimeRecombDown=crtLifetime_dataRecombDown[1]
    crtLifetimeErrRecombDown=crtLifetime_dataRecombDown[2]
    crt_runNumberRecombDown=crtLifetime_dataRecombDown[3]
    qaqcCRT, qaqcCRTHigh, qaqcCRTLow=eLifetimeToQaQc(crtLifetime,crtLifetimeAlt,  crtLifetimeErr,crtLifetimeRecombUp,crtLifetimeRecombDown,cathode=0)


    dates_crt = matplotlib.dates.date2num(date_crt)
    f, ax = plt.subplots()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())

    purityMatchedDatesTPC,purityMatchedDates,deltaT=closest_dates(date_cathode, date_crt)
    qaqcCRTFilter=[qaqcCRT[i] for i in purityMatchedDates]
    qaqcCRTLowFilter=[qaqcCRT[i]-qaqcCRTLow[i] for i in purityMatchedDates]
    qaqcCRTHighFilter=[qaqcCRTHigh[i]-qaqcCRT[i] for i in purityMatchedDates]
    print("Length:",len(qaqcCRT),len(qaqcCRTFilter))
    qaqcCRT=qaqcCRTFilter
    qaqcCRTHigh=qaqcCRTHighFilter
    qaqcCRTLow=qaqcCRTLowFilter
    qaqcCathodeBLFilter=[qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
    qaqcCathodeBLHighFilter=[qaqcCathodeBLHigh[i]-qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
    qaqcCathodeBLLowFilter=[qaqcCathodeBL[i]-qaqcCathodeBLLow[i] for i in purityMatchedDatesTPC]
    qaqcCathodeBL=qaqcCathodeBLFilter
    qaqcCathodeBLHigh=qaqcCathodeBLHighFilter
    qaqcCathodeBLLow=qaqcCathodeBLLowFilter
    dates_cathode=[dates_cathode[i] for i in purityMatchedDatesTPC]
    dates_crt=[dates_crt[i] for i in purityMatchedDates]
    chi2=chi2_calcuator(qaqcCathodeBL,qaqcCRT, qaqcCathodeBLHigh, qaqcCathodeBLLow, qaqcCRTHigh, qaqcCRTLow,deltaT)
    plt.errorbar(dates_crt,qaqcCRT,yerr=[qaqcCRTLow,qaqcCRTHigh],fmt="v", label="CRT-TPC Matched Muons", capsize=3)
    plt.errorbar(dates_cathode,qaqcCathodeBL,yerr=[qaqcCathodeBLLow,qaqcCathodeBLHigh],fmt="o",label="Middle Cathode-Crossers Muons")



    plt.xlim(dates_cathode[0],dates_cathode[-1])
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,12,00,00)), 3.0))
    
    
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_cathode]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_cathode]
    plt.plot(dates_cathode,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_cathode,tau100ms,":",color="green",label=r"$\tau$=100 ms")

    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathode[-1])
        plt.xticks(np.arange(novDate, max(dates_cathode), 3.0))

    plt.title("")
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,"Center of the Detector BL",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    plt.savefig("qaqcBLMiddleCRT_"+str(int(novDates))+".pdf")
    plt.show()
    
def plotMidLifetimeBL(novDates=True):
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2.txt")          
    date_cathode=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_alt.txt")       
    date_cathodeAlt=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return
    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=tauUncertainty(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp,lifetime_cathodeBLRecombDown,cathode=1)
    dates_cathode = matplotlib.dates.date2num(date_cathode)


    crtLifetime_data=getCRTData("crt_lifetimeresults_Nov2018.txt")          
    date_crt=crtLifetime_data[0]
    dates_crt = matplotlib.dates.date2num(date_crt)
    crtLifetime=crtLifetime_data[1]
    crtLifetimeErr=crtLifetime_data[2]
    crt_runNumber=crtLifetime_data[3]
    
    crtLifetime_dataAltSCE=getCRTData("crt_lifetimeresults_Nov2018_alt.txt")          
    date_crtAlt=crtLifetime_dataAltSCE[0]
    crtLifetimeAlt=crtLifetime_dataAltSCE[1]
    crtLifetimeErrAlt=crtLifetime_dataAltSCE[2]
    crt_runNumberAlt=crtLifetime_dataAltSCE[3]
    crtLifetime_dataRecombUp=getCRTData("crt_lifetimeresults_Nov2018_recombUp.txt")          
    date_crtRecombUp=crtLifetime_dataRecombUp[0]
    crtLifetimeRecombUp=crtLifetime_dataRecombUp[1]
    crtLifetimeErrRecombUp=crtLifetime_dataRecombUp[2]
    crt_runNumberRecombUp=crtLifetime_dataRecombUp[3]

    crtLifetime_dataRecombDown=getCRTData("crt_lifetimeresults_Nov2018_recombDown.txt")          
    date_crtRecombDown=crtLifetime_dataRecombDown[0]
    crtLifetimeRecombDown=crtLifetime_dataRecombDown[1]
    crtLifetimeErrRecombDown=crtLifetime_dataRecombDown[2]
    crt_runNumberRecombDown=crtLifetime_dataRecombDown[3]
    qaqcCRT, qaqcCRTHigh, qaqcCRTLow=tauUncertainty(crtLifetime,crtLifetimeAlt,  crtLifetimeErr,crtLifetimeRecombUp,crtLifetimeRecombDown,cathode=0)
    
    
    purityMatchedDatesTPC,purityMatchedDates,deltaT=closest_dates(date_cathode, date_crt)
    qaqcCRTFilter=[qaqcCRT[i] for i in purityMatchedDates]
    qaqcCRTLowFilter=[qaqcCRT[i]-qaqcCRTLow[i] for i in purityMatchedDates]
    qaqcCRTHighFilter=[qaqcCRTHigh[i]-qaqcCRT[i] for i in purityMatchedDates]
    qaqcCRT=qaqcCRTFilter
    qaqcCRTHigh=qaqcCRTHighFilter
    qaqcCRTLow=qaqcCRTLowFilter
    qaqcCathodeBLFilter=[qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
    qaqcCathodeBLHighFilter=[qaqcCathodeBLHigh[i]-qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
    qaqcCathodeBLLowFilter=[qaqcCathodeBL[i]-qaqcCathodeBLLow[i] for i in purityMatchedDatesTPC]
    qaqcCathodeBL=qaqcCathodeBLFilter
    qaqcCathodeBLHigh=qaqcCathodeBLHighFilter
    qaqcCathodeBLLow=qaqcCathodeBLLowFilter
    dates_cathode=[dates_cathode[i] for i in purityMatchedDatesTPC]
    dates_crt=[dates_crt[i] for i in purityMatchedDates]
    chi2=chi2_calcuator(qaqcCathodeBL,qaqcCRT, qaqcCathodeBLHigh, qaqcCathodeBLLow, qaqcCRTHigh, qaqcCRTLow,deltaT)
    plt.errorbar(dates_crt,qaqcCRT,yerr=[qaqcCRTLow,qaqcCRTHigh],fmt="v", label="CRT-TPC Matched Muons", capsize=3)
    plt.errorbar(dates_cathode,qaqcCathodeBL,yerr=[qaqcCathodeBLLow,qaqcCathodeBLHigh],fmt="o",label="Middle Cathode-Crossers Muons")
    print("Chi2: ", chi2)
    legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))

    dates_crt = matplotlib.dates.date2num(date_crt)
    f, ax = plt.subplots()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())

    plt.errorbar(dates_cathode,qaqcCathodeBL,yerr=[qaqcCathodeBL-qaqcCathodeBLLow,qaqcCathodeBLHigh-qaqcCathodeBL],fmt="o",label="Cathode-Crossing Muons BL")
    plt.errorbar(dates_crt,qaqcCRT,yerr=[qaqcCRT-qaqcCRTLow,qaqcCRTHigh-qaqcCRT],fmt="v", label="CRT-TPC Mathced Muons")

    #plt.fill_between(dates_cathode, qaqcCathodeBLLow,qaqcCathodeBLHigh, alpha=0.3,label="Cathode-Crossing Sample BL")
    #plt.fill_between(dates_crt,qaqcCRTLow,qaqcCRTHigh, alpha=0.3, label="CRT-TPC Matched Tracks")
    plt.xlim(dates_cathode[0],dates_cathode[-1])
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,12,00,00)), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathode[-1])
        plt.xticks(np.arange(novDate, max(dates_cathode), 3.0))

    plt.title("")
    plt.legend(loc='lower right',fontsize="small", title=legendTitle)
    plt.gcf().autofmt_xdate()
    plt.ylabel("$\tau$ [ms]")
    plt.yscale("log")
    plt.ylim(0,1E2)
    plt.xticks(rotation=25)
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,"Center of the Detector BL",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    plt.savefig("lifetimeBLMiddleCRT_"+str(int(novDates))+".pdf")
    plt.show()


    
    
def plotQaQcBL(locationY,locationZ,plane, withPur=0,whichPur="None",withCRT=0,novDates=True,statOnly=0):
    if (locationY!="mid" and locationY!="high" and locationY!="low"): locationY="mid"
    if (locationZ!="first" and locationZ!="second" and locationZ!="third"): locationZ="second"
    if (plane!=0 and plane!=1 and plane!=2): plane=2
    planeLabel="U"
    if (plane==1): planeLabel="V"
    if (plane==2): planeLabel="Z"
    locationYLabel="Top"
    if (locationY=="mid"): locationYLabel="Middle"
    if (locationY=="low"): locationYLabel="Bottom"
    locationZLabel="First"
    if (locationZ=="second"): locationZLabel="2nd"
    elif (locationZ=="third"): locationZLabel="3rd"
    else: locationZLabel="1st"
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathode=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAlt=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return
    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    dates_cathode = matplotlib.dates.date2num(date_cathode)
    if (withPur==True):
        if (whichPur!="Top" and whichPur!="Mid" and whichPur!="Bot"): whichPur="Mid"
        invertDates=0
        if (whichPur=="Mid" or whichPur=="Bot"): invertDates=True
        pur=getPurMonData("prm_"+str(whichPur)+"_EFieldCorlifetime_ProtoDUNE-SP-I_data20210120.csv",0, invertDates)
        date_pur=pur[0]
        purQaQc=pur[1]
        purQaQcLow=pur[2]
        purQaQcTop=pur[3]
        dates_pur= matplotlib.dates.date2num(date_pur)
    if (withCRT==True):
        crtLifetime_data=getCRTData("crt_lifetimeresults_Nov2018.txt")          
        date_crt=crtLifetime_data[0]
        crtLifetime=crtLifetime_data[1]
        crtLifetimeErr=crtLifetime_data[2]
        crt_runNumber=crtLifetime_data[3]
        
        crtLifetime_dataAltSCE=getCRTData("crt_lifetimeresults_Nov2018_alt.txt")          
        date_crtAlt=crtLifetime_dataAltSCE[0]
        crtLifetimeAlt=crtLifetime_dataAltSCE[1]
        crtLifetimeErrAlt=crtLifetime_dataAltSCE[2]
        crt_runNumberAlt=crtLifetime_dataAltSCE[3]
        
        crtLifetime_dataRecombUp=getCRTData("crt_lifetimeresults_Nov2018_recombUp.txt")          
        date_crtRecombUp=crtLifetime_dataRecombUp[0]
        crtLifetimeRecombUp=crtLifetime_dataRecombUp[1]
        crtLifetimeErrRecombUp=crtLifetime_dataRecombUp[2]
        crt_runNumberRecombUp=crtLifetime_dataRecombUp[3]

        crtLifetime_dataRecombDown=getCRTData("crt_lifetimeresults_Nov2018_recombDown.txt")          
        date_crtRecombDown=crtLifetime_dataRecombDown[0]
        crtLifetimeRecombDown=crtLifetime_dataRecombDown[1]
        crtLifetimeErrRecombDown=crtLifetime_dataRecombDown[2]
        crt_runNumberRecombDown=crtLifetime_dataRecombDown[3] 
        crtLifetime_dataRecombDown=crtLifetime
        crtLifetime_dataRecombUp=crtLifetime
        crtLifetimeAlt=crtLifetime
        qaqcCRT, qaqcCRTHigh, qaqcCRTLow=eLifetimeToQaQc(crtLifetime,crtLifetimeAlt,  crtLifetimeErr,crtLifetimeRecombUp,crtLifetimeRecombDown,cathode=0,statOnly=0)
        qaqcCRTStatOnly, qaqcCRTHighStatOnly, qaqcCRTLowStatOnly=eLifetimeToQaQc(crtLifetime,crtLifetimeAlt,  crtLifetimeErr,crtLifetimeRecombUp,crtLifetimeRecombDown,cathode=0,statOnly=statOnly)

        lifetime_cathodeBLAlt=lifetime_cathodeBL
        lifetime_cathodeBLRecombUp=lifetime_cathodeBL
        lifetime_cathodeBLRecombDown=lifetime_cathodeBL
        qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1,statOnly=0)
        qaqcCathodeBLStatOnly, qaqcCathodeBLHighStatOnly, qaqcCathodeBLLowStatOnly=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1,statOnly=statOnly)

        dates_crt= matplotlib.dates.date2num(date_crt)
    f, ax = plt.subplots()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    legendTitle=""
    if (withPur==True):
        yErrLow=[abs(a_i - b_i) for a_i, b_i in zip(purQaQc, purQaQcLow)]
        yErrHigh=[abs(a_i - b_i) for a_i, b_i in zip(purQaQcTop, purQaQc)]
        purityMatchedDatesTPC,purityMatchedDates,deltaT=closest_dates(date_cathode, date_pur)
        qaqcPurFilter=[purQaQc[i] for i in purityMatchedDates]
        qaqcPurLowFilter=[yErrLow[i] for i in purityMatchedDates]
        qaqcPurHighFilter=[yErrHigh[i] for i in purityMatchedDates]
        purQaqcLowFilter=[purQaQcLow[i] for i in purityMatchedDates]
        qaqcPurHigh=qaqcPurHighFilter
        qaqcPurLow=qaqcPurLowFilter
        qaqcPur=qaqcPurFilter
        dates_pur=[dates_pur[i] for i in purityMatchedDates]
        qaqcCathodeBLFilter=[qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBLHighFilter=[qaqcCathodeBLHigh[i]-qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBLLowFilter=[qaqcCathodeBL[i]-qaqcCathodeBLLow[i] for i in purityMatchedDatesTPC]

        qaqcCathodeBL=qaqcCathodeBLFilter
        qaqcCathodeBLHigh=qaqcCathodeBLHighFilter
        qaqcCathodeBLLow=qaqcCathodeBLLowFilter
        dates_cathode=[dates_cathode[i] for i in purityMatchedDatesTPC]

        plt.errorbar(dates_pur,qaqcPur,yerr=[qaqcPurLow,qaqcPurHigh],fmt="v", label=locationYLabel+" Purity Monitor")
        plt.errorbar(dates_cathode,qaqcCathodeBL,yerr=[qaqcCathodeBLLow,qaqcCathodeBLHigh],fmt="o",label=locationYLabel+" Cathode-Crossing Muons")

        chi2=chi2_calcuator(qaqcCathodeBL,qaqcPur, qaqcCathodeBLHigh, qaqcCathodeBLLow, qaqcPurHigh, qaqcPurLow,deltaT)
        print("Chi2: ", chi2)
        legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
    if (novDates==True):
        end=matplotlib.dates.date2num(datetime(2018,11,13,00,00))
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,end)
        plt.xticks(np.arange(novDate, end, 3.0))
        dates_temp=[novDate,end]
    if (withCRT==True):
        
        purityMatchedDatesTPC,purityMatchedDates,deltaT=closest_dates(date_cathode, date_crt)
        qaqcCRTFilter=[qaqcCRT[i] for i in purityMatchedDates]
        qaqcCRTLowFilter=[qaqcCRT[i]-qaqcCRTLow[i] for i in purityMatchedDates]
        qaqcCRTHighFilter=[qaqcCRTHigh[i]-qaqcCRT[i] for i in purityMatchedDates]
        qaqcCRTLowFilterStat=[qaqcCRTStatOnly[i]-qaqcCRTLowStatOnly[i] for i in purityMatchedDates]
        qaqcCRTHighFilterStat=[qaqcCRTHighStatOnly[i]-qaqcCRTStatOnly[i] for i in purityMatchedDates]
        qaqcCRT=qaqcCRTFilter
        qaqcCRTHigh=qaqcCRTHighFilter
        qaqcCRTLow=qaqcCRTLowFilter
        qaqcCathodeBLFilter=[qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBLHighFilter=[qaqcCathodeBLHigh[i]-qaqcCathodeBL[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBLLowFilter=[qaqcCathodeBL[i]-qaqcCathodeBLLow[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBLLowFilterStat=[qaqcCathodeBLStatOnly[i]-qaqcCathodeBLLowStatOnly[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBLHighFilterStat=[qaqcCathodeBLHighStatOnly[i]-qaqcCathodeBLStatOnly[i] for i in purityMatchedDatesTPC]
        qaqcCathodeBL=qaqcCathodeBLFilter
        qaqcCathodeBLHigh=qaqcCathodeBLHighFilter
        qaqcCathodeBLLow=qaqcCathodeBLLowFilter
        dates_cathode=[dates_cathode[i] for i in purityMatchedDatesTPC]
        dates_crt=[dates_crt[i] for i in purityMatchedDates]
        chi2=chi2_calcuator(qaqcCathodeBL,qaqcCRT, qaqcCathodeBLHigh, qaqcCathodeBLLow, qaqcCRTHigh, qaqcCRTLow,deltaT)
        print("Chi2: ", chi2)
        print(purityMatchedDates,purityMatchedDatesTPC,len(purityMatchedDates),len(purityMatchedDatesTPC))
        legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
        if (statOnly==1):
            chi2=chi2_calcuator(qaqcCathodeBL,qaqcCRT, qaqcCathodeBLHighFilterStat, qaqcCathodeBLLowFilterStat, qaqcCRTHighFilterStat, qaqcCRTLowFilterStat,deltaT)
            print("Chi2: ", chi2)
            print(purityMatchedDates,purityMatchedDatesTPC,len(purityMatchedDates),len(purityMatchedDatesTPC))
            legendTitle='$\chi_{stat}^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
        tau10ms=np.exp(-tmax/10)
        tau10ms=[tau10ms for i in dates_temp]
        plt.errorbar(dates_crt,qaqcCRT,yerr=[qaqcCRTLow,qaqcCRTHigh],fmt="v", label="CRT-TPC Matched Muons", capsize=3)
        print("Last Uncertainty: ", qaqcCRTHigh[-10:],qaqcCRT[-10:],qaqcCRTLow[-10:])
        print("Last Uncertainty2: ", qaqcCRTHigh[-10:],qaqcCRT[-10:],qaqcCRTLow[-10:])

        plt.errorbar(dates_cathode,qaqcCathodeBL,yerr=[qaqcCathodeBLLow,qaqcCathodeBLHigh],fmt="o",label="Cathode-Crossing Muons")
    dates_temp=[matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00))]

    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")

    plt.xlim(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)))
    plt.yticks(np.linspace(0.5,1.0,11))
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,12,00,00)), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_temp[-1])
        plt.xticks(np.arange(novDate, max(dates_temp), 1.0))
        dates_temp=[novDate,max(dates_temp)]
    plt.title("")
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationZLabel)+" BL APA",verticalalignment="bottom", horizontalalignment='center',transform=ax.transAxes,fontsize=13)
    plt.legend(loc='lower right', fontsize="small", title=legendTitle, title_fontsize="small")

    plt.gcf().autofmt_xdate()
    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    if (withCRT==True): plt.savefig("qaqcBL_"+locationY+"_"+locationZ+"_plane"+str(plane)+"_withCRT_NovDates"+str(int(novDates))+".pdf", bbox_inches="tight")
    else: plt.savefig("qaqcBL_"+locationY+"_"+locationZ+"_plane"+str(plane)+"_"+str(whichPur)+"NovDates"+str(int(novDates))+".pdf", bbox_inches="tight")
    plt.show()
def plotQaQcBLvsBR(locationY,locationZ,plane, novDates=True,statOnly=0):
    if (locationY!="mid" and locationY!="high" and locationY!="low"): locationY="mid"
    if (locationZ!="first" and locationZ!="second" and locationZ!="third"): locationZ="second"
    if (plane!=0 and plane!=1 and plane!=2): plane=2
    planeLabel="U"
    if (plane==1): planeLabel="V"
    if (plane==2): planeLabel="Z"
    locationYLabel="Top"
    if (locationY=="mid"): locationYLabel="Middle"
    if (locationY=="low"): locationYLabel="Bottom"
    locationZLabel="First"
    if (locationZ=="second"): locationZLabel="2nd"
    elif (locationZ=="third"): locationZLabel="3rd"
    else: locationZLabel="1st"
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathode=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    runNumber=cathodeData[-1]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAlt=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return

    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1,statOnly=0)
    qaqcCathodeBR, qaqcCathodeBRHigh, qaqcCathodeBRLow=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1,statOnly=0)
    qaqcCathodeBLStatOnly, qaqcCathodeBLHighStatOnly, qaqcCathodeBLLowStatOnly=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1,statOnly=statOnly)
    qaqcCathodeBRStatOnly, qaqcCathodeBRHighStatOnly, qaqcCathodeBRLowStatOnly=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1,statOnly=statOnly)
    
    dates_cathode = matplotlib.dates.date2num(date_cathode)
    f, ax = plt.subplots()
    
    
    
    
    purityMatchedDates,purityMatchedDatesBR,deltaT=closest_dates(date_cathode, date_cathode)



    qaqcCathodeBLFilter=[qaqcCathodeBL[i] for i in purityMatchedDates]
    qaqcCathodeBLHighFilter=[qaqcCathodeBLHigh[i]-qaqcCathodeBL[i] for i in purityMatchedDates]
    qaqcCathodeBLLowFilter=[qaqcCathodeBL[i]-qaqcCathodeBLLow[i] for i in purityMatchedDates]
    qaqcCathodeBLHighFilterStat=[qaqcCathodeBLHighStatOnly[i]-qaqcCathodeBLStatOnly[i] for i in purityMatchedDates]
    qaqcCathodeBLLowFilterStat=[qaqcCathodeBLStatOnly[i]-qaqcCathodeBLLowStatOnly[i] for i in purityMatchedDates]
    qaqcCathodeBL=qaqcCathodeBLFilter
    qaqcCathodeBLHigh=qaqcCathodeBLHighFilter
    qaqcCathodeBLLow=qaqcCathodeBLLowFilter
    
    qaqcCathodeBRFilter=[qaqcCathodeBR[i] for i in purityMatchedDatesBR]
    qaqcCathodeBRLowFilter=[qaqcCathodeBRHigh[i]-qaqcCathodeBR[i] for i in purityMatchedDatesBR]
    qaqcCathodeBRHighFilter=[qaqcCathodeBR[i]-qaqcCathodeBRLow[i] for i in purityMatchedDatesBR]
    qaqcCathodeBRLowFilterStat=[qaqcCathodeBRHighStatOnly[i]-qaqcCathodeBRStatOnly[i] for i in purityMatchedDatesBR]
    qaqcCathodeBRHighFilterStat=[qaqcCathodeBRStatOnly[i]-qaqcCathodeBRLowStatOnly[i] for i in purityMatchedDatesBR]
    dates_cathode=[dates_cathode[i] for i in purityMatchedDates]
    
    qaqcCathodeBR=qaqcCathodeBRFilter
    qaqcCathodeBRHigh=qaqcCathodeBRHighFilter
    qaqcCathodeBRLow=qaqcCathodeBRLowFilter
    chi2=chi2_calcuator(qaqcCathodeBL,qaqcCathodeBR, qaqcCathodeBLHigh, qaqcCathodeBLLow, qaqcCathodeBRHigh, qaqcCathodeBRLow,deltaT)
    print("Chi2: ", chi2)
    legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
    if (statOnly==1):
        chi2=chi2_calcuator(qaqcCathodeBL,qaqcCathodeBR, qaqcCathodeBLHighFilterStat, qaqcCathodeBLLowFilterStat, qaqcCathodeBRHighFilterStat, qaqcCathodeBRLowFilterStat,deltaT)
        print("Chi2: ", chi2)
        legendTitle='$\chi_{stat}^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    
   
    plt.errorbar(dates_cathode,qaqcCathodeBL,yerr=[qaqcCathodeBLLow,qaqcCathodeBLHigh],fmt="o", color="#1f77b4",label="Cathode-Crossing Muons BL at " +str(locationYLabel)+" Height",)
    plt.errorbar(dates_cathode,qaqcCathodeBR,yerr=[qaqcCathodeBRLow,qaqcCathodeBRHigh],fmt="v", color="r",label="Cathode-Crossing Muons BR at "+str(locationYLabel)+" Height" )
    end=matplotlib.dates.date2num(datetime(2018,11,13,00,00))
    start=matplotlib.dates.date2num(datetime(2018,10,9,0,00))
    plt.xlim(start,end)
    dates_temp=[start,end]
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathode[-1])
        plt.xticks(np.arange(novDate, max(dates_cathode), 3.0))


    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationZLabel)+" Set of APAs",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    plt.legend(loc='lower right',fontsize="small", title=legendTitle,title_fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.yticks(np.linspace(0.5,1.0,11))

    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcBLvsBR_"+locationY+"_"+locationZ+"_plane"+str(plane)+"_NovDates"+str(int(novDates))+".pdf", bbox_inches="tight")
    plt.show()
def plotQaQcDiffPlanes(locationY,locationZ, novDates=True):
    if (locationY!="mid" and locationY!="high" and locationY!="low"): locationY="mid"
    if (locationZ!="first" and locationZ!="second" and locationZ!="third"): locationZ="second"

    locationYLabel="Top"
    if (locationY=="mid"): locationYLabel="Middle"
    if (locationY=="low"): locationYLabel="Bottom"
    locationZLabel="First"
    if (locationZ=="second"): locationZLabel="2nd"
    elif (locationZ=="third"): locationZLabel="3rd"
    else: locationZLabel="1st"
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane2.txt")          
    date_cathodePlane2=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane2_alt.txt")          
    date_cathodeAltPlane2=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane2_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane2_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4]
    qaqcCathodeBLPlane2, qaqcCathodeBLHighPlane2, qaqcCathodeBLLowPlane2=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane1.txt")          
    date_cathodePlane1=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane1_alt.txt")          
    date_cathodeAltPlane1=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4] 
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane1_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane1_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4]
    qaqcCathodeBLPlane1, qaqcCathodeBLHighPlane1, qaqcCathodeBLLowPlane1=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane0.txt")          
    date_cathodePlane0=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane0_alt.txt")          
    date_cathodeAltPlane0=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]  
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane0_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane0_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4]
    qaqcCathodeBLPlane0, qaqcCathodeBLHighPlane0, qaqcCathodeBLLowPlane0=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    f, ax = plt.subplots()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    
    end=matplotlib.dates.date2num(datetime(2018,11,13,00,00))
    start=matplotlib.dates.date2num(datetime(2018,10,9,0,00))
    plt.xlim(start,end)
    tau30ms=np.exp(-tmax/30)
    dates_temp=[start,end]
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)), 3.0))

    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")


    plt.plot(dates_cathodePlane2,qaqcCathodeBLPlane2,"o",label="Plane Z")

    plt.plot(dates_cathodePlane1,qaqcCathodeBLPlane1,"v",label="Plane V")
    plt.plot(dates_cathodePlane0,qaqcCathodeBLPlane0,"s",label="Plane U")    
    plt.fill_between(dates_cathodePlane2, qaqcCathodeBLLowPlane2,qaqcCathodeBLHighPlane2, alpha=0.3)

    plt.fill_between(dates_cathodePlane1, qaqcCathodeBLLowPlane1,qaqcCathodeBLHighPlane1, alpha=0.3)

    plt.fill_between(dates_cathodePlane0, qaqcCathodeBLLowPlane0,qaqcCathodeBLHighPlane0, alpha=0.3)

    
    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,12,00,00)), 3.0))


    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,12,00,00)), 3.0))

    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))


    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationZLabel)+" BL APA at "+str(locationYLabel)+" Height",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.yticks(np.linspace(0.5,1.0,11))

    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffPlanes_"+locationY+"_"+locationZ+"_NovDates"+str(int(novDates))+".pdf", bbox_inches="tight")
    plt.show()

def plotQaQcDiffHeights(locationZ,plane, novDates=False):
    if (locationZ!="first" and locationZ!="second" and locationZ!="third"): locationZ="second"
    if (plane!=0 and plane!=1 and plane!=2): plane=2

    planeLabel="U"
    if (plane==1): planeLabel="V"
    if (plane==2): planeLabel="Z"
    locationZLabel="First"
    if (locationZ=="second"): locationZLabel="Second"
    if (locationZ=="third"): locationZLabel="Third"
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathodePlane2=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane2=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBLPlane2, qaqcCathodeBLHighPlane2, qaqcCathodeBLLowPlane2=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathodePlane1=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane1=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4] 
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBLPlane1, qaqcCathodeBLHighPlane1, qaqcCathodeBLLowPlane1=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathodePlane0=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane0=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    
    qaqcCathodeBLPlane0, qaqcCathodeBLHighPlane0, qaqcCathodeBLLowPlane0=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    f, ax = plt.subplots()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    end=matplotlib.dates.date2num(datetime(2018,11,13,00,00))
    start=matplotlib.dates.date2num(datetime(2018,10,9,0,00))
    plt.xlim(start,end)
    tau30ms=np.exp(-tmax/30)
    dates_temp=[start,end]
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)), 3.0))
    tau30ms=[tau30ms for i in dates_temp]
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")    
    purityMatchedDatesPlane2,purityMatchedDatesPlane0,deltaT=closest_dates(date_cathodePlane2, date_cathodePlane0)

    qaqcCathodePlane2Filter=[qaqcCathodeBLPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2HighFilter=[qaqcCathodeBLHighPlane2[i]-qaqcCathodeBLPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2LowFilter=[qaqcCathodeBLPlane2[i]-qaqcCathodeBLLowPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2=qaqcCathodePlane2Filter
    qaqcCathodePlane2High=qaqcCathodePlane2HighFilter
    qaqcCathodePlane2Low=qaqcCathodePlane2LowFilter
    
    qaqcCathodePlane0Filter=[qaqcCathodeBLPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0LowFilter=[qaqcCathodeBLHighPlane0[i]-qaqcCathodeBLPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0HighFilter=[qaqcCathodeBLPlane0[i]-qaqcCathodeBLLowPlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane0=[dates_cathodePlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane2=[dates_cathodePlane2[i] for i in purityMatchedDatesPlane2]

    qaqcCathodePlane0=qaqcCathodePlane0Filter
    qaqcCathodePlane0High=qaqcCathodePlane0HighFilter
    qaqcCathodePlane0Low=qaqcCathodePlane0LowFilter
    chi2=chi2_calcuator(qaqcCathodePlane2,qaqcCathodePlane0, qaqcCathodePlane2High, qaqcCathodePlane2Low, qaqcCathodePlane0High, qaqcCathodePlane0Low,deltaT)
    print("Chi2: ", chi2)
    legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
    




    plt.errorbar(dates_cathodePlane2,qaqcCathodePlane2,yerr=[qaqcCathodePlane2Low,qaqcCathodePlane2High],fmt="o",label="Top")
    plt.errorbar(dates_cathodePlane0,qaqcCathodePlane0,yerr=[qaqcCathodePlane0Low,qaqcCathodePlane0High],fmt="s",label="Bottom")

    plt.xlim(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)))
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)), 3.0))

    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))

    plt.title("")
    
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationZLabel)+" BL APA",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    
    plt.legend(loc='lower right',fontsize="small",title=legendTitle)

    plt.gcf().autofmt_xdate()
    plt.yticks(np.linspace(0.5,1.0,11))

    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffHeights_"+locationZ+"_Plane"+str(plane)+"NovDates"+str(int(novDates))+".pdf", bbox_inches="tight")
    plt.show()

def plotQaQcDiffHeightsBR(locationZ,plane, novDates=False):
    if (locationZ!="first" and locationZ!="second" and locationZ!="third"): locationZ="second"
    if (plane!=0 and plane!=1 and plane!=2): plane=2

    planeLabel="U"
    if (plane==1): planeLabel="V"
    if (plane==2): planeLabel="Z"
    locationZLabel="First"
    if (locationZ=="second"): locationZLabel="Second"
    if (locationZ=="third"): locationZLabel="Third"
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathodePlane2=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane2=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBRPlane2, qaqcCathodeBRHighPlane2, qaqcCathodeBRLowPlane2=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathodePlane1=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane1=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4] 
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBRPlane1, qaqcCathodeBRHighPlane1, qaqcCathodeBRLowPlane1=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+".txt")          
    date_cathodePlane0=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane0=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(locationZ)+"_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBRPlane0, qaqcCathodeBRHighPlane0, qaqcCathodeBRLowPlane0=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    f, ax = plt.subplots()
    dates_temp=[matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00))]

    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    purityMatchedDatesPlane2,purityMatchedDatesPlane0,deltaT=closest_dates(date_cathodePlane2, date_cathodePlane0)

    qaqcCathodePlane2Filter=[qaqcCathodeBRPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2HighFilter=[qaqcCathodeBRHighPlane2[i]-qaqcCathodeBRPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2LowFilter=[qaqcCathodeBRPlane2[i]-qaqcCathodeBRLowPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2=qaqcCathodePlane2Filter
    qaqcCathodePlane2High=qaqcCathodePlane2HighFilter
    qaqcCathodePlane2Low=qaqcCathodePlane2LowFilter
    
    qaqcCathodePlane0Filter=[qaqcCathodeBRPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0LowFilter=[qaqcCathodeBRHighPlane0[i]-qaqcCathodeBRPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0HighFilter=[qaqcCathodeBRPlane0[i]-qaqcCathodeBRLowPlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane0=[dates_cathodePlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane2=[dates_cathodePlane2[i] for i in purityMatchedDatesPlane2]

    qaqcCathodePlane0=qaqcCathodePlane0Filter
    qaqcCathodePlane0High=qaqcCathodePlane0HighFilter
    qaqcCathodePlane0Low=qaqcCathodePlane0LowFilter
    chi2=chi2_calcuator(qaqcCathodePlane2,qaqcCathodePlane0, qaqcCathodePlane2High, qaqcCathodePlane2Low, qaqcCathodePlane0High, qaqcCathodePlane0Low,deltaT)
    print("Chi2: ", chi2)
    legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
    




    plt.errorbar(dates_cathodePlane2,qaqcCathodePlane2,yerr=[qaqcCathodePlane2Low,qaqcCathodePlane2High],fmt="o",label="Top")
    plt.errorbar(dates_cathodePlane0,qaqcCathodePlane0,yerr=[qaqcCathodePlane0Low,qaqcCathodePlane0High],fmt="s",label="Bottom")

    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,12,00,00)), 3.0))
    plt.xlim(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))
    plt.title("")
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationZLabel)+" BR APA",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    plt.legend(loc='lower right',fontsize="small", title=legendTitle)
    plt.gcf().autofmt_xdate()
    plt.yticks(np.linspace(0.5,1.0,11))

    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffHeights_"+locationZ+"_Plane"+str(plane)+"NovDates"+str(int(novDates))+"_BR.pdf", bbox_inches="tight")
    plt.show()

    
def plotQaQcDiffAPAs(locationY,plane, novDates=False):
    #if (locationZ!="first" and locationZ!="second" and locationZ!="third"): locationZ="second"
    if (plane!=0 and plane!=1 and plane!=2): plane=2

    planeLabel="U"
    if (plane==1): planeLabel="V"
    if (plane==2): planeLabel="Z"
    locationYLabel="Top"
    if (locationY=="mid"): locationYLabel="Middle"
    if (locationY=="low"): locationYLabel="Bottom"
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+".txt")          
    date_cathodePlane2=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane2=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBLPlane2, qaqcCathodeBLHighPlane2, qaqcCathodeBLLowPlane2=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+".txt")          
    date_cathodePlane1=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane1=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4] 
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBLPlane1, qaqcCathodeBLHighPlane1, qaqcCathodeBLLowPlane1=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+".txt")          
    date_cathodePlane0=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane0=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    
    qaqcCathodeBLPlane0, qaqcCathodeBLHighPlane0, qaqcCathodeBLLowPlane0=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp, lifetime_cathodeBLRecombDown,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    f, ax = plt.subplots()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    dates_temp=[matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00))]
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")
    purityMatchedDatesPlane2,purityMatchedDatesPlane0,deltaT=closest_dates(date_cathodePlane2, date_cathodePlane0)

    qaqcCathodePlane2Filter=[qaqcCathodeBLPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2HighFilter=[qaqcCathodeBLHighPlane2[i]-qaqcCathodeBLPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2LowFilter=[qaqcCathodeBLPlane2[i]-qaqcCathodeBLLowPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2=qaqcCathodePlane2Filter
    qaqcCathodePlane2High=qaqcCathodePlane2HighFilter
    qaqcCathodePlane2Low=qaqcCathodePlane2LowFilter
    
    qaqcCathodePlane0Filter=[qaqcCathodeBLPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0LowFilter=[qaqcCathodeBLHighPlane0[i]-qaqcCathodeBLPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0HighFilter=[qaqcCathodeBLPlane0[i]-qaqcCathodeBLLowPlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane0=[dates_cathodePlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane2=[dates_cathodePlane2[i] for i in purityMatchedDatesPlane2]

    qaqcCathodePlane0=qaqcCathodePlane0Filter
    qaqcCathodePlane0High=qaqcCathodePlane0HighFilter
    qaqcCathodePlane0Low=qaqcCathodePlane0LowFilter
    chi2=chi2_calcuator(qaqcCathodePlane2,qaqcCathodePlane0, qaqcCathodePlane2High, qaqcCathodePlane2Low, qaqcCathodePlane0High, qaqcCathodePlane0Low,deltaT)
    print("Chi2: ", chi2)
    legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
    




    plt.errorbar(dates_cathodePlane2,qaqcCathodePlane2,yerr=[qaqcCathodePlane2Low,qaqcCathodePlane2High],fmt="o",label="Third APA")
    plt.errorbar(dates_cathodePlane0,qaqcCathodePlane0,yerr=[qaqcCathodePlane0Low,qaqcCathodePlane0High],fmt="s",label="First APA")

    
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)), 3.0))
    plt.xlim(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))

    plt.title("")
    
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationYLabel)+" BL APA",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    
    plt.legend(loc='lower right',fontsize="small",title=legendTitle)
    plt.gcf().autofmt_xdate()
    plt.yticks(np.linspace(0.5,1.0,11))

    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffAPAs_"+locationY+"_Plane"+str(plane)+"NovDates"+str(int(novDates))+".pdf", bbox_inches="tight")
    plt.show()
    
def plotQaQcDiffAPAsBR(locationY,plane, novDates=False):
    if (plane!=0 and plane!=1 and plane!=2): plane=2
    locationYLabel="Top"
    if (locationY=="mid"): locationYLabel="Middle"
    if (locationY=="low"): locationYLabel="Bottom"
    planeLabel="U"
    if (plane==1): planeLabel="V"
    if (plane==2): planeLabel="Z"

    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+".txt")          
    date_cathodePlane2=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane2=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_third_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBRPlane2, qaqcCathodeBRHighPlane2, qaqcCathodeBRLowPlane2=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+".txt")          
    date_cathodePlane1=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane1=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4] 
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_second_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBRPlane1, qaqcCathodeBRHighPlane1, qaqcCathodeBRLowPlane1=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+".txt")          
    date_cathodePlane0=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+"_alt.txt")          
    date_cathodeAltPlane0=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+"_recombUp.txt")          
    date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
    lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
    lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
    lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
    lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
    cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_first_plane"+str(plane)+"_recombDown.txt")          
    date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
    lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
    lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
    lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
    lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
    qaqcCathodeBRPlane0, qaqcCathodeBRHighPlane0, qaqcCathodeBRLowPlane0=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,lifetime_cathodeBRRecombUp, lifetime_cathodeBRRecombDown,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    f, ax = plt.subplots()
    dates_temp=[matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00))]
    tau30ms=np.exp(-2.25/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_temp]
    tau100ms=np.exp(-tmax/100)
    tau100ms=[tau100ms for i in dates_temp]
    plt.plot(dates_temp,tau30ms,"--",color="red",label=r"$\tau$=30 ms")    
    plt.plot(dates_temp,tau100ms,":",color="green",label=r"$\tau$=100 ms")

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())

    purityMatchedDatesPlane2,purityMatchedDatesPlane0,deltaT=closest_dates(date_cathodePlane2, date_cathodePlane0)

    qaqcCathodePlane2Filter=[qaqcCathodeBRPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2HighFilter=[qaqcCathodeBRHighPlane2[i]-qaqcCathodeBRPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2LowFilter=[qaqcCathodeBRPlane2[i]-qaqcCathodeBRLowPlane2[i] for i in purityMatchedDatesPlane2]
    qaqcCathodePlane2=qaqcCathodePlane2Filter
    qaqcCathodePlane2High=qaqcCathodePlane2HighFilter
    qaqcCathodePlane2Low=qaqcCathodePlane2LowFilter
    
    qaqcCathodePlane0Filter=[qaqcCathodeBRPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0LowFilter=[qaqcCathodeBRHighPlane0[i]-qaqcCathodeBRPlane0[i] for i in purityMatchedDatesPlane0]
    qaqcCathodePlane0HighFilter=[qaqcCathodeBRPlane0[i]-qaqcCathodeBRLowPlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane0=[dates_cathodePlane0[i] for i in purityMatchedDatesPlane0]
    dates_cathodePlane2=[dates_cathodePlane2[i] for i in purityMatchedDatesPlane2]

    qaqcCathodePlane0=qaqcCathodePlane0Filter
    qaqcCathodePlane0High=qaqcCathodePlane0HighFilter
    qaqcCathodePlane0Low=qaqcCathodePlane0LowFilter
    chi2=chi2_calcuator(qaqcCathodePlane2,qaqcCathodePlane0, qaqcCathodePlane2High, qaqcCathodePlane2Low, qaqcCathodePlane0High, qaqcCathodePlane0Low,deltaT)
    print("Chi2: ", chi2)
    legendTitle='$\chi^{2}/df$='+str(int(chi2[0]))+"/"+str(int(chi2[1]))
    




    plt.errorbar(dates_cathodePlane2,qaqcCathodePlane2,yerr=[qaqcCathodePlane2Low,qaqcCathodePlane2High],fmt="o",label="Third APA")
    plt.errorbar(dates_cathodePlane0,qaqcCathodePlane0,yerr=[qaqcCathodePlane0Low,qaqcCathodePlane0High],fmt="s",label="First APA")






    
    plt.xticks(np.arange(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)), 3.0))
    plt.xlim(matplotlib.dates.date2num(datetime(2018,10,9,0,00)), matplotlib.dates.date2num(datetime(2018,11,13,00,00)))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))
    plt.title("")
    ax.text(0.22, 1.0, textstring,
        verticalalignment="bottom", horizontalalignment='center',
        transform=ax.transAxes,fontsize=13,fontweight='bold')
    ax.text(0.75,1.0,str(locationYLabel)+" BR APA",verticalalignment="bottom", horizontalalignment='center',
            transform=ax.transAxes,fontsize=13)
    plt.legend(loc='lower right',fontsize="small", title=legendTitle)
    plt.gcf().autofmt_xdate()
    plt.yticks(np.linspace(0.5,1.0,11))

    plt.ylim(0.58,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffAPAs_"+locationY+"_Plane"+str(plane)+"NovDates"+str(int(novDates))+"_BR.pdf", bbox_inches="tight")
    plt.show()

#plotMidBL(novDates=True)
#plotMidBL(novDates=False)  
#plotMidLifetimeBL(novDates=True)
#plotMidLifetimeBL(novDates=False)  

locationsY=["mid","high","low"]
locationsZ=["first","second","third"]
planes=[2]
for i in planes:
    for j in locationsY:
        whichPur="Mid"
        if (j=="high"): whichPur="Top"
        if (j=="low"): whichPur="Bot"
        for k in locationsZ:
            plotQaQcBL(j,k,i, withPur=1,whichPur=whichPur,withCRT=0,novDates=False)
            plotQaQcBLvsBR(j,k,i,False,1)

whichPur="Mid"
plotQaQcBL("mid","second",2, withPur=0,whichPur=whichPur,withCRT=1,novDates=True,statOnly=1)

plotQaQcDiffHeights("second",2, novDates=False)
plotQaQcDiffHeightsBR("second",2, novDates=False)
plotQaQcDiffHeights("third",2, novDates=False)
plotQaQcDiffHeightsBR("third",2, novDates=False)
plotQaQcDiffHeights("first",2, novDates=False)
plotQaQcDiffHeightsBR("first",2, novDates=False)
plotQaQcDiffHeights("second",2,novDates=True)


plotQaQcDiffAPAsBR("high",2, novDates=False)
plotQaQcDiffAPAs("high",2, novDates=False)
plotQaQcDiffAPAsBR("mid",2, novDates=False)
plotQaQcDiffAPAs("mid",2, novDates=False)
plotQaQcDiffAPAsBR("low",2, novDates=False)
plotQaQcDiffAPAs("low",2,novDates=False)

plotQaQcDiffPlanes("mid","second",False)
plotQaQcBLvsBR("high","first",2,False)
x=[-360+15*i for i in range(0,46)]
y=[1.01101,1.01093,1.01065,1.01065,1.01016,1.00939,1.00712,1.00712,1.00618,1.00576,1.00339,1.00203,1.00203,1.00068,0.998181,0.996836,0.996836,0.99538,0.993998,0.992707,0.991846,0.991846,0.991164,0.990913,0.993319,0.993319,0.994223,0.995122,0.996017,0.996017,0.996744,0.997667,0.998714,1.00021,1.00021,1.00043,1.00173,1.0026,1.0026,1.00463,1.00539,1.00597,1.00642,1.00642,1.00678,1.00695]
plt.xlabel("Drift distance position (cm)")
plt.ylabel("Recombination correction")
plt.plot(x,y,"o")
plt.show()


cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2.txt")          
date_cathode=cathodeData[0]
lifetime_cathodeBL=cathodeData[1]
lifetime_cathodeErrBL=cathodeData[2]
lifetime_cathodeBR=cathodeData[3]
lifetime_cathodeErrBR=cathodeData[4]
runNumber=cathodeData[-1]
cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_alt.txt")       
date_cathodeAlt=cathodeDataAlt[0]
lifetime_cathodeBLAlt=cathodeDataAlt[1]
lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
lifetime_cathodeBRAlt=cathodeDataAlt[3]
lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
cathodeDataRecombUp=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_recombUp.txt")          
date_cathodeRecombUpPlane0=cathodeDataRecombUp[0]
lifetime_cathodeBLRecombUp=cathodeDataRecombUp[1]
lifetime_cathodeErrBLRecombUp=cathodeDataRecombUp[2]
lifetime_cathodeBRRecombUp=cathodeDataRecombUp[3]
lifetime_cathodeErrBRRecombUp=cathodeDataRecombUp[4] 
cathodeDataRecombDown=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_second_plane2_recombDown.txt")          
date_cathodeRecombDownPlane0=cathodeDataRecombDown[0]
lifetime_cathodeBLRecombDown=cathodeDataRecombDown[1]
lifetime_cathodeErrBLRecombDown=cathodeDataRecombDown[2]
lifetime_cathodeBRRecombDown=cathodeDataRecombDown[3]
lifetime_cathodeErrBRRecombDown=cathodeDataRecombDown[4] 
qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,lifetime_cathodeBLRecombUp,lifetime_cathodeBLRecombDown,cathode=1,printTable=True,runNumbers=runNumber)
