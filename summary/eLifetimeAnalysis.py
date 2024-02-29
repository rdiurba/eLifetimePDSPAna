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

def combineErrr(lifetime,lifetimeAlt,statErr):
    sceErr=[]
    i=0
    while i<len(statErr):
        sceErr.append(abs(lifetime[i]-lifetimeAlt[i]))

        i=i+1
    return np.array(sceErr)


def eLifetimeToQaQc(lifetimes,lifetimesAlt,statErr,cathode=1):
    qaqcHigh=[]
    qaqcLow=[]
    qaqcMid=[]
    i=0
    diffConst=0.024
    if (cathode==0): diffConst=0.036
    #print(len(statErr),statErr[len(statErr)-1],len(lifetimes))
    while i<len(lifetimes):
        qaqcMidTemp=np.exp(-tmax/lifetimes[i])
        qaqcMidAltTemp=np.exp(-tmax/lifetimesAlt[i])

        diffQaQcAlt=abs(qaqcMidTemp-qaqcMidAltTemp)
        if (abs(lifetimes[i])<abs(statErr[i])): statErr[i]=abs(lifetimes[i])-0.01
        if (qaqcMidTemp>1): qaqcMidTemp=1 
        qaqcMid.append(qaqcMidTemp)
        qaqcHighTemp=np.exp(-tmax/(lifetimes[i]+statErr[i]))+diffConst+diffQaQcAlt
        qaqcLowTemp=(np.exp(-tmax/(lifetimes[i]-statErr[i])))-diffConst-diffQaQcAlt
        if (qaqcHighTemp>1): qaqcHighTemp=1
        qaqcHigh.append(qaqcHighTemp)
        if (qaqcLowTemp<0): qaqcLowTemp=0 
        qaqcLow.append(qaqcLowTemp)
        if (qaqcLowTemp<0.1): 
            print("Oops bad error low<0: "+str(i)+","+str(lifetimes[i])+","+str(lifetimesAlt[i])+","+str(statErr[i]))
            print(qaqcHighTemp,qaqcLowTemp)
        i=i+1    
    return np.array(qaqcMid),np.array(qaqcHigh),np.array(qaqcLow)


def chi2_calcuator(x,y,xerr,yerr):
    i=0
    chi2=0
    while i<len(x):
        num=(x[i]-y[i])**2
        denom=(xerr[i]**2+yerr[i]**2)
        chi2=chi2+(num/denom)
        i=i+1
    return(chi2,len(x))
def closest_dates(date_TPC, date_Pur, lifetime_Pur,lifetime_PurErr):
    matchedLifetime=[]
    matchedLifetimeErr=[]
    for d in date_TPC:
        lowest_delta=timedelta(hours=36)
        lowest_deltaArray=0
        j=0
        for c in date_Pur:
            j=j+1
            if (abs(d-c)<lowest_delta):
                lowest_delta=abs(d-c)
                lowest_deltaArray=j-1
        matchedLifetime.append(lifetime_Pur[lowest_deltaArray])
        matchedLifetimeErr.append(lifetime_PurErr[lowest_deltaArray])
    return(matchedLifetime,matchedLifetimeErr)

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
            if (float(row[2])>200):
                lifetime_pur.append(200)
            else:
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
    return date_pur,qaqc_purAdj,qaqc_purAdjLow,qaqc_purAdjTop, 
def getCRTData(fileName):
    lifetime_CRT=[]
    lifetime_CRTErr=[]
    date_CRT=[]
    runNumber_CRT=[]
    with open(fileName, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:

            if (row[0][0]=="T" or row[0][0]=="M" or int(row[0])==5768):
                continue
            if (float(row[5])<0 or np.isnan(float(row[4]))==1):
                continue
            if (float(row[5])<0 or float(row[5])>200):
                lifetime_CRT.append(200)
                lifetime_CRTErr.append(70)
            else:
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
            if (float(row[-8])<0 or float(row[-8])>200):
                lifetime_cathodeBR.append(200)
                lifetime_cathodeErrBR.append(70)
            else:
                lifetime_cathodeBR.append(float(row[-8]))
                lifetime_cathodeErrBR.append(float(row[-7]))
            if (float(row[3])<0 or float(row[3])>200):
                lifetime_cathodeAverage.append(200)
                lifetime_cathodeErrAvg.append(70)
            else:
                lifetime_cathodeAverage.append(float(row[3]))
                lifetime_cathodeErrAvg.append(float(row[4]))
            if (float(row[-4])<0 or float(row[-4])>200):
                lifetime_cathodeBL.append(200)
                lifetime_cathodeErrBL.append(70)
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
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return
    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
    dates_cathode = matplotlib.dates.date2num(date_cathode)


    crtLifetime_data=getCRTData("crt_lifetimeresults_Nov2018.txt")          
    date_crt=crtLifetime_data[0]
    crtLifetime=crtLifetime_data[1]
    crtLifetimeErr=crtLifetime_data[2]
    crt_runNumber=crtLifetime_data[3]
    
    crtLifetime_dataAltSCE=getCRTData("crt_lifetimeresults_Nov2018_alt.txt")          
    date_crtAlt=crtLifetime_data[0]
    crtLifetimeAlt=crtLifetime_data[1]
    crtLifetimeErrAlt=crtLifetime_data[2]
    crt_runNumberAlt=crtLifetime_data[3]

    qaqcCRTBL, qaqcCRTBLHigh, qaqcCRTBLLow=eLifetimeToQaQc(crtLifetime,crtLifetimeAlt,crtLifetimeErr,cathode=0)
    dates_crt = matplotlib.dates.date2num(date_crt)

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.plot(dates_cathode,qaqcCathodeBL,"o")
    plt.plot(dates_crt,qaqcCRTBL,"o")
    plt.fill_between(dates_cathode, qaqcCathodeBLLow,qaqcCathodeBLHigh, alpha=0.3,label="Cathode-Crossing Sample BL")
    plt.fill_between(dates_crt,qaqcCRTBLLow,qaqcCRTBLHigh, alpha=0.3, label="CRT-TPC Matched Tracks")
    plt.xlim(dates_cathode[0],dates_cathode[-1])
    plt.xticks(np.arange(min(dates_cathode), max(dates_cathode), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathode[-1])
        plt.xticks(np.arange(novDate, max(dates_cathode), 3.0))
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_cathode]
    plt.plot(dates_cathode,tau30ms,"--",color="red")
    tau100ms=np.exp(-2.25/100)
    tau100ms=[tau100ms for i in dates_cathode]
    plt.plot(dates_cathode,tau100ms,"--",color="black",label=r"$\tau$=100 ms")
    plt.title("Minimum Fraction of Charge Recovered for Center of the Detector BL")
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.6,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcBLMiddleCRT_"+str(int(novDates))+".pdf")
    plt.show()
def plotQaQcBL(locationY,locationZ,plane, withPur=0,whichPur="None",withCRT=0,novDates=True):
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
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_NomErr.txt")          
    date_cathode=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_NomErr_alt.txt")          
    date_cathodeAlt=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return
    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
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
        date_crtAlt=crtLifetime_data[0]
        crtLifetimeAlt=crtLifetime_data[1]
        crtLifetimeErrAlt=crtLifetime_data[2]
        crt_runNumberAlt=crtLifetime_data[3]
        
        qaqcCRT, qaqcCRTHigh, qaqcCRTLow=eLifetimeToQaQc(crtLifetime,crtLifetimeAlt,  crtLifetimeErr,cathode=0)
        dates_crt= matplotlib.dates.date2num(date_crt)

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.plot(dates_cathode,qaqcCathodeBL,"o",label=locationYLabel+" Cathode-Crossers")
    plt.fill_between(dates_cathode, qaqcCathodeBLLow,qaqcCathodeBLHigh, alpha=0.3)
    if (withPur==True):
        plt.plot(dates_pur,purQaQc,"v", label=locationYLabel+" Purity Monitor")
        plt.fill_between(dates_pur,purQaQcLow,purQaQcTop, alpha=0.3)    
    if (withCRT==True):
        plt.plot(dates_crt,qaqcCRT,"v", label="CRT-TPC Sample")
        plt.fill_between(dates_crt,qaqcCRTLow,qaqcCRTHigh, alpha=0.3)  
    plt.xlim(dates_cathode[0],dates_cathode[-1])
    plt.xticks(np.arange(min(dates_cathode), max(dates_cathode), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathode[-1])
        plt.xticks(np.arange(novDate, max(dates_cathode), 3.0))
    tau30ms=np.exp(-tmax/30)
    tau30ms=[tau30ms for i in dates_cathode]
    plt.plot(dates_cathode,tau30ms,"--",color="red",label=r"$\tau$=30 ms")
    tau100ms=np.exp(-2.25/100)
    tau100ms=[tau100ms for i in dates_cathode]
    plt.plot(dates_cathode,tau100ms,"--",color="black",label=r"$\tau$=100 ms")
    plt.title("Minimum Fraction of Charge Recovered for the "+str(locationZLabel)+" BL APA")
    plt.legend(loc='lower right', fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.6,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    if (withCRT==True): plt.savefig("qaqcBL_"+locationY+"_"+locationZ+"_plane"+str(plane)+"_withCRT_NovDates"+str(int(novDates))+".pdf")
    else: plt.savefig("qaqcBL_"+locationY+"_"+locationZ+"_plane"+str(plane)+"_"+str(whichPur)+"NovDates"+str(int(novDates))+".pdf")
    plt.show()
def plotQaQcBLvsBR(locationY,locationZ,plane, novDates=True):
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
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_NomErr.txt")          
    date_cathode=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    runNumber=cathodeData[-1]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane"+str(plane)+"_NomErr_alt.txt")          
    date_cathodeAlt=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    if (len(date_cathode)<2 or len(date_cathodeAlt)<2): return

    qaqcCathodeBL, qaqcCathodeBLHigh, qaqcCathodeBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
    qaqcCathodeBR, qaqcCathodeBRHigh, qaqcCathodeBRLow=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,cathode=1)
    dates_cathode = matplotlib.dates.date2num(date_cathode)
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.plot(dates_cathode,qaqcCathodeBL,"o", color="#1f77b4",label="Cathode-Crossing Sample BL at " +str(locationYLabel)+" Height",)
    plt.plot(dates_cathode,qaqcCathodeBR,"v", color="r",label="Cathode-Crossing Sample BR at "+str(locationYLabel)+" Height" )
    plt.fill_between(dates_cathode, qaqcCathodeBLLow,qaqcCathodeBLHigh, alpha=0.3,color="#1f77b4")
    plt.fill_between(dates_cathode, qaqcCathodeBRLow,qaqcCathodeBRHigh, alpha=0.3,color="r")
    plt.xlim(date_cathode[0],dates_cathode[-1])
    plt.xticks(np.arange(min(dates_cathode), max(dates_cathode), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathode[-1])
        plt.xticks(np.arange(novDate, max(dates_cathode), 3.0))
    tau30ms=np.exp(-2.25/30)
    tau30ms=[tau30ms for i in dates_cathode]
    plt.plot(dates_cathode,tau30ms,"--",color="red",label=r"$\tau$=30 ms")
    tau100ms=np.exp(-2.25/100)
    tau100ms=[tau100ms for i in dates_cathode]
    plt.plot(dates_cathode,tau100ms,"--",color="black",label=r"$\tau$=100 ms")
    plt.title("Minimum Fraction of Charge Recovered Measured for the "+str(locationZLabel)+" Set of APAs")
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.6,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcBLvsBR_"+locationY+"_"+locationZ+"_plane"+str(plane)+"_NovDates"+str(int(novDates))+".pdf")
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
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane2_NomErr.txt")          
    date_cathodePlane2=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane2_NomErr_alt.txt")          
    date_cathodeAltPlane2=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]
    qaqcCathodeBLPlane2, qaqcCathodeBLHighPlane2, qaqcCathodeBLLowPlane2=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane1_NomErr.txt")          
    date_cathodePlane1=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane1_NomErr_alt.txt")          
    date_cathodeAltPlane1=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]  
    qaqcCathodeBLPlane1, qaqcCathodeBLHighPlane1, qaqcCathodeBLLowPlane1=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
    cathodeData=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane0_NomErr.txt")          
    date_cathodePlane0=cathodeData[0]
    lifetime_cathodeBL=cathodeData[1]
    lifetime_cathodeErrBL=cathodeData[2]
    lifetime_cathodeBR=cathodeData[3]
    lifetime_cathodeErrBR=cathodeData[4]
    cathodeDataAlt=getDataCathodeCrosser("mich1_lifetimeresults_All_"+str(locationY)+"_"+str(locationZ)+"_plane0_NomErr_alt.txt")          
    date_cathodeAltPlane0=cathodeDataAlt[0]
    lifetime_cathodeBLAlt=cathodeDataAlt[1]
    lifetime_cathodeErrBLAlt=cathodeDataAlt[2]
    lifetime_cathodeBRAlt=cathodeDataAlt[3]
    lifetime_cathodeErrBRAlt=cathodeDataAlt[4]  
    qaqcCathodeBLPlane0, qaqcCathodeBLHighPlane0, qaqcCathodeBLLowPlane0=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())

    plt.plot(dates_cathodePlane2,qaqcCathodeBLPlane2,"o",label="Plane Z")

    plt.plot(dates_cathodePlane1,qaqcCathodeBLPlane1,"v",label="Plane V")
    plt.plot(dates_cathodePlane0,qaqcCathodeBLPlane0,"s",label="Plane U")    
    plt.fill_between(dates_cathodePlane2, qaqcCathodeBLLowPlane2,qaqcCathodeBLHighPlane2, alpha=0.3)

    plt.fill_between(dates_cathodePlane1, qaqcCathodeBLLowPlane1,qaqcCathodeBLHighPlane1, alpha=0.3)

    plt.fill_between(dates_cathodePlane0, qaqcCathodeBLLowPlane0,qaqcCathodeBLHighPlane0, alpha=0.3)

    
    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(min(dates_cathodePlane2), max(dates_cathodePlane2), 3.0))

    tau30ms=np.exp(-2.25/30)
    tau30ms=[tau30ms for i in dates_cathodePlane2]
    plt.plot(dates_cathodePlane2,tau30ms,"--",color="red",label=r"$\tau$=30 ms")
    tau100ms=np.exp(-2.25/100)
    tau100ms=[tau100ms for i in dates_cathodePlane2]
    plt.plot(dates_cathodePlane2,tau100ms,"--",color="black",label=r"$\tau$=100 ms")
    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(min(dates_cathodePlane2), max(dates_cathodePlane2), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))

    plt.title("Minimum Fraction of Charge Recovered for "+str(locationZLabel)+" BL APA at "+str(locationYLabel)+" Height")
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.6,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffPlanes_"+locationY+"_"+locationZ+"_NovDates"+str(int(novDates))+".pdf")
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
    qaqcCathodeBLPlane2, qaqcCathodeBLHighPlane2, qaqcCathodeBLLowPlane2=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
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
    qaqcCathodeBLPlane1, qaqcCathodeBLHighPlane1, qaqcCathodeBLLowPlane1=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)
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
    qaqcCathodeBLPlane0, qaqcCathodeBLHighPlane0, qaqcCathodeBLLowPlane0=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt,  lifetime_cathodeErrBL,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())


    plt.plot(dates_cathodePlane0,qaqcCathodeBLPlane0,"o")
    plt.plot(dates_cathodePlane1,qaqcCathodeBLPlane1,"o")    
    plt.plot(dates_cathodePlane2,qaqcCathodeBLPlane2,"o")
    plt.fill_between(dates_cathodePlane2, qaqcCathodeBLLowPlane2,qaqcCathodeBLHighPlane2, alpha=0.3,label="Top")
    plt.fill_between(dates_cathodePlane1, qaqcCathodeBLLowPlane1,qaqcCathodeBLHighPlane1, alpha=0.3,label="Middle")
    plt.fill_between(dates_cathodePlane0, qaqcCathodeBLLowPlane0,qaqcCathodeBLHighPlane0, alpha=0.3,label="Bottom")

    
    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(min(dates_cathodePlane2), max(dates_cathodePlane2), 3.0))

    tau30ms=np.exp(-2.25/30)
    tau30ms=[tau30ms for i in dates_cathodePlane2]
    plt.plot(dates_cathodePlane2,tau30ms,"--",color="red",label=r"$\tau$=30 ms")
    tau100ms=np.exp(-2.25/100)
    tau100ms=[tau100ms for i in dates_cathodePlane2]
    plt.plot(dates_cathodePlane2,tau100ms,"--",color="black",label=r"$\tau$=100 ms")

    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(min(dates_cathodePlane2), max(dates_cathodePlane2), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))

    plt.title("Minimum Fraction of Charge Recovered Measured for the "+str(locationZLabel)+" BL APA")
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.6,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffHeights_"+locationZ+"_Plane"+str(plane)+"NovDates"+str(int(novDates))+".pdf")
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
    qaqcCathodeBRPlane2, qaqcCathodeBRHighPlane2, qaqcCathodeBRLowPlane2=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,cathode=1)
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
    qaqcCathodeBRPlane1, qaqcCathodeBRHighPlane1, qaqcCathodeBRLowPlane1=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,cathode=1)
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
    qaqcCathodeBRPlane0, qaqcCathodeBRHighPlane0, qaqcCathodeBRLowPlane0=eLifetimeToQaQc(lifetime_cathodeBR,lifetime_cathodeBRAlt,  lifetime_cathodeErrBR,cathode=1)

    dates_cathodePlane2 = matplotlib.dates.date2num(date_cathodePlane2)
    dates_cathodePlane1 = matplotlib.dates.date2num(date_cathodePlane1)
    dates_cathodePlane0 = matplotlib.dates.date2num(date_cathodePlane0)
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())


    plt.plot(dates_cathodePlane0,qaqcCathodeBRPlane0,"o")
    plt.plot(dates_cathodePlane1,qaqcCathodeBRPlane1,"o")    
    plt.plot(dates_cathodePlane2,qaqcCathodeBRPlane2,"o")
    plt.fill_between(dates_cathodePlane2, qaqcCathodeBRLowPlane2,qaqcCathodeBRHighPlane2, alpha=0.3,label="Top")
    plt.fill_between(dates_cathodePlane1, qaqcCathodeBRLowPlane1,qaqcCathodeBRHighPlane1, alpha=0.3,label="Middle")
    plt.fill_between(dates_cathodePlane0, qaqcCathodeBRLowPlane0,qaqcCathodeBRHighPlane0, alpha=0.3,label="Bottom")

    
    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(min(dates_cathodePlane2), max(dates_cathodePlane2), 3.0))

    tau30ms=np.exp(-2.25/30)
    tau30ms=[tau30ms for i in dates_cathodePlane2]
    plt.plot(dates_cathodePlane2,tau30ms,"--",color="red",label=r"$\tau$=30 ms")
    tau100ms=np.exp(-2.25/100)
    tau100ms=[tau100ms for i in dates_cathodePlane2]
    plt.plot(dates_cathodePlane2,tau100ms,"--",color="black",label=r"$\tau$=100 ms")

    plt.xlim(dates_cathodePlane2[0],dates_cathodePlane2[-1])
    plt.xticks(np.arange(min(dates_cathodePlane2), max(dates_cathodePlane2), 3.0))
    if (novDates==True):
        novDate=matplotlib.dates.date2num(datetime(2018,11,1,00,00))
        plt.xlim(novDate,dates_cathodePlane2[-1])
        plt.xticks(np.arange(novDate, max(dates_cathodePlane2), 3.0))

    plt.title("Minimum Fraction of Charge Recovered Measured for the "+str(locationZLabel)+" BR APA")
    plt.legend(loc='lower right',fontsize="small")
    plt.gcf().autofmt_xdate()
    plt.ylim(0.6,1.02)
    plt.ylabel("$Q_{a}/Q_{c}$")
    plt.xticks(rotation=25)
    plt.savefig("qaqcDiffHeights_"+locationZ+"_Plane"+str(plane)+"NovDates"+str(int(novDates))+"_BR.pdf")
    plt.show()
plotMidBL(novDates=True)
plotMidBL(novDates=False)   
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
            plotQaQcBLvsBR(j,k,i,False)
plotQaQcBL("mid","second",2, withPur=0,whichPur=whichPur,withCRT=1,novDates=True)
plotQaQcDiffHeights("second",2, novDates=False)
plotQaQcDiffHeightsBR("second",2, novDates=False)

plotQaQcDiffPlanes("mid","second",False)
plotQaQcBLvsBR("high","first",2,False)
x=[-360+15*i for i in range(0,46)]
y=[1.01101,1.01093,1.01065,1.01065,1.01016,1.00939,1.00712,1.00712,1.00618,1.00576,1.00339,1.00203,1.00203,1.00068,0.998181,0.996836,0.996836,0.99538,0.993998,0.992707,0.991846,0.991846,0.991164,0.990913,0.993319,0.993319,0.994223,0.995122,0.996017,0.996017,0.996744,0.997667,0.998714,1.00021,1.00021,1.00043,1.00173,1.0026,1.0026,1.00463,1.00539,1.00597,1.00642,1.00642,1.00678,1.00695]
plt.xlabel("Drift distance position (cm)")
plt.ylabel("Recombination correction")
plt.plot(x,y,"o")
plt.show()


CRTOnly=1   
purTop=getPurMonData("prm_Top_EFieldCorlifetime_ProtoDUNE-SP-I_data20210120.csv",0)
date_pur=purTop[0]
purQaQcTop=purTop[1]
purQaQcTopLow=purTop[2]
purQaQcTopTop=purTop[3]

purMid=getPurMonData("prm_Mid_EFieldCorlifetime_ProtoDUNE-SP-I_data20210120.csv",0, True)
date_purMid=purMid[0]
purQaQcMid=purMid[1]
purQaQcMidLow=purMid[2]
purQaQcMidTop=purMid[3]

purLow=getPurMonData("prm_Bot_EFieldCorlifetime_ProtoDUNE-SP-I_data20210120.csv",0, True)
date_purLow=purLow[0]
purQaQcLow=purLow[1]
purQaQcLowLow=purLow[2]
purQaQcLowTop=purLow[3]


crtLifetime_data=getCRTData("crt_lifetimeresults_Nov2018.txt")          
date_crt=crtLifetime_data[0]
crtLifetime=crtLifetime_data[1]
crtLifetimeErr=crtLifetime_data[2]
crt_runNumber=crtLifetime_data[3]

crtLifetime_dataAltSCE=getCRTData("crt_lifetimeresults_Nov2018_alt.txt")          
date_crtAlt=crtLifetime_data[0]
crtLifetimeAlt=crtLifetime_data[1]
crtLifetimeErrAlt=crtLifetime_data[2]
crt_runNumberAlt=crtLifetime_data[3]

qaqcCRT, qaqcCRTHigh, qaqcCRTLow=eLifetimeToQaQc(crtLifetime,crtLifetimeAlt,  crtLifetimeErr,cathode=0)
dates_crt= matplotlib.dates.date2num(date_crt)


tpc="second"
plane="2"

cathodeMid=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(tpc)+"_plane"+str(plane)+".txt")          
date_cathodeTop=cathodeMid[0]
lifetime_cathodeBL=cathodeMid[1]
lifetime_cathodeErrBL=cathodeMid[2]
lifetime_cathodeBR=cathodeMid[3]
lifetime_cathodeErrBR=cathodeMid[4]

cathodeMid=getDataCathodeCrosser("mich1_lifetimeresults_All_high_"+str(tpc)+"_plane"+str(plane)+"_alt.txt")          
date_cathodeTop=cathodeMid[0]
lifetime_cathodeBLAlt=cathodeMid[1]
lifetime_cathodeErrBLAlt=cathodeMid[2]
lifetime_cathodeBRAlt=cathodeMid[3]
lifetime_cathodeErrBRAlt=cathodeMid[4]

qaqcBLTop, qaqcBLTopHigh, qaqcBLTopLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt, lifetime_cathodeErrBL,cathode=1)



cathodeMid=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(tpc)+"_plane"+str(plane)+".txt")          
date_cathodeLow=cathodeMid[0]
lifetime_cathodeBL=cathodeMid[1]
lifetime_cathodeErrBL=cathodeMid[2]
lifetime_cathodeBR=cathodeMid[3]
lifetime_cathodeErrBR=cathodeMid[4]

cathodeMid=getDataCathodeCrosser("mich1_lifetimeresults_All_low_"+str(tpc)+"_plane"+str(plane)+"_alt.txt")          
date_cathodeLow=cathodeMid[0]
lifetime_cathodeBLAlt=cathodeMid[1]
lifetime_cathodeErrBLAlt=cathodeMid[2]
lifetime_cathodeBRAlt=cathodeMid[3]
lifetime_cathodeErrBRAlt=cathodeMid[4]

qaqcBLBottom, qaqcBLBottomHigh, qaqcBLBottomLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt, lifetime_cathodeErrBL,cathode=1)



cathodeMid=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(tpc)+"_plane"+str(plane)+".txt")          
date_cathode=cathodeMid[0]
lifetime_cathodeBL=cathodeMid[1]
lifetime_cathodeErrBL=cathodeMid[2]
lifetime_cathodeBR=cathodeMid[3]
lifetime_cathodeErrBR=cathodeMid[4]

cathodeMid=getDataCathodeCrosser("mich1_lifetimeresults_All_mid_"+str(tpc)+"_plane"+str(plane)+"_alt.txt")          
date_cathode=cathodeMid[0]
lifetime_cathodeBLAlt=cathodeMid[1]
lifetime_cathodeErrBLAlt=cathodeMid[2]
lifetime_cathodeBRAlt=cathodeMid[3]
lifetime_cathodeErrBRAlt=cathodeMid[4]

qaqcBL, qaqcBLHigh, qaqcBLLow=eLifetimeToQaQc(lifetime_cathodeBL,lifetime_cathodeBLAlt, lifetime_cathodeErrBL,cathode=1)



cathodeHigh=getDataCathodeCrosser("mich1_lifetimeresults_All_high_second_plane2.txt")          
date_cathodeHigh=cathodeHigh[0]
lifetime_cathodeBLHigh=cathodeHigh[1]
lifetime_cathodeErrBLHigh=cathodeHigh[2]
lifetime_cathodeBRHigh=cathodeHigh[3]
lifetime_cathodeErrBRHigh=cathodeHigh[4]
lifetime_cathodeHigh=cathodeHigh[5]
lifetime_cathodeErrHigh=cathodeHigh[6]

cathodeLow=getDataCathodeCrosser("mich1_lifetimeresults_All_low_second_plane2.txt")          
date_cathodeLow=cathodeLow[0]
lifetime_cathodeBLLow=cathodeLow[1]
lifetime_cathodeErrBLLow=cathodeLow[2]
lifetime_cathodeBRLow=cathodeLow[3]
lifetime_cathodeErrBRLow=cathodeLow[4]
lifetime_cathodeLow=cathodeLow[5]
lifetime_cathodeErrLow=cathodeLow[6]





dates_pur= matplotlib.dates.date2num(date_pur)
dates_purMid= matplotlib.dates.date2num(date_purMid)

dates_purLow= matplotlib.dates.date2num(date_purLow)

dates_cathode = matplotlib.dates.date2num(date_cathode)
dates_cathodeHigh = matplotlib.dates.date2num(date_cathodeHigh)
dates_cathodeTop = matplotlib.dates.date2num(date_cathodeTop)

dates_cathodeLow = matplotlib.dates.date2num(date_cathodeLow)
dates_crt=matplotlib.dates.date2num(date_crt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())

#plt.errorbar(dates_cathodeHigh,lifetime_cathodeBLHigh,lifetime_cathodeErrBLHigh,marker="o", ls="none",label="Cathode-crossers Top")
#plt.errorbar(dates_cathodeLow,lifetime_cathodeBLLow,lifetime_cathodeErrBLLow, marker="v",ls='none', label="Cathode-crossers Low")
plt.errorbar(dates_cathode,lifetime_cathodeBL,lifetime_cathodeErrBL,marker="o",ls='none', label="Cathode-crossers Middle")
plt.errorbar(dates_crt,crtLifetime,crtLifetimeErr,marker="v",ls='none', label="CRT-TPC Matched Tracks")
plt.legend(loc='upper left', title="Statistical Uncertainties")
plt.ylim(0,150)
plt.ylabel("Drift Electron Lifetime (ms)")
plt.xlim(dates_crt[0],dates_crt[-1])
plt.xticks(np.arange(min(dates_crt), max(dates_crt)+1, 3.0))
plt.xticks(rotation=25)
plt.gcf().autofmt_xdate()
plt.savefig("lifetimeMiddle.pdf")
plt.show()


