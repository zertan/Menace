#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import numpy as np
#from scipy import signal
from lmfit.models import Model
#from lmfit import conf_interval
import sys
import pandas
from pandas import DataFrame
import matplotlib.pyplot as plt
import xmltodict
import os
from scipy.fftpack import fft, ifft
from scipy.interpolate import UnivariateSpline

def mv_filt(L,omega):
    return (1/float(L))*(1-np.exp(-omega*1j*L))/(1-np.exp(-omega*1j))

def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[i:i+n]

def chunks2(x, b, s):
	l=len(x)
	b2=int(np.floor(float(b)/float(2)))
	#out=np.zeros(int(np.floor(float(l)/float(s))))
	for i in range(0, int(np.floor(float(l)/float(s)))):
		if (i*s+b2>l):
			r=i*s+b2-l
			x2=x[i*s-b2:i*s+b2-r]
			x2.extend(x[0:r])
		elif (i*s-b2<0):
			r=-i*s+b2
			x2=x[0:i*s+b2]
			x2[:0]=x[l-r:l]
		else:
			x2=x[i*s-b2:i*s+b2]
		yield x2
		#print(str(r))
		#out[i]=np.sum(x2)
#	return out

def chunks3(x, b, s):
	l=len(x)
	#b2=int(np.floor(float(b)/float(2)))
	#out=np.zeros(int(np.floor(float(l)/float(s))))
	for i in range(0, int(np.floor(float(l)/float(s)))):
		if (i*s+b>l):
			r=i*s+b-l
			x2=x[i*s:i*s+b-r]
			x2.extend(x[0:r])
		else:
			x2=x[i*s:i*s+b]
		yield x2

def binData(x,binSize):
	l=np.ceil((len(x)/binSize))
	out=np.zeros(np.int32(l)+1)
	tmp=chunks(x,binSize)
	for i,val in enumerate(tmp):
		out[i]=np.sum(val)
	return out

def movingAverageC(x,binSize,slide):
	l=len(x)
	iEnd=np.floor(l/slide)
	out=np.zeros(int(iEnd));

	for i,val in enumerate(chunks2(x,binSize,slide)):
		out[i]=np.sum(val)
	return out

def movingMedianC(x,binSize,slide):
    out=np.zeros(int(np.floor(len(x)/slide)));

    for i,val in enumerate(chunks2(x,binSize,slide)):
        out[i]=np.median(val)
    return out

#def movingMedianC(x,binSize,slide):
#	l=len(x)
#	iEnd=np.floor(l/slide)
#	out=chunks2(x,binSize,slide);
#
#	for i,val in enumerate(out):
#		r=np.mod(i*slide+binSize,l)
#		tmp=x[i*slide:i*slide+binSize-r]
#		out[i]=np.median(tmp)
#		#out[i]=np.median(np.concatenate((tmp,r)))
#	return out


def movingAverage(x,binSize,slide):
	l=np.floor((len(x)-binSize)/slide)
	out=np.zeros(l);
		
	for i,val in enumerate(out):
		out[i]=np.sum(x[i*slide:i*slide+binSize])
	return out

def movingMedian(x,binSize,slide):
	l=np.floor((len(x)-binSize)/slide)
	out=np.zeros(l);
		
	for i,val in enumerate(out):
		out[i]=np.median(x[i*slide:i*slide+binSize])
	return out

def R2(bestFit,data):
	SSres=np.sum(np.power(data-bestFit,2))
	SStot=np.sum(np.power(data-np.mean(data),2))
	return 1-SSres/SStot
	

def piecewiseLinear(x,Tc,Oc,Ol,delta,genLen):
	Tl=np.mod(Ol+delta,genLen)
	#if (Tl-Ol<0):
	#	a=(Tc-Oc)/-delta
	#else:
	a=(Tc-Oc)/(Tl-Ol)
	x1=min(Tl,Ol)
	x2=max(Tl,Ol)
	y1=Tc if x1==Tl else Oc
	y2=Tc if x2==Tl else Oc
	conds=[x < x1,(x>=x1) & (x<=x2),x>x2]
	funcs=[lambda x: -a*x+y1+a*x1,lambda x: a*x+y1-a*x1,lambda x: -a*x+y2+a*x2]
	return np.piecewise(x,conds,funcs)

def piecewiseLinear2(x,Tc,Oc,Ol,genLen,delta):
	#def periodically_continued(a, b):
	#  	interval = b - a
    #	return lambda f: lambda x: f((x - a) % interval + a)

	#@periodically_continued(-1, 1)
	#def f(x):
    #	return x

	#g = periodically_continued(0, 1)(lambda x: -x)

	def absFun(x,Ol,Oc,Tc):
		a=np.fabs(np.float64(Tc-Oc)/np.float64(genLen*delta))
		return np.abs(-a*x+Oc-Tc)+Tc-(0.5-delta)

	def periodic_function(func, period, offset):
		return func((x - offset) % period)
	#print(str(a))

	#conds=[absFun(x,Ol,Oc,Tc)>Oc,absFun(x,Ol,Oc,Tc)<Tc,(x>=Tc) & (x<=Oc)]
	#funcs=[lambda x: Oc,lambda x: Tc,lambda x: periodic_function(absFun(x,Ol,Oc,Tc),genLen,Ol)]

	#return np.piecewise(x,conds,funcs)
	#return np.abs(-a*x+Oc-Tc)+Oc
	#return -x
	#1 if i<100 else 2 if i>100 else 0
	return periodic_function(lambda x: Tc if absFun(x,Ol,Oc,Tc)<Tc else lambda x: Oc if absFun(x,Ol,Oc,Tc)>Oc else lambda x: absFun(x,Ol,Oc,Tc),genLen,Ol)

def piecewiseLinear3(x,Tc,Tl,Oc,Ol,genLen):
	a=np.fabs(np.float64(Tc-Oc))/np.fabs(np.float64(genLen/2))
	#a=np.fabs(np.float64(Tc-Oc))/np.fabs(np.float64(Tl-Ol))
	#def absFun(x,Ol,Oc,Tl,Tc):
		#a=np.fabs(np.float64(Tc-Oc)/np.fabs(np.float64(Tl-Ol)))
		#a=np.fabs(np.float64(Tc-Oc)/np.float64(genLen/2))
		#return np.abs(-a*x+Oc-Tc)+Tc

#	def periodic_function(func, period, offset):
#		return func((x - offset) % period)

	return np.abs(-a*((x-Ol)%genLen)+Oc-Tc)+Tc

def mfreqz(b,a=1):
	w,h = signal.freqz(b,a)
	h_dB = 20 * np.log10 (abs(h))
	plt.subplot(211)
	plt.plot(w/max(w),h_dB)
	plt.ylim(-150, 5)
	plt.ylabel('Magnitude (db)')
	plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
	plt.title(r'Frequency response')
	plt.subplot(212)
	h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
	plt.plot(w/max(w),h_Phase)
	plt.ylabel('Phase (radians)')
	plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
	plt.title(r'Phase response')
	plt.subplots_adjust(hspace=0.5)
	plt.savefig(str(sys.argv[1])+".mfreqz.png")
	plt.clf()

def rejectOutliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

def pxn(x, C, l):
    return (2**(1 + C - (2*C*x)/float(l))*C*np.log(2))/(float(l)*(-1 + 2**C));

def fit_signal(signal,l):
    x1=np.linspace(0,1,len(signal))
    
    piecewiseModel=Model(piecewise_prob)
    piecewiseModel.set_param_hint('l',     value=1,vary=False)
    piecewiseModel.set_param_hint('C',   vary=True,  value=.1,min=0,max=1)
    piecewiseModel.make_params()
    
    res=piecewiseModel.fit(signal,x=x1,weights=np.sin(1.5*x1)+1.5)
    return res

def smooth_signal(signal,length):
    xi=np.linspace(0,1,len(signal))
    xi2=np.linspace(0,1,length)
    iss=UnivariateSpline(xi,signal,k=3)
    iss.set_smoothing_factor(.00001)
    return iss(xi2)

def sharpen_filter(signal):
    x = np.linspace(-0.00001*np.pi/float(1), 1*np.pi/float(1),10**3)
    return ifft(fft(signal,10**3)/mv_filt(1.27,x))[0:len(signal)]

acc,fileEnding=os.path.splitext(sys.argv[1])

#print acc

try:
	if fileEnding==".depth":
		data=pandas.read_csv(sys.argv[1],delimiter=" ")
		data=np.array(data).astype('int32')
	elif fileEnding==".npy":
		data=np.load(sys.argv[1])
		#acc,fileEnding=os.path.splitext(acc)
	else:
		print("No file found: "+sys.argv[1])
		sys.exit()
except IOError as e:
	print("I/O error "+sys.argv[1])
	sys.exit()

with open(os.path.join(sys.argv[2],acc+".xml")) as fd:
    obj = xmltodict.parse(fd.read())

genomeLen=int(obj['DocSum']['Item'][8]['#text'])
bacteriaName=obj['DocSum']['Item'][1]['#text']
### oriC

try:
	oriData=DataFrame.from_csv(os.path.join(sys.argv[3],"bacteria_record.dat"),  sep='	',index_col=1)
except:
	pass

try:
	if(data.shape[1]==1):
		x = np.arange(1,genomeLen)-1
		yr = data[:, 0].astype(float)
	elif(data.shape[1]==2):
		x = data[:, 0]-1
		y = data[:, 1].astype(float)
		yr = np.zeros(genomeLen)
		yr[x] = y
except:
	genomeLen=data.shape[0]
	x = np.arange(1,genomeLen)-1
	#yr = yr = np.zeros(genomeLen)
	yr = data

###########

#print(sys.argv[1])

#try:
#	print(oriData[repr(sys.argv[1][:-6]])
#except:
#	print('Non existent oric.')

#print(oriData.to_string())
#s=oriData.loc['Organism']
#for column in list(oriData.columns.values):
#print(str(column))

#print(oriData.index)

coveragePercentage=np.sum(yr)/genomeLen

if (coveragePercentage<0.05):
	print(sys.argv[1]+": Coverage too low, exiting.")
	sys.exit()

#doric=np.loadtxt(sys.argv[1],delimiter=" ")
#doric=np.array(data).astype('int32')

#yr=np.zeros(genomeLen)
#yr=yr*np.mean(y)
#yr[x]=y

# Hermansson filter outliers in initial y vector
# exit if too many bins are empty
binLen=2000
bins=binData(yr,binLen)
bins2=binData(yr,binLen)

tmp=0
med=np.median(bins)
std=np.std(bins, dtype=np.float64)
medYr=np.median(yr)

for i,val in enumerate(bins):
	if (val>0):
		tmp=tmp+1

	if (val>2.5*std+med or val<-2.5*std+med):
		yr[i*binLen:i*binLen+binLen]=medYr
		bins2[i]=med


#print(str(3*std+med))

#yr=rejectOutliers(yr)

# exit if too many bins are empty
#bins=binData(yr,10000)

#tmp=0
#for val in bins:
#	if (val>0):
#		tmp=tmp+1

if (float(tmp)/float(np.ceil(len(yr)/binLen))<0.6):
	print(sys.argv[1]+": Bin coverage too low, exiting.")
	sys.exit()

#plt.clf()
#plt.plot(bins,'b')
#plt.savefig(str(sys.argv[1])+".bins.png")

#plt.plot(yr[100000:500000],         'k-')
#plt.savefig(str(sys.argv[1])+".yr.png")
#plt.clf()

#med=np.median(bins)
#bins=[value for value in bins if (value > (1/2)*med) & (value<2*med)]

#med=np.median(yr)
#print(str(med))
#yr=[value for value in yr if (value<10*med)]

y1=yr/coveragePercentage
#med=np.median(y)
#y1=[value for value in yr if (value > (1/4)*med) & (value<4*med)]

#b,a=signal.iirdesign(0.01, 0.02, gstop= 60, gpass=1, ftype='ellip')

#y1=signal.filtfilt( b,a, yr)

y1=movingAverageC(y1.tolist(),10**4,100)

#y1=rejectOutliers(y1,m=4)

med=np.median(y1)
y1=[value for value in y1 if (value > (1./6.)*med) & (value<6*med)]

#med=np.median(y1)
#y1=[value for value in y1 if (value > (1/5)*med) & (value<5*med)]

plt.clf()
plt.plot(bins2,'b')
plt.savefig(str(sys.argv[1])+".bins2.png")

plt.clf()
plt.plot(y1,'b')
plt.savefig(str(sys.argv[1])+".y1.png")

#var=np.var(y1, dtype=np.float64)
#print("\nvar y1: "+str(np.sqrt(var)))

#from scipy.fftpack import fft
# Number of samplepoints
# sample spacing
#N=len(y1)
#T = 1.0
#yf = fft(y1)
#xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
#plt.plot(xf, 2.0/N * np.abs(yf[0:N/2]))
#plt.grid()
#plt.savefig(str(sys.argv[1])+".fft.png")
#plt.clf()

# provide them to firwin
#h=signal.firwin(numtaps=len(y1), cutoff=Fc, nyq=Fs/2)
y1=movingMedianC(y1,10**4,100)

#med=np.median(y1)
#y1=[value for value in y1 if (value > (1./6.)*med) & (value<8.*med)]

#var=np.var(np.log2(y1), dtype=np.float64)
#print("\nvar y1_2: "+str(np.sqrt(var))+"\n")


#b,a=signal.iirdesign(0.01, 0.02, gstop= 60, gpass=1, ftype='ellip')

#y1=signal.filtfilt( b,a, y1)
#mfreqz(b,a)

#y1=movingMedian(y1,7*10**3,100)

#med=np.median(y1)
#y1=[value for value in y1 if (value > (1.25)*med) & (value<.75*med)]

#y1=movingMedian(y1,3*100,2)

#y1=signal.lfilter( h, 1.0, y1)

#med=np.median(y1)
#y1=[value for value in y1 if (value > (1/4)*med) & (value<4*med)]

#print(sys.argv[1]+": var="+str(np.var(np.log2(y1), dtype=np.float64)))

#med=np.median(y1)
#y1=[value for value in y1 if (value > (1/2)*med) & (value<2*med)]

#var=np.var(np.log2(y1), dtype=np.float64)
#med=np.median(np.log2(y1))

#y1=[value for value in y1 if (np.log2(value) < med+2*np.sqrt(var)) & (np.log2(value)>med-2*np.sqrt(var))]

if (len(y1)<25):
	print(sys.argv[1]+": Smoothed length too low ("+str(len(y1))+"), exiting.")
	sys.exit()

#x1=np.linspace(0,2000000-1,1000)

#print(str(piecewiseLinear2(1000,40,100,100,3*10^6,genomeLen)))

#plt.plot(x1,piecewiseLinear2(x1,40,100,1000000,2000000,0.5))

#plt.plot(np.linspace(0,genomeLen,len(bins)),bins/10000,         'bo')
#plt.plot(x1, y1,         'bo')
#plt.savefig(str(sys.argv[1])+".bins.png")
#plt.clf()

y1=smooth_signal(y1,10**3)
y1=np.real(sharpen_filter(y1))
plt.clf()
plt.plot(y1,'b')
plt.savefig("".join([sys.argv[1],"blak"])+".png")


print(sys.argv[1]+": Coverage OK ("+str(coveragePercentage)+" x).")
x1=np.linspace(0,genomeLen-1,len(y1))

### Korem fit

# fit piecewise using least squares (LM)

piecewiseModel=Model(piecewiseLinear)
piecewiseModel.set_param_hint('Ol',     value=0.5*genomeLen,vary=True,min=0,max=genomeLen)
#piecewiseModel.set_param_hint('Tl',     value=0.6*genomeLen,vary=True,min=0,max=genomeLen)
piecewiseModel.set_param_hint('Tc',     value = np.median(np.log2(y1)), vary=True)
piecewiseModel.set_param_hint('Oc',     value = np.median(np.log2(y1)), vary=True)
piecewiseModel.set_param_hint('delta',  value = 0.545*genomeLen, min=0.45*genomeLen, max=0.55*genomeLen, vary=True)#,expr='Tl-Ol if Tl-Ol > 0 else Ol-Tl')

#piecewiseModel.set_param_hint('Tl',     value=0.6*genomeLen,vary=True,min=0,max=genomeLen,expr='mod(Ol+delta,genLen)')
piecewiseModel.set_param_hint('genLen',   vary=False,  value=genomeLen)
piecewiseModel.make_params()

result=piecewiseModel.fit(np.log2(y1),x=x1)
resR2=R2(result.best_fit,y1)
# var=np.var(abs(np.log2(y1)-result.best_fit), dtype=np.float64)
#med=np.median(np.log2(y1))

# y1=[value for value in y1 if (np.log2(value) < med+1.45*np.sqrt(var)) & (np.log2(value)>med-1.45*np.sqrt(var))]

# x1=np.linspace(0,genomeLen-1,len(y1))

# piecewiseModel=Model(piecewiseLinear)
# piecewiseModel.set_param_hint('Ol',     value=0.3*genomeLen,vary=True,min=0,max=genomeLen)
# piecewiseModel.set_param_hint('Tl',     value=0.6*genomeLen,vary=True,min=0,max=genomeLen)
# piecewiseModel.set_param_hint('Tc',     value = np.median(np.log2(y1)), vary=True)
# piecewiseModel.set_param_hint('Oc',     value = np.median(np.log2(y1)), vary=True)
# #piecewiseModel.set_param_hint('delta',  value = 0.5*genomeLen, min=0.45*genomeLen, max=0.55*genomeLen, vary=True,expr='fabs(Tl-Ol)')
# piecewiseModel.set_param_hint('delta',  value = 0.5*genomeLen,expr='fabs(Tl-Ol)')
# piecewiseModel.make_params()

# result=piecewiseModel.fit(np.log2(y1),x=x1)

### fit special

# piecewiseModel=Model(piecewiseLinear3)
# piecewiseModel.set_param_hint('Ol',     value=0.5*genomeLen,vary=True,min=0,max=genomeLen)
# piecewiseModel.set_param_hint('Tl',     value=0.5*genomeLen,vary=True,min=0,max=genomeLen)
# #expr='fmod(Tl-delta*genLen/2, genLen)')
# piecewiseModel.set_param_hint('Tc',     value = np.median(np.log2(y1)), vary=True)
# piecewiseModel.set_param_hint('Oc',     value = np.median(np.log2(y1)), vary=True)
# #piecewiseModel.set_param_hint('delta',  value = 0.5, min=0.45, max=0.55, vary=True)
# piecewiseModel.set_param_hint('genLen',  value = genomeLen)
# piecewiseModel.make_params()

# result=piecewiseModel.fit(np.log2(y1),x=x1)

bestVal=result.best_values
bestVal['Tl']=np.mod(bestVal['Ol']+bestVal['delta'],bestVal['genLen'])

if(bestVal['Tc']>bestVal['Oc']):
	bestVal['Tc'],bestVal['Oc'] = bestVal['Oc'],bestVal['Tc']
	tmpT=bestVal['Ol']
	bestVal['Ol']=np.mod(bestVal['Ol']+bestVal['delta'],bestVal['genLen'])
	bestVal['Tl']=tmpT
	#bestVal['Tl'],bestVal['Ol'] = bestVal['Ol'],bestVal['Tl']

PTR=2**bestVal['Oc']/2**bestVal['Tc']

if(PTR<1.1):
	print(sys.argv[1]+": PTR < 1.1 ("+str(PTR)+"), exiting.")
	sys.exit()

#### p-value test
counter=0
nIter=2000
cCoeff=np.power(np.corrcoef(result.best_fit,y1),2)
for i in range(1,nIter):
	y1Tmp=np.random.permutation(y1)
	#print("\n"+np.array_str(y1Tmp))
	resTmp=piecewiseModel.fit(np.log2(y1Tmp),x=x1)
	#bvTmp=resTmp.best_values
	
	#bvTmp['Tl']=np.mod(bvTmp['Ol']+bvTmp['delta'],bvTmp['genLen'])
	#if(bvTmp['Tc']>bvTmp['Oc']):
	#	bvTmp['Tc'],bvTmp['Oc'] = bvTmp['Oc'],bvTmp['Tc']
	#	tmpT=bvTmp['Ol']
#		bvTmp['Ol']=np.mod(bvTmp['Ol']+bvTmp['delta'],bvTmp['genLen'])
#		bvTmp['Tl']=tmpT


	#if (abs((bvTmp['Oc']-bvTmp['Tc'])/(bvTmp['Ol']-bvTmp['Tl']))>=abs((bestVal['Oc']-bestVal['Tc'])/(bestVal['Ol']-bestVal['Tl']))):
	cCoeffTmp=np.power(np.corrcoef(resTmp.best_fit,y1),2)
#	print("\n"+str(R2(resTmp.best_fit,y1)))
#	if(R2(resTmp.best_fit,y1)>=resR2):
#	print("\n"+str(cCoeffTmp[0,1]))
	if(cCoeffTmp[0,1]>=cCoeff[0,1]):
		counter+=1

#print("\nres: "+str(resR2))
#print("\nres: "+str(cCoeff[0,1]))
print("\ncntr: "+str(counter))

if (float(counter)/float(nIter)>0.03):
	print("p-value: "+str(float(counter)/float(nIter))+" > 0.03, exiting.")
	sys.exit()

####

np.save(sys.argv[1],np.log2(y1))
np.save(sys.argv[1]+".best",np.array([bestVal['Ol'],bestVal['Tl'],bestVal['Oc'],bestVal['Tc']]))

f = open(sys.argv[1]+'.log', 'w')

f.write(result.fit_report())

f.write("\n\nName: "+bacteriaName)
f.write("\n\nFit success: "+str(result.success))

f.write("\nPTR:\t"+str(PTR))
f.write("\nMedian:\t"+str(np.median(y1)))
for key,val in bestVal.items():
	f.write("\n"+key+":\t"+str(val))

#f.write("\n"+result.message)

f.write("\n####\n\n")

#ci = conf_interval(result, sigmas=[0.68,0.9], trace=False, verbose=False)
#f.write("\nConf. interval\t.95\t.68\t0\t.68\t.95")
#for key,val in ci.items():
#	f.write("\n"+key+":\t")
#	for val2 in val:
#		f.write(str(val2[1])+"\t")
f.close()

#x2=np.linspace(0,1000000,100)

#plt.plot(x2,piecewiseLinear3(x2,40,10,100,10**6/4,1000000))

#plt.plot(np.linspace(0,genomeLen,len(bins)),bins/10000,         'bo')
#plt.plot(x1, y1,         'bo')
#plt.ylim([0,max(y1)])
#plt.savefig(str(sys.argv[1])+".bins.png")
#plt.clf()
plt.clf()

result.plot()

#plt.plot(x1, np.log2(y1),         'bo')
try:
	oriLoc=oriData.loc[acc,'oriC location']
	oriLoc=oriLoc.split('..',1)
	oriLoc=int(oriLoc[0])
	med=np.median(result.best_fit)
	plt.plot(oriLoc, med,         'g*')
	print(oriData.loc[acc,'Organism'])
except:
	pass
#plt.plot(x1, result.init_fit, 'k--')
#plt.plot(x1, result.best_fit, 'r-')


plt.savefig(str(sys.argv[1])+".png")
plt.clf()
plt.plot(x1,y1,'b')
plt.savefig(str(sys.argv[1])+".coverage.png")
plt.clf()
plt.plot(bins,'b')
plt.savefig(str(sys.argv[1])+".bins.png")

