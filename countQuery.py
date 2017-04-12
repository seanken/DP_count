######################################
##Implement both naive exponential and our version 
######################################
import math;
from numpy.random import multinomial as mult;
import matplotlib.pyplot as plt;
import time;
##
##Implements the naive exponential mechanism
##Return -1 if parameters break rules
##
def NaiveExp(y,epsilon,betap,betan,alphap,alphan,rmin,rmax,getZ=False,getOmega=False):
	if alphap>1 or alphan>1 or alphap<0 or alphan<0:
		return -1;
	if betap<0 or betan<0:
		return -1;
	if rmax<rmin:
		return -1;

	scores=[0.0 for i in range(rmin,rmax+1)];
	chang=max(betan,betap);	


	for i in range(rmin,rmax+1):
		if i<y:
			scores[i-rmin]=betan*(y-i)**alphan;
		else:
			scores[i-rmax]=betap*(i-y)**alphap;


	mx=min(scores);
	scores=[s-mx for s in scores];


	scores=[epsilon*s/(2*chang) for s in scores]

	scores=[math.exp(-s) for s in scores];

	Z=sum(scores);

	if getZ:
		return Z;

	scores=[s/Z for s in scores];


	pick=mult(1,scores);

	i=min([i for i in range(0,len(pick)) if pick[i]>0])+rmin;
	if getOmega:
		return epsilon/(2.0*chang)
	return i;


##
##Calculates the real epsilon for a given omega
##
def RealEps(omega,betap,betan,alphap,alphan,rmin,rmax):
	epsilon=omega*2*max(betap,betan);
	Z1=NaiveExp(rmin,epsilon,betap,betan,alphap,alphan,rmin,rmax,getZ=True);

	Z2=NaiveExp(rmin+1,epsilon,betap,betan,alphap,alphan,rmin,rmax,getZ=True);
	Z3=NaiveExp(rmax,epsilon,betap,betan,alphap,alphan,rmin,rmax,getZ=True);
	Z4=NaiveExp(rmax-1,epsilon,betap,betan,alphap,alphan,rmin,rmax,getZ=True);

	a=math.log(Z4/Z3)+betap*omega;
	b=math.log(Z2/Z1)+betan*omega;
	return max(a,b);


##
##calulcates the risk for a given omega
##
def Risk(y,omega,betap,betan,alphap,alphan,rmin,rmax):
	scores=[0.0 for i in range(rmin,rmax+1)];
	chang=max(betan,betap);	

	for i in range(rmin,rmax+1):
		if i<y:
			scores[i-rmin]=betan*(y-i)**alphan;
		else:
			scores[i-rmax]=betap*(i-y)**alphap;

	mx=min(scores);
	scores2=[s-mx for s in scores];

	scores2=[omega*s for s in scores2]

	scores2=[math.exp(-s) for s in scores2];

	Z=sum(scores2);


	scores2=[s/Z for s in scores2];
	
	return sum([scores[i]*scores2[i] for i in range(0,len(scores))]);

##
##Implements our version
##
def SmartExp(y,epsilon,betap,betan,alphap,alphan,rmin,rmax,k=10,getOmega=False):

	chang=max(betap,betan);
	poss=[(1.0+i/float(k))*epsilon/(2.0*chang) for i in range(0,k+1)]
	
	realeps=[2*epsilon for i in range(0,k+1)];
	a=0;
	b=k;
	while b-a>1:
		c=int(math.floor((a+b)/2.0));
		realeps[c]=RealEps(poss[c],betap,betan,alphap,alphan,rmin,rmax) 
		if realeps[c]>epsilon:
			b=c;
		else:
			a=c;

	omega=poss[a];

	return NaiveExp(y,omega*2*chang,betap,betan,alphap,alphan,rmin,rmax,getOmega=getOmega)


##
##Plots Risk vs epsilon for different ks and naive method. If k==1 use naive
##
def RiskVsEpsilon(alphan,alphap,betap,betan,rmin,rmax,y,eps=[.1*l for l in range(1,21)],ks=[1,10,100]):
	graphs=[[0.0 for e in eps] for k in ks];
	for i in range(0,len(ks)):
		for j in range(0,len(eps)):
			e=eps[j];
			k=ks[i];
			if k==0:
				se=NaiveExp(y,e,betap,betan,alphap,alphan,rmin,rmax,getOmega=True);
			else:
				se=SmartExp(y,e,betap,betan,alphap,alphan,rmin,rmax,k=k,getOmega=True);
			graphs[i][j]=Risk(y,se,betap,betan,alphap,alphan,rmin,rmax);

	fig=plt.figure();
	ax=fig.add_subplot(111)
	for g in graphs:
		ax.plot(eps,g);
	ax.spines['top'].set_visible(False);
	ax.spines['right'].set_visible(False);
	ax.xaxis.set_ticks_position('none');

	ax.yaxis.set_ticks_position('none');

	plt.xlabel("Epsilon");
	plt.ylabel("Risk");

	plt.yscale("log");
	fig.savefig("pics/Risk.png")
	print "Ratio:"
	#print graph[0][-1]/float(graph[2][-1]);


##
##Risk versus y
##
def RiskVsY(alphan,alphap,betap,betan,rmin,rmax,epsilon,ys=[10**i for i in range(1,8)],ks=[0,10,100]):
	graphs=[[0.0 for y in ys] for k in ks];
	e=epsilon;
	for i in range(0,len(ks)):
		for j in range(0,len(ys)):
			y=ys[j];
			k=ks[i];
			print y;
			if k==0:
				se=NaiveExp(y,e,betap,betan,alphap,alphan,rmin,rmax,getOmega=True);
			else:
				se=SmartExp(y,e,betap,betan,alphap,alphan,rmin,rmax,k=k,getOmega=True);
			graphs[i][j]=Risk(y,se,betap,betan,alphap,alphan,rmin,rmax);

	fig=plt.figure();
	ax=fig.add_subplot(111)
	for g in graphs:
		ax.plot([i for i in ys],[math.log(i) for i in g]);
	ax.spines['top'].set_visible(False);
	ax.spines['right'].set_visible(False);
	ax.xaxis.set_ticks_position('none');

	ax.yaxis.set_ticks_position('none');

	plt.xlabel("Log of Ymax");
	plt.ylabel("Log of Risk");

	plt.yscale("log");
	fig.savefig("pics/RiskVsY.png")


##
##Table used for figure 1a-d
##
def makeTable():
    betan=1
    rmin=3
    rmax=1000000
    res=[]
    for e in [.1*float(i) for i in range(1,21)]:
        print(e)
        for k in [0,10,100]:
            for betap in [1,2,3,4,5]:
                for alpha in [.1*float(j) for j in range(1,11)]:
                    for y in 50*range(1,11):
                        print(y)
                        if k==0:
                            se=NaiveExp(y,e,betap,betan,alpha,alpha,rmin,rmax,getOmega=True);
                        else:
                            se=SmartExp(y,e,betap,betan,alpha,alpha,rmin,rmax,k=k,getOmega=True);
                        risk=Risk(y,se,betap,betan,alpha,alpha,rmin,rmax);
                        res.append([e,k,betap,alpha,y])
    fil=open("table.privcount.txt","w")
    for r in res:
        toSave=" ".join(map(str,r))
        toSave=toSave+"\n"
        fil.write(toSave)
    fil.close()


##
##Runtime in seconds versys Ymax
##
def RuntimeVsYmax(alphan,alphap,betap,betan,rmin,y,epsilon,rMaxs=[10**i for i in range(1,8)],ks=[0,10,100]):
	graphs=[[0.0 for ymax in rMaxs] for k in ks];
	e=epsilon;
	for i in range(0,len(ks)):
		for j in range(0,len(rMaxs)):
			rmax=rMaxs[j];
			k=ks[i];
			
			start=time.clock();
			if k==0:
				se=NaiveExp(y,e,betap,betan,alphap,alphan,rmin,rmax);
			else:
				se=SmartExp(y,e,betap,betan,alphap,alphan,rmin,rmax,k=k);
			end=time.clock();
			graphs[i][j]=end-start;

	fig=plt.figure();
	ax=fig.add_subplot(111)
	for g in graphs:
		ax.plot([i for i in rMaxs],[i for i in g]);
	ax.spines['top'].set_visible(False);
	ax.spines['right'].set_visible(False);
	ax.xaxis.set_ticks_position('none');

	ax.yaxis.set_ticks_position('none');

	plt.xlabel(r"$r_{max}$");
	plt.ylabel("Runtime in Seconds");

	plt.xscale("log");
	plt.yscale("log");
	fig.savefig("pics/Runtime.png")


##
##Mu versus eps
##
def muVsEps(alphan,alphap,betap,betan,rmin,rmax,showIt=False,eps=[.1*i for i in range(1,21)],ks=[0,10,100]):
	graphs=[[0.0 for e in eps] for k in ks];
	#for i in range(0,len(ks)):
	for j in range(0,len(eps)):
		for i in range(0,len(ks)):
			e=eps[j];
			k=ks[i];
		
			start=time.clock();
			if k==0:
				se=NaiveExp(y,e,betap,betan,alphap,alphan,rmin,rmax,getOmega=True);
			else:
				se=SmartExp(y,e,betap,betan,alphap,alphan,rmin,rmax,k=k,getOmega=True);
			print se
			end=time.clock();
			graphs[i][j]=se;
			if graphs[i][j]/graphs[0][j]>2:
				print "ouch"
		print "\n\n"
	fig=plt.figure();
	ax=fig.add_subplot(111)
	for g in graphs:
		ax.plot(eps,g);
	ax.spines['top'].set_visible(False);
	ax.spines['right'].set_visible(False);
	ax.xaxis.set_ticks_position('none');

	ax.yaxis.set_ticks_position('none');

	plt.xlabel(r"$\epsilon$");
	plt.ylabel(r"$\omega$");
	if showIt:
		fig.show();
		return;
	fig.savefig("pics/MuVsEps.png")


if __name__=="__main__":
        makeTable()
"""
    print("hi!")
	Fig1Only=True;
	y=100;
	epsilon=1.0;
	alphan=1.0;
	alphap=1.0;
	betap=2.0;
	betan=1.0;
	rmin=3;
	rmax=10**6;

	print "Fig 1"
	#muVsEps(alphan,alphap,betap,betan,rmin,rmax,showIt=False,eps=[.1*i for i in range(1,20)],ks=[0,10,100])
	#print "Fig 3"
	#RiskVsY(alphan,alphap,betap,betan,rmin,rmax,epsilon);
	if Fig1Only:
		print "Fig 2"
		RuntimeVsYmax(alphan,alphap,betap,betan,rmin,y,epsilon,rMaxs=[10**i for i in range(1,8)],ks=[0,10,100]);
		print "Fig 4"
		#RiskVsEpsilon(alphan,alphap,betap,betan,rmin,rmax,y,eps=[.1*l for l in range(1,21)],ks=[0,10,100]);

	print "Naive Version:"
	start=time.clock();
	ne=NaiveExp(y,epsilon,betap,betan,alphap,alphan,rmin,rmax,getOmega=True);
	end=time.clock();
	print str(end-start);
	print Risk(ne,betap,betan,alphap,alphan,rmin,rmax);
	print ne;
	for k in [10,20,30,100]:
		print "\n";
		print "For k="+str(k);
		
		start=time.clock();
		se=SmartExp(y,epsilon,betap,betan,alphap,alphan,rmin,rmax,k=k,getOmega=True);
		end=time.clock();
		print str(end-start);
		print Risk(se,betap,betan,alphap,alphan,rmin,rmax);
		print se;
	"""
