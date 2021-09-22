import numpy as np
import sys

rp=float(sys.argv[1])
rnp=float(sys.argv[2])
aeff=0.50*(rp**2+rnp**2)
rtot=rp+rnp

M2P=np.zeros((6,6))
R2P=np.zeros((6,6))
rij=np.zeros(3)
udr=np.zeros(3)
rmat=np.zeros((3,3))
umat=np.identity(3)

for i in range(3):
        M2P[i,i]=1.0/rp
        M2P[i+3,i+3]=1.0/rnp

dd=0.1

n=9000; m=14; m1=int(m/2)
A1=np.zeros((n,m))
B1=np.zeros(n)
A2=np.zeros((n,m))
B2=np.zeros(n)
k=0
for i in range(100):
	for j in range(100):
		rij = np.asarray([float(i)*dd,float(j)*dd,0.0])
		dr1 = np.sqrt(np.dot(rij,rij))
		if(dr1 >= rtot):
			idr1=1.0/dr1; idr3=idr1**3
			udr=rij*idr1
			rmat=np.outer(udr,udr)
			M2P[0:3,3:6] = 0.75*idr1*(umat+rmat)+0.5*aeff*idr3*(umat-3.0*rmat)
			M2P[3:6,0:3] = M2P[0:3,3:6]
			R2P=np.linalg.inv(M2P)
			for p in range(m1):
				A1[k,p] = idr1**p
				A2[k,p] = A1[k,p]
			for p in range(m1,m):
				A1[k,p] = rmat[0,0] * idr1**(p-m1)
				A2[k,p] = A1[k,p]
			B1[k] = R2P[0,0]
			B2[k] = R2P[0,3]
			k+=1
		if(k>=n):
			break
	else:
		continue
	break

a,b,c,d=np.linalg.lstsq(A1,B1,rcond=None)
f=open('r11_coeff.txt','w')
for i in range(m):
	f.write('%20.8f \n'%a[i])
f.close()

a,b,c,d=np.linalg.lstsq(A2,B2,rcond=None)
f=open('r12_coeff.txt','w')
for i in range(m):
	f.write('%20.8f \n'%a[i])
f.close()
