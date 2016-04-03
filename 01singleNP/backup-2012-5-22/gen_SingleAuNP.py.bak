import math
from math import cos, sin, sqrt
from numpy import matrix
from random import randint
from numpy import array
from numpy import cross
from numpy import dot
from numpy import abs
import re
import string

import sys
sys.path.append('./') 
import PDB_module

def pointsOnSphere(N):
    N = float(N) # in case we got an int which we surely got
    pts = []
 
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / N
    for k in range(0, int(N)):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r, y, math.sin(phi)*r])
 
    return pts

def getX(x,scale,shift=0):
    return x*scale+shift

#number of atom of a big ball
N=40
#number of mao per big ball
L=20
#number of atom per mao
lenspacer=5
lenlinker=3

#radius of big ball
s = 8.0
ssp =1.0
sfl =1.0
mol = []
bnd = []
ang = []
ndx = 0
#generate a big ball
for n,pt in enumerate(pointsOnSphere(N)):
	ndx = ndx+1
	#mol formate (ndx, numform-atom-type, string-form-of-atom-type, residue-type, x, y z)
	mol.append([ndx, 1, 'NP', 'DNP', getX(pt[0],s), getX(pt[1],s), getX(pt[2],s)])

s = s+ssp
roots = pointsOnSphere(L)
nroots = []
for n,pt in enumerate(roots):
	nroots.append(ndx)
	i=0
	for i in range(lenspacer):
		ndx = ndx+1
		mol.append([ndx, 2, 'SP', 'DNP', getX(pt[0],s+i*ssp), getX(pt[1],s+i*ssp), getX(pt[2],s+i*ssp)])
		
		#GET bonds, angles 
		if (i==0):
			p1=array([pt[0], pt[1], pt[2]])
			pointer=0
			pCos=-1
			for m in range(N):
				p2=array([mol[m][4], mol[m][5], mol[m][6]])
				if(dot(p1, p2)> pCos):
					pCos=dot(p1, p2)
					pointer=m

			bnd.append([pointer+1, ndx])
		elif (i>0):
			bnd.append([ndx-1, ndx])
			if(i>1):
				ang.append([ndx-2, ndx-1, ndx])
		
				

	#generate FL vector and FLs
	n1=array([pt[0], pt[1], pt[2]])
	#print n1
	while (True):
		Randint=randint(0, L-1)
		n2=array([roots[Randint][0], roots[Randint][1], roots[Randint][2]])
		#print n2
		if (abs(dot(n1, n2))<0.9):
			n3=cross(n1, n2)
			#print n3
			break
	nf=cross(n1, n3)	
	#print n1, n2, n3, nf
	for j in range(lenlinker):
		ndx = ndx+1
		mol.append([ndx, 3, 'C1', 'DNP', getX(pt[0],s+(i+j+1)*ssp), getX(pt[1],s+(i+j+1)*ssp), getX(pt[2],s+(i+j+1)*ssp)])
		
		#Get bonds, angles
		if (j==0):
			if (lenspacer>0):
				bnd.append([ndx-1, ndx])
			if (lenspacer>=2):
				ang.append([ndx-2, ndx-1, ndx])
		elif (j>0):
			bnd.append([ndx-4, ndx])
			if (j==1):
				if (lenspacer>=1):
					ang.append([ndx-5, ndx-4, ndx])
			elif (j>=2):
				ang.append([ndx-8, ndx-4, ndx])	
		
		#generate ct2 
		ndx = ndx+1
		#AN='%s%d'%('C',j+1)
		ct2=list([getX(n3[0], sfl, getX(pt[0],s+(i+j+1)*ssp)), \
                         getX(n3[1], sfl, getX(pt[1],s+(i+j+1)*ssp)), \
                         getX(n3[2], sfl, getX(pt[2],s+(i+j+1)*ssp))])
		mol.append([ndx, 3+j+2, 'C2', 'DNP', ct2[0], ct2[1], ct2[2] ])
		
		#Get bonds, angles
		bnd.append([ndx-1, ndx])
		
		#generate FLs
		#FL-1
		ndx = ndx+1
		mol.append([ndx, 4, 'FL', 'DNP', getX(nf[0], sfl, ct2[0]), \
                                              getX(nf[1], sfl, ct2[1]), \
                                              getX(nf[2], sfl, ct2[2])])
		
		#Get bonds, angles
		bnd.append([ndx-1, ndx])

		#FL-2
		ndx = ndx+1
                mol.append([ndx, 4, 'FL', 'DNP', getX(nf[0]*(-1.0), sfl, ct2[0]), \
                                              getX(nf[1]*(-1.0), sfl, ct2[1]), \
                                              getX(nf[2]*(-1.0), sfl, ct2[2])])
		
		#Get bonds, angles
		bnd.append([ndx-2, ndx])
		ang.append([ndx-1, ndx-2, ndx])

		#the ct2 at the end of a branch
		if (j==(lenlinker-1)):
			ndx = ndx+1
			mol.append([ndx, 4, 'FL', 'DNP', getX(n1[0], sfl, ct2[0]),\
                                                        getX(n1[1], sfl, ct2[1]), \
                                                        getX(n1[2], sfl, ct2[2])])
		#Get bonds, angle for ct2-ct2-ct2
		#bond
		if (j>=1 and j<(lenlinker-1)):
			bnd.append([ndx-6, ndx-2]) 
		elif (j==(lenlinker-1)):
			bnd.append([ndx-3, ndx])
			bnd.append([ndx-7, ndx-3])
		#angle
		if ( j>1 and j< (lenlinker-1) ):
			ang.append([ndx-10, ndx-6, ndx-2])
		elif (j==(lenlinker-1)):
			ang.append([ndx-11, ndx-7, ndx-3])
			ang.append([ndx-7, ndx-3, ndx])

#PRINT PDB
filename="DNP.pdb"
try:
        fp = open(filename , 'w')
except:
        print "Error: No such file: "+filename
        exit(1)
for j,p in enumerate(mol):

	atom = PDB_module.Atom_class(AtomSerial=p[0], AtomName=p[1], ResidueName=p[3], ResidueSerial=1, AtomCoorX=p[4], \
              AtomCoorY=p[5], AtomCoorZ=p[6])	
	fp.write("%s\n" %(atom.atom_2_PDBformat()))
fp.close()



#PSF

#bond
bondtype = []
#//bondList: BondIndex, BondTypeIndex, BondAtomIndex1, BondAtomIndex2, BondAtomIndex1Type, BondAtomIndex2Type
bondList = []

bondtypetmp = []
bondtypeindex = 0
for j, p in enumerate(bnd):
	if (j==0):
		bondtypeindex =bondtypeindex +1
		bondtypetmp = [mol[p[0]-1][1], mol[p[1]-1][1]]
		bondtypetmp = sorted (bondtypetmp)
		bondtype.append([bondtypeindex, bondtypetmp[0], bondtypetmp[1]])
		#print bondtype[0][1], bondtype[0][2]
		bondList.append([j+1, bondtypeindex, p[0], p[1], bondtypetmp[0], bondtypetmp[1]])
	else:
		bondtypetmp=[mol[p[0]-1][1], mol[p[1]-1][1]]
		bondtypetmp = sorted (bondtypetmp)
		#print bondtypetmp
		pointer=True
		for b, bt in enumerate(bondtype):
			if([bt[1], bt[2]] == [bondtypetmp[0], bondtypetmp[1]]):
				bondList.append([j+1, bt[0], p[0], p[1], bt[1], bt[2]])
				break
			#elif(([bt[1], bt[2]] != [bondtypetmp[0], bondtypetmp[1]]) and \
			elif(b == (len(bondtype)-1)):
				pointer=False	
				bondtypeindex =bondtypeindex +1
				bondtype.append([bondtypeindex, bondtypetmp[0], bondtypetmp[1]])
				bondList.append([j+1, bondtypeindex, p[0], p[1], bondtypetmp[0], bondtypetmp[1]])
				#print bondtypeindex, bondtypetmp[0], bondtypetmp[1]
				break
for j, p in enumerate(bondList):
	continue
	#print j, p, bnd[j]

print "========BOND TYPE========"
if(len(bnd) != len(bondList)):
	print "Error in Bond Type Statistic"
	print "Num of bnd is %d" %(len(bnd))
	print "Num of bondList is %d" %(len(bondList))
	exit(1)
for j, p in enumerate(bondtype):
	print p[0], p[1], p[2]
				
#bond
angleType = []
#AngleList: AngleIndex, AngleTypeIndex, AngleAtomIndex1, AngleAtomIndex2, AngleAtomIndex3, 
#                                       AngleAtomIndex1Type, AngleAtomIndex2Type, AngleAtomIndex3Type,
angleList = []
angletypetmp = []
angletypetmpEnd = []
angletypeindex = 0
for j, p in enumerate(ang):
	if (j==0):
		angletypeindex = angletypeindex +1
		angletypetmp = [mol[p[0]-1][1], mol[p[1]-1][1], mol[p[2]-1][1]]
		
		#check Angle's End atoms
		#the End of an angle is equal, i.e., (3 5 7) == (7, 5, 3)
		angletypetmpEnd =[angletypetmp[0], angletypetmp[2]]
		angletypetmpEnd = sorted(angletypetmpEnd)
		angletypetmp = [angletypetmpEnd[0], angletypetmp[1], angletypetmpEnd[1]]
		#

		angleType.append([angletypeindex, angletypetmp[0], angletypetmp[1], angletypetmp[2]])
		angleList.append([j+1, angletypeindex, p[0], p[1], p[2], angletypetmp[0], angletypetmp[1], angletypetmp[2]])
	else:
		angletypetmp = [mol[p[0]-1][1], mol[p[1]-1][1], mol[p[2]-1][1]]
               
		#Check Angle's End atoms 
		angletypetmpEnd =[angletypetmp[0], angletypetmp[2]]
                angletypetmpEnd = sorted(angletypetmpEnd)
                angletypetmp = [angletypetmpEnd[0], angletypetmp[1], angletypetmpEnd[1]]
		#

		pointer = True
		for an, ant in enumerate(angleType):
			if([ant[1], ant[2], ant[3]] == [angletypetmp[0], angletypetmp[1], angletypetmp[2]]):
				angleList.append([j+1, ant[0], p[0], p[1], p[2], ant[1], ant[2], ant[3]])
				break
			elif(([ant[1], ant[2], ant[3]] != [angletypetmp[0], angletypetmp[1], angletypetmp[2]]) and \
                              (an == len(angleType)-1)):
				pointer = False
				angletypeindex = angletypeindex +1
				angleType.append([angletypeindex, angletypetmp[0], angletypetmp[1], angletypetmp[2]])
				angleList.append([j+1, angletypeindex, p[0], p[1], p[2], \
                                                  angletypetmp[0], angletypetmp[1], angletypetmp[2]])
				break
		
print "=======ANGLE TYPE========"
if(len(bnd) != len(bondList)):
        print "Error in Angle Type Statistic"
        print "Num of angle is %d" %(len(ang))
        print "Num of anlgeList is %d" %(len(angleList))
        exit(1)

for j, p in enumerate(angleType):
	print p[0], p[1], p[2], p[3]





#PRINIT PSF
filename="DNP.psf"
try:
        fp = open(filename , 'w')
except:
        print "Error: No such file: "+filename
        exit(1)

#tile
fp.write("RESIDNAME  DNP\n")
fp.write("%s %3d\n" %('NUMATOM', len(mol)))
fp.write("%s %3d\n" %('NUMBOND', len(bnd)))
fp.write("%s %3d\n" %('NUMANGLE', len(ang)))
fp.write("%s %3d\n" %('NUMIMPR', 0))

line=''

#ATOM part
for j,p in enumerate(mol):
        line = PDB_module.Atom2PSF(atomHead='ATOM', atomSerial=p[0], atomName=p[1], atomType=p[1], \
                        Mass=1.0, Charge=0.0, Unset=0.0)
        fp.write("%s\n" %(line))

#BOND part
for j,p in enumerate(bondList):
	line = PDB_module.Bond2PSF(bondHead='BOND', bondSerial=p[0], bondTypeSerial=p[1], bondIndex1=p[2], bondIndex2=p[3], \
                                  bondIndex1Type=p[4], bondIndex2Type=p[5])
	fp.write("%s\n" %(line))

#ANGLE part 
for j,p in enumerate(angleList):
	line = PDB_module.Angle2PSF(angleHead='ANGLE', angleSerial=p[0], angleTypeSerial=p[1], \
                                    angleIndex1=p[2], angleIndex2=p[3], angleIndex3=p[4], \
                                    angleIndex1Type=p[5], angleIndex2Type=p[6], angleIndex3Type=p[7])
	fp.write("%s\n" %(line))

fp.close()

for n,pt in enumerate(roots):
	for i in range(1,lenspacer):
		ndx = ndx+1
		mol.append([2, getX(pt[0],s+i), getX(pt[1],s+i), getX(pt[2],s+i)])
		if (i==1):
			bl=nroots[n]
		else:
			bl=ndx-1
		bnd.append([bl,ndx])

nmol=8
type = ['C','O']
center = [ 
[0,0,0],
[0,0,1],
[0,1,0],
[0,1,1],
[1,0,0],
[1,0,1],
[1,1,0],
[1,1,1],
]
#print len(mol)
for i in range(0,nmol):
	x=center[i][0]*48
	y=center[i][1]*48
	z=center[i][2]*48
	ndx0 = int(len(mol))*i+1
	#for j,p in enumerate(mol):
#print "ATOM%7d  %c   SPH A   1    %8.3f%8.3f%8.3f  1.00  0.00" % (j+ndx0, type[p[0]-1], p[1],p[2],p[3])
		
		#print "%d %d %d %8.3f%8.3f%8.3f" % (j+ndx0, i+1, p[0], p[1]+x,p[2]+y,p[3]+z) 

nb = 0
for i in range(0, nmol):
	ndx0 = int(len(mol))*i
	for n,b in enumerate(bnd):
		nb = nb+1
		#print "%d %d %d %d" % (nb, 1, b[0]+ndx0, b[1]+ndx0)

