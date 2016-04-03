import math
 
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
N=960
#number of mao per big ball
L=20
#number of atom per mao
lenb=16

#radius of big ball
s = 8.0
mol = []
bnd = []
ndx = 0
#generate a big ball
for n,pt in enumerate(pointsOnSphere(N)):
	ndx = ndx+1
	mol.append([1, getX(pt[0],s), getX(pt[1],s), getX(pt[2],s)])

s = s+1.0
roots = pointsOnSphere(L)
nroots = []
for n,pt in enumerate(roots):
	ndx = ndx+1
	nroots.append(ndx)
	i=0
	mol.append([1, getX(pt[0],s+i), getX(pt[1],s+i), getX(pt[2],s+i)])

for n,pt in enumerate(roots):
	for i in range(1,lenb):
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
	for j,p in enumerate(mol):
#print "ATOM%7d  %c   SPH A   1    %8.3f%8.3f%8.3f  1.00  0.00" % (j+ndx0, type[p[0]-1], p[1],p[2],p[3])
		print "%d %d %d %8.3f%8.3f%8.3f" % (j+ndx0, i+1, p[0], p[1]+x,p[2]+y,p[3]+z) 

nb = 0
for i in range(0, nmol):
	ndx0 = int(len(mol))*i
	for n,b in enumerate(bnd):
		nb = nb+1
		print "%d %d %d %d" % (nb, 1, b[0]+ndx0, b[1]+ndx0)

