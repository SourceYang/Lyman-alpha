#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:01:36 2022

@author: Emily
"""
import py21cmfast as p21c
import numpy as np
#import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import scipy

def contS(box,x,y,size):
    right=x+1
    left=x-1
    up=y+1
    down=y-1
    if x == size-1:
        right = 0
    if x ==0:
        left = size-1
    if y == size-1:
        up = 0
    if y ==0:
        down = size-1
    if x == size:
        x=0
        left=size-1
        right = 1
       
    if y == size:
        y=0
        down=size-1
        up=1
    fraca = box[right][y]
    fracb = box[left][y]
    fracc = box[x][up]
    fracd = box[x][down]
    if fraca<.5 or fracb<.5 or fracc<.5 or fracd<.5:
        return True
    else:
        return False
def computeDistance(x,y,z,x2,y2,z2,size):
    if z==size-1:
        if z2==0:
            z2=size
    return np.sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)

def TwoDContours(box1):
    size=len(box1[0])
    deltax=[0,1,2,3,4]
    deltay=[0,1,-1,2,-2,3,-3,4,-4]
    zlists2 = []
    xlist2=[]
    ylist2=[]
    xylists2 = []
    countedz= []
    pairslist =[]
    countz=0
    xylists3 =[]
    for z in range(0,size,1):
        slc = np.take(box1, z, axis=2)
        smoothed = gaussian_filter(slc.T,2.5, mode='wrap')  
        smoothed=smoothed.T
        countz=0
        for y in range(0,size,1):
            for x in range(0,size,1):
                 frac1 = smoothed[x][y]
                 check = 0
                 count=0
                 if frac1>0.5 and check<.5:
                     proc = contS(smoothed,x,y,size)
                     if proc==True:
                         if (count<2):
                             for a in range(0,2,1):
                                 deltaX = deltax[a]
                                 for b in range(0,3,1):
                                     deltaY= deltay[b]
                                     if deltaX+np.abs(deltaY)<5 and deltaX+np.abs(deltaY)>0 and x+deltaX<101 and y+deltaY<101:
                                         xcheck=x+deltaX
                                         ycheck = y+deltaY
                                         if x+deltaX==size:
                                             xcheck= 0
                                         if y+deltaY==size:
                                             ycheck= 0
                                         frac2 = smoothed[xcheck][ycheck]
                                         if frac2>.5:
                                             ad = contS(smoothed,x+deltaX,y+deltaY,size)
                                             if ad==True:
                                                 xylists2.append((x,y,z,x+deltaX,y+deltaY))
                                                 xylists3.append((x,y,z,x+deltaX,y+deltaY))
                                                 xlist2.append(x)
                                                 xlist2.append(x+deltaX)
                                                 xlist2.append(None)
                                                 ylist2.append(y)
                                                 ylist2.append(y+deltaY)
                                                 ylist2.append(None)
                                                 zlists2.append(z)
                                                 zlists2.append(z)
                                                 zlists2.append(None)
                                                 count +=1
                                                 countz+=1
        countedz.append(countz)
        pairslist.append(xylists2)
        xylists2 =[]
    return pairslist

def ThreeDlineslist(pairslist,size):
    complist=[]
    lines =[]
    for t in range(0,size,1): #z index
        linestemp=[]
        for p in range(0,len(pairslist[t]),1):
            x1=pairslist[t][p][0]
            y1=pairslist[t][p][1]
            z1=pairslist[t][p][2]
            x2=pairslist[t][p][3]
            y2=pairslist[t][p][4]
            xtemp=0
            ytemp=0
            ztemp=0
            zindex=t
            comp=100000000.0
            if zindex ==size-1:
                zindex =-1
            for p2 in range(0,len(pairslist[zindex+1]),1):
                x2=pairslist[zindex+1][p2][0]
                y2=pairslist[zindex+1][p2][1]
                z2=pairslist[zindex+1][p2][2]
                temp=computeDistance(x1,y1,z1,x2,y2,z2,size)
                if temp<8:
                    if temp<comp:
                        comp=temp 
                        xtemp=x2
                        ytemp=y2
                        ztemp=z2
            complist.append(comp)
            if comp<8 :
                linestemp.append((x1,y1,z1,xtemp,ytemp,ztemp))
            else:
                linestemp.append((x1,y1,z1,x1,y1,z1))
            if temp<=1:
                p2=len(pairslist[zindex+1])
        lines.append(linestemp)
    return lines
def frontTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    z1a=z1
    z2a=z2
    z3a=z3
    if z1-z2>2:
        z2a = z1+z2
    if z2-z1>2:
        z1a = z1+z2
    if z1-z3>2:
        z3a = z1+z3
    if z3-z1>2:
        z1a = z1+z3    
    if z3-z2>2:
        z2a = z3+z2
    if z2-z3>2:
        z3a = z3+z2    
        
   # triangle0 = triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    center0 = triangle.center(x1,y1,z1a,x2,y2,z2a,x3,y3,z3a)
    area0 = triangle.area(x1,y1,z1a,x2,y2,z2a,x3,y3,z3a)
    #print(area0)
    nhat0 = triangle.nhat(area0,x1,y1,z1a,x2,y2,z2a,x3,y3,z3a)
    triangle1= triangles(area0,center0,(x1,y1,z1a),(x2,y2,z2a),(x3,y3,z3a),nhat0)
    return triangle1
    
def callsTriangle(zslice,pairslist, upleft, upright, bottomright, internalTriangleList,countCalls,size):
    if countCalls<10:
        for ind in range(zslice+1,zslice+2,1):
            if ind<size-1:  #is this doing the top right? 
                for h in range(len(pairslist[ind])):
                    #(x,y,z,x+deltaX,y+deltaY)
                    if((pairslist[ind][h][0],pairslist[ind][h][1],pairslist[ind][h][2])==upleft):
                        upmid = pairslist[ind][h][3],pairslist[ind][h][4],pairslist[ind][h][2]
                        if upmid==upright:
                            triangle1= frontTriangle(bottomright[0], bottomright[1], bottomright[2], upleft[0], upleft[1], upleft[2], upmid[0], upmid[1], upmid[2])
                            triangle2= frontTriangle(bottomright[0], bottomright[1], bottomright[2], upright[0], upright[1], upright[2], upmid[0], upmid[1], upmid[2])
                            internalTriangleList.append(triangle1)
                            internalTriangleList.append(triangle2)
                            
                        else:
                            
                            triangle0 = frontTriangle(bottomright[0], bottomright[1], bottomright[2], upleft[0], upleft[1], upleft[2], upmid[0], upmid[1], upmid[2])
                            internalTriangleList.append((triangle0))
                            callsTriangle(zslice,pairslist,upmid,upright,bottomright,internalTriangleList,countCalls+1,size)

                            
def flatTop2(zslice,pairslist, baseleft, baseright, internalTrianglesTop,w):
    w+=1
    for pa in range(len(pairslist[zslice])):
        if((pairslist[zslice][pa][3],pairslist[zslice][pa][4],pairslist[zslice][pa][2])==baseleft):
            up=pairslist[zslice][pa][0],pairslist[zslice][pa][1],pairslist[zslice][pa][2]
            triangle0 = frontTriangle(baseleft[0],baseleft[1],baseleft[2],baseright[0],baseright[1],baseright[2],up[0],up[1],up[2])
            if (triangle0.vertex1[2]==zslice and triangle0.vertex2[2]==zslice and triangle0.vertex3[2]==zslice):
                internalTrianglesTop.append(triangle0)
                
            if ((baseleft[0],baseleft[1],baseleft[2],baseright[0],baseright[1]) in pairslist[zslice]):
                pairslist[zslice].remove((baseleft[0],baseleft[1],baseleft[2],baseright[0],baseright[1]))
                pairslist[zslice].append((10000000,100000000,100000000,10000000,10000000))
    
            if ((up[0],up[1],up[2],baseleft[0],baseleft[1]) in pairslist[zslice]):
                pairslist[zslice].remove((up[0],up[1],up[2],baseleft[0],baseleft[1]))
                pairslist[zslice].append((100000,100000,100000,100000,100000))
            if ((baseleft[0],baseleft[1],up[2],up[0],up[1]) in pairslist[zslice]):
                pairslist[zslice].remove((baseleft[0],baseleft[1],up[2],up[0],up[1]))
                pairslist[zslice].append((100000,100000,100000,100000,100000))   

            temp=up
            flatTop2(zslice,pairslist,temp,baseright,internalTrianglesTop,w)
            
        elif((pairslist[zslice][pa][3],pairslist[zslice][pa][4],pairslist[zslice][pa][2])==baseright and (pairslist[zslice][pa][0],pairslist[zslice][pa][1],pairslist[zslice][pa][2]) != baseleft):
           up=pairslist[zslice][pa][0],pairslist[zslice][pa][1],pairslist[zslice][pa][2]
           triangle0 = frontTriangle(baseleft[0],baseleft[1],baseleft[2],baseright[0],baseright[1],baseright[2],up[0],up[1],up[2])
           if (triangle0.vertex1[2]==zslice and triangle0.vertex2[2]==zslice and triangle0.vertex3[2]==zslice):
               internalTrianglesTop.append(triangle0)
              
           if ((baseleft[0],baseleft[1],baseleft[2],baseright[0],baseright[1]) in pairslist[zslice]):
               pairslist[zslice].remove((baseleft[0],baseleft[1],baseleft[2],baseright[0],baseright[1]))
               pairslist[zslice].append((1000000,10000000,1000000,1000000,1000000))   
           if((baseright[0],baseright[1],baseright[2],up[0],up[1]) in pairslist[zslice]):
                pairslist[zslice].remove((baseright[0],baseright[1],baseright[2],up[0],up[1]))
                pairslist[zslice].append((100000,100000,100000,100000,100000)) 
           if((up[0],up[1],up[2],baseright[0],baseright[1]) in pairslist[zslice]):
                pairslist[zslice].remove((up[0],up[1],up[2],baseright[0],baseright[1]))
                pairslist[zslice].append((100000,100000,100000,100000,100000)) 
                   
           temp=up
           flatTop2(zslice,pairslist,temp,baseleft,internalTrianglesTop,w)
           
def findFronts(pairslist,lines,size):
    trianglelist=[]            
              
    for i in range(0,size,1): #i indexes zslice
        for j in range(0,len(pairslist[i]),1): #j loops over the pairs in zslice i
            x1=pairslist[i][j][0]
            y1=pairslist[i][j][1]
            z1=pairslist[i][j][2]
            x2=pairslist[i][j][3]
            y2=pairslist[i][j][4]
            #bottomleft = (x1,y1,z1)
            bottomright = (x2,y2,z1)
            zindex = i
            triangleadd =0
            if zindex ==size-1:
                zindex =-1
            st =0 
            while st ==0:
                for l in range(len(lines[i])):
                    if lines[i][l][0]==x1 and lines[i][l][1]==y1 and lines[i][l][2]==z1:
                        up1 = (lines[i][l][3],lines[i][l][4],lines[i][l][5])
                        cont1=0
                        if up1 == (x1,y1,z1):
                            w=0
                            flatTop2(zindex, pairslist,(x1,y1,z1),(x2,y2,z1),trianglelist,w)
                            
                            st=1
                           
                        else:
                            for ll in range(0,len(lines[i]),1):
                                if (lines[i][ll][0],lines[i][ll][1],lines[i][ll][2])==bottomright:
                                    up2 = (lines[i][ll][3],lines[i][ll][4],lines[i][ll][5])
                                    if up1==up2 or up2==bottomright:
                                        triangle= frontTriangle(x1, y1, z1, x2, y2, z1,up1[0],up1[1], up1[2])
                                        trianglelist.append((triangle))
                                        triangleadd+=1
                                        st=1  
                                    else :
                                        cont1+=1
                                        triangle0=frontTriangle(x1,y1,z1,x2,y2,z1,up1[0],up1[1],up1[2])
                                        triangle1=frontTriangle(x2,y2,z1,up2[0],up2[1],up2[2],up1[0],up1[1],up1[2])
                                        trianglelist.append((triangle0))
                                        trianglelist.append((triangle1))
                                        midind = up1[2]
                                        
                                        callsTriangle(midind,pairslist,up1,up2,bottomright, trianglelist,0,size)
                 
                st=1
    return trianglelist
def Periodic(x,y,z,size):
    return [(x+size)%size, (y+size)%size, (z+size)%size]
    # if x==size:
    #     x=0
    # if y==size:
    #     y=0
    # if z==size:
    #     z=0
    # if x==-1:
    #     x=size-1
    # if y==-1:
    #     y=size-1
    # if z==-1:
    #     z=size-1
    # return [x,y,z] 

def findGamma12Loc(trianglesfinal,box1,size,gammabox):
    gamma12pointfortriangleloc=np.zeros((len(trianglesfinal),3),dtype=int)
    l_step = np.sqrt(3)
    for r in range(len(trianglesfinal)):
        if trianglesfinal[r].area >.00001:
            (cx,cy,cz) = trianglesfinal[r].center
            (nx,ny,nz) =  trianglesfinal[r].nhat
            
            
            loc = np.rint([cx-l_step * nx,cy-l_step *ny,cz-l_step *nz])
            #print(loc)
            #print(r)
            loc = [int(loc[0]),int(loc[1]),int(loc[2])]
            loc = Periodic(loc[0],loc[1],loc[2],size)
            #print(loc)
            if (box1[loc[0]][loc[1]][loc[2]] > .5):
                trianglesfinal[r].nhat= [-1*nx,-1 *ny,-1*nz]
                #trianglesfinal[r][6]= -1*ny
                #trianglesfinal[r][7]= -1*nz #need to deal with this eventually

                loc =  np.rint([cx+l_step *nx,cy+l_step *ny,cz+l_step *nz])
                loc = Periodic(loc[0],loc[1],loc[2],size)
                loc = [int(loc[0]),int(loc[1]),int(loc[2])]

            # If gamma12 value is 0, point further along the nhat direction
            # while (gammabox[loc[0]][loc[1]][loc[2]] == 0):
            #     trianglesfinal[r].nhat= [-1*nx,-1*ny,-1*nz]
            #     #trianglesfinal[r][6]= -1*ny
            #     #trianglesfinal[r][7]= -1*nz #need to deal with this eventually 
            #     loc =  np.rint([loc[0] - 0.2*trianglesfinal[r].nhat[0], loc[1] - 0.2 * trianglesfinal[r].nhat[1], loc[2] - 0.2 * trianglesfinal[r].nhat[2]])
            #     loc = Periodic(loc[0],loc[1],loc[2],size)
            #     loc = [int(loc[0]),int(loc[1]),int(loc[2])]

            gamma12pointfortriangleloc[r][0],gamma12pointfortriangleloc[r][1],gamma12pointfortriangleloc[r][2] = loc[0],loc[1],loc[2]
    return gamma12pointfortriangleloc

def findGamma12Val(gamma12loc,trianglesfinal,gammabox):
    gamma12pointfortriangleval=np.zeros((len(trianglesfinal)))
    gammabox_max_filtered = scipy.ndimage.maximum_filter(gammabox,size=(len(gammabox),len(gammabox),len(gammabox)), mode='wrap')
    for r in range(len(trianglesfinal)):
        if trianglesfinal[r].area >.00001:
            gamma12pointfortriangleval[r] = gammabox[int(gamma12loc[r][0])][int(gamma12loc[r][1])][int(gamma12loc[r][2])]    
        if gamma12pointfortriangleval[r]==0:
            gamma12pointfortriangleval[r] = gammabox_max_filtered[int(gamma12loc[r][0])][int(gamma12loc[r][1])][int(gamma12loc[r][2])]

    return gamma12pointfortriangleval                              
            
def CleanUpList(trianglelist):
    trianglesfinal = [ ele for ele in trianglelist if ele is not None ]
    tyTest=np.full((len(trianglesfinal)),False,dtype=bool)   
    tyTest = [True for ele in trianglesfinal if type(ele) is list]
    exists =  True in tyTest

    while exists == True:
        for h in range(len(trianglesfinal)):
            if type(trianglesfinal[h]) == list:
                for f in range(len(trianglesfinal[h])-1,-1,-1):
                    trianglesfinal.append(trianglesfinal[h][f])
                trianglesfinal.pop(h)
            if [] in trianglesfinal==True:
                trianglesfinal.remove([])
            tyTest=np.full((len(trianglesfinal)),False,dtype=bool)    
            tyTest = [True for ele in trianglesfinal if type(ele) is list]
            exists =  True in tyTest
                
    trianglesfinal = [ ele for ele in trianglesfinal if ele is not None ]    
    return trianglesfinal
    
                

class triangle:
    
    def __init__(self,x1,y1,z1,x2,y2,z2,x3,y3,z3):
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.y1 = y1
        self.y2 = y2
        self.y3 = y3
        self.z1 = z1
        self.z2 = z2
        self.z3 = z3
        
        #@staticmethod
    def area(x1,y1,z1,x2,y2,z2,x3,y3,z3):
        a= np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        b= np.sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
        c= np.sqrt((x1-x3)**2+(y1-y3)**2+(z1-z3)**2)
        s = (a+b+c)/2.0
        area1 = np.sqrt(s*(s-a)*(s-b)*(s-c))
        return area1
    #@staticmethod
    def center(x1,y1,z1,x2,y2,z2,x3,y3,z3):
        center1 = ((x1+x2+x3)/3,(y1+y2+y3)/3,(z1+z2+z3)/3)
        return center1
    #@staticmethod
    def nhat(area1,x1,y1,z1,x2,y2,z2,x3,y3,z3):
        leg1 = [x1-x2,y1-y2,z1-z2]
        leg2 = [x2-x3,y2-y3,z2-z3]
        nfull = np.cross(leg1,leg2)
        nhat = nfull/np.sqrt(nfull[0]**2 + nfull[1]**2 +nfull[2]**2)
        #area1 = self.area(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        if area1==0.0:
            nhat = [1000,1000,1000]
        return nhat
            
        

class triangles:
    
    
    def __init__(self,area,center,vertex1,vertex2,vertex3,nhat):
        self.area = area
        self.center = center
        self.vertex1 = vertex1
        self.vertex2= vertex2
        self.vertex3= vertex3
        self.nhat = nhat
        
    @property
    def vertex1(self):
        return  self._vertex1
    @property
    def vertex2(self):
        return  self._vertex2
    @property
    def vertex3(self):
         return  self._vertex3
    """
    @property
    def area(self):
        return self.area
    @property
    def center(self):
        return self.center
    @property
    def nhat(self):
        return self.nhat
    """
    @vertex1.setter
    def vertex1(self,a):
        self._vertex1 =(a[0],a[1],a[2])
    @vertex2.setter
    def vertex2(self,a):
        self._vertex2 =(a[0],a[1],a[2])
    @vertex3.setter
    def vertex3(self,a):
        self._vertex3 =(a[0],a[1],a[2])
    """
    @area.setter
    def area(self):
        (x1,y1,z1) = self.vertex1
        (x2,y2,z2) = self.vertex2
        (x3,y3,z3) = self.vertex2
        a= np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        b= np.sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
        c= np.sqrt((x1-x3)**2+(y1-y3)**2+(z1-z3)**2)
        s = (a+b+c)/2.0
        area1 = np.sqrt(s*(s-a)*(s-b)*(s-c))
        self._area = area1
    @center.setter
    def center(self):
        (x1,y1,z1) = self.vertex1
        (x2,y2,z2) = self.vertex2
        (x3,y3,z3) = self.vertex2
        center1 = ((x1+x2+x3)/3,(y1+y2+y3)/3,(z1+z2+z3)/3)
        self._center= center1
    @nhat.setter
    def nhat(self):
        (x1,y1,z1) = self.vertex1
        (x2,y2,z2) = self.vertex2
        (x3,y3,z3) = self.vertex2
        leg1 = [x1-x2,y1-y2,z1-z2]
        leg2 = [x2-x3,y2-y3,z2-z3]
        nfull = np.cross(leg1,leg2)
        nhat = nfull/np.sqrt(nfull[0]**2 + nfull[1]**2 +nfull[2]**2)
        
        if self.area==0.0:
            nhat = [1000,1000,1000]
        self._nhat = nhat
    """    
class Edgefinder:
    
    def __init__(self,box,gammabox): #fronts,gamma12loc,gamma12val):
        self.box=box
        #self.fronts = fronts
        self.gammabox = gammabox
        #self.gamma12loc = gamma12loc
        #self.gamma12val = gamma12val
        
        
    @property
    def box(self):
        return self._box
    #@property
    #def fronts(self):
    #    return self._fronts
    
    #@property
    #def gamma12loc(self):
    #    return self._gamma12loc
    #@property
    #def gamma12val(self):
    #    return self._gamma12val
    @property
    def gammabox(self):
        return self._gammabox
    
    
    @box.setter
    def box(self,box):
        self._box = box
        

    @gammabox.setter
    def gammabox(self,gammabox1):
        self._gammabox = gammabox1    
    
    #@fronts.setter
    def pairs(self):
        return TwoDContours(self.box)
    def lines(self,pairslist):
        #pairslist = TwoDContours(self.box)
        linesup = ThreeDlineslist(pairslist,len(self.box[0]))
        return linesup
    def fronts(self,pairslist,linesup):
        #pairslist = TwoDContours(self.box)
        #linesup = ThreeDlineslist(pairslist,len(self.box[0]))
        #return pairslist
        inbetweentriangles =findFronts(pairslist,linesup,len(self.box[0]))
        trianglesfinal =CleanUpList(inbetweentriangles)
        return trianglesfinal
        #self._fronts = trianglesfinal
        
    #@gamma12loc.setter
    def gamma12loc(self,trianglelist):
        return findGamma12Loc(trianglelist, self.box, len(self.box[0]), self.gammabox)
        #self._gamma12loc = findGamma12Loc(self.fronts,self.box,len(self.box[0]))
    
    #@gamma12val.setter
    def gamma12val(self,trianglelist,gammaloc):
        return findGamma12Val(gammaloc, trianglelist, self.gammabox)
        #self._gamma12val = findGamma12Val(self.gamma12loc,self.fronts,self.gammabox)
        
  
"""    This is what I used to test the above work and it matches what i was getting before      
    
ionized_box = p21c.cache_tools.readbox(direc ="/Users/Emily/21cmFAST-cache", fname="IonizedBox_165367d07d17980e2c10fd22abfeb0be_r54321.h5")    
test = Edgefinder(getattr(ionized_box,'xH_box'),getattr(ionized_box,'Gamma12_box')) 
pairs1 = test.pairs() 
lines1= test.lines(pairs1)  
fronts = test.fronts(pairs1,lines1)  #fronts made of triangles class objects (correct number)
gamma12locations= test.gamma12loc(fronts)
gamma12values = test.gamma12val(fronts,gamma12locations)
           
       
"""       
        
        
    
        
    