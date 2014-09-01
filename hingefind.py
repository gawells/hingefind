#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import prody as pd
import transformations as tf
import os 
import sys
import argparse

class HingeFind:
    """docstring for HingeFind"""
    def __init__(self, mobile,reference,initrad=20,maxcounter=20,ndomains=999,cutdom=10):
        # super(HingeFind, self).__init__()
        self.mobile = mobile
        self.reference = reference
        self.initrad = initrad
        self.maxcounter = maxcounter
        self.ndomains = ndomains
        self.domain_count = 0
        self.cutdom = cutdom

        self.mobFullCA = self.mobile.select("name CA and protein ")
        self.refFullCA = self.reference.select("name CA and protein ")
        self.mobCA = self.mobile.select("name CA and protein ")
        self.refCA = self.reference.select("name CA and protein ")

        self.mobDomList = {}
        self.refDomList = {}

        self.dindex = []
        self.sup_rmsd = []

        try :     
            self.refCA.numAtoms() == self.mobCA.numAtoms()
        except :
            print "init> two selections need same number of atoms!"

        try:
            self.cutdom <= 5
        except:
            print "init> cutdom needs to be larger than 4!"

    def lsubstract(self,list1,list2):
        return [l for l in list1 if l not in list2]

    def arraystr(self,a):
        return " ".join([str(x) for x in a])

    def indexSel(self,structure,indices):        
        if len(indices) > 0:
            return structure.select("index "+self.arraystr(indices))
        else :
            return None

    def resStr(self,structure):                
        if structure:
            return self.arraystr(structure.getResnums())
        else :
            return None
    
    def criterion (self,mob_sel,ref_sel,eps):
    # def criterion (self,mobile,reference,eps):
        mob_coord = mob_sel.getCoords()
        ref_coord = ref_sel.getCoords()

        mob_ids =mob_sel.getIndices()
        ref_ids =ref_sel.getIndices()

        mob_list = []
        ref_list = []

        i = 0

        while i < len(mob_coord):
        # for (c1,c2) in [(c1,c2) for c1 in mob_coord for c2 in ref_coord]:            
            if pd.calcDistance(mob_coord[i],ref_coord[i]) < eps:
                mob_list.append(mob_ids[i])
                ref_list.append(ref_ids[i])
            i += 1  
        
        mob_surv = self.indexSel(mob_sel,mob_list)
        ref_surv = self.indexSel(ref_sel,ref_list)

        return [mob_surv,ref_surv]
        # return [mob_list,ref_list]

    def superimpose(self,eps):
    	survivor = self.criterion(self.mobCA,self.refCA,eps)
    	mob_subset = survivor[0]
    	ref_subset = survivor[1]
    	if mob_subset:
            if mob_subset.numAtoms() < 5:
    		  return survivor

            t = pd.calcTransformation(mob_subset,ref_subset)
            pd.applyTransformation(t,self.mobile)
        else:
            return [None,None]

    	return survivor    	

    def seed(self,mob_sel,ref_sel,radius):
    	mobcoord = mob_sel.getCoords()
    	mob_ids = mob_sel.getIndices()
    	ref_ids = ref_sel.getIndices()

    	mob_cd1 = mobcoord[0]
    	mob_id1 = mob_ids[0]

    	mob_list = []
    	ref_list = []

    	i = 0
    	for cd1 in mobcoord:
    		mob_id = mob_ids[i]
    		ref_id = ref_ids[i]
    		if pd.calcDistance(cd1,mob_cd1) < radius:
    			mob_list.append(mob_id)
    			ref_list.append(ref_id)
    		i +=1

        mob_subset = self.indexSel(self.mobile, mob_list)
        ref_subset = self.indexSel(self.reference, ref_list)
    	return [mob_subset,ref_subset]

    def convergence(self,eps,radius):
    	startsel = self.seed(self.mobCA,self.refCA,radius)
    	mob_start_sel = startsel[0]
    	ref_start_sel = startsel[1]
        print "Seed: %s"%mob_start_sel.getResnums()

    	if mob_start_sel.numAtoms() < 5:
    		return startsel
        t = pd.calcTransformation(mob_start_sel,ref_start_sel,weights=mob_start_sel.getMasses())
        pd.applyTransformation(t,self.mobile)

        count = 0
        prev = self.superimpose(eps)
        curr = self.superimpose(eps)

        # while self.arraystr(prev[0].getResnums()) != self.arraystr(curr[0].getResnums()):
        while self.resStr(prev[0]) != self.resStr(curr[0]):
            prev = curr
            curr = self.superimpose(eps)

            count += 1

            if count == self.maxcounter:
                print "convergence - warning: a domain did not converge."
                break

            if curr[0].numAtoms() < 5 or not curr[0]:
                print "Less than 5"
                return [[None], [None]]

        return curr

    def update_boundaries(self,eps,number):
        mob_coord = self.mobFullCA.getCoords()
        ref_coord = self.refFullCA.getCoords()

        newsup_rmsd = []
        newindex = []

        i = 0

        while i < len(mob_coord):
            mcrd = mob_coord[i]
            rcrd = ref_coord[i]
            dist = pd.calcDistance(mcrd,rcrd)
            if dist < self.sup_rmsd[i]:
                newsup_rmsd.append(dist)
                if dist < eps:
                    newindex.append(number)
                else:
                    newindex.append(self.dindex[i])
            else:
                newsup_rmsd.append(self.sup_rmsd[i])
                newindex.append(self.dindex[i])

            i += 1

        self.sup_rmsd = newsup_rmsd
        self.dindex = newindex

    def partition(self,eps):
        # initialize
        # self.mobCA = self.mobFullCA
        # self.refCA = self.refFullCA
        self.sup_rmsd = [999]*self.mobFullCA.numAtoms()
        self.dindex = [-1]*self.mobFullCA.numAtoms()
        self.domain_count = 0
        seedradius = self.initrad

        # partition the protein
        while self.domain_count < self.ndomains:
            domain = self.convergence(eps, seedradius)

            if domain[0]:
                if domain[0].numAtoms() > 0:
                    print "Domain: %s Size: %d"%(domain[0].getResnums(),len(domain[0].getResnums()))

            mob_domain = domain[0]
            ref_domain = domain[1]

            if not mob_domain or mob_domain.numAtoms() == 0:
                print "partition> trying 20% smaller seed radius..."
                seedradius *= 0.8
                continue

            # remove from ca1 and ca2 those atoms which have been found
            mob_dom_list = mob_domain.getIndices()
            ref_dom_list = ref_domain.getIndices()
            mobCA_list = self.mobCA.getIndices()
            refCA_list = self.refCA.getIndices()
            mob_disjoint_l = [i for i in mobCA_list if not i in mob_dom_list]            
            ref_disjoint_l = [i for i in refCA_list if not i in ref_dom_list]

            if len(mob_disjoint_l) > 0:
                self.mobCA = self.indexSel(self.mobile, mob_disjoint_l)
                self.refCA = self.indexSel(self.reference, ref_disjoint_l)

            if mob_domain.numAtoms() > 4:
                seedradius = self.initrad
                self.update_boundaries(eps, self.domain_count)

            self.domain_count += 1

            if self.domain_count == 1: 
                print "partition> now have 1 domain"
            else:
                print "partition> now have %d domains"%self.domain_count

            if len(mob_disjoint_l) < self.cutdom:
                print "partittion> protein now fully partitioned"
                break

            self.mobDomList = {}
            self.refDomList = {}
            currentDomains = set(self.dindex)        
            for d in currentDomains:
                self.mobDomList[d] = []
                self.refDomList[d] = []

            mob_ids = self.mobFullCA.getIndices()
            ref_ids = self.refFullCA.getIndices()

            i = 0
            for d in self.dindex:
                self.mobDomList[d].append(mob_ids[i])
                self.refDomList[d].append(ref_ids[i])
                i += 1

            # print self.mobCA.getResnums()
            # for d in self.mobDomList.keys():
            #     print self.indexSel(self.mobile, self.mobDomList[d]).getResnums()


    def hinge (beforeM, afterM, beforeR, afterR, mobile, rotmat): 
        '''
        hinge from hingefind.tcl recapitulated in python
        
        '''
        
        com1 = calcCenter(beforeM,beforeM.getMasses())
        # print "COM beforeM: %s"%com1
        com2 = calcCenter(afterM,afterM.getMasses())
        # print "COM afterM: %s"%com2
        cdis = np.linalg.norm(com2-com1)
        print "Distance moved: %f Å"%cdis

        # calculate bisecting point and normalvector of bisecting plane
        bi = (com1 + com2) / 2
        pl = (com1 - com2) / np.abs(np.linalg.norm(com1 - com2))

        # compute rotation matrix for rotations about com2
        t21 = tf.translation_matrix(com2 - com1)
        rot2 = np.dot(rotmat,t21)

        # rotate com1 twice about com2
        # (arbitrary choice of rotation to generate new points to make a plane,
        # rotation axis is perpindicular to this?)
        p1 = np.dot(rot2,np.append(com2,1))
        p2 = np.dot(rot2,p1)

        # compute the rotation axis from the 3 points
        rideal =  np.cross(p2[:3] - com2, p1[:3] - com2)
        rideal /= np.linalg.norm(rideal)

        # project com2 onto the 3-point-plane
        new = com1 - (rideal * np.dot(rideal,(com1 - com2)))

        # compute rotation angle (least squares fit)
        cosine = np.dot((new - com2)/np.abs(np.linalg.norm(new - com2)), \
            (new - p1[:3])/np.abs(np.linalg.norm(new - p1[:3])))
        angl = arccos(cosine)
        angl_deg = angl * 180 / np.pi
        print "Rotation: %f°"%angl_deg

        # compute projection of rot axis on bisecting plane
        perp = np.dot(rideal, pl)
        angp = np.abs(arcsin(perp))
        pro = (rideal - (pl*perp))/np.abs(np.linalg.norm(rideal - (pl*perp)))

        # compute decomposition angle
        tang = cos(angp) * tan(angl*0.5)
        angle = 2*arctan(tang)

        deg_angle = angle * 180 / np.pi
        print "Hingefind effective rotation: %f°"%deg_angle

        # compute pivot point
        hi = bi + np.cross(pl, pro)*0.5*cdis/tang

        # # translate by effective rotation and reset
        # (rotate by $angle around $hi in the plane of $com1,$com2 and $hi?)
        t = Transformation(tf.rotation_matrix(angle, np.cross((com2 - hi),(com1 - hi)),hi))
        com3 = calcCenter(beforeM,weights=beforeM.getMasses())
        applyTransformation(t,mobile)
        com3 = calcCenter(beforeM,weights=beforeM.getMasses())
        rmspro = calcRMSD(beforeM,afterM,weights=beforeM.getMasses())
        t = calcTransformation(afterR,beforeR)
        applyTransformation(t,mobile)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--reference",type=str,help="Reference structure")
    parser.add_argument("-m","--mobile",type=str,help="Mobile structure")
    parser.add_argument("-d1","--domain1",type=str,help="Atom selection for domain 1")
    parser.add_argument("-d2","--domain2",type=str,help="Atom selection for domain 2")

    args = parser.parse_args()
    kwargs = {}

    if not args.reference:
        args.reference = "1lfh.pdb"
    
    if not args.mobile:
        args.mobile = "1lfg.pdb"

    if not args.domain1:
        args.domain1 = "resid 5 to 85 87 88 90 91 251 to 293 295 to 302 304 to 332 685 687 to 689 691 and name CA"

    if not args.domain2:
        args.domain2 = "resid 92 to 100 104 to 140 143 to 218 220 to 250 and name CA"

    reference = pd.parsePDB(args.reference)
    mobile = pd.parsePDB(args.mobile)

    hf = HingeFind(mobile,reference)
    
    # test criterion
    # selstr = "resid 92 to 100 104 to 140 143 to 218 220 to 250 and name CA"
    # moveBy = pd.calcTransformation(mobile.select(selstr),reference.select(selstr))
    # pd.applyTransformation(moveBy,mobile)
    # print hf.criterion(mobile.select(selstr),reference.select(selstr),1.0)

    hf.partition(1.0)

    # D1str = args.domain1
    # D1_0 = structure1.select(D1str)
    # D1_now = structure2.select(D1str)
    
    # t = calcTransformation(D1_now,D1_0,weights=D1_now.getMasses())
    # applyTransformation(t,structure2)   
    # writePDB('sup1.pdb',structure2)

    # D2str = args.domain2
    # D2_0 = structure1.select(D2str)    
    # D2_now = structure2.select(D2str)
    # t = calcTransformation(D2_now,D2_0,weights=D2_now.getMasses())
    
    # # from hingefind
    # print "**Hingefind**"
    # hinge(D2_0, D2_now,D1_0,D1_now,structure2,t.getMatrix())
    
    # print "**Prody**"
    # # print t.getMatrix()
    # q = tf.quaternion_from_matrix(t.getMatrix()) ##
    # # print t.getRotation()
    # print "Quaternion from matrix: %s"%q
    # rotation = (2*arccos(q[0]))/np.pi*180
    # print  "Rotation from quaternion: %f°"%rotation
    
    # applyTransformation(t,structure2)   
    # writePDB('sup2.pdb',structure2)

    # From hingefind : D0 D1 = 18.83°
    # From hingefind : D1 D0 = 19.80°

if __name__ == '__main__':
    main()
