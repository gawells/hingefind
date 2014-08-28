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
        self.cutdom = cutdom

        self.mobFullCA = self.mobile.select("name CA and protein")
        self.refFullCA = self.reference.select("name CA and protein")
        self.mobCA = self.mobile.select("name CA and protein")
        self.refCA = self.reference.select("name CA and protein")

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

    def criterion (self,mob_sel,ref_sel,eps):
    # def criterion (self,mobile,reference,eps):
        mob_coord = mob_sel.getCoords()
        ref_coord = ref_sel.getCoords()

        mob_ids =mob_sel.getIndices()
        ref_ids =ref_sel.getIndices()

        mob_list = []
        ref_list = []

        i = 0

        for (c1,c2) in [(c1,c2) for c1 in mob_coord for c2 in ref_coord]:
            if pd.calcDistance(c1,c2) < eps:
                mob_list.append(mob_ids[i])
                ref_list.append(ref_ids[i])
                i += 1
        
        return [mob_list,ref_list]

    def arraystr(a):
    	return " ".join([str(x) for x in a])

    def superimpose(self,mob_sel,ref_sel,eps):
    	survivor = criterion(mob_sel,ref_sel,eps)
    	mob_subset = survivor[0]
    	ref_subset = survivor[1]
    	if len(mob_subset) < 5:
    		return survivor

    	mob_subset_sel = self.mobile.select("index "+arraystr(mob_subset))
    	ref_subset_sel = self.reference.select("index "+arraystr(ref_subset))
    	t = pd.calcTransformation(mob_subset_sel,ref_subset_sel)
    	pd.applyTransformation(t,self.mobile)
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
    	return [mob_list,ref_list]

    def convergence(self,eps,radius):
    	startids = seed(self.mobCA,self.refCA,radius)
    	startids_mob = startids[0]
    	startids_ref = startids[1]
    	if len(startids_mob) < 5:
    		return startids

    	mob_start_sel = self.mobile.select("index "+arraystr(startids_mob))
    	ref_start_sel = self.reference.select("index "+arraystr(startids_ref))
    	t = pd.calcTransformation(mob_start_sel,ref_start_sel,weights=mob_start_sel.getMasses())

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
    
    # # test criterion
    # selstr = "resid 92 to 100 104 to 140 143 to 218 220 to 250 and name CA"
    # moveBy = pd.calcTransformation(mobile.select(selstr),reference.select(selstr))
    # pd.applyTransformation(moveBy,mobile)
    # print hf.criterion(selstr,selstr,1.0)

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
