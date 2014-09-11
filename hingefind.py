#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import prody as pd              #  https://github.com/prody/ProDy.git
import transformations as tf    #  https://github.com/malcolmreynolds/transformations.git
import argparse
import re

class HingeFind:
    """docstring for HingeFind"""
    def __init__(self, mobile,reference,initrad=20,maxcounter=20,ndomains=999,cutdom=10,\
        trajectory=None,structure=None,output="hingefind_domains"):
        self.mobFileName = mobile
        self.mobile = pd.parsePDB(mobile)
        self.refFileName = reference
        self.reference = pd.parsePDB(reference)
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

        self.trajectory = trajectory
        self.structure = structure
        self.output = output

        # if structure and trajectory:
        #     self.structure = pd.parsePDB(structure)
        #     self.trajectory = pd.Trajectory(trajectory)
        #     self.trajectory.link(self.structure)

        try :     
            self.refCA.numAtoms() == self.mobCA.numAtoms()
        except :
            print "init> two selections need same number of atoms!"

        try:
            self.cutdom <= 5
        except:
            print "init> cutdom needs to be larger than 4!"


    def lsubstract(self,list1,list2):
        '''
        Return disjoint list of elements in list1 but not in list2

        '''
        return [l for l in list1 if l not in list2]


    def arraystr(self,a):
        '''
        Return list as space delimited string
        '''
        return " ".join([str(x) for x in a])


    def indexSel(self,structure,indices):        
        '''
        Return atom selection correpsoding to index list
        '''
        if len(indices) > 0:
            return structure.select("index "+self.arraystr(indices))
        else :
            return None


    def resStr(self,structure):                   
        '''
        Return space delimited string of residue numbers in input structure

        '''         
        if structure:
            return self.arraystr(structure.getResnums())
        else :
            return None

    
    def criterion (self,mob_sel,ref_sel,eps):
        '''
        - find those atoms which have rmsd < eps in the two selections
        - return the info as a list of prody Cα atom selections 
        - called by procedure 'superimpose'.
        '''

        mob_coord = mob_sel.getCoords()
        ref_coord = ref_sel.getCoords()
        mob_ids =mob_sel.getIndices()
        ref_ids =ref_sel.getIndices()
        mob_list = []
        ref_list = []

        i = 0
        while i < len(mob_coord):
            if pd.calcDistance(mob_coord[i],ref_coord[i]) < eps:
                mob_list.append(mob_ids[i])
                ref_list.append(ref_ids[i])
            i += 1  
        
        mob_surv = self.indexSel(mob_sel,mob_list)
        ref_surv = self.indexSel(ref_sel,ref_list)

        return [mob_surv,ref_surv]


    def superimpose(self,eps):
        '''
        - superposition based on the atoms which match the best.
        - return the best matching subsets as a list of prody Cα atom selections.
        - called by procedure 'convergence'.
        '''

        # get atoms that are within 'eps'
    	survivor = self.criterion(self.mobCA,self.refCA,eps)
    	mob_subset = survivor[0]
    	ref_subset = survivor[1]

        # and fit by the best matching set
        if mob_subset:
            if mob_subset.numAtoms() < 5:
    		  return survivor

            t = pd.calcTransformation(mob_subset,ref_subset)
            pd.applyTransformation(t,self.mobile)
        else:
            return [None,None]

    	return survivor    	

    def seed(self,mob_sel,ref_sel,radius):
        '''
        - find a seed subset (all atoms of 'mob_sel' within 'radius' of the
          first atom in 'mob_sel').
        - return the info as a list of prody Cα atom selections.
        - called by procedure 'convergence'.
        '''

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
        '''
          - iterate until a domain is found
.          - return domain as a list of prody Cα atom selections.
          - called by procedure 'partition'.
        '''

        # initiate search in seed-subset
    	startsel = self.seed(self.mobCA,self.refCA,radius)
    	mob_start_sel = startsel[0]
    	ref_start_sel = startsel[1]
    	if mob_start_sel.numAtoms() < 5:
    		return startsel
            # nothing to find from seed set, return formally as a domain

        # if seed-subset sufficiently large, proceed with initial fit
        t = pd.calcTransformation(mob_start_sel,ref_start_sel,weights=mob_start_sel.getMasses())
        pd.applyTransformation(t,self.mobile)

        # search iteratively for the domain until residue list converges
        count = 0
        prev = self.superimpose(eps)
        curr = self.superimpose(eps)
        while self.resStr(prev[0]) != self.resStr(curr[0]):
            prev = curr
            curr = self.superimpose(eps)

            count += 1

            if count == self.maxcounter:
                print "convergence - warning: a domain did not converge."                
                break

        if not curr[0] or curr[0].numAtoms() < 5:
            return [None, None]
            # even though seed set was sufficiently large, 
            # the found domain is too small. return nothing and 
            # let 'partition' try again with smaller radius

        return curr


    def update_boundaries(self,eps,number):
        '''
        - update the 'sup_rmsd' and 'dindex' lists with the 'number' domain.
        - actually "save" the domain and update it's boundaries with earlier
          found domains.
        - called by procedure 'partition'.
        '''

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
        '''
        - the main partitioning procedure, generates the domains.
        - stores the domains of molecules 1 and 2 in 'domindex' arrays. 
        '''

        # initialize
        self.sup_rmsd = [999]*self.mobFullCA.numAtoms()
        self.dindex = [-1]*self.mobFullCA.numAtoms()
        self.domain_count = 0
        seedradius = self.initrad

        # partition the protein
        while self.domain_count < self.ndomains:
            domain = self.convergence(eps, seedradius)

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

            # and reset the definitions of ca1 and ca2
            if len(mob_disjoint_l) > 0:
                self.mobCA = self.indexSel(self.mobile, mob_disjoint_l)
                self.refCA = self.indexSel(self.reference, ref_disjoint_l)

            if mob_domain.numAtoms() > 4:
                # convergence found a domain after normal iteration
                # reset seed radius and update domain boundaries
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

            # store the unsorted domains in the 'domindex' arrays
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

    def sortDomList(self,domlist):
        '''
        Sort domains by size (number of residues)
        '''

        i = 0
        tmplist = {}
        l = sorted(domlist.items(), key=lambda x: len(x[1]),reverse=True)
        for e in l:
            if e[0] != -1:
                tmplist[i] = e[1]
                i += 1

        return tmplist

    def sortDomains(self):
        mob_sorted = self.sortDomList(self.mobDomList)
        ref_sorted = self.sortDomList(self.refDomList)
        uncoverged = len(self.mobDomList[-1])
        converged = 0

        self.mobDomList = mob_sorted
        self.refDomList = ref_sorted
        dcount = 0

        for k in sorted(self.mobDomList.keys()):
            reslist = self.indexSel(self.mobile,self.mobDomList[k]).getResnums()
            print "sort> domain number %d - number of atoms %d - resids: %s"\
            %(k,len(reslist),self.arraystr(reslist))
            converged += len(reslist)
            dcount += 1

        self.domain_count = dcount

        print "sort> total number of domains larger than %d atoms: %d"%(self.cutdom,self.domain_count)
        print "sort> total number of residues involved: %d"%converged
        print "sort> uncoverged rest: %d"%uncoverged


    def calcHinge (self,RefDom=0, MobDom=1): 
        '''
        - generate and render the effective rotation axis of the movement 
          of 'domid2' (moving domain) relative to 'domid1' (reference
          domain).
        - called by the user.
        - (don't confuse reference/mobile domain with reference/mobile structure.)
        '''

        beforeR = self.indexSel(self.reference, self.refDomList[RefDom])
        afterR = self.indexSel(self.mobile, self.mobDomList[RefDom])
        beforeM = self.indexSel(self.reference, self.refDomList[MobDom])
        afterM = self.indexSel(self.mobile, self.mobDomList[MobDom])

        t = pd.calcTransformation(afterR,beforeR,weights=afterR.getMasses())
        pd.applyTransformation(t,self.mobile)

        com1 = pd.calcCenter(beforeM,beforeM.getMasses())
        com2 = pd.calcCenter(afterM,afterM.getMasses())
        cdis = np.linalg.norm(com2-com1)
                
        # calculate the best fit rmsd of the moving domain and reset
        trans_mat = pd.calcTransformation(afterM,beforeM,weights=afterR.getMasses())
        rotmat = trans_mat.getMatrix()
        pd.applyTransformation(trans_mat,self.mobile)
        rmsid = pd.calcRMSD(afterM,beforeM)
        t = pd.calcTransformation(afterR,beforeR)
        pd.applyTransformation(t,self.mobile)

        # calculate bisecting point and normalvector of bisecting plane
        bi = (com1 + com2) / 2
        pl = (com1 - com2) / np.abs(np.linalg.norm(com1 - com2))

        # compute rotation matrix for rotations about com2
        t21 = tf.translation_matrix(com2 - com1)
        rot2 = np.dot(rotmat,t21)

        # rotate com1 twice about com2
        # (arbitrary choice of rotation to generate new points to make a plane,
        # rotation axis is perpendicular to this?)
        p1 = np.dot(rot2,np.append(com2,1))
        p2 = np.dot(rot2,p1)

        # compute the rotation axis from the 3 points
        rideal =  np.cross(p2[:3] - com2, p1[:3] - com2)
        rideal /= np.linalg.norm(rideal)

        # project com2 onto the 3-point-plane
        new = com1 - (rideal * np.dot(rideal,(com1 - com2)))

        # compute rotation angle (least squares fit)
        ## print 2*np.arccos(tf.quaternion_from_matrix(rotmat)[0])/np.pi*180 
        angl = tf.rotation_from_matrix(rotmat)[0]
        angl_deg_o = angl * 180 / np.pi
        
        # compute projection of rot axis on bisecting plane
        perp = np.dot(rideal, pl)
        angp = np.abs(np.arcsin(perp))
        pro = (rideal - (pl*perp))/np.abs(np.linalg.norm(rideal - (pl*perp)))

        # compute decomposition angle
        tang = np.cos(angp) * np.tan(angl*0.5)
        angle = 2*np.arctan(tang)
        deg_angle = angle * 180 / np.pi

        # compute pivot point
        hi = bi + np.cross(pl, pro)*0.5*cdis/tang

        # # translate by effective rotation and reset
        # (rotate by $angle around $hi in the plane of $com1,$com2 and $hi?)
        t = pd.Transformation(tf.rotation_matrix(angle, np.cross((com2 - hi),(com1 - hi)),hi))
        pd.applyTransformation(t,self.mobile)
        rmspro = pd.calcRMSD(afterM,beforeM,weights=beforeM.getMasses()) 
        relerr = 100 *(rmspro - rmsid)/cdis
        deg_angp = angp*180/np.pi

        return [RefDom,MobDom,cdis,angl_deg_o,hi,pro,deg_angle,rmsid,rmspro,relerr,deg_angp,com1,com2]
        #       0      1      2    3          4  5   6         7     8      9      10       11   12

    def hingeOutput(self,refdom,mobdom,cdis,angl_deg_o,pivot,axis,eff_rot,rmsid,rmspro,relerr,proj):
        print "hinge> results:"
        print "hinge> Structures superimposed by domain %d (reference domain)"%refdom
        print "hinge> Distance moved of domain %d: %f Å"%(mobdom,cdis)
        print "Hinge> \"Overall\" Rotation: %f°"%angl_deg_o
        print "hinge> pivot point: %s"%self.arraystr(pivot)
        print "hinge> effective rotation axis: %s (left-handed rotation)"%self.arraystr(axis)
        print "hinge> effective rotation angle: %f°"%eff_rot
        print "hinge> accuracy:"
        print "hinge> rmsd (least squares): %f"%rmsid
        print "hinge> rmsd (effective rotation): %f"%rmspro
        print "hinge> relative error: %f %%"%relerr
        print "hinge> projection (deviation) angle: %f°"%proj

    def hingeTraj(self,domain1,domain2):
        '''
        Calculate hinge bending angle over trajectory

        '''
        overall = []
        effective_rot = []

        if self.structure and self.trajectory:
            trajfiles=self.trajectory.split(",")        

            self.mobile = pd.parsePDB(self.structure)
            print("Adding trajectory: %s\n"%trajfiles[0])
            self.trajectory = pd.Trajectory(trajfiles[0])
            for trj in trajfiles[1:]:
                print("Adding trajectory: %s\n"%trj)
                self.trajectory.addFile(trj)

            self.trajectory.link(self.mobile)

        ltraj_tenth = int(len(self.trajectory)/10)
        for i, frame in enumerate(self.trajectory):    
            prog = i%ltraj_tenth
            if prog == 0:
                print str(i/ltraj_tenth*10)+"%"

            results = self.calcHinge(domain1,domain2)
            overall.append([i,results[3]])
            effective_rot.append([i,results[6]])

        np.savetxt(self.output+".overall_rotation.dat",overall,delimiter=' ' )
        np.savetxt(self.output+".effective_rotation.dat",effective_rot,delimiter=' ' )

    def vmd_render(self,ref_domain,mob_domain,pivot,pro,com1,com2):
        '''
        Write .vmd script
        '''

        vmdfile = open(self.output+".vmd",'w')
        vmdfile.write(
'''
proc newrep {{mol top} {selection "all"} {style "newcartoon"} {color "chain"} {mat "AOChalky"}} {
    mol addrep $mol
    set lastrep [expr [molinfo $mol get numreps]-1]
    mol modselect $lastrep $mol $selection
    mol modstyle $lastrep $mol $style
    mol modcolor $lastrep $mol $color
    mol modmaterial $lastrep $mol $mat
    mol selupdate [lastrep] top on
}

mol new %s 
mol new %s

newrep 0 "protein" newcartoon element
newrep 1 "protein" newcartoon element BrushedMetal
'''%(self.refFileName,self.mobFileName)
            )


        colour = 0
        for d in sorted(self.mobDomList.keys()):
            indexstrM = self.arraystr(self.mobDomList[d])
            indexstrR = self.arraystr(self.refDomList[d])
            vmdfile.write(
'''
newrep 0 "same residue as index %s" newcartoon "ColorID %d"
set ref_dom_%d [atomselect 0 "index %s"]
newrep 1 "same residue as index %s" newcartoon "ColorID %d" BrushedMetal
set mob_dom_%d [atomselect 1 "index %s"]
'''%(indexstrR,colour,colour,indexstrR,\
    indexstrM,colour,colour,indexstrM)
                )
            colour += 1

        if type(ref_domain) is int:
            vmdfile.write('''
set fit [measure fit $mob_dom_%d $ref_dom_%d]
[atomselect 1 "all"] move $fit

'''%(ref_domain,ref_domain))

        com1 = self.arraystr(com1)
        com2 = self.arraystr(com2)
        pivot = self.arraystr(pivot)
        pro = self.arraystr(pro)

        vmdfile.write(
'''
set laxis 80 
set raxis 0.5

set ax1 [vecadd {%s} [vecscale {%s} [expr 0.50 * $laxis]]]
set ax2 [vecsub {%s} [vecscale {%s} [expr 0.42 * $laxis]]]
set ax3 [vecsub {%s} [vecscale {%s} [expr 0.50 * $laxis]]]
draw color [expr %d %% 17]
draw cylinder $ax1 $ax2 radius $raxis
draw cone $ax2 $ax3 radius [expr 2.5 * $raxis]
draw sphere $ax1 radius $raxis
draw cylinder {%s} {%s} radius $raxis
draw sphere {%s} radius $raxis
draw cylinder {%s} {%s} radius $raxis
draw sphere {%s} radius $raxis
'''%(pivot,pro,\
    pivot,pro,\
    pivot,pro,\
    mob_domain,\
    com1,pivot,\
    com1,\
    com2,pivot,\
    com2))
        
        vmdfile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--reference",type=str,help="Reference structure (default = 1lfh)",default="1lfh")
    parser.add_argument("-m","--mobile",type=str,help="Mobile structure (default = 1lfg)",default="1lfg")
    parser.add_argument("-d1","--domain1",type=str,help="Atom selection for domain 1 (default = domain 1 for 1lfg/h)",\
        default="resid 5 to 85 87 88 90 91 251 to 293 295 to 302 304 to 332 685 687 to 689 691 and name CA")
    parser.add_argument("-d2","--domain2",type=str,help="Atom selection for domain 2 (default = domain 2 for 1lfg/h)",\
        default="resid 92 to 100 104 to 140 143 to 218 220 to 250 and name CA")
    
    parser.add_argument("-t","--trajectory",type=str,help="Plot angle for trajectory",default=None)
    parser.add_argument("-s","--structure",type=str,help="Reference structure for trajectory",default=None)
    parser.add_argument("-e","--eps",type=float,help="Tolerance (default = 1.0 Å)",default=1.0)
    parser.add_argument("-c","--cutdom",type=int,help="Minimum number of atoms per domain (default = 10)",default=10)
    parser.add_argument("-o","--ouptut",type=str,help="Output file",default="hingefind_domains")

    args = parser.parse_args()
    kwargs = {}

    hf = HingeFind(args.mobile,args.reference,cutdom=args.cutdom,trajectory=args.trajectory,structure=args.structure)    

    if re.search("^\d+$", args.domain1) and re.search("^\d+$", args.domain2) :
        hf.partition(args.eps)
        hf.sortDomains()
        results = hf.calcHinge(int(args.domain1), int(args.domain2))
        hf.vmd_render(int(args.domain1),int(args.domain2),results[4],results[5],results[11],results[12])
        hf.hingeOutput(results[0],results[1],results[2],results[3],results[4],\
            results[5],results[6],results[7],results[8],results[9],results[10])

        if args.structure and args.trajectory:
            hf.hingeTraj(int(args.domain1), int(args.domain2))
    
    else: # if using predetermined vmd atom selections
        mobDom1 = hf.mobile.select(args.domain1).getIndices()
        mobDom2 = hf.mobile.select(args.domain2).getIndices()
        refDom1 = hf.reference.select(args.domain1).getIndices()
        refDom2 = hf.reference.select(args.domain2).getIndices()

        hf.mobDomList[0] = mobDom1
        hf.mobDomList[1] = mobDom2
        hf.refDomList[0] = refDom1
        hf.refDomList[1] = refDom2

        results = hf.calcHinge(0, 1)

        hf.hingeOutput(results[0],results[1],results[2],results[3],results[4],\
            results[5],results[6],results[7],results[8],results[9],results[10])

        hf.vmd_render(0,1,results[4],results[5],results[11],results[12])

        if args.structure and args.trajectory:
            hf.hingeTraj(0, 1)
        

if __name__ == '__main__':
    main()
