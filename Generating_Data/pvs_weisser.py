import os, sys, pythia8, ROOT
import numpy as np
pre = '/Users/weisser/MIT_Dropbox/LbVMWeisser_shared/Tracking/Simulated_Velo/fastsim'
save_pre = '/Users/weisser/MIT_Dropbox/LbVMWeisser_shared/Tracking/Simulated_Velo/LHCbPVFinding_DataSets'
save_path = save_pre + '/Data_ROOT/'
cwd = os.getcwd()
sys.path.insert(0, pre)
os.chdir(pre)
import velo
os.chdir(cwd)

# Writer class.
class Writer():
    def __init__(self):
        from collections import OrderedDict
        self.vars = OrderedDict()
        self.null = ROOT.vector('double')(1, 0)
    def init(self, tree):
        for key, val in self.vars.iteritems(): tree.Branch(key, val)
    def add(self, var): self.vars[var] = ROOT.vector('double')()
    def var(self, var, val = None, idx = -2):
        if not var in self.vars: return self.null.back()
        var = self.vars[var]
        if idx < -1: var.push_back(0 if val == None else val)
        if idx < 0: idx = var.size() - 1
        elif idx >= var.size(): idx = -1
        if idx < 0: return self.null[0]
        if val != None: var[idx] = val
        return var[idx]
    def size(self, var): return self.vars[var].size()
    def clear(self):
        for key, val in self.vars.iteritems(): val.clear()

def Hits(mtr, prt):
    hits = []
    p = prt.pAbs()
    if p == 0: return hits
    vx, vy, vz = prt.xProd(), prt.yProd(), prt.zProd(),
    px, py, pz = prt.px()/p, prt.py()/p, prt.pz()/p
    hit = mtr.intersect(vx, vy, vz, px, py, pz)
    while hit.T() != 0:
        vx, vy, vz = hit.X(), hit.Y(), hit.Z()
        hits += [[vx, vy, vz]]
        vx, vy, vz = vx + px*0.1, vy + py*0.1, vz + pz*0.1
        hit = mtr.intersect(vx, vy, vz, px, py, pz)
    return hits

# Initialize Pythia.
random = ROOT.TRandom3()
pythia = pythia8.Pythia('', False)
pythia.readString('Print:quiet = on')
pythia.readString('SoftQCD:all = on')
pythia.init()
material = velo.ModuleMaterial('/Users/weisser/MIT_Dropbox/LbVMWeisser_shared/'
                               'Tracking/Simulated_Velo/fastsim/dat/run3.root')

# Create the output TFile and TTree.
tfile = ROOT.TFile(save_path+'pvs_weisser.root', 'RECREATE')
#tfile = ROOT.TFile(save_path+'pvs_weisser_test2.root', 'RECREATE')
#tfile = ROOT.TFile(save_path+'pvs_weisser_easy.root', 'RECREATE')
#tfile = ROOT.TFile(save_path+'pvs_weisser_1M.root', 'RECREATE')
ttree = ROOT.TTree('data', '')

# Create the writer handler, add branches, and initialize.
writer = Writer()
writer.add('pvr_x')
writer.add('pvr_y')
writer.add('pvr_z')
writer.add('hit_x')
writer.add('hit_y')
writer.add('hit_z')
writer.add('hit_prt')
writer.add('prt_pid')
writer.add('prt_px')
writer.add('prt_py')
writer.add('prt_pz')
writer.add('prt_e')
writer.add('prt_x')
writer.add('prt_y')
writer.add('prt_z')
writer.init(ttree)

number_rejected_events = 0

# Fill the events.
iEvt, tEvt = 0, 1e4
#iEvt, tEvt = 0, 1e6
while iEvt < tEvt:
    if not pythia.next(): continue
    else: iEvt += 1; writer.clear()
    if (iEvt%1000 ==0): print "Event : ", iEvt/1000, "k / ", tEvt/1000, "k"
    # All distance measurements are in units of mm
    #xPv, yPv, zPv = 0, 0, random.Gaus(100, 63) # normal LHCb operation without pv x and y spread
    xPv, yPv, zPv = random.Gaus(0, 0.055), random.Gaus(0, 0.055), random.Gaus(100, 63) # normal LHCb operation
    #xPv, yPv, zPv = 0, 0, np.random.choice([0, 200], 1, p=[0.5, 0.5])[0]

    #pvr x and y spead can be found https://arxiv.org/pdf/1410.0149.pdf page 42. z dependent
    [-1000,-750, -500, -250] # mm  

    writer.var('pvr_x', xPv)
    writer.var('pvr_y', yPv) 
    writer.var('pvr_z', zPv) 
    number_of_detected_particles = 0
    for prt in pythia.event:
        if not prt.isFinal or prt.charge() == 0: continue
        prt.xProd(prt.xProd() + xPv) # Need to change the origin of the event before getting the hits
        prt.yProd(prt.yProd() + yPv)
        prt.zProd(prt.zProd() + zPv)
        hits = Hits(material, prt)
        if len(hits) == 0: continue
        if len(hits) > 2: number_of_detected_particles += 1
        #prt.xProd(prt.xProd() + xPv)
        #prt.yProd(prt.yProd() + yPv)
        #prt.zProd(prt.zProd() + zPv)
        writer.var('prt_pid', prt.id()) 
        writer.var('prt_px',  prt.px()) 
        writer.var('prt_py',  prt.py()) 
        writer.var('prt_pz',  prt.pz()) 
        writer.var('prt_e',   prt.e()) 
        writer.var('prt_x',   prt.xProd()) 
        writer.var('prt_y',   prt.yProd()) 
        writer.var('prt_z',   prt.zProd()) 
        for xHit, yHit, zHit in hits:
            #xHit_recorded, yHit_recorded, zHit_recorded = random.Gaus(xHit, 0.025), random.Gaus(yHit, 0.025), random.Gaus(zHit, 0.025) # _old_smearing
            xHit_recorded, yHit_recorded, zHit_recorded = np.random.uniform(-0.0275,0.0275)+xHit, np.random.uniform(-0.0275,0.0275)+yHit, zHit # normal
            #xHit_recorded, yHit_recorded, zHit_recorded = np.random.uniform(-0.0275,0.0275)+xHit, np.random.uniform(-0.0275,0.0275)+yHit, random.Gaus(zHit, 0.025) #test1
            #xHit_recorded, yHit_recorded, zHit_recorded = random.Gaus(xHit, 0.025), random.Gaus(yHit, 0.025), random.Gaus(zHit, 0.025) #test2 

            #writer.var('hit_x', xHit) 
            #writer.var('hit_y', yHit) 
            #writer.var('hit_z', zHit)
            writer.var('hit_x', xHit_recorded) 
            writer.var('hit_y', yHit_recorded) 
            writer.var('hit_z', zHit_recorded) 
            writer.var('hit_prt', writer.size('prt_e') - 1)
    if number_of_detected_particles < 5: iEvt -= 1; number_rejected_events+=1; continue
    # Fill the TTree.
    ttree.Fill()

print "Number of rejected events : ", number_rejected_events, "\t for ", tEvt, " accepted"
  
# Write and close the TTree and TFile.
ttree.Write(ttree.GetName(), ROOT.TObject.kOverwrite)
tfile.Close()
