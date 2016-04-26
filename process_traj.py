import mdtraj
import numpy

tt_trr = mdtraj.load( "dna_models/dna17bp_sol_neu-MD_1.trr", top="dna_models/dna17bp_sol_neu-MD_1.pdb")
atoms_to_keep = []

vs = 0
hh = 0
na = 0
for atom in tt_trr.topology.atoms:
    e = atom.element.symbol
    if e in ['O','P','C','N']:
        atoms_to_keep.append(atom.index)
    elif e in ['VS']:
        vs +=1
    elif e in ["NA"]:
        na+=1
    else:
        hh+=1

        
print tt_trr
print "Found %d virtual sites"%vs
print "Found %d hydrogen sites"%hh
print "Found %d sodium sites"%na
print "atoms kept %d"%len(atoms_to_keep)
print atoms_to_keep[0:30]

tt_trr.atom_slice(atoms_to_keep,inplace=True)


print tt_trr
output = "dna_models/dna17bp_sol_neu-MD_1_sliced.trr"
top = "dna_models/dna17bp_sol_neu-MD_1_sliced.pdb"

tt_trr.save_trr(output)
tt_trr[0].save_pdb(top)

tt_sliced=mdtraj.load("dna_models/dna17bp_sol_neu-MD_1_sliced.trr",top="dna_models/dna17bp_sol_neu-MD_1_sliced.pdb")
print tt_sliced
