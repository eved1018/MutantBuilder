import argparse
import itertools
import os
import shutil

from modeller import *
from modeller.automodel import *
from modeller.automodel import autosched
from modeller.optimizers import ConjugateGradients, MolecularDynamics

"""
Modeller and python3 are needed to run. 

PARAMETERS:
    -pdb: PDB file 
    -c: Chain to mutate
    -f: Config file containing which positions to mutate (defualt 'config.txt')
    -o: Output folder (defualt 'mutants')
    -r: add flag to not use homology modeling and use restricted model instead 

Example use: 
python mutantBuilder.py -pdb 6waq.pdb -c A

CONFIG FILE:
- The config.txt is used to guide the mutagensis each line contains a residue postion and amino acid to mutate to. 
- A '*' can be used to mutate to all resdiues. 
- Multiple mutants can be done at the same time using a comma as a seperator. (on the same line)

        EX 1: 2 seperate single mutants. 
        1 ALA
        4 GLN 

        EX 2: 1 double mutant.
        1 ALA, 4 GLN 

        EX 3: All mutants at postion 5
        5 *
"""

def userinterface():
    """
        User Interface:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", "--PDBFileName", help = "Name of saved PDB file.", required = True)
    parser.add_argument("-c", "--Chain", help = "Chain indicator", required = True)
    parser.add_argument("-f", "--File", help = "File containing mutant infromation see README for details ", default="config.txt")
    parser.add_argument("-o", "--OutFolder", help = "output folder for mutants", default="mutants")
    parser.add_argument("-r", "--Userestricted", help = "Dont use homology modeling", action='store_true', default=False)


    args = parser.parse_args()
    pdb = args.PDBFileName
    chain = args.Chain
    outfolder = args.OutFolder
    Userestricted = args.Userestricted
    os.makedirs(outfolder, exist_ok=True)
    os.makedirs("tmp", exist_ok=True)

    config_file = args.File

    mcontainer = MutantContainer()
    mcontainer.readFromFile(config_file)
    
    return pdb, mcontainer, chain, Userestricted, outfolder

class MutantContainer:
    """
        Use to hold mutants, parse config file 
    """
    def __init__(self) -> None:
        self.mutations = []
        return

    def addMutation(self, mutant):
        self.mutations.append(mutant)
        return

    def getMutations(self):
        return self.mutations

    def readFromFile(self, filename):
        with open(filename, "r") as fh:
            for line in fh.readlines():
                if "*" in line:
                    subs = line.split(",")
                    positions = [i.split()[0] for i in subs]
                    combos = itertools.product(self.allAminoAcids(), repeat = len(positions))
                    for i in combos:
                        mutant = Mutant()
                        for pos, aa in zip(positions, i):
                            mutant.addMutant(pos, aa)
                        self.mutations.append(mutant)
                else:
                    mutant = Mutant()
                    for submutant in line.split(","):
                        position, mutantAA  = submutant.split()
                        mutantAA = self.mutationConfiger(mutantAA)
                        mutant.addMutant(position, mutantAA)
                    self.mutations.append(mutant)
        
        return
    
    def allAminoAcids(self):
        aa3 = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG', 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']
        return aa3

    def mutationConfiger(self, mutantAA):
        three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
                    'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
                    'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
                    'G':'GLY', 'P':'PRO', 'C':'CYS'}
        aa3 = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG', 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']
        aa1 = ['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
        
        mutantAA =  mutantAA.upper()

        if len(mutantAA) == 3 and  mutantAA in aa3:
            return mutantAA

        elif len(mutantAA) == 1 and mutantAA in aa1:
            return three_letter[mutantAA]

        else:
            print(f"Unrecognizable Amino Acid in input: {mutantAA}")
            return "nan"
            
class Mutant:
    """
        define object for each round of mutation, ie single mutant, double mutant 
    """
    def __init__(self) -> None:
        self.mutants = []
        return 

    def addMutant(self, position, mutantAA):
        self.mutants.append([position, mutantAA])
        return

    def getMutant(self):
        return self.mutants

####################################################################################

# Homology Modeling 
def HomologyModel(pdb, mutant, chain, mutant_name, outfolder):
    """
        build homology model and refine using automodel 
    """
    log.none()
    # Set Modeller environment (including search patch for model.read())
    env = Environ()
    env.io.atom_files_directory = "./:../atom_files/"

     # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')


    # Create a new empty alignment and model:
    aln = Alignment(env)
    mdl = Model(env)

    mdl.read(file=pdb, model_segment=('FIRST:@', 'END:'))

    aln.append_model(mdl, align_codes=pdb)
    
    s2 =  Selection()
    #set up the mutate residue selection segment
    for respos, restype in mutant:
        s = Selection()
        s.add(mdl.chains[chain].residues[respos])
        s.mutate(residue_type=restype)
        s2.add(mdl.chains[chain].residues[respos])
    
    aln.append_model(mdl, align_codes="target")
    # aln.write(file="wt_ali.ali")

    # Align them by sequence
    aln.malign(gap_penalties_1d=(-500, -300))
    # aln.write(file= tmp_file_name )

    # check the alignment for its suitability for modeling
    aln.check()
    # aln.write(file= tmp_file_name + ".ali" )
    mutant_pdb = autoModelPDB(pdb, aln )
    # mutant_pdb = "tmp/" + mutant_pdb
    shutil.move( "tmp/" + mutant_pdb, f"{outfolder}/homolog-" + mutant_name )
    return

def autoModelPDB(pdb, ali_file ):
    """
        Automodel the homology model
    """
    env = environ(rand_seed = -8123)
    log.none()
    env.io.hetatm = False
    os.chdir("tmp/")
    a = automodel(env, alnfile=ali_file, knowns= pdb ,sequence= "target", assess_methods=(assess.DOPE))
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    os.chdir("../")
    return a.outputs[0]['name']

####################################################################################

# Restricted Modeling 
def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False

def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = Selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)

def restrictedModel(modelname, mutant, chain, mutant_name, outfolder):
#first argument
    # log.verbose()
    log.none()

# Set a different value for rand_seed to get a different final model
    env = Environ(rand_seed=-49837)

    env.io.hetatm = True
#soft sphere potential
    env.edat.dynamic_sphere=False
#lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

# Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')


# Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = Model(env, file=modelname)
    ali = Alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

# #set up the mutate residue selection segment
#     s = Selection(mdl1.chains[chain].residues[respos])

# #perform the mutate residue operation
#     s.mutate(residue_type=restyp)

    s2 =  Selection()
    #set up the mutate residue selection segment
    for respos, restype in mutant:
        s = Selection()
        # print(respos,restype)
        s.add(mdl1.chains[chain].residues[respos])
        s.mutate(residue_type=restype)
        s2.add(mdl1.chains[chain].residues[respos])
    

#get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=modelname)

# Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])


# Transfer all the coordinates you can from the template native structure
# to the mutant (this works even if the order of atoms in the native PDB
# file is not standard):
#here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

# Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

#yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = Model(env, file=modelname)

    #required to do a transfer_res_numb
    #ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
    #transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2,ali)

#It is usually necessary to write the mutated sequence out and read it in
#before proceeding, because not all sequence related information about MODEL
#is changed by this command (e.g., internal coordinates, charges, and atom
#types and radii are not updated).

    mdl1.write(file=outfolder + "/" + mutant_name+'.tmp')
    mdl1.read(file=outfolder + "/" + mutant_name+'.tmp')

#set up restraints before computing energy
#we do this a second time because the model has been written out and read in,
#clearing the previously set restraints
    make_restraints(mdl1, ali)

#a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms=1

    sched = autosched.loop.make_for_model(mdl1)

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s2.energy()

    s2.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)

#feels environment (energy computed on pairs that have at least one member
#in the selected)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)

    s2.energy()

#give a proper name
    mdl1.write(file=outfolder + "/" + mutant_name)

#delete the temporary file
    os.remove(outfolder + "/"+ mutant_name+'.tmp')
    return 

####################################################################################

def main():
    
    pdb, mcontainer, chain, Userestricted, outfolder = userinterface()
    
    mutations  = mcontainer.getMutations()
    model  = restrictedModel if Userestricted else HomologyModel

    for mutants in mutations:
        mutants  = mutants.getMutant()
        mutant_name =  "_".join([j for i in mutants for j in i]) + ".pdb"
        model(pdb, mutants, chain, mutant_name, outfolder)
    return
    
if __name__ == '__main__':
    main()
