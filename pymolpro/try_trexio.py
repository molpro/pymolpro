import pymolpro
p=pymolpro.Project("Zn")
p.write_input('geometry={Zn};df-hf;show,energy')
p.run(wait=True)
print(p.out)
print(p.variable('ENERGY'))
p.orbitals_to_molden('thing',0)
p.orbitals_to_trexio('thing',0)