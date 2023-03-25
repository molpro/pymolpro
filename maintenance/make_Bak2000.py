import pymolpro
import os

molecules = {}
molecules["F2"] = {"De_exp": 163.35, "De_calc": 161.04, "geometry": "F;F,F,1.41134168", "stoichiometry": "FF"}
molecules["H2"] = {"De_exp": 458.04, "De_calc": 458.13, "geometry": "H;H,H,0.74184515", "stoichiometry": "HH"}
molecules["HF"] = {"De_exp": 593.16, "De_calc": 593.20, "geometry": "H;F,H,0.91576748", "stoichiometry": "HF"}
molecules["CH2"] = {"De_exp": 757.06, "De_calc": 575.87, "geometry": "C;H,C,1.10678917;H,C,1.10678917,H,102.02249921",
                    "stoichiometry": "CHH"}
molecules["HNO"] = {"De_exp": 861.50, "De_calc": 860.12, "geometry": "H;N,H,1.05238269;O,N,1.20847128,H,108.08480747",
                    "stoichiometry": "HNO"}
molecules["N2"] = {"De_exp": 956.28, "De_calc": 954.68, "geometry": "N;N,N,1.09811947", "stoichiometry": "NN"}
molecules["H2O"] = {"De_exp": 975.28, "De_calc": 975.32, "geometry": "H;O,H,0.95711647;H,O,0.95711647,H,104.22359713",
                    "stoichiometry": "HHO"}
molecules["CO"] = {"De_exp": 1086.70, "De_calc": 1086.70, "geometry": "C;O,C,1.12887680", "stoichiometry": "CO"}
molecules["H2O2"] = {"De_exp": 1126.91, "De_calc": 1125.82,
                     "geometry": "H1;O1,H1,0.9619;O2,O1,1.4497,H1,99.99;H2,O2,0.9619,O1,99.99,H1,112.7",
                     "stoichiometry": "HHOO"}
molecules["NH3"] = {"De_exp": 1247.88, "De_calc": 1247.15,
                    "geometry": "N1;X,N1,10;H1,N1,1.01116937,X,112.42635522;H2,N1,1.01116937,X,112.42635522,H1,120;H3,N1,1.01116937,X,112.42635522,H2,120",
                    "stoichiometry": "NHHH"}
molecules["HCN"] = {"De_exp": 1312.75, "De_calc": 1310.72, "geometry": "H;C,H,1.06553056;N,C,1.15383982,H,180",
                    "stoichiometry": "HCN"}
molecules["CH2O"] = {"De_exp": 1566.58, "De_calc": 1567.78,
                     "geometry": "C;O,C,1.20430342;H1,C,1.10077162,O,121.78046616;H2,C,1.10077162,O,121.78046616,H1,180",
                     "stoichiometry": "CHHO"}
molecules["CO2"] = {"De_exp": 1632.46, "De_calc": 1632.77, "geometry": "O1;C,O1,1.16039505;O2,C,1.16039505,O1,180",
                    "stoichiometry": "COO"}
molecules["C2H2"] = {"De_exp": 1697.84, "De_calc": 1696.81,
                     "geometry": "H1;C1,H1,1.06209557;C2,C1,1.20366063,H1,180;H2,C2,1.06209557,C1,180,H1,0",
                     "stoichiometry": "CCHH"}
molecules["CH4"] = {"De_exp": 1759.33, "De_calc": 1759.30,
                    "geometry": "C1;H4,C1,1.08637390;H1,C1,1.08637390,H4,109.47122063;H2,C1,1.08637390,H4,109.47122063,H1,120;H3,C1,1.08637390,H4,109.47122063,H2,120",
                    "stoichiometry": "CHHHH"}
molecules["C2H4"] = {"De_exp": 2359.82, "De_calc": 2360.55,
                     "geometry": "C1;C2,C1,1.33117445;H11,C1,1.08089487,C2,121.44926360;H12,C1,1.08089487,C2,121.44926360,H11,180;H21,C2,1.08089487,C1,121.44926360,H11,0;H22,C2,1.08089487,C1,121.44926360,H12,0",
                     "stoichiometry": "CCHHHH"}

reactions = {}
reactions["1"] = {"CO": -1, "H2": -1, "CH2O": +1}
reactions["2"] = {"N2": -1, "H2": -3, "NH3": +2}
reactions["3"] = {"C2H2": -1, "H2": -1, "C2H4": +1}
reactions["4"] = {"CO2": -1, "H2": -4, "CH4": +1, "H2O": +2}
reactions["5"] = {"CH2O": -1, "H2": -2, "CH4": +1, "H2O": +1}
reactions["6"] = {"CO": -1, "H2": -3, "CH4": +1, "H2O": +1}
reactions["7"] = {"HCN": -1, "H2": -3, "CH4": +1, "NH3": +1}
reactions["8"] = {"H2O2": -1, "H2": -1, "H2O": +2}
reactions["9"] = {"HNO": -1, "H2": -2, "H2O": +1, "NH3": +1}
reactions["10"] = {"C2H2": -1, "H2": -3, "CH4": +2}
reactions["11"] = {"CH2": -1, "H2": -1, "CH4": +1}
reactions["12"] = {"F2": -1, "H2": -1, "HF": +2}
reactions["13"] = {"CH2": -2, "C2H4": +1}

reactions_db = pymolpro.database.Database()
for name, molecule in molecules.items():
    reactions_db.add_molecule(name, molecule['geometry'],
                              reference_energy=molecule['De_exp'] / 2625.49963948, description=name)

atomisations_db = reactions_db

for name, reaction in reactions.items():
    reactions_db.add_reaction(name, reaction, description=name)

reactions_db.dump(
    os.path.realpath(os.path.join(__file__, '..', '..', 'share', 'database', 'Bak2000_reactions' + '.json')))

for element in ['H', 'C', 'N', 'O', 'F']:
    atomisations_db.add_molecule(element, element, reference_energy=0.0, description=element)

for name in molecules:
    stoi = {name: -1}
    string = molecules[name]['stoichiometry']
    while string != "":
        if string[0] in stoi:
            stoi[string[0]] += 1  # assumes only 1-character element names
        else:
            stoi[string[0]] = 1
        string = string[1:]
    atomisations_db.add_reaction(name, stoichiometry=stoi, description=name)
atomisations_db.dump(
    os.path.realpath(os.path.join(__file__, '..', '..', 'share', 'database', 'Bak2000_atomisations' + '.json')))
