from pyteomics import fasta, mass

def calculate_mass_charge(faa):
    with fasta.read(faa) as f:
        for record in f:
            sequence = record[1]
            
            mw = mass.calculate_mass(sequence=sequence, monoisotopic=True)
            print(f"Molecular Weight: {mw} Da")
