from flask import Flask, render_template, request, jsonify
import os
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors  # Gets used to get info from RDKit even if not highlighted here
from rdkit.Chem import Draw
import math

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('main.html')

@app.route('/molecules', methods=['GET'])
def get_info():
    inpt_mol = request.args.get('molecule_id')

    if inpt_mol is None:
        return "Please enter a valid molecule ID."
    
    smiles = get_smiles(inpt_mol)

    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        name = get_pubchem_name(smiles)
        id = get_chembl_id(smiles)
    else:
        # Add caption if there is an image, else don't
        return render_template('molecule.html', image='', name=None, id=None, formula=None, weight=None, caption="")
    
    if name.lower() == "formaldehyde": name = "Methanal"  # PubChem returns Formaldehyde for IUPAC name
    elif name.lower() == "acetaldehyde": name = "Ethanal" # PubChem returns Acetaldehyde for IUPAC name
    elif name.lower() == "benzaldehyde": name = "Benzenecarbaldehyde"  # PubChem returns Benzaldehyde for IUPAC name

    # Put image in static folder so it's usable by Flask
    img_path = os.path.join('static', 'mol.png')
    Draw.MolToFile(mol, img_path)
    # Get info from RDKit
    mw = Chem.Descriptors.MolWt(mol)
    mw = str(math.trunc(mw * 100) / 100)  + " amu"  # Truncate to 2 decimal places
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

    return render_template('molecule.html', image='mol.png', name=name, id=id, formula=formula, weight=mw, caption="RDKit")

def get_smiles(inpt_mol):
    # Only spend time on ChEMBL ID if a ChEMBL ID is entered
    if len(inpt_mol) > 6 and inpt_mol[:6].lower() == "chembl":
        print("Trying ChEMBL ID...")
        smiles = try_chembl(inpt_mol)
        if smiles: return smiles

    print("Trying SMILES code...")
    if try_smiles(inpt_mol):
        return inpt_mol

    print("Trying IUPAC name...")
    smiles = try_iupac(inpt_mol)
    if smiles: return smiles

    return None
    
def try_smiles(code):
    mol = Chem.MolFromSmiles(code)
    return mol != None

def try_iupac(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    try:
        return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except (KeyError, IndexError):
        return None

def try_chembl(id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{id}.json"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Error fetching data from ChEMBL: {response.status_code}")
        return None

    try:
        data = response.json()
        smiles = data['molecule_structures']['canonical_smiles']
        return smiles
    except (KeyError, TypeError):
        print("SMILES not found for this ChEMBL ID.")
        return None

def get_chembl_id(smiles):
    # Go to EBI database to get ChEMBL  data through the SMILES code
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?smiles={smiles}"
    response = requests.get(url)

    # Check if the response is successful, much better success rate when user enters a ChEMBL ID
    if response.status_code != 200:
        print(f"Error fetching data from ChEMBL: {response.status_code}")
        return None

    try:
        data = response.json()
        chembl_id = data['molecules'][0]['molecule_chembl_id']
        return chembl_id
    except (KeyError, IndexError):
        return "No ChEMBL ID found"

def get_pubchem_name(smiles):
    # Gets the IUPAC name by first getting the PubChem CID by using the SMILES code
    # and then uses the PubChem CID to get the IUPAC name
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    response = requests.get(url)

    if response.status_code != 200:
        print("Error fetching CID:", response.status_code, response.text)
        return "No PubChem SMILES code found"

    try:
        cid = response.json()["IdentifierList"]["CID"][0]
    except (KeyError, IndexError):
        return "No PubChem CID found"

    url_name = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
    response_name = requests.get(url_name)

    if response_name.status_code != 200:
        print("Error getting name:", response_name.status_code, response_name.text)
        return "Name not available"

    # Even with successful response, many output names are not IUPAC (i. e. aldehydes)
    try:
        name = response_name.json()["PropertyTable"]["Properties"][0]["IUPACName"]
        return name.capitalize()
    except (KeyError, IndexError):
        return "Name not found"
    
if __name__ == '__main__':
    app.run(debug=True)