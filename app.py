from flask import Flask, render_template, request, jsonify
import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors  # Gets used to get info from RDKit even if not highlighted here
from rdkit.Chem import Draw
from rdkit.Chem import inchi
from rdkit.Chem import rdMolDescriptors
import math
import urllib.parse

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('main.html')

@app.route('/molecules', methods=['GET'])
def get_info():
    inpt_mol = request.args.get('molecule_id')
    property = request.args.get('property')

    if inpt_mol is None:
        return "Please enter a valid molecule ID."
    
    smiles = get_smiles(inpt_mol)

    if smiles:
        mol = Chem.MolFromSmiles(smiles)
    else:
        # Add caption if there is an image, else don't
        return render_template('molecule.html', image='', name=None, smiles=None, chembl_id=None, cid=None, \
                               formula=None, weight=None, logP=None, tpsa=None, mech=None, caption="", mol3D=None, \
                                WIP=True)
    
    cid, name = get_pubchem_name(smiles)
    chembl_id = get_chembl_id(smiles)

    # Put image in static folder so it's usable by Flask
    img_path = os.path.join('static', 'mol.png')
    Draw.MolToFile(mol, img_path)
    # Get info from RDKit
    mw = Chem.Descriptors.MolWt(mol)
    mw = str(math.trunc(mw * 100) / 100)  + " amu"  # Truncate to 2 decimal places
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    logP = Crippen.MolLogP(mol)
    logP = str(math.trunc(logP * 100) / 100)  # Round LogP to 2 decimal places
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    tpsa = str(math.trunc(tpsa * 100) / 100)  # Round TPSA to 2 decimal places
    mech = get_mechanism(chembl_id)
    mol3D = get3D(mol)

    # WIP code for getting similar molecules
    if property != None:
        print(f"Finding similar molecules of property {property}...")
        new_smiles = find_similar_molecules(mol, property, formula, smiles, chembl_id)
        new_names = []
        new_imgs = []
        for i in range(len(new_smiles)):
            new_cid, new_name = get_pubchem_name(new_smiles[i])
            new_names.append(new_name)
            new_mol = Chem.MolFromSmiles(new_smiles[i])
            new_img_path = os.path.join('static', f'mol{i}.png')
            new_imgs.append(f'mol{i}.png')
            Draw.MolToFile(new_mol, new_img_path)
        return render_template('expanded_molecule.html', image='mol.png', name=name, smiles=smiles, chembl_id=chembl_id, \
                               cid=cid, formula=formula, weight=mw, logP=logP, tpsa=tpsa, mech=mech, caption="go", \
                                mol3D=mol3D, new_smiles=new_smiles, new_names=new_names, new_imgs=new_imgs, \
                                    property=property, len=len(new_smiles))

    return render_template('molecule.html', image='mol.png', name=name, smiles=smiles, chembl_id=chembl_id, cid=cid, \
                           formula=formula, weight=mw, logP=logP, tpsa=tpsa, mech=mech, caption="go", \
                            mol3D=mol3D, WIP=True)

def get_smiles(inpt_mol):
    # Only spend time on ChEMBL ID if a ChEMBL ID is entered
    if len(inpt_mol) > 6 and inpt_mol[:6].lower() == "chembl":
        smiles = try_chembl(inpt_mol)
        if smiles: return smiles

    if try_smiles(inpt_mol):
        return inpt_mol

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
    
def cid_to_smiles(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Error fetching data from PubChem: {response.status_code}")
        return None

    try:
        return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except (KeyError, IndexError):
        print("SMILES not found for this CID.")
        return None

def get_chembl_id(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid golecule")
        return "Invalid molecule"
    
    inchiKey = inchi.MolToInchiKey(mol)
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={inchiKey}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    try:
        return response.json()['molecules'][0]['molecule_chembl_id']
    except (KeyError, IndexError):
        return "No ChEMBL ID found"

def get_pubchem_name(smiles):
    formatted_smiles = urllib.parse.quote(smiles)
    # Gets the IUPAC name by first getting the PubChem CID by using the SMILES code
    # and then uses the PubChem CID to get the IUPAC name
    cidURL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{formatted_smiles}/cids/JSON"
    cidResponse = requests.get(cidURL)

    if cidResponse.status_code != 200:
        print("Error fetching CID:", cidResponse.text)
        return None, "Name not found"

    try:
        cid = cidResponse.json()["IdentifierList"]["CID"][0]
    except (KeyError, IndexError):
        return None, "Name not found"

    nameURL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
    nameResponse = requests.get(nameURL)

    if nameResponse.status_code != 200:
        print("Error getting name:", nameResponse.status_code, nameResponse.text)
        return cid, "Unnamed compound"
    # Check if the response contains the expected data

    # Even with successful response, many output names are not IUPAC (i. e. aldehydes)
    try:
        name = nameResponse.json()["PropertyTable"]["Properties"][0]["IUPACName"]
        return cid, name.capitalize()
    except (KeyError, IndexError):
        return cid, "Unnamed compound"
    
def get_mechanism(id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={id}"
    response = requests.get(url)
    if response.status_code != 200:
        return "Mechanism retrieval failed"

    try:
        data = response.json()
        return data['mechanisms'][0]['mechanism_of_action']
    except (KeyError, IndexError):
        return "Not Pharmacological"

def get3D(mol):
    # Generate 3D coordinates for the molecule
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol)

    try:
        AllChem.UFFOptimizeMolecule(mol)
    except:
        print("UFF optimization failed")

    return Chem.MolToMolBlock(mol)

def find_similar_molecules(mol, property, formula, smiles, chembl_id):
    smiles_list = []
    chembl_list = []
    cid_list = []
    # Re-get values for most properties because they were truncated in main app.route() function
    if property == "formula":
        chembl_list = get_similar_formula(formula)
        # Use PubChem or ChEMBL to search by formula
    elif property == "structure":
        cid_list = get_similar_structure(smiles, 85, 15)
    elif property == "mw":
        print("Finding MW")
        mw = Chem.Descriptors.MolWt(mol)
        chembl_list = get_similar_mw(mw, 10)
        print(len(chembl_list))
    elif property == "logp":
        logp = Crippen.MolLogP(mol)
        chembl_list = get_similar_logp(logp, 0.3)
        # Query based on logP range
    elif property == "tpsa":
        tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
        chembl_list = get_similar_tpsa(mol, 10)
        # Query based on TPSA range
    elif property == "mechanism":
        chembl_list = get_similar_mechanism(chembl_id)
        
    if chembl_list != [] and smiles_list == []:
        print("Found ChEMBL IDs")
        for id in chembl_list:
            smile = try_chembl(id)
            if smile and smile != smiles: smiles_list.append(smile)
        return smiles_list[:8]  # Limit to first 8 similar molecules
    elif cid_list != [] and smiles_list == []:
        for id in cid_list:
            smile = cid_to_smiles(id)
            if smile and smile != smiles: smiles_list.append(smile)
        return smiles_list[:8]
    else:
        print("No similar molecules found")
        return ["No similar molecules found"]

def get_similar_formula(formula):
    print(f"Searching ChEMBL DB for formula: {formula}")

    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties.full_molformula={formula}"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"ChEMBL query failed with status {response.status_code}")
        return None

    try:
        data = response.json()
        chembl_ids = [entry['molecule_chembl_id'] for entry in data['molecules']]
        return chembl_ids
    except (KeyError, IndexError):
        return None
    
def get_similar_structure(smiles, min_similarity, cap):
    print(f"Searching PubChem DB for structure: {smiles}")
    formatted_smiles = urllib.parse.quote(smiles)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{formatted_smiles}/JSON?Threshold={min_similarity}&MaxRecords={cap}"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"PubChem similarity search failed: {response.status_code}")
        return None

    try:
        data = response.json()
        cids = data['IdentifierList']['CID']
        return cids  # PubChem CIDs
    except (KeyError, IndexError):
        return None
    
def get_similar_mw(mw, delta):
    print(f"Searching ChEMBL DB for MW: {mw}")
    min_mw = mw - delta
    max_mw = mw + delta

    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties.mw_freebase__gte={min_mw}&molecule_properties.mw_freebase__lte={max_mw}"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"ChEMBL MW search failed: {response.status_code}")
        return None

    try:
        data = response.json()
        chembl_ids = [entry['molecule_chembl_id'] for entry in data['molecules']]
        return chembl_ids
    except (KeyError, IndexError):
        return None

def get_similar_logp(logp, delta):
    print(f"Searching ChEMBL DB for LogP: {logp}")
    min_logp = logp - delta
    max_logp = logp + delta

    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties.alogp__gte={min_logp}&molecule_properties.alogp__lte={max_logp}"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"ChEMBL LogP search failed: {response.status_code}")
        return None

    try:
        data = response.json()
        chembl_ids = [entry['molecule_chembl_id'] for entry in data['molecules']]
        return chembl_ids
    except (KeyError, IndexError):
        return None
    
def get_similar_tpsa(mol, delta):
    print(f"Searching ChEMBL DB for TPSA: {mol}")
    tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
    min_tpsa = tpsa - delta
    max_tpsa = tpsa + delta

    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties.topologicalpolar_surfacearea__gte={min_tpsa}&molecule_properties.topologicalpolar_surfacearea__lte={max_tpsa}"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"ChEMBL TPSA search failed: {response.status_code}")
        return None

    try:
        data = response.json()
        chembl_ids = [entry['molecule_chembl_id'] for entry in data['molecules']]
        return chembl_ids
    except (KeyError, IndexError):
        return None

def get_similar_mechanism(chembl_id):
    print(f"Searching ChEMBL DB for mechanism of action: {chembl_id}")
    url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={chembl_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None

    data = response.json()['mechanisms'][0]['mechanism_of_action']
    mechanism = data
    # Find all molecules with same mechanism in ChEMBL
    search_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?mechanism_of_action={mechanism}"
    result = requests.get(search_url)

    if result.status_code != 200: return None
    else:
        sim_mols = [m['molecule_chembl_id'] for m in result.json().get('mechanisms', [])]
        return sim_mols # Limit to first 5 similar molecules
    
if __name__ == '__main__':
    app.run(host="127.0.0.1", port=8080, debug=True)