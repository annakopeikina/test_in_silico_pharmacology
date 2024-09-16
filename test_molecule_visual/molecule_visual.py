import requests
from rdkit import Chem
from rdkit.Chem import Draw
import os
import sys
import io

# Устанавливаем кодировку UTF-8 для вывода
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def fetch_smiles(drug_name):
    # Здесь можно использовать API для получения SMILES по названию препарата
    # Пример запроса к PubChem API
    api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/CanonicalSMILES/JSON"
    response = requests.get(api_url)
    if response.status_code == 200:
        data = response.json()
        if 'PropertyTable' in data:
            properties = data['PropertyTable']['Properties'][0]
            return properties.get("CanonicalSMILES")
    return None

def visualize_molecule(smiles, output_path):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    img.save(output_path)

def main():
    if len(sys.argv) != 3:
        print("Использование: python molecule_visual.py <drug_name> <folder_path>")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    folder_path = sys.argv[2]

    smiles = fetch_smiles(drug_name)
    if smiles:
        output_path = os.path.join(folder_path, "molecule.png")
        visualize_molecule(smiles, output_path)
        print(f"Изображение молекулы сохранено в {output_path}")
    else:
        print("SMILES для данного препарата не найден.")
        sys.exit(1)

if __name__ == "__main__":
    main()
