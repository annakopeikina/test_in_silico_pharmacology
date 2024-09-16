import sys
import os
import requests
import shutil
import io

if hasattr(sys.stdout, 'buffer'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def search_by_trade_name(trade_name):
    api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{trade_name}/property/MolecularFormula,InChIKey,Title/JSON"
    response = requests.get(api_url)
    if response.status_code == 200:
        data = response.json()
        if 'PropertyTable' in data:
            properties = data['PropertyTable']['Properties'][0]
            return properties.get("Title"), properties.get("MolecularFormula"), properties.get("InChIKey")
    return None, None, None

def search_drug(drug_name):
    latin_name, brutto_formula, chemical_name = search_by_trade_name(drug_name)
    if latin_name:
        return latin_name, brutto_formula, chemical_name
    return None, None, None

def main():
    if len(sys.argv) != 2:
        print("Использование: python search.py <drug_name>")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    latin_name, brutto_formula, chemical_name = search_drug(drug_name)
    
    if latin_name:
        project_root = r"C:\Users\annav\OneDrive\Desktop\test_in_silico_pharmacology"
        drug_folder = os.path.join(project_root, "drug_data", drug_name.replace(" ", "_"))
        
        if not os.path.exists(drug_folder):
            os.makedirs(drug_folder)

        with open(os.path.join(drug_folder, "info.txt"), "w", encoding="utf-8") as info_file:
            info_file.write(f"Латинское название: {latin_name}\n")
            info_file.write(f"Брутто формула: {brutto_formula}\n" if brutto_formula else "Брутто формула не найдена.\n")
            info_file.write(f"Химическое название (InChIKey): {chemical_name}\n" if chemical_name else "Химическое название не найдено.\n")
        
        print(drug_folder)
    else:
        print("Препарат не найден.")
        sys.exit(1)

if __name__ == "__main__":
    main()


# def search_by_trade_name(trade_name):
    # '''Функция принимает торговое название препарата trade_name в качестве аргумента,
    # отправляет запрос к API для поиска информации о препарате по его торговому названию.
    # Если препарат найден, функция возвращает его латинское название, брутто формулу и 
    # химическое название. В случае, если препарат не найден, функция возвращает None 
    # для каждого из трех параметров.'''
    # api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/search?trade_name={trade_name}"
    # response = requests.get(api_url)
    # if response.status_code == 200:
        # data = response.json()
        # if data.get("found"):
            # return data["latin_name"], data.get("brutto_formula"), data.get("chemical_name"
    # return None, None, None
# 
# 
# def search_by_latin_name(latin_name):
    # '''Функция принимает латинское название препарата latin_name в качестве аргумента и 
    # отправляет запрос к API для поиска информации о препарате по его латинскому названию. 
    # Если препарат найден, функция возвращает его брутто формулу и химическое название. 
    # Если препарат не найден, функция возвращает None для обоих параметров.'''
    # api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/search?latin_name={latin_name}"
    # response = requests.get(api_url)
    # if response.status_code == 200:
        # data = response.json()
        # if data.get("found"):
            # return data.get("brutto_formula"), data.get("chemical_name")
    # return None, None
# 
# 
# def search_by_brutto_formula(brutto_formula):
    # '''Функция принимает брутто формулу препарата brutto_formula в качестве аргумента и 
    # отправляет запрос к API для поиска информации о препарате по его брутто формуле. 
    # Если препарат найден, функция возвращает его химическое название. Если препарат не найд
    # функция возвращает None.'''
    # api_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/search?brutto_formula={brutto_formula}"
    # response = requests.get(api_url)
    # if response.status_code == 200:
        # data = response.json()
        # if data.get("found"):
            # return data.get("chemical_name")
    # return None