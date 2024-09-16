import xml.etree.ElementTree as ET
import os
import sys
import requests
import io

# Устанавливаем кодировку UTF-8 для вывода
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_chembl_id(drug_name):
    search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?pref_name__iexact={drug_name}"
    response = requests.get(search_url)

    if response.status_code == 200:
        try:
            root = ET.fromstring(response.content)
            molecules = root.findall(".//molecule")
            if molecules:
                drug_id = molecules[0].find("molecule_chembl_id").text
                return drug_id
            else:
                print("Препарат не найден.")
                return None
        except ET.ParseError:
            print("Ошибка парсинга XML. Ответ сервера:")
            print(response.text)
            return None
    else:
        print(f"Ошибка при запросе препарата: {response.status_code}")
        return None

def get_drug_mechanism(drug_id):
    mechanism_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id={drug_id}"
    response = requests.get(mechanism_url)

    if response.status_code == 200:
        try:
            return ET.fromstring(response.content)
        except ET.ParseError:
            print("Ошибка парсинга XML.")
            return None
    else:
        print(f"Ошибка при запросе информации о механизме действия: {response.status_code}")
        return None

def save_readable_info(drug_info, file_path):
    with open(file_path, "w", encoding="utf-8") as file:
        for element in drug_info.iter():
            if element.text:
                file.write(f"{element.tag}: {element.text.strip()}\n")

def main():
    if len(sys.argv) != 2:
        print("Использование: python pharmacodynamics.py <drug_name>")
        sys.exit(1)

    drug_name = sys.argv[1]
    drug_id = get_chembl_id(drug_name)

    if drug_id:
        drug_mechanism = get_drug_mechanism(drug_id)

        if drug_mechanism:
            project_root = r"C:\Users\annav\OneDrive\Desktop\test_in_silico_pharmacology"
            drug_folder = os.path.join(project_root, "drug_data", drug_name.replace(" ", "_"))

            if not os.path.exists(drug_folder):
                os.makedirs(drug_folder)

            file_path = os.path.join(drug_folder, "mechanism.txt")
            save_readable_info(drug_mechanism, file_path)

            print(f"Информация о механизме действия сохранена в {file_path}")
        else:
            print("Информация о механизме действия не найдена.")
            sys.exit(1)
    else:
        print("Препарат не найден.")
        sys.exit(1)

if __name__ == "__main__":
    main()
