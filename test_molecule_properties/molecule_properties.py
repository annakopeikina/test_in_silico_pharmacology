import json
import os
import sys
from chembl_webresource_client.new_client import new_client
import io

# Устанавливаем кодировку UTF-8 для вывода
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def fetch_drug_data(drug_name):
    molecule = new_client.molecule
    try:
        result = molecule.filter(preferred_name=drug_name).only(['molecule_chembl_id', 'molecule_properties'])
        if result:
            return result[0]
        return {'error': 'Препарат не найден'}
    except Exception as e:
        return {'error': str(e)}

def save_data(drug_name, data, folder_path):
    try:
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, 'molecule_properties.json')
        
        with open(file_path, 'w', encoding='utf-8') as file:
            json.dump(data, file, indent=4, ensure_ascii=False)
        
        print(f"Данные о препарате '{drug_name}' сохранены в '{file_path}'.")
    except Exception as e:
        print(f"Ошибка при сохранении данных: {e}")

def main():
    if len(sys.argv) != 3:
        print("Использование: python molecule_properties.py <drug_name> <folder_path>")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    folder_path = sys.argv[2]

    data = fetch_drug_data(drug_name)
    print(f"Полученные данные: {data}")  # Отладочное сообщение
    save_data(drug_name, data, folder_path)

if __name__ == "__main__":
    main()
