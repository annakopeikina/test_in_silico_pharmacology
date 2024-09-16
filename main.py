import subprocess
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("Использование: python main.py <drug_name>")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    
    # Указание корневой папки проекта
    project_root = r"C:\Users\annav\OneDrive\Desktop\test_in_silico_pharmacology"
    drug_folder = os.path.join(project_root, "drug_data", drug_name.replace(" ", "_"))

    # Создание папки, если она не существует
    if not os.path.exists(drug_folder):
        os.makedirs(drug_folder)

    try:
        # Путь к molecule_properties.py
        molecule_properties_path = 'C:/Users/annav/OneDrive/Desktop/test_in_silico_pharmacology/test_molecule_properties/molecule_properties.py'

        # Путь к pharmacodynamics.py
        pharmacodynamics_path = 'C:/Users/annav/OneDrive/Desktop/test_in_silico_pharmacology/test_pharmacodynamics/pharmacodynamics.py'

        # Путь к molecule_visual.py
        molecule_visual_path = 'C:/Users/annav/OneDrive/Desktop/test_in_silico_pharmacology/test_molecule_visual/molecule_visual.py'

        # Выполнение molecule_properties.py
        result_molecule = subprocess.run(
            ['python', molecule_properties_path, drug_name, drug_folder],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )

        if result_molecule.returncode != 0:
            print(f"Ошибка выполнения molecule_properties.py: {result_molecule.stderr}")
            sys.exit(1)
        
        print(result_molecule.stdout)

        # Выполнение pharmacodynamics.py
        result_pharmacodynamics = subprocess.run(
            ['python', pharmacodynamics_path, drug_name],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )

        if result_pharmacodynamics.returncode != 0:
            print(f"Ошибка выполнения pharmacodynamics.py: {result_pharmacodynamics.stderr}")
            sys.exit(1)
        
        print(result_pharmacodynamics.stdout)

        # Выполнение molecule_visual.py
        result_visual = subprocess.run(
            ['python', molecule_visual_path, drug_name, drug_folder],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )

        if result_visual.returncode != 0:
            print(f"Ошибка выполнения molecule_visual.py: {result_visual.stderr}")
            sys.exit(1)
        
        print(result_visual.stdout)

    except Exception as e:
        print(f"Ошибка: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
