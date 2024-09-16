# In Silico Pharmacology Project

## Описание

Этот проект предназначен для автоматизации поиска данных о препаратах, анализа их свойств, фармакодинамики и визуализации. Включает несколько ключевых компонентов:

1. **main.py**: Основной скрипт для запуска проекта. Выполняет анализ молекулярных свойств, фармакодинамических характеристик и визуализацию молекул. Этот скрипт объединяет функциональность всех модулей и управляет процессом обработки данных.

2. **search.py**: Модуль поиска информации о препаратах. Использует PubChem API для поиска препаратов по их торговому названию и общему названию. Реализует функции для поиска и получения данных о препарате.

3. **molecule_properties.py**: Модуль получения молекулярных свойств из внешних API. Содержит функции для получения данных о препарате и сохранения этих данных в указанную папку. Обрабатывает информацию о молекуле и ее свойствах.

4. **pharmacodynamics.py**: Модуль для анализа механизма действия препаратов. Включает функции для получения ID ChEMBL для препарата, извлечения информации о механизме действия по этому ID и сохранения информации в читаемом формате.

5. **molecule_visuals.py**: Модуль для визуализации структуры молекул. Реализует функции для получения SMILES строки для препарата и конвертации этой строки в изображение молекулы, которое сохраняется в указанном формате.

## Структура проекта

**in_silico_pharmacology**/
│
├── main.py                                              # Основной скрипт для запуска проекта
├── README.md                                            # Документация по проекту
├── drug_data/                                           # Папка для хранения данных о препаратах
│   └── Atenolol/                                        # Пример папки для препарата (заменить на нужные препараты)
│       ├── info.txt                                     # Общая информация о препарате
│       ├── mechanism.json                               # Данные о механизме действия (формат JSON)
│       ├── mechanism.txt                                # Механизм действия в текстовом формате
│       ├── molecule_properties.json                     # Молекулярные свойства препарата (формат JSON)
│       └── molecule.png                                 # Изображение молекулы препарата
│
├── search/                                              # Модуль поиска информации о препаратах
│   └── search.py
│       ├── def search_by_trade_name(trade_name)         # Поиск препарата по торговому названию
│       ├── def search_drug(drug_name)                   # Поиск препарата по его названию
│
├── molecule_properties/                                 # Модуль получения молекулярных свойств из внешних API
│   └── molecule_properties.py
│       ├── def fetch_drug_data(drug_name)               # Получение данных о препарате
│       ├── def save_data(drug_name, data, folder_path)  # Сохранение данных о препарате в папку
│
├── pharmacodynamics/                                    # Модуль для механизма действия препарата
│   └── pharmacodynamics.py
│       ├── def get_chembl_id(drug_name)                 # Получение ID ChEMBL для препарата
│       ├── def get_drug_mechanism(drug_id)              # Получение механизма действия по ID ChEMBL
│       ├── def save_readable_info(drug_info, file_path) # Сохранение информации о механизме действия
│
├── molecule_visuals/                                    # Модуль для визуализации структуры молекул
    └── molecule_visuals.py
        ├── def fetch_smiles(drug_name)                  # Получение SMILES строки для препарата
        ├── def visualize_molecule(smiles, output_path)  # Конвертация SMILES в изображение молекулы

drug_data/: Папка для хранения данных о препаратах. Включает поддиректории для каждого препарата с файлами информации о препарате, механизме действия, молекулярных свойствах и изображении молекулы.

search/: Модуль для поиска информации о препаратах.

molecule_properties/: Модуль для получения молекулярных свойств.

pharmacodynamics/: Модуль для анализа механизма действия.

molecule_visuals/: Модуль для визуализации молекул.

Установка и запуск
Установите необходимые зависимости, выполнив команду:

bash
Copy code
pip install -r requirements.txt
Запустите основной скрипт:

bash
Copy code
python main.py
Для поиска данных о препарате используйте функции из search.py.

Для получения и сохранения молекулярных свойств используйте функции из molecule_properties.py.

Для анализа механизма действия используйте функции из pharmacodynamics.py.

Для визуализации молекул используйте функции из molecule_visuals.py.

Пример использования
В main.py реализован пример использования всех функций, объединяющий их в один рабочий процесс.

### Требования

Для работы проекта необходим Python версии 3.10 или 3.11.

### Установка зависимостей

Проект требует установки нескольких библиотек, включая requests для HTTP-запросов к внешним API, rdkit для работы с молекулярными данными и matplotlib для визуализации молекул.
Создайте файл `requirements.txt`, куда добавьте необходимые библиотеки:

```txt
requests
requests>=2.18.4
requests-cache~=0.7.0
easydict
urllib3
lxml
chembl_webresource_client
rdkit


# pip install -r requirements.txt
# source myenv/Scripts/activate

###создание .gitignore 
Виртуальное окружение
myenv/
venv/
env/
ENV/

PyCharm
.idea/

Модули Python
__pycache__/
*.pyc
*.pyo
*.pyd
*.pdb

Jupyter Notebook
.ipynb_checkpoints

Директории и файлы для систем сборки
build/
dist/
*.egg-info/
*.egg
*.log

Конфиденциальные файлы
*.env
*.secret

Макеты и временные файлы
*.swp
*.swo
*.tmp
*.bak

Файлы кэша и состояния
*.cache
*.coverage
*.coverage.*
nosetests.xml
coverage.xml

 Другие
*.DS_Store
*.directory
*.vscode/


Добавление пустого файла __init__.py в корневую папку проекта и в каждую подпапку, содержащую модули, пустой файл __init__.py сделает папки распознаваемыми как пакеты Python.

Активация виртуального окружения:
bash -
myenv\Scripts\activate
pshell-
source /c/Users/annav/OneDrive/Desktop/test_in_silico_pharmacology/myenv/Scripts/activate

Регистрация виртуального окружения для Jupyter Notebook
pip install ipykernel
python -m ipykernel install --user --name=myenv --display-name "MyEnv"

## Запуск программы: python main.py <Drug_name>

    Например: python main.py Atenolol

## для .ipynb Jupyter Notebook

%pip install pandas
%pip install matplotlib
%pip install seaborn
%pip install requests
%pip install chardet
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import chardet
from rdkit import Chem
from rdkit.Chem import Draw
import os
import chardet
