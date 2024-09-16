# Модуль `molecule_properties`

Этот модуль отвечает за получение данных о молекулярных свойствах лекарственных препаратов и сохранение их в формате JSON.

## Функции

- `fetch_drug_data(drug_name, folder_path)`: Получает данные о молекулярных свойствах препарата и сохраняет их в JSON файл.

## Использование

1. Импортируйте модуль:
python
    fetch_drug_data("Атенолол", "./drug_data/Atenolol")
    from molecule_properties.molecule_properties import fetch_drug_data
    fetch_drug_data("Атенолол", "./drug_data/Atenolol")
Это сохранит данные о молекулярных свойствах в файл molecule_properties.json в указанной папке.(Atenolol используется как пример)
