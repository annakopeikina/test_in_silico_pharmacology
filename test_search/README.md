# Описание

Для поиска лекарственных препаратов по торговому названию, латинскому наименованию или брутто-формуле могут быть использованы специализированные API. Существуют несколько ресурсов, которые предоставляют необходимые данные:

DrugBank (<https://www.drugbank.ca>): Предоставляет доступ к API для поиска информации о лекарственных препаратах по их названиям, химической структуре и другим характеристикам. API ключи являются коммерческими, ими целесообразно воспользоваться, если проект финансируется компанией или государствомю.

PubChem (<https://pubchem.ncbi.nlm.nih.gov/>): API этого ресурса позволяет искать данные о химических соединениях, веществах и их свойствах по названию, формуле или идентификатору PubChem. Пример запроса:
<https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{название}/property/CanonicalSMILES,InChIKey,Formula/JSON>. PubChem предоставляет доступ к данным через PUG REST API (Power User Gateway), который не требует генерации индивидуального API-ключа. Этот API бесплатен и предоставляет программный доступ к различным типам данных о химических соединениях, веществах и их свойствах.

ChEMBL (<https://www.ebi.ac.uk/chembl/>): Этот API ориентирован на биологически активные соединения и позволяет искать препараты по названию, молекуле и другим параметрам. Пример запроса:
<https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}>.

RxNorm (<https://rxnav.nlm.nih.gov/>): API ресурса предназначено для получения информации о лекарствах и их торговых названиях. Пример запроса:
<https://rxnav.nlm.nih.gov/REST/rxcui?name={название_лекарства}>.

Каждый из этих ресурсов содержит подробную документацию с описанием параметров запросов и форматов ответов. Выбор конкретного API зависит от типа данных, которые требуется получить для решения задачи.

Для наиболее полной информации о фармакологических данных выберем DrugBank <https://www.drugbank.com/> — это один из лучших вариантов. Он предоставляет обширную базу данных, которая включает как структурные, так и фармакологические данные о лекарственных веществах, их механизмах действия, метаболизме, взаимодействиях, побочных эффектах и многом другом.

Почему стоит выбрать DrugBank:
Обширные фармакологические данные: DrugBank содержит информацию о лекарствах, их механизмах действия, побочных эффектах, взаимодействиях с другими препаратами и их метаболических путях.
API-доступ: API-ключ, который позволит автоматизировать запросы к базе данных.
Информация о клинических исследованиях: DrugBank также предоставляет данные о клиническом использовании лекарств, что полезно для in silico исследований.
Поддержка различных форматов данных: поиск лекарства по различным параметрам — торговому названию, химической структуре, терапевтическим свойствам и другим критериям.
Download data <https://go.drugbank.com/releases/latest#open-data>
Поскольку DrugBank предоставляет только коммерческие API (годовую подписку могут позволить только крупные корпорации, фармкомпании или государство), этот вариант можно рассматривать только в случае такой коммерческой заинтересованности.

Для реализации проекта был выбран бесплатный ресурс <https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access>

## Drug Information Search Tool

Этот скрипт предназначен для поиска информации о препаратах с использованием PubChem API. При вводе названия препарата, скрипт отправляет запросы к API и собирает данные, такие как латинское название, брутто формула и химическое название (InChIKey). Информация сохраняется в специально созданных папках в директории проекта.

## Функционал

- Поиск по торговому названию препарата.
- Поиск по латинскому названию, если по торговому названию не найдено.
- Поиск по брутто формуле, если по торговому и латинскому названиям не найдено.
- Удаление старых данных и создание новой папки для хранения информации о каждом новом запросе.

## Использование

1. **Установка зависимостей**:
   - Убедитесь, что у вас установлен Python версии 3.10.11 или выше.
   - Создайте виртуальное окружение и активируйте его:

     ```bash
     python -m venv .venv
     source .venv/bin/activate  # Для Windows используйте: .venv\Scripts\activate
     ```

   - Установите необходимые библиотеки:

     ```bash
     pip install requests
     ```

2. **Запуск скрипта**:
   - Поместите файл `search.py` в корень проекта.
   - Запустите скрипт с помощью команды:

     ```bash
     python search.py
     ```

   - Введите название препарата латинскими буквами, когда будет предложено.

## Пример использования

1. Запустите скрипт:

   ```bash
   python search.py

### 2. Установка зависимостей

После активации виртуального окружения установите необходимые зависимости:
pip install -r requirements.txt

Пример:
Введите название препарата, например: Amiodarone.

Скрипт создаст папку drug_data (если она уже существует, она будет удалена) и внутри нее папку Amiodarone.

В папке Amiodarone будет создан файл info.txt с найденной информацией о препарате.

Примечания
При каждом новом запросе папка drug_data в корне проекта будет удаляться и заново создаваться, чтобы хранить только данные последнего запроса.
Папка для конкретного препарата будет создана внутри drug_data с именем, соответствующим названию препарата (пробелы заменены на подчеркивания).
Собранные данные не сохраняются, так как предполагается работа с большими объемами данных.
Лицензия
Пока эта разработка находится в стадии тестирования, и её использование запрещено. Лицензия будет предоставлена позднее.
