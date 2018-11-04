## Метаданные, генерируемые при анализе кардиосигнала

Функция *blobapi_detect_qrs* выполняет обнаружение QRS-комплексов и сегментация отдельных волн.

qrs_metadata = blobapi_detect_qrs(inbuf, min_qrs_ms=20, postprocessing=True)

    Предварительные условия:
    blobapi_detect_qrs желательно вызывать после функций подавления
    дрейфа изолинии blobapi_fix_baseline и подавления сетевой помехи
    blobapi_mains_correction (порядок вызова не имеет значения)

    :param inbuf: входной буфер (остается неизменным)
    :param min_qrs_ms: минимальная длительность QRS-комплекса
    :param postprocessing: расчет вторичных параметров (ритм, ST и др.)
    :return: qrs_metadata (список найденных комплексов)


Возвращает метаданные в виде списка словарей.

```
>>> qrs_metadata
[qrs0, qrs1, qrs2, ...]
```
Каждый словарь в списке содержит информацию по одному QRS-комплексу в виде пар "ключ-значение"

Внутри blobapi_detect_qrs выполняется последовательный вызов 3 функций:
*qrs_detection* - первичное выделение QRS-комплексов
*find_points* - поиск характерных точек
*metadata_postprocessing* - расчет вторичных параметров

Для повторного расчета вторичных параметров (например,
после ручной коррекции характерных точек), необходимо вызвать функцию
*blobapi_postprocessing_qrs*

Примечания:

Параметры, которые невозможно определить, записываются как None.
Параметры, которые никогда не могут быть None после вызова соответствуюей функции,
отмечены **жирным шрифтом**.
Тип данных "array" указывает на то, что параметр индивидуально рассчитывается в каждом отведении.


| Ключ | Расшифровка | Тип данных | Размерность или диапазон | Какая процедура рассчитывает |
| ---- |:---------- | :--------- | :---------- | ---------------------------: |
| **qrs_start** | Начало QRS | float | [с] от начала записи | qrs_detection |
| **qrs_end** | Конец QRS | float | [с] от начала записи | qrs_detection |
| **qrs_center** | Середина QRS | float | [с] от начала записи | qrs_detection |
| **complex_type** | тип комплекса | char | 'N'\|'S'\|'V'\|'U' | metadata_postprocessing |
| qrs_class_id | № класса QRS | int | - | incremental_classifier |
| **flags** | флаги "артефакт" или "экстрасистола" | string | ''\|'A'\|'E' | metadata_postprocessing, incremental_classifier |
| p_start | начало P-зубца | int array | № отсчета | find_points |
| p_end | конец P-зубца | int array | № отсчета | find_points |
| p_pos | вершина P-зубца | int array | № отсчета | find_points |
| p_height | высота P-зубца над/под изолинией| float array | мВ | metadata_postprocessing |
| q_pos | вершина Q-зубца | int array | № отсчета | find_points |
| q_height | высота Q-зубца над/под изолинией| float array | мВ | metadata_postprocessing |
| r_start | начало R-зубца | int array | № отсчета | find_points |
| r_end | конец R-зубца | int array | № отсчета | find_points |
| r_pos | вершина R-зубца | int array | № отсчета | find_points |
| r_height | высота R-зубца над/под изолинией | float array | мВ | metadata_postprocessing |
| s_pos | вершина S-зубца | int array | № отсчета | find_points |
| s_height | высота S-зубца над/под изолинией | float array | мВ | metadata_postprocessing |
| t_start | начало T-зубца | int array | № отсчета | find_points |
| t_end | конец T-зубца | int array | № отсчета | find_points |
| t_pos | вершина T-зубца | int array | № отсчета | find_points |
| t_height | высота T-зубца над/под изолинией | float array | мВ | metadata_postprocessing |
| **параметры ритма** |
| RR | RR-интервал | float | [с] | metadata_postprocessing |
| heartrate | ЧСС = 60000/RR| float | уд/мин | metadata_postprocessing |
| isolevel | уровень изолинии | float array | как у сигнала | find_points |
| **ST-сегмент** |
| st_start | Начало ST (J) | int array | № отсчета | metadata_postprocessing |
| st_plus | Точка J+0.08 | int array | № отсчета | metadata_postprocessing |
| st_end | Конец ST = начало T | int array | № отсчета | metadata_postprocessing |
| st_start_level | Смещение J от изолинии | float array | мВ | metadata_postprocessing |
| st_plus_level | Смещение J+ от изолинии | float array | мВ | metadata_postprocessing |
| st_end_level | Смещение конца ST от изолинии | float array | мВ | metadata_postprocessing |
| st_offset | Среднее смещение ST от изолинии | float array | мВ | metadata_postprocessing |
| st_duration | Длительность ST | float array | [мс] | metadata_postprocessing |
| st_slope | Наклон ST | float array | - | metadata_postprocessing |
| **QT-интервал** |
| qt_duration | Длительность QT интервала| float array | [с] | metadata_postprocessing |
| qtc_duration | Длительность корригированного QT интервала | float array | [с] | metadata_postprocessing |

Особенности оценки параметров
*qt_duration*
характеризует атрио-вентрикулярную проводимость, то есть
проведение электрического импульса через соединение между предсердиями
и желудочками (через АВ-узел). В случае, если в комплексе отсутствует
выраженный зубец Q, началом интервала QT считается начало R-зубца.