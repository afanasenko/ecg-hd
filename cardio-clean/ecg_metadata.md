**Описание метаданных, генерируемых при анализе кардиосигнала**

Функция *blobapi_detect_qrs*: обнаружение QRS-комплексов и сегментация отдельных волн.

qrs_metadata = blobapi_detect_qrs(inbuf, min_qrs_ms=20, postprocessing=True)

    Предварительные условия:
    blobapi_detect_qrs желательно вызывать после функций подавления
    дрейфа изолинии blobapi_fix_baseline и подавления сетевой помехи
    blobapi_mains_correction (порядок вызова не имеет значения)

    :param inbuf: входной буфер (остается неизменным)
    :param min_qrs_ms: минимальная длительность QRS-комплекса
    :param delineate: выполнить сегментацию найденных комплексов
    :return: qrs_metadata (список найденных комплексов)


Возвращает метаданные в виде списка словарей.

```
>>> qrs_metadata
[qrs0, qrs1, qrs2, ...]
```
Каждый словарь в списке содержит информацию по одному QRS-комплексу в виде пар "ключ-значение"

Примечания:

Параметры, которые невозможно определить, записываются как None.
Параметры, которые не могут быть None после вызова соответствуюей функции,
отмечены **жирным шрифтом**.
Параметры, имеюие тип данных array, индивидуально рассчитываются в каждом отведдении.


| Ключ | Назначение | Тип данных | Размерность | Какая процедура рассчитывает |
| ---- |:---------- | :--------- | :---------- | ---------------------------: |
| **qrs_start** | Начало QRS | float | [с] от начала записи | qrs_detection |
| **qrs_end** | Конец QRS | float | [с] от начала записи | qrs_detection |
| **qrs_center** | Середина QRS | float | [с] от начала записи | qrs_detection |
| qrs_class_id | № класса QRS | int | - | incremental_classifier |
| **artifact** | артефакт | bool | - | incremental_classifier |
| qrsType | форма комплекса | string | строковый код | find_points |
| p_start | начало P-зубца | int array | № отсчета | find_points |
| p_end | конец P-зубца | int array | № отсчета | find_points |
| p_pos | вершина P-зубца | int array | № отсчета | find_points |
| p_height | высота P-зубца | float array | как у сигнала | metadata_postprocessing |
| q_pos | вершина Q-зубца | int array | № отсчета | find_points |
| q_height | высота Q-зубца | float array | как у сигнала | metadata_postprocessing |
| r_pos | вершина R-зубца | int array | № отсчета | find_points |
| r_height | высота R-зубца | float array | как у сигнала | metadata_postprocessing |
| s_pos | вершина S-зубца | int array | № отсчета | find_points |
| s_height | высота S-зубца | float array | как у сигнала | metadata_postprocessing |
| t_start | начало T-зубца | int array | № отсчета | find_points |
| t_end | конец T-зубца | int array | № отсчета | find_points |
| t_pos | вершина T-зубца | int array | № отсчета | find_points |
| t_height | высота T-зубца | float array | как у сигнала | metadata_postprocessing |
| параметры ритма |
| RR | RR-интервал | float | [мс] | metadata_postprocessing |
| heartrate | ЧСС = 60000/RR| float | уд/мин | metadata_postprocessing |
| isolevel | уровень изолинии | float array | как у сигнала | find_points |
| ST-сегмент |
| st_start | Начало ST (J) | int array | № отсчета | metadata_postprocessing |
| st_plus | Точка J+0.08 | int array | № отсчета | metadata_postprocessing |
| st_end | Конец ST = начало T | int array | № отсчета | metadata_postprocessing |
| st_start_level | Смещение J от изолинии | float array | как у сигнала | metadata_postprocessing |
| st_plus_level | Смещение J+ от изолинии | float array | как у сигнала | metadata_postprocessing |
| st_end_level | Смещение конца ST от изолинии | float array | как у сигнала | metadata_postprocessing |
| st_offset | Средниее смещение ST от изолинии | float array | как у сигнала | metadata_postprocessing |
| st_duration | Длительность ST | float array | [мс] | metadata_postprocessing |
| st_slope | Наклон ST | float array | отн. ед. | metadata_postprocessing |
