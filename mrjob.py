#!/usr/bin/env python
# coding: utf-8

import sys
import os

"""
    Mapper, if we need to change input format
    Output: class_mark | features
"""
def mapper_main(delim):

    for line in sys.stdin:

        p = [x.strip() for x in line.split(delim)]
        features = p[:-1]
        class_mark = p[-1]

        sys.stdout.write("\t".join([class_mark, features]))
"""
    Reducers to prepare data for means counting
    Output: class_mark | n_obj | feature_1_sum | ... | feature_n_sum
"""
def means_reducer():
    mem = []

    for line in sys.stdin:

        p = [x.strip() for x in line.split("\t")]

        class_mark = p[0]
        features = p[1:]

        if mem:
            if class_mark == mem[0]:
                # objects count
                mem[1] += 1
                # sum of respective features
                for i, f in enumerate(features):
                    mem[i+2] += f
            else:
                sys.stdout.write("\t".join([str(x) for x in mem]))
                mem[0] = class_mark
                for i, f in enumerate(features):
                    mem[i+2] = f
                mem[1] = 1
        else:
            mem[0] = class_mark
            for i, f in enumerate(features):
                mem[i+1] = f
            mem[1] = 1

    sys.stdout.write("\t".join([str(x) for x in mem]))

""" 
    Merge reducers output to compute mean vectors, overall mean 
    Output: sample_data = {"class_mark": [n_obj, feature_1_sum, ..., feature_n_sum]}
"""
def read_mean_files():
    sample_data = {}

    for line in sys.stdin:
        p = [x.strip() for x in line.split('\t')]
        class_mark = p[0]
        if class_mark not in sample_data:
            sample_data[class_mark] = np.array(p[1:])
        else:
            sample_data[class_mark] += np.array(p[1:])

    return sample_data

""" 
    Reducer to compute scatter within matrix 
    Output: flattened scatter within matrix
    !!! Needs to run on a single reducer.
"""
# def scatter_within_one_reducer(mean_v, n_features):
#     mem = ""
#     s_i = None
#     s_w = np.zeros((n_features, n_features))

#     for line in sys.stdin:

#         p = [x.strip() for x in line.split("\t")]

#         class_mark = p[0]
#         features = np.asarray(map(float, p[1:]), dtype=np.float64)

#         if mem:
#             if class_mark == mem:
#                 s_i += (features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).dot((features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).T) 
#             else:
#                 s_w += s_i
#                 mem = class_mark
#                 s_i = np.zeros((n_features, n_features))
#         else:
#             mem = class_mark
#             s_i = np.zeros((n_features, n_features))

#     s_w += s_i
#     print "\t".join(map(str, np.ravel(s_w)))
"""
    Переделала функцию для нескольких reducer-ов
    Output: flattened scatter matrix for each class
"""
def scatter_within_reducer(mean_v, n_features):
    mem = ""
    s_i = None

    for line in sys.stdin:

        p = [x.strip() for x in line.split("\t")]

        class_mark = p[0]
        features = np.asarray(map(float, p[1:]), dtype=np.float64)

        if mem:
            if class_mark == mem:
                s_i += (features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).dot((features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).T) 
            else:
                sys.stdout.write("\t".join([mem, map(str, np.ravel(s_i))]))
                mem = class_mark
                s_i = np.zeros((n_features, n_features))
                s_i += (features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).dot((features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).T) 
        else:
            mem = class_mark
            s_i = np.zeros((n_features, n_features))
            s_i += (features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).dot((features.reshape(n_features, 1) - mean_v[class_mark].reshape(n_features, 1)).T) 

    sys.stdout.write("\t".join([mem, map(str, np.ravel(s_i))]))

"""
    функция считывает выходные данные reducer-ов, 
    вычисляет сумму и преобразует вектор обратно в матрицу для дальнейших расчетов
"""
def read_scatter_files(n_features):
    s_w = np.zeros((n_features, n_features))
    s_w_flatten = np.zeros((n_features * n_features))

    for line in sys.stdin:
        p = [x.strip() for x in line.split('\t')]
        s_w_flatten += np.array(map(float, p[1:]))
    s_w_flatten = np.ravel(s_w_flatten)
    s_w = np.array([s_w_flatten[i:i + n_features] for i in range(0, len(s_w_flatten), n_features)])

    return s_w

"""  
    Вектор проекций.
    Сначала рассчитываем собственные значения и векторы для матрицы s_w^(-1)*s_b
    Найденные значения и векторы попарно помещаем в список, 
    при этом собственные значения берем по модулю, а затем сортируем по ним список.
    Теперь собственные значения расположены в порядке уменьшения.
    Вытаскиваем вектор, соответствующий первому(самому большому) собственному значению.
    Это и есть искомый вектор для классификации.
"""
def projection_vector(s_w, s_b):
    w, v = np.linalg.eig(np.linalg.inv(s_w).dot(s_b))
    eig_list = [(np.abs(w[i]), v[:,i]) for i in range(len(w))]
    eig_list = sorted(eig_list,key=lambda tup:tup[0], reverse=True)
    p_v = eig_list[0][1].reshape(n,1)
    return p_v

def build_options():
    parser = argparse.ArgumentParser()
    nowtime = time.time()

    parser.add_argument(
        '-d', '--delimiter',
        default='\t',
        type=str,
        help='Delimiter of input file'
    )
    parser.add_option(
        '-c', '--command',
        type='string',
        dest='command',
        default='mapper',
        help='Command to launch mapper, reducer for means or for scatter matrix'
    )
    parser.add_argument(
        '-i', '--input-dir',
        default='',
        type=str,
        help='Root input folder'
    )
    parser.add_argument(
        '-o', '--output-dir',
        default='',
        type=str,
        help='Root output folder'
    )

    options, other_params = parser.parse_known_args()
    return options

def main():
    options = build_options()
    delim = options.delimiter
    cmd = options.command

    if cmd == "mapper":
        mapper_main(delim)
    elif cmd == "reducer_m":
        means_reducer()
    elif cmd == "reducer_s":
        scatter_within_reducer(mean_v, n_features)

    sample_data = read_mean_files(sys.stdin)   
    n_features = len(sample_data[sample_data.keys()[0]][1:])
    mean_v = {class_mark: [] for class_mark in sample_data}
    all_obj = 0
    overallmean = np.zeros((n_features, 1)) 
    s_b = np.zeros((n_features, n_features))

    # compute mean vectors for each class and overall mean
    for class_mark, data in sample_data.items():
        n_obj = data[0]
        mean_v[class_mark] = data[1:] / n_obj
        all_obj += n_obj
        overallmean += mean_v[class_mark].reshape(n_features,1)
    overallmean = overallmean / all_obj
    
    # compute between class scatter matrix
    for class_mark, data in sample_data.items():
        s_b += data[0] * (mean_v[class_mark].reshape(n_features, 1) - overallmean).dot((mean_v[class_mark].reshape(n_features,1) - overallmean).T)

    p_v = projection_vector(s_w, s_b)
    
if __name__ == "__main__":
    main()
