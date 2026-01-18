from sklearn.model_selection import StratifiedKFold, train_test_split
import argparse
import json
import numpy as np
import os
import pandas as pd
import time
import warnings

from training_func.utils import *
import training_func.environment as env


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-name', default='RDDL', type=case_name,
                        help='name of this deep learning project')
    parser.add_argument('-f', default=5, type=int, help=argparse.SUPPRESS, dest='k_fold')
    parser.add_argument('-e', default=0.1, type=data_size,
                        help='the ensemble set percentage of the whole dataset', dest='ensemble_size')
    parser.add_argument('-t', default=0.1, type=data_size,
                        help='the test set percentage of the whole dataset', dest='test_size')
    parser.add_argument('-w', default='sample_weights.csv', help=argparse.SUPPRESS, dest='sw_file')
    parser.add_argument('-r', default=10, type=int,
                        help='same random seed lead to same results', dest='random_seed')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='do not display the program messages on screen')
    a = parser.parse_args()

    env.environment_setup(a.name, quiet_mode=a.quiet)

    if a.ensemble_size + a.test_size > 1.0:
        parser.error('argument -e, -t: the sum of the splitting percentages for the ensemble set and the test set should be less than 1.0.')

    data = pd.DataFrame(columns=['sample', 'label', 'weight'])

    for label, source_dir in [(1, env.WORKSPACE['USER_pos']), (0, env.WORKSPACE['USER_neg'])]:

        sw_file_path = source_dir / a.sw_file
        sample_list = list_files(source_dir, cond=lambda x: x != a.sw_file)
        index = 'positive' if label == 1 else 'negative'

        if sw_file_path.exists():
            weight_map = {x: w for x, w in pd.read_csv(sw_file_path).values}
        else:
            weight_map = {sample: 1.0 for sample in sample_list}
            warnings.warn('%s is not provided, the weights of all %s data will be set to 1.0.' % (sw_file_path, index))
            time.sleep(0.75)

        missing_weights = list()

        for sample in sample_list:

            if sample in weight_map:
                weight = weight_map[sample]
            else:
                weight = 1.0
                missing_weights.append(sample)

            data.loc[len(data)] = [str(source_dir / sample), label, weight]

        if missing_weights:
            if len(missing_weights) > 3:
                missing_weights = missing_weights[:3] + ['%d more' % (len(missing_weights) - 3)]
            msg = ' and '.join(', '.join(missing_weights).rsplit(', ', 1))
            warnings.warn('The %s data %s are not specified in %s, the weights will be set to 1.0.' % (index, msg, sw_file_path))
            time.sleep(0.75)

    other_indices, test_indices = train_test_split(np.arange(len(data)), test_size=a.test_size,
                                                   random_state=a.random_seed, stratify=data['label'])
    training_indices, ensemble_indices = train_test_split(other_indices, test_size=a.ensemble_size / (1 - a.test_size),
                                                          random_state=a.random_seed, stratify=data.loc[other_indices]['label'])

    splitting_info_path = env.WORKSPACE['RDDL_splitting_info']
    cv_info = dict()
    sfolder = StratifiedKFold(n_splits=a.k_fold, random_state=a.random_seed, shuffle=True)

    for file_name in list_files(splitting_info_path, cond=lambda x: x.endswith('.csv')):
        (splitting_info_path / file_name).unlink()

    for cv, (_, val) in enumerate(sfolder.split(training_indices, data.loc[training_indices]['label'])):
        val_indices = training_indices[val]
        file_name = splitting_info_path / ('training_%d.csv' % (cv + 1))
        data.loc[val_indices].to_csv(file_name, index=False)
        cv_info[file_name] = np.bincount(data.loc[val_indices]['label'], minlength=2)
   
    file_name = splitting_info_path / 'ensemble.csv'
    data.loc[ensemble_indices].to_csv(file_name, index=False)
    cv_info[file_name] = np.bincount(data.loc[ensemble_indices]['label'], minlength=2)

    file_name = splitting_info_path / 'test.csv'
    data.loc[test_indices].to_csv(file_name, index=False)
    cv_info[file_name] = np.bincount(data.loc[test_indices]['label'], minlength=2)
    
    with open(splitting_info_path / ('data_splitting_%s.log' % get_nowtime()), 'w') as f:
        f.write('Case name: %s\n' % a.name)
        f.write('K-fold: %d\n' % a.k_fold)
        f.write('Size of ensemble set: %f\n' % a.ensemble_size)
        f.write('Size of test set: %f\n' % a.test_size)
        f.write('Random seed: %d' % a.random_seed)
    
    for file_name, (neg_count, pos_count) in cv_info.items():
        total_count = neg_count + pos_count
        print(file_name)
        print('=' * 50)
        print('Positive  %5d  (%4.1f%% of this split)' % (pos_count, pos_count / total_count * 100))
        print('Negative  %5d  (%4.1f%% of this split)' % (neg_count, neg_count / total_count * 100))
        print('-' * 50)
        print('Total     %5d  (%4.1f%% of all data)\n' % (total_count, total_count / len(data) * 100))

    print('Total data: %d' % len(data))
    print('Positive data: %d' % (data['label'] == 1).sum())
    print('Negative data: %d' % (data['label'] == 0).sum())
