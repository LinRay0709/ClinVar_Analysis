from tensorflow.keras.preprocessing.image import ImageDataGenerator, img_to_array, load_img
import csv
import imghdr
import math
import numpy as np
import os
import pandas as pd
import pickle
import sys
import tensorflow as tf
import warnings

from training_func.utils import balancing_methods, list_files
import training_func.environment as env


def cal_z_stats(samples, input_num, nor_mode=0):
    
    z_mean = [0.0] * input_num

    if nor_mode == 0:
        z_std = [1.0] * input_num

    elif nor_mode == 1:
        z_std = [255.0] * input_num

    elif nor_mode == 2:

        ss_mean = [0.0] * input_num

        for sample in samples:
            x = make_data_batch([sample], input_num)
            for i in range(input_num):
                z_mean[i] += np.mean(x[i]) / len(samples)
                ss_mean[i] += np.mean(np.power(x[i], 2)) / len(samples)
        
        z_std = np.sqrt(ss_mean - np.power(z_mean, 2)).tolist()

    if input_num == 1:
        return z_mean[0], z_std[0]
    else:
        return z_mean, z_std


# Return generator
def load_data(split, fold_num, input_num, batch_size=1, method='Baseline', nor_mode=0, z_mean=None, z_std=None):

    data_dir = env.WORKSPACE['RDDL_splitting_info']

    if split == 'training':
        data = pd.DataFrame()
        for i in range(1, env.K_FOLD + 1):
            if i == fold_num:
                continue
            data = pd.concat([data, pd.read_csv(data_dir / ('training_%d.csv' % i))],
                             axis=0, ignore_index=True)
        data = resample(data, split, method)
        z_mean, z_std = cal_z_stats(data['sample'], input_num, nor_mode=nor_mode)
        alpha = 1 - np.mean(data['label'])
    
    elif split == 'validation':
        data = pd.read_csv(data_dir / ('training_%d.csv' % fold_num))
        data = resample(data, split, method, shuffle=False)
    
    elif split == 'ensemble' or split == 'test':
        data = pd.read_csv(data_dir / ('%s.csv' % split))

    else:
        raise ValueError('the name of data split can only be '
                         '\"training\", \"validation\", \"ensemble\" or \"test\"')

    generator = create_generator(data, input_num, batch_size, z_mean, z_std,
                                 balanced_batch=(split == 'training' and method == 'BB'),
                                 specific_weight=(split == 'training' and method == 'SW'))
    steps = int(math.ceil(len(data) / batch_size))
    labels = data['label'].tolist()

    if split == 'training':
        return generator, steps, labels, z_mean, z_std, alpha
    else:
        return generator, steps, labels


def create_generator(data, input_num, batch_size, z_mean, z_std, balanced_batch=False, specific_weight=False):

    data = data.reset_index(drop=True)

    if balanced_batch:
        sampler = np.reshape([data.index[data['label'] == _] for _ in (0, 1)], -1, order='F')
    else:
        sampler = np.arange(len(data))

    while True:

        for i in range(int(math.ceil(len(data) / batch_size))):

            batch_data = data.loc[sampler[i * batch_size:(i + 1) * batch_size]]
            x = make_data_batch(batch_data['sample'], input_num, z_mean=z_mean, z_std=z_std)
            y = tf.one_hot(batch_data['label'], 2)
            
            if specific_weight:
                yield x, y, batch_data['weight'].values
            else:
                yield x, y

        if balanced_batch:
            np.random.shuffle(sampler[0::2])
            np.random.shuffle(sampler[1::2])
        else:
            np.random.shuffle(sampler)


def resample(data, split, method, shuffle=True):
   
    data = data.reset_index(drop=True)
    label_counts = data['label'].value_counts()
    
    many_label, few_label = label_counts.index
    many_count, few_count = label_counts.values
    more_name, few_name = ('positive', 'negative') if many_label == 1 else ('negative', 'positive')
    many_samples, few_samples = data[data['label'] == many_label], data[data['label'] == few_label]

    if method == 'US1-1':
        data = pd.concat([few_samples, many_samples.sample(n=few_count)], axis=0)

    elif method == 'US1-2':
        sampled_size = few_count * 2
        if sampled_size <= many_count:
            data = pd.concat([few_samples, many_samples.sample(n=sampled_size)], axis=0)
        else:
            raise RuntimeError('The dataset ratio of %s to %s are 1:%g, which can not be undersampled to 1:2. '
                               'The US1-2 method will be skipped.' % (few_name, more_name, many_count / few_count))

    elif method == 'OS1-1' or method == 'BB':
        data = pd.concat([data, few_samples.sample(n=many_count - few_count, replace=True)], axis=0)

    elif method == 'OS1-2':
        sampled_size = many_count // 2 - few_count
        if sampled_size >= 0:
            data = pd.concat([data, few_samples.sample(n=sampled_size, replace=True)], axis=0)
        else:
            raise RuntimeError('The dataset ratio of %s to %s are 1:%g, which can not be oversampled to 1:2. '
                               'The OS1-2 method will be skipped.' % (few_name, more_name, many_count / few_count))

    elif method == 'DA':

        if not all(imghdr.what(_) is not None for _ in data['sample']):
            raise RuntimeError('The DA method can not apply non-image classification case.')

        save_dir = env.WORKSPACE['DA_images']

        for da_file in list_files(save_dir, cond=lambda x: x.startswith(split)):
            (save_dir / da_file).unlink()

        datagen = ImageDataGenerator(rotation_range=30,
                                     width_shift_range=0.2,
                                     height_shift_range=0.2,
                                     brightness_range=[0.8, 1.2],
                                     shear_range=0.15,
                                     zoom_range=[0.75, 1.25],
                                     horizontal_flip=True,
                                     vertical_flip=True)
        generator = datagen.flow_from_dataframe(dataframe=few_samples, x_col='sample',
                                                target_size=load_img(data['sample'].values[0]).size[::-1],
                                                batch_size=1, class_mode=None,
                                                save_to_dir=save_dir, save_prefix=split)
       
        print('Selected balancing method: data augmentation. '
              'Now processing the image augmentation. This may take several minutes...')
        for i in range(many_count - few_count):
            next(generator)

        for da_file in list_files(save_dir, cond=lambda x: x.startswith(split)):
            data.loc[len(data)] = [str(save_dir / da_file), few_label, np.nan]
   
    if shuffle:
        data = data.sample(frac=1).reset_index(drop=True)
   
    return data


def make_data_batch(sample_paths, input_num, z_mean=None, z_std=None):

    if z_mean is None:
        z_mean = 0.0 if input_num == 1 else [0.0] * input_num
    if z_std is None:
        z_std = 1.0 if input_num == 1 else [1.0] * input_num

    batch = list() if input_num == 1 else [list() for _ in range(input_num)]

    for path in sample_paths:

        if imghdr.what(path) is not None:
            arr = img_to_array(load_img(path))
        elif str(path).endswith(('.pk', '.pickle')):
            with open(path, 'rb') as f:
                arr = pickle.load(f)
        else:
            raise IOError('the file %s cannot be processed' % path)

        if input_num == 1:
            batch.append(arr)
        else:
            for i in range(input_num):
                batch[i].append(arr[i])

    if input_num == 1:
        batch = (np.array(batch, dtype=float) - z_mean) / z_std
    else:
        for i in range(input_num):
            batch[i] = (np.array(batch[i], dtype=float) - z_mean[i]) / z_std[i]

    return batch
