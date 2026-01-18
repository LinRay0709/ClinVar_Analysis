import math
import numpy as np
import tensorflow as tf

from training_func.environment import memory_recovery


class Evaluater(tf.keras.callbacks.Callback):

    def __init__(self, prefix, generator, steps, batch_size, verbose=1):

        self.prefix = prefix
        self.generator = generator
        self.steps = steps
        self.batch_size = batch_size
        self.verbose = verbose

    def on_epoch_end(self, epoch, logs=None):
       
        stat = self.model.evaluate(self.generator, batch_size=self.batch_size,
                                   verbose=self.verbose, steps=self.steps)
        logs.update({'%s_%s' % (self.prefix, self.model.metrics_names[i]): s for i, s in enumerate(stat)})

        memory_recovery()


class LearningRateScheduler():

    def __init__(self, epochs, initial_lr, decay_rate, bottom_lr):

        self._initial_lr = initial_lr
        self._decay_rate = decay_rate
        self._bottom_lr = bottom_lr
        self._change_pt = math.ceil(epochs * 0.2)
        self._slow_pt = math.ceil(epochs * 0.9)

    def _cosine_decay(self, epoch):

        if epoch < self._change_pt:
            lr = self._initial_lr * np.cos(1 - (epoch / self._change_pt))
        elif epoch < self._slow_pt:
            lr = self._initial_lr * (self._decay_rate ** (epoch - self._change_pt))
        else:
            lr = self._bottom_lr

        print('LR: %.12f' % lr)

        return lr

    def scheduler(self):

        return tf.keras.callbacks.LearningRateScheduler(self._cosine_decay)


class SamplingLoss():

    def __init__(self, method, alpha=0.5):

        self._method = method
        self._alpha = alpha

        if method == 'CW':
            self._gamma = 0
        elif method == 'FL-gamma-1':
            self._gamma = 1
        elif method == 'FL-gamma-2':
            self._gamma = 2

    def class_weight(self, y_true, y_pred):

        y_true = tf.cast(y_true, tf.float32)
        y_pred = tf.clip_by_value(y_pred, tf.keras.backend.epsilon(), 1.0 - tf.keras.backend.epsilon())

        y_t = tf.reduce_sum(y_true * y_pred, axis=1)
        ce = -tf.keras.backend.log(y_t)

        alpha_t = tf.reduce_sum(y_true * [1.0 - self._alpha, self._alpha], axis=1)
        weight = tf.pow(1.0 - y_t, self._gamma) * alpha_t

        if self._method == 'CW':
            loss = tf.reduce_sum(weight * ce) / tf.reduce_sum(weight)
        else:
            loss = tf.reduce_mean(weight * ce)
        
        return loss

    # for MFE
    def mean_false_error(self, y_true, y_pred):

        y_true = tf.cast(y_true, tf.float32)
        y_pred = tf.cast(y_pred, tf.float32)

        # get label index & neural == 2
        index_1 = tf.where(tf.equal(0.0, y_true[:, 0]))
        index_0 = tf.where(tf.equal(1.0, y_true[:, 0]))

        # split positive set and negative set
        y_t_1 = tf.gather(y_true, indices=index_1)
        y_p_1 = tf.gather(y_pred, indices=index_1)
        y_t_0 = tf.gather(y_true, indices=index_0)
        y_p_0 = tf.gather(y_pred, indices=index_0)
        loss = 0.5 * self.calculate_error(y_t_1, y_p_1) + 0.5 * self.calculate_error(y_t_0, y_p_0)

        return loss

    # for mfe calculate error
    def calculate_error(self, y_true, y_pred):

        y_true = tf.cast(y_true, tf.float32)
        y_pred = tf.cast(y_pred, tf.float32)
        
        sqr_error = tf.keras.backend.square(y_true - y_pred)
        sum_error_in_neural = tf.reduce_sum(tf.keras.backend.mean(sqr_error, axis=0))
        sum_error_in_neural = tf.where(tf.math.is_nan(sum_error_in_neural),
                                       tf.zeros_like(sum_error_in_neural),
                                       sum_error_in_neural)

        return sum_error_in_neural

    def sampling_loss(self):

        if self._method in ['CW', 'FL-gamma-1', 'FL-gamma-2']:
            return self.class_weight
        elif self._method == 'MFE':
            return self.mean_false_error
        else:
            return tf.keras.metrics.categorical_crossentropy
