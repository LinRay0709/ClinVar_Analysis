from tensorflow.keras.metrics import *
import tensorflow as tf


def load_model(file_path):

    return tf.keras.models.load_model(file_path,
                                      custom_objects={'CustomModel': CustomModel},
                                      compile=False)


class CustomModel(tf.keras.Model):

    def __init__(self, **kwargs):

        super(CustomModel, self).__init__(**kwargs)

        self.dr_loss_metric = Mean(name='loss_dropout')
        self.dr_acc_metric = BinaryAccuracy(name='accuracy_dropout')

        self.loss_metric = Mean(name='loss')
        self.acc_metric = BinaryAccuracy(name='accuracy')
        self.spe_metric = Recall(name='specificity')
        self.prec_metric = Precision(name='precision')
        self.rec_metric = Recall(name='recall')
        self.auc_metric = AUC(name='auROC', curve='ROC')

        self.tp_metric = TruePositives(name='TP')
        self.fp_metric = FalsePositives(name='FP')
        self.fn_metric = FalseNegatives(name='FN')
        self.tn_metric = TrueNegatives(name='TN')

    def test_step(self, data):

        if len(data) == 3:
            x, y, w = data
            w = tf.divide(w, tf.reduce_mean(w))
        else:
            x, y = data
            w = None

        # Compute predictions
        y_pred = self(x, training=False)
        loss = self.compiled_loss(y, y_pred, sample_weight=w,
                                  regularization_losses=self.losses)

        # Update the metrics.
        self.dr_loss_metric.update_state(loss)
        self.dr_acc_metric.update_state(y, y_pred)
        self.loss_metric.update_state(loss)

        self._custom_eval(y, y_pred)

        return {m.name: m.result() for m in self.metrics}

    def train_step(self, data):

        if len(data) == 3:
            x, y, w = data
            w = tf.divide(w, tf.reduce_mean(w))
        else:
            x, y = data
            w = None

        with tf.GradientTape() as tape:
            y_pred = self(x, training=True)
            loss = self.compiled_loss(y, y_pred, sample_weight=w,
                                      regularization_losses=self.losses)

        # Compute gradients & update weights
        trainable_vars = self.trainable_variables
        gradients = tape.gradient(loss, trainable_vars)
        self.optimizer.apply_gradients(zip(gradients, trainable_vars))
        
        # Forward pass with dropout
        self.dr_loss_metric.update_state(loss)
        self.dr_acc_metric.update_state(y, y_pred)

        # Forward pass without dropout
        y_pred_no_reg = self(x, training=False)
        loss_no_reg = self.compiled_loss(y, y_pred_no_reg, sample_weight=w,
                                         regularization_losses=self.losses)
        self.loss_metric.update_state(loss_no_reg)

        self._custom_eval(y, y_pred_no_reg)

        return {m.name: m.result() for m in self.metrics}

    @property
    def metrics(self):

        return [self.dr_loss_metric, self.dr_acc_metric, self.loss_metric, self.acc_metric,
                self.spe_metric, self.prec_metric, self.rec_metric, self.auc_metric,
                self.tp_metric, self.fp_metric, self.fn_metric, self.tn_metric]
     
    #### evaluate specified metrics
    def _custom_eval(self, y_true_one_hot, y_pred_one_hot):

        y_true = tf.math.argmax(y_true_one_hot, axis=-1)
        y_pred = tf.math.argmax(y_pred_one_hot, axis=-1)
        y_label = tf.gather(y_true_one_hot, 1, axis=-1)
        y_prob = tf.gather(y_pred_one_hot, 1, axis=-1)

        self.acc_metric.update_state(y_true, y_pred)
        self.spe_metric.update_state(1 - y_true, 1 - y_pred)
        self.prec_metric.update_state(y_true, y_pred)
        self.rec_metric.update_state(y_true, y_pred)
        self.auc_metric.update_state(y_label, y_prob)
        self.tp_metric.update_state(y_true, y_pred)
        self.fp_metric.update_state(y_true, y_pred)
        self.fn_metric.update_state(y_true, y_pred)
        self.tn_metric.update_state(y_true, y_pred)
