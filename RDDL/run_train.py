import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import argparse
import importlib
import json
import numpy as np
import pickle
import random
import tensorflow as tf

from training_func.callbacks import *
from training_func.data_processor import load_data
from training_func.plotting import *
from training_func.utils import *
import training_func.environment as env


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-name', default='RDDL', type=case_name,
                        help='name of this deep learning project')
    parser.add_argument('-m', choices=balancing_methods(), required=True,
                        help='the balanced method', metavar='METHOD', dest='method')
    parser.add_argument('-nor', default=0, type=int, choices=[0, 1, 2],
                        help='the specified normalization option of this balancing scheme', dest='nor_mode')
    parser.add_argument('-g', default=0, type=int,
                        help='the specific GPU number', dest='gpu_number')
    parser.add_argument('-v', '--verbose', default=1, type=int, choices=[0, 1, 2],
                        help='verbosity mode')
    a = parser.parse_args()

    env.environment_setup(a.name, gpu_device='CUDA', gpu_number=a.gpu_number, quiet_mode=(a.verbose == 0))
    user_model = importlib.import_module('%s.USER_model' % a.name)

    with open(env.WORKSPACE['USER_hyperparameters'] / ('%s.json' % a.method), 'r') as f:
        p = argparse.Namespace(**json.load(f))

    # save the history in order to draw the learning curve
    z_stats, histories, predictions = list(), list(), list()

    # start training five folds
    for fold_num in range(1, env.K_FOLD + 1):

        tf.random.set_seed(p.random_seed)
        np.random.seed(p.random_seed)
        random.seed(p.random_seed)

        model = user_model.self_defined_model(p.dropout_rate, model_name='%s_%s_%d' % (a.name, a.method, fold_num))
        model.summary()
        input_num = len(model.inputs)
        
        train_gen, train_steps, _, z_mean, z_std, alpha = load_data('training', fold_num, input_num,
                                                                    batch_size=p.batch_size, method=a.method,
                                                                    nor_mode=a.nor_mode)
        val_gen, val_steps, _ = load_data('validation', fold_num, input_num,
                                          batch_size=p.batch_size, z_mean=z_mean, z_std=z_std)
        val_evaluater = Evaluater('val', val_gen, val_steps, p.batch_size, verbose=a.verbose)

        scheduler = LearningRateScheduler(p.epochs, p.initial_lr, p.decay_rate, p.bottom_lr)
        loss = SamplingLoss(a.method, alpha=alpha)

        if a.method in balancing_methods(resampled=True):
            samp_val_gen, samp_val_steps, _ = load_data('validation', fold_num, input_num,
                                                        batch_size=p.batch_size, method=a.method,
                                                        z_mean=z_mean, z_std=z_std)
            samp_val_evaluater = Evaluater('samp_val', samp_val_gen, samp_val_steps,
                                           p.batch_size, verbose=a.verbose)
            callbacks = [val_evaluater, samp_val_evaluater, scheduler.scheduler()]
        else:
            callbacks = [val_evaluater, scheduler.scheduler()]

        print('\n%s %s (Fold %d/%d)' % (a.name, a.method, fold_num, env.K_FOLD))
        model.compile(optimizer=tf.keras.optimizers.Adam(), loss=loss.sampling_loss())
        history = model.fit(train_gen, epochs=p.epochs, batch_size=p.batch_size,
                            steps_per_epoch=train_steps, verbose=a.verbose, callbacks=callbacks)

        # save history & the weight of the model
        z_stats.append((z_mean, z_std))
        histories.append(refine_history(history.history))
        model.save(env.WORKSPACE['models'] / ('%s_model_%d.h5' % (a.method, fold_num)))

        # start predicting five folds
        pred_gen, pred_steps, pred_labels = load_data('validation', fold_num, input_num,
                                                      z_mean=z_mean, z_std=z_std)
        pred_result = model.predict(pred_gen, steps=pred_steps, verbose=a.verbose)
        predictions.append({'y_true': pred_labels, 'y_pred': pred_result[:, 1].tolist()})
       
        del model
        env.memory_recovery()

    with open(env.WORKSPACE['tmp'] / ('%s_histories.pickle' % a.method), 'wb') as f:
        pickle.dump(histories, f)
    with open(env.WORKSPACE['tmp'] / ('%s_predictions.pickle' % a.method), 'wb') as f:
        pickle.dump(predictions, f)
    
    with open(env.WORKSPACE['models'] / 'ensemble_info.json', 'r+') as f:
        
        try:
            en_info = json.load(f)
        except:
            en_info = dict()

        f.truncate(0)
        f.seek(0)
        
        if a.method not in en_info:
            en_info[a.method] = dict()

        en_info[a.method]['data_time'] = get_nowtime()
       
        for method in en_info:
            en_info[method]['en_val_rank'] = None
            en_info[method]['selected'] = None

        en_info[a.method]['cv_folds'] = list()
        val_ranks = performance_rank(predictions)
        
        for i in range(env.K_FOLD):
            en_info[a.method]['cv_folds'].append({'val_rank': val_ranks.index(i) + 1,
                                                  'z_mean': z_stats[i][0],
                                                  'z_std': z_stats[i][1]})
        f.write(json.dumps(en_info, sort_keys=False, indent=4))
    
    draw_metric(env.WORKSPACE['plots'] / ('%s_learning_curve_f1-score.png' % a.method),
                histories, 'f1-score',
                title='%s %s %d-Fold Cross-Validation' % (a.name, a.method, env.K_FOLD))
    draw_curve(env.WORKSPACE['plots'] / ('%s_validation_ROC.png' % a.method),
               predictions, 'ROC',
               title='%s %s %d-Fold Cross-Validation ROC' % (a.name, a.method, env.K_FOLD))
    draw_curve(env.WORKSPACE['plots'] / ('%s_validation_PRC.png' % a.method),
               predictions, 'PRC',
               title='%s %s %d-Fold Cross-Validation PRC' % (a.name, a.method, env.K_FOLD))
    
    for metric in ['loss', 'accuracy', 'specificity', 'precision', 'recall', 'auROC']:
        draw_metric(env.WORKSPACE['other_metrics'] / ('%s_learning_curve_%s.png' % (a.method, metric)),
                    histories, metric,
                    title='%s %s %d-Fold Cross-Validation' % (a.name, a.method, env.K_FOLD))

    print()
    performance = [metric_scores(**_) for _ in predictions]
    tabulate([_.values() for _ in performance],
             headers=[to_formal(_) for _ in performance[0]],
             v_headers=['Fold %d' % (_ + 1) for _ in range(env.K_FOLD)],
             title='%d-FOLD CROSS-VALIDATION PERFORMANCE METRICS' % env.K_FOLD, float_fmt='%.3f')
