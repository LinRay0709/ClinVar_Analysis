import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import argparse
import json
import pickle
import warnings

from training_func.custom_model import load_model
from training_func.data_processor import load_data
from training_func.plotting import draw_curve
from training_func.utils import *
import training_func.environment as env


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-name', default='RDDL', type=case_name,
                        help='name of this deep learning project')
    parser.add_argument('-t', default=3, type=int, help=argparse.SUPPRESS, dest='top')
    parser.add_argument('-n', action='store_true', help=argparse.SUPPRESS, dest='with_baseline')
    parser.add_argument('-g', default=0, type=int,
                        help='the specific GPU number', dest='gpu_number')
    parser.add_argument('-v', '--verbose', default=1, type=int, choices=[0, 1, 2],
                        help='verbosity mode')
    a = parser.parse_args()

    env.environment_setup(a.name, gpu_device='CUDA', gpu_number=a.gpu_number, quiet_mode=(a.verbose == 0))
    en_file = env.WORKSPACE['models'] / 'ensemble_info.json'

    with open(en_file, 'r') as f:
        en_info = json.load(f)

    ensemble_methods = [_ for _ in balancing_methods(baseline=a.with_baseline) if _ in en_info]
    model_paths = dict()

    for method in ensemble_methods.copy():
        cv_folds = en_info[method]['cv_folds']
        best_fold_num = np.argmin([_['val_rank'] for _ in cv_folds]) + 1
        model_path = env.WORKSPACE['models'] / ('%s_model_%d.h5' % (method, best_fold_num))
        if model_path.exists():
            model_paths[method] = (model_path, cv_folds[best_fold_num - 1]['z_mean'], cv_folds[best_fold_num - 1]['z_std'])
        else:
            warnings.warn('%s model does not exist.' % model_path)
            ensemble_methods.remove(method)

    if a.with_baseline:
        if 'Baseline' in ensemble_methods:
            ensemble_methods.remove('Baseline')
        else:
            parser.error('argument -n: specify the comparison with baseline but it does not exist.')

    if ensemble_methods:
        if len(ensemble_methods) < a.top:
            warnings.warn('Less than %d models were given; only ensemble the given models.' % a.top)
            a.top = len(ensemble_methods)
    else:
        parser.error('argument -name: the balancing model %s was not trained.' % a.name)

    # ensemble evaluation
    ensemble_predictions = list()
    print('\nValidly implemented model(s): %s' % ', '.join(ensemble_methods))
    print('Evaluate on ensemble validation set.\n')

    for method in ensemble_methods:

        model = load_model(model_paths[method][0])
        input_num = len(model.inputs)

        print('%s %s' % (a.name, method))
        pred_gen, pred_steps, pred_labels = load_data('ensemble', None, input_num,
                                                      z_mean=model_paths[method][1], z_std=model_paths[method][2])
        pred_result = model.predict(pred_gen, steps=pred_steps, verbose=a.verbose)
        ensemble_predictions.append({'y_true': pred_labels, 'y_pred': pred_result[:, 1].tolist()})

        del model
        env.memory_recovery()

    en_val_ranks = performance_rank(ensemble_predictions)
    ensemble_methods = [ensemble_methods[_] for _ in en_val_ranks]
    ensemble_predictions = [ensemble_predictions[_] for _ in en_val_ranks]

    with open(en_file, 'w') as f:
        for i, method in enumerate(ensemble_methods):
            en_info[method]['en_val_rank'] = i + 1
            en_info[method]['selected'] = (i + 1 <= a.top)
        f.write(json.dumps(en_info, sort_keys=False, indent=4))

    # test evaluation
    test_predictions = list()
    test_methods = ensemble_methods[:a.top]
    print('\nTop %d models: %s' % (a.top, ', '.join(test_methods)))
    print('Evaluate on test set.\n')

    if a.with_baseline:
        test_methods += ['Baseline']

    for method in test_methods:

        model = load_model(model_paths[method][0])
        input_num = len(model.inputs)

        print('%s %s' % (a.name, method))
        pred_gen, pred_steps, pred_labels = load_data('test', None, input_num,
                                                      z_mean=model_paths[method][1], z_std=model_paths[method][2])
        pred_result = model.predict(pred_gen, steps=pred_steps, verbose=a.verbose)
        test_predictions.append({'y_true': pred_labels, 'y_pred': pred_result[:, 1].tolist()})

        del model
        env.memory_recovery()

    if a.with_baseline:
        combined_predictions = [test_predictions.pop()]
        combined_labels = ['Ensemble Model', 'Baseline']
    else:
        combined_predictions = list()
        combined_labels = ['Ensemble Model']

    combined_predictions.insert(0, {'y_true': test_predictions[0]['y_true'],
                                    'y_pred': np.mean([_['y_pred'] for _ in test_predictions], axis=0)})

    print()
    performance = [metric_scores(**_) for _ in ensemble_predictions]
    tabulate([_.values() for _ in performance],
             headers=[to_formal(_) for _ in performance[0].keys()],
             v_headers=ensemble_methods,
             title='ENSEMBLE VALIDATION PERFORMANCE METRICS', float_fmt='%.3f')

    print()
    performance = [metric_scores(**_) for _ in combined_predictions]
    tabulate([_.values() for _ in performance],
             headers=[to_formal(_) for _ in performance[0].keys()],
             v_headers=combined_labels,
             title='TEST PERFORMANCE METRICS', float_fmt='%.3f')

    for curve_type in ['ROC', 'PRC']:
        draw_curve(env.WORKSPACE['other_metrics'] / ('ensemble_validation_%s.png' % curve_type),
                   ensemble_predictions, curve_type, labels=ensemble_methods,
                   title='%s Ensemble Validation %s' % (a.name, curve_type))
        draw_curve(env.WORKSPACE['plots'] / ('test_%s.png' % curve_type),
                   combined_predictions, curve_type, labels=combined_labels,
                   title='%s Test %s' % (a.name, curve_type))
