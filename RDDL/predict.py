import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

from training_func.custom_model import load_model
from training_func.data_processor import make_data_batch
from training_func.utils import case_name, list_files
import training_func.environment as env


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-name', default='RDDL', type=case_name,
                        help='name of this deep learning project')
    parser.add_argument('-i', type=Path, required=True,
                        help='the input directory path', dest='in_dir')
    parser.add_argument('-o', required=True,
                        help='the output file path', dest='out_file')
    parser.add_argument('-g', default=0, type=int,
                        help='the specific GPU number', dest='gpu_number')
    parser.add_argument('-v', '--verbose', default=1, type=int, choices=[0, 1, 2],
                        help='verbosity mode')
    a = parser.parse_args()

    env.environment_setup(a.name, gpu_device='CUDA', gpu_number=a.gpu_number, quiet_mode=(a.verbose == 0))

    samples = list_files(a.in_dir)
    model_paths = dict()
    y_pred = list()

    with open(env.WORKSPACE['models'] / 'ensemble_info.json', 'r') as f:
        
        en_info = json.load(f)
        test_methods = sorted([_ for _ in en_info if en_info[_]['selected']],
                              key=lambda x: en_info[x]['en_val_rank'])
        
        if not test_methods:
            parser.error('argument -name: case %s has not yet run the ensemble step' % a.name)
        
        for method in test_methods:
            cv_folds = en_info[method]['cv_folds']
            best_fold_num = np.argmin([_['val_rank'] for _ in cv_folds]) + 1
            model_paths[method] = (env.WORKSPACE['models'] / ('%s_model_%d.h5' % (method, best_fold_num)),
                                   cv_folds[best_fold_num - 1]['z_mean'],
                                   cv_folds[best_fold_num - 1]['z_std'])

    for i, (method, model_path) in enumerate(model_paths.items()):

        model = load_model(model_path[0])
        input_num = len(model.inputs)

        print('Prediction - %d/%d' % (i + 1, len(model_paths)))
        x = (make_data_batch([a.in_dir / _], input_num, z_mean=model_path[1], z_std=model_path[2]) for _ in samples)
        pred_result = model.predict(x, steps=len(samples), verbose=a.verbose)
        y_pred.append(pred_result[:, 1].tolist())

        del model
        env.memory_recovery()

    out_df = pd.DataFrame({'sample': samples, 'probability': np.mean(y_pred, axis=0)})
    out_df.to_csv(a.out_file, index=False)
