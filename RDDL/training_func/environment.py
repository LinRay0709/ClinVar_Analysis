import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from pathlib import Path
import gc
import numpy as np
import sys
import tensorflow as tf
import warnings

from training_func.utils import list_files

WORKSPACE = None
K_FOLD = 0

def environment_setup(name, gpu_device=None, gpu_number=-1, quiet_mode=False):

    global WORKSPACE, K_FOLD

    WORKSPACE = dict()
    directories = [('USER_pos', ),
                   ('USER_neg', ),
                   ('USER_hyperparameters', ),
                   ('RDDL_splitting_info', ),
                   ('RDDL_outputs', ),
                   ('RDDL_outputs', 'models'),
                   ('RDDL_outputs', 'plots'),
                   ('RDDL_outputs', 'plots', 'other_metrics'),
                   ('RDDL_outputs', 'tmp'),
                   ('RDDL_outputs', 'tmp', 'DA_images')]
    
    for d in directories:
        path = Path(name).joinpath(*d)
        path.mkdir(exist_ok=True)
        WORKSPACE[d[-1]] = path

    (WORKSPACE['models'] / 'ensemble_info.json').touch()

    K_FOLD = len(list_files(WORKSPACE['RDDL_splitting_info'], cond=lambda x: x.startswith('training')))
    
    if quiet_mode:
        warnings.filterwarnings('ignore')
        sys.stdout = open(os.devnull, 'a')
        # sys.stderr = open(os.devnull, 'a')

    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    np.seterr(invalid='ignore')
    # tf.config.run_functions_eagerly(True)

    if gpu_number == -1:
        tf.config.threading.set_intra_op_parallelism_threads(8)
        tf.config.threading.set_inter_op_parallelism_threads(8)
    elif gpu_device is not None:
        if gpu_device == 'CUDA':
            os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_number)
        else:
            raise ValueError('unsupported GPU device name: %s' % gpu_device)
        config = tf.compat.v1.ConfigProto(gpu_options=tf.compat.v1.GPUOptions(allow_growth=True))
        sess = tf.compat.v1.Session(config=config)


def memory_recovery():
        
    gc.collect()
    tf.keras.backend.clear_session()

    # TF graph isn't same as Keras graph
    tf.compat.v1.reset_default_graph()
