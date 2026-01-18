from datetime import datetime
from email import encoders
from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
from sklearn.metrics import auc, precision_recall_curve, roc_curve
import argparse
import email
import numpy as np
import os
import smtplib


def case_name(name):

    name = os.path.basename(os.path.abspath(name))

    if os.path.isdir(name):
        return name
    else:
        raise argparse.ArgumentTypeError('name of working directory %s does not exist.' % name)


def data_size(size):

    try:
        size = float(size)
    except ValueError as e:
        raise argparse.ArgumentTypeError(e)

    if 0 < size < 1:
        return size
    else:
        raise argparse.ArgumentTypeError('the percentage for data set splitting must be between 0.0 and 1.0.')


def get_nowtime():

    return datetime.now().strftime('%Y%m%d%H%M%S')


def list_files(dir_path, cond=None):

    return list(filter(cond, next(os.walk(dir_path))[2]))


def balancing_methods(resampled=None, baseline=True):
    
    resampled_methods = ['US1-1', 'US1-2', 'OS1-1', 'OS1-2', 'BB', 'DA']
    unresampled_methods = ['CW', 'FL-gamma-1', 'FL-gamma-2', 'MFE', 'SW']

    if resampled is None:
        selection = ['Baseline'] + resampled_methods + unresampled_methods
    else:
        selection = resampled_methods if resampled else ['Baseline'] + unresampled_methods

    if not baseline and 'Baseline' in selection:
        selection.remove('Baseline')

    return selection


def refine_history(history):

    history = history.copy()
    prefixes = ['', 'val_', 'samp_val_'] if 'samp_val_loss' in history else ['', 'val_']

    for prefix in prefixes:
        precision = np.array(history[prefix + 'precision'])
        recall = np.array(history[prefix + 'recall'])
        history[prefix + 'f1-score'] = np.nan_to_num((2 * precision * recall) / (precision + recall)).tolist()

    return history


def metric_scores(y_true, y_pred):

    y_pred_class = np.around(y_pred)
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    precision, recall, _ = precision_recall_curve(y_true, y_pred)

    return {'accuracy': accuracy_score(y_true, y_pred_class),
            'specificity': recall_score(y_true, y_pred_class, pos_label=0),
            'precision': precision_score(y_true, y_pred_class, zero_division=0),
            'recall': recall_score(y_true, y_pred_class),
            'f1-score': f1_score(y_true, y_pred_class, zero_division=0),
            'auROC': auc(fpr, tpr),
            'auPRC': auc(recall, precision)}


def performance_rank(predictions, keys=['auROC', 'auPRC']):

    criterion = list()
    
    for prediction in predictions:
        performance = metric_scores(**prediction)
        criterion.append([performance[_] for _ in keys])
    
    keys = ['sum'] + keys
    criterion = np.argsort(criterion, axis=0).argsort(axis=0)
    criterion = [tuple(_) for _ in np.append(criterion.sum(axis=1)[:, None], criterion, axis=1)]
    criterion = np.array(criterion, dtype=[(_, int) for _ in keys])
    
    return np.argsort(criterion, order=keys)[::-1].tolist()


def to_formal(name):

    formal_names = {'loss': 'Loss',
                    'accuracy': 'Accuracy',
                    'specificity': 'Specificity',
                    'precision': 'Precision',
                    'recall': 'Recall',
                    'f1-score': 'F1-Score'}

    return name if name not in formal_names else formal_names[name]


def tabulate(values, headers=None, v_headers=None, title=None, float_fmt='%.3f'):

    values = [list(_) for _ in values]

    if headers is not None:
        headers = list(headers)
        item_widths = [len(_) for _ in headers]
    else:
        item_widths = [0 for _ in range(len(values[0]))]

    if v_headers is not None:
        if headers is not None:
            headers.insert(0, '')
        for i, row in enumerate(values):
            row.insert(0, v_headers[i])
        item_widths.insert(0, 0)
    
    for row in values:
        for i, v in enumerate(row):
            row[i] = float_fmt % v if isinstance(v, float) else str(v)
            item_widths[i] = max(item_widths[i], len(row[i]))

    sep_line = '+%s+' % '+'.join('-' * (w + 2) for w in item_widths)

    if title is not None:
        print('+%s+' % ('-' * (len(sep_line) - 2)))
        print('| %%-%ds |' % (len(sep_line) - 4) % title)
    
    print(sep_line)

    if headers is not None:
        for i, h in enumerate(headers):
            print('| %%%ds ' % item_widths[i] % h, end='')
        print('|')
        print(sep_line)

    for row in values:
        for i, v in enumerate(row):
            print('| %%%ds ' % item_widths[i] % v, end='')
        print('|')
    print(sep_line)


def send_email(account, password, recievers, attached_files, title='No Title', content=''):

    msg = MIMEMultipart()

    for i, attached_file in enumerate(attached_files):

        with open(attached_file, 'rb') as f:

            file_name = os.path.basename(attached_file)

            mime = MIMEBase('image', 'png', filename=file_name)
            mime.add_header('Content-Disposition', 'attachment', filename=file_name)
            mime.add_header('Content-ID', '<%d>' % i)
            mime.add_header('X-Attachment-Id', '%d' % i)
            mime.set_payload(f.read())

        email.encoders.encode_base64(mime)
        msg.attach(mime)

    msg.attach(MIMEText('<html><body><div>%s</div></body></html>' % content.replace('\n', '<br>'), 'html', 'utf-8'))

    msg['Subject'] = title
    msg['From'] = account
    msg['To'] = ','.join(recievers)

    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.ehlo()
    server.starttls()
    server.login(account, password)
    server.send_message(msg, account, recievers)
    server.quit()
