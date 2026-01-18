from sklearn.metrics import auc, precision_recall_curve, roc_curve
import math
import matplotlib.pyplot as plt
import numpy as np

from training_func.utils import to_formal
import training_func.environment as env


# Draw the picturer of k-fold trainning results and output the filename of the picture
def draw_metric(image_path, histories, metric, title='Learning Curve'):

    epochs = len(histories[0][metric])
    x_ticks = np.arange(1, epochs + 1)

    top_value = 0.0
    properties = {metric: ('Training', 'lightskyblue', 'teal'),
                  'val_%s' % metric: ('Validation', 'pink', 'salmon'),
                  'samp_val_%s' % metric: ('Sampled Validation', 'gold', 'orange')}

    plt.figure(figsize=(10, 5))

    for metric_key, (label_name, color, markercolor) in properties.items():

        if metric_key in histories[0]:

            stats = np.array([h[metric_key] for h in histories])
            min_stats = np.min(stats, axis=0)
            max_stats = np.max(stats, axis=0)
            mean_stats = np.mean(stats, axis=0)
    
            plt.plot(x_ticks, min_stats, '--', color=color, linewidth=1.5, alpha=0.6)
            plt.plot(x_ticks, max_stats, '--', color=color, linewidth=1.5, alpha=0.6)
            plt.fill_between(x_ticks, min_stats, max_stats, color=color, alpha=0.2)
            plt.plot(x_ticks, mean_stats, '-o', color=color, linewidth=2,
                     markeredgecolor=markercolor, markeredgewidth=2, markersize=4,
                     label=r'%d-Fold %s Average (%.3f $\pm$ %.3f)' % (env.K_FOLD, label_name, mean_stats[-1], stats.std(axis=0)[-1]))

            top_value = max(top_value, np.max(max_stats) * 1.5)

    plt.grid(color='gainsboro', linestyle='--', linewidth=1)
    plt.xticks(np.r_[np.arange(1, epochs, math.ceil(epochs / 10)), epochs], fontsize=14)

    if metric == 'loss':
        plt.ylim(0.0, top_value)
        plt.yticks(np.around(np.linspace(0.0, top_value, 6), 3), fontsize=14)
    else:
        plt.ylim(0.0, 1.0)
        plt.yticks(np.linspace(0.0, 1.0, 6), fontsize=14)
        plt.gca().set_yticklabels(['%d%%' % (_ * 100) for _ in plt.gca().get_yticks()])

    plt.title(title, fontsize=16, pad=10)
    plt.xlabel('Epoch', fontsize=16, labelpad=10)
    plt.ylabel(to_formal(metric), fontsize=16, labelpad=10)
    plt.legend(fontsize=12)
    plt.tight_layout(pad=2)

    plt.savefig(image_path, dpi=300)
    plt.close()


# Draw ROC or PR curve and calculate area
def draw_curve(image_path, predictions, curve_type, labels=None, title='Learning Curve'):

    y_interps = list()
    areas = list()
    baseline = np.linspace(0, 1, 3500)
    plt.figure(figsize=(8, 8))

    if curve_type == 'ROC':
        plt.plot([0, 1], [0, 1], '--', color='gainsboro', linewidth=1)

    # Calculate & interpolate
    for prediction in predictions:

        if curve_type == 'ROC':
            fpr, tpr, _ = roc_curve(prediction['y_true'], prediction['y_pred'])
            y_interp = np.interp(baseline, fpr, tpr)
        elif curve_type == 'PRC':
            precision, recall, _ = precision_recall_curve(prediction['y_true'], prediction['y_pred'])
            y_interp = np.interp(baseline, recall[::-1], precision[::-1])

        y_interps.append(y_interp)
        areas.append(auc(baseline, y_interp))
    
    if labels is None:

        lower_y_interp = np.min(y_interps, axis=0)
        upper_y_interp = np.max(y_interps, axis=0)
        mean_y_interp = np.mean(y_interps, axis=0) 

        plt.plot(baseline, lower_y_interp, '--', color='orange', linewidth=1.5, alpha=0.6)
        plt.plot(baseline, upper_y_interp, '--', color='orange', linewidth=1.5, alpha=0.6)
        plt.fill_between(baseline, lower_y_interp, upper_y_interp, color='orange', alpha=0.2)
        plt.plot(baseline, mean_y_interp, color='orangered', linewidth=2,
                 label=r'%d-Fold Validation Average (au%s = %.3f $\pm$ %.3f)' % (env.K_FOLD, curve_type, auc(baseline, mean_y_interp), np.std(areas)))

    else:

        colors = plt.cm.get_cmap('tab10').colors + plt.cm.get_cmap('Accent').colors
        for i in range(len(predictions)):
            plt.plot(baseline, y_interps[i], color=colors[i], zorder=-i,
                     label='%s (au%s = %.3f)' % (labels[i], curve_type, areas[i]))
    
    # Some specific edit about the roc picture
    if curve_type == 'ROC':
        plt.xlabel('1 - Specificity', fontsize=16, labelpad=10)
        plt.ylabel('Sensitivity', fontsize=16, labelpad=10)
    elif curve_type == 'PRC':
        plt.xlabel('Recall', fontsize=16, labelpad=10)
        plt.ylabel('Precision', fontsize=16, labelpad=10)

    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.xticks(np.arange(0.0, 1.2, 0.2), fontsize=14)
    plt.yticks(np.arange(0.2, 1.2, 0.2), fontsize=14)
    plt.gca().set_xticklabels(['%d%%' % (_ * 100) for _ in plt.gca().get_xticks()])
    plt.gca().set_yticklabels(['%d%%' % (_ * 100) for _ in plt.gca().get_yticks()])

    plt.title(title, fontsize=16, pad=10)
    plt.legend(fontsize=12)
    plt.tight_layout(pad=2)

    plt.savefig(image_path, dpi=300)
    plt.close()
