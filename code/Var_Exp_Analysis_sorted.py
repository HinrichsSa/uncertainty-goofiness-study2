import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import os as os
from scipy.optimize import differential_evolution, minimize
from scipy.stats import logistic
from scipy import stats
import multiprocessing as mp
import sys as sys
import seaborn as sbn


# NOTE: Adapt this path accordingly
os.chdir('D:\Studium\Auslandsstudium\TuitionWaver_Master\Masterthesis\Analysis\Exp_Variance\Modelling\sorted\code')


# NOTE: Create Figure for illustration of hypotheses
def plot_hypothesis():
    f_name = '../fit_input/master_data.csv'

    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]

    d = pd.read_csv(f_name)
    rot0 = d[d['sub'] == 9]['Appl_Perturb'].values
    rot1 = d[d['sub'] == 1]['Appl_Perturb'].values
    rot2 = d[d['sub'] == 2]['Appl_Perturb'].values

    p1 = [0.2, 0.99, 25, 0.1, 0.3, 10]
    x_pred1 = simulate_state_space_with_g_func_2_state(p1, rot0)
    p2 = [0.05, 0.99, 25, 0.1, 0.3, 10]
    x_pred2 = simulate_state_space_with_g_func_2_state(p2, rot1)
    p3 = [0.01, 0.99, 25, 0.1, 0.3, 10]
    x_pred3 = simulate_state_space_with_g_func_2_state(p3, rot2)

    plt.figure(figsize=(10, 4))
    legend0 = mpatches.Patch(color=colors[0], label='0°sd group')
    legend1 = mpatches.Patch(color=colors[1], label='4°sd group')
    legend2 = mpatches.Patch(color=colors[2], label='12°sd group')
    plt.legend(handles=[legend0, legend1, legend2], loc="upper left")
    plt.plot(x_pred1[:, 5], alpha=0.4, c=colors[0])
    plt.plot(x_pred2[:, 5], alpha=0.4, c=colors[1])
    plt.plot(x_pred3[:, 5], alpha=0.4, c=colors[2])

    plt.ylabel('Error (in °)')
    plt.xlabel('Trial')
    plt.suptitle('Illustration Prediction Experimental Outcome', fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/illustrate_exp_design/predcition_exp_outcome.png')
    plt.close()


# NOTE: Simulate and fit the 2 state model
def simulate_state_space_with_g_func_2_state(p, rot):
    alpha_s = p[0]
    beta_s = p[1]
    g_sigma_s = p[2]
    alpha_f = p[3]
    beta_f = p[4]
    g_sigma_f = p[5]

    num_trials = rot.shape[0]

    delta = np.zeros(num_trials)

    theta_values = np.linspace(0.0, 330.0, 12) - 150.0
    theta_train_ind = np.where(theta_values == 0.0)[0][0]
    theta_ind = theta_train_ind * np.ones(num_trials, dtype=np.int8)

    def g_func(theta, theta_mu, sigma):
        if sigma != 0:
            G = np.exp(-(theta - theta_mu)**2 / (2 * sigma**2))
        else:
            G = np.zeros(12)
        return (G)

    x = np.zeros((12, num_trials))
    xs = np.zeros((12, num_trials))
    xf = np.zeros((12, num_trials))
    for i in range(0, num_trials - 1):
        if i > 1000 and i <= 1072:
            delta[i] = 0.0
        elif i > 1172:
            delta[i] = 0.0
        else:
            delta[i] = x[theta_ind[i], i] - rot[i]

        Gs = g_func(theta_values, theta_values[theta_ind[i]], g_sigma_s)
        Gf = g_func(theta_values, theta_values[theta_ind[i]], g_sigma_f)
        if np.isnan(rot[i]):
            xs[:, i + 1] = beta_s * xs[:, i]
            xf[:, i + 1] = beta_f * xf[:, i]
        else:
            xs[:, i + 1] = beta_s * xs[:, i] - alpha_s * delta[i] * Gs
            xf[:, i + 1] = beta_f * xf[:, i] - alpha_f * delta[i] * Gf

        x[:, i + 1] = xs[:, i + 1] + xf[:, i + 1]

    return (xs.T + xf.T)


def fit_obj_func_sse_2_state(params, *args):
    alpha_s = params[0]
    beta_s = params[1]
    g_sigma_s = params[2]
    alpha_f = params[3]
    beta_f = params[4]
    g_sigma_f = params[5]

    if alpha_s >= alpha_f:
        sse = 10**100
        r_squared = 2
    elif beta_s <= beta_f:
        sse = 10**100
        r_squared = 2
    else:
        x_obs = args[0]
        rot = args[1]
        x_pred = simulate_state_space_with_g_func_2_state(params, rot)

        sst_rec = np.zeros(12)
        sse_rec = np.zeros(12)
        for i in range(12):
            # # pick trial indices for cost function
            # fit_inds = np.concatenate(
            #     (np.arange(600, 700, 1), np.arange(1000, 1072, 1),
            #      np.arange(1172, 1272, 1)))
            # sse_rec[i] = (np.nansum((x_obs[fit_inds, i] - x_pred[fit_inds, i])**2))
            # sst_rec[i] = (np.nansum((x_obs[fit_inds, i] - np.nanmean(x_obs[fit_inds, i]))**2))
            sse_rec[i] = (np.nansum((x_obs[:, i] - x_pred[:, i])**2))
            sst_rec[i] = (np.nansum((x_obs[:, i] - np.nanmean(x_obs[:, i]))**2))

            sst = np.nansum(sst_rec)
            sse = np.nansum(sse_rec)

            r_squared = 1 - (sse / sst)

    return (sse, r_squared)
    # NOTE: Activate the return function below for the bootstrap fit function
    # return sse


def fit_state_space_with_g_func_2_state_grp_bootstrap_ee():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)

    n_boot_samp = 1000

    for i in d['cnd'].unique():
        for j in d['rot_dir'].unique():
            p_rec = -1 * np.ones((n_boot_samp, 6))
            for b in range(n_boot_samp):
                print(i, j, b)
                dd = d[(d['cnd'] == i) & (d['rot_dir'] == j)]
                rot = dd['Appl_Perturb'].values[0:1272]
                subs = dd['sub'].unique()
                boot_subs = np.random.choice(subs,
                                             size=subs.shape[0],
                                             replace=True)
                x_boot_rec = []
                for k in boot_subs:
                    x_boot_rec.append(d[d['sub'] == k])
                    x_boot = pd.concat(x_boot_rec)

                x_obs = x_boot.groupby(['cnd', 'rot_dir', 'Target',
                                        'trial']).mean()
                x_obs.reset_index(inplace=True)

                x_obs = x_obs[["Endpoint_Error", "target_deg", "trial"]]
                x_obs = x_obs.pivot(index="trial",
                                    columns="target_deg",
                                    values="Endpoint_Error")

                x_obs = x_obs.values

                args = (x_obs, rot)
                bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60))
                results = differential_evolution(func=fit_obj_func_sse_2_state,
                                                 bounds=bounds,
                                                 args=args,
                                                 maxiter=300,
                                                 disp=False,
                                                 polish=False,
                                                 updating="deferred",
                                                 workers=1)
                p = results["x"]
                p_rec[b, :] = p

                f_name_p = "../fits/fit_grp_2state_bootstrap_ee" + str(
                    i) + '_' + j
                with open(f_name_p, "w") as f:
                    np.savetxt(f, p_rec, "%0.4f", ",")

    return p_rec


def fit_state_space_with_g_func_2_state_grp_bootstrap_bcee():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)

    n_boot_samp = 1000

    for i in d['cnd'].unique():
        for j in d['rot_dir'].unique():
            p_rec = -1 * np.ones((n_boot_samp, 6))
            for b in range(n_boot_samp):
                print(i, j, b)
                dd = d[(d['cnd'] == i) & (d['rot_dir'] == j)]
                rot = dd['Appl_Perturb'].values[0:1272]
                subs = dd['sub'].unique()
                boot_subs = np.random.choice(subs,
                                             size=subs.shape[0],
                                             replace=True)
                x_boot_rec = []
                for k in boot_subs:
                    x_boot_rec.append(d[d['sub'] == k])
                    x_boot = pd.concat(x_boot_rec)

                x_obs = x_boot.groupby(['cnd', 'rot_dir', 'Target',
                                        'trial']).mean()
                x_obs.reset_index(inplace=True)

                x_obs = x_obs[["bcee", "target_deg", "trial"]]
                x_obs = x_obs.pivot(index="trial",
                                    columns="target_deg",
                                    values="bcee")

                x_obs = x_obs.values

                args = (x_obs, rot)
                bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60))
                results = differential_evolution(func=fit_obj_func_sse_2_state,
                                                 bounds=bounds,
                                                 args=args,
                                                 maxiter=300,
                                                 disp=False,
                                                 polish=False,
                                                 updating="deferred",
                                                 workers=1)
                p = results["x"]
                p_rec[b, :] = p

                f_name_p = "../fits/fit_grp_2state_bootstrap_bcee" + str(
                    i) + '_' + j
                with open(f_name_p, "w") as f:
                    np.savetxt(f, p_rec, "%0.4f", ",")

    return p_rec


def fit_state_space_with_g_func_2_state_grp_bcee():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)
    p_rec = np.empty((0, 6))

    rot = d[d['sub'] == 1]['Appl_Perturb'].values

    x_obs_all = d.groupby(['cnd', 'Target', 'trial', 'rot_dir']).mean()
    x_obs_all.reset_index(inplace=True)

    for i in range(x_obs_all['cnd'].unique().shape[0]):
        for j in ('cw', 'ccw'):
            x_obs = x_obs_all[x_obs_all['cnd'] == i]
            x_obs = x_obs[x_obs['rot_dir'] == j]
            x_obs = x_obs[["bcee", "Target", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="bcee")

            x_obs = x_obs.values

            args = (x_obs, rot)
            bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60))
            results = differential_evolution(func=fit_obj_func_sse_2_state,
                                             bounds=bounds,
                                             args=args,
                                             maxiter=300,
                                             disp=False,
                                             polish=True,
                                             updating="deferred",
                                             workers=1)
            p = results["x"]

            p_rec = np.append(p_rec, [p], axis=0)

            f_name_p = "../fits/fit_2state_grp_bcee" + str(i) + '_' + j
            with open(f_name_p, "w") as f:
                np.savetxt(f, p, "%0.4f", ",")

    return p_rec


def fit_state_space_with_g_func_2_state_grp_ee():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)
    p_rec = np.empty((0, 6))

    rot = d[d['sub'] == 1]['Appl_Perturb'].values

    x_obs_all = d.groupby(['cnd', 'Target', 'trial', 'rot_dir']).mean()
    x_obs_all.reset_index(inplace=True)

    for i in range(x_obs_all['cnd'].unique().shape[0]):
        for j in ('cw', 'ccw'):
            x_obs = x_obs_all[x_obs_all['cnd'] == i]
            x_obs = x_obs[x_obs['rot_dir'] == j]
            x_obs = x_obs[["Endpoint_Error", "Target", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="Endpoint_Error")

            x_obs = x_obs.values

            args = (x_obs, rot)
            bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60))
            results = differential_evolution(func=fit_obj_func_sse_2_state,
                                             bounds=bounds,
                                             args=args,
                                             maxiter=300,
                                             disp=False,
                                             polish=True,
                                             updating="deferred",
                                             workers=1)
            p = results["x"]

            p_rec = np.append(p_rec, [p], axis=0)

            f_name_p = "../fits/fit_2state_grp_ee" + str(i) + '_' + j
            with open(f_name_p, "w") as f:
                np.savetxt(f, p, "%0.4f", ",")

    return p_rec


def fit_state_space_with_g_func_2_state_ee():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)
    sub = d['sub'].unique()
    length_names = sub.shape[0]

    p_rec = np.empty((0, 6))
    for i in range(length_names):
        print(sub[i])
        x_obs = d[d['sub'] == sub[i]]
        rot = x_obs["Appl_Perturb"].values
        x_obs = x_obs[["Endpoint_Error", "Target", "target_deg", "trial"]]
        x_obs = x_obs.pivot(index="trial",
                            columns="target_deg",
                            values="Endpoint_Error")

        x_obs = x_obs.values

        args = (x_obs, rot)
        bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60))
        results = differential_evolution(func=fit_obj_func_sse_2_state,
                                         bounds=bounds,
                                         args=args,
                                         maxiter=300,
                                         disp=False,
                                         polish=True,
                                         updating="deferred",
                                         workers=1)
        p = results["x"]

        p_rec = np.append(p_rec, [p], axis=0)

        f_name_p = "../fits/fit_2state_ee" + str(sub[i])
        with open(f_name_p, "w") as f:
            np.savetxt(f, p, "%0.4f", ",")

    return p_rec


def fit_state_space_with_g_func_2_state_bcee():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)
    sub = d['sub'].unique()
    length_names = sub.shape[0]

    p_rec = np.empty((0, 6))
    for i in range(length_names):
        print(sub[i])
        x_obs = d[d['sub'] == sub[i]]
        rot = x_obs["Appl_Perturb"].values
        x_obs = x_obs[["bcee", "Target", "target_deg", "trial"]]
        x_obs = x_obs.pivot(index="trial",
                            columns="target_deg",
                            values="bcee")

        x_obs = x_obs.values

        args = (x_obs, rot)
        bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60))
        results = differential_evolution(func=fit_obj_func_sse_2_state,
                                         bounds=bounds,
                                         args=args,
                                         maxiter=300,
                                         disp=False,
                                         polish=True,
                                         updating="deferred",
                                         workers=1)
        p = results["x"]

        p_rec = np.append(p_rec, [p], axis=0)

        f_name_p = "../fits/fit_2state_bcee" + str(sub[i])
        with open(f_name_p, "w") as f:
            np.savetxt(f, p, "%0.4f", ",")

    return p_rec


# NOTE: Inspect and plot the fits
def inspect_fits_individual():
    f_name = '../fit_input/master_data.csv'
    d = pd.read_csv(f_name)
    sub = d['sub'].unique()
    length_names = sub.shape[0]

    for i in range(length_names):
        x_obs = d[d['sub'] == sub[i]]
        rot = x_obs["Appl_Perturb"].values
        x_obs = x_obs[["bcee", "Target", "target_deg", "trial"]]
        x_obs = x_obs.pivot(index="trial",
                            columns="target_deg",
                            values="bcee")
        x_obs = x_obs.values

        params = np.loadtxt("../fits/fit_2state_bcee" + str(sub[i]))
        x_pred = simulate_state_space_with_g_func_2_state(params, rot)

        for ii in range(12):
            plt.scatter(np.arange(0, x_obs.shape[0]),
                        x_obs[:, ii],
                        alpha=0.2)
            plt.plot(x_pred)

        plt.savefig('../figures/fits_individual_subjects/fit_combined_' +
                    str(sub[i]) + ".png")
        plt.close()


def inspect_fits_grp():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)

    theta_values = np.linspace(0.0, 330.0, 12) - 150.0
    rot = d[d['sub'] == 1]['Appl_Perturb'].values
    x_obs_all = d.groupby(['cnd', 'trial', 'Target', 'rot_dir']).mean()
    x_obs_all.reset_index(inplace=True)

    for i in range(x_obs_all['cnd'].unique().shape[0]):
        for j in ('cw', 'ccw'):
            x_obs = x_obs_all[x_obs_all['cnd'] == i]
            x_obs = x_obs[x_obs['rot_dir'] == j]
            x_obs = x_obs[["bcee", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="bcee")

            x_obs = x_obs.values

            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

            params = np.loadtxt("../fits/fit_2state_grp_bcee" + str(i) + '_' + j)
            x_pred = simulate_state_space_with_g_func_2_state(params, rot)

            for ii in range(12):
                ax[0].scatter(np.arange(0, x_obs.shape[0]),
                              x_obs[:, ii],
                              alpha=0.05)
                ax[0].plot(x_pred)

            xg_obs = np.nanmean(x_obs[1000:1072, :], 0)
            xg_pred = np.mean(x_pred[1000:1072, :], 0)
            ax[1].plot(xg_obs, '-')
            ax[1].plot(xg_pred, '-')
            ax[1].set_xticks(ticks=np.arange(0, 12, 1))
            ax[1].set_xticklabels(labels=theta_values)

            plt.savefig('../figures/fits_groups_rot-cnd/fit_combined_grp_' +
                        str(i) + '_' + str(j) + ".png")
            plt.close()


def inspect_fits_grp_overview():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)

    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]
    colors_obs = ['black'] * 12
    colors_obs[5] = colors[0]

    colors_pred = ['black'] * 12
    colors_pred[5] = colors[1]

    a = 0.05 * np.ones(12)
    a[5] = 1.0

    fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12, 8))

    for i in d['cnd'].unique():
        for j in range(len(d['rot_dir'].unique())):
            rd = d['rot_dir'].unique()[j]
            dd = d[(d['cnd'] == i) & (d['rot_dir'] == rd)]
            rot = dd['Appl_Perturb'].values[0:1272]
            x_obs = dd.groupby(['cnd', 'rot_dir', 'Target', 'trial']).mean()
            x_obs.reset_index(inplace=True)
            x_obs = x_obs[["bcee", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="bcee")

            x_obs = x_obs.values

            params = np.loadtxt('../fits/fit_grp_2state_bootstrap_bcee' + str(i) + '_' + rd, delimiter=',')
            params = params.mean(axis=0)
            x_pred = simulate_state_space_with_g_func_2_state(params, rot)

            num_trials = rot.shape[0]
            theta_values = np.linspace(0.0, 330.0, 12) - 150.0
            theta_train_ind = np.where(theta_values == 0.0)[0][0]
            theta_ind = theta_train_ind * np.ones(num_trials, dtype=np.int8)

            title = str('cnd = ') + str(i) + ', rot_dir = ' + rd
            for k in range(12):
                ax[j, i].plot(x_obs[:, k], '-', alpha=a[k], c=colors_obs[k])
                ax[j, i].plot(x_pred[:, k], '-', alpha=a[k], c=colors_pred[k])
                ax[j, i].set_ylim([-10, 30])
                ax[j, i].set_xlim([550, 1272])
                ax[j, i].set_title(title)

            xg_obs = np.nanmean(x_obs[1000:1072, :], 0)
            xg_pred = np.nanmean(x_pred[1000:1072, :], 0)
            ax[j + 2, i].plot(theta_values, xg_obs, '-', c=colors_obs[5])
            ax[j + 2, i].plot(theta_values, xg_pred, '-', c=colors_pred[5])
            ax[j + 2, i].plot(theta_values, xg_obs, '.', c=colors_obs[5])
            ax[j + 2, i].plot(theta_values, xg_pred, '.', c=colors_pred[5])
            ax[j + 2, i].set_title(title)
            ax[j + 2, i].set_ylim([-10, 20])

    plt.suptitle('Group Fits bcEE', fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/bootstrap_results/group_fits_overview_bcee' + ".png")

    fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12, 8))

    for i in d['cnd'].unique():
        for j in range(len(d['rot_dir'].unique())):
            rd = d['rot_dir'].unique()[j]
            dd = d[(d['cnd'] == i) & (d['rot_dir'] == rd)]
            rot = dd['Appl_Perturb'].values[0:1272]
            x_obs = dd.groupby(['cnd', 'rot_dir', 'Target', 'trial']).mean()
            x_obs.reset_index(inplace=True)
            x_obs = x_obs[["Endpoint_Error", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="Endpoint_Error")

            x_obs = x_obs.values

            params = np.loadtxt('../fits/fit_grp_2state_bootstrap_ee' + str(i) + '_' + rd, delimiter=',')
            params = params.mean(axis=0)
            x_pred = simulate_state_space_with_g_func_2_state(params, rot)

            num_trials = rot.shape[0]
            theta_values = np.linspace(0.0, 330.0, 12) - 150.0
            theta_train_ind = np.where(theta_values == 0.0)[0][0]
            theta_ind = theta_train_ind * np.ones(num_trials, dtype=np.int8)

            title = str('cnd = ') + str(i) + ', rot_dir = ' + rd
            for k in range(12):
                ax[j, i].plot(x_obs[:, k], '-', alpha=a[k], c=colors_obs[k])
                ax[j, i].plot(x_pred[:, k], '-', alpha=a[k], c=colors_pred[k])
                ax[j, i].set_ylim([-10, 30])
                ax[j, i].set_xlim([550, 1272])
                ax[j, i].set_title(title)

            xg_obs = np.nanmean(x_obs[1000:1072, :], 0)
            xg_pred = np.nanmean(x_pred[1000:1072, :], 0)
            ax[j + 2, i].plot(theta_values, xg_obs, '-', c=colors_obs[5])
            ax[j + 2, i].plot(theta_values, xg_pred, '-', c=colors_pred[5])
            ax[j + 2, i].plot(theta_values, xg_obs, '.', c=colors_obs[5])
            ax[j + 2, i].plot(theta_values, xg_pred, '.', c=colors_pred[5])
            ax[j + 2, i].set_title(title)
            ax[j + 2, i].set_ylim([-10, 20])

    plt.suptitle('Group Fits EE', fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/bootstrap_results/group_fits_overview_ee' + ".png")


# NOTE: Compute stats, test for significant parameter differences
def bootstrap_t(x_obs, y_obs, x_samp_dist, y_samp_dist, n):
    d_obs = x_obs - y_obs

    d_boot = np.zeros(n)
    xs = np.random.choice(x_samp_dist, n, replace=True)
    ys = np.random.choice(y_samp_dist, n, replace=True)
    d_boot = xs - ys
    d_boot = d_boot - d_boot.mean()

    p_null = (1 + np.sum(np.abs(d_boot) > np.abs(d_obs))) / (n + 1)
    return (p_null)


def inspect_boot_bcee():

    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]
    patterns = ["", "\\\\"]

    r_squared_rec = np.empty(0)

    f_name = '../fit_input/master_data.csv'

    data = pd.read_csv(f_name)
    rot = data[data['sub'] == 1]['Appl_Perturb'].values

    x_obs_all = data.groupby(['cnd', 'Target', 'trial', 'rot_dir']).mean()
    x_obs_all.reset_index(inplace=True)

    dd = []
    cnd = [0, 1, 2]
    rot_dir = ['cw', 'ccw']
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
    sys.stdout = open("../stats/bcee_r_squared_means_ci.txt", "w")
    for j in range(len(rot_dir)):
        for i in cnd:
            d = pd.read_csv('../fits/fit_grp_2state_bootstrap_bcee' + str(i) +
                            '_' + rot_dir[j],
                            header=None)
            d.columns = [
                'alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f'
            ]
            d['rot_dir'] = rot_dir[j]
            d['cnd'] = i
            dd.append(d)

            alpha_s = d['alpha_s'].values
            beta_s = d['beta_s'].values
            sigma_s = d['sigma_s'].values
            alpha_f = d['alpha_f'].values
            beta_f = d['beta_f'].values
            sigma_f = d['sigma_f'].values

            # NOTE: Compute Goodness of fit heatmap - r squared
            x_obs = x_obs_all[x_obs_all['cnd'] == i]
            x_obs = x_obs[x_obs['rot_dir'] == rot_dir[j]]
            x_obs = x_obs[["bcee", "Target", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="bcee")

            x_obs = x_obs.values

            for b in range(len(alpha_s)):
                params = [alpha_s[b], beta_s[b], sigma_s[b],
                          alpha_f[b], beta_f[b], sigma_f[b]]
                result = fit_obj_func_sse_2_state(params, x_obs, rot)
                r_squared = result[1]
                r_squared_rec = np.append(r_squared_rec, r_squared)

            f_name_p = "../fits/r_squared_bootstrap_fits_bcee" + str(i) + '_' + str(rot_dir[j])
            with open(f_name_p, "w") as f:
                np.savetxt(f, r_squared_rec, "%0.4f", ",")

            # NOTE: Save R2 means and CIs in a file
            alpha = 0.05

            ci = np.reshape(
                np.percentile(r_squared_rec, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))

            print(j, i, r_squared_rec.mean(), ci)

            # NOTE: Create heatmap plots
            r_squared_rec = ((np.asarray(r_squared_rec)).reshape(40, 25))

            ax[j, i].set_xticks([])
            ax[j, i].set_yticks([])
            ax[j, i].axis('off')
            ax[j, i].set_title('cnd ' + str(i) + ' ' + str(rot_dir[j]) + ' rotation')

            sbn.heatmap(r_squared_rec, vmin=0, vmax=0.5, ax=ax[j, i])

            r_squared_rec = np.empty(0)

    plt.suptitle('Goodness of Bootstrap Fits - R squared values', fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/bootstrap_results/heatmap_bootstrap_fits' + '.png')
    plt.close()
    sys.stdout.close()

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
    sys.stdout = open("../stats/bcee_params_means_ci.txt", "w")  # to print output in file
    print('rot', 'cnd',
          'alpha_s_mean', 'alpha_s_ci', 'alpha_f_mean', 'alpha_f_ci',
          'beta_s_mean', 'beta_s_ci', 'beta_f_mean', 'beta_f_ci',
          'sigma_s_mean', 'sigma_s_ci', 'sigma_f_mean', 'sigma_f_ci')
    for j in range(len(rot_dir)):
        for i in cnd:
            d = pd.read_csv('../fits/fit_grp_2state_bootstrap_bcee' + str(i) +
                            '_' + rot_dir[j],
                            header=None)
            d.columns = [
                'alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f'
            ]
            d['rot_dir'] = rot_dir[j]
            d['cnd'] = i
            dd.append(d)

            alpha_s = d['alpha_s'].values
            beta_s = d['beta_s'].values
            sigma_s = d['sigma_s'].values
            alpha_f = d['alpha_f'].values
            beta_f = d['beta_f'].values
            sigma_f = d['sigma_f'].values

            alpha = 0.05

            ci1 = np.reshape(
                np.percentile(alpha_s, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            cw_legend = mpatches.Patch(facecolor='#FFFFFF', hatch='', edgecolor='k', label='CW')
            ccw_legend = mpatches.Patch(facecolor='#FFFFFF', hatch='\\\\', edgecolor='k', label='CCW')
            ax[0, 0].legend(handles=[cw_legend, ccw_legend], loc="upper right")
            ax[0, 0].bar(0 + i * len(cnd) + j,
                         alpha_s.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[0, 0].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci1, 'k')
            ax[0, 0].set_title('alpha_s')

            ci2 = np.reshape(
                np.percentile(beta_s, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[0, 1].bar(0 + i * len(cnd) + j,
                         beta_s.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[0, 1].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci2, 'k')
            ax[0, 1].set_title('beta_s')

            ci3 = np.reshape(
                np.percentile(sigma_s, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[0, 2].bar(0 + i * len(cnd) + j,
                         sigma_s.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[0, 2].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci3, 'k')
            ax[0, 2].set_title('sigma_s')

            ci4 = np.reshape(
                np.percentile(alpha_f, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[1, 0].bar(0 + i * len(cnd) + j,
                         alpha_f.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[1, 0].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci4, 'k')
            ax[1, 0].set_title('alpha_f')

            ci5 = np.reshape(
                np.percentile(beta_f, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[1, 1].bar(0 + i * len(cnd) + j,
                         beta_f.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[1, 1].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci5, 'k')
            ax[1, 1].set_title('beta_f')

            ci6 = np.reshape(
                np.percentile(sigma_f, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[1, 2].bar(0 + i * len(cnd) + j,
                         sigma_f.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[1, 2].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci6, 'k')
            ax[1, 2].set_title('sigma_f')

            print(j, i, alpha_s.mean(), ci1, alpha_f.mean(), ci4,
                  beta_s.mean(), ci2, beta_f.mean(), ci5,
                  sigma_s.mean(), ci3, sigma_f.mean(), ci6)

    plt.suptitle('Baseline Corrected Estimates', fontweight='bold')
    plt.setp(ax, xticks=[0.5, 3.5, 6.5], xticklabels=['0°\u03C3', '4°\u03C3', '12°\u03C3'])
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/bootstrap_results/fit_estimates_bcee' + ".png")
    plt.close()
    sys.stdout.close()

    # NOTE: Report stats
    sys.stdout = open("../stats/p_values_bcee.txt", "w")  # to print output in file
    dd = pd.concat(dd)
    for j in ['alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f']:
        for i in dd['cnd'].unique():
            xx = dd[(dd['cnd'] == i) & (dd['rot_dir'] == 'cw')][j].values
            yy = dd[(dd['cnd'] == i) & (dd['rot_dir'] == 'ccw')][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, x)

    print('\n')
    print('\n')

    for j in ['alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f']:
        print('\n')
        for i in dd['rot_dir'].unique():
            print('\n')
            xx = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '00', x)

            xx = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '11', x)

            xx = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '22', x)

            xx = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '01', x)

            xx = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '02', x)

            xx = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '12', x)

        # x = bootstrap_t(xx, yy, 1000)
        # x = stats.ks_2samp(xx, yy)
        # x = x[1]

    sys.stdout.close()


def inspect_boot_ee():

    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]
    patterns = ["", "\\\\"]

    dd = []
    cnd = [0, 1, 2]
    rot_dir = ['cw', 'ccw']
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
    for j in range(len(rot_dir)):
        for i in cnd:
            d = pd.read_csv('../fits/fit_grp_2state_bootstrap_ee' + str(i) +
                            '_' + rot_dir[j],
                            header=None)
            d.columns = [
                'alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f'
            ]
            d['rot_dir'] = rot_dir[j]
            d['cnd'] = i
            dd.append(d)

            alpha_s = d['alpha_s'].values
            beta_s = d['beta_s'].values
            sigma_s = d['sigma_s'].values
            alpha_f = d['alpha_f'].values
            beta_f = d['beta_f'].values
            sigma_f = d['sigma_f'].values

            alpha = 0.05

            ci = np.reshape(
                np.percentile(alpha_s, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            cw_legend = mpatches.Patch(facecolor='#FFFFFF', hatch='', edgecolor='k', label='CW')
            ccw_legend = mpatches.Patch(facecolor='#FFFFFF', hatch='\\\\', edgecolor='k', label='CCW')
            ax[0, 0].legend(handles=[cw_legend, ccw_legend], loc="upper right")
            ax[0, 0].bar(0 + i * len(cnd) + j,
                         alpha_s.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[0, 0].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci, 'k')
            ax[0, 0].set_title('alpha_s')

            ci = np.reshape(
                np.percentile(beta_s, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[0, 1].bar(0 + i * len(cnd) + j,
                         beta_s.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[0, 1].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci, 'k')
            ax[0, 1].set_title('beta_s')

            ci = np.reshape(
                np.percentile(sigma_s, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[0, 2].bar(0 + i * len(cnd) + j,
                         sigma_s.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[0, 2].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci, 'k')
            ax[0, 2].set_title('sigma_s')

            ci = np.reshape(
                np.percentile(alpha_f, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[1, 0].bar(0 + i * len(cnd) + j,
                         alpha_f.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[1, 0].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci, 'k')
            ax[1, 0].set_title('alpha_f')

            ci = np.reshape(
                np.percentile(beta_f, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[1, 1].bar(0 + i * len(cnd) + j,
                         beta_f.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[1, 1].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci, 'k')
            ax[1, 1].set_title('beta_f')

            ci = np.reshape(
                np.percentile(sigma_f, [100 * (alpha / 2), 100 * (1.0 - alpha / 2)]),
                (2, 1))
            ax[1, 2].bar(0 + i * len(cnd) + j,
                         sigma_f.mean(),
                         color=colors[i],
                         hatch=patterns[j])
            ax[1, 2].plot([0 + i * len(cnd) + j, 0 + i * len(cnd) + j], ci, 'k')
            ax[1, 2].set_title('sigma_f')

    plt.suptitle('Not Baseline Corrected Estimates', fontweight='bold')
    plt.setp(ax, xticks=[0.5, 3.5, 6.5], xticklabels=['0°\u03C3', '4°\u03C3', '12°\u03C3'])
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/bootstrap_results/fit_estimates_ee' + ".png")
    plt.close()

    # NOTE: Report stats
    sys.stdout = open("../stats/p_values_ee.txt", "w")  # to print output in file
    dd = pd.concat(dd)
    for j in ['alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f']:
        for i in dd['cnd'].unique():
            xx = dd[(dd['cnd'] == i) & (dd['rot_dir'] == 'cw')][j].values
            yy = dd[(dd['cnd'] == i) & (dd['rot_dir'] == 'ccw')][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, x)

    print('\n')
    print('\n')

    for j in ['alpha_s', 'beta_s', 'sigma_s', 'alpha_f', 'beta_f', 'sigma_f']:
        print('\n')
        for i in dd['rot_dir'].unique():
            print('\n')
            xx = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '00', x)

            xx = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '11', x)

            xx = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '22', x)

            xx = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '01', x)

            xx = dd[(dd['cnd'] == 0) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '02', x)

            xx = dd[(dd['cnd'] == 1) & (dd['rot_dir'] == i)][j].values
            yy = dd[(dd['cnd'] == 2) & (dd['rot_dir'] == i)][j].values
            x = bootstrap_t(xx.mean(), yy.mean(), xx, yy, 1000)
            print(j, i, '12', x)

        # x = bootstrap_t(xx, yy, 1000)
        # x = stats.ks_2samp(xx, yy)
        # x = x[1]

    sys.stdout.close()


# NOTE: Analsying Usedependency as possible explenation for rot_dir effect
def simulate_two_state_space_with_g_func_ud(p, rot):
    alpha_s = p[0]
    beta_s = p[1]
    g_sigma_s = p[2]
    alpha_f = p[3]
    beta_f = p[4]
    g_sigma_f = p[5]
    alpha_ud = p[6]
    beta_ud = p[7]
    sigma_ud = p[8]

    num_trials = rot.shape[0]

    theta_values = np.linspace(0.0, 330.0, 12) - 150.0
    theta_train_ind = np.where(theta_values == 0.0)[0][0]
    theta_ind = theta_train_ind * np.ones(num_trials, dtype=np.int8)

    delta = np.zeros(num_trials)
    x = np.zeros((12, num_trials))
    xs = np.zeros((12, num_trials))
    xf = np.zeros((12, num_trials))
    xud = np.zeros((12, num_trials))
    for i in range(0, num_trials - 1):
        if i > 1000 and i <= 1072:
            delta[i] = 0.0
        elif i > 1172:
            delta[i] = 0.0
        else:
            delta[i] = x[theta_ind[i], i] - rot[i]

        def g_func(theta, theta_mu, sigma):
            if sigma != 0:
                G = np.exp(-(theta - theta_mu)**2 / (2 * sigma**2))
            else:
                G = np.zeros(12)
            return G

        def ud_func(theta, theta_mu, sigma):
            if sigma != 0:
                G = -2 * (logistic.cdf(theta, theta_mu, sigma) - 0.5)
                # G = np.exp(-(theta - theta_mu)**2 / (2 * sigma**2))
                # G[np.where(theta > theta_mu)] = -G[np.where(theta > theta_mu)]
                # G[np.where(theta == theta_mu)] = 0
            else:
                G = np.zeros(12)
            return G

        G_s = g_func(theta_values, theta_values[theta_ind[i]], g_sigma_s)
        G_f = g_func(theta_values, theta_values[theta_ind[i]], g_sigma_f)
        G_ud = ud_func(theta_values, x[theta_ind[i], i], sigma_ud)

        if np.isnan(rot[i]):
            xs[:, i + 1] = beta_s * xs[:, i]
            xf[:, i + 1] = beta_f * xf[:, i]
            xud[:, i + 1] = beta_ud * xud[:, i]
        else:
            xs[:, i + 1] = beta_s * xs[:, i] - alpha_s * delta[i] * G_s
            xf[:, i + 1] = beta_f * xf[:, i] - alpha_f * delta[i] * G_f
            xud[:, i + 1] = beta_ud * xud[:, i] + alpha_ud * G_ud

        x[:, i + 1] = xs[:, i + 1] + xf[:, i + 1] + xud[:, i + 1]

    return x.T


def fit_obj_func_sse_2_state_ud(params, *args):
    x_obs = args[0]
    rot = args[1]

    # hack to impose parameter constraints
    if len(params) > 3:
        alpha_s = params[0]
        beta_s = params[1]
        g_sigma_s = params[2]
        alpha_f = params[3]
        beta_f = params[4]
        g_sigma_f = params[5]
        if (alpha_s >= alpha_f) or (beta_s <= beta_f):
            sse = 10**100
            return sse

    x_pred = simulate_two_state_space_with_g_func_ud(params, rot)

    sse_rec = np.zeros(12)
    for i in range(12):
        # pick trial indices for cost function
        fit_inds = np.concatenate(
            (np.arange(600, 700, 1), np.arange(1000, 1072,
                                               1), np.arange(1172, 1272, 1)))
        # fit_inds = np.arange(0, 1272, 1)
        sse_rec[i] = (np.nansum((x_obs[fit_inds, i] - x_pred[fit_inds, i])**2))
        sse = np.nansum(sse_rec)

    return sse


def fit_state_space_with_g_func_2_state_grp_bcee_ud():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)
    p_rec = np.empty((0, 9))

    rot = d[d['sub'] == 1]['Appl_Perturb'].values

    x_obs_all = d.groupby(['cnd', 'Target', 'trial', 'rot_dir']).mean()
    x_obs_all.reset_index(inplace=True)

    for i in range(x_obs_all['cnd'].unique().shape[0]):
        for j in ('cw', 'ccw'):
            x_obs = x_obs_all[x_obs_all['cnd'] == i]
            x_obs = x_obs[x_obs['rot_dir'] == j]
            x_obs = x_obs[["bcee", "Target", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="bcee")

            x_obs = x_obs.values

            args = (x_obs, rot)
            bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60),
                      (0, 1), (0, 1), (0, 60))
            results = differential_evolution(func=fit_obj_func_sse_2_state_ud,
                                             bounds=bounds,
                                             args=args,
                                             maxiter=300,
                                             disp=False,
                                             polish=True,
                                             updating="deferred",
                                             workers=1)
            p = results["x"]

            p_rec = np.append(p_rec, [p], axis=0)

            f_name_p = "../fits/fit_grp_2state_bcee_ud" + str(i) + '_' + j
            with open(f_name_p, "w") as f:
                np.savetxt(f, p, "%0.4f", ",")

    return p_rec


def fit_state_space_with_g_func_2_state_grp_boot_bcee_ud():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)

    n_boot_samp = 1000

    for i in d['cnd'].unique():
        for j in d['rot_dir'].unique():
            p_rec = -1 * np.ones((n_boot_samp, 9))
            for b in range(n_boot_samp):
                print(i, j, b)
                dd = d[(d['cnd'] == i) & (d['rot_dir'] == j)]
                rot = dd['Appl_Perturb'].values[0:1272]
                subs = dd['sub'].unique()
                boot_subs = np.random.choice(subs,
                                             size=subs.shape[0],
                                             replace=True)
                x_boot_rec = []
                for k in boot_subs:
                    x_boot_rec.append(d[d['sub'] == k])
                    x_boot = pd.concat(x_boot_rec)

                x_obs = x_boot.groupby(['cnd', 'rot_dir', 'Target',
                                        'trial']).mean()
                x_obs.reset_index(inplace=True)

                x_obs = x_obs[["bcee", "target_deg", "trial"]]
                x_obs = x_obs.pivot(index="trial",
                                    columns="target_deg",
                                    values="bcee")

                x_obs = x_obs.values

                args = (x_obs, rot)
                bounds = ((0, 1), (0, 1), (0, 60), (0, 1), (0, 1), (0, 60),
                          (0, 1), (0, 1), (0, 60))
                results = differential_evolution(func=fit_obj_func_sse_2_state_ud,
                                                 bounds=bounds,
                                                 args=args,
                                                 maxiter=300,
                                                 disp=False,
                                                 polish=False,
                                                 updating="deferred",
                                                 workers=1)
                p = results["x"]
                p_rec[b, :] = p

                f_name_p = "../fits/fit_grp_2state_bootstrap_bcee_ud" + str(
                    i) + '_' + j
                with open(f_name_p, "w") as f:
                    np.savetxt(f, p_rec, "%0.4f", ",")


def inspect_fits_grp_overview_ud():
    f_name = '../fit_input/master_data.csv'

    d = pd.read_csv(f_name)

    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]
    colors_obs = ['black'] * 12
    colors_obs[5] = colors[0]

    colors_pred = ['black'] * 12
    colors_pred[5] = colors[1]

    a = 0.05 * np.ones(12)
    a[5] = 1.0

    fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12, 8))

    for i in d['cnd'].unique():
        for j in range(len(d['rot_dir'].unique())):
            rd = d['rot_dir'].unique()[j]
            dd = d[(d['cnd'] == i) & (d['rot_dir'] == rd)]
            rot = dd['Appl_Perturb'].values[0:1272]
            x_obs = dd.groupby(['cnd', 'rot_dir', 'Target', 'trial']).mean()
            x_obs.reset_index(inplace=True)
            x_obs = x_obs[["bcee", "target_deg", "trial"]]
            x_obs = x_obs.pivot(index="trial",
                                columns="target_deg",
                                values="bcee")

            x_obs = x_obs.values

            # params = np.loadtxt('../fits/fit_grp_2state_bootstrap_bcee_ud' +
            #                     str(i) + '_' + rd, delimiter=',')
            # params = params.mean(axis=0) -> REACTIVATE!
            params = np.loadtxt('../fits/fit_grp_2state_bcee_ud' +
                                str(i) + '_' + rd, delimiter=',')  # -> DEACTIVATE!
            x_pred = simulate_two_state_space_with_g_func_ud(params, rot)

            num_trials = rot.shape[0]
            theta_values = np.linspace(0.0, 330.0, 12) - 150.0
            theta_train_ind = np.where(theta_values == 0.0)[0][0]
            theta_ind = theta_train_ind * np.ones(num_trials, dtype=np.int8)

            title = str('cnd = ') + str(i) + ', rot_dir = ' + rd
            for k in range(12):
                ax[j, i].plot(x_obs[:, k], '-', alpha=a[k], c=colors_obs[k])
                ax[j, i].plot(x_pred[:, k], '-', alpha=a[k], c=colors_pred[k])
                ax[j, i].set_ylim([-10, 30])
                ax[j, i].set_xlim([550, 1272])
                ax[j, i].set_title(title)

            xg_obs = np.nanmean(x_obs[1000:1072, :], 0)
            xg_pred = np.nanmean(x_pred[1000:1072, :], 0)
            ax[j + 2, i].plot(theta_values, xg_obs, '-', c=colors_obs[5])
            ax[j + 2, i].plot(theta_values, xg_pred, '-', c=colors_pred[5])
            ax[j + 2, i].plot(theta_values, xg_obs, '.', c=colors_obs[5])
            ax[j + 2, i].plot(theta_values, xg_pred, '.', c=colors_pred[5])
            ax[j + 2, i].set_title(title)
            ax[j + 2, i].set_ylim([-10, 25])

    plt.suptitle('Group Fits bcEE Use-Dependence Model', fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../figures/bootstrap_results/group_fits_overview_ud' + ".png")


# NOTE: Execute functions
# NOTE: Change return function in fit_obj_func_sse_2_state
#       before running the bootstrap functions
# fit_state_space_with_g_func_2_state_grp_bootstrap_bcee()
# fit_state_space_with_g_func_2_state_grp_bcee()
# fit_state_space_with_g_func_2_state_bcee()

# fit_state_space_with_g_func_2_state_grp_bootstrap_ee()
# fit_state_space_with_g_func_2_state_grp_ee()
# fit_state_space_with_g_func_2_state_ee()

# inspect_boot_bcee()
# inspect_boot_ee()

# inspect_fits_individual()
# inspect_fits_grp()
# inspect_fits_grp_overview()

# fit_state_space_with_g_func_2_state_grp_bcee_ud()
# fit_state_space_with_g_func_2_state_grp_boot_bcee_ud()
# inspect_fits_grp_overview_ud()

plot_hypothesis()
