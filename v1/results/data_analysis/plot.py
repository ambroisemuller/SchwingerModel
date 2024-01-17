import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys;

if __name__=="__main__":

    args = sys.argv

    T = int(args[1])
    L = int(args[2])
    kappa = float(args[3])
    beta = float(args[4])

    plot_dH = ('dH' in args) or ('all' in args)
    plot_acc = ('acc' in args) or ('all' in args)
    plot_plaq = ('plaq' in args) or ('all' in args)
    plot_qtop = ('qtop' in args) or ('all' in args)
    plot_corr = ('corr' in args) or ('all' in args)

    ntherm = 0
    if 'all' in args: 
        ntherm = int(args[-1])
    
    if plot_dH:
        df_dH = pd.read_csv("../observables/dH.csv")
        data_dH = {col: df_dH[col].values for col in df_dH.columns}
    if plot_acc:
        df_acc = pd.read_csv("../observables/acc.csv")
        data_acc = {col: df_acc[col].values for col in df_acc.columns}
    if plot_plaq:
        df_plaq = pd.read_csv("../observables/plaq.csv")
        data_plaq = {col: df_plaq[col].values for col in df_plaq.columns}
    if plot_qtop:
        df_qtop = pd.read_csv("../observables/qtop.csv")
        data_qtop = {col: df_qtop[col].values for col in df_qtop.columns}

    if plot_dH:
        t = data_dH['time']
        dH = data_dH['dH']
    if plot_acc:
        t = data_acc['time']
        acc = data_acc['acc']
    if plot_plaq:
        t = data_plaq['time']
        plaq = data_plaq['plaq']
    if plot_qtop:
        t = data_qtop['time']
        qtop = data_qtop['qtop']

    fig, ((ax_dH, ax_acc), (ax_plaq, ax_qtop)) = plt.subplots(2, 2, figsize=(12, 6))
    fig.suptitle(r"$T = $%d, $L = $%d, $\kappa = $%f, $\beta = $%f" % (T, L, kappa, beta))
    if plot_dH:
        ax_dH.plot(t, dH, '.-', label='HMC measurements')
        ax_dH.plot(t[ntherm:], [np.mean(dH[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(dH[ntherm:]), np.std(dH[ntherm:])))
        ax_dH.set_xlabel(r"Simulation time $\tau$")
        ax_dH.set_ylabel(r"Energy change $dH$")
        ax_dH.legend()
    if plot_acc:
        ax_acc.plot(t, acc, '.-', label='HMC measurements')
        ax_acc.plot(t[ntherm:], [np.mean(acc[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(acc[ntherm:]), np.std(acc[ntherm:])))
        ax_acc.set_xlabel(r"Simulation time $\tau$")
        ax_acc.set_ylabel(r"Acceptance")
        ax_acc.legend()
    if plot_plaq:
        ax_plaq.plot(t, plaq, '.-', label='HMC measurements')
        ax_plaq.plot(t[ntherm:], [np.mean(plaq[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(plaq[ntherm:]), np.std(plaq[ntherm:])))
        ax_plaq.set_xlabel(r"Simulation time $\tau$")
        ax_plaq.set_ylabel(r"Plaquette")
        ax_plaq.legend()
    if plot_qtop:
        ax_qtop.plot(t, qtop, '.-', label='HMC measurements')
        ax_qtop.plot(t[ntherm:], [np.mean(qtop[ntherm:]) for _ in t[ntherm:]], '--', label="avg %f\nstd %f"%(np.mean(qtop[ntherm:]), np.std(qtop[ntherm:])))
        ax_qtop.set_xlabel(r"Simulation time $\tau$")
        ax_qtop.set_ylabel(r"Topological charge")
        ax_qtop.legend()

    fig.tight_layout()

    plot_file = "../plots/general.png"
    fig.savefig(plot_file)
    print("General observable plots saved to "+plot_file)

    if plot_corr:

        df_PP = pd.read_csv("../observables/PP.csv")
        time = df_PP.iloc[:, 0].values
        number_of_time_values = len(time)
        df_PP = df_PP.drop(df_PP.columns[0], axis=1)
        if df_PP.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_PP.values.reshape(number_of_time_values, T, L)
        average_values = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        fig_PP = plt.figure(figsize=(12,5))
        fig_PP.suptitle(r"$T = $%d, $L = $%d, $\kappa = $%f, $\beta = $%f" % (T, L, kappa, beta))
        ax = fig_PP.add_subplot(122)
        heatmap = ax.imshow(np.log(average_values), cmap='viridis', aspect='equal')
        contours = ax.contour(np.log(average_values), colors='black', linewidths=0.5)
        ax.clabel(contours, inline=True, fontsize=8)
        ax.set_xlabel('x')
        ax.set_ylabel('t')
        ax.set_title(r'Pseudoscalar current correlator $\langle P(t, x) P(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax, fraction=0.026, pad=0.04)
        cbar.set_label(r'$\log \langle P(t, x) P(0, 0) \rangle$')
        ax1 = fig_PP.add_subplot(121)
        ax1.errorbar(np.arange(L), average_summed_values, std_summed_values)
        ax1.set_yscale('log')
        ax1.set_xlabel('x')
        ax1.set_ylabel(r'$\log \sum_t \langle P(t, x) P(0, 0) \rangle$')
        ax1.set_title("Summed pseudoscalar current correlator")
        fig_PP.tight_layout()
        fig_PP.savefig("../plots/corr_PP.png")

        pp_values = average_values[1:-1, 1:-1]


        fig_VV, ((ax00, ax00s, ax01, ax01s), (ax10, ax10s, ax11, ax11s)) = plt.subplots(2, 4, figsize=(18, 12*(T/L)-1))

        df_00 = pd.read_csv("../observables/V0V0.csv")
        time = df_00.iloc[:, 0].values
        number_of_time_values = len(time)
        df_00 = df_00.drop(df_00.columns[0], axis=1)
        if df_00.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_00.values.reshape(number_of_time_values, T, L)
        average_values = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        heatmap = ax00.imshow(np.log(average_values), cmap='viridis', aspect='equal')
        contours = ax00.contour(np.log(average_values), colors='black', linewidths=0.5)
        ax00.clabel(contours, inline=True, fontsize=8)
        ax00.set_xlabel('x')
        ax00.set_ylabel('t')
        ax00.set_title(r'Vector current correlator $\langle V_0(t, x) V_0(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax00, fraction=0.026, pad=0.04)
        cbar.set_label(r'$\log \langle V_0(t, x) V_0(0, 0) \rangle$')
        ax00s.errorbar(np.arange(L), average_summed_values, std_summed_values)
        ax00s.set_yscale('log')
        ax00s.set_xlabel('x')

        df_01 = pd.read_csv("../observables/V0V1.csv")
        time = df_01.iloc[:, 0].values
        number_of_time_values = len(time)
        df_01 = df_01.drop(df_01.columns[0], axis=1)
        if df_01.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_01.values.reshape(number_of_time_values, T, L)
        average_values = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        heatmap = ax01.imshow(((average_values)), cmap='viridis', aspect='equal')
        contours = ax01.contour(((average_values)), colors='black', linewidths=0.5)
        ax01.clabel(contours, inline=True, fontsize=8)
        ax01.set_xlabel('x')
        ax01.set_ylabel('t')
        ax01.set_title(r'Vector current correlator $\langle V_0(t, x) V_1(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax01, fraction=0.026, pad=0.04)
        cbar.set_label(r'$\langle V_0(t, x) V_1(0, 0) \rangle$')
        ax01s.errorbar(np.arange(L), average_summed_values, std_summed_values)
        # ax01s.set_yscale('log')
        ax01s.set_xlabel('x')

        df_10 = pd.read_csv("../observables/V1V0.csv")
        time = df_10.iloc[:, 0].values
        number_of_time_values = len(time)
        df_10 = df_10.drop(df_10.columns[0], axis=1)
        if df_10.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_10.values.reshape(number_of_time_values, T, L)
        average_values = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        heatmap = ax10.imshow(((average_values)), cmap='viridis', aspect='equal')
        contours = ax10.contour(((average_values)), colors='black', linewidths=0.5)
        ax10.clabel(contours, inline=True, fontsize=8)
        ax10.set_xlabel('x')
        ax10.set_ylabel('t')
        ax10.set_title(r'Vector current correlator $\langle V_1(t, x) V_0(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax10, fraction=0.026, pad=0.04)
        cbar.set_label(r'$ \langle V_1(t, x) V_0(0, 0) \rangle$')
        ax10s.errorbar(np.arange(L), average_summed_values, std_summed_values)
        # ax10s.set_yscale('log')
        ax10s.set_xlabel('x')

        df_11 = pd.read_csv("../observables/V1V1.csv")
        time = df_11.iloc[:, 0].values
        number_of_time_values = len(time)
        df_11 = df_11.drop(df_11.columns[0], axis=1)
        if df_11.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_11.values.reshape(number_of_time_values, T, L)
        average_values = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        heatmap = ax11.imshow(np.log(average_values), cmap='viridis', aspect='equal')
        contours = ax11.contour(np.log(average_values), colors='black', linewidths=0.5)
        ax11.clabel(contours, inline=True, fontsize=8)
        ax11.set_xlabel('x')
        ax11.set_ylabel('t')
        ax11.set_title(r'Vector current correlator $\langle V_1(t, x) V_1(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax11, fraction=0.026, pad=0.04)
        cbar.set_label(r'$\log \langle V_1(t, x) V_1(0, 0) \rangle$')
        ax11s.errorbar(np.arange(L), average_summed_values, std_summed_values)
        ax11s.set_yscale('log')
        ax11s.set_xlabel('x')

        fig_VV.tight_layout()
        fig_VV.savefig("../plots/corr_VV.png")




        fig_AP, ((ax0, ax0s, ax1, ax1s), (ax_div, ax_divs, ax_m, ax_ms)) = plt.subplots(2, 4, figsize=(18, 12*(T/L)-1))

        df_0 = pd.read_csv("../observables/A0P.csv")
        time = df_0.iloc[:, 0].values
        number_of_time_values = len(time)
        df_0 = df_0.drop(df_0.columns[0], axis=1)
        if df_0.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_0.values.reshape(number_of_time_values, T, L)
        average_values0 = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        heatmap = ax0.imshow((((average_values0))), cmap='viridis', aspect='equal')
        contours = ax0.contour((((average_values0))), colors='black', linewidths=0.5)
        ax0.clabel(contours, inline=True, fontsize=8)
        ax0.set_xlabel('x')
        ax0.set_ylabel('t')
        ax0.set_title(r'Correlator $\langle A_0(t, x) P(0, 0) \rangle + \langle P(x, t) A_0(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax0, fraction=0.026, pad=0.04)
        cbar.set_label(r'$\log \langle A_0(t, x) P(0, 0) \rangle + \langle P(x, t) A_0(0, 0) \rangle$')
        ax0s.errorbar(np.arange(L), average_summed_values, std_summed_values)
        # ax0s.set_yscale('log')
        ax0s.set_xlabel('x')

        df_1 = pd.read_csv("../observables/A1P.csv")
        time = df_1.iloc[:, 0].values
        number_of_time_values = len(time)
        df_1 = df_1.drop(df_1.columns[0], axis=1)
        if df_1.shape[1] != T * L:
            raise ValueError("The number of elements in the CSV does not match number_of_time_values * T * L")
        reshaped_array = df_1.values.reshape(number_of_time_values, T, L)
        average_values1 = np.mean(reshaped_array[ntherm:, :, :], axis=0)
        summed_values = np.sum(reshaped_array, axis=1)
        average_summed_values = np.mean(summed_values[ntherm:, :], axis=0)
        std_summed_values = np.std(summed_values[ntherm:, :], axis=0)
        heatmap = ax1.imshow((((average_values1))), cmap='viridis', aspect='equal')
        contours = ax1.contour((((average_values1))), colors='black', linewidths=0.5)
        ax1.clabel(contours, inline=True, fontsize=8)
        ax1.set_xlabel('x')
        ax1.set_ylabel('t')
        ax1.set_title(r'Correlator $\langle A_1(t, x) P(0, 0) \rangle + \langle P(x, t) A_1(0, 0) \rangle$')
        cbar = fig.colorbar(heatmap, ax=ax1, fraction=0.026, pad=0.04)
        cbar.set_label(r'$\langle A_1(t, x) P(0, 0) \rangle + \langle P(x, t) A_1(0, 0) \rangle$')
        ax1s.errorbar(np.arange(L), average_summed_values, std_summed_values)
        # ax1s.set_yscale('log')
        ax1s.set_xlabel('x')

        print(average_values0.shape)
        d0_c0 = 0.5*(average_values0[2:, 1:-1] - average_values0[:-2, 1:-1])
        print(average_values1.shape)
        d1_c1 = 0.5*(average_values1[1:-1, 2:] - average_values1[1:-1, :-2])

        div_arr = d0_c0 + d1_c1
        heatmap = ax_div.imshow(div_arr, cmap='viridis', aspect='equal')
        contours = ax_div.contour(div_arr, colors='black', linewidths=0.5)
        ax_div.clabel(contours, inline=True, fontsize=8)
        ax_div.clabel(contours, inline=True, fontsize=8)
        ax_div.set_xlabel('x')
        ax_div.set_ylabel('t')
        ax_div.set_title(r'$\partial_\mu [ \langle A_\mu(t, x) P(0, 0) \rangle + \langle P(x, t) A_\mu(0, 0) \rangle]$')
        cbar = fig.colorbar(heatmap, ax=ax_div, fraction=0.026, pad=0.04)
        cbar.set_label(r'divergence')
        div_sum = np.sum(div_arr, axis=0)
        ax_divs.plot(div_sum, '.-')
        # ax_divs.set_yscale('log')
        ax_divs.set_xlabel('x')

        m_arr = div_arr/pp_values
        heatmap = ax_m.imshow(m_arr, cmap='viridis', aspect='equal')
        contours = ax_m.contour(m_arr, colors='black', linewidths=0.5)
        ax_m.clabel(contours, inline=True, fontsize=8)
        ax_m.clabel(contours, inline=True, fontsize=8)
        ax_m.set_xlabel('x')
        ax_m.set_ylabel('t')
        ax_m.set_title(r' $\frac{\partial_\mu [ \langle A_\mu(t, x) P(0, 0) \rangle + \langle P(x, t) A_\mu(0, 0) \rangle]}{\langle P(0, 0) \rangle}$')
        cbar = fig.colorbar(heatmap, ax=ax_m, fraction=0.026, pad=0.04)
        cbar.set_label(r'PCAC mass')
        m_sum = np.sum(m_arr, axis=0)
        ax_ms.plot(m_sum, '.-')
        # ax_ms.set_yscale('log')
        ax_ms.set_xlabel('x')

        fig_AP.tight_layout()
        fig_AP.savefig("../plots/corr_AP.png")


    print("re-run and exclude thermalization trajectories in averages with\npython plot.py %d %d %f %f all $\{N_THERM\}" % (T, L, kappa, beta))