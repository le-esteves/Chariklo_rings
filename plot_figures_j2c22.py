import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Circle
import os, sys, glob, pickle
from datetime import datetime
import pandas as pd
import numpy as np

def plot_snapshot_w_tp(ax, T, i, df_disk, df_central, size, point_size_factor):
    if i >= 3:
        ax.set_xlabel('x (R$_{Chariklo}$)')
    else:
        ax.xaxis.set_tick_params(labelbottom=False)
    if i % 3 == 0:
        ax.set_ylabel("y (R$_{Chariklo}$)")
    else:
        ax.yaxis.set_tick_params(labelleft=False)
    ax.set_ylim(-size, size)
    ax.set_xlim(-size, size)
    ax.set_xticks(np.arange(-size,size+0.1,2))
    ax.set_yticks(np.arange(-size,size+0.1,2))
    ax.set_aspect(1)

    # Adiciona o número do plot no canto superior direito
    annotation = ax.annotate(
        aux[i], 
        xy=(0.97, 0.97), 
        xycoords='axes fraction',
        fontsize=13, 
        #fontweight='bold',
        ha='right', 
        va='top',
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7)
    )
    #print (df_central)

    # Dados do corpo central
    e_x = float(df_central["x"].iloc[0])
    e_y = float(df_central["y"].iloc[0])
    e_r = float(df_central["radius"].iloc[0])
    e_m = float(df_central["mass"].iloc[0])

    # Colormap
    cmap = plt.get_cmap('gnuplot')
    color_test = cmap(0.0)
    color_central = cmap(1.0)

    # Plotando o corpo central
    ax.add_patch(Circle((e_x, e_y), e_r, color=color_central))
    
    # raio das partículas
    par_radius = 0.0001

    # Plotando as partículas do disco
    scatter = ax.scatter(df_disk["x"], df_disk["y"], s=par_radius * point_size_factor, 
                         color=color_test, label="Test particles")

    # Plotando as partículas dos clusteres
    #unique_labels = df_clusters[df_clusters["time"] == T]['label'].unique()
    # unique_labels = df_clusters['label'].unique()
    # for cluster_idx, cluster_label in enumerate(unique_labels):
    #     cluster_particles = df_clusters_par[df_clusters_par["label"] == cluster_label]
    #     #print (cluster_particles)
    #     cluster_mass = np.sum(cluster_particles["mass"])
    #     cluster_color = cmap(norm(cluster_mass))

    #     # Desenha as partículas do cluster com uma cor diferente
    #     ax.scatter(cluster_particles['x'], cluster_particles['y'], 
    #                s=(np.pi * par_radius ** 2) * point_size_factor, 
    #                color=cluster_color, label=f'Cluster {cluster_label}')

    #     # Desenha um círculo ao redor do cluster (opcional)
    #     if not cluster_particles.empty:
    #         continue
    #         cluster_center_x = cluster_particles['x'].mean()
    #         cluster_center_y = cluster_particles['y'].mean()
    #         cluster_radius = np.sqrt((cluster_particles['x'] - cluster_center_x)**2 + 
    #                                  (cluster_particles['y'] - cluster_center_y)**2).max()
    #         #print (cluster_radius)

    #         circle = Circle((cluster_center_x, cluster_center_y), cluster_radius*2, 
    #                         fill=False, edgecolor='cyan', linewidth=.3, alpha=1)
    #         ax.add_patch(circle)

    #ax.text(0.9 * size, 0.9 * size, f"orbit {T/2/np.pi:.0f}", ha="right")
    time = T*0.59/24 # T x 0.59 = tempo em horas para Rchariklo = 124 km e Mchariklo = 6.3e18 kg
    ax.set_title(f"{time:.0f} days")

    # Plotando o limite do fluido de Roche
    #ax.add_patch(Circle((e_x, e_y), 2.9, edgecolor='green', facecolor='None'))

    # Plotando as ressonancias orbital/rotação 1:1, 2:1, 3:1
    # ax.add_patch(Circle((e_x, e_y), 1.526, edgecolor='green', facecolor='None', lw=0.6))
    # ax.add_patch(Circle((e_x, e_y), 2.422, edgecolor='blue', facecolor='None', lw=0.6))
    # ax.add_patch(Circle((e_x, e_y), 3.174, edgecolor='red', facecolor='None', lw=0.6))
    ax.add_patch(Circle((e_x, e_y), 1.526, edgecolor='green', facecolor='None', lw=0.6))
    ax.add_patch(Circle((e_x, e_y), 2.422, edgecolor='blue', facecolor='None', lw=0.6))
    ax.add_patch(Circle((e_x, e_y), 2.8, edgecolor='red', facecolor='None', lw=0.6))
    return scatter

if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        if arg.startswith("path="):
            path = arg.split("=")[1]

plt.rcParams["axes.labelsize"] = 13
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 13

namefile = path + f'Analysis/plot_xy_snapshots.png'
clusters_file = path + 'Outputs/clusters_data.pkl'
clusters_par_file = path + 'Outputs/clusters_particles_data.pkl'
disk_file = path + 'Outputs/disk_data.pkl'
central_file = path + 'Outputs/central_data.pkl'

df_clusters = pd.read_pickle(clusters_file)
df_clusters_par = pd.read_pickle(clusters_par_file)
df_disk = pd.read_pickle(disk_file)
df_central = pd.read_pickle(central_file)

# Encontrar os valores mínimo e máximo de massa em todos os dados
# all_masses = pd.concat([
#     df_clusters_par.groupby('label')['mass'].sum(), 
#     df_disk['mass'], 
#     df_central['mass']
# ])
# massmin = all_masses.min()
# massmax = all_masses.max()

# Criar a normalização e colormap que será usado por todos os subplots
# norm = LogNorm(vmin=10**int(np.log10(massmin)), vmax=1.05)
# cmap = plt.get_cmap('gnuplot')

Time = df_central["time"].to_list()
T_ = [0, 400, 8130, 20000, 40000, 80000]
aux = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
size = 4
point_size_factor = 300

# Criar figura com espaço para a colorbar
fig = plt.figure(figsize=(9, 6))
# Ajustar o layout para deixar espaço para a colorbar à direita
fig.subplots_adjust(bottom=0.12, hspace=0.13, wspace=0.13)
# Criar gridspec para os subplots
gs = gridspec.GridSpec(nrows=2, ncols=3, figure=fig)

scatters = []  # Para armazenar os scatter plots para a colorbar

for i, T in enumerate(T_):
    ax = fig.add_subplot(gs[i])
    
    # selecionar apenas os dados onde T = T_[i]
    df_central['time_diff'] = (df_central['time'] - T).abs()
    min_diff_central = df_central['time_diff'].min()
    df_centralT = df_central[df_central['time_diff'] == min_diff_central].copy()
    true_time = df_centralT["time"].to_list()[0]

    # Repetir o processo para df_disk
    df_diskT = df_disk[df_disk['time'] == true_time].copy()

    T_dump = 40000
    if (T == T_dump):
        with open(f"{path}/Outputs/snap_{true_time:.0f}.pkl", "wb") as f:
            pickle.dump(df_diskT, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Repetir o processo para df_clusters
    #df_clustersT = df_clusters[df_clusters['time'] == true_time].copy()

    # Repetir o processo para df_clusters_par
    #df_clusters_parT = df_clusters_par[df_clusters_par['time'] == true_time].copy()
    
    print(f"Plotting snapshot: {T:.2f} | {true_time:.2f}")
    #print(f"numeber of particles sats: {len(df_clusters_parT)} | disk: {len(df_diskT)}")

    scatter = plot_snapshot_w_tp(ax, true_time, i, df_diskT, df_centralT, size, point_size_factor)
    scatters.append(scatter)

# Adicionar uma única colorbar horizontal na parte inferior do plot
# cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.03])  # [left, bottom, width, height]
# cbar = fig.colorbar(scatters[0], cax=cbar_ax, orientation='horizontal')
# cbar.set_label('mass (M$_{\\oplus}$)')

# Salva e fecha o gráfico
fig.savefig(namefile, dpi=400, bbox_inches='tight')
plt.close(fig)