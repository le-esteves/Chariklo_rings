import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Circle
import os, sys, glob, pickle
from datetime import datetime
import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN

import matplotlib.pyplot as plt

def analyze_rings_snapshot(df_rings_par, df_rings, snapshot_time,
                                    n_sectors=18, save_prefix="rings_analysis"):
    """
    Analisa anéis em um snapshot: densidade angular/latitudinal físicas e versão polar.
    
    df_rings_par : DataFrame com partículas (x,y,z,...) e 'label'
    df_rings     : DataFrame com estatísticas dos anéis
    snapshot_time: tempo real do snapshot (usado no título/arquivo)
    n_sectors    : número de setores angulares (18=20° cada)
    save_prefix  : prefixo dos arquivos de saída
    """

    # -------------------
    # Constantes físicas
    # -------------------
    Mchariklo = 6.3e18   # kg
    Rchariklo = 124.0    # km
    rho_p = 1.0          # g/cm3 = 1000 kg/m3
    mp = Mchariklo/1000 / 50000  # massa individual das partículas [g]
    
    # converter unidades: assumimos coordenadas em unidades de Rchariklo
    unit_length_cm = Rchariklo * 100000  # 1 unidade = 124 km = 12400000 cm
    unit_area = (unit_length_cm**2)  # cm²

    results = []
    n_rings = df_rings["label"].nunique()
    ring_labels = sorted(df_rings["label"].unique())
    ring_colors_ = ['tab:blue', 'tab:orange', 'tab:purple']

    # Preparar figuras e eixos
    fig, axes = plt.subplots(nrows=n_rings, ncols=4, figsize=(15, 4*n_rings))
    if n_rings == 1:
        axes = np.array([axes])  # garante array 2D

    # Pré-calcular escalas globais
    # 1) Longitude
    sector_edges = np.linspace(0, 360, n_sectors+1)
    # 2) Z
    all_z_vals = df_rings_par["z"] * Rchariklo  # km
    z_min, z_max = -5, 5 #all_z_vals.min(), all_z_vals.max()
    bins_z = np.linspace(z_min, z_max, 81)
    dz = bins_z[1] - bins_z[0]
    # 3) r (km)
    all_r_km = np.sqrt(df_rings_par["x"]**2 + df_rings_par["y"]**2) * Rchariklo
    r_min, r_max = 0, 500 #all_r_km.min(), all_r_km.max()

    for idx, lbl in enumerate(ring_labels):
        sub = df_rings_par[df_rings_par["label"] == lbl]
        if len(sub) == 0:
            continue

        # Coordenadas polares
        r = np.sqrt(sub["x"]**2 + sub["y"]**2) * Rchariklo * unit_length_cm  # cm
        r_km = np.sqrt(sub["x"]**2 + sub["y"]**2) * Rchariklo  # km
        theta = np.arctan2(sub["y"], sub["x"])              # rad
        theta_deg = (np.degrees(theta) + 360) % 360         # [0,360)

        # ---- 1) Densidade angular (massa por área em setores) ----
        sector_mass_density_all = []
        for i in range(n_sectors):
            mask = (theta_deg >= sector_edges[i]) & (theta_deg < sector_edges[i+1])
            r_sector = r[mask]
            n_part = len(r_sector)
            if n_part == 0:
                sector_mass_density_all.append(0)
                continue

            r_mean = r_sector.mean()
            dr = r_sector.max() - r_sector.min() if n_part > 1 else 1e-10
            dtheta_rad = np.radians(sector_edges[i+1]-sector_edges[i])
            area_sector = r_mean * dr * dtheta_rad  # cm² (aprox)

            mass_sector = n_part * mp  # g
            sigma_sector = mass_sector / area_sector if area_sector > 0 else 0  # g/cm²
            sector_mass_density_all.append(sigma_sector)

            # salvar estatísticas
            results.append({
                "snapshot_time": snapshot_time,
                "label": lbl,
                "sector_deg": f"{sector_edges[i]:.0f}-{sector_edges[i+1]:.0f}",
                "n_particles": n_part,
                "r_mean": r_mean,
                "width": dr,
                "sigma": sigma_sector
            })

        # Plotar todos os anéis em cinza claro
        ax1 = axes[idx, 0]
        for j, lbl2 in enumerate(ring_labels):
            sub2 = df_rings_par[df_rings_par["label"] == lbl2]
            r2 = np.sqrt(sub2["x"]**2 + sub2["y"]**2) * Rchariklo * unit_length_cm
            theta2 = np.arctan2(sub2["y"], sub2["x"])
            theta2_deg = (np.degrees(theta2) + 360) % 360
            sector_mass_density2 = []
            for i in range(n_sectors):
                mask2 = (theta2_deg >= sector_edges[i]) & (theta2_deg < sector_edges[i+1])
                r_sector2 = r2[mask2]
                n_part2 = len(r_sector2)
                if n_part2 == 0:
                    sector_mass_density2.append(0)
                    continue
                r_mean2 = r_sector2.mean()
                dr2 = r_sector2.max() - r_sector2.min() if n_part2 > 1 else 1e-10
                dtheta_rad2 = np.radians(sector_edges[i+1]-sector_edges[i])
                area_sector2 = r_mean2 * dr2 * dtheta_rad2
                mass_sector2 = n_part2 * mp
                sigma_sector2 = mass_sector2 / area_sector2 if area_sector2 > 0 else 0
                sector_mass_density2.append(sigma_sector2)
            color = ring_colors_[j] if lbl2 == lbl else "lightgray"
            alpha = 0.7 if lbl2 == lbl else 0.3
            ax1.bar((sector_edges[:-1]+sector_edges[1:])/2, sector_mass_density2,
                    width=360/n_sectors, color=color, alpha=alpha, label=f"Anel {lbl2}" if lbl2 == lbl else None)
        ax1.set_xlabel("Longitude (°)")
        ax1.set_ylabel(r"$\Sigma$ (g/cm$^2$)")
        #ax1.set_yscale("log")
        ax1.set_title(f"Anel {lbl} - densidade angular")
        ax1.set_xlim(0, 360)
        #ax1.legend()

        # ---- 2) Densidade vertical (massa/altura) ----
        ax2 = axes[idx, 1]
        for j, lbl2 in enumerate(ring_labels):
            sub2 = df_rings_par[df_rings_par["label"] == lbl2]
            z_vals2 = sub2["z"] * Rchariklo  # km
            hist2, _ = np.histogram(z_vals2, bins=bins_z)
            mass_bins2 = hist2 * mp  # g
            color = ring_colors_[j] if lbl2 == lbl else "lightgray"
            alpha = 0.7 if lbl2 == lbl else 0.3
            ax2.bar((bins_z[:-1]+bins_z[1:])/2, mass_bins2, width=dz, color=color, alpha=alpha, label=f"Anel {lbl2}" if lbl2 == lbl else None)
        ax2.set_xlabel("z (km)")
        ax2.set_ylabel("Massa (g)")
        #ax2.set_yscale("log")
        ax2.set_title(f"Anel {lbl} - densidade vertical")
        ax2.set_xlim(z_min, z_max)
        #ax2.legend()

        # ---- 3) r vs longitude ----
        ax3 = axes[idx, 2]
        for j, lbl2 in enumerate(ring_labels):
            sub2 = df_rings_par[df_rings_par["label"] == lbl2]
            r_km2 = np.sqrt(sub2["x"]**2 + sub2["y"]**2) * Rchariklo
            theta2 = np.arctan2(sub2["y"], sub2["x"])
            theta2_deg = (np.degrees(theta2) + 360) % 360
            color = ring_colors_[j] if lbl2 == lbl else "lightgray"
            alpha = 0.7 if lbl2 == lbl else 0.3
            ax3.scatter(theta2_deg, r_km2, s=2, alpha=alpha, color=color, label=f"Anel {lbl2}" if lbl2 == lbl else None)
        ax3.set_xlabel("Longitude (°)")
        ax3.set_ylabel("r (km)")
        ax3.set_xlim(0, 360)
        ax3.set_ylim(300, r_max)
        ax3.set_title(f"Anel {lbl} - r vs longitude")

        # ---- 4) Plot polar (r vs theta) ----
        # Remove o eixo padrão do grid para o plot polar
        fig.delaxes(axes[idx, 3])
        # Cria eixo polar na mesma posição do grid
        ax4 = fig.add_subplot(axes[idx, 3].get_subplotspec(), projection='polar')
        for j, lbl2 in enumerate(ring_labels):
            sub2 = df_rings_par[df_rings_par["label"] == lbl2]
            r_km2 = np.sqrt(sub2["x"]**2 + sub2["y"]**2) * Rchariklo
            theta2 = np.arctan2(sub2["y"], sub2["x"])
            color = ring_colors_[j] if lbl2 == lbl else "lightgray"
            alpha = 0.7 if lbl2 == lbl else 0.3
            ax4.scatter(theta2, r_km2, s=2, alpha=alpha, color=color, label=f"Anel {lbl2}" if lbl2 == lbl else None)
        #ax3.set_title(f"Anel {lbl} - polar")
        ax4.set_ylim(r_min, r_max)
        ax4.set_rticks([100, 200, 300, 400, 500])
        ax4.set_yticklabels(['100', '', '300', '', '500'])
        ax4.set_rlabel_position(90)
        ax4.text(np.radians(90), 0, "Distance (km)", rotation=0, ha='center', va='bottom', fontsize=12)
        ax4.set_theta_zero_location("E")
        ax4.set_theta_direction(-1)
        #ax3.legend()

    plt.suptitle(f"Snapshot {snapshot_time:.0f}", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{save_prefix}_T{snapshot_time:.0f}.png", dpi=300)
    plt.close(fig)

    df_sector = pd.DataFrame(results)
    df_sector.to_csv(f"{save_prefix}_T{snapshot_time:.0f}_sector_table.csv", index=False)
    return df_sector

def identify_rings(df_diskT, T, eps=0.05, min_samples=30):
    """
    Identifica anéis a partir da distribuição radial das partículas
    usando DBSCAN em 1D no raio.
    
    df_diskT: dataframe de partículas em um snapshot
    eps: largura máxima (em unidades de raio) para considerar partículas no mesmo anel
    min_samples: mínimo de partículas para formar um anel
    
    Retorna:
        df_rings_par: dataframe com partículas etiquetadas com label do anel
        df_rings: propriedades médias dos anéis
    """
    # Raio orbital de cada partícula
    r = np.sqrt(df_diskT["x"]**2 + df_diskT["y"]**2).values.reshape(-1, 1)
    
    # DBSCAN em 1D (raio)
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(r)
    labels = db.labels_
    
    # Salva partículas com label
    df_rings_par = df_diskT.copy()
    df_rings_par["label"] = labels  # -1 = ruído
    
    # Agora calcula propriedades de cada anel
    rings = []
    for lbl in np.unique(labels):
        if lbl == -1:  # ignorar ruído
            continue
        sub = df_rings_par[df_rings_par["label"] == lbl]
        r_vals = np.sqrt(sub["x"]**2 + sub["y"]**2)
        theta_vals = np.arctan2(sub["y"], sub["x"])
        
        rings.append({
            "label": lbl,
            "n_particles": len(sub),
            "r_mean": r_vals.mean(),
            "r_std": r_vals.std(),
            "r_min": r_vals.min(),
            "r_max": r_vals.max(),
            "angular_extent": theta_vals.max() - theta_vals.min()
        })
    
    df_rings = pd.DataFrame(rings)
    # Save each ring as a pickle file
    for ring in rings:
        ring_label = ring["label"]
        ring_particles = df_rings_par[df_rings_par["label"] == ring_label]
        ring_particles.to_pickle(f"{path}Analysis/ring_{ring_label}_particles_T{T:.0f}.pkl")
    return df_rings_par, df_rings

def plot_snapshot_w_tp(ax, T, i, df_disk, df_rings, df_rings_par, df_central, size, point_size_factor):
    #if i >= 3:
    ax.set_xlabel('x (R$_{Chariklo}$)')
    #else:
    #    ax.xaxis.set_tick_params(labelbottom=False)
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
    ring_colors_ = ['tab:blue', 'tab:orange', 'tab:purple']
    # Classify unique labels of rings by their mean distance to the central body (r_mean)
    unique_labels = df_rings.sort_values('r_mean')['label'].values
    for ring_idx, rings_label in enumerate(unique_labels):
        ring_particles = df_rings_par[df_rings_par["label"] == rings_label]
        #print (cluster_particles)

        # Desenha as partículas do cluster com uma cor diferente
        ax.scatter(ring_particles['x'], ring_particles['y'], 
                   s=par_radius * point_size_factor, 
                   color=ring_colors_[ring_idx], label=f'Ring {rings_label}')

    #ax.text(0.9 * size, 0.9 * size, f"orbit {T/2/np.pi:.0f}", ha="right")
    time = T*0.59/24 # T x 0.59 = tempo em horas para Rchariklo = 124 km e Mchariklo = 6.3e18 kg
    ax.set_title(f"{time:.0f} days")

    # Plotando o limite do fluido de Roche
    #ax.add_patch(Circle((e_x, e_y), 2.9, edgecolor='green', facecolor='None'))

    # Plotando as ressonancias orbital/rotação 1:1, 2:1, 3:1
    ax.add_patch(Circle((e_x, e_y), 1.526, edgecolor='green', facecolor='None', lw=0.6, ls='--'))
    ax.add_patch(Circle((e_x, e_y), 2.422, edgecolor='blue', facecolor='None', lw=0.6, ls='--'))
    ax.add_patch(Circle((e_x, e_y), 3.174, edgecolor='red', facecolor='None', lw=0.6, ls='--'))
    return scatter

if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        if arg.startswith("path="):
            path = arg.split("=")[1]

plt.rcParams["axes.labelsize"] = 13
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 13

namefile = path + f'Analysis/plot_rings_snapshots.png'
clusters_file = path + 'Outputs/clusters_data.pkl'
clusters_par_file = path + 'Outputs/clusters_particles_data.pkl'
disk_file = path + 'Outputs/disk_data.pkl'
central_file = path + 'Outputs/central_data.pkl'

df_clusters = pd.read_pickle(clusters_file)
df_clusters_par = pd.read_pickle(clusters_par_file)
df_disk = pd.read_pickle(disk_file)
df_central = pd.read_pickle(central_file)

Time = df_central["time"].to_list()
T_ = [45000, 50000, 60000]
aux = ["(a)", "(b)", "(c)"]
size = 4
point_size_factor = 300

# Criar figura com espaço para a colorbar
fig = plt.figure(figsize=(9, 4))
# Ajustar o layout para deixar espaço para a colorbar à direita
fig.subplots_adjust(bottom=0.12, hspace=0.13, wspace=0.13)
# Criar gridspec para os subplots
gs = gridspec.GridSpec(nrows=1, ncols=3, figure=fig)

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

    # detectar os anéis e salvar em df_rings e df_rings_par
    df_rings_par, df_rings = identify_rings(df_diskT, true_time, eps=0.03, min_samples=100)

    # plotar a densidade longitudinal e angular de partículas nos anéis detectados
    df_sector = analyze_rings_snapshot(df_rings_par, df_rings, true_time,
                                   n_sectors=18,
                                   save_prefix=f"{path}Analysis/rings")
    
    print(f"Plotting snapshot: {T:.2f} | {true_time:.2f} | n_part: {len(df_diskT)} | detected rings: {len(df_rings)}")
    #print(f"numeber of particles sats: {len(df_clusters_parT)} | disk: {len(df_diskT)}")

    scatter = plot_snapshot_w_tp(ax, true_time, i, df_diskT, df_rings, df_rings_par, df_centralT, size, point_size_factor)
    scatters.append(scatter)

# Salva e fecha o gráfico
fig.savefig(namefile, dpi=400, bbox_inches='tight')
plt.close(fig)