import os
import shutil

# Função para modificar o conteúdo do arquivo .c
def change_c_vars(file_path, vars):
    with open(file_path, 'r') as file:
        content = file.readlines()

    # Modificar as variáveis no conteúdo do arquivo
    for i, line in enumerate(content):
        for var, val in vars.items():
            if line.lstrip().startswith(var):
                content[i] = f"\t{var} = {val}; // altered by create_sim \n"

    with open(file_path, 'w') as file:
        file.writelines(content)

# Caminho do arquivo .c original
common_path = "example/openmp_j2_c22/"

set_ = ['chariklo_tp', 'chariklo_1e-3_e.1', 'chariklo_1e-3_e.5', 'chariklo_1e-4_e.1', 'chariklo_1e-4_e.5', 'chariklo_1e-3_e.1_ag', 'chariklo_1e-3_e.5_ag', 'chariklo_1e-4_e.1_ag', 'chariklo_1e-4_e.5_ag']
disc_ = ['chariklo_tp', 'chariklo_1e-3', 'chariklo_1e-3', 'chariklo_1e-4', 'chariklo_1e-4', 'chariklo_1e-3', 'chariklo_1e-3', 'chariklo_1e-4', 'chariklo_1e-4']
type_ = ['reb_collision_resolve_hardsphere_merge_central' for _ in range(9)]
threads_ = [16 for _ in range(9)]
cr_ = [0, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5]
tstep_ = [5e-3 for _ in range(9)]
out_int_ = [6 for _ in range(9)]
softening_ = [2.9256932e-03*0.9, 2.9256932e-03*0.9, 9.35990162e-03*0.9, 1.357986495e-03*0.9, 1.357986495e-03*0.9, 2.9256932e-03*0.9, 9.35990162e-03*0.9, 1.357986495e-03*0.9, 1.357986495e-03*0.9]
N_active_ = [1 for _ in range(5)]+['r->N' for _ in range(4)]
gravity_ = ['REB_GRAVITY_BASIC' for _ in range(5)]+['REB_GRAVITY_TREE' for _ in range(4)]
change_central_ = [1 for _ in range(5)]+[0 for _ in range(4)]
remove_tp_flag_ = [1]+[0 for _ in range(8)] # 1 to remove test particles with q inside CB, 0 to non-test particles

# Copiar o arquivo e modificar as variáveis em cada pasta de destino
for set, disc, cr, type, threads, dt, out_int, soft, N_active, gravity, change_central, tp_flag in zip(set_, disc_, cr_, type_, threads_, tstep_, 
                                                           out_int_, softening_, N_active_, gravity_, change_central_, remove_tp_flag_):
    # Criar a pasta se não existir
    if not os.path.exists(set):
        os.makedirs(set)
        print(f"{set} created!")
    else:
        print(f"{set} already exists!")

    vars_in_c = {
    'entrada': f'fopen("../../Inputs/disc_{disc}.in", "r")',
    'r->gravity': gravity,
    'const double tmax': '30000',
    'r->softening': soft,
    'double cr': cr,
    'r->collision_resolve': type,
    'int np': threads,
    'double output_interval': out_int,
    'r->dt': dt,
    'r->N_active': N_active,
    'r->central_to_zero': change_central,
    'int test_particle_flag': tp_flag
    }

    for i in range(1):
        sim_name = f"{set}/sim_{i:02d}"
        if not os.path.exists(sim_name):
            os.makedirs(sim_name)
            print(f"\t{sim_name} created!")
        else:
            print(f"\t{sim_name} already exists!")
        
        # Copiando arquivos
        c_file_copy = os.path.join(sim_name, 'problem.c')
        #c_file_copy2 = os.path.join(sim_name, 'problem_restart.c')

        shutil.copyfile(os.path.join(common_path, 'problem.c'), c_file_copy)
        shutil.copyfile(os.path.join(common_path, 'makefile'), os.path.join(sim_name, 'makefile'))

        # Modificar as variáveis no arquivo copiado
        change_c_vars(c_file_copy, vars_in_c)
        #change_c_vars(c_file_copy2, {k: v for k, v in vars_in_c.items() 
        #                             if k != 'entrada'})

print("\nEnd!")