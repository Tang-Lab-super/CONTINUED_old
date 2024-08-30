import os
import numpy as np
from scipy.stats import gaussian_kde


def Print_Clustered_Mass_By_Sample(mass_clustered, mass_sample, lipid, small_mol, output_prefix):
    import os

    # get sample names
    samples = {}
    for m in mass_sample:
        for s in mass_sample[m]:
            samples[s] = samples.get(s, 0) + 1

    samples = sorted(samples.keys())

    output_tab = f'{output_prefix}.clustered_mass.table.with.anno.txt'
    with open(output_tab, 'w') as fh_out:
        header_part_1 = samples
        header_part_2 = [f'{s}_mass' for s in samples]
        header_part_3 = ['anno_lipid', 'anno_small_mol']
        fh_out.write('\t'.join(['Index'] + header_part_1 + ['Cluster_type', 'Mass_num'] + header_part_2 + header_part_3) + '\n')

        indices = mass_clustered.keys()

        for index in indices:
            sample_c_mass = {s: [] for s in samples}
            sample_c_mass_flag = {s: 0 for s in samples}
            anno_lipid = []
            anno_small_mol = []

            tag = mass_clustered[index]['tag']
            mass = mass_clustered[index]['mass']

            for m in mass:
                for s in samples:
                    if m in mass_sample and s in mass_sample[m]:
                        sample_c_mass_flag[s] += 1
                        sample_c_mass[s].append(m)

                if m in lipid:
                    anno_lipid.append(';'.join([str(m),
                                                ','.join(lipid[m]['IonFormula']),
                                                ','.join(lipid[m]['LipidIon'])]))

                if m in small_mol:
                    anno_small_mol.append(';'.join([str(m),
                                                    ','.join(small_mol[m]['Formula']),
                                                    ','.join(small_mol[m]['Tag']),
                                                    ','.join(small_mol[m]['Name'])]))

            output_flag = [1 if sample_c_mass_flag[s] > 0 else 0 for s in samples]
            output_mass = []
            # [','.join(sample_c_mass[s]) if sample_c_mass[s] else 'None' for s in samples]
            for s in samples:
                if len(sample_c_mass[s]) == 0:
                    output_mass.append('None')
                else:
                    output_mass.append(','.join([str(x) for x in sample_c_mass[s]]))

            output_anno_lipid = '|'.join(anno_lipid) if anno_lipid else 'None'
            output_anno_small_mol = '|'.join(anno_small_mol) if anno_small_mol else 'None'

            num = sum(output_flag)

            fh_out.write('\t'.join([str(index)] + list(map(str, output_flag)) + [tag, str(num)] + output_mass + [output_anno_lipid, output_anno_small_mol]) + '\n')

    return 0

def Clustering_Mass_by_KDE(mass_index_group, lipid, small_mol, mass_cutoff):
    from copy import deepcopy
    clustered_mass = {}

    weight_lipid = 1
    weight_small_mol = 1
    weight_mass = 1

    for index in sorted(mass_index_group.keys()):
        mass = sorted(mass_index_group[index])

        mass_diff = mass[-1] - mass[0]
        mass_num = len(mass)

        if mass_num == 1:
            if mass[0] in lipid:
                tag = 'Solo_Lipid'
            elif mass[0] in small_mol:
                tag = 'Solo_Small_Mol'
            else:
                tag = 'Solo_Mass'
            clustered_mass[index] = {'mass': mass, 'num': mass_num, 'tag': tag}
        
        elif 1 < mass_num <= 3 and mass_diff <= mass_cutoff:
            clustered_mass[index] = {'mass': mass, 'num': mass_num, 'tag': 'Cluster_No_KDE'}
        
        elif 1 < mass_num <= 3 and mass_diff > mass_cutoff:
            mass_split = []
            index_split = 0
            for i in range(mass_num):
                if i == 0:
                    mass_split.append(mass[i])
                elif mass[i] - mass[i - 1] <= mass_cutoff:
                    mass_split.append(mass[i])
                elif mass[i] - mass[i - 1] > mass_cutoff:
                    clustered_mass[f"{index}_{index_split}"] = {'mass': mass_split.copy(), 'num': len(mass_split), 'tag': 'Cluster_No_KDE'}
                    mass_split = [mass[i]]
                    index_split += 1
                else:
                    print('unknown case:', mass)
                    exit(1)
            clustered_mass[f"{index}_{index_split}"] = {'mass': mass_split.copy(), 'num': len(mass_split), 'tag': 'Cluster_No_KDE'}
        
        elif mass_num > 3 and mass_diff <= mass_cutoff:
            clustered_mass[index] = {'mass': mass, 'num': mass_num, 'tag': 'Cluster_No_KDE'}
        
        elif mass_num > 3 and mass_diff > mass_cutoff:
            pd = Kernel_Density_Estimate(mass, weight_lipid, weight_small_mol, weight_mass, lipid, small_mol)
            pd_bak = deepcopy(pd)

            clusters = Find_Probability_Density_Regions(pd_bak)
            Print_KDE_and_Cluster(pd, clusters, index, mass)

            mass_grouped = []
            for cluster in clusters:
                mass_start = cluster[0]
                mass_end = cluster[1]
                mm = [m for m in mass if mass_start <= m <= mass_end]
                if mm:
                    mass_grouped.append(mm)

            for i, m_grouped in enumerate(mass_grouped):
                m_diff = m_grouped[-1] - m_grouped[0]
                m_num = len(m_grouped)
                if m_num == 1:
                    if m_grouped[0] in lipid:
                        clustering_tag = 'KDE_Solo_Lipid'
                    elif m_grouped[0] in small_mol:
                        clustering_tag = 'KDE_Solo_Small_Mol'
                    else:
                        clustering_tag = 'KDE_Solo_Mass'
                elif m_num > 1 and m_diff <= mass_cutoff:
                    clustering_tag = 'KDE_No_Exceed'
                elif m_num > 1 and m_diff > mass_cutoff:
                    clustering_tag = 'KDE_Exceed'
                else:
                    print('unknown case: KDE clustering!')
                    print(f'mass_num: {m_num}, mass_diff: {m_diff}')
                    exit(1)
                clustered_mass[f"{index}_{i}"] = {'mass': m_grouped, 'num': m_num, 'tag': clustering_tag}
        else:
            print('unknown case: grouping!')
            print(f'mass_num: {mass_num}, mass_diff: {mass_diff}')
            exit(1)

    return clustered_mass

def Find_Probability_Density_Regions(pd):
    from copy import deepcopy

    # more than 20% of summit ( pd_max )
    peak_change_cutoff = 0.2

    # fixed bug here and clean the codes
    # up tendency
    if pd[1][1] >= pd[0][1]:
        regions = []

        pd_min_left = pd.pop(0)
        pd_min_right = pd_min_left
        pd_max = pd_min_left

        while pd:
            e = pd.pop(0)
            if pd:
                if not pd_min_left and not pd_min_right and not pd_max:
                    pd_min_left = e
                    pd_min_right = e
                    pd_max = e
                elif e[1] >= pd_min_right[1] and e[1] <= pd[0][1]:
                    pd_min_right = e
                    pd_max = e
                elif e[1] >= pd_min_right[1] and e[1] > pd[0][1]:
                    pd_min_right = e
                    pd_max = e
                elif e[1] <= pd_min_right[1] and e[1] >= pd[0][1]:
                    pd_min_right = e
                elif e[1] <= pd_min_right[1] and e[1] <= pd[0][1]:
                    pd_min_right = e
                    regions.append([pd_min_left, pd_min_right, pd_max])
                    pd_min_left = None
                    pd_min_right = None
                    pd_max = None
                else:
                    print(f"unknown case: up tendency\n{pd_min_left}; {pd_max}; {pd_min_right}")
                    print(pd)
                    exit(1)
            else:
                # fixed the bug here
                if pd_min_left and pd_max:
                    pd_min_right = e
                    regions.append([pd_min_left, pd_min_right, pd_max])
                elif not pd_min_left and not pd_max:
                    pd_min_left = e
                    pd_min_right = e
                    pd_max = e
                    regions.append([pd_min_left, pd_min_right, pd_max])
                else:
                    print('up tendency: unknown case')
                    exit(1)

        region_merged = Merge_Regios(regions, peak_change_cutoff)
        return region_merged

    elif pd[1][1] < pd[0][1]:
        regions = []

        pd_min_left = pd.pop(0)
        pd_min_right = pd_min_left
        pd_max = pd_min_left

        while pd:
            e = pd.pop(0)
            if pd:
                if not pd_min_left and not pd_min_right and not pd_max:
                    pd_min_left = e
                    pd_min_right = e
                    pd_max = e
                elif e[1] <= pd_min_right[1] and e[1] >= pd[0][1]:
                    pd_min_right = e
                elif e[1] <= pd_min_right[1] and e[1] < pd[0][1]:
                    pd_min_right = e
                    regions.append([pd_min_left, pd_min_right, pd_max])
                    pd_min_left = None
                    pd_min_right = None
                    pd_max = None
                elif e[1] >= pd_min_right[1] and e[1] <= pd[0][1]:
                    pd_min_right = e
                    pd_max = e
                elif e[1] >= pd_min_right[1] and e[1] >= pd[0][1]:
                    pd_min_right = e
                    pd_max = e
                else:
                    print(f"unknown case: down tendency\n{pd_min_left}; {pd_max}; {pd_min_right}")
                    print(pd)
                    exit(1)
            else:
                # fixed the bug here
                if pd_min_left and pd_max:
                    pd_min_right = e
                    regions.append([pd_min_left, pd_min_right, pd_max])
                elif not pd_min_left and not pd_max:
                    pd_min_left = e
                    pd_min_right = e
                    pd_max = e
                    regions.append([pd_min_left, pd_min_right, pd_max])
                else:
                    print('down tendency: unknown case')
                    exit(1)

        region_merged = Merge_Regios(regions, peak_change_cutoff)
        return region_merged

    else:
        print("unconsidered case:")
        print('NOT: pd[1][1] > pd[0][1]')
        print('print the 1st and 2nd intensity values:')
        print(f'{pd[1][1]}\t{pd[0][1]}')
        exit(1)

def Merge_Regios(regions, peak_change_cutoff):
    region_merged = []
    peak_1st = regions.pop(0)

    summit_density = peak_1st[2][1]
    right_side_density = peak_1st[1][1]

    region_merge_left = peak_1st[0][0]
    region_merge_right = peak_1st[1][0]

    while regions:
        r = regions.pop(0)
        if (summit_density - right_side_density) <= summit_density * peak_change_cutoff:
            region_merge_right = r[1][0]
            summit_density = r[2][1]
            right_side_density = r[1][1]
        else:
            region_merged.append([region_merge_left, region_merge_right])
            region_merge_left = r[0][0]
            region_merge_right = r[1][0]
            summit_density = r[2][1]
            right_side_density = r[1][1]

    region_merged.append([region_merge_left, region_merge_right])

    return region_merged

def Kernel_Density_Estimate(mass, weight_lipid, weight_small_mol, weight_mass, lipid, small_mol):
    bin_num = 200
    weights = []

    # Create the data array and the corresponding weights array
    data = np.array(mass)
    for m in mass:
        if m in lipid:
            weights.append(weight_lipid)
        elif m in small_mol:
            weights.append(weight_small_mol)
        else:
            weights.append(weight_mass)
    
    weights = np.array(weights)

    # Create the KDE with weights
    kde = gaussian_kde(data, weights=weights, bw_method='scott')

    min_val = data.min()
    max_val = data.max()
    x = np.linspace(min_val, max_val, bin_num)
    pdf = kde(x)

    return list(zip(x, pdf))

def Group_Mass(mass, lipid, small_mol, mass_cutoff):
    mass_all = list(mass.keys()) + list(lipid.keys()) + list(small_mol.keys())
    mass_all = sorted(mass_all)

    mass_index_group = {}
    mass_group = []
    index = 0
    mass_1st = None

    for i in range(len(mass_all)):
        if mass_1st is None and not mass_group:
            mass_1st = mass_all[i]
            mass_group.append(mass_all[i])
        else:
            if (mass_all[i] - mass_all[i-1]) <= mass_cutoff:
                mass_group.append(mass_all[i])
            elif (mass_all[i] - mass_all[i-1]) > mass_cutoff:
                index += 1
                mass_index_group[index] = mass_group.copy()
                mass_1st = None
                mass_group = []
                mass_1st = mass_all[i]
                mass_group.append(mass_all[i])
            else:
                raise ValueError(f"Unknown case: {mass_all[i-1]}, {mass_all[i]}, {mass_all[i+1]}")
    
    index += 1
    mass_index_group[index] = mass_group.copy()

    return mass_index_group

def Print_Mass_Diff_By_Samples(sample_mass, output):
    with open(output, 'w') as fh_out:
        fh_out.write("\t".join(['mass_dis', 'sample']) + "\n")
        
        for sample in sorted(sample_mass.keys()):
            mass = sorted(sample_mass[sample].keys())
            for i in range(len(mass) - 1):
                dis = mass[i + 1] - mass[i]
                fh_out.write("\t".join([str(dis), sample]) + "\n")
    
    return 0

def Parsing_Small_Molecule(input_file):
    small_mol = {}
    
    mass_H = 1.0078
    mass_HCOO = 44.9977
    
    with open(input_file, 'r') as fh_in:
        for line in fh_in:
            if line.startswith('Name'):
                continue
            line = line.strip()
            name, formula, mass = line.split('\t')
            mass = float(mass)
            
            # Add H
            if (mass - mass_H) not in small_mol:
                small_mol[mass - mass_H] = {'Name': [], 'Formula': [], 'Tag': ''}
            small_mol[mass - mass_H]['Name'].append(name)
            small_mol[mass - mass_H]['Formula'].append(formula)
            small_mol[mass - mass_H]['Tag'] = 'H'
            
            # Add HCOO
            if (mass + mass_HCOO) not in small_mol:
                small_mol[mass + mass_HCOO] = {'Name': [], 'Formula': [], 'Tag': ''}
            small_mol[mass + mass_HCOO]['Name'].append(name)
            small_mol[mass + mass_HCOO]['Formula'].append(formula)
            small_mol[mass + mass_HCOO]['Tag'] = 'HCOO'
    
    return small_mol

def Parsing_Lipid(input_file):
    lipid = {}
    
    with open(input_file, 'r') as fh_in:
        for line in fh_in:
            if line.startswith('LipidIon'):
                continue
            line = line.strip()
            lipidIon, mass, formula = line.split('\t')
            mass = float(mass)
            
            if mass not in lipid:
                lipid[mass] = {'LipidIon': [], 'IonFormula': []}
            lipid[mass]['LipidIon'].append(lipidIon)
            lipid[mass]['IonFormula'].append(formula)
    
    return lipid

def Parsing_Mass_Table(input_list, dir_sample):
    sample_mass = {}
    mass_sample = {}
    mass = {}

    with open(input_list, 'r') as fh_in_1:
        for sample_name in fh_in_1:
            sample_name = sample_name.strip()
            input_mass = os.path.join(dir_sample, sample_name, f'{sample_name}_3k_signal.lock_mass.txt')
            
            with open(input_mass, 'r') as fh_in_2:
                # Read only the first line
                line = fh_in_2.readline().strip()
                list_items = line.split('\t')

                for item in list_items[3:]:
                    item = float(item)
                    if sample_name not in sample_mass:
                        sample_mass[sample_name] = {}
                    if item not in sample_mass[sample_name]:
                        sample_mass[sample_name][item] = 0
                    sample_mass[sample_name][item] += 1

                    if item not in mass_sample:
                        mass_sample[item] = {}
                    if sample_name not in mass_sample[item]:
                        mass_sample[item][sample_name] = 0
                    mass_sample[item][sample_name] += 1

                    if item not in mass:
                        mass[item] = 0
                    mass[item] += 1

    return sample_mass, mass_sample, mass

def Print_KDE_and_Cluster(pd, clusters, index, mass):
    output_file = f'plot.kde.index_{index}.txt'
    
    with open(output_file, 'w') as fh_out:
        kde_mass = []
        for item in pd:
            kde_mass.append(item[1])
            fh_out.write(f"{item[0]}\t{item[1]}\tKDE\n")
        
        kde_mass.sort()
        
        for cluster in clusters:
            fh_out.write(f"{cluster[0]}\t{kde_mass[0]}\tCluster\n")
            fh_out.write(f"{cluster[1]}\t{kde_mass[0]}\tCluster\n")
        
        for m in mass:
            fh_out.write(f"{m}\t{kde_mass[0]}\tMass\n")
    
    return 0
