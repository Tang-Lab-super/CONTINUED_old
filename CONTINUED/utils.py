#################################################### Continued pipeline ########################################
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import cv2
from tqdm import tqdm

def check_makedir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)

def preprocess(df_desi):
    all_mz = df_desi.columns[~df_desi.columns.str.startswith('Unnamed')]
    new_columns = pd.Index(['x', 'y']).append(all_mz)
    df_desi = df_desi.iloc[:, 1:]
    df_desi.columns = new_columns
    df_desi.index = [f'Bin_{x}' for x in df_desi.index]
    return df_desi

def load_raw_desi_data(desi_lm_file, desi_unlm_file, num_mz=3000):
    df_desi_lm = pd.read_csv(desi_lm_file, sep='\t', skiprows=3, index_col=False)
    df_desi_unlm = pd.read_csv(desi_unlm_file, sep='\t', skiprows=3, index_col=False)

    df_desi_lm = preprocess(df_desi_lm)
    df_desi_unlm = preprocess(df_desi_unlm)

    df_desi_lm = df_desi_lm.iloc[:, :(num_mz+2)]
    df_desi_unlm = df_desi_unlm.iloc[:, :(num_mz+2)]
    return df_desi_lm, df_desi_unlm

def check_tissue_mz(mz_tissue_type, mz_tissue, thresh, df_desi_unlm, df_desi_lm, cmap, to_dir):
    # visualization the tissue mz
    if mz_tissue_type == 'unlm':
        df_plot = df_desi_unlm.copy()
    else:
        df_plot = df_desi_lm.copy()

    dic_x = dict(zip(np.sort(df_plot['x'].unique()), range(len(np.sort(df_plot['x'].unique())))))
    dic_y = dict(zip(np.sort(df_plot['y'].unique()), range(len(np.sort(df_plot['y'].unique())))))
    df_plot['x_array'] = df_plot['x'].map(dic_x)
    df_plot['y_array'] = df_plot['y'].map(dic_y)

    spot_row = np.max(df_plot['x_array']) - np.min(df_plot['x_array'])
    spot_col = np.max(df_plot['y_array']) - np.min(df_plot['y_array'])

    values = df_plot[mz_tissue]

    fig, ax = plt.subplots(figsize=(7 * spot_row / spot_col, 7))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    ax.scatter(x=df_plot['x_array']+1, y=df_plot['y_array']+1, c=values, s=marker_size, edgecolors='black', linewidth=0, cmap=cmap)
    # ax.scatter(x=tissue_edge[1]+1, y=tissue_edge[0]+1, c="#bdbdbd", s=marker_size/2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{to_dir}/1.1.check.tissue_mz.pdf')


    df_plot[mz_tissue][df_plot[mz_tissue] < thresh] = 0
    df_plot[mz_tissue][df_plot[mz_tissue] >= thresh] = 1
    values = df_plot[mz_tissue]
    fig, ax = plt.subplots(figsize=(7 * spot_row / spot_col, 7))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    ax.scatter(x=df_plot['x_array']+1, y=df_plot['y_array']+1, c=values, s=marker_size, edgecolors='black', linewidth=0, cmap=cmap)
    # ax.scatter(x=tissue_edge[1]+1, y=tissue_edge[0]+1, c="#bdbdbd", s=marker_size/2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{to_dir}/1.2.check.tissue_mz.thresh.pdf')

def tissue_detection(df_desi, mz_tissue, to_dir, thresh=500, otsu=False, dilate_size=2, tissue_erode_size=1): 
    df_tissue = df_desi.pivot(index='y', columns='x', values=mz_tissue)
    if otsu:
        df_tissue = 255 * df_tissue / df_tissue.max().max()
        df_tissue = df_tissue.astype(np.uint8)
        cv2.imwrite(f'{to_dir}/2.1.tissue_mz.png', df_tissue.values)
        # 使用cv2.THRESH_OTSU 自动寻找最佳阈值，并检查二值化结果
        ret, tissue_thresh = cv2.threshold(df_tissue.values, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
        cv2.imwrite(f'{to_dir}/2.2.tissue.otsu.png', tissue_thresh)
    else:
        df_tissue[df_tissue > thresh] = 100000000
        df_tissue[df_tissue <= thresh] = 0
        df_tissue = 255 * df_tissue / df_tissue.max().max()
        df_tissue = df_tissue.astype(np.uint8)
        cv2.imwrite(f'{to_dir}/2.1.tissue_mz.png', df_tissue.values)
        # 使用cv2.THRESH_OTSU 自动寻找最佳阈值，并检查二值化结果
        tissue_thresh = df_tissue.values.astype(np.uint8)
        cv2.imwrite(f'{to_dir}/2.2.tissue.thresh.png', tissue_thresh)

    # 膨胀收缩后找最大连通域定于为组织区域
    kernel_dilate_erode = cv2.getStructuringElement(cv2.MORPH_RECT, (dilate_size, dilate_size))
    dilated = cv2.dilate(tissue_thresh, kernel_dilate_erode)
    eroded = cv2.erode(dilated, kernel_dilate_erode)
    cv2.imwrite(f'{to_dir}/2.3.tissue.eroded.png', eroded)

    contours, hierarchy_ = cv2.findContours(eroded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    max_area = (0, 0)
    for i in range(len(contours)):
        area = cv2.contourArea(contours[i])
        if area > max_area[1]:
            max_area = (i, area)
    max_contour = contours[max_area[0]]
    mask = np.zeros((tissue_thresh.shape[0], tissue_thresh.shape[1], 3))
    tissue = cv2.drawContours(mask, [max_contour], -1, (255, 255, 255), -1)
    tissue = tissue.astype(np.uint8)
    gray = cv2.cvtColor(tissue, cv2.COLOR_BGR2GRAY)
    cv2.imwrite(f'{to_dir}/2.4.tissue.detection.raw.png', gray)

    df_tissue_raw = pd.DataFrame(gray, index=range(gray.shape[0]), columns=range(gray.shape[1]))
    df_tissue_raw.to_csv(f'{to_dir}/2.result.tissue.detection.raw.txt', sep='\t')

    # 最终将组织收缩5个单位，以去除边界的污染信号

    kernel_erode_tissue = cv2.getStructuringElement(cv2.MORPH_RECT, (tissue_erode_size, tissue_erode_size))
    eroded_tissue = cv2.erode(gray, kernel_erode_tissue)
    cv2.imwrite(f'{to_dir}/2.5.tissue.eroded.png', eroded_tissue)
    df_tissue_eroded = pd.DataFrame(eroded_tissue, index=range(eroded_tissue.shape[0]), columns=range(eroded_tissue.shape[1]))
    df_tissue_eroded.to_csv(f'{to_dir}/2.result.tissue.detection.eroded.txt', sep='\t')
    return eroded_tissue

def remove_noise_signal(df_desi_unlm, tissue_mask, num=90, bg_percent=0.5): 
    true_signal = []
    noise_signal = []
    with tqdm(total = df_desi_unlm.shape[1] - 2) as pbar:
        for j in range(2, df_desi_unlm.shape[1]):
            tis = df_desi_unlm.iloc[:,[0,1,j]]
            per = np.percentile(tis.iloc[:,2],num)
            hmat = tis.pivot(index='y', columns='x', values=tis.columns[2]) 
            hmat.index = range(np.shape(hmat)[0])
            hmat.columns = range(np.shape(hmat)[1])
            hmat.iloc[:,0:5] = 0 
            count90 = np.sum(np.array(hmat>=per))  
            hmat[hmat<per] = 0
            hmat = tissue_mask*hmat 
            count_tis = np.sum(np.array(hmat>0))
            if count_tis/count90 > bg_percent: 
                true_signal.append(df_desi_unlm.columns[j])
            else:
                noise_signal.append(df_desi_unlm.columns[j])
            pbar.update(1)
    final_cols = ['x', 'y'] + true_signal
    df_desi_unlm = df_desi_unlm[final_cols]
    return noise_signal, true_signal, df_desi_unlm

def create_desi_obj(df_final):
    df_meta = df_final.iloc[:, :4]
    df_meta['x_array_plot'] = df_meta['x_array'] - df_meta['x_array'].min() + 1
    df_meta['y_array_plot'] = df_meta['y_array'] - df_meta['y_array'].min() + 1
    df_value = df_final.iloc[:, 4:]
    adata = anndata.AnnData(X=df_value, obs=df_meta)
    adata.obsm["spatial"] = df_meta[['x_array', 'y_array']].values
    return adata

def desi_clustering(adata, resolution, prefix, col16):
    sc.pp.scale(adata)
    sc.pp.highly_variable_genes(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    # sc.tl.umap(adata)
    sc.tl.leiden(adata, n_iterations=2, resolution=resolution)
    adata.obs['clustering.scanpy'] = [int(x) for x in adata.obs['leiden']]
    spatialdimplot_desi(adata, 'clustering.scanpy', prefix, col16)
    return adata

def get_final_matrix(df_desi_unlm_filter, tissue_mask):
    dic_x = dict(zip(np.sort(df_desi_unlm_filter['x'].unique()), range(len(np.sort(df_desi_unlm_filter['x'].unique())))))
    dic_y = dict(zip(np.sort(df_desi_unlm_filter['y'].unique()), range(len(np.sort(df_desi_unlm_filter['y'].unique())))))
    df_desi_unlm_filter['x_array'] = df_desi_unlm_filter['x'].map(dic_x)
    df_desi_unlm_filter['y_array'] = df_desi_unlm_filter['y'].map(dic_y)
    df_desi_unlm_filter['pos_tmp'] = df_desi_unlm_filter['x_array'].astype('str').str.cat(df_desi_unlm_filter['y_array'].astype('str'), sep='_')
    tissue_pos_tmp = np.where(tissue_mask!=0)
    tissue_pos = set([f'{tissue_pos_tmp[1][x]}_{tissue_pos_tmp[0][x]}' for x in range(len(tissue_pos_tmp[0]))])
    df_desi_unlm_final= df_desi_unlm_filter[df_desi_unlm_filter['pos_tmp'].isin(tissue_pos)]
    df_desi_unlm_final.drop('pos_tmp', axis=1, inplace=True)
    new_columns = ['x', 'y', 'x_array', 'y_array'] + list(df_desi_unlm_final.columns[df_desi_unlm_final.columns.str.startswith('mz_')])
    df_desi_unlm_final = df_desi_unlm_final[new_columns]
    return df_desi_unlm_final

def spatialfeatureplot_desi(adata, gene, prefix, cmap):
    spot_row = np.max(adata.obs['x_array_plot']) - np.min(adata.obs['x_array_plot'])
    spot_col = np.max(adata.obs['y_array_plot']) - np.min(adata.obs['y_array_plot'])

    values = adata.X[:, adata.var.index == gene][:, 0]

    fig, ax = plt.subplots(figsize=(7 * spot_row / spot_col, 7))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    ax.scatter(x=adata.obs['x_array_plot']+1, y=adata.obs['y_array_plot']+1, c=values, s=marker_size, edgecolors='black', linewidth=0, cmap=cmap)
    # ax.scatter(x=tissue_edge[1]+1, y=tissue_edge[0]+1, c="#bdbdbd", s=marker_size/2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}.pdf')

def spatialdimplot_desi(adata, label, prefix, col16):

    spot_row = np.max(adata.obs['x_array_plot']) - np.min(adata.obs['x_array_plot'])
    spot_col = np.max(adata.obs['y_array_plot']) - np.min(adata.obs['y_array_plot'])

    values = adata.obs[label]
    colors = [col16[x] for x in values]

    fig, ax = plt.subplots(figsize=(7 * spot_row / spot_col, 7))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    ax.scatter(x=adata.obs['x_array_plot']+1, y=adata.obs['y_array_plot']+1, c=colors, s=marker_size, edgecolors='black', linewidth=0)
    # ax.scatter(x=tissue_edge[1]+1, y=tissue_edge[0]+1, c="#bdbdbd", s=marker_size/2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}.pdf')  

def get_curcontour_pixel_num(raw_img, cur_contour):
    # image,contours,contourIdx,color,thickness
    # contourldx为负值绘制所有轮廓
    # thickness<0填充连通域，thickness>0绘制轮廓
    mask = np.zeros(raw_img.shape)
    cv2.drawContours(mask, cur_contour, -1, 255, -1)
    return(np.sum(mask != 0))

def select_contour(raw_img, contours, threshold):
    selected_countours = []
    for cur_counter in contours:
        if get_curcontour_pixel_num(raw_img, cur_counter) > threshold:
            selected_countours.append(cur_counter)
    return selected_countours

def plot_cluster_border(adata, prefix, threshold_border, col16):
    spot_row = np.max(adata.obs['x_array_plot']) - np.min(adata.obs['x_array_plot'])
    spot_col = np.max(adata.obs['y_array_plot']) - np.min(adata.obs['y_array_plot'])

    #### 单个cluster
    for cluster in adata.obs['clustering.scanpy'].unique():
        fig, ax = plt.subplots(figsize=(7 * spot_row / spot_col, 7))
        ax.set_xlim((0 ,spot_row+3))
        ax.set_ylim((0, spot_col+3))

        r = 0.55
        r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
        marker_size = np.pi * r_**2

        df_plot = adata.obs[['x_array_plot', 'y_array_plot', 'clustering.scanpy']]
        df_plot['clustering.scanpy'][df_plot['clustering.scanpy']==cluster] = 255
        df_plot['clustering.scanpy'][df_plot['clustering.scanpy']!=255] = 0

        img = df_plot.pivot_table(columns='x_array_plot', index='y_array_plot', values='clustering.scanpy')
        img.fillna(0, inplace=True)
        img = img.astype(np.uint8)
        img = img.values

        contours, _ = cv2.findContours(img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        selected_contours = select_contour(np.zeros(img.shape), contours, threshold_border)

        mask = np.zeros(img.shape)
        cv2.drawContours(mask, selected_contours, -1, 255, 1)

        df_border = pd.DataFrame(mask)
        df_border = df_border.stack().reset_index()
        df_border.columns = ['y', 'x', 'value']

        df_border['x'] = df_border['x'] + 1
        df_border['y'] = df_border['y'] + 1
        df_border = df_border[df_border['value'] != 0]

        ax.scatter(x=df_border['x'], y=df_border['y'], c=col16[cluster], s=marker_size, edgecolors='black', linewidth=0)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)

        ax.figure.savefig(f'{prefix}.cluster.{cluster}.pdf')


    #### 所有cluster
    fig, ax = plt.subplots(figsize=(7 * spot_row / spot_col, 7))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2

    for cluster in adata.obs['clustering.scanpy'].unique():
        df_plot = adata.obs[['x_array_plot', 'y_array_plot', 'clustering.scanpy']]
        df_plot['clustering.scanpy'][df_plot['clustering.scanpy']==cluster] = 255
        df_plot['clustering.scanpy'][df_plot['clustering.scanpy']!=255] = 0

        img = df_plot.pivot_table(columns='x_array_plot', index='y_array_plot', values='clustering.scanpy')
        img.fillna(0, inplace=True)
        img = img.astype(np.uint8)
        img = img.values

        contours, _ = cv2.findContours(img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        selected_contours = select_contour(np.zeros(img.shape), contours, threshold_border)

        mask = np.zeros(img.shape)
        cv2.drawContours(mask, selected_contours, -1, 255, 1)

        df_border = pd.DataFrame(mask)
        df_border = df_border.stack().reset_index()
        df_border.columns = ['y', 'x', 'value']

        df_border['x'] = df_border['x'] + 1
        df_border['y'] = df_border['y'] + 1
        df_border = df_border[df_border['value'] != 0]

        ax.scatter(x=df_border['x'], y=df_border['y'], c=col16[cluster], s=marker_size, edgecolors='black', linewidth=0)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}.cluster.merge.pdf')