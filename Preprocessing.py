import pandas as pd
import cv2 as cv
import cv2
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def preprocessing(data15):
    """
    :param df:raw DataFrame
    """
    data15 = data15.drop([0])
    data15 = data15.drop(columns='Unnamed: 0')
    newcname = data15.columns[3:]  
    last_column = data15.columns[-1]
    data15 = data15.drop(columns=last_column)
    newcname = newcname.insert(0, 'y')
    newcname = newcname.insert(0, 'x')
    data15.columns = newcname
    data15.index = range(data15.shape[0])
    return data15

def get_raw_tissue(data15, mz, threshold):
    tis = data15.loc[:,('x','y',mz)]
    tis.loc[tis[mz]<=threshold, mz]=0    #背景调0
    tis.loc[tis[mz]>threshold, mz]=1000
    hmat = tis.pivot(index='y', columns='x', values=mz)
    hmat.index = range(np.shape(hmat)[0])
    hmat.columns = range(np.shape(hmat)[1])    
    hmat.to_csv('tis{}.txt'.format(threshold), sep='\t')
    return hmat

    
def get_max_contour(binarized_tissue, dilate_size):
    thresh = binarized_tissue.values
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (dilate_size, dilate_size))
    dilated = cv2.dilate(thresh, kernel)
    eroded = cv2.erode(dilated, kernel)
    eroded = eroded.astype(np.uint8)
    contours, hierarchy_ = cv2.findContours(eroded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    max_area = (0, 0)
    for i in range(len(contours)):
        area = cv2.contourArea(contours[i])
        if area > max_area[1]:
            max_area = (i, area)
    max_contour = contours[max_area[0]]
    mask = np.zeros((thresh.shape[0], thresh.shape[1], 3))
    tissue = cv.drawContours(mask, [max_contour], -1, (255, 255, 255), -1)
    tissue = tissue.astype(np.uint8)
    gray = cv.cvtColor(tissue, cv.COLOR_BGR2GRAY)
    cv.imwrite('gray.png', gray)
    df = pd.DataFrame(gray, index=range(gray.shape[0]), columns=range(gray.shape[1]))
    df.to_csv('raw_tissue.txt', sep='\t')
    print('Step 2: binarized tissue dilated and eroded has been done!')
    print('used dilate size = {}'.format(str(dilate_size)))
    return df


def erode_tissue(kernel_size):
    img = cv.imread('gray.png')
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
    eroded = cv.erode(img, kernel)
    eroded_gray = cv.cvtColor(eroded, cv.COLOR_BGR2GRAY)
    cv.imwrite('eroded_gray.png', eroded_gray)
    df = pd.DataFrame(eroded_gray, index=range(eroded_gray.shape[0]), columns=range(eroded_gray.shape[1]))
    df.to_csv('eroded_tissue.txt', sep='\t')
    print('Step 2: enroded tissue has been done!')
    print('used enrod size = {}'.format(str(kernel_size)))
    return df


def remove_bl(data15,df,num,bg_percent):    #背景信号和边缘模糊信号都去掉
    df[df==255]=1
    mz_tissue = []
    for j in range(2,np.shape(data15)[1]):   #1000个m/z
        tis = data15.iloc[:,[0,1,j]]
        tis.index = range(1,np.shape(tis)[0]+1)
        per = np.percentile(tis.iloc[:,2],num)
        hmat = tis.pivot(index='y', columns='x', values=tis.columns[2]) 
        hmat.index = range(np.shape(hmat)[0])
        hmat.columns = range(np.shape(hmat)[1])
        hmat.iloc[:,0:5]=0                       #先去除可能的污染区域
        count90 = np.sum(np.array(hmat>=per))   #去除亮条之后切片表达值高于90%的元素有4626个，
        hmat[hmat<per]=0                        #保留90%以上的元素位置，总共有4628个位置非0
        hmat = df*hmat                          #这个是保守组织区域
        count_tis = np.sum(np.array(hmat>0))       #寻找90%以上的元素落在组织区域的数量
        if count_tis/count90>bg_percent:        #80%以上最强信号都落在组织区域认为是真实信号
            mz_tissue.append(tis.columns[2])
    data15_signal = data15.loc[:,mz_tissue]
    data15_signal.insert(0, 'x', data15['x'])
    data15_signal.insert(1, 'y', data15['y'])
    data15_bl = data15.drop(columns=mz_tissue)   #开始提取blacklist
    mzbg = list(data15_bl)     #保留bl和模糊地带
    del mzbg[0]
    del mzbg[0]
    mzbg = list(map(float, mzbg))
    mz_tissue = list(map(float, mz_tissue))
    return mzbg, mz_tissue, data15_signal

def preprocessing_signal(Rscores, ST87_20210331, threshold):
    Rscores = Rscores.drop(columns=[' DriftTime'])
    Rscores = Rscores.drop(columns=[' Analyte'])
    Rscores = Rscores.sort_values(by = ' R',ascending=False)
    Rscores['tmp'] = list(map(lambda x:'%1.4f' % x, Rscores['M/z']))   #index会自动舍去0，强制保留格式
    ST87_20210331_signal = ST87_20210331[Rscores[Rscores[' R']>=threshold]['tmp']]
    ST87_20210331_signal.insert(0,'x', ST87_20210331['x'])
    ST87_20210331_signal.insert(1,'y', ST87_20210331['y'])
    ST87_20210331_signal.to_csv('{}_signal.txt'.format(sample), sep='\t', index=False)
    return ST87_20210331_signal


def get_tissue(data15_removebl, df):
    rowcount = np.shape(data15_removebl)[1]
    tis_df = np.zeros((rowcount-2,np.shape(data15_removebl)[0]))
    tis_df = pd.DataFrame(tis_df)
    for j in range(2,rowcount):   #1000个m/z
        tis = data15_removebl.iloc[:,[0,1,j]]
        tis.index = range(1,np.shape(tis)[0]+1)
        hmat = tis.pivot(index='y', columns='x', values=tis.columns[2]) 
        hmat.index = range(np.shape(hmat)[0])
        hmat.columns = range(np.shape(hmat)[1])
        hmat = df*hmat
        tis_df.loc[j-2] = np.array(hmat).ravel()
    tis_df = pd.DataFrame(tis_df)  
    newcname = data15_removebl.columns[2:rowcount]
    tis_df.index = newcname
    #准备列名，zip函数将两个集合对应元素拼接
    xcount = len(set(data15_removebl['x'].values.flatten()))  
    ycount = len(set(data15_removebl['y'].values.flatten()))  
    z = list(zip(list(range(1,xcount+1))*ycount,[val for val in range(1,ycount+1) for i in range(xcount)]))
    colname = [str(a) + "x" + str(b) for a, b in z]
    tis_df.columns = colname  
    tis_df = tis_df.loc[:,~(tis_df==0).all(axis=0)]  
    return tis_df
    
def find_coor_selected(ST87B_20210331un_3k, df_0331):
    rowcount = np.shape(ST87B_20210331un_3k)[1]
    tis_df = np.zeros((rowcount-2,np.shape(ST87B_20210331un_3k)[0]))
    tis_df = pd.DataFrame(tis_df)
    tis = ST87B_20210331un_3k.iloc[:,[0,1,2]]
    tis.index = range(1,np.shape(tis)[0]+1)
    hmat = tis.pivot(index='y', columns='x', values=tis.columns[2]) 
    df_0331.index = hmat.index
    df_0331.columns = hmat.columns
    hmat = df_0331*hmat
    where = np.where(df_0331 == 1)
    y = df_0331.index[np.where(df_0331 == 1)[0]]
    x = df_0331.columns[np.where(df_0331 == 1)[1]]
    coor = {'x':x,'y':y}
    coor = pd.DataFrame(coor)
    return coor
    
def max_intensity(ST06_20211019un_signal):   
    col = ST06_20211019un_signal.iloc[:,[0,1,2]]
    ST06_20211019un_signal = ST06_20211019un_signal[ST06_20211019un_signal.columns[3:][np.max(ST06_20211019un_signal.iloc[:,3:], axis = 0)>=400]]  
    ST06_20211019un_signal = pd.concat([col, ST06_20211019un_signal],axis = 1)
    return ST06_20211019un_signal

    
def mznum_400(ST06_20211019un_signal):
    mz = list(ST06_20211019un_signal.columns[3:])
    mz = list(map(float, mz))
    mz_lip = len([x for x in mz if x>=400])
    mz_meta = len([x for x in mz if x<400])
    return mz_meta, mz_lip
    
    
def get_mz(ST132_20210420):
    coln0420 = list(ST132_20210420)
    del coln0420[0]   
    del coln0420[0]   
    coln0420 = list(map(float, coln0420))
    return coln0420
    
    