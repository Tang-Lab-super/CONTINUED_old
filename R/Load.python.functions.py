import pandas as pd
import cv2 as cv
import cv2
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def preprocessing(rawdata):
    """
    param rawdata: Raw dataFrame with the first and last columns are NaNs, the first row are NaNs. 
    The first three columns have 'Unnamed' column names.
    
    output predata: The processed dataFrame with no NaNs. The first and second column names are 'x' and 'y'.
    The remaining column names represent m/z values.
    """
    predata = rawdata.drop([0])
    predata = predata.drop(columns='Unnamed: 0')
    newcname = predata.columns[3:]  
    last_column = predata.columns[-1]
    predata = predata.drop(columns=last_column)
    newcname = newcname.insert(0, 'y')
    newcname = newcname.insert(0, 'x')
    predata.columns = newcname
    predata.index = range(predata.shape[0])
    return predata

def get_raw_tissue(predata, mz, threshold):
    """
    param predata: The ouput of 'preprocessing'.
    param mz: Choose a m/z that has a high signal throughout the tissue area and a low background area signal.
    param threshold: A threshold is selected between the lowest signal in the tissue region and the highest
    signal in the background region. The m/z and threshold were selected according to the shape of the HE staining slice.
    
    output hmat: A matrix with the same dimensions as the desi scan area. The elements are just 0 and 1000. 0 represents
    background and 1000 represents the tissue area.
    """
    tis = predata.loc[:,('x','y',mz)]
    tis.loc[tis[mz]<=threshold, mz]=0    
    tis.loc[tis[mz]>threshold, mz]=1000
    hmat = tis.pivot(index='y', columns='x', values=mz)
    hmat.index = range(np.shape(hmat)[0])
    hmat.columns = range(np.shape(hmat)[1])    
    hmat.to_csv('tis{}.txt'.format(threshold), sep='\t')
    return hmat

    
def get_max_contour(hmat, dilate_size):
    """
    param hmat: The ouput of 'get_raw_tissue'.
    param dilate_size: Default is 2.
    
    output: 'gray.png' and 'raw_tissue.txt' will generated in the setting file path. Check 'gray.png' to make 
    sure you get a satisfactory shape for the tissue area.
    """
    thresh = hmat.values
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
    print('Step 1: binarized tissue dilated and eroded has been done!')
    print('used dilate size = {}'.format(str(dilate_size)))
    return df


def erode_tissue(kernel_size):
    """
    param kernel_size: Default is 5. The larger the value, the more inward the edge of the selected area shrinks.

    output: 'eroded_gray.png' and 'eroded_tissue.txt' will generated in the setting file path. Check 'eroded_gray.png' to make 
    sure you get a satisfactory shape for the tissue area.
    output df: A matrix with 0 and 255. 0 represents background and 1000 represents the tissue area.
    """
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


def remove_bl(predata,df,num,bg_percent):    
    """
    param predata: The ouput of 'preprocessing'.
    param df: The output of 'erode_tissue'.
    param num: Range from 0 to 100. The default is 90, which means the coordinate points in the whole 
    scan area whose signal is higher than 90% are used as the judgment object.  
    param bg_percent: Range from 0 to 1. The default is 0.5. Calculate the ratio between the tissue 
    region and the background region where the signal is higher than 90% in the whole scanning region. 
    If the ratio is higher than 0.5, the m/z is considered to be A true metabolite. The m/z is stored in 
    'mz_tissue'. If the ratio is less than 0.5, the m/z is a pollutant and is stored in 'mzbg'.
    
    output mzbg: m/z of pollutants.
    output mz_tissue: m/z of metabolites.
    output predata_signal: Compared with 'predata', 'predata_signal' removes contaminants.
    """
    df[df==255]=1
    mz_tissue = []
    for j in range(2,np.shape(predata)[1]):   
        tis = predata.iloc[:,[0,1,j]]
        tis.index = range(1,np.shape(tis)[0]+1)
        per = np.percentile(tis.iloc[:,2],num)
        hmat = tis.pivot(index='y', columns='x', values=tis.columns[2]) 
        hmat.index = range(np.shape(hmat)[0])
        hmat.columns = range(np.shape(hmat)[1])
        hmat.iloc[:,0:5]=0                      
        count90 = np.sum(np.array(hmat>=per))   
        hmat[hmat<per]=0                        
        hmat = df*hmat                         
        count_tis = np.sum(np.array(hmat>0))      
        if count_tis/count90>bg_percent:       
            mz_tissue.append(tis.columns[2])
    predata_signal = predata.loc[:,mz_tissue]
    predata_signal.insert(0, 'x', predata['x'])
    predata_signal.insert(1, 'y', predata['y'])
    predata_bl = predata.drop(columns=mz_tissue)  
    mzbg = list(predata_bl)     
    del mzbg[0]
    del mzbg[0]
    mzbg = list(map(float, mzbg))
    mz_tissue = list(map(float, mz_tissue))
    return mzbg, mz_tissue, predata_signal


    
def find_coor_selected(predata_signal, df):
    """
    param predata_signal: The output of 'remove_bl'.
    param df: The output of 'erode_tissue'.
    
    output coor: The x,y coordinates of the tissue region.
    """
    
    rowcount = np.shape(predata_signal)[1]
    tis_df = np.zeros((rowcount-2,np.shape(predata_signal)[0]))
    tis_df = pd.DataFrame(tis_df)
    tis = predata_signal.iloc[:,[0,1,2]]
    tis.index = range(1,np.shape(tis)[0]+1)
    hmat = tis.pivot(index='y', columns='x', values=tis.columns[2]) 
    df.index = hmat.index
    df.columns = hmat.columns
    hmat = df*hmat
    where = np.where(df == 1)
    y = df.index[np.where(df == 1)[0]]
    x = df.columns[np.where(df == 1)[1]]
    coor = {'x':x,'y':y}
    coor = pd.DataFrame(coor)
    return coor
    


def find_interest_factor(predata, mz_from_path, mz, threshold):
    """
    param predata: The ouput of 'preprocessing'.
    param mz_from_path: File path from which parameter 'mz' comes.
    param mz: Choose a m/z that has a high signal throughout the tissue area and a low background area 
    signal. 
    param threshold: A threshold is selected between the lowest signal in the tissue region and the 
    highest signal in the background region. The m/z and threshold were selected according to the 
    shape of the HE staining slice.
    
    output: You will get a dataframe 'df_data' and two .png files under mz_from_path, the images can help
    you adjust the parameters.     
    """
    os.chdir(mz_from_path)     
    mz = mz
    threshold = threshold
    hmat_0331 = get_raw_tissue(predata, mz, threshold)
    df_data = get_max_contour(hmat_0331, 2)
    df_data = erode_tissue(5)   
    return df_data


def get_output_csv_file(predata, predata_lm, df_data, csv_output_dir_path, sample):
    """
    param predata: The ouput of 'preprocessing'. Rawdata is generated under parameter 'unlockmass'.
    param predata_lm: The ouput of 'preprocessing'. Rawdata is generated under parameter 'lockmass'.
    param df_data: The ouput of 'find_interest_factor'.
    param csv_output_dir_path: Set any dir path you want to store the results.
    param sample: Set a name that distinguishes different tissues.
    
    output: Three .txt files.
    output predata_lm_signal: Compared with 'predata_lm', 'predata_lm_signal' removes contaminants.
    output predata_signal: Compared with 'predata', 'predata_signal' removes contaminants.
    output coor: The x,y coordinates of the tissue region.
    """
    mzbg, signal, predata_signal = remove_bl(predata, df_data, 90, 0.5)
    print(len(signal))
    predata_signal.insert(0,'Index',range(1,np.shape(predata_signal)[0]+1))
    csv_output_file = csv_output_dir_path + '/' + sample + '_signal.unlock_mass.txt'
    predata_signal.to_csv(csv_output_file, sep='\t', index = False) 
    
    mzbg, signal_lm, predata_lm_signal = remove_bl(predata_lm, df_data, 90, 0.5)
    print(len(signal_lm))
    predata_lm_signal.insert(0,'Index',range(1,np.shape(predata_lm_signal)[0]+1))
    csv_output_file_lm = csv_output_dir_path + '/' + sample + '_signal.lock_mass.txt'
    predata_lm_signal.to_csv(csv_output_file_lm, sep='\t', index = False) 
    
    coor = find_coor_selected(predata_signal, df_data)
    csv_output_coor = csv_output_dir_path + '/' + sample + '_signal.selected.coor.txt'
    coor.to_csv(csv_output_coor, sep='\t', index = False) 

    return predata_signal, predata_lm_signal, coor


    
