"""
Load function.py.
Create a new folder named 'csv_output_dir_path'. Then create some folders under it. One sample corresponds to one folder. Each folder is named after the 
corresponding sample name. The Preprocessing.py output files are stored in these folders.
"""

# Find a suitable sieve that matches the HE image.
# Look for your m/z and threshold from the data under mz_from_path. mz_from_path may be the same as data_file_path, depending on where you store the data.
# Two image files gray.png and eroded_gray.png are generated under the mz_from_path. You can check to see if the sieve fits your HE's shape. If yes, proceed to 
# the next step. If no, modify your mz and threshold until satisfied.
def find_interest_sieve(data_file_path, mz_from_path, mz, threshold, sample):
        data = pd.read_csv(data_file_path, sep = '\t')
        os.chdir(mz_from_path) 
        mz = mz
        threshold = threshold
        sample = sample
        hmat = get_raw_tissue(data, mz, threshold)
        df_data = get_max_contour(hmat, 2)
        df_data = erode_tissue(5)   
        return df_data

# Enter df_data. You will get three csv files as input to the filter.mass.R.
def get_output_csv_file(data_file_path, data_file_path_lm, df_data, csv_output_dir_path, sample):
        data = pd.read_csv(data_file_path, sep = '\t')
        data = preprocessing(data)
        data = data.iloc[:,:3002]    #choose the number of m/z
        data_lm = pd.read_csv(data_file_path_lm, sep = '\t')   #lockmass
        data_lm = preprocessing(data_lm)
        data_lm = data_lm.iloc[:,:3002] 

        coln_mzbg, coln_signal, data_signal = remove_bl(data, df_data, 90, 0.5)
        print(len(coln_signal))
        data_signal.insert(0,'Index',range(1,np.shape(data_signal)[0]+1))
        csv_output_file = paste0(csv_output_dir_path, '/', sample,'_3k_signal.unlock_mass.txt')
        data_signal.to_csv(csv_output_file, sep='\t', index = False) 

        coln_mzbg_lm, coln_signal_lm, data_signal_lm = remove_bl(data_lm, df_data, 90, 0.5)
        print(len(coln_signal_lm))
        data_signal_lm.insert(0,'Index',range(1,np.shape(data_signal_lm)[0]+1))
        csv_output_file_lm = paste0(csv_output_dir_path, '/', sample,'_3k_signal.lock_mass.txt')
        data_signal_lm.to_csv(csv_output_file_lm, sep='\t', index = False) 

        coor = find_coor_selected(data_signal_lm, df_data)
        csv_output_coor = paste0(csv_output_dir_path, '/', sample,'_3k_signal.selected.coor.txt')
        coor.to_csv(csv_output_coor, sep='\t', index = False) 

        return data_signal, data_signal_lm, coor









