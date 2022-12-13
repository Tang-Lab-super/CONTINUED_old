

def find_interest_factor(data_file_path, mz_from_path, mz, threshold, sample):
        """
        找出符合要求的筛子因子
        """
        ST87_20210331un_5k = pd.read_csv(data_file_path, sep = '\t')
        os.chdir(mz_from_path)     #mz_from_path路径下必须有数据，从中寻找阈值m/z
        mz = mz
        threshold = threshold
        sample = sample
        hmat_0331 = get_raw_tissue(ST87_20210331un_5k, mz, threshold)
        df_data = get_max_contour(hmat_0331, 2)
        """
        路径mz_from_path下会生成两个图像文件gray.png和erode.gray.png，可以查看筛子是否符合形状，
        若符合，继续往下，若不符合，适当修改mz和threshold
        """
        df_data = erode_tissue(5)   
        return df_data


def get_output_csv_file(data_file_path, data_file_path_lm, df_data, csv_output_dir_path, sample):
        """
        将find_interest_factor找到的df_data作为参数输入，输出三个csv文件，作为下一步**.R的输入
        """
        ST87_20210331un_5k = pd.read_csv(data_file_path, sep = '\t')
        ST87_20210331un_5k = preprocessing(ST87_20210331un_5k)
        ST87_20210331un_3k = ST87_20210331un_5k.iloc[:,:3002]    #取前3000个m/z
        ST87_20210331_5k = pd.read_csv(data_file_path_lm, sep = '\t')   #lockmass
        ST87_20210331_5k = preprocessing(ST87_20210331_5k)
        ST87_20210331_3k = ST87_20210331_5k.iloc[:,:3002] 

        coln0331_mzbg, coln0331_signal, ST87_20210331un_signal = remove_bl(ST87_20210331un_3k, df_data, 90, 0.5)
        print(len(coln0331_signal))
        ST87_20210331un_signal.insert(0,'Index',range(1,np.shape(ST87_20210331un_signal)[0]+1))
        csv_output_file = paste0(csv_output_dir_path, '/', sample,'_3k_signal.unlock_mass.txt')
        ST87_20210331un_signal.to_csv(csv_output_file, sep='\t', index = False) 

        coln0331_mzbg, coln0331un_signal, ST87_20210331_signal = remove_bl(ST87_20210331_3k, df_data, 90, 0.5)
        print(len(coln0331un_signal))
        ST87_20210331_signal.insert(0,'Index',range(1,np.shape(ST87_20210331_signal)[0]+1))
        csv_output_file_lm = paste0(csv_output_dir_path, '/', sample,'_3k_signal.lock_mass.txt')
        ST87_20210331_signal.to_csv(csv_output_file_lm, sep='\t', index = False) 
        #找坐标
        coor = find_coor_selected(ST87_20210331_signal, df_data)
        csv_output_coor = paste0(csv_output_dir_path, '/', sample,'_3k_signal.selected.coor.txt')
        coor.to_csv(csv_output_coor, sep='\t', index = False) 

        return ST87_20210331un_signal, ST87_20210331_signal, coor









