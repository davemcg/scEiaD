# requires --gres=gpu:v100x:1, and module load CUDA/10.1 && module load cuDNN/7.6.5/CUDA-10.1 
# %%
import pandas as pd 
import numpy as np
import os 
from sklearn.model_selection import train_test_split, cross_val_predict
from xgboost import XGBClassifier
from sklearn.metrics import classification_report
import time 
import argparse
import warnings
import pickle
warnings.filterwarnings("ignore")

class SklDataObj:
    def __init__(self, X_data, feature_cols,  lab2id,name=None, sampling_method = None, k=None):
        self.label2id = lab2id
        if sampling_method is None:
            X_arr = X_data.filter(feature_cols).to_numpy()# first test on scvi only 
            Y_arr = X_data['cell_type_id'].to_numpy()
            X_train, X_test, Y_train, Y_test =train_test_split(X_arr,Y_arr,
                                                             test_size=.22222, random_state=42,stratify=Y_arr)

        elif sampling_method == 'downsample': 
            X_data_ds, removed_data = self.down_sampler(X_data, k)
            X_arr = X_data_ds.filter(feature_cols).to_numpy()# first test on scvi only 
            Y_arr = X_data_ds['cell_type_id'].to_numpy()
            remove_x_arr = removed_data.filter(feature_cols).to_numpy()
            remove_y_arr = removed_data['cell_type_id'].to_numpy()
            X_train, X_tmp, Y_train, Y_tmp =train_test_split(X_arr,Y_arr,
                                                             test_size=.22222, random_state=42,stratify=Y_arr)
            X_test = np.concatenate([X_tmp, remove_x_arr])
            Y_test  = np.concatenate([Y_tmp, remove_y_arr])

        
        self.X_train=X_train
        self.Y_train=Y_train
        self.X_test=X_test
        self.Y_test=Y_test
        self.name=name
        self.model = None
    def down_sampler(self, df,k):
        # down sample all celltypes that have total count above k
        cell_counts = df.CellType.value_counts()
        cells_to_ds = list(cell_counts[ cell_counts > k ].index)
        downsampled_data = df[df['CellType'].isin(cells_to_ds)].groupby('CellType').apply(pd.DataFrame.sample, n=k)
        clean_data_downsampled = pd.concat( [df[~df['CellType'].isin(cells_to_ds)], downsampled_data]) 
        removed_samples = df[~df.Barcode.isin(clean_data_downsampled['Barcode'])]
        return clean_data_downsampled, removed_samples
    
    def train(self, model, model_name=None):
        self.model = model
        start = time.time()
        self.model.fit(self.X_train, self.Y_train)
        end = time.time()
        diff = end - start
        print(f'\nTraining Time: {int(diff/60)} min {diff%60} seconds\n')
        self.model_name = model_name
    def test(self, test_ext_data=False, X=None, Y=None):
        if test_ext_data:
            self.X_test = X
            self.Y_test = Y
        start = time.time()
        self.Y_test_predicted =  self.model.predict(self.X_test)
        end = time.time()
        diff = end - start
        print(f'Test Time: {int(diff/60)} min {diff%60} seconds\n')
        
        cr = classification_report(self.Y_test, self.Y_test_predicted, output_dict = True)
        i="accuracy"
        print(f'Total Accuracy: {round(cr[i], 4)}')
        for i in ['macro avg', 'weighted avg']:
            ol=f'{i}:\t'
            for j in ['precision', 'recall', 'f1-score', 'support']:
                ol+= f'{j}: {round(cr[i][j], 4)}\t'
            print(ol)
            del cr[i]
        print('\n')
        num_celltype = self.label2id.shape[0]
        ids = [str(i) for i in range(num_celltype)]
        class_df = pd.DataFrame([ cr[i] for i in ids  ])
        class_df['cell_type_id'] = pd.Series(ids).astype('int64')
        class_df = class_df.merge(self.label2id, on = 'cell_type_id')
        self.class_rep_df = class_df.sort_values(by= 'f1-score')
        print(class_df.sort_values(by= 'f1-score').to_string() )
    def predict(X):
        start = time.time()
        pred_probs = self.model.predict_proba(X)
        preds = np.asarray([np.argmax(i) for i in pred_probs] )
        end = time.time()
        diff = end - start
        print(f'Prediction Time: {int(diff/60)} min {diff%60} seconds\n')
        return preds, pred_probs

def encode_age(stld):
    stld['Age'].fillna(10000, inplace = True)
    conditions = [  (stld['Age'] < -165) & (stld['organism'] ==   'Homo sapiens')  ,
                   (stld['Age'] >= -165) & (stld['Age'] < 0 ) & (stld['organism'] == 'Homo sapiens')  ,
                    #(stld['Age'] >= 0) & (stld['organism'] ==   'Homo sapiens')  ,
                    (stld['Age'] <= -3) & (stld['organism'] =='Mus musculus'  ) ,
                    (stld['Age'] > -3) & (stld['Age'] <=8) & (stld['organism'] ==   'Mus musculus'  ),
                  #(stld['Age'] > 8) & (stld['organism'] ==   'Mus musculus')  
                 ]
    values = ['early_dev',
                'late_dev',
             'early_dev', 
             'late_dev']
    stld.loc[:,'age_group'] = np.select(conditions, values, default = 'adult')
    stld = stld.reset_index(drop = True)
    categorical_columns= [ 'age_group']
    cat_col_df = pd.get_dummies(stld[categorical_columns])
    stld = stld.drop(columns = categorical_columns).join(cat_col_df)
    return stld

parser = argparse.ArgumentParser(description = 'train cell type predictor or predict cell types')
parser.add_argument('mode',action = 'store', choices = ['train', 'predict'])
parser.add_argument('--workingDir', action = 'store', default = None, help = 'change directory to here (must be an absolute path)')
parser.add_argument('--inputMatrix', action = 'store', default=None, help = 'input .tsv matrix to train / predict on ')
parser.add_argument('--featureCols', action = 'store', default = None, help = 'Optional file with newline seperated feature names to train on')
parser.add_argument('--trainedModelFile', action = 'store', default = 'cell_type_predictor.xgb', help = 'file to save or load trained model')
parser.add_argument('--predictions', action = 'store', default='predictions.tsv', help = 'output .tsv file to save predictions')
args = parser.parse_args()
if args.workingDir is not None:
    os.chdir(args.workingDir)

non_feature_cols = ['cluster', 'batch', 'cluster', 'subcluster', 'sample_accession', 'library_layout', 'Platform', 'UMI', 'Covariate', 'integration_group', 
                        'Paper', 'biosample_title', 'sample', 'barcode', 'donor', 'region', 'Method', 
                      'study_accession', 'SubCellType','TissueNote', 'orig.ident', 'Tissue', 'organism', 'Age']

bad_cell_types = ["RPE/Margin/Periocular Mesenchyme/Lens Epithelial Cells", "Droplet", "Droplets", 'Doublet', 'Doublets', 'Mast']
def train(args, non_feature_cols, bad_cell_types):
    all_data = pd.read_csv(args.inputMatrix, sep = '\t')
    is_bad_celltype = (all_data.CellType.isin(bad_cell_types)) |(all_data.CellType.isnull())
    print(f'Removing {sum(is_bad_celltype)} cells with missing or excluded cell types\n')
    cell_types =  list(all_data[~is_bad_celltype].CellType.value_counts().index)
    cell_type2id =pd.DataFrame({ 'CellType': cell_types , 'cell_type_id': list(range(len(cell_types))) } )
    clean_labeled_data = encode_age(all_data[~is_bad_celltype])
    clean_labeled_data = clean_labeled_data.merge(cell_type2id, on = 'CellType', how = 'left')
    id_cols = ['Barcode', 'CellType', 'cell_type_id']

    if args.featureCols is not None:
        with open(args.featureCols) as infile:
            feature_cols = [line.strip('\n') for line in infile]
    else:
        feature_cols = list(clean_labeled_data.drop(columns = non_feature_cols + id_cols).columns)
    print('Training on Following Features:')
    print(feature_cols)
    data_obj = SklDataObj(clean_labeled_data, feature_cols, cell_type2id)
    xgbc_gpu = XGBClassifier(tree_method = 'gpu_hist', gpu_id = 0)
    data_obj.train(xgbc_gpu)
    data_obj.test()
    trained_model = data_obj.model
    with open(args.trainedModelFile, 'wb+') as modelfile:
        pickle.dump((trained_model, feature_cols, cell_type2id ),modelfile )
    
def predict(args, non_feature_cols, bad_cell_types):
    print('Loading Data...\n')
    with open(args.trainedModelFile, 'rb') as modelfile:
        model_info= pickle.load(modelfile )
        trained_model = model_info[0]
        feature_cols = model_info[1]
        cell_type2id = model_info[2]
        cell_type2id = cell_type2id.sort_values('cell_type_id')
        all_data = pd.read_csv(args.inputMatrix, sep = '\t').pipe(encode_age)
    barcodes = all_data.loc[:,'Barcode']
    if args.featureCols is not None:
        with open(args.featureCols) as infile:
            feature_cols = [line.strip('\n') for line in infile]
    X = all_data.filter(feature_cols).to_numpy()
    print('Predicting Data...\n')
    pred_probs = trained_model.predict_proba(X)
    preds_class = np.asarray([np.argmax(i) for i in pred_probs] )
    pred_probs_df = pd.DataFrame(pred_probs, columns = cell_type2id['CellType'])
    pred_probs_df['Barcode'] = barcodes
    full_pred_df = pd.DataFrame({'Barcode':barcodes,  'cell_type_id' : preds_class }).merge(cell_type2id).merge(pred_probs_df)
    full_pred_df.to_csv(args.predictions, sep = '\t', index = False)

if args.mode == 'train':
    train(args,  non_feature_cols, bad_cell_types)
elif args.mode == 'predict':
    predict(args,  non_feature_cols, bad_cell_types)

    