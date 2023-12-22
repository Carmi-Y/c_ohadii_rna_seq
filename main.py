import os
import glob
import scipy
import pathlib
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import scienceplots

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--input_path', help='The input directory path containing raw text files to be processed', required=True)
    parser.add_argument('-o', '--output_path', help='The output directory', required=True)
    
    args = parser.parse_args()
    input_path = os.path.normpath(args.input_path)
    output_path = os.path.normpath(args.output_path)

    # Grab all txt files in the directory
    files = glob.glob(f"{input_path}/*.xlsx")
    
    organized_files_df = get_files_data_in_dict(files)

    joined_expression_df = join_expression_df_with_GO_terms(organized_files_df['expression_phase_I_vs_phase_II'], organized_files_df['GO_term_keys'])
    
    # export the joined dataframe to a csv file
    joined_expression_df.to_csv(os.path.join(output_path, 'joined_expression_df.csv'), index=False)

    # Create a directory for the plots
    plots_dir = create_directory(output_path, 'plots')

        
        
def get_files_data_in_dict(files):
    '''
    Description
    -----------
    Get the dataframes from each file under a specific key in a dictionary

    Parameters
    ----------
    files : list
        A list of file paths
    
    Returns
    -------
    df_dict : dict
        A dictionary of dataframes with the following keys:

        'GO_term_keys'
        'expression_phase_I_vs_phase_II'
        'plastid_expression_phase_I_vs_phase_II'
        'mito_expression_phase_I_vs_phase_II'

        Each key containing the dataframe of the data from the corresponding file
    '''
    df_dict = {
        'GO_term_keys': '', 'expression_phase_I_vs_phase_II': '',
        'plastid_expression_phase_I_vs_phase_II': '', 'mito_expression_phase_I_vs_phase_II': ''
    }

    for file in files:
        file_name = pathlib.Path(file).stem
        # Check if the file was already read and added to the dictionary
        if df_dict[file_name] != '':
            raise ValueError(f'The file {file_name} was already read, please check the input directory for duplicates')

        # Read the excel file
        df = pd.read_excel(file)
        # Add the dataframe to the dictionary
        df_dict[file_name] = df

    return df_dict


def join_expression_df_with_GO_terms(expression_df, go_term_df):
    '''
    Description
    -----------
    Join the expression dataframe with the GO term dataframe

    Parameters
    ----------
    expression_df : pandas.DataFrame
        The expression dataframe. Must contain a column named 'gene'
    go_term_df : pandas.DataFrame
        The GO term dataframe. Must contain a column named 'gene'

    Returns
    -------
    joined_df : pandas.DataFrame
        exression_df with two extra columns: 'name' and 'description'
        name : The name of the GO term
        description : The description of the GO term
    '''
    # Check if the expression dataframe contains a column named 'gene'
    if 'gene' not in list(expression_df.columns):
        raise ValueError('The expression dataframe must contain a column named "gene"')
    # Check if the GO term dataframe contains a column named 'gene'
    if 'gene' not in list(go_term_df.columns):
        raise ValueError('The GO term dataframe must contain a column named "gene"')
    
    # If expression_df contains a name column, drop it
    if 'name' in expression_df.columns:
        expression_df = expression_df.drop(columns=['name'])

    # Drop the 'significant' column from expression_df
    expression_df = expression_df.drop(columns=['significant'])

    # Drop rows where q_value is NaN or grate than 0.05
    expression_df = expression_df.dropna(subset=['q_value'])
    expression_df = expression_df[expression_df['q_value'] <= 0.05]

    # Drop rows where either T1 or T2 is NaN or 0
    expression_df = expression_df.dropna(subset=['T1', 'T2'])
    expression_df = expression_df[(expression_df['T1'] != 0) & (expression_df['T2'] != 0)]

    # remove the text after '.' in the gene column
    expression_df['gene'] = expression_df['gene'].str.split('.').str[0]
    
    # Join the two dataframes
    expression_df = expression_df.merge(go_term_df, on='gene', how='left')

    # Drop rows where GO is 0
    expression_df = expression_df[expression_df['GO'] != 0]

    # Drop rows where name is '---NA---'
    expression_df = expression_df[expression_df['name'] != '---NA---']

    return expression_df
    

def create_directory(parent_directory, nested_directory_name):
    '''
    Description
    -----------
    Create a directory if it does not exist
    
    Parameters
    ----------
    parent_directory : str
        The path to the directory under which the new directory will be created
    nested_directory_name : str
        The name of the nested directory to be created
    '''
    # Create the output directory path
    new_dir_path = os.path.join(parent_directory, nested_directory_name)
    # Create the directory if it does not exist
    if not os.path.isdir(new_dir_path):
        os.mkdir(new_dir_path)
    return new_dir_path


if __name__ == "__main__":
    main()