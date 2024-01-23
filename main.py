import os
import glob
import scipy
import pathlib
import argparse
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
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


    GO_abundace_summary_df = get_GO_data(joined_expression_df)

    # Save the dataframe to a csv file
    GO_abundace_summary_df.to_csv(os.path.join(output_path, 'summary_go_abundances.csv'), index=False)
        

    # Create a directory for the plots
    plots_dir = create_directory(output_path, 'plots')

    make_GO_plot('C', ['chloroplast envelope', 'chloroplast inner membrane', 'chloroplast stroma', 'chloroplast thylakoid',
                       'photosystem I', 'thylakoid lumen', 'thylakoid membrane'], GO_abundace_summary_df, plots_dir,
                       'Chloroplast', sns.light_palette("seagreen", as_cmap=True))
    

    make_GO_plot('C', ['mitochondrial inner membrane', 'mitochondrial respiratory chain complex I', 'mitochondrial ribosome'],
                    GO_abundace_summary_df, plots_dir, 'Mitochondrion', 'Reds')
    
    make_GO_plot('C', ['endosome', 'extracellular exosome', 'extrinsic component of membrane','extracellular region'],
                 GO_abundace_summary_df, plots_dir, 'Endo and Exo cytosis', 'Blues')
    
    make_GO_plot('C', ['endomembrane system', 'endoplasmic reticulum', 'endoplasmic reticulum lumen', 
                       'endoplasmic reticulum membrane', 'cytoskeleton', 'cell cortex', 'myosin complex',
                        'nuclear envelope', 'nucleoplasm' 'nucleus', 'vacuole'], GO_abundace_summary_df, plots_dir, 'Endomembrane systems', 'Purples')
        
    make_GO_plot('C', ['anaphase-promoting complex', 'condensin complex', 'apoplast', 'BRCA1-A complex', 
                       'cell wall', 'MCM complex', 'microtubule', 'phragmoplast', 'spindle microtubule', 'U7 snRNP', 'nucleolus'], GO_abundace_summary_df, plots_dir,
                       'Cell cycle and DNA repair', sns.color_palette("flare", as_cmap=True))
    
    make_GO_plot('C', ['DNA-directed RNA polymerase II, core complex', 'Elongator holoenzyme complex',
                       'mRNA cleavage and polyadenylation specificity factor complex', 'signal peptidase complex', 'ribosome',
                       'small-subunit processome', 'transcription factor TFIID complex'], GO_abundace_summary_df, plots_dir,
                       'Transcription and translation', sns.diverging_palette(220, 20, as_cmap=True))

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


def get_GO_data(joined_expression_df):
    '''
    Description
    -----------
    Get the GO data from the joined dataframe

    Parameters
    ----------
    joined_expression_df : pandas.DataFrame

    Returns
    -------
    df : pandas.DataFrame
        A dataframe containing the GO data with the following columns:
        phase : The phase of the gene expression
        domain : The domain of the GO term (C, P, F)
        description : The description of the GO term
        abundance : The abundance of the GO term - sum of all reads (in rpkm) that fell under the GO term
    '''
    # In a GO string, the first letter represents the ontology (C, P, F)
    # C - cellular component
    cellular_component_phase_I_keys_and_proportion = {}
    cellular_component_phase_II_keys_and_proportion = {}
    # P - biological process
    biological_process_phase_I_keys_and_proportion = {}
    biological_process_phase_II_keys_and_proportion = {}
    # F - molecular function
    molecular_function_phase_I_keys_and_proportion = {}
    molecular_function_phase_II_keys_and_proportion = {}

    # Make the abundance calculation using dictionaries for efficiency
    for index ,row in joined_expression_df.iterrows():
        expression_in_phae_I = row['T1']
        expression_in_phae_II = row['T2']

        for go_term_inner_item in row['GO'].split(";"):
            go_term_inner_item = go_term_inner_item.strip().split(":")
            go_term_inner_item_key = go_term_inner_item[0]
            go_term_inner_item_value = go_term_inner_item[1]

            if go_term_inner_item_key == 'C':
                if go_term_inner_item_value not in cellular_component_phase_I_keys_and_proportion.keys():
                    cellular_component_phase_I_keys_and_proportion[go_term_inner_item_value] = expression_in_phae_I
                    cellular_component_phase_II_keys_and_proportion[go_term_inner_item_value] = expression_in_phae_II
                else:
                    cellular_component_phase_I_keys_and_proportion[go_term_inner_item_value] += expression_in_phae_I
                    cellular_component_phase_II_keys_and_proportion[go_term_inner_item_value] += expression_in_phae_II
            elif go_term_inner_item_key == 'P':
                if go_term_inner_item_value not in biological_process_phase_I_keys_and_proportion.keys():
                    biological_process_phase_I_keys_and_proportion[go_term_inner_item_value] = expression_in_phae_I
                    biological_process_phase_II_keys_and_proportion[go_term_inner_item_value] = expression_in_phae_II
                else:
                    biological_process_phase_I_keys_and_proportion[go_term_inner_item_value] += expression_in_phae_I
                    biological_process_phase_II_keys_and_proportion[go_term_inner_item_value] += expression_in_phae_II
            elif go_term_inner_item_key == 'F':
                if go_term_inner_item_value not in molecular_function_phase_I_keys_and_proportion.keys():
                    molecular_function_phase_I_keys_and_proportion[go_term_inner_item_value] = expression_in_phae_I
                    molecular_function_phase_II_keys_and_proportion[go_term_inner_item_value] = expression_in_phae_II
                else:
                    molecular_function_phase_I_keys_and_proportion[go_term_inner_item_value] += expression_in_phae_I
                    molecular_function_phase_II_keys_and_proportion[go_term_inner_item_value] += expression_in_phae_II
            else:
                raise ValueError(f'Unexpected GO term key in: {row["GO"]} at row {index}')
            
    
    # Put the contents of the dictionaries in a dataframe
    phases = []
    domains = []
    descriptions = []
    abundances = []

    domains_letters = ['C', 'C', 'P', 'P', 'F', 'F']
    phases_letters = ['I', 'II', 'I', 'II', 'I', 'II']
    for index, dict in enumerate([
                    cellular_component_phase_I_keys_and_proportion, cellular_component_phase_II_keys_and_proportion,
                    biological_process_phase_I_keys_and_proportion, biological_process_phase_II_keys_and_proportion,
                    molecular_function_phase_I_keys_and_proportion, molecular_function_phase_II_keys_and_proportion,
                ]):
        curr_domain = domains_letters[index]
        curr_phase = phases_letters[index]
        for key, value in dict.items():
            phases.append(curr_phase)
            domains.append(curr_domain)
            descriptions.append(key)
            abundances.append(value)


    # Create a dataframe from the lists
    return pd.DataFrame({
        'phase': phases,
        'domain': domains,
        'description': descriptions,
        'abundance': abundances
    })


def make_GO_plot(domain, descriptions, GO_abundace_summary_df, plots_dir, save_name, colormap_for_plot):
    '''
    Description
    -----------
    Make a plot of the GO terms
    
    Parameters
    ----------
    domain : str
        The domain of the GO term (C, P, F)
    descriptions : list of str
        A list of GO term descriptions
    GO_abundace_summary_df : pandas.DataFrame
        A dataframe containing the GO data with the following columns:
        phase : The phase of the gene expression
        domain : The domain of the GO term (C, P, F)
        description : The description of the GO term
        abundance : The abundance of the GO term - sum of all reads (in rpkm) that fell under the GO term
    plots_dir : str
        The path to the plots directory
    
    Returns
    -------
    None
    '''
    domain_data = GO_abundace_summary_df.loc[(GO_abundace_summary_df['domain'] == domain)]
    phase_I_data = domain_data.loc[((domain_data['phase'] == 'I') & (domain_data['description'].isin(descriptions)))]
    phase_II_data = domain_data.loc[((domain_data['phase'] == 'II') & (domain_data['description'].isin(descriptions)))]

    joined_df = phase_I_data.merge(phase_II_data, on='description', how='left', suffixes=('_I', '_II'))
    joined_df = joined_df.drop(columns=['domain_I', 'domain_II', 'phase_I', 'phase_II'])
    joined_df = joined_df.rename(columns={'description': 'description', 'abundance_I': 'abundance_phase_I', 'abundance_II': 'abundance_phase_II'}) 
    joined_df['log2_ratio_phaseII_over_phaseI'] = np.log2(joined_df['abundance_phase_II'] / joined_df['abundance_phase_I'])
    joined_df = joined_df.sort_values(by=['log2_ratio_phaseII_over_phaseI'], ascending=False)


    # Plot the data
    colormap = colormap_for_plot

    # Normalize values to map to the colormap
    norm = plt.Normalize(joined_df['log2_ratio_phaseII_over_phaseI'].min(), joined_df['log2_ratio_phaseII_over_phaseI'].max())

    # Create a colorbar scalar map
    colors = plt.cm.ScalarMappable(norm=norm, cmap=colormap)

    fig, ax = plt.subplots(figsize=(10, 10))
    description_for_sidp =  joined_df['description'].str.replace(' ', '\n')
    bars = ax.bar(description_for_sidp, joined_df['log2_ratio_phaseII_over_phaseI'], color=colors.to_rgba(joined_df['log2_ratio_phaseII_over_phaseI']), edgecolor='black', linewidth=1)

    ax.axline((0, 0), slope=0, color='black', linewidth=1)

    ax.set_title(r'$log_2$' + r'$(\frac{T2}{T1})$' + f' for {save_name} GO terms', fontsize=20)
    ax.set_xlabel('GO term', fontsize=16)
    ax.set_ylabel(r'$log_2$' + r'$(\frac{T2}{T1})$', fontsize=16)

    plt.tight_layout()
    
    plt.savefig(os.path.join(plots_dir, f'log2_ratio_of_phase_II_over_phase_I_for_{save_name}_GO_terms.png'))

if __name__ == "__main__":
    main()