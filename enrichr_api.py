import json
import requests
import pandas
from os import remove
from os.path import join
from tempfile import _get_candidate_names, gettempdir



###################
#### FUNCTIONS ####
###################
# After: https://stackoverflow.com/questions/26541416/generate-temporary-file-names-without-creating-actual-file-in-python
def get_tempfile_name(some_id : str = 'TEMP.txt'):
    return join(next(_get_candidate_names()) + some_id)



def get_or_exception(url_ : str, userListId_ : str = None, query_string_ : str = None, gene_set_library_ : str = None, gene_ : str = None, filename_ : str = None):
    if url_ == 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s':
        response_ = requests.get(url_ % userListId_)
        if not response_.ok:
            raise Exception('Error getting gene list')

        data_ = json.loads(response_.text)

    elif url_ == 'http://amp.pharm.mssm.edu/Enrichr/enrich':
        response_ = requests.get(
            url_ + query_string_ % (userListId_, gene_set_library_)
         )
        if not response_.ok:
            raise Exception('Error fetching enrichment results')

        data_ = json.loads(response_.text)

    elif url_ == 'http://amp.pharm.mssm.edu/Enrichr/genemap':
        response_ = requests.get(url_ + query_string_ % gene_)
        if not response_.ok:
            raise Exception('Error searching for terms')

        data_ = json.loads(response_.text)

    elif url_ == 'http://amp.pharm.mssm.edu/Enrichr/export':
        url = url_ + query_string_ % (userListId_, filename_, gene_set_library_)
        response_ = requests.get(url, stream=True)

        temp_filename = f'{gene_set_library_}__{get_tempfile_name()}'

        with open(temp_filename, mode='w') as fp:
             fp.write(response_.text)

        data_ = pandas.read_csv(temp_filename, sep='\t')

        data_['Database'] = gene_set_library_

        remove(temp_filename)

    else:
        raise Exception('Bad url provided')

    return  data_



def create_merged_df(userListId__ : str, databases_to_query_array, url__ : str = 'http://amp.pharm.mssm.edu/Enrichr/export',  query_string__ = '?userListId=%s&filename=%s&backgroundType=%s', ):
    temp1 = pandas.DataFrame({})
    for value in databases_to_query_array:
        if temp1.empty:
            temp1 = get_or_exception(url_ = url__, userListId_ = userListId__, query_string_ = query_string__, gene_set_library_ = value, filename_ = value)
        else:
            temp2 = get_or_exception(url_ = url__, userListId_ = userListId__, query_string_ = query_string__, gene_set_library_ = value, filename_ = value)

            temp1 = temp1.append(temp2, sort=False)

    return temp1



def filter_results(df_, col_overlap_count : str = 'Overlap', col_genes_in_term : str = 'Genes_captured', regex_to_extract_number_from_overlap : str = '(.*)/', genes_in_term_cutoff : int = 2, col_p_val : str = 'P-value', p_val_cutoff = 0.05):
    # Filter based on number of genes in given term
    df_[col_genes_in_term] = df_[col_overlap_count].str.extract(regex_to_extract_number_from_overlap)
    df_[col_genes_in_term] = pandas.to_numeric(df_[col_genes_in_term])
    df_ = df_.loc[df_[col_genes_in_term] > genes_in_term_cutoff]
    df_.pop(col_genes_in_term)

    # Filter based on p-value
    df_[col_p_val] = pandas.to_numeric(df_[col_p_val])
    df_ = df_.loc[df_[col_p_val] < p_val_cutoff]

    return df_

###################
#### FUNCTIONS ####
###################



##################
#### METADATA ####
##################

#### ADRESSES ####
enrichr_addlist = 'http://amp.pharm.mssm.edu/Enrichr/addList'
enrichr_view = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s' #Needs arguments: url_, userListId_
enrichr_enrich = 'http://amp.pharm.mssm.edu/Enrichr/enrich' #Needs arguments: url_, userListId_, query_string_, gene_set_library_
enrichr_genemap = 'http://amp.pharm.mssm.edu/Enrichr/genemap' #Needs arguments:url_, userListId_, query_string_, gene_
enrichr_export = 'http://amp.pharm.mssm.edu/Enrichr/export' #Needs arguments:
#### ADRESSES ####

#### INPUT ####
# input needs to be in a format 'HNF1A\nPCYT2\nSIGMAR1'
with open('enrichr_input_P___YB_YC1.txt', 'r') as f:
    input = f.read()

description = 'silvio_P___YB_YC1'
output_name = 'annotation' + description + '.tsv'

payload = {
    'list': (None, input),
    'description': (None, description)
}
#### INPUT ####

#### DATABASES ####
with open('enrichr_databases.txt', 'r') as f:
    databases = f.read()

databases = databases.split('\n')

# test_databases = {'KEGG_2015', 'Reactome_2016'}
#### DATABASES ####

#### PREPARE A LIST ####
response_addlist = requests.post(enrichr_addlist, files = payload)
if not response_addlist.ok:
    raise Exception('Error analyzing gene list')

data_response_addlist = json.loads(response_addlist.text)

userListId = data_response_addlist['userListId']
#### PREPARE A LIST ####

#### POSSIBLE GETS ####
# response_test = get_or_exception(url_ = enrichr_view, userListId_ = userListId)
# response_test = get_or_exception(url_ = enrichr_enrich, userListId_ = userListId, query_string_ = '?userListId=%s&backgroundType=%s', gene_set_library_ = 'KEGG_2015')
# response_test = get_or_exception(url_ = enrichr_genemap, userListId_ = userListId, query_string_ = '?json=true&setup=true&gene=%s', gene_ = 'ALKBH2')
# response_test = get_or_exception(url_ = enrichr_export, userListId_ = userListId, query_string_ = '?userListId=%s&filename=%s&backgroundType=%s', gene_set_library_ = 'KEGG_2015', filename_ = 'example_enrichment')
#### POSSIBLE GETS ####

##################
#### METADATA ####
##################



#################################
#### CREATE MERGED DATAFRAME ####
#################################
merged_df = create_merged_df(userListId__ = userListId, databases_to_query_array = databases)

#df_ = df_.loc[df_[col_p_val] < p_val_cutoff]
filtered_merged_df = filter_results(df_ = merged_df)

filtered_merged_df.loc[:,'Overlap'] = filtered_merged_df['Overlap'].str.replace(pat = '/', repl = ' out of ')

filtered_merged_df.to_csv(output_name, sep = '\t', decimal = ',', index = False)
#################################
#### CREATE MERGED DATAFRAME ####
#################################
