import json
import requests
import pandas
from os import remove
from os.path import join
from tempfile import _get_candidate_names, gettempdir

#### FUNCTIONS ####
# def get_or_exception(url_ : str, userListId_ : str):
#     response_ = requests.get(url_ % userListId_)
#     if not response_.ok:
#         raise Exception('Error getting gene list')
#
#     data_ = json.loads(response_.text)
#
#     return  data_
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

response_test = get_or_exception(url_ = enrichr_export, userListId_ = userListId, query_string_ = '?userListId=%s&filename=%s&backgroundType=%s', gene_set_library_ = 'KEGG_2015', filename_ = 'example_enrichment')


# def get_or_exception(url_ : str, userListId_ : str = None, query_string_ : str = None, gene_set_library_ : str = None, gene_ : str = None, filename_ : str = None):
#     if url_ == 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s':
#         response_ = requests.get(url_ % userListId_)
#         if not response_.ok:
#             raise Exception('Error getting gene list')
#
#         data_ = json.loads(response_.text)
#
#     elif url_ == 'http://amp.pharm.mssm.edu/Enrichr/enrich':
#         response_ = requests.get(
#             url_ + query_string_ % (userListId_, gene_set_library_)
#          )
#         if not response_.ok:
#             raise Exception('Error fetching enrichment results')
#
#         data_ = json.loads(response_.text)
#
#     elif url_ == 'http://amp.pharm.mssm.edu/Enrichr/genemap':
#         response_ = requests.get(url_ + query_string_ % gene_)
#         if not response_.ok:
#             raise Exception('Error searching for terms')
#
#         data_ = json.loads(response_.text)
#
#     elif url_ == 'http://amp.pharm.mssm.edu/Enrichr/export':
#         url = url_ + query_string_ % (userListId_, filename_, gene_set_library_)
#         response_ = requests.get(url, stream=True)
#
#         with open(filename_ + '.txt', 'wb') as f:
#             for chunk in response_.iter_content(chunk_size=1024):
#                 if chunk:
#                     f.write(chunk)
#
#         data_ = 'Results exported'
#
#     else:
#         raise Exception('Bad url provided')
#
#     return  data_

#### FUNCTIONS ####


#### METADATA ####
enrichr_addlist = 'http://amp.pharm.mssm.edu/Enrichr/addList'
enrichr_view = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s' #Needs arguments: url_, userListId_
enrichr_enrich = 'http://amp.pharm.mssm.edu/Enrichr/enrich' #Needs arguments: url_, userListId_, query_string_, gene_set_library_
enrichr_genemap = 'http://amp.pharm.mssm.edu/Enrichr/genemap' #Needs arguments:url_, userListId_, query_string_, gene_
enrichr_export = 'http://amp.pharm.mssm.edu/Enrichr/export' #Needs arguments:
#### METADATA ####

#### INPUT ####
# input needs to be in a format 'HNF1A\nPCYT2\nSIGMAR1'
with open('enrichr_input.txt', 'r') as f:
    input = f.read()

description = 'hiroaki'

payload = {
    'list': (None, input),
    'description': (None, description)
}
#### INPUT ####

#### ADD A LIST ####
response_addlist = requests.post(enrichr_addlist, files=payload)
if not response_addlist.ok:
    raise Exception('Error analyzing gene list')

data_response_addlist = json.loads(response_addlist.text)

userListId = data_response_addlist['userListId']
#### ADD A LIST ####

response_test = get_or_exception(url_ = enrichr_view, userListId_ = userListId)

response_test = get_or_exception(url_ = enrichr_enrich, userListId_ = userListId, query_string_ = '?userListId=%s&backgroundType=%s', gene_set_library_ = 'KEGG_2015')

response_test = get_or_exception(url_ = enrichr_genemap, userListId_ = userListId, query_string_ = '?json=true&setup=true&gene=%s', gene_ = 'ALKBH2')

# This writes actuall data table
response_test = get_or_exception(url_ = enrichr_export, userListId_ = userListId, query_string_ = '?userListId=%s&filename=%s&backgroundType=%s', gene_set_library_ = 'KEGG_2015', filename_ = 'example_enrichment')

test = pandas.read_csv('KEGG_2015.txt', sep='\t')

test = pandas.read_csv(response_test.text, sep='\t')
