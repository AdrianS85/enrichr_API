import json
import requests
import pandas

#### FUNCTIONS ####
def get_or_exception(url_ : str, userListId_ : str):
    response_ = requests.get(url_ % userListId_)
    if not response_.ok:
        raise Exception('Error getting gene list')

    data_ = json.loads(response_.text)
    
    return  data_


# def get_or_exception(url_ : str, userListId_ : str):
#     response_ = requests.get(url_ % userListId_)
#     if not response_.ok:
#         raise Exception('Error getting gene list')
# 
#     data_ = json.loads(response_.text)
# 
#     return  data_
    
#### FUNCTIONS ####


#### METADATA ####
enrichr_addlist = 'http://amp.pharm.mssm.edu/Enrichr/addList'
enrichr_view = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'
enrichr_enrich = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
enrichr_genemap = 'http://amp.pharm.mssm.edu/Enrichr/genemap'
enrichr_export = 'http://amp.pharm.mssm.edu/Enrichr/export'
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




def get_or_exception(url_ : str, userListId_ : str, ):
    response_ = requests.get(url_ % userListId_)
    if not response_.ok:
        raise Exception('Error getting gene list')

    data_ = json.loads(response_.text)
    
    return  data_


ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = userListId
gene_set_library = 'KEGG_2015'
response = requests.get(
    ENRICHR_URL + query_string % (user_list_id, gene_set_library)
 )
if not response.ok:
    raise Exception('Error fetching enrichment results')



data = json.loads(response.text)

app_json = json.dumps(data)

test = pandas.DataFrame.from_dict(data)
test2 = pandas.read_json(response.text, orient='table')


print(data)
type(response.text)


with open('test.json', 'w') as f:
  json.dump(response.text, f)


ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/genemap'
query_string = '?json=true&setup=true&gene=%s'
gene = 'AKT1'
response = requests.get(ENRICHR_URL + query_string % gene)
if not response.ok:
    raise Exception('Error searching for terms')

data = json.loads(response.text)
print(data)




ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
query_string = '?userListId=%s&filename=%s&backgroundType=%s'
user_list_id = 363320
filename = 'example_enrichment'
gene_set_library = 'KEGG_2015'

url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
response = requests.get(url, stream=True)

with open(filename + '.txt', 'wb') as f:
    for chunk in response.iter_content(chunk_size=1024):
        if chunk:
            f.write(chunk)
