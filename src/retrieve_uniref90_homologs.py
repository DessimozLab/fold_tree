# importing libraries
import json, pandas as pd,re,time,zlib,requests,os
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry
import requests, tqdm, dill
from collections import deque
import numpy as np
from Bio import PDB

# define auxiliary functions
def getting_protein_ids(id_string):
    if '(' in id_string:
        return id_string.rpartition(' ')[2].replace('(', '').replace(')', '')
    else:
        return id_string

def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise
        
def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]

def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
        
def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])
        
def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)
        
def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results

def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]

def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text

def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""

def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")
    
def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results

def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

def split_list(input_list, chunk_size):
    # Create a deque object from the input list
    deque_obj = deque(input_list)
    # While the deque object is not empty
    chunks = []
    while deque_obj:
        # Pop chunk_size elements from the left side of the deque object
        # and append them to the chunk list
        chunk = []
        for _ in range(chunk_size):
            if deque_obj:
                chunk.append(deque_obj.popleft())
        # Yield the chunk
        chunks.append(chunk)
    return chunks

def get_uniref90_cluster_id(ID_list):
    try:
        if len(ID_list) > 1:
            id_request = ''
            for accession in ID_list:
                id_request = id_request + 'uniprot_id:' + accession + '%20OR%20'
            # getting rid of last separator
            id_request = id_request.rpartition('%20OR%20')[0]
            requestURL = "https://rest.uniprot.org/uniref/search?query={0}&format=json&size=500".format(id_request)
        else:
            requestURL = "https://rest.uniprot.org/uniref/search?query=uniprot_id:{0}&format=json".format(ID_list[0])
        #print(requestURL)
        r = requests.get(requestURL, headers={ "Accept" : "application/json"})
        if not r.ok:
          r.raise_for_status()
          sys.exit()
        results_json = r.json()
        return [result['id'] for result in results_json['results'] if result['entryType'] == 'UniRef90']
    except:
        return {'results': {'primaryAccession': None}}

def splitting_dots_and_commas(x):
    if '.' in x:
        return x.rpartition('.')[0]
    elif ',' in x:
        return x.rpartition(',')[0]
    else:
        return x

def get_uniprotkbid_uniref90(cluster_ids):
    try:
        # mapping from UniRef90 to UniRef90 in order to get members and bypass the 500 seqs restriction for API requests
        # if cluster_ids has more than 400 clusters, then group
        if len (cluster_ids) >=  400:
            cluster_id_groups = split_list(input_list = cluster_ids, chunk_size = 400)
        else:
            cluster_id_groups = [cluster_ids]
        # create a dictionary to allocate results
        res_dictionary = {'results': {}}
        # retrieve results
        for cluster_group in cluster_id_groups:
            job_id = submit_id_mapping(from_db='UniRef90', to_db="UniRef90", ids=cluster_group)
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                results = get_id_mapping_results_search(link) # the name of the group is in results['to']['id']
                # getting a dictionary of the form {id: [members], ..., id:[members]}
                groups2members_dict = {result_dictionary['to']['id']: result_dictionary['to']['members'] for result_dictionary in results['results']}
                # now: get into each group, get rid of deprecated IDs
                # then update the dictionary to get that list
                # then update the results dict to include each of this polished dictionaries
                # the following part is done better if its done in a batch first, creating a dictionary that allocates either the good ID or a None
                # after getting this in batch, groups can be polished pretty easly. This is better than doing many requests
                # get the whole set of IDs
                whole_set_ids = [ID for group_name,id_list in groups2members_dict.items() for ID in id_list]
                # now perform the conversion in batch
                # create empty dict
                memberids2uniprotkb = {}
                # split into chunks
                member_ids_chunks = split_list(input_list = whole_set_ids, chunk_size = 400) # spliting into chunks
                for member_ids_group in member_ids_chunks:
                    # submit batch
                    job_id = submit_id_mapping(from_db='UniProtKB', to_db="UniProtKB", ids=member_ids_group)
                    if check_id_mapping_results_ready(job_id):
                        link = get_id_mapping_results_link(job_id)
                        # get results
                        results = get_id_mapping_results_search(link)
                        # update the dictionary for each result
                        memberids2uniprotkb.update({x['from']: x['to']['primaryAccession'] for x in results['results']}) # appending
                # now we can convert IDs in each case and add them
                for uniref90id,member_ids in groups2members_dict.items():
                    # convert ids
                    member_ids_final = [memberids2uniprotkb.get(splitting_dots_and_commas(x), None) for x in member_ids if splitting_dots_and_commas(x) in memberids2uniprotkb.keys()]
                    # updating the results dict
                    res_dictionary['results'].update({uniref90id: member_ids_final})
        # returning dictionary of primary accessions
        return res_dictionary
    except Exception as e:
        print(e)
        return {'results': [{'primaryAccession': None}]}

def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
        
def get_uniprotkbid_length(ID_list):
    if isinstance(ID_list, list) and len(ID_list) > 1:
        id_request = ''
        for accession in ID_list:
            id_request = id_request + 'accession:' + accession + '%20OR%20'
        # getting rid of last separator
        id_request = id_request.rpartition('%20OR%20')[0]
        requestURL = "https://rest.uniprot.org/uniprotkb/search?query={0}&format=json&fields=accession,length&size=500".format(id_request)
    else:
        if isinstance(ID_list,list):
            ID_list = ID_list[0]
        requestURL = "https://rest.uniprot.org/uniprotkb/search?query=accession:{0}&format=json&fields=accession,length".format(ID_list)
    #print(requestURL)
    time.sleep(1)
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    return r.json()
        
def prune_group_by_length(uniprot_ids):
    # retrieve sequences lengths by batch of 400
    try:
        if len(uniprot_ids) == 1:
            return uniprot_ids
        elif len(uniprot_ids) >=  400:
            id_groups = split_list(input_list = uniprot_ids, chunk_size = 400)
        else:
            id_groups = [uniprot_ids]
        # go over the list of chunks
        # create a dictionary to allocate id: length
        dictionary_lengths_result = {}
        for id_group in id_groups:
            dictionary_lengths_result.update({x['primaryAccession']: x['sequence']['length'] for x 
                                              in get_uniprotkbid_length(id_group)['results']})
        # get mean and std deviation
        group_length_mean = np.mean([item for key,item in dictionary_lengths_result.items()])
        group_length_std = np.std([item for key,item in dictionary_lengths_result.items()])
        # now get rid of those sequences with lengths strongly deviating from the mean
        return [key for key,item in dictionary_lengths_result.items() if 
                    (item >= group_length_mean-group_length_std) and (item <= group_length_mean+group_length_std)]
    except Exception as e:
        print(e)

def param_to_num(param):
    if float(param).is_integer():
        return int(param)
    else:
        return float(param)

    
# get the API into work
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# performing pipeline
for foldseek_result in snakemake.input:
    # get PDB id
    pdb_id = foldseek_result.rpartition('/')[2].rpartition('.')[0].rpartition('foldseek_search_')[2]
    # load parameters
    prob_threshold = param_to_num(snakemake.params.prob_threshold)
    qcov_threshold = param_to_num(snakemake.params.qcov_threshold)
    scov_threshold = param_to_num(snakemake.params.scov_threshold)
    evalue_threshold = param_to_num(snakemake.params.evalue_threshold)
    # create target folder
    seed_folder = foldseek_result.rpartition('/')[0].rpartition('/')[0]
    target_dir = '{5}/{0}_prob_{1}_qcov_{2}_scov_{3}_eval_{4}_uniref90_homologs_structs/'.format(pdb_id,str(prob_threshold),str(qcov_threshold),str(scov_threshold),str(evalue_threshold), seed_folder)
    create_dir(target_dir)
    # load results table
    results_table = pd.read_csv(foldseek_result, sep = '\t')
    # filter the results based on the parameters specified
    results_table_filtered = (results_table.query("db == 'afdb50'").
                       assign(qcov = lambda df: abs(df.qEndPos-df.qStartPos)*100/df.qLen,
                              scov = lambda df: abs(df.dbEndPos-df.dbStartPos)*100/df.dbLen).
                       query("prob >= @prob_threshold and qcov >= @qcov_threshold and scov >= @scov_threshold and eval <= @evalue_threshold"))
    # get UniRef90 clusters
    foldseek_hit_ids = [x.split(' ')[0].split('-')[1] for x in results_table_filtered['target'].to_list()]
    # get sequences belonging to the clusters
    # chunk into lists of 400 members
    print('                                          '+' '*(len('[seed structure: {0}]'.format(pdb_id))+1))
    print('##########################################'+'#'*(len('[seed structure: {0}]'.format(pdb_id))+1))
    print('Retrieving target UniRef90 clusters IDs... [seed structure: {0}]'.format(pdb_id))
    print('##########################################'+'#'*(len('[seed structure: {0}]'.format(pdb_id))+1))
    print('                                          '+' '*(len('[seed structure: {0}]'.format(pdb_id))+1))
    foldseek_hit_ids_chunks = split_list(input_list= list(set(foldseek_hit_ids)), 
                                         chunk_size= 400)
    foldseek_clusterids_uni90_groups = [get_uniref90_cluster_id(hit_id_list) for hit_id_list in tqdm.tqdm(foldseek_hit_ids_chunks) if 
                                           not isinstance(get_uniref90_cluster_id(hit_id_list), dict)]
    # converting into a flat list and de-duplicate
    foldseek_clusterids_uni90_groups_full = [hit_id for uniref_list in foldseek_clusterids_uni90_groups for hit_id in uniref_list]
    foldseek_clusterids_uni90_groups_full = list(set(foldseek_clusterids_uni90_groups_full))
    # split the list into chunks of 400 proteins
    foldseek_clusterids_uni90_chunks = split_list(input_list= list(set(foldseek_clusterids_uni90_groups_full)), 
                                                  chunk_size= 400)
    # get the members
    print('                                              ')
    print('##############################################')
    print('Retrieving target UniRef90 clusters members...')
    print('##############################################')
    print('                                              ')
    foldseek_uniref90_members = [get_uniprotkbid_uniref90(x) for x in tqdm.tqdm(foldseek_clusterids_uni90_chunks)] # get_uniprotkbid_uniref50(chunk)
    
    # get length of sequences for each cluster
    selected_groups = {}
    for group,member_list in tqdm.tqdm(foldseek_uniref90_members[0]['results'].items()):
        # prune the groups
        if len(member_list) > 1:
            selected_seqs = prune_group_by_length(uniprot_ids = member_list)
            selected_groups.update({group: selected_seqs})
        else:
            # update the dictionary
            selected_groups.update({group: member_list})
    # filter sequences by length
    # get the whole list
    selected_pdb_ids = [pdb_id for key,item in selected_groups.items() for pdb_id in item]
    # save ID list as identifiers.txt file
    print('                                ')
    print('################################')
    print('Creating identifiers.txt file...')
    print('################################')
    print('                                ')
    with open(target_dir+'/identifiers.txt', 'w') as output_file:
        for ID in selected_pdb_ids:
            output_file.write(ID)
            output_file.write('\n')
    print('     ')
    print('#####')
    print('Done!')
    print('#####')
    print('     ')