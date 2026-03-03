



def prepare_astral_input(uniprot_df , output):
    finalset = pd.read_csv(uniprot_df)
    finalset['species'] = finalset['Taxonomic lineage (Ids)'].map(lambda x: x.split(',')[-1].split('(')[0].strip())
    finalset['species'] = finalset['species'].map(lambda x: x.split('_')[0])
    mapper = dict(zip(finalset['Entry'], finalset['species']))
    with open(uniprot_df.replace('.csv','_speciesmap.txt' ) , 'w') as f:
        for k,v in mapper.items():
            f.write('{}\t{}\n'.format(k,int(v)))
    
    #get ncbi tree of species set
    species_set = list(finalset.species.unique())
    species_set = [int(s) for s in species_set]
    st = ncbi.get_topology(species_set, intermediate_nodes=False)
    st.name = 'root'
    #write with internal node names

    st.write(outfile=uniprot_df.replace('.csv','_ncbi_tree.nwk' ) , format= 1 )
    return mapper , uniprot_df.replace('.csv','_speciesmap.txt' )  , uniprot_df.replace('.csv','_ncbi_tree.nwk' ) 


