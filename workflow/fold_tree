
#once snakemake is installed use the following command to test the struct tree
import snakemake.utils
snakemake.utils.min_version("7.8.0")
snake_dir = workflow.basedir

rootdir = ''.join([ sub + '/' for sub in snake_dir.split('/')[:-1] ] )
print(rootdir)
mattypes = ['foldtree', 'alntmscore', 'lddt']

configfile: rootdir+ "workflow/config/config_vars.yaml"
# remote homologues search parameters

foldseekpath = config["foldseek_path"]
if foldseekpath == 'provided':
	foldseekpath = rootdir + "foldseek/foldseek"

rule all:
	input:
	#get all treescore and rf distance files for all alntypes
		#expand( "{folder}/RFdistances_{exp}_.json" , folder = folders , exp = exp) ,
		expand( "{folder}/plddt.json" , folder = config["folder"] ) ,
		expand("{folder}/{mattype}_struct_tree.PP.nwk.rooted.final", folder = config["folder"], mattype = mattypes),

rule mad_root_post:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/{mattype}_struct_tree.PP.nwk"
	output:
		"{folder}/{mattype}_struct_tree.PP.nwk.rooted.final"
	log:
		"{folder}/logs/{mattype}_struct_madroot_post.log"
	script:
		'../src/process_madroot.py'
	
rule mad_root_struct:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/{mattype}_struct_tree.PP.nwk"
	output:
		"{folder}/{mattype}_struct_tree.PP.nwk.rooted"
	log:
		"{folder}/logs/{mattype}_struct_madroot.log"
	shell:
		rootdir+'madroot/mad {wildcards.folder}/{wildcards.mattype}_struct_tree.PP.nwk'

rule postprocess:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/{mattype}_struct_tree.nwk"
	output:
		"{folder}/{mattype}_struct_tree.PP.nwk"
	log:
		"{folder}/logs/{mattype}_struct_postprocess.log"
	script:
		'../src/postprocess.py'

rule quicktree:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/{mattype}_fastmemat.txt"
	output:
		"{folder}/{mattype}_struct_tree.nwk"
	log:
		"{folder}/logs/{mattype}_quicktree.log"
	shell:
		'quicktree -i m {wildcards.folder}/{wildcards.mattype}_fastmemat.txt > {wildcards.folder}/{wildcards.mattype}_struct_tree.nwk '

rule foldseek2distmat:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/allvall_1.csv"
	output:
		"{folder}/foldtree_fastmemat.txt",
		"{folder}/alntmscore_fastmemat.txt",
		"{folder}/lddt_fastmemat.txt",
	params:
		fmt = None,
	log:
		"{folder}/logs/foldseek2distmat.log"
	script:
		"../src/foldseekres2distmat_simple.py"

rule foldseek_allvall_1:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/finalset.csv"
	output:
		"{folder}/allvall_1.csv"
	log:
		"{folder}/logs/foldseekallvall.log"
	shell:
		foldseekpath + " easy-search {wildcards.folder}/structs/ {wildcards.folder}/structs/ {wildcards.folder}/allvall_1.csv {wildcards.folder}/tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore' --exhaustive-search --alignment-type 2 -e inf --threads " + str(config['foldseek_cores']) 

rule dl_ids_sequences:
	conda: 
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		ids="{folder}/identifiers.txt",
	output:
		"{folder}/sequence_dataset.csv",
	log:
		"{folder}/logs/dlsequences.log"
	params:
		cath=config["cath"],
		custom_structs=config["custom_structs"],
		clean_folder=config["clean_folder"],
	script:
		"../src/dl_sequences.py"

rule plddt:
	conda:
		#"config/fold_tree.yaml"
		"foldtree"
	input:
		"{folder}/finalset.csv",
	output:
		"{folder}/plddt.json",
	log:
		"{folder}/logs/plddt.log"
	script:
		'../src/grabplddt.py'

rule dl_ids_structs:
	input:
		"{folder}/sequence_dataset.csv",
	output:
		#"{folder}/sequences.fst",
		"{folder}/finalset.csv",
	conda: 
		#"config/fold_tree.yaml"
		"foldtree",
	log:
		"{folder}/logs/dlstructs.log",
	params:
		filtervar=config["filter"],
		cath=False,
		filtervar_min=config["filter_min"],
		filtervar_avg=config["filter_avg"],
		custom_structs=config["custom_structs"],
		clean_folder=config["clean_folder"],
	script:
		"../src/dl_structs.py"