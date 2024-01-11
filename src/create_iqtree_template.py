

template_str = '''
#nexus


begin sets;
    charset part1 = aa.final.fasta: *;
    charset part2 = 3di.final.fasta: *;
    charpartition mine = LG+I+G:part1, ./3di_substmat.txt+I+F:part2;
end;
'''

template_str = template_str.replace('aa.final.fasta', snakemake.input[0])
template_str = template_str.replace('3di.final.fasta', snakemake.input[1])
template_str = template_str.replace('./3di_substmat.txt', snakemake.params.submat)
print(template_str)
with open(snakemake.output[0], 'w') as f:
    f.write(template_str)


template_str = '''
#nexus


begin sets;
    charset part1 = aa.final.fasta: *;
    charset part2 = 3di.final.fasta: *;
    charpartition mine = LG+I+G:part1, ./3di_substmat.txt+I+F:part2;
end;
'''

template_str = template_str.replace('aa.final.fasta', snakemake.input[2])
template_str = template_str.replace('3di.final.fasta', snakemake.input[3])
template_str = template_str.replace('./3di_substmat.txt', snakemake.params.submat)
print(template_str)
with open(snakemake.output[1], 'w') as f:
    f.write(template_str)


