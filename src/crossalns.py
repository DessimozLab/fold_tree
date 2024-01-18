
#copy the alignment logic from structal to sequence and vice versa
import Bio.SeqIO
import Bio.AlignIO

def crossaln(aln1f, aln2f, outaln):
    #copy aln1 logic to aln2
    #aln1 and al2 are both file names

    #read in the alignments
    aln1 = [ rec for rec in Bio.SeqIO.parse(aln1f, "fasta")]
    aln2 = [ rec for rec in  Bio.SeqIO.parse(aln2f, "fasta") ]

    #get the sequences by removing the gaps
    seqs1= { rec.id: str(rec.seq) for rec in Bio.SeqIO.parse(aln1f, "fasta") }
    seqs2= { rec.id: str(rec.seq).replace("-","") for rec in Bio.SeqIO.parse(aln2f, "fasta") }

    print('seqs1',seqs1.keys())
    print('seqs2',seqs2.keys())
    #copy the gap positions from aln1 to aln2
    for rec in aln2:
        seq = seqs2[rec.id]
        seq = list(seq)
        i = 0
        for c in seqs1[rec.id]:
            if c == "-":
                seq.insert(i, "-")
            else:
                i += 1
        seqs2[rec.id] = "".join(seq)
    #write out the new alignment
    with open(outaln, "w") as f:
        for rec in aln2:
            f.write(">%s\n%s\n" % (rec.id, seqs2[rec.id]))
    return outaln

if __name__ == "__main__":

    #snakemake input and output files
    aln1 = snakemake.input[0]
    aln2 = snakemake.input[1]
    outaln = snakemake.output[0]


    crossaln(aln1, aln2, outaln)

    aln1 = snakemake.input[1]
    aln2 = snakemake.input[0]
    outaln = snakemake.output[1]

    crossaln(aln1, aln2, outaln)
