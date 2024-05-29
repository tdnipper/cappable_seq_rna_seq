import pandas as pd

def annotate_genes(file):
    f = pd.read_csv(file)
    basename = file.split('.')[0]

    # Parse GTF file and create a dictionary mapping gene IDs to gene names
    gene_dict = {}
    with open('genome/hybrid_annotated_cat.gtf', 'r') as gtf:
        for line in gtf:
            if "gene_name" in line:
                gene_id = line.split()[9].strip('"').strip('";')
                gene_name = line.split("gene_name")[1].split(';')[0].strip().strip('"')
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = gene_name
                else:
                    continue

    # Annotate genes in DataFrame
    try:
        f["gene_name"] = f.iloc[:,0].map(gene_dict)
    except: f["gene_name"] = 'unnamed'
    
    f.rename(columns={f.columns[0]: 'gene_id'}, inplace=True)
    return f.to_csv(f'{basename}_annotated.csv', index=False)

annotate_genes('res05_infection.csv')
annotate_genes('res05_condition.csv')