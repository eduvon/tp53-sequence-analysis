info = "ALLELEID=347484;CLNDN=Li-Fraumeni_syndrome_1;CLNHGVS=NC_000017.11:g.7668539G>A;CLNSIG=Benign/Likely_benign;GENEINFO=TP53:7157;MC=SO:0001624|3_prime_UTR_variant;ORIGIN=1;RS=114831472"

data =info.split(";")

data_dict = dict(item.split("=", 1) for item in data)


# print(f"Clinical significance:", data_dict['CLNSIG'])
# print(f"Variant type:", data_dict['MC'].split('|')[1])
# print("dbSNP:", f"rs{data_dict['RS']}")

line = "17	7668434	35555	T	G	.	.	AF_TGP=0.00260;ALLELEID=44228;CLNDISDB=.|MONDO:MONDO:0015356,MeSH:D009386,MedGen:C0027672,Orphanet:140162|MedGen:C3661900|MONDO:MONDO:0013876,MedGen:C3553606,OMIM:614740|Gene:553989,MedGen:C1835398,OMIM:151623,Orphanet:524|MONDO:MONDO:0018875,MedGen:C0085390,OMIM:PS151623,Orphanet:524;CLNDN=TP53-related_disorder|Hereditary_cancer-predisposing_syndrome|not_provided|Basal_cell_carcinoma,_susceptibility_to,_7|Li-Fraumeni_syndrome_1|Li-Fraumeni_syndrome;CLNHGVS=NC_000017.11:g.7668434T>G;CLNREVSTAT=reviewed_by_expert_panel;CLNSIG=Benign;CLNVC=single_nucleotide_variant;GENEINFO=TP53:7157;MC=SO:0001624|3_prime_UTR_variant;ORIGIN=1;RS=78378222"

def parse_vcf_record(line):
    fields = line.split("\t")

    chrom = fields[0]
    position = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    info = fields[7]

    data =info.split(";")
    data_dict = dict(item.split("=", 1) for item in data)

    return {
        "chrom": chrom,
        "position": position,
        "ref": ref,
        "alt": alt,
        "clinical_significance": data_dict['CLNSIG'],
        "variant_type": data_dict['MC'].split('|')[1],
        "rs": f"rs{data_dict['RS']}"
    }

variant = parse_vcf_record(line)
print(variant)