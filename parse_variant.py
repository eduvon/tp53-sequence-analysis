def parse_vcf_record(line):
    fields = line.split("\t")

    chrom = fields[0]
    position = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    info = fields[7]

    data = info.split(";")
    data_dict = dict(item.split("=", 1) for item in data)

    consequences = []

    if "MC" in data_dict:
        for consequence in data_dict["MC"].split(","):
            description = consequence.split("|")[1]
            consequences.append(description)
    else:
        consequences = [data_dict["CLNVC"]]

    if "RS" in data_dict:
        rs = f"rs{data_dict['RS']}"
    else:
        rs = None

    return {
        "chrom": chrom,
        "position": position,
        "ref": ref,
        "alt": alt,
        "clinical_significance": data_dict.get("CLNSIG"),
        "consequences": consequences,
        "rs": rs,
        "genomic_hgvs": data_dict.get("CLNHGVS")
    }