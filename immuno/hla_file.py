import logging
from mhc_common import normalize_hla_allele_name

def read_hla_file(path, permissive_parsing=True):
    """
    Read in HLA alleles and normalize them, returning a list of HLA allele
    names.
    """
    assert path.endswith(".hla"), \
        "Expected HLA file %s to end with suffix .hla" % path

    logging.info("Reading HLA file %s", path)
    alleles = []
    with open(path, 'r') as f:
        contents = f.read()
        for line in contents.split("\n"):
            for raw_allele in line.split(","):
                if permissive_parsing:
                    # get rid of surrounding whitespace
                    raw_allele = raw_allele.strip()
                    # sometimes we get extra columns with scores,
                    # ignore those
                    raw_allele = raw_allele.split(" ")[0]
                    raw_allele = raw_allele.split("\t")[0]
                    raw_allele = raw_allele.split("'")[0]
                if len(raw_allele) > 0:
                    alleles.append(
                        normalize_hla_allele_name(
                            raw_allele))
    return alleles

