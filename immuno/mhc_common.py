import re

def seq_to_str(obj):
    """
    Given a sequence convert it to a comma separated string. 
    If, however, the argument is a single object, return its string representation.
    """
    if isinstance(obj, (unicode, str)):
        return obj
    elif isinstance(obj, (list, tuple)):
        return  ",".join([str(x) for x in obj])
    else:
        return str(obj)

def convert_str(obj):
    """
    Given a string, convert it to an int or float if possible.
    """
    if obj is None: 
        return obj 
    try:
        try:
            return int(obj)
        except:
            return float(obj)
    except:
        return str(obj)

def normalize_hla_allele_name(hla):
    """
    HLA allele names can look like:
        - HLA-A*03:02
        - HLA-A02:03
        - HLA-A:02:03
        - HLA-A2
        - A2 
        - A*03:02
        - A02:02
        - A:02:03
    ...should all be normalized to:
        HLA-A*03:02:03
    """
    hla = hla.strip().upper()
    match = re.match('(HLA\-)?([A-Z])(\*|:)?([0-9][0-9]?):?([0-9][0-9]?)$', hla)
    assert match, "Malformed HLA type %s" % hla 
    (_, gene, _, family, protein) = match.groups()
    if len(family) == 1:
        family = "0" + family 
    if len(protein) == 1:
        protein = "0" + protein 
    return "HLA-%s*%s:%s" % (gene, family, protein )
