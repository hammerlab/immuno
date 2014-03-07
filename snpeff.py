class SnpEffFields():

  SNPEFF_INFO_FIELD_KEY = "EFF"

  EFFECT_KEY            = "SNPEFF_EFFECT"
  IMPACT_KEY            = "SNPEFF_IMPACT"
  FUNCTIONAL_CLASS_KEY  = "SNPEFF_FUNCTIONAL_CLASS"
  CODON_CHANGE_KEY      = "SNPEFF_CODON_CHANGE"
  AMINO_ACID_CHANGE_KEY = "SNPEFF_AMINO_ACID_CHANGE"
  GENE_NAME_KEY         = "SNPEFF_GENE_NAME"
  GENE_BIOTYPE_KEY      = "SNPEFF_GENE_BIOTYPE"
  TRANSCRIPT_ID_KEY     = "SNPEFF_TRANSCRIPT_ID"
  EXON_ID_KEY           = "SNPEFF_EXON_ID"

class SnpEffEffectType():
  SPLICE_SITE_ACCEPTOR = "SPLICE_SITE_ACCEPTOR"
  SPLICE_SITE_DONOR = "SPLICE_SITE_DONOR"
  START_LOST = "START_LOST"
  EXON_DELETED = "EXON_DELETED"
  FRAME_SHIFT = "FRAME_SHIFT"
  STOP_GAINED = "STOP_GAINED"
  STOP_LOST = "STOP_LOST"

  NON_SYNONYMOUS_CODING = "NON_SYNONYMOUS_CODING"
  CODON_CHANGE = "CODON_CHANGE"
  CODON_INSERTION = "CODON_INSERTION"
  CODON_CHANGE_PLUS_CODON_INSERTION = "CODON_CHANGE_PLUS_CODON_INSERTION"
  CODON_DELETION = "CODON_DELETION"
  CODON_CHANGE_PLUS_CODON_DELETION = "CODON_CHANGE_PLUS_CODON_DELETION"
  UTR_5_DELETED = "UTR_5_DELETED"
  UTR_3_DELETED = "UTR_3_DELETED"

  SYNONYMOUS_START = "SYNONYMOUS_START"
  NON_SYNONYMOUS_START = "NON_SYNONYMOUS_START"
  START_GAINED = "START_GAINED"
  SYNONYMOUS_CODING = "SYNONYMOUS_CODING"
  SYNONYMOUS_STOP = "SYNONYMOUS_STOP"
  NON_SYNONYMOUS_STOP = "NON_SYNONYMOUS_STOP"

  NONE = "NONE"
  CHROMOSOME = "CHROMOSOME"
  CUSTOM = "CUSTOM"
  CDS = "CDS"
  GENE = "GENE"
  TRANSCRIPT = "TRANSCRIPT"
  EXON = "EXON"
  INTRON_CONSERVED = "INTRON_CONSERVED"
  UTR_5_PRIME = "UTR_5_PRIME"
  UTR_3_PRIME = "UTR_3_PRIME"
  DOWNSTREAM = "DOWNSTREAM"
  INTRAGENIC = "INTRAGENIC"
  INTERGENIC = "INTERGENIC"
  INTERGENIC_CONSERVED = "INTERGENIC_CONSERVED"
  UPSTREAM = "UPSTREAM"
  REGULATION = "REGULATION"
  INTRON = "INTRON"



class SnpEffEffect():
  """

  EFF=EXON
  (MODIFIER|||||FTCD|processed_transcript|CODING|ENST00000498355|6|1),
  FRAME_SHIFT(HIGH||-|-221|495|FTCD|protein_coding|CODING|ENST00000355384|6|1)

  """
  EFFECT_PATTERN = re.compile(r'(\w+)\(.+\)')
  EFFECT_FIELDS = re.compile(r'\w+\((.+)\)')


  def __init__(self, line = None):
    # self.effect 0
    # self.impact 1
    # self.functionalClass
    # self.codonChange
    # self.aminoAcidChange
    # self.geneName
    # self.geneBiotype
    # self.coding
    # self.transcriptID
    # self.exonID

    if line is not None:
      self.parse(line)


  def parse( self, record ):
    self.effect = SnpEffEffect.EFFECT_PATTERN.match(record).group(1)
    fields = SnpEffEffect.EFFECT_FIELDS.match(record).group(1)
    self.full_line = fields
    fields = fields.split("|")

    self.impact = fields[0]
    self.functionalClass = fields[1]

    self.transcriptId = fields[8]

  def get_transcript(self):
    return Ensembl.get_transcripts(self.transcriptID)


if __name__ == '__main__':
  SnpEffEffect("EXON(MODIFIER|||||FTCD|processed_transcript|CODING|ENST00000498355|6|1)")

