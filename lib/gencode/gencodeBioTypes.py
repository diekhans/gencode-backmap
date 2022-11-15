"""
Declarations and operations on biotypes.  Used in ensuring that
all biotypes are can be mapped to a reduced set for coloring and
selection in the browser.
"""
from pycbio.sys.symEnum import SymEnum, SymEnumValue, auto
from pycbio.gencode import GencodeGenesException

class BioType(SymEnum):
    overlapping_ncRNA_3prime = SymEnumValue(auto(), "3prime_overlapping_ncRNA")
    antisense = auto()
    antisense_RNA = antisense
    bidirectional_promoter_lncRNA = auto()
    IG_C_gene = auto()
    IG_C_pseudogene = auto()
    IG_D_gene = auto()
    IG_D_pseudogene = auto()
    IG_J_gene = auto()
    IG_J_pseudogene = auto()
    IG_LV_gene = auto()
    IG_pseudogene = auto()
    IG_V_gene = auto()
    IG_V_pseudogene = auto()
    lincRNA = auto()
    lncRNA = auto()
    macro_lncRNA = auto()
    miRNA = auto()
    misc_RNA = auto()
    Mt_rRNA = auto()
    Mt_tRNA = auto()
    non_coding = auto()
    nonsense_mediated_decay = auto()  # transcripts only
    non_stop_decay = auto()  # transcripts only
    polymorphic_pseudogene = auto()
    processed_pseudogene = auto()
    processed_transcript = auto()
    protein_coding = auto()
    pseudogene = auto()
    retained_intron = auto()  # transcripts only
    ribozyme = auto()
    rRNA = auto()
    rRNA_pseudogene = auto()
    scaRNA = auto()
    scRNA = auto()
    sense_intronic = auto()
    sense_overlapping = auto()
    snoRNA = auto()
    snRNA = auto()
    sRNA = auto()
    TEC = auto()
    transcribed_processed_pseudogene = auto()
    transcribed_unitary_pseudogene = auto()
    transcribed_unprocessed_pseudogene = auto()
    translated_processed_pseudogene = auto()
    translated_unprocessed_pseudogene = auto()
    TR_C_gene = auto()
    TR_D_gene = auto()
    TR_J_gene = auto()
    TR_J_pseudogene = auto()
    TR_V_gene = auto()
    TR_V_pseudogene = auto()
    unitary_pseudogene = auto()
    unprocessed_pseudogene = auto()
    vaultRNA = auto()
    vault_RNA = auto()
    protein_coding_LoF = auto()
    artifact = auto()
    artifactual_duplication = auto()
    protein_coding_CDS_not_defined = auto()

class GencodeFunction(SymEnum):
    """category used to color UCSC browser tracks"""
    coding = auto()
    nonCoding = auto()
    pseudo = auto()
    other = auto()


class GnecodeGeneCategory(SymEnum):
    """gene categories uses in GENCODE stats reporting, from
     http://ftp.ebi.ac.uk/pub/databases/gencode/_README_stats.txt"""
    protein_coding = auto()
    lncRNA = auto()
    smallRNA = auto()
    pseudoGene = auto()
    immunoSegment = auto()


bioTypesTranscriptOnly = frozenset([BioType.nonsense_mediated_decay,
                                    BioType.non_stop_decay,
                                    BioType.retained_intron])

bioTypesCoding = frozenset([BioType.IG_C_gene,
                            BioType.IG_D_gene,
                            BioType.IG_J_gene,
                            BioType.IG_V_gene,
                            BioType.IG_LV_gene,
                            BioType.polymorphic_pseudogene,
                            BioType.protein_coding_LoF,
                            BioType.protein_coding_CDS_not_defined,
                            BioType.protein_coding,
                            BioType.nonsense_mediated_decay,
                            BioType.TR_C_gene,
                            BioType.TR_D_gene,
                            BioType.TR_J_gene,
                            BioType.TR_V_gene,
                            BioType.non_stop_decay])
bioTypesNonCoding = frozenset([BioType.antisense,
                               BioType.lincRNA,
                               BioType.lncRNA,
                               BioType.miRNA,
                               BioType.misc_RNA,
                               BioType.Mt_rRNA,
                               BioType.Mt_tRNA,
                               BioType.non_coding,
                               BioType.processed_transcript,
                               BioType.rRNA,
                               BioType.snoRNA,
                               BioType.scRNA,
                               BioType.snRNA,
                               BioType.overlapping_ncRNA_3prime,
                               BioType.sense_intronic,
                               BioType.sense_overlapping,
                               BioType.macro_lncRNA,
                               BioType.ribozyme,
                               BioType.scaRNA,
                               BioType.sRNA,
                               BioType.vaultRNA,
                               BioType.vault_RNA,
                               BioType.bidirectional_promoter_lncRNA])
bioTypesLncRna = frozenset([BioType.lincRNA,
                            BioType.lncRNA])
bioTypesOther = frozenset([BioType.retained_intron,
                           BioType.TEC,
                           BioType.artifact,
                           BioType.artifactual_duplication])
bioTypesPseudo = frozenset([BioType.IG_J_pseudogene,
                            BioType.IG_pseudogene,
                            BioType.IG_V_pseudogene,
                            BioType.processed_pseudogene,
                            BioType.pseudogene,
                            BioType.rRNA_pseudogene,
                            BioType.transcribed_processed_pseudogene,
                            BioType.transcribed_unprocessed_pseudogene,
                            BioType.unitary_pseudogene,
                            BioType.transcribed_unitary_pseudogene,
                            BioType.unprocessed_pseudogene,
                            BioType.IG_C_pseudogene,
                            BioType.IG_D_pseudogene,
                            BioType.TR_J_pseudogene,
                            BioType.TR_V_pseudogene,
                            BioType.translated_processed_pseudogene,
                            BioType.translated_unprocessed_pseudogene,
                            ])
bioTypesTranscribedPseudo = frozenset([BioType.transcribed_processed_pseudogene,
                                       BioType.transcribed_unprocessed_pseudogene,
                                       BioType.transcribed_unitary_pseudogene])

assert (len(bioTypesCoding) + len(bioTypesNonCoding) + len(bioTypesOther) + len(bioTypesPseudo)) == len(list(BioType))

bioTypesTR = frozenset((BioType.TR_C_gene,
                        BioType.TR_D_gene,
                        BioType.TR_J_gene,
                        BioType.TR_V_gene))
bioTypesIG = frozenset((BioType.IG_C_gene,
                        BioType.IG_D_gene,
                        BioType.IG_J_gene,
                        BioType.IG_V_gene))


# small, non-coding RNAs from
bioTypesSmallNonCoding = frozenset([BioType.miRNA,
                                    BioType.misc_RNA,
                                    BioType.Mt_rRNA,
                                    BioType.Mt_tRNA,
                                    BioType.ribozyme,
                                    BioType.rRNA,
                                    BioType.snoRNA,
                                    BioType.scRNA,
                                    BioType.snRNA,
                                    BioType.sRNA,
                                    BioType.scaRNA,
                                    BioType.vaultRNA,
                                    BioType.vault_RNA])

# imported from external databases
bioTypesNonCodingExternalDb = frozenset([BioType.miRNA,
                                         BioType.misc_RNA,
                                         BioType.Mt_rRNA,
                                         BioType.Mt_tRNA,
                                         BioType.rRNA,
                                         BioType.snoRNA,
                                         BioType.snRNA])

# GENE category mappings
gencodeCategoryProteinCoding = frozenset([BioType.protein_coding])
gencodeCategoryLncRNA = frozenset([BioType.processed_transcript,
                                   BioType.lincRNA,
                                   BioType.overlapping_ncRNA_3prime,
                                   BioType.antisense,
                                   BioType.non_coding,
                                   BioType.sense_intronic,
                                   BioType.sense_overlapping,
                                   BioType.TEC,
                                   BioType.macro_lncRNA,
                                   BioType.bidirectional_promoter_lncRNA,
                                   BioType.lncRNA])
gencodeCategorySmallRNA = frozenset([BioType.snRNA,
                                     BioType.snoRNA,
                                     BioType.rRNA,
                                     BioType.Mt_tRNA,
                                     BioType.Mt_rRNA,
                                     BioType.misc_RNA,
                                     BioType.miRNA,
                                     BioType.ribozyme,
                                     BioType.sRNA,
                                     BioType.scaRNA,
                                     BioType.vaultRNA])
gencodeCategoryPseudoGene = bioTypesPseudo
gencodeCategoryImmunoSegment = frozenset([BioType.IG_C_gene,
                                          BioType.IG_D_gene,
                                          BioType.IG_J_gene,
                                          BioType.IG_V_gene,
                                          BioType.IG_LV_gene,
                                          BioType.TR_C_gene,
                                          BioType.TR_D_gene,
                                          BioType.TR_J_gene,
                                          BioType.TR_V_gene])


def getFunctionForBioType(bt):
    """map a raw biotype to a function.  Note that transcript
    function isn't as simple as just translating this type,
    as gene biotype must be considered as well.
    """
    if bt in bioTypesCoding:
        return GencodeFunction.coding
    elif bt in bioTypesNonCoding:
        return GencodeFunction.nonCoding
    elif bt in bioTypesPseudo:
        return GencodeFunction.pseudo
    elif bt in bioTypesOther:
        return GencodeFunction.other
    else:
        raise GencodeGenesException("unknown biotype: " + str(bt))

def getTranscriptFunction(geneBioType, transcriptBioType):
    # all transcripts in pseudogenes are psuedogene transcripts
    if geneBioType in bioTypesPseudo:
        return GencodeFunction.pseudo
    else:
        return getFunctionForBioType(transcriptBioType)
