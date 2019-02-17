"""
Objects to store gencode gene annotations
"""
from pycbio.sys import PycbioException
from pycbio.sys.enumeration import Enumeration

BioType = Enumeration("BioType", ("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "IG_LV_gene",
                                  "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
                                  "TR_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
                                  "ambiguous_orf", "antisense", "lincRNA",
                                  "Mt_rRNA", "Mt_tRNA", "Mt_tRNA_pseudogene",
                                  "TEC",
                                  "TR_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
                                  "miRNA", "miRNA_pseudogene", "misc_RNA",
                                  "ncrna_host", "misc_RNA_pseudogene",
                                  "non_coding", "nonsense_mediated_decay",
                                  "polymorphic_pseudogene",
                                  "processed_pseudogene",
                                  "processed_transcript", "protein_coding",
                                  "pseudogene", "rRNA", "rRNA_pseudogene",
                                  "retained_intron", "retrotransposed",
                                  "scRNA_pseudogene", "snRNA",
                                  "snRNA_pseudogene", "snoRNA",
                                  "snoRNA_pseudogene", "tRNA_pseudogene",
                                  "transcribed_processed_pseudogene",
                                  "transcribed_unprocessed_pseudogene",
                                  "unitary_pseudogene",
                                  "transcribed_unitary_pseudogene",
                                  "unprocessed_pseudogene", 
                                  ("overlapping_ncrna_3prime", "3prime_overlapping_ncrna"),
                                  "disrupted_domain",
                                  "sense_intronic", "sense_overlapping",
                                  "non_stop_decay", "translated_processed_pseudogene",
                                  "translated_unprocessed_pseudogene",
                                  "known_ncrna",
                                  "macro_lncRNA",
                                  "ribozyme",
                                  "scaRNA",
                                  "sRNA",
                                  "vaultRNA"))
BioStatus = Enumeration("BioStatus", ("KNOWN", "PUTATIVE", "NOVEL", "UNKNOWN", "KNOWN_BY_PROJECTION"))
BioTag = Enumeration("BioTag", ("2way_pseudo_cons", "CCDS", "PAR", "exp_conf",
                                "cds_end_NF", "cds_start_NF", "mRNA_end_NF",
                                "mRNA_start_NF",
                                ("not_organism_supported", "not_organism_supported", ("not_organism-supported",)), # alias is for typo in V7
                                "pseudo_consens", "seleno", "non_ATG_start",
                                "alternative_5_UTR", "alternative_3_UTR",
                                "non_canonical_other", "non_canonical_U12",
                                "non_canonical_TEC",
                                "non_canonical_conserved",
                                "non_canonical_polymorphism",
                                "non_canonical_genome_sequence_error",
                                ("non_submitted_evidence", "non_submitted_evidence", ("non-submitted_evidence",)),  # alias is old value
                                "NAGNAG_splice_site", "upstream_ATG",
                                "downstream_ATG", "NMD_exception",
                                "upstream_uORF", "overlapping_uORF",
                                ("not_best_in_genome_evidence", "not_best_in_genome_evidence", ("not_best-in-genome_evidence",)),  # alias is old value
                                "basic",
                                "readthrough_transcript",
                                "dotter_confirmed",
                                "RNA_Seq_supported_only",
                                "NMD_likely_if_extended",
                                "retained_intron_final",
                                "CAGE_supported_TSS",
                                "retained_intron_first",
                                "sequence_error",
                                "low_sequence_quality",
                                "bicistronic",
                                "retained_intron_CDS",
                                "RP_supported_TIS",
                                "inferred_exon_combination",
                                "appris_principal",
                                "appris_principal_1",
                                "appris_principal_2",
                                "appris_principal_3",
                                "appris_principal_4",
                                "appris_principal_5",
                                "appris_alternative_1",
                                "appris_alternative_2",
                                "appris_candidate_longest",
                                "appris_candidate",
                                "appris_candidate_ccds",
                                "appris_candidate_longest_ccds",
                                "appris_candidate_longest_seq",
                                "appris_candidate_highest_score",))
GencodeMethod = Enumeration("GencodeMethod", ("manual", "automatic"))
GencodeExtendedMethod = Enumeration("GencodeExtendedMethod", ("manualOnly", "automaticOnly", "merged"))

GencodeFunction = Enumeration("GencodeFunction", ("pseudo", "coding", "nonCoding", "problem"))
GencodeCompleteness = Enumeration("GencodeCompleteness", ("partial", "full"))

bioTypesCoding = frozenset([BioType.IG_C_gene,
                            BioType.IG_D_gene,
                            BioType.IG_J_gene,
                            BioType.IG_V_gene,
                            BioType.IG_LV_gene,
                            BioType.polymorphic_pseudogene,
                            BioType.protein_coding,
                            BioType.nonsense_mediated_decay,
                            BioType.TR_gene,
                            BioType.TR_C_gene,
                            BioType.TR_D_gene,
                            BioType.TR_J_gene,
                            BioType.TR_V_gene,
                            BioType.non_stop_decay])
bioTypesNonCoding = frozenset([BioType.antisense,
                               BioType.lincRNA,
                               BioType.miRNA,
                               BioType.misc_RNA,
                               BioType.ncrna_host,
                               BioType.Mt_rRNA,
                               BioType.Mt_tRNA,
                               BioType.non_coding,
                               BioType.processed_transcript,
                               BioType.rRNA,
                               BioType.snoRNA,
                               BioType.snRNA,
                               BioType.overlapping_ncrna_3prime,
                               BioType.sense_intronic,
                               BioType.sense_overlapping,
                               BioType.known_ncrna,
                               BioType.macro_lncRNA,
                               BioType.ribozyme,
                               BioType.scaRNA,
                               BioType.sRNA,
                               BioType.vaultRNA])
bioTypesProblem = frozenset([BioType.retained_intron,
                             BioType.TEC,
                             BioType.disrupted_domain,
                             BioType.ambiguous_orf])
bioTypesPseudo = frozenset([BioType.IG_J_pseudogene,
                            BioType.IG_pseudogene,
                            BioType.IG_V_pseudogene,
                            BioType.miRNA_pseudogene,
                            BioType.misc_RNA_pseudogene,
                            BioType.Mt_tRNA_pseudogene,
                            BioType.processed_pseudogene,
                            BioType.pseudogene,
                            BioType.rRNA_pseudogene,
                            BioType.scRNA_pseudogene,
                            BioType.snoRNA_pseudogene,
                            BioType.snRNA_pseudogene,
                            BioType.transcribed_processed_pseudogene,
                            BioType.transcribed_unprocessed_pseudogene,
                            BioType.tRNA_pseudogene,
                            BioType.TR_pseudogene,
                            BioType.unitary_pseudogene,
                            BioType.transcribed_unitary_pseudogene,
                            BioType.unprocessed_pseudogene,
                            BioType.IG_C_pseudogene,
                            BioType.IG_D_pseudogene,
                            BioType.TR_J_pseudogene,
                            BioType.TR_V_pseudogene,
                            BioType.retrotransposed,
                            BioType.translated_processed_pseudogene,
                            BioType.translated_unprocessed_pseudogene,
                            ])

assert((len(bioTypesCoding)+len(bioTypesNonCoding)+len(bioTypesProblem)+len(bioTypesPseudo)) == len(BioType.values))

bioTypesTR = frozenset((BioType.TR_gene,
                        BioType.TR_C_gene,
                        BioType.TR_D_gene,
                        BioType.TR_J_gene,
                        BioType.TR_V_gene))
bioTypesIG = frozenset((BioType.IG_C_gene,
                        BioType.IG_D_gene,
                        BioType.IG_J_gene,
                        BioType.IG_V_gene))

# characterized/uncharacterized ncRNAs for basic gencode set definition
bioTypesNonCodingCharacterized = frozenset([BioType.antisense,
                                            BioType.miRNA,
                                            BioType.Mt_rRNA,
                                            BioType.Mt_tRNA,
                                            BioType.rRNA,
                                            BioType.snoRNA,
                                            BioType.snRNA,
                                            BioType.ncrna_host,
                                            BioType.ribozyme,
                                            BioType.scaRNA,
                                            BioType.sRNA,
                                            BioType.vaultRNA,])
bioTypesNonCodingUncharacterized = frozenset([BioType.non_coding,
                                              BioType.lincRNA,
                                              BioType.processed_transcript,
                                              BioType.misc_RNA,
                                              BioType.overlapping_ncrna_3prime,
                                              BioType.sense_intronic,
                                              BioType.sense_overlapping,
                                              BioType.known_ncrna,
                                              BioType.macro_lncRNA,])
assert(len(bioTypesNonCodingCharacterized)+len(bioTypesNonCodingUncharacterized) == len(bioTypesNonCoding))

# imported from external databases
bioTypesNonCodingExternalDb = frozenset([BioType.miRNA,
                                         BioType.Mt_rRNA,
                                         BioType.Mt_tRNA,
                                         BioType.rRNA,
                                         BioType.snoRNA,
                                         BioType.snRNA,
                                         BioType.misc_RNA])
def getFunctionForBioType(bt):
    """map a raw biotype to a function.  Note that transcript
    function isn't as simple as just translating this type,
    as gene biotype must be considered as well.
    """
    if bt in bioTypesCoding:
        return GencodeFunction.coding
    if bt in bioTypesNonCoding:
        return GencodeFunction.nonCoding
    if bt in bioTypesPseudo:
        return GencodeFunction.pseudo
    if bt in bioTypesProblem:
        return GencodeFunction.problem
    raise GencodeGenesException("unknown biotype: " + str(bt))

class GencodeGenesException(PycbioException):
    "exception associated with Gencode Genes objects"
    def __init__(self, msg, cause=None):
        PycbioException.__init__(self, msg, cause)

