"""
Objects to store gencode gene annotations
"""
import sys, re
from pycbio.hgdata.genePred import GenePredReader, GenePredDbReader
from pycbio.tsv import TsvReader, TabFileReader
from pycbio.sys.enumeration import Enumeration
from pycbio.sys import PycbioException
from pycbio.hgdata.rangeFinder import RangeFinder

# FIXME: Currently tested by programs in:
#   modules/gencode/src/progs/gencodeMakeTracks/
# having separate unit tests would be helpful
# FIXME: no test for database access

# FIXME: would be better to have constructor table both files
# rather than load and finish methods onlyExisting=True is a hack,
# but need due to chrom filter.

# FIXME  FIXME FIXME: should figure out how to make this table with all the
# definitions, then create the enumerations from that.

# FIXME: sql select stuff is a bit hacky

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
                                  "scRNA", "scRNA_pseudogene",
                                  "snRNA", "snRNA_pseudogene",
                                  "snoRNA", "snoRNA_pseudogene",
                                  "tRNA_pseudogene",
                                  "transcribed_processed_pseudogene",
                                  "transcribed_unprocessed_pseudogene",
                                  "unitary_pseudogene",
                                  "transcribed_unitary_pseudogene",
                                  "unprocessed_pseudogene", 
                                  ("overlapping_ncrna_3prime", "3prime_overlapping_ncrna"),
                                  ("overlapping_ncRNA_3prime", "3prime_overlapping_ncRNA"),
                                  "disrupted_domain",
                                  "sense_intronic", "sense_overlapping",
                                  "non_stop_decay", "translated_processed_pseudogene",
                                  "translated_unprocessed_pseudogene",
                                  "known_ncrna",
                                  "macro_lncRNA",
                                  "ribozyme",
                                  "scaRNA",
                                  "sRNA",
                                  "vaultRNA",
                                  "bidirectional_promoter_lncrna",
                                  "bidirectional_promoter_lncRNA"))
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
                                "appris_candidate_highest_score",
                                "inferred_transcript_model",
                                "RNA_Seq_supported_partial"))
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
                               BioType.scRNA,
                               BioType.snRNA,
                               BioType.overlapping_ncrna_3prime,
                               BioType.overlapping_ncRNA_3prime,
                               BioType.sense_intronic,
                               BioType.sense_overlapping,
                               BioType.known_ncrna,
                               BioType.macro_lncRNA,
                               BioType.ribozyme,
                               BioType.scaRNA,
                               BioType.sRNA,
                               BioType.vaultRNA,
                               BioType.bidirectional_promoter_lncrna,
                               BioType.bidirectional_promoter_lncRNA])
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

bioTagNotFullCds = frozenset([BioTag.cds_start_NF,
                              BioTag.cds_end_NF])
bioTagNotFullMRna = frozenset([BioTag.mRNA_start_NF,
                              BioTag.mRNA_end_NF])
bioTagNotFull = frozenset(bioTagNotFullCds | bioTagNotFullMRna)

# characterized/uncharacterized ncRNAs for basic gencode set definition
bioTypesNonCodingCharacterized = frozenset([BioType.antisense,
                                            BioType.miRNA,
                                            BioType.Mt_rRNA,
                                            BioType.Mt_tRNA,
                                            BioType.rRNA,
                                            BioType.snoRNA,
                                            BioType.scRNA,
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
                                              BioType.overlapping_ncRNA_3prime,
                                              BioType.sense_intronic,
                                              BioType.sense_overlapping,
                                              BioType.known_ncrna,
                                              BioType.macro_lncRNA,
                                              BioType.bidirectional_promoter_lncrna,
                                              BioType.bidirectional_promoter_lncRNA])
assert(len(bioTypesNonCodingCharacterized)+len(bioTypesNonCodingUncharacterized) == len(bioTypesNonCoding))

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
                                    BioType.sRNA])

# imported from external databases
bioTypesNonCodingExternalDb = frozenset([BioType.miRNA,
                                         BioType.misc_RNA,
                                         BioType.Mt_rRNA,
                                         BioType.Mt_tRNA,
                                         BioType.rRNA,
                                         BioType.snoRNA,
                                         BioType.snRNA])
def _getFunctionForBioType(bt):
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

def _noneIfEmpty(s):
    return s if len(s) > 0 else None

def _emptyIfNone(s):
    return s if s is not None else ""

def _sourceToExtendedMethod(src):
    hasHav = src.find("havana") >= 0
    hasEns = src.find("ensembl") >= 0
    if hasHav and hasEns:
        return GencodeExtendedMethod.merged
    elif hasHav:
        return GencodeExtendedMethod.manualOnly
    else:
        return GencodeExtendedMethod.automaticOnly

class GencodeGenesException(PycbioException):
    "exception associated with Gencode Genes objects"
    def __init__(self, msg, cause=None):
        PycbioException.__init__(self, msg, cause)

class BioTags(set):
    "set of tags of type BioTag"
    def __init__(self, tagsInit=None):
        "tagsInit can be a string that is parsed or a list or tuple to initialize the set"
        if tagsInit is not None:
            if isinstance(tagsInit, str):
                self.__initFromStr(tagsInit)
            else:
                set.__init__(self, tagsInit)

    def __initFromStr(self, tagsStr):
        tagsStr = tagsStr.strip()
        if len(tagsStr) > 0:
            for t in tagsStr.split(","):
                t = t.strip()
                if t.startswith("seleno_"):  # deal with seleno_62
                    t = "seleno"
                self.add(BioTag(t))

    def __str__(self):
        t = [str(t) for t in self]
        t.sort()
        return ",".join(t)

    def getSorted(self):
        l = list(self)
        l.sort()
        return l

    def isFullLength(self):
        for tag in self:
            if tag in bioTagNotFull:
                return False
        return True

class GencodeTranscriptLocus(object):
    """object for a transcript at a given locus.  Need due to PAR"""
    __slots__ = ("transcript", "gp", "geneLocus")

    def __init__(self, transcript, gp):
        self.transcript = transcript
        self.geneLocus = None
        self.gp = gp

    def isSingleExon(self):
        return (len(self.gp.exons) == 1)

    def __str__(self):
        return self.gp.name + "\t" + self.gp.chrom + ":" + str(self.gp.txStart) + "-" + str(self.gp.txEnd)

    @property
    def id(self):
        return self.transcript.id
    @property
    def chrom(self):
        return self.gp.chrom
    @property
    def chromStart(self):
        return self.gp.txStart
    @property
    def chromEnd(self):
        return self.gp.txEnd
    @property
    def strand(self):
        return self.gp.strand

    def hasMultiLoci(self):
        "is the transcript associated with multiple loci?"
        return len(self.transcript.transcriptLoci) > 1
    
    def hasCds(self):
        return (self.gp.cdsStart < self.gp.cdsEnd)

    def __cmp__(self, other):
        "sort by genomic location"
        diff = cmp(self.chrom, other.chrom)
        if diff != 0:
            return diff
        diff = cmp(self.chromStart, other.chromStart)
        if diff != 0:
            return diff
        diff = cmp(self.chromEnd, other.chromEnd)
        if diff != 0:
            return diff
        return 0

class GencodeTranscript(object):
    "object for a single transcript"
    __slots__ = ("id", "name", "bioType", "havanaId", "ccdsId", "tags", "gene", "level", "extendedMethod", "transcriptSupportLevel", "transcriptLoci")

    def __init__(self, transcriptId):
        self.id = transcriptId
        self.name = None
        self.bioType = None
        self.havanaId = None
        self.ccdsId = None
        self.tags = BioTags()
        self.level = None
        self.gene = None
        self.extendedMethod = None
        self.transcriptSupportLevel = None   # only set from GTF/GFF3
        self.transcriptLoci = []  # GencodeTranscriptLocus objects

    def finish(self):
        "finish construction, sorting into predictable order"
        self.transcriptLoci.sort(key=lambda l: (l.gp.chrom, l.gp.txStart))

    def hasCds(self):
        return self.transcriptLoci[0].hasCds()

    def isCoding(self):
        return self.bioType in bioTypesCoding

    def isNonCoding(self):
        return self.bioType in bioTypesNonCoding

    def isPseudo(self):
        return self.bioType in bioTypesPseudo

    def isProblem(self):
        return self.bioType in bioTypesProblem

    def isHLA(self):
        "is this an HLA transcript?"
        return self.gene.isHLA()

    def isOlfactoryReceptor(self):
        "is this an olfactory receptor transcript?"
        return self.gene.isOlfactoryReceptor()

    def isImmunoglobin(self):
        return self.bioType in bioTypesIG

    def isTCellReceptor(self):
        return self.bioType in bioTypesTR

    def isFullLength(self):
        return self.tags.isFullLength()

    def isSingleExon(self):
        return self.transcriptLoci[0].isSingleExon()

    def getFunction(self):
        # all transcripts in pseudogenes are psuedogene transcripts
        if self.gene.getFunction() == GencodeFunction.pseudo:
            return GencodeFunction.pseudo
        else:
            return _getFunctionForBioType(self.bioType)

    def getCompleteness(self):
        return GencodeCompleteness.full if self.isFullLength() else GencodeCompleteness.partial

    def getMethod(self):
        return GencodeMethod.manual if self.havanaId is not None else GencodeMethod.automatic

    def getExtendedMethod(self):
        """get extended method that indicates merged transcripts.  Must have 
        loaded the Transcript_source metadata"""
        if self.extendedMethod is None:
            raise GencodeGenesException("Transcript_source metadata must be loaded to get transcript extended method: " + self.id)
        return self.extendedMethod

    __infoTsvBasicHeader = ("geneId", "geneName", "geneType", "geneStatus", "transcriptId", "transcriptName", "transcriptType", "transcriptStatus", "havanaGene", "havanaTranscript", "ccdsId", "level")
    __infoTsvWithTagsHeader = __infoTsvBasicHeader + ("tags",)

    @staticmethod
    def getInfoTsvHeader(inclTags=True):
        "get header for use toInfoRow"
        if inclTags:
            return GencodeTranscript.__infoTsvWithTagsHeader
        else:
            return GencodeTranscript.__infoTsvBasicHeader
    
    def toInfoRow(self, inclTags=True):
        # biostatus no longer supported
        row = [self.gene.id, self.gene.name, self.gene.bioType, "", self.id, self.name, self.bioType, "", _emptyIfNone(self.gene.havanaId),  _emptyIfNone(self.havanaId),  _emptyIfNone(self.ccdsId), self.level]
        if inclTags:
            row.append(str(self.tags))
        return row

class GencodeGeneLocus(object):
    """object to group all of the GencodeTranscriptLocus objects for a given
    GencodeGene.  Doesn't consider strand, since anti-sense transcripts are part of
    the same gene."""
    __slots__ = ("gene", "chrom", "chromStart", "chromEnd", "strand", "transcriptLoci")

    def __init__(self, gene, chrom):
        self.gene = gene
        self.chrom = chrom
        self.chromStart = None
        self.chromEnd = None
        self.strand = None
        self.transcriptLoci = []

    def finish(self):
        "finish construction, sorting into predictable order"
        self.transcriptLoci.sort(key=lambda l: (l.gp.chrom, l.gp.txStart))

    @property
    def id(self):
        return self.gene.id

    def isSingleExonGene(self):
        "is this a single-exon gene (all transcripts are single exon)"
        for tl in self.transcriptLoci:
            if not tl.isSingleExon():
                return False
        return True

    def add(self, transcriptLocus):
        assert(transcriptLocus.gp.chrom == self.chrom)
        self.transcriptLoci.append(transcriptLocus)
        if self.chromStart is None:
            self.chromStart = transcriptLocus.gp.txStart
            self.chromEnd = transcriptLocus.gp.txEnd
            self.strand = transcriptLocus.gp.strand
        else:
            self.chromStart = min(self.chromStart, transcriptLocus.gp.txStart)
            self.chromEnd = max(self.chromEnd, transcriptLocus.gp.txEnd)
            if self.strand != transcriptLocus.gp.strand:
                raise GencodeGenesException("gene has transcripts with different strands: " + self.gene.id)

class GencodeGene(object):
    "object for a gene"
    __slots__ = ("id", "name", "bioType", "havanaId", "transcripts", "extendedMethod", "geneLoci")

    def __init__(self, geneId):
        self.id = geneId
        self.name = None
        self.bioType = None
        self.havanaId = None
        self.extendedMethod = None
        self.transcripts = [] # GencodeTranscript objects
        self.geneLoci = []    # GencodeGeneLocus objects

    def finish(self):
        "finish construction, sorting into predictable order"
        self.transcripts.sort(key=lambda t: (t.id,))
        self.geneLoci.sort(key=lambda l: (l.chrom, l.chromStart))

    def hasCds(self):
        for trans in self.transcripts:
            if trans.hasCds():
                return True
        return False

    def isCoding(self):
        return self.bioType in bioTypesCoding

    def isNonCoding(self):
        return self.bioType in bioTypesNonCoding

    def isPseudo(self):
        return self.bioType in bioTypesPseudo

    def isProblem(self):
        return self.bioType in bioTypesProblem

    def isHLA(self):
        "is this an HLA gene?"
        return (self.name.find("HLA-") == 0)

    __orGeneSymRe = re.compile("^OR[0-9].+")
    def isOlfactoryReceptor(self):
        "is this an olfactory receptor gene?"
        return self.__orGeneSymRe.match(self.name) is not None

    def getFunction(self):
        return _getFunctionForBioType(self.bioType)

    def getMethod(self):
        return GencodeMethod.manual if self.havanaId is not None else GencodeMethod.automatic

    def getExtendedMethod(self):
        """get extended method that indicates merged genes.  Must have 
        loaded the Gene_source metadata"""
        if self.extendedMethod is None:
            raise GencodeGenesException("Gene_source metadata must be loaded to get gene extended method: " + self.id)
        return self.extendedMethod

    def obtainLocus(self, chrom):
        """obtain the GencodeGeneLocus object for the specified chromosome, or
        create one"""
        for locus in self.geneLoci:
            if locus.chrom == chrom:
                return locus
        locus = GencodeGeneLocus(self, chrom)
        self.geneLoci.append(locus)
        return locus

class GencodeGenes(object):
    "object to information and optional genePreds of all gencode genes"

    def __init__(self):
        self.transcriptsById = dict()
        self.genesById = dict()
        # build in a lazy manner
        self.transcriptLociByRange = None
        self.geneLociByRange = None

    def getTranscript(self, transcriptId):
        "get a transcript object, or none if it doesn't exist"
        return self.transcriptsById.get(transcriptId)

    def getTranscriptRequired(self, transcriptId):
        "get a transcript object, or exception if it doesn't exist"
        trans = self.transcriptsById.get(transcriptId)
        if trans is None:
            raise GencodeGenesException("transcript not found in attributes: " + transcriptId)
        return trans

    def getTranscripts(self):
        "get an iterator over all transcript objects"
        return self.transcriptsById.itervalues()

    def obtainTranscript(self, transcriptId):
        "get or create a transcript object"
        trans = self.getTranscript(transcriptId)
        if trans is None:
            trans = self.transcriptsById[transcriptId] = GencodeTranscript(transcriptId)
        return trans

    def getGene(self, geneId):
        "get a gene object, or none if it doesn't exist"
        return self.genesById.get(geneId)

    def getGeneRequired(self, geneId):
        "get a gene object, or none if it doesn't exist"
        gene = self.genesById.get(geneId)
        if gene is None:
            raise GencodeGenesException("gene not found in attributes: " + geneId)
        return gene

    def obtainGene(self, geneId):
        "get or create a gene object"
        gene = self.getGene(geneId)
        if gene is None:
            gene = self.genesById[geneId] = GencodeGene(geneId)
        return gene

    def getGeneLoci(self):
        "generator over all GencodeGeneLocus objects"
        for gene in self.genesById.itervalues():
            for geneLocus in gene.geneLoci:
                yield geneLocus
        
    def getOverlappingGeneLoci(self, chrom, chromStart, chromEnd, strand=None):
        """Get GeneLoci overlapping range.  This covers entry range of gene,
        including introns"""
        if self.geneLociByRange is None:
            self.__buildGeneLociByRange()
        return self.geneLociByRange.overlapping(chrom, chromStart, chromEnd, strand)

    def getTranscriptLoci(self):
        "generator over all GencodeTranscriptLocus objects"
        for trans in self.transcriptsById.itervalues():
            for transLocus in trans.transcriptLoci:
                yield transLocus

    def getOverlappingTranscriptLoci(self, chrom, chromStart, chromEnd, strand=None):
        """"get TranscriptLoci overlapping range. This only overlap with
        exons"""
        if self.transcriptLociByRange is None:
            self.__buildTranscriptLociByRange()
        return self.transcriptLociByRange.overlapping(chrom, chromStart, chromEnd, strand)

    def getTranscriptsSortByLocus(self):
        """get a list of transcript objects sorted to optimize access by
        location and test consistency."""
        transes = list(self.transcriptsById.itervalues())
        transes.sort(key=lambda t: (t.transcriptLoci[0].gp.chrom, t.transcriptLoci[0].gp.txStart, t.transcriptLoci[0].gp.txEnd, t.transcriptLoci[0].gp.name))
        return transes

    def getTranscriptsSortById(self):
        """get a list of transcript objects sorted by id."""
        transes = list(self.transcriptsById.itervalues())
        transes.sort(key=lambda t: (t.id,))
        return transes

    def __buildTranscriptLociByRange(self):
        "construct the transcriptLociByRange index when needed"
        self.transcriptLociByRange = RangeFinder()
        for transLocus in self.getTranscriptLoci():
            gp = transLocus.gp
            for exon in gp.exons:
                self.transcriptLociByRange.add(gp.chrom, exon.start, exon.end, transLocus, gp.strand)

    def __buildGeneLociByRange(self):
        "construct the geneLociByRange index when needed"
        self.geneLociByRange = RangeFinder()
        for geneLocus in self.getGeneLoci():
            self.geneLociByRange.add(geneLocus.chrom, geneLocus.chromStart, geneLocus.chromEnd, geneLocus, strand=geneLocus.strand)

    def __setInfoRowGene(self, info, gene):
        gene.name = info.geneName
        gene.bioType = info.geneType
        gene.havanaId = info.havanaGene
        
    def __loadInfoRowGene(self, info):
        gene = self.obtainGene(info.geneId)
        if gene.name is None:
            self.__setInfoRowGene(info, gene)
        elif not ((gene.name == info.geneName) and (gene.bioType == info.geneType)):
            raise GencodeGenesException("inconsistent gene meta-data from multiple info rows: " + info.geneId)
        return gene

    def __setInfoRowTranscript(self, info, gene, trans):
        trans.name = info.transcriptName
        trans.bioType = info.transcriptType
        trans.havanaId = info.havanaTranscript
        trans.ccdsId = info.ccdsId
        trans.level = info.level

    def __linkTranscriptToGene(self, gene, trans):
        trans.gene = gene
        gene.transcripts.append(trans)
        for transLocus in trans.transcriptLoci:
            geneLocus = gene.obtainLocus(transLocus.gp.chrom)
            geneLocus.add(transLocus)
            transLocus.geneLocus = geneLocus

    def __loadInfoRowTranscript(self, info, gene):
        # allows for rows being accidental duplicated
        trans = self.obtainTranscript(info.transcriptId)
        if trans.name is None:
            self.__setInfoRowTranscript(info, gene, trans)
        elif not ((trans.name == info.transcriptName) and (trans.bioType == info.transcriptType)):
            raise GencodeGenesException("inconsistent transcript meta-data from multiple info rows: " + info.transcriptId)
        # accumulate tags (PAR tag not on all copies).  May not have it tag if loading from database dump
        tags = getattr(info, "tags", None)
        if tags is not None:
            trans.tags |= tags
        if trans.gene is None:
            self.__linkTranscriptToGene(gene, trans)

    def __loadInfoRow(self, info, onlyExisting):
        if (not onlyExisting) or (info.transcriptId in self.transcriptsById):
            try:
                gene = self.__loadInfoRowGene(info)
                self.__loadInfoRowTranscript(info, gene)
            except Exception,ex:
                raise GencodeGenesException("error loading gencode info for %s" % (info.transcriptId,), ex), None, sys.exc_info()[2]

    def loadInfoFromFile(self, gencodeInfo, onlyExisting=True):
        """load file produced by the gencodeGtfToAttrs script.  If onlyExisting
        is true, only entries for transcripts already loaded from genePreds
        area included.
        """
        infoTypeMap = {"geneType": BioType, "transcriptType": BioType, "level": int, "tags": BioTags, "havanaGene": _noneIfEmpty, "havanaTranscript": _noneIfEmpty, "ccdsId": _noneIfEmpty}
        for info in TsvReader(gencodeInfo, typeMap=infoTypeMap):
            self.__loadInfoRow(info, onlyExisting)

    class InfoDbRow(object):
        "used to hold a row of attrs table from database query"
        def __init__(self, row):
            self.geneId = row[0]
            self.geneName = row[1]
            self.geneType = BioType(row[2])
            # row[3] - geneStatus no longer supported
            self.transcriptId = row[4]
            self.transcriptName = row[5]
            self.transcriptType = BioType(row[6])
            # row[7] transcriptStatus no longer supported
            self.havanaGene = _noneIfEmpty(row[8])
            self.havanaTranscript = _noneIfEmpty(row[9])
            self.ccdsId = _noneIfEmpty(row[10])
            self.level = int(row[11])
            self.transcriptClass = row[12]

    @staticmethod
    def __makeChromSubSelectClauses(chromSubSelects):
        clauses = ["(transcriptId in ({}))".format(chromSubSelect) for chromSubSelect in chromSubSelects]
        return "({})".format(" or ".join(clauses))
        
    @staticmethod
    def __makeTranscriptIdSelectWhere(chromSubSelects=None, transcriptIdSubSelect=None):
        """generate, possibly empty, where clause for tables keyed by tables,
        multiple chromSubSelects are or-ed due to pseudo genes being in
        another tables"""
        clauses = []
        if chromSubSelects is not None:
            clauses.append(GencodeGenes.__makeChromSubSelectClauses(chromSubSelects))
        if transcriptIdSubSelect is not None:
            clauses.append("(transcriptId in (" + transcriptIdSubSelect + "))")
        if len(clauses) > 0:
            return  " where " + " and ".join(clauses)
        else:
            return ""

    def loadInfoFromDb(self, conn, table, onlyExisting=True, chromSubSelects=None, transcriptIdSubSelect=None):
        """load meta information from attributes table.  If onlyExisting
        is true, only entries for transcripts already loaded from genePreds
        area included.  If transcriptIdSubSelect is specified, it should be
        a select to return just the transcriptIds of interest.  Useful
        for per-chrom filtering.
        """
        sql = "select geneId, geneName, geneType, geneStatus, transcriptId, transcriptName, transcriptType, transcriptStatus, havanaGeneId, havanaTranscriptId, ccdsId, level, transcriptClass from " + table \
            + self.__makeTranscriptIdSelectWhere(chromSubSelects, transcriptIdSubSelect)
        cursor = conn.cursor()
        try:
            cursor.execute(sql)
            while True:
                row = cursor.fetchone()
                if row is None:
                    break
                self.__loadInfoRow(self.InfoDbRow(row), onlyExisting)
        finally:
            cursor.close()

    def __checkLociOverlap(self, gp, transcript):
        for locus in transcript.transcriptLoci:
            if gp.chrom == locus.gp.chrom:
                raise GencodeGenesException("multiple annotations for \"%s\" on \"%s\"" % (gp.name, gp.chrom))

    def __loadGenePred(self, gp):
        "add a genePred"
        try:
            trans = self.obtainTranscript(gp.name)
            self.__checkLociOverlap(gp, trans)
            transLocus = GencodeTranscriptLocus(trans, gp)
            trans.transcriptLoci.append(transLocus)
        except Exception,ex:
            raise GencodeGenesException("error loading genePred for %s at %s:%d-%d" % (gp.name, gp.chrom, gp.txStart, gp.txEnd), ex), None, sys.exc_info()[2]

    def loadGenePredsFromFile(self, gpFile, chroms=None):
        "load genePreds from a file, optionally filtering for a list or set of chromosome"
        for gp in GenePredReader(gpFile):
            if (chroms is None) or (gp.chrom in chroms):
                self.__loadGenePred(gp)

    @staticmethod
    def __makeChromSqlWhere(chroms):
        return "(chrom in (" + ",".join(['"'+c+'"' for c in chroms]) + "))"

    @staticmethod
    def __makeGenePredSelectWhere(chroms=None, transcriptIdSubSelect=None):
        "generate, possibly empty, where clause for genePred tables"
        clauses = []
        if chroms is not None:
            clauses.append("(" + GencodeGenes.__makeChromSqlWhere(chroms) + ")")
        if transcriptIdSubSelect is not None:
            clauses.append("(name in ("+ transcriptIdSubSelect + "))")
        if len(clauses) > 0:
            return  " where " + " and ".join(clauses)
        else:
            return ""

    def loadGenePredsFromDb(self, conn, table, chroms=None, transcriptIdSubSelect=None):
        """load genePreds from a database, optionally filtering for list or set of chromosomes
        or names returned by a subselect"""
        sql = "select * from " + table + self.__makeGenePredSelectWhere(chroms, transcriptIdSubSelect)
        for gp in GenePredDbReader(conn, sql):
            self.__loadGenePred(gp)

    def __checkTransHasGene(self, transcript):
        "a transcript would not have a gene if it was in the genePred, but not info"
        if transcript.gene is None:
            raise GencodeGenesException("no gene id associated with transcript \"" + transcript.id + "\"")

    def loadTranscriptSourceFromFile(self, transcriptSource):
        "load transcript source metadata from a file"
        for row in TabFileReader(transcriptSource):
            trans = self.getTranscript(row[0])
            if trans is not None:
                trans.extendedMethod = _sourceToExtendedMethod(row[1])

    def loadGeneSourceFromFile(self, geneSource):
        "load gene source metadata from a file"
        for row in TabFileReader(geneSource):
            gene = self.getGene(row[0])
            if gene is not None:
                gene.extendedMethod = _sourceToExtendedMethod(row[1])

    def loadTagsFromFile(self, tags):
        "load tags file.  Used with table dumps when tags are not in info file. "
        for row in TabFileReader(tags):
            trans = self.getTranscript(row[0])
            if trans is not None:
                trans.tags.add(BioTag(row[1]))

    def loadTagsFromDb(self, conn, table, chromSubSelects=None, transcriptIdSubSelect=None):
        """load tags file.  Used with table dumps when tags are not in info
        file.  If transcriptIdSubSelect is specified, it should be a select to
        return just the transcriptIds of interest.  Useful for per-chrom
        filtering."""
        sql = "select transcriptId, tag from " + table \
            + self.__makeTranscriptIdSelectWhere(chromSubSelects, transcriptIdSubSelect)

        cursor = conn.cursor()
        try:
            cursor.execute(sql)
            while True:
                row = cursor.fetchone()
                if row is None:
                    break
                trans = self.getTranscript(row[0])
                if trans is not None:
                    trans.tags.add(BioTag(row[1]))
        finally:
            cursor.close()

    def check(self, errFh=None):
        """check the that all loci have associated transcripts/gene information.
        If errFh is supplied, then output all problems before raising an exception."""
        errCnt = 0
        for transcript in self.transcriptsById.itervalues():
            try:
                self.__checkTransHasGene(transcript)
            except GencodeGenesException, ex:
                errCnt += 1
                if errFh is not None:
                    errFh.write("Error: " + str(ex) + "\n")
                else:
                    raise ex
        if errCnt != 0: 
            raise GencodeGenesException("missing gene information for " + str(errCnt) + " loci, check info file")

    def finish(self, errFh=None):
        """run checking and sort internal lists so that results are deterministic"""
        for gene in  self.genesById.itervalues():
            gene.finish()
        self.check(errFh)

    @staticmethod
    def loadFromFiles(gpFile=None, gencodeInfo=None, geneSource=None, transcriptSource=None, tags=None, errFh=None, chroms=None, onlyExisting=True):
        "factory method to load from files (see load functions)"
        gencode = GencodeGenes()
        if gpFile is not None:
            gencode.loadGenePredsFromFile(gpFile, chroms=chroms)
        else:
            onlyExisting = False  # force loading only from info
        if gencodeInfo is not None:
            gencode.loadInfoFromFile(gencodeInfo, onlyExisting=onlyExisting)
        if geneSource is not None:
            gencode.loadGeneSourceFromFile(geneSource)
        if transcriptSource is not None:
            gencode.loadTranscriptSourceFromFile(transcriptSource)
        if tags is not None:
            gencode.loadTagsFromFile(tags)
        gencode.finish(errFh=errFh)
        return gencode

    @staticmethod
    def __makeChromSubSelects(gpTables, chroms=None):
        "generate subselect for chroms, or None"
        if chroms is not None:
            return ["select name from {}  where {}".format(gpTable, GencodeGenes.__makeChromSqlWhere(chroms)) for gpTable in gpTables]
        else:
            return None

    @staticmethod
    def loadFromDb(conn, gpTables, attrsTable, tagTable, errFh=None, chroms=None, transcriptIdSubSelect=None):
        """Factory method to load from a databases (see load functions). If transcriptIdSubSelect is specified, it should be a select to
        return just the transcriptIds of interest.  gpTables can be a string or a list of tables"""
        if isinstance(gpTables, str):
            gpTables = [gpTables]

        # chrom sub-select speeds up loads when doing limit testing.
        chromSubSelects = GencodeGenes.__makeChromSubSelects(gpTables, chroms)

        gencode = GencodeGenes()
        for gpTable in gpTables:
            gencode.loadGenePredsFromDb(conn, gpTable, chroms=chroms, transcriptIdSubSelect=transcriptIdSubSelect)
        gencode.loadInfoFromDb(conn, attrsTable, chromSubSelects, transcriptIdSubSelect=transcriptIdSubSelect)
        gencode.loadTagsFromDb(conn, tagTable, chromSubSelects, transcriptIdSubSelect=transcriptIdSubSelect)
        gencode.finish(errFh=errFh)
        return gencode

    @staticmethod
    def loadFromDbByVersion(conn, gencodeVersion, errFh=None, basicSet=False, pseudo=True, chroms=None, transcriptIdSubSelect=None):
        """Factory method to load from a databases assuming gencode table
        naming conventions. gencodeVersion should be in the form `VNN'.  If
        transcriptIdSubSelect is specified, it should be a select to return
        just the transcriptIds of interest."""

        gpTables = []
        if basicSet:
            gpTables.append("wgEncodeGencodeBasic" + gencodeVersion)
        else:
            gpTables.append("wgEncodeGencodeComp" + gencodeVersion)
        if pseudo:
            gpTables.append("wgEncodeGencodePseudoGene" + gencodeVersion)
        attrsTable = "wgEncodeGencodeAttrs" + gencodeVersion
        tagTable = "wgEncodeGencodeTag" + gencodeVersion

        return GencodeGenes.loadFromDb(conn, gpTables, attrsTable, tagTable, errFh=errFh, chroms=chroms, transcriptIdSubSelect=transcriptIdSubSelect)
