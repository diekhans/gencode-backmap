#include "gxfRecord.hh"

const string GxfFeature::GENE = "gene";
const string GxfFeature::TRANSCRIPT = "transcript";
const string GxfFeature::EXON = "exon";
const string GxfFeature::CDS = "CDS";
const string GxfFeature::START_CODON = "start_codon";
const string GxfFeature::UTR = "UTR";
const string GxfFeature::STOP_CODON = "stop_codon";
const string GxfFeature::STOP_CODON_REDEFINED_AS_SELENOCYSTEINE = "stop_codon_redefined_as_selenocysteine";


// standard attribute names
const string GxfFeature::ID_ATTR = "ID";
const string GxfFeature::PARENT_ATTR = "Parent";
const string GxfFeature::GENE_ID_ATTR = "gene_id";
const string GxfFeature::GENE_NAME_ATTR = "gene_name";
const string GxfFeature::GENE_TYPE_ATTR = "gene_type";
const string GxfFeature::GENE_STATUS_ATTR = "gene_status";
const string GxfFeature::GENE_HAVANA_ATTR = "havana_gene";
const string GxfFeature::TRANSCRIPT_ID_ATTR = "transcript_id";
const string GxfFeature::TRANSCRIPT_NAME_ATTR = "transcript_name";
const string GxfFeature::TRANSCRIPT_TYPE_ATTR = "transcript_type";
const string GxfFeature::TRANSCRIPT_STATUS_ATTR = "transcript_status";
const string GxfFeature::TRANSCRIPT_HAVANA_ATTR = "havana_transcript";
const string GxfFeature::EXON_ID_ATTR = "exon_id";
const string GxfFeature::EXON_NUMBER_ATTR = "exon_number";
const string GxfFeature::TAG_ATTR = "tag";

const string GxfFeature::SOURCE_HAVANA = "HAVANA";
const string GxfFeature::SOURCE_ENSEMBL = "ENSEMBL";

/* return base columns (excluding attributes) as a string */
string GxfFeature::baseColumnsAsString() const {
    return fSeqid + "\t" + fSource + "\t" + fType + "\t" + to_string(fStart) + "\t"
        + to_string(fEnd) + "\t" + fScore + "\t" + fStrand + "\t" + fPhase + "\t";
}

/* get the id based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getTypeId() const {
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_ID_ATTR, emptyString);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_ID_ATTR, emptyString);
    } else if (fType == GxfFeature::EXON) {
        return getAttrValue(GxfFeature::EXON_ID_ATTR, emptyString);
    } else {
        return emptyString;
    }
}

/* get the havana id based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getHavanaTypeId() const {
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_HAVANA_ATTR, emptyString);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_HAVANA_ATTR, emptyString);
    } else {
        return emptyString;
    }
}

/* get the name based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getTypeName() const {
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_NAME_ATTR, emptyString);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_NAME_ATTR, emptyString);
    } else {
        return emptyString;
    }
}

/* get the biotype based on feature type, or empty string if it doesn't have an
 * id */
const string& GxfFeature::getTypeBiotype() const {
    static const string emptyString;
    if (fType == GxfFeature::GENE) {
        return getAttrValue(GxfFeature::GENE_TYPE_ATTR, emptyString);
    } else if (fType == GxfFeature::TRANSCRIPT) {
        return getAttrValue(GxfFeature::TRANSCRIPT_TYPE_ATTR, emptyString);
    } else {
        return emptyString;
    }
}
