/*
 * Tree structure use to store genes
 */
#include "featureIO.hh"
#include <iostream>
#include <algorithm>
#include "gxfIO.hh"


/* Constructor */
FeatureParser::FeatureParser(const string& gxfFile):
    fGxfParser(GxfParser::factory(gxfFile, featureFactory)),
    fNextGene(NULL) {
}

/* Destructor */
FeatureParser::~FeatureParser() {
    delete fGxfParser;
}

/* Get the next feature, skipping other records */
Feature* FeatureParser::nextFeature() {
    GxfRecord* gxfRecord;
    while ((gxfRecord = fGxfParser->next()) != NULL) {
        if (instanceOf(gxfRecord, Feature)) {
            return static_cast<Feature*>(gxfRecord);
        } else {
            delete gxfRecord;
        }
    }
    return NULL;
}

/* Get the next gene feature, possible queued */
Feature* FeatureParser::nextGeneFeature() {
    Feature* gene;
    if (fNextGene != NULL) {
        gene = fNextGene;
        fNextGene = NULL;
    } else {
        gene = nextFeature();
    }
    if ((gene != NULL) and not gene->isGene()) {
        throw invalid_argument("expected gene feature, found: " + gene->toString());
    }
    return gene;
}

/*
 * Find the parent for GFF3.
 */
Feature* FeatureParser::findGff3Parent(Feature* geneLeaf,
                                       const Feature* feature) {
    const string& parentId = feature->getAttrValue(Feature::PARENT_ATTR);
    Feature* parent = geneLeaf;
    while ((parent != NULL) and (parent->getAttrValue(Feature::ID_ATTR) != parentId)) {
        parent = parent->getParent();
    }
    if (parent == NULL) {
        throw invalid_argument("parent node " +  parentId + " for " + feature->getAttrValue(GxfFeature::ID_ATTR) + " not found");
    }
    return parent;
}

/* Get the desired type the parent of GTF feature
 * WARNING: this assumed that the hierarchy is
 * gene->transcript->{everything else}
 */
const string& FeatureParser::getGtfParentType(const string& featureType) {
    assert(featureType != GxfFeature::GENE);
    if (featureType == GxfFeature::TRANSCRIPT) {
        return GxfFeature::GENE;
    } else {
        return GxfFeature::TRANSCRIPT;
    }
}

/*
 * Find the parent for a GTF record.  This is painful guess based on the
 * GENCODE file order and know how GENCODE is structures.
 */
Feature* FeatureParser::findGtfParent(Feature* geneLeaf,
                                     const Feature* feature) {
    const string& parentType = getGtfParentType(feature->getType());
    Feature* parent = geneLeaf;
    while ((parent != NULL) and (parent->getType() != parentType)) {
        parent = parent->getParent();
    }
    if (parent == NULL) {
        throw invalid_argument("parent node of type " + parentType + "  not found for type " + parent->getType());
    }
    return parent;
}

/* find the parent */
Feature* FeatureParser::findParent(Feature* geneLeaf,
                                   const Feature* feature) {
    if (fGxfParser->getFormat() == GFF3_FORMAT) {
        return findGff3Parent(geneLeaf, feature);
    } else {
        return findGtfParent(geneLeaf, feature);
    }
}

/*
 * Process a feature for a gene.
 */
void FeatureParser::addGeneFeature(Feature* gene,
                                   Feature*& geneLeaf,
                                   Feature* feature) {
    Feature* parent = findParent(geneLeaf, feature);
    parent->addChild(feature);
    geneLeaf = feature;
}

/*
 * Load all child records associated with a given gene.  This whole thing is
 * annoying due to the lack of explicit structure in GTF.
 */
void FeatureParser::loadGeneChildren(Feature* gene) {
    assert(gene->getType() == GxfFeature::GENE);

    Feature* geneLeaf = gene;  // were we are currently working
    Feature* feature;
    while ((feature = nextFeature()) != NULL) {
        if (feature->isGene()) {
            fNextGene = feature;
            break;
        } else {
            addGeneFeature(gene, geneLeaf, feature);
        }
    }
}

/*
 * Remove transcript attributes that were accidentally left on genes
 * in some gencode versions.
 */
void FeatureParser::removeTransAttrsOnGenes(Feature* gene) {
    AttrVals& attrs = gene->getAttrs();
    attrs.remove(GxfFeature::TRANSCRIPT_ID_ATTR);
    attrs.remove(GxfFeature::TRANSCRIPT_TYPE_ATTR);
    attrs.remove(GxfFeature::TRANSCRIPT_STATUS_ATTR);
    attrs.remove(GxfFeature::TRANSCRIPT_NAME_ATTR);
}

/* factory function to create a feature record */
GxfFeature* FeatureParser::featureFactory(const string& seqid, const string& source, const string& type,
                                          int start, int end, const string& score, const string& strand,
                                          const string& phase, const AttrVals& attrs) {
    Feature* feature = new Feature(seqid, source, type, start, end, score, strand, phase, attrs);
    if (type == Feature::GENE) {
        removeTransAttrsOnGenes(feature);
    }
    return feature;
}
    
/* load next gene */
Feature* FeatureParser::nextGene() {
    Feature* gene = nextGeneFeature();
    if (gene != NULL) {
        loadGeneChildren(gene);
    }
    return gene;;
}

