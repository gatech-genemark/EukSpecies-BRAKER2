#!/usr/bin/env python
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Flag partial CDS in Ensembl annotation. Split the enriched
# annotation file into complete/incomplete genes and transcripts.
# ==============================================================


import argparse
import csv
import re


def getSignature(row):
    return row[0] + "_" + row[3] + "_" + row[4] + "_" + row[6]


def extractFeature(text, feature):
    regex = feature + ' "([^;]+)";'
    search = re.search(regex, text)
    if search:
        return search.groups()[0]
    else:
        return None


def loadOriginalAnnotation(annotFile):
    starts = set()
    stops = set()
    for row in csv.reader(open(annotFile), delimiter='\t'):
        if len(row) < 8:
            continue
        if row[2] == "start_codon":
            starts.add(getSignature(row))
        elif row[2] == "stop_codon":
            stops.add(getSignature(row))

    return starts, stops


def findIncomplete(enrichedAnnotFile, starts, stops):
    incompleteStartsTranscripts = set()
    incompleteStopsTranscripts = set()
    incompleteGenes = set()
    for row in csv.reader(open(enrichedAnnotFile), delimiter='\t'):
        if row[2] == "start_codon":
            if not getSignature(row) in starts:
                incompleteStartsTranscripts.add(extractFeature(row[8], "transcript_id"))
                incompleteGenes.add(extractFeature(row[8], "gene_id"))
        elif row[2] == "stop_codon":
            if not getSignature(row) in stops:
                incompleteStopsTranscripts.add(extractFeature(row[8], "transcript_id"))
                incompleteGenes.add(extractFeature(row[8], "gene_id"))
    return incompleteStartsTranscripts, incompleteStopsTranscripts, incompleteGenes


def saveOutputs(args, incompleteStartsTranscripts, incompleteStopsTranscripts, incompleteGenes):
    incompleteTrOut = open(args.incompleteTranscriptsOutput, "w")
    completeTrOut = open(args.completeTranscriptsOutput, "w")
    fullOut = open(args.fullOutput, "w")
    incompleteGeneOut = open(args.incompleteGenesOutput, "w")
    completeGeneOut = open(args.completeGenesOutput, "w")

    for row in csv.reader(open(args.enrichedAnnot), delimiter='\t'):
        transcript = extractFeature(row[8], "transcript_id")
        gene = extractFeature(row[8], "gene_id")
        incomplete = False
        type = extractFeature(row[8], "cds_type")

        if transcript in incompleteStartsTranscripts:
            incomplete = True
            if row[2] == "start_codon":
                continue
            elif row[2] == "CDS" and (type == "Initial" or type == "Single"):
                row[2] = "CDS_partial"

        if transcript in incompleteStopsTranscripts:
            incomplete = True
            if row[2] == "stop_codon":
                continue
            elif row[2] == "CDS" and (type == "Terminal" or type == "Single"):
                row[2] = "CDS_partial"

        if incomplete:
            incompleteTrOut.write("\t".join(row) + "\n")
        else:
            completeTrOut.write("\t".join(row) + "\n")
        fullOut.write("\t".join(row) + "\n")

        if gene in incompleteGenes:
            incompleteGeneOut.write("\t".join(row) + "\n")
        else:
            completeGeneOut.write("\t".join(row) + "\n")

    incompleteTrOut.close()
    completeTrOut.close()
    fullOut.close()
    incompleteGeneOut.close()
    completeGeneOut.close()


def main():
    args = parseCmd()
    originalStarts, originalStops = loadOriginalAnnotation(args.originalAnnot)
    incompleteStartsTranscripts, incompleteStopsTranscripts, incompleteGenes = findIncomplete(args.enrichedAnnot, originalStarts, originalStops)
    saveOutputs(args, incompleteStartsTranscripts, incompleteStopsTranscripts, incompleteGenes)


def parseCmd():

    parser = argparse.ArgumentParser(description='Flag partial CDS in Ensembl annotation.')

    parser.add_argument('originalAnnot', metavar='originalAnnot.gtf', type=str,
                        help='Original Ensembl annotation file')

    parser.add_argument('enrichedAnnot', metavar='enrichedAnnot.gtf', type=str,
                        help='Annotation processed with the enrich_gff.pl script, converted to gtf')

    parser.add_argument('--incompleteTranscriptsOutput',  required=True,
                        help='Save incomplete transcripts into this file')

    parser.add_argument('--completeTranscriptsOutput',  required=True,
                        help='Save complete transcripts into this file')

    parser.add_argument('--incompleteGenesOutput',  required=True,
                        help='Save incomplete genes (at least one transcript '
                        'in the gene is incomplete) into this file')

    parser.add_argument('--completeGenesOutput',  required=True,
                        help='Save complete genes (All transcripts in the gene'
                        ' must be complete) into this file')

    parser.add_argument('--fullOutput',  required=True,
                        help='Save both incomplete and complete transcripts into this file')

    return parser.parse_args()


if __name__ == '__main__':
    main()
