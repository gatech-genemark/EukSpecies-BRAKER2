#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# This script categorizes genes in gtf annotation (previously enriched with
# enrich_gff.pl script) into 4 groups: completeTranscripts,
# incompleteTranscripts, completeGenes, incompleteGenes. Complete transcript
# starts with "ATG" and ends with "TAG", "TAA", or "TGA". Complete gene has all
# its transcripts complete.
# ==============================================================


import argparse
import csv
import re
import sys


class Transcript(object):
    """docstring for Transcript"""
    def __init__(self, ID, geneID):
        self.ID = ID
        self.geneID = geneID
        self.start = ""
        self.stop = ""

    def addStart(self, row):
        seq = extractFeatureGtf(row[8], "site_seq")
        if row[6] == '+':
            self.start += seq
        else:
            self.start = seq + self.start


    def addStop(self, row):
        seq = extractFeatureGtf(row[8], "site_seq")
        if row[6] == '+':
            self.stop += seq
        else:
            self.stop = seq + self.stop

def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def loadTranscripts(gffFile):
    transcripts = {}
    for row in csv.reader(open(gffFile), delimiter='\t'):
        ID = extractFeatureGtf(row[8], "transcript_id")
        geneID = extractFeatureGtf(row[8], "gene_id")
        if ID not in transcripts:
            transcripts[ID] = Transcript(ID, geneID)

        if row[2] == "start_codon":
            transcripts[ID].addStart(row)

        if row[2] == "stop_codon":
            transcripts[ID].addStop(row)
    return transcripts


def sortTranscripts(transcripts, args):
    completeTranscripts = set()
    incompleteGenes = set()
    for key in transcripts:
        transcript = transcripts[key]
        if ((transcript.start == "ATG" or
             (transcript.start == "CTG" and args.allowCTG)) and 
           (transcript.stop == "TGA" or
            transcript.stop == "TAA" or
            transcript.stop == "TAG")):
            completeTranscripts.add(transcript.ID)
        else:
            incompleteGenes.add(transcript.geneID)
    return completeTranscripts, incompleteGenes


def printTranscripts(completeTranscripts, incompleteGenes, args):
    completeTrOut = open(args.completeTranscripts, "w")
    incompleteTrOut = open(args.incompleteTranscripts, "w")
    completeGenesOut = open(args.completeGenes, "w")
    incompleteGenesOut = open(args.incompleteGenes, "w")

    for row in csv.reader(open(args.annot), delimiter='\t'):
        ID = extractFeatureGtf(row[8], "transcript_id")
        geneID = extractFeatureGtf(row[8], "gene_id")
        if ID in completeTranscripts:
            completeTrOut.write("\t".join(row) + "\n")
        else:
            incompleteTrOut.write("\t".join(row) + "\n")

        if geneID in incompleteGenes:
            incompleteGenesOut.write("\t".join(row) + "\n")
        else:
            completeGenesOut.write("\t".join(row) + "\n")

    completeTrOut.close()
    incompleteTrOut.close()
    completeGenesOut.close()
    incompleteGenesOut.close()

def main():
    args = parseCmd()
    transcripts = loadTranscripts(args.annot)
    completeTranscripts, incompleteGenes = sortTranscripts(transcripts, args)
    printTranscripts(completeTranscripts, incompleteGenes, args)


def parseCmd():

    parser = argparse.ArgumentParser(description='This script categorizes genes in \
        gtf annotation (previously enriched with enrich_gff.pl script) into 4 groups:\
        completeTranscripts, incompleteTranscripts, completeGenes, incompleteGenes.\
        Complete transcript starts with "ATG" and ends with "TAG", "TAA", or "TGA".\
        Complete gene has all its transcripts complete.')

    parser.add_argument('annot', metavar='annot.gtf', type=str,
                        help='Annotation file in gtf format processed by \
                        enrich_gff.pl script.')

    parser.add_argument('--completeTranscripts', required=True,
                        help='Output with complete transcripts')
    
    parser.add_argument('--incompleteTranscripts', required=True,
                        help='Output with incomplete transcripts')

    parser.add_argument('--completeGenes', required=True,
                        help='Output with complete genes. All transcripts \
                        in a gene need to be complete.')

    parser.add_argument('--incompleteGenes', required=True,
                        help='Output with incomplete genes')

    parser.add_argument('--allowCTG',  default=False, action='store_true',
                        help='Allow CTG as a valid start codon.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
