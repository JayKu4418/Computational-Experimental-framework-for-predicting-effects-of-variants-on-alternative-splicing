#!/bin/bash

shapemapper_profileFile=$1
nameFile=$2
scriptFolder="/home/jkumar/Projects/Model_MAPTsplicing/scripts/Processing_Scripts/"

# Renormalize values per nucleotide
python ${scriptFolder}/Renormalize_DMSdata_MapFileShapeFile.py -p ${shapemapper_profileFile} -w ${nameFile}_RenormalizedPerNucleotide

#Scale Gs and Us so that they are much lower in value than As and Cs
python ${scriptFolder}/ScalingGsandUsAfterNormalization.py -m ${nameFile}_RenormalizedPerNucleotide.map -w ${nameFile}_RenormalizedPerNucleotide_ScaledGUs

#convert Gs and Us to NAs so that they can be used for predicting sstructure 
python ${scriptFolder}/convertGsAndUsToLargeNegativeValuesForShapeMapperMapFile.py ${nameFile}_RenormalizedPerNucleotide_ScaledGUs.map

#Rewrite the ShapeMapper file to have different coordinates
python ${scriptFolder}/reWriteMapFile_WithDifferentCoordinates.py -m ${nameFile}_RenormalizedPerNucleotide_ScaledGUs_GsAndUsSetToLargeNegativeValue.map -s 33 -e 0 -a RewrittenCoordinatesIncludeAllExon10