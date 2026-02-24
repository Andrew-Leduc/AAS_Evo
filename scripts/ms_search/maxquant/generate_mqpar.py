# -*- coding: utf-8 -*-
"""
generate_mqpar.py

Generate MaxQuant mqpar.xml parameter files for per-plex TMT-11 MS searches.

MaxQuant handles decoys internally (reversed sequences, decoyMode=revert) —
do NOT add REV_ entries to the FASTA. With minPeptides=1, single-peptide
protein identifications are reported, which is critical for our mutant entries.

Results land in {output_dir}/{plex_id}/combined/
RAW files are symlinked into results/{plex_id}/ so MaxQuant places combined/
there automatically (MaxQuant creates combined/ in the parent of the RAW files).

Usage:
    python3 generate_mqpar.py \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        --raw-dir /path/to/RAW/ \
        --fasta-dir /path/to/FASTA/per_plex_mq/ \
        -o /path/to/MQ_SEARCH/ \
        [--threads 16]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from xml.etree import ElementTree as ET
from xml.dom import minidom


TMT_CHANNEL_MAP = {
    "tmt_126":  "126C",
    "tmt_127n": "127N",
    "tmt_127c": "127C",
    "tmt_128n": "128N",
    "tmt_128c": "128C",
    "tmt_129n": "129N",
    "tmt_129c": "129C",
    "tmt_130n": "130N",
    "tmt_130c": "130C",
    "tmt_131":  "131N",
    "tmt_131c": "131C",
}

# MaxQuant internal modification names for TMT11plex.
# Names verified against modifications.xml in MaxQuant 2.7.5.0.
# Key naming rules:
#   - Channel 126:  no N/C suffix  → "TMT10plex-Lys126"  / "TMT10plex-Nter126"
#   - Channels 127–130: N/C suffix → "TMT10plex-Lys127N" / "TMT10plex-Nter127N" etc.
#   - Channel 131N: no N/C suffix  → "TMT10plex-Lys131"  / "TMT10plex-Nter131"
#   - Channel 131C: TMT11plex      → "TMT11plex-Lys131C" / "TMT11plex-Nter131C"
# Tuple: (internalLabel for Lys, terminalLabel for N-term, channel_key)
TMT11_CHANNELS = [
    ("TMT10plex-Lys126",   "TMT10plex-Nter126",   "126C"),
    ("TMT10plex-Lys127N",  "TMT10plex-Nter127N",  "127N"),
    ("TMT10plex-Lys127C",  "TMT10plex-Nter127C",  "127C"),
    ("TMT10plex-Lys128N",  "TMT10plex-Nter128N",  "128N"),
    ("TMT10plex-Lys128C",  "TMT10plex-Nter128C",  "128C"),
    ("TMT10plex-Lys129N",  "TMT10plex-Nter129N",  "129N"),
    ("TMT10plex-Lys129C",  "TMT10plex-Nter129C",  "129C"),
    ("TMT10plex-Lys130N",  "TMT10plex-Nter130N",  "130N"),
    ("TMT10plex-Lys130C",  "TMT10plex-Nter130C",  "130C"),
    ("TMT10plex-Lys131",   "TMT10plex-Nter131",   "131N"),
    ("TMT11plex-Lys131C",  "TMT11plex-Nter131C",  "131C"),
]


def load_tmt_map(tmt_map_path):
    plex_to_files    = defaultdict(set)
    plex_to_channels = defaultdict(list)
    with open(tmt_map_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run_id      = row["run_metadata_id"].strip()
            file_name   = row["file_name"].strip()
            channel     = row.get("tmt_channel", "").strip().lower()
            case_id     = row.get("case_submitter_id", "").strip()
            sample_type = row.get("sample_type", "").strip()
            plex_to_files[run_id].add(file_name)
            plex_to_channels[run_id].append({
                "file_name":   file_name,
                "channel":     channel,
                "case_id":     case_id,
                "sample_type": sample_type,
            })
    return plex_to_files, plex_to_channels


def sub(parent, tag, text=None):
    el = ET.SubElement(parent, tag)
    if text is not None:
        el.text = str(text)
    return el


def build_mqpar(raw_files, fasta_path, out_dir, plex_id,
                channel_entries, n_threads):
    """
    Build and return a MaxQuant mqpar.xml ElementTree for one plex.

    Structure mirrors sample_mqpar.txt (MaxQuant 2.7.5.0 GUI-generated).

    Custom choices vs sample:
      - identifierParseRule: >([^\\s]*)   works for sp|/tr|/mut|/comp| headers
        (sample uses >[^|]*\\|(.*?)\\| which only captures UniProt accession)
      - taxonomyParseRule: ""             our headers have no OX= field
      - RAW files: symlinked into results/{plex_id}/ so combined/ lands there
      - Carbamidomethyl (C): kept as fixed mod (matches sample_mqpar)
      - writeAllPeptidesTable: True       needed for SLURM completion detection
    """
    root = ET.Element("MaxQuantParams",
                      attrib={
                          "xmlns:xsd": "http://www.w3.org/2001/XMLSchema",
                          "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
                      })

    # ── FASTA ─────────────────────────────────────────────────────────────────
    fastas = sub(root, "fastaFiles")
    fi = sub(fastas, "FastaFileInfo")
    sub(fi, "fastaFilePath",         fasta_path)
    sub(fi, "identifierParseRule",   r">([^\s]*)")   # custom: handles mut|/comp| prefixes
    sub(fi, "descriptionParseRule",  r">(.*)")
    sub(fi, "taxonomyParseRule",     "")              # no OX= in our headers
    sub(fi, "variationParseRule",    "")
    sub(fi, "modificationParseRule", "")
    sub(fi, "taxonomyId",            "")

    sub(root, "fastaFilesProteogenomics")  # empty — required
    sub(root, "fastaFilesFirstSearch")     # empty — required
    sub(root, "fixedSearchFolder",         "")

    # ── GLOBAL SETTINGS (field order matches sample_mqpar.txt / MQ 2.7.5.0) ──
    sub(root, "andromedaCacheSize",          350000)
    sub(root, "advancedRatios",              "True")
    sub(root, "pvalThres",                   0.005)
    sub(root, "rtShift",                     "False")
    sub(root, "separateLfq",                 "False")
    sub(root, "lfqStabilizeLargeRatios",     "True")
    sub(root, "lfqRequireMsms",              "True")
    sub(root, "lfqBayesQuant",               "False")
    sub(root, "decoyMode",                   "revert")  # MQ generates reversed decoys
    sub(root, "includeContaminants",         "True")
    sub(root, "maxPeptideMass",              4600)
    sub(root, "epsilonMutationScore",        "True")
    sub(root, "mutatedPeptidesSeparately",   "True")
    sub(root, "proteogenomicPeptidesSeparately", "True")
    sub(root, "minDeltaScoreUnmodifiedPeptides", 0)
    sub(root, "minDeltaScoreModifiedPeptides",   6)
    sub(root, "minScoreUnmodifiedPeptides",      0)
    sub(root, "minScoreModifiedPeptides",        40)
    sub(root, "secondPeptide",               "True")
    sub(root, "matchBetweenRuns",            "False")
    sub(root, "matchUnidentifiedFeatures",   "False")
    sub(root, "matchBetweenRunsFdr",         "False")
    sub(root, "dependentPeptides",           "False")
    sub(root, "dependentPeptideFdr",         0)
    sub(root, "dependentPeptideMassBin",     0)
    sub(root, "dependentPeptidesBetweenRuns",              "False")
    sub(root, "dependentPeptidesWithinExperiment",         "False")
    sub(root, "dependentPeptidesWithinParameterGroup",     "False")
    sub(root, "dependentPeptidesRestrictFractions",        "False")
    sub(root, "dependentPeptidesFractionDifference",       0)
    sub(root, "ibaq",                        "False")
    sub(root, "top3",                        "False")
    sub(root, "independentEnzymes",          "False")
    sub(root, "useDeltaScore",               "False")
    sub(root, "splitProteinGroupsByTaxonomy","False")
    sub(root, "taxonomyLevel",               "Species")
    sub(root, "avalon",                      "False")
    sub(root, "nModColumns",                 3)
    sub(root, "ibaqLogFit",                  "False")
    sub(root, "ibaqChargeNormalization",     "False")
    sub(root, "razorProteinFdr",             "True")
    sub(root, "deNovoSequencing",            "False")
    sub(root, "deNovoVarMods",               "False")
    sub(root, "deNovoCompleteSequence",      "False")
    sub(root, "deNovoCalibratedMasses",      "False")
    sub(root, "deNovoMaxIterations",         0)
    sub(root, "deNovoProteaseReward",        0)
    sub(root, "deNovoProteaseRewardTof",     0)
    sub(root, "deNovoAgPenalty",             0)
    sub(root, "deNovoGgPenalty",             0)
    sub(root, "deNovoUseComplementScore",    "True")
    sub(root, "deNovoUseProteaseScore",      "True")
    sub(root, "deNovoUseWaterLossScore",     "True")
    sub(root, "deNovoUseAmmoniaLossScore",   "True")
    sub(root, "deNovoUseA2Score",            "True")
    sub(root, "deNovoMassClusterTolDa",      0)
    sub(root, "deNovoScalingFactor",         0)
    sub(root, "massDifferenceSearch",        "False")
    sub(root, "isotopeCalc",                 "False")
    sub(root, "minPeptideLength",            7)       # ← correct 2.7.5.0 field name
    sub(root, "psmFdrCrosslink",             0.01)
    sub(root, "peptideFdr",                  0.01)
    sub(root, "proteinFdr",                  0.01)
    sub(root, "siteFdr",                     0.01)
    sub(root, "minPeptideLengthForUnspecificSearch", 8)
    sub(root, "maxPeptideLengthForUnspecificSearch", 25)
    sub(root, "useNormRatiosForOccupancy",   "True")
    sub(root, "minPeptides",                 1)    # ← key: single-peptide proteins reported
    sub(root, "minRazorPeptides",            1)
    sub(root, "minUniquePeptides",           0)
    sub(root, "useCounterparts",             "False")
    sub(root, "advancedSiteIntensities",     "True")
    sub(root, "minRatioCount",               2)
    sub(root, "restrictProteinQuantification", "True")
    restrict_mods = sub(root, "restrictMods")
    sub(restrict_mods, "string", "Oxidation (M)")
    sub(restrict_mods, "string", "Acetyl (Protein N-term)")
    sub(root, "matchingTimeWindow",          0)
    sub(root, "matchingIonMobilityWindow",   0)
    sub(root, "alignmentTimeWindow",         0)
    sub(root, "alignmentIonMobilityWindow",  0)
    sub(root, "numberOfCandidatesMsms",      15)
    sub(root, "compositionPrediction",       0)
    sub(root, "quantMode",                   1)    # 1 = reporter ion (TMT)
    sub(root, "massDifferenceMods")                 # empty
    sub(root, "mainSearchMaxCombinations",   200)
    sub(root, "writeMsScansTable",           "False")
    sub(root, "writeMsmsScansTable",         "True")
    sub(root, "writePasefMsmsScansTable",    "True")
    sub(root, "writeAccumulatedMsmsScansTable", "True")
    sub(root, "writeMs3ScansTable",          "True")
    sub(root, "writeAllPeptidesTable",       "True")  # needed for SLURM completion check
    sub(root, "writeMzRangeTable",           "True")
    sub(root, "writeDiaFragmentTable",       "False")
    sub(root, "writeDiaFragmentQuantTable",  "False")
    sub(root, "writeMzTab",                  "False")
    sub(root, "disableMd5",                  "False")
    sub(root, "cacheBinInds",                "True")
    sub(root, "etdIncludeB",                 "False")
    sub(root, "ms2PrecursorShift",           0)
    sub(root, "complementaryIonPpm",         20)
    sub(root, "variationParseRule",          "")
    sub(root, "variationMode",               "none")
    sub(root, "useSeriesReporters",          "False")
    sub(root, "name",                        plex_id)
    sub(root, "maxQuantVersion",             "2.7.5.0")
    sub(root, "pluginFolder",                "")
    sub(root, "numThreads",                  n_threads)
    sub(root, "emailAddress",                "")
    sub(root, "smtpHost",                    "")
    sub(root, "emailFromAddress",            "")
    # fixedCombinedFolder removed in MQ 2.7.5.0; output dir controlled by symlinks
    sub(root, "fullMinMz",   -1.7976931348623157e+308)
    sub(root, "fullMaxMz",    1.7976931348623157e+308)
    sub(root, "sendEmail",                   "False")
    sub(root, "ionCountIntensities",         "False")
    sub(root, "verboseColumnHeaders",        "False")
    sub(root, "calcPeakProperties",          "False")
    sub(root, "showCentroidMassDifferences", "False")
    sub(root, "showIsotopeMassDifferences",  "False")
    # useDotNetCore removed in MQ 2.7.5.0
    sub(root, "profilePerformance",          "False")

    # ── RAW FILES ─────────────────────────────────────────────────────────────
    fp_el   = sub(root, "filePaths")
    exp_el  = sub(root, "experiments")
    frac_el = sub(root, "fractions")
    ptm_el  = sub(root, "ptms")
    pg_el   = sub(root, "paramGroupIndices")

    for raw in sorted(raw_files):
        basename_no_ext = os.path.splitext(os.path.basename(raw))[0]
        sub(fp_el,   "string",  raw)
        sub(exp_el,  "string",  basename_no_ext)  # per-file experiment name (matches MQ CLI template)
        sub(frac_el, "short",   32767)             # 32767 = no fractionation (matches MQ CLI template)
        sub(ptm_el,  "boolean", "False")
        sub(pg_el,   "int",     0)

    ref_ch = sub(root, "referenceChannel")
    for _ in sorted(raw_files):
        sub(ref_ch, "string", "")

    # Post-referenceChannel global fields
    sub(root, "lfqTopNPeptides",             0)
    sub(root, "diaJoinPrecChargesForLfq",    "False")
    sub(root, "diaFragChargesForQuant",      1)
    sub(root, "gridSpacing",                 0.1)
    sub(root, "proteinGroupingFile",         "")
    sub(root, "simplePepCalculation",        "False")
    sub(root, "useAndromeda20",              "False")
    sub(root, "useAndromeda20DefaultModel",  "False")
    sub(root, "andromeda20AltModelPath",     "")
    sub(root, "intensityPredictionFolder",   "")
    sub(root, "encoding",                    0)
    sub(root, "customTxtFolder",             "")

    # ── PARAMETER GROUP ───────────────────────────────────────────────────────
    pgs = sub(root, "parameterGroups")
    pg  = sub(pgs,  "parameterGroup")

    sub(pg, "andromeda20AltModelPath",   "")
    sub(pg, "andromeda20DefaultModel",   "False")
    sub(pg, "useAndromeda20",            "False")
    sub(pg, "msInstrument",              0)     # 0 = Orbitrap
    sub(pg, "maxCharge",                 7)
    sub(pg, "minPeakLen",                2)
    sub(pg, "diaMinPeakLen",             1)
    sub(pg, "useMs1Centroids",           "False")
    sub(pg, "useMs2Centroids",           "False")
    sub(pg, "cutPeaks",                  "True")
    sub(pg, "gapScans",                  1)
    sub(pg, "minTime",                   "NaN")
    sub(pg, "maxTime",                   "NaN")
    sub(pg, "matchType",                 "MatchFromAndTo")
    sub(pg, "intensityDetermination",    0)
    sub(pg, "centroidMatchTol",          8)
    sub(pg, "centroidMatchTolInPpm",     "True")
    sub(pg, "centroidHalfWidth",         35)
    sub(pg, "centroidHalfWidthInPpm",    "True")
    sub(pg, "valleyFactor",              1.4)
    sub(pg, "isotopeValleyFactor",       1.2)
    sub(pg, "advancedPeakSplitting",     "False")
    sub(pg, "intensityThresholdMs1Dda",  0)
    sub(pg, "intensityThresholdMs1Dia",  0)
    sub(pg, "intensityThresholdMs2",     0)

    label_mods = sub(pg, "labelMods")
    sub(label_mods, "string", "")    # one empty string child — required

    sub(pg, "lcmsRunType",               "Reporter MS2")
    sub(pg, "reQuantify",                "False")
    sub(pg, "lfqMode",                   0)
    sub(pg, "lfqNormClusterSize",        80)
    sub(pg, "lfqMinEdgesPerNode",        3)
    sub(pg, "lfqAvEdgesPerNode",         6)
    sub(pg, "lfqMaxFeatures",            100000)
    sub(pg, "neucodeMaxPpm",             0)
    sub(pg, "neucodeResolution",         0)
    sub(pg, "neucodeResolutionInMda",    "False")
    sub(pg, "neucodeInSilicoLowRes",     "False")
    sub(pg, "fastLfq",                   "True")
    sub(pg, "lfqRestrictFeatures",       "False")
    sub(pg, "lfqMinRatioCount",          1)
    sub(pg, "lfqMinRatioCountDia",       2)
    sub(pg, "lfqPrioritizeMs1Dia",       "True")
    sub(pg, "lfqLamdaLimit",             20)
    sub(pg, "lfqMaxIter",                10)
    sub(pg, "maxLabeledAa",              0)
    sub(pg, "maxNmods",                  5)
    sub(pg, "maxMissedCleavages",        2)
    sub(pg, "maxMissedCleavagesDia",     1)
    sub(pg, "multiplicity",              1)
    sub(pg, "enzymeMode",                0)
    sub(pg, "complementaryReporterType", 0)
    sub(pg, "reporterNormalization",     0)
    sub(pg, "neucodeIntensityMode",      0)

    fixed = sub(pg, "fixedModifications")
    sub(fixed, "string", "Carbamidomethyl (C)")    # standard for IAA-alkylated TMT

    enzymes = sub(pg, "enzymes")
    sub(enzymes, "string", "Trypsin/P")
    sub(pg, "enzymesFirstSearch")
    sub(pg, "enzymeModeFirstSearch",               0)
    sub(pg, "useEnzymeFirstSearch",                "False")
    sub(pg, "useVariableModificationsFirstSearch",  "False")

    var = sub(pg, "variableModifications")
    sub(var, "string", "Oxidation (M)")
    sub(var, "string", "Acetyl (Protein N-term)")

    sub(pg, "useMultiModification",    "False")
    sub(pg, "multiModifications")      # empty

    # TMT11 isobaric labels — verified against MQ 2.7.5.0 modifications.xml
    iso_el = sub(pg, "isobaricLabels")
    for internal_lys, terminal_nter, _ch_key in TMT11_CHANNELS:
        ili = sub(iso_el, "IsobaricLabelInfo")
        sub(ili, "internalLabel",      internal_lys)
        sub(ili, "terminalLabel",      terminal_nter)
        sub(ili, "correctionFactorM2", 0)
        sub(ili, "correctionFactorM1", 0)
        sub(ili, "correctionFactorP1", 0)
        sub(ili, "correctionFactorP2", 0)
        sub(ili, "tmtLike",            "True")

    sub(pg, "neucodeLabels")
    sub(pg, "variableModificationsFirstSearch")
    sub(pg, "hasAdditionalVariableModifications",  "False")
    sub(pg, "additionalVariableModifications")
    sub(pg, "additionalVariableModificationProteins")

    sub(pg, "doMassFiltering",             "True")
    sub(pg, "firstSearchTol",              20)     # MS1 first pass (ppm)
    sub(pg, "mainSearchTol",               4.5)    # MS1 main (ppm)
    sub(pg, "searchTolInPpm",              "True")
    sub(pg, "isotopeMatchTol",             2)
    sub(pg, "isotopeMatchTolInPpm",        "True")
    sub(pg, "isotopeTimeCorrelation",      0.6)
    sub(pg, "theorIsotopeCorrelation",     0.6)
    sub(pg, "checkMassDeficit",            "True")
    sub(pg, "recalibrationInPpm",          "True")
    sub(pg, "intensityDependentCalibration","False")
    sub(pg, "minScoreForCalibration",      70)
    sub(pg, "matchLibraryFile",            "False")
    sub(pg, "libraryFile",                 "")
    sub(pg, "matchLibraryMassTolPpm",      0)
    sub(pg, "matchLibraryTimeTolMin",      0)
    sub(pg, "matchLabelTimeTolMin",        0)
    sub(pg, "reporterMassTolerance",       0.003)  # 3 mDa — high-res MS2
    sub(pg, "reporterPif",                 0)
    sub(pg, "filterPif",                   "False")
    sub(pg, "reporterFraction",            0)
    sub(pg, "reporterBasePeakRatio",       0)

    # TIMS settings — all zero/default for Orbitrap data
    sub(pg, "timsHalfWidth",               0)
    sub(pg, "timsStep",                    0)
    sub(pg, "timsResolution",              0)
    sub(pg, "timsMinMsmsIntensity",        0)
    sub(pg, "timsRemovePrecursor",         "True")
    sub(pg, "timsIsobaricLabels",          "False")
    sub(pg, "timsCollapseMsms",            "True")

    # Cross-linking settings — disabled
    sub(pg, "crossLinkingType",            0)
    sub(pg, "crossLinker",                 "")
    sub(pg, "minMatchXl",                  3)
    sub(pg, "minPairedPepLenXl",           6)
    sub(pg, "minScoreDipeptide",           40)
    sub(pg, "minScoreMonopeptide",         0)
    sub(pg, "minScorePartialCross",        10)
    sub(pg, "crosslinkOnlyIntraProtein",   "False")
    sub(pg, "crosslinkIntensityBasedPrecursor", "True")
    sub(pg, "isHybridPrecDetermination",   "False")
    sub(pg, "topXcross",                   3)
    sub(pg, "doesSeparateInterIntraProteinCross", "False")
    sub(pg, "crosslinkMaxMonoUnsaturated", 0)
    sub(pg, "crosslinkMaxMonoSaturated",   0)
    sub(pg, "crosslinkMaxDiUnsaturated",   0)
    sub(pg, "crosslinkMaxDiSaturated",     0)
    sub(pg, "crosslinkModifications")
    sub(pg, "crosslinkFastaFiles")
    sub(pg, "crosslinkSites")
    sub(pg, "crosslinkNetworkFiles")
    sub(pg, "crosslinkMode",               "")

    sub(pg, "peakRefinement",              "False")
    sub(pg, "peakRefinementCrosslinking",  "False")
    sub(pg, "isobaricSumOverWindow",       "True")
    sub(pg, "isobaricWeightExponent",      "0.75")
    sub(pg, "collapseMsmsOnIsotopePatterns", "False")

    # DIA settings — all defaults (we run DDA/Reporter MS2)
    sub(pg, "diaLibraryType",              0)
    sub(pg, "diaLibraryPaths")
    sub(pg, "diaEvidencePaths")
    sub(pg, "diaMsmsPaths")
    sub(pg, "diaLabelIndsForLibraryMatch")
    sub(pg, "diaInitialPrecMassTolPpm",    10)
    sub(pg, "diaInitialFragMassTolPpm",    20)
    sub(pg, "diaCorrThresholdFeatureClustering",  0.85)
    sub(pg, "diaPrecTolPpmFeatureClustering",     2)
    sub(pg, "diaFragTolPpmFeatureClustering",     2)
    sub(pg, "diaScoreN",                   12)
    sub(pg, "diaScoreNAdditional",         5)
    sub(pg, "diaMinScore",                 1.99)
    sub(pg, "diaXgBoostBaseScore",         0.4)
    sub(pg, "diaXgBoostSubSample",         0.9)
    sub(pg, "centroidPosition",            0)
    sub(pg, "diaQuantMethod",              7)
    sub(pg, "diaFeatureQuantMethod",       0)
    sub(pg, "lfqNormType",                 1)
    sub(pg, "diaTopNForQuant",             6)
    sub(pg, "diaTopNCorrelatingFragmentsForQuant", 6)
    sub(pg, "diaFragmentCorrelationForQuant",      0.78)
    sub(pg, "lfqUsePeptideCorrelation",    "True")
    sub(pg, "lfqTopNCorrelatingPeptides",  100)
    sub(pg, "lfqPeptideCorrelation",       0)
    sub(pg, "diaMinMsmsIntensityForQuant", 0)
    sub(pg, "diaTopMsmsIntensityQuantileForQuant", 0.85)
    sub(pg, "diaMinFragmentOverlapScore",  0)
    sub(pg, "diaMinPrecursorScore",        0)
    sub(pg, "diaUseProfileCorrelation",    "False")
    sub(pg, "diaMinPrecProfileCorrelation",0)
    sub(pg, "diaMinFragProfileCorrelation",0)
    sub(pg, "diaClassificationModel",      1)
    sub(pg, "diaXgBoostMinChildWeight",    9)
    sub(pg, "diaXgBoostMaximumTreeDepth",  4)
    sub(pg, "diaXgBoostEstimators",        580)
    sub(pg, "diaXgBoostGamma",             0.9)
    sub(pg, "diaXgBoostMaxDeltaStep",      3)
    sub(pg, "diaGlobalMl",                 "True")
    sub(pg, "diaAdaptiveMassAccuracy",     "False")
    sub(pg, "diaMassWindowFactor",         3.3)
    sub(pg, "diaPermuteRt",               "False")
    sub(pg, "diaPermuteCcs",              "False")
    sub(pg, "diaBackgroundSubtraction",   "False")
    sub(pg, "diaBackgroundSubtractionQuantile", 0.5)
    sub(pg, "diaBackgroundSubtractionFactor",   4)
    sub(pg, "diaLfqRatioType",             0)
    sub(pg, "diaTransferQvalue",           0.9)
    sub(pg, "diaTransferQvalueBetweenLabels",    0.01)
    sub(pg, "diaTransferQvalueBetweenFractions", 0.01)
    sub(pg, "diaTransferQvalueBetweenFaims",     0.01)
    sub(pg, "diaOnlyIsosForRecal",         "True")
    sub(pg, "diaMinPeaks",                 5)
    sub(pg, "diaUseFragIntensForMl",       "False")
    sub(pg, "diaUseFragMassesForMl",       "False")
    sub(pg, "diaMaxTrainInstances",        500000)
    sub(pg, "diaMaxFragmentCharge",        3)
    sub(pg, "diaAdaptiveMlScoring",        "False")
    sub(pg, "diaDynamicScoringMaxInstances", 25000)
    sub(pg, "diaMaxPrecursorMz",           0)
    sub(pg, "diaHardRtFilter",             "True")
    sub(pg, "diaHardCcsFilter",            "True")
    sub(pg, "diaConvertLibraryCharge2Fragments", "False")
    sub(pg, "diaChargeNormalizationLibrary",     "True")
    sub(pg, "diaChargeNormalizationSample",      "True")
    sub(pg, "diaDeleteIntermediateResults",      "True")
    sub(pg, "diaDeleteRawFilesAfterSearch",      "False")
    sub(pg, "diaScoreWeightScanIndex",     -1)
    sub(pg, "diaScoreWeightScanValue",     1)
    sub(pg, "diaNumNonleadingMatches",     1)
    sub(pg, "diaFragmentModelType",        0)
    sub(pg, "diaAltFragmentModelPath",     "")
    sub(pg, "diaUseDefaultRtModel",        "True")
    sub(pg, "diaAltRtModelPath",           "")
    sub(pg, "diaUseDefaultCcsModel",       "True")
    sub(pg, "diaAltCcsModelPath",          "")
    sub(pg, "diaBatchProcessing",          "True")
    sub(pg, "diaBatchSize",                0)
    sub(pg, "diaDecoyMode",                0)
    sub(pg, "diaFirstBatch",               -1)
    sub(pg, "diaLastBatch",                -1)
    sub(pg, "diaOnlyPreprocess",           "False")
    sub(pg, "diaMultiplexQuantMethod",     1)
    sub(pg, "diaOnlyPostprocess",          "False")
    sub(pg, "diaSkipLibraryGeneration",    "False")
    sub(pg, "diaRequirePrecursor",         "False")
    sub(pg, "diaFuturePeptides",           "False")
    sub(pg, "diaOverrideRtWithPrediction", "False")
    sub(pg, "diaMaxModifications",         2)
    sub(pg, "diaMaxPositionings",          20)
    sub(pg, "diaPredictCharges",           "False")
    sub(pg, "diaUseProbScore",             "True")
    sub(pg, "diaProbScoreP",               0.055)
    sub(pg, "diaProbScoreG",               0.5)
    sub(pg, "diaProbScoreStep",            0.1)
    sub(pg, "isM2FragTypeOverride",        "False")
    sub(pg, "ms2FragTypeOverride",         0)
    sub(pg, "classicLfqForSingleShots",    "True")
    sub(pg, "subsamplePeptidesForLfqNorm", "True")
    sub(pg, "lfqSubsampleNormConstant",    400)
    sub(pg, "sequenceBasedModifier",       "False")
    sub(pg, "diaRtFromSamplesForExport",   "True")
    sub(pg, "diaCcsFromSamplesForExport",  "True")
    sub(pg, "diaLibraryExport",            0)
    sub(pg, "diaUseApexWeightsForPtmLoc",  "False")
    sub(pg, "diaSecondScoreForMultiplex",   "True")
    sub(pg, "diaBestPrecursorChargeForQuant",      "False")
    sub(pg, "diaBestPrecursorMassRangeForQuant",   "False")
    sub(pg, "diaBestPrecursorIntensityForQuant",   "False")
    sub(pg, "diaPtmLocMethod",             1)
    sub(pg, "diaSecondScoring",            "True")
    sub(pg, "diaIndexedSearch",            "True")
    sub(pg, "diaIncreaseSequenceCoverage", "True")
    sub(pg, "diaSequenceCoverageFdr",      0.01)
    sub(pg, "diaHasIntensitiesForMl",      "True")
    sub(pg, "diaFdrCvIncludeCharge",       "True")
    sub(pg, "diaMinCharge",                1)
    sub(pg, "diaMaxCharge",                4)
    sub(pg, "diaDeltaMzLower",             0.55)
    sub(pg, "diaDeltaMzUpper",             0.05)
    sub(pg, "diaUseNeutralLossesForQuant", "True")
    sub(pg, "diaSecondBestCorrelationThreshold", 0)
    sub(pg, "diaUseSecondMatchAsDecoy",    "False")
    sub(pg, "mobilityPredictionTwoConformations", "False")
    sub(pg, "diaHasLibrarySupplement",     "False")
    sub(pg, "diaSupplementalFile",         "")
    sub(pg, "diaSupplementalEvidenceFile", "")
    sub(pg, "diaNeuralNetNumIter",         0)
    sub(pg, "diaNeuralNetDropout",         0)
    sub(pg, "diaNeuralNetHiddenLayers",    0)
    sub(pg, "diaNeuralNetBatchSize",       0)
    sub(pg, "diaUseRnnRtModel",            "True")
    sub(pg, "diaUseRnnCcsModel",           "False")
    sub(pg, "diaUseCharge2Fragments",      "True")
    sub(pg, "diaUseWaterLosses",           "False")
    sub(pg, "diaUseAmmoniaLosses",         "False")
    sub(pg, "diaUseRtFeaturesInMlScoring", "True")
    sub(pg, "diaCountSamples",             "False")
    sub(pg, "diaCountSamplesSecondScore",  "True")
    sub(pg, "diaCountSamplesUpperPercentile",    50)
    sub(pg, "diaCountSamplesNumLowerPercentiles", 5)
    sub(pg, "diaSeparateCalibrationPerLib","False")

    # ── MSMS PARAMS ARRAY ─────────────────────────────────────────────────────
    msms_arr = sub(root, "msmsParamsArray")

    _msms_configs = [
        # (Name,      MatchTol, MatchPpm,  DeisoTol, DeisoPpm,  DeNovoTol, DeNovoPpm, Deiso,   Topx)
        ("FTMS",     20,    "True",  7,     "True",  25,   "True",  "True",  12),
        ("ITMS",     0.5,   "False", 0.15,  "False", 0.5,  "False", "False", 8),
        ("TOF",      25,    "True",  0.01,  "False", 25,   "True",  "True",  16),
        ("ASTRAL",   25,    "True",  0.01,  "False", 25,   "True",  "True",  16),
        ("UNKNOWN",  20,    "True",  7,     "True",  25,   "True",  "True",  12),
    ]
    for (name, match_tol, match_ppm, deiso_tol, deiso_ppm,
         denovo_tol, denovo_ppm, deiso, topx) in _msms_configs:
        mp = sub(msms_arr, "msmsParams")
        sub(mp, "Name",                    name)
        sub(mp, "MatchTolerance",          match_tol)
        sub(mp, "MatchToleranceInPpm",     match_ppm)
        sub(mp, "DeisotopeTolerance",      deiso_tol)
        sub(mp, "DeisotopeToleranceInPpm", deiso_ppm)
        sub(mp, "DeNovoTolerance",         denovo_tol)
        sub(mp, "DeNovoToleranceInPpm",    denovo_ppm)
        sub(mp, "Deisotope",               deiso)
        sub(mp, "Topx",                    topx)
        sub(mp, "TopxInterval",            100)
        sub(mp, "HigherCharges",           "True")
        sub(mp, "IncludeWater",            "True")
        sub(mp, "IncludeAmmonia",          "True")
        sub(mp, "IncludeWaterCross",       "False")
        sub(mp, "IncludeAmmoniaCross",     "False")
        sub(mp, "DependentLosses",         "True")
        sub(mp, "Recalibration",           "False")

    # ── FRAGMENTATION PARAMS ARRAY ────────────────────────────────────────────
    frag_arr = sub(root, "fragmentationParamsArray")
    for frag_name in ("CID", "HCD", "ETD", "PQD", "ETHCD", "ETCID", "UVPD", "Unknown"):
        fp = sub(frag_arr, "fragmentationParams")
        sub(fp, "Name",                     frag_name)
        sub(fp, "UseIntensityPrediction",   "False")
        sub(fp, "UseSequenceBasedModifier", "False")
        sub(fp, "InternalFragments",        "False")
        sub(fp, "InternalFragmentWeight",   1)
        sub(fp, "InternalFragmentAas",      "KRH")

    return root


def prettify(root):
    raw = ET.tostring(root, encoding="unicode")
    parsed = minidom.parseString(raw)
    return parsed.toprettyxml(indent="  ", encoding=None)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MaxQuant mqpar.xml files for per-plex TMT searches."
    )
    parser.add_argument("--tmt-map", required=True,
                        help="TMT mapping file (pdc_file_tmt_map.tsv)")
    parser.add_argument("--raw-dir", required=True,
                        help="Directory containing RAW files")
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory with per-plex FASTAs (from combine_fastas.py)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for MQ_SEARCH setup")
    parser.add_argument("--threads", type=int, default=16,
                        help="CPU threads per search (default: 16)")
    args = parser.parse_args()

    for path, label in [(args.tmt_map,  "TMT map"),
                        (args.raw_dir,  "RAW directory"),
                        (args.fasta_dir,"FASTA directory")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    mqpar_dir      = os.path.join(args.output, "mqpar")
    plex_list_path = os.path.join(args.output, "plex_list.txt")
    os.makedirs(mqpar_dir, exist_ok=True)

    print(f"Loading TMT map: {args.tmt_map}")
    plex_to_files, plex_to_channels = load_tmt_map(args.tmt_map)
    print(f"  {len(plex_to_files)} plexes")

    available_raws = {f for f in os.listdir(args.raw_dir)
                      if f.lower().endswith(".raw")}
    available_fastas = {f.replace(".fasta", "")
                        for f in os.listdir(args.fasta_dir)
                        if f.endswith(".fasta")}
    print(f"  {len(available_raws)} RAW files")
    print(f"  {len(available_fastas)} per-plex FASTAs")

    plex_list = []
    n_skip_fasta = n_skip_raw = 0

    for plex_id in sorted(plex_to_files.keys()):
        if plex_id not in available_fastas:
            n_skip_fasta += 1
            continue

        found_raws = [
            os.path.join(args.raw_dir, f)
            for f in sorted(plex_to_files[plex_id])
            if f in available_raws
        ]
        if not found_raws:
            n_skip_raw += 1
            continue

        fasta_path = os.path.abspath(
            os.path.join(args.fasta_dir, f"{plex_id}.fasta"))
        out_dir = os.path.join(args.output, "results", plex_id)
        os.makedirs(out_dir, exist_ok=True)

        # Symlink RAW files into results/{plex_id}/.
        # MaxQuant places combined/ in the parent directory of the RAW files,
        # so symlinking here causes combined/ to appear at
        # MQ_SEARCH/results/{plex_id}/combined/ as desired.
        linked_raws = []
        for raw in found_raws:
            raw_name = os.path.basename(raw)
            link = os.path.join(out_dir, raw_name)
            if not os.path.islink(link) and not os.path.exists(link):
                os.symlink(raw, link)
            linked_raws.append(link)

        root = build_mqpar(
            raw_files=linked_raws,
            fasta_path=fasta_path,
            out_dir=os.path.abspath(out_dir),
            plex_id=plex_id,
            channel_entries=plex_to_channels[plex_id],
            n_threads=args.threads,
        )

        mqpar_path = os.path.join(mqpar_dir, f"{plex_id}.xml")
        with open(mqpar_path, "w") as f:
            f.write(prettify(root))

        plex_list.append(plex_id)
        print(f"  {plex_id}: {len(found_raws)} RAW files → {mqpar_path}")

    with open(plex_list_path, "w") as f:
        for p in plex_list:
            f.write(p + "\n")

    print(f"\nPlexes ready:             {len(plex_list)}")
    print(f"Skipped (no FASTA):       {n_skip_fasta}")
    print(f"Skipped (no RAW files):   {n_skip_raw}")
    print(f"mqpar files:              {mqpar_dir}")
    print(f"Plex list:                {plex_list_path}")
    print(f"\nNext: bash scripts/ms_search/maxquant/run_maxquant.sh")


if __name__ == "__main__":
    main()
