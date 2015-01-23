package com.compomics.util.experiment.identification.identification_parameters;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.IdentificationAlgorithmParameter;

/**
 * The Tide specific parameters.
 *
 * @author Harald Barsnes
 */
public class TideParameters implements IdentificationAlgorithmParameter {

    /**
     * Version number for deserialization.
     */
    static final long serialVersionUID = 1049197890002802776L;
    /**
     * The name of the FASTA index folder.
     */
    private String fastIndexFolderName = "fasta-index";
    /**
     * The maximum number of variable modifications allowed on a single peptide.
     * The default is no limit (set to null).
     */
    private Integer maxVariablePtmsPerPeptide = null;
    /**
     * The maximum number of variable modifications of each type allowed on a
     * single peptide.
     */
    private Integer maxVariablePtmsPerTypePerPeptide = 2; // @TODO: make this modification specific?
    /**
     * The minimum peptide length.
     */
    private Integer minPeptideLength = 6;
    /**
     * The maximal peptide length.
     */
    private Integer maxPeptideLength = 30; // note that for Tide default is 50
    /**
     * The minimum precursor mass considered.
     */
    private Double minPrecursorMass = 200.0;
    /**
     * The maximum precursor mass considered.
     */
    private Double maxPrecursorMass = 7200.0;
    /**
     * The decoy format.
     */
    private String decoyFormat = "none"; // none|shuffle|peptide-reverse|protein-reverse
    /**
     * Keep terminal amino acids when creating decoys.
     */
    private String keepTerminalAminoAcids = "NC"; // N|C|NC|none
    /**
     * The seeed of the random number generator with the given unsigned integer
     * when generating the decoy decoys. When given the string "time," the seed
     * is set with the system time.
     */
    private Integer decoySeed = 1;
    /**
     * The name of the output folder (relative to the Tide working folder).
     */
    private String outputFolderName = "crux-output";
    /**
     * If true, a list of all the peptides will be printed to the output folder.
     */
    private Boolean printPeptides = false;
    /**
     * The verbosity of the progress output: 0-fatal errors, 10-non-fatal
     * errors, 20-warnings, 30-information on the progress of execution, 40-more
     * progress information, 50-debug info, 60-detailed debug info.
     */
    private Integer verbosity = 30;
    /**
     * If true, a monoisotopic precursor mass is used, false uses average.
     */
    private Boolean monoisotopicPrecursor = true;
    /**
     * If true, starting methionine peptides will be included both with and
     * without the starting M.
     */
    private Boolean clipNtermMethionine = false;
    /**
     * The enzyme digestion type: full-digest or partial-digest.
     */
    private String digestionType = "full-digest";
    /**
     * If true, the SP score will be computed.
     */
    private Boolean computeSpScore = false;
    /**
     * The maximum number of spectrum matches per peptide.
     */
    private Integer numberOfSpectrumMatches = 10;
    /**
     * If true, the exact p-values will be computed.
     */
    private Boolean computeExactPValues = false;
    /**
     * The minimum spectrum m/z to search for.
     */
    private Double minSpectrumMz = 0.0;
    /**
     * The maximum spectrum m/z to search for. Null if not maximum.
     */
    private Double maxSpectrumMz = null;
    /**
     * The minimum number of peaks in a spectrum for it to be searched.
     */
    private Integer minSpectrumPeaks = 20;
    /**
     * The spectrum charges to search for: 1,2,3 or all. With 'all' every
     * spectrum will be searched and spectra with multiple charge states will be
     * searched once at each charge state. With 1, 2, or 3 only spectra with
     * that charge will be searched.
     */
    private String spectrumCharges = "all";
    /**
     * If true, the precursor peak will be removed. The range removed is
     * specified by
     */
    private Boolean removePrecursor = false;
    /**
     * The tolerance in (Th) used when removing the precursor using
     * removePrecursor.
     */
    private Double removePrecursorTolerance = 1.5;
    /**
     * Show search progress by printing every n spectra searched. Set to 0 to
     * show no search progress.
     */
    private Integer printProgressIndicatorSize = 1000;
    /**
     * If true, the search includes flanking peaks around singly charged b and y
     * theoretical ions. Each flanking peak occurs in the adjacent m/z bin and
     * has half the intensity of the primary peak.
     */
    private Boolean useFlankingPeaks = false;
    /**
     * Controls whether neutral loss ions are considered in the search. Two
     * types of neutral losses are included and are applied only to singly
     * charged b- and y-ions: loss of ammonia (NH3, 17.0086343 Da) and H2O
     * (18.0091422). Each neutral loss peak has intensity 1/5 of the primary
     * peak.
     */
    private Boolean useNeutralLossPeaks = false;
    /**
     * Before calculation of the XCorr score, the m/z axes of the observed and
     * theoretical spectra are discretized. This parameter specifies the size of
     * each bin. The exact formula is floor((x/mz-bin-width) + 1.0 -
     * mz-bin-offset), where x is the observed m/z value. For low resolution ion
     * trap ms/ms data 1.0005079 and for high resolution ms/ms 0.02 is
     * recommended.
     */
    private Double mzBinWidth = 0.02;
    /**
     * In the discretization of the m/z axes of the observed and theoretical
     * spectra, this parameter specifies the location of the left edge of the
     * first bin, relative to mass = 0 (i.e., mz-bin-offset = 0.xx means the
     * left edge of the first bin will be located at +0.xx Da). The parameter
     * must lie in the range 0 ≤ mz-bin-offset ≤ 1.
     */
    private Double mzBinOffset = 0.0;
    /**
     * If true, target and decoy search results are reported in a single file
     * named "tide-search.txt," and only the top-scoring N matches (as specified
     * via --top-match) are reported for each spectrum, irrespective of whether
     * the matches involve target or decoy peptides.
     */
    private Boolean concatenateTargetDecoy = false;
    /**
     * Specify the name of the file where the binarized fragmentation spectra
     * will be stored. Subsequent runs of crux tide-search will execute more
     * quickly if provided with the spectra in binary format. The filename is
     * specified relative to the current working directory, not the Crux output
     * directory (as specified by --output-dir).
     */
    private String storeSpectraFileName = null;
    /**
     * If true, tab delimited text file output is generated.
     */
    private Boolean textOutput = true;
    /**
     * If true, SQT output is generated.
     */
    private Boolean sqtOutput = false;
    /**
     * If true, pepxml output is generated.
     */
    private Boolean pepXmlOutput = false;
    /**
     * If true, mzid output is generated.
     */
    private Boolean mzidOutput = false;
    /**
     * If true, Percolator input file is generated.
     */
    private Boolean pinOutput = false;
    /**
     * If true, the tide output and index folders are removed when the search
     * has completed.
     */
    private Boolean removeTempFolders = true;

    /**
     * Constructor.
     */
    public TideParameters() {
        // @TODO: add --peptide-centric-search?
    }

    @Override
    public Advocate getAlgorithm() {
        return Advocate.tide;
    }

    @Override
    public boolean equals(IdentificationAlgorithmParameter identificationAlgorithmParameter) {

        if (identificationAlgorithmParameter instanceof TideParameters) {
            TideParameters tideParameters = (TideParameters) identificationAlgorithmParameter;

            if (minPeptideLength != tideParameters.getMinPeptideLength()) {
                return false;
            }
            if (maxPeptideLength != tideParameters.getMaxPeptideLength()) {
                return false;
            }
            double diff = Math.abs(minPrecursorMass - tideParameters.getMinPrecursorMass());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(maxPrecursorMass - tideParameters.getMaxPrecursorMass());
            if (diff > 0.0000000000001) {
                return false;
            }
            if (maxVariablePtmsPerPeptide != null && tideParameters.getMaxVariablePtmsPerPeptide() == null) {
                return false;
            }
            if (maxVariablePtmsPerPeptide == null && tideParameters.getMaxVariablePtmsPerPeptide() != null) {
                return false;
            }
            if (maxVariablePtmsPerPeptide != null && tideParameters.getMaxVariablePtmsPerPeptide() != null) {
                if (maxVariablePtmsPerPeptide.intValue() != tideParameters.getMaxVariablePtmsPerPeptide().intValue()) {
                    return false;
                }
            }
            if (maxVariablePtmsPerTypePerPeptide != null && tideParameters.getMaxVariablePtmsPerTypePerPeptide() == null) {
                return false;
            }
            if (maxVariablePtmsPerTypePerPeptide == null && tideParameters.getMaxVariablePtmsPerTypePerPeptide() != null) {
                return false;
            }
            if (maxVariablePtmsPerTypePerPeptide != null && tideParameters.getMaxVariablePtmsPerTypePerPeptide() != null) {
                if (maxVariablePtmsPerTypePerPeptide.intValue() != tideParameters.getMaxVariablePtmsPerTypePerPeptide().intValue()) {
                    return false;
                }
            }
            if (!decoyFormat.equalsIgnoreCase(tideParameters.getDecoyFormat())) {
                return false;
            }
            if (!keepTerminalAminoAcids.equalsIgnoreCase(tideParameters.getKeepTerminalAminoAcids())) {
                return false;
            }
            if (decoySeed != tideParameters.getDecoySeed().intValue()) {
                return false;
            }
            if (outputFolderName.equalsIgnoreCase(tideParameters.getOutputFolderName())) {
                return false;
            }
            if (printPeptides != tideParameters.getPrintPeptides().booleanValue()) {
                return false;
            }
            if (verbosity != tideParameters.getVerbosity().intValue()) {
                return false;
            }
            if (monoisotopicPrecursor != tideParameters.getMonoisotopicPrecursor().booleanValue()) {
                return false;
            }
            if (clipNtermMethionine != tideParameters.getClipNtermMethionine().booleanValue()) {
                return false;
            }
            if (digestionType.equalsIgnoreCase(tideParameters.getDigestionType())) {
                return false;
            }
            if (computeSpScore != tideParameters.getComputeSpScore().booleanValue()) {
                return false;
            }
            if (numberOfSpectrumMatches != tideParameters.getNumberOfSpectrumMatches().intValue()) {
                return false;
            }
            if (computeExactPValues != tideParameters.getComputeExactPValues().booleanValue()) {
                return false;
            }
            if (minSpectrumMz != tideParameters.getMinSpectrumMz().doubleValue()) {
                return false;
            }
            if (maxSpectrumMz != tideParameters.getMaxSpectrumMz().doubleValue()) {
                return false;
            }
            if (minSpectrumPeaks != tideParameters.getMinSpectrumPeaks().intValue()) {
                return false;
            }
            if (spectrumCharges.equalsIgnoreCase(tideParameters.getSpectrumCharges())) {
                return false;
            }
            if (removePrecursor != tideParameters.getRemovePrecursor().booleanValue()) {
                return false;
            }
            if (removePrecursorTolerance != tideParameters.getRemovePrecursorTolerance().doubleValue()) {
                return false;
            }
            if (printProgressIndicatorSize != tideParameters.getPrintProgressIndicatorSize().intValue()) {
                return false;
            }
            if (useFlankingPeaks != tideParameters.getUseFlankingPeaks().booleanValue()) {
                return false;
            }
            if (useNeutralLossPeaks != tideParameters.getUseNeutralLossPeaks().booleanValue()) {
                return false;
            }
            if (mzBinWidth != tideParameters.getMzBinWidth().doubleValue()) {
                return false;
            }
            if (mzBinOffset != tideParameters.getMzBinOffset().doubleValue()) {
                return false;
            }
            if (concatenateTargetDecoy != tideParameters.getConcatenatTargetDecoy().booleanValue()) {
                return false;
            }
            if (storeSpectraFileName.equalsIgnoreCase(tideParameters.getStoreSpectraFileName())) {
                return false;
            }
            if (textOutput != tideParameters.getTextOutput().booleanValue()) {
                return false;
            }
            if (sqtOutput != tideParameters.getSqtOutput().booleanValue()) {
                return false;
            }
            if (pepXmlOutput != tideParameters.getPepXmlOutput().booleanValue()) {
                return false;
            }
            if (mzidOutput != tideParameters.getMzidOutput().booleanValue()) {
                return false;
            }
            if (pinOutput != tideParameters.getPinOutput().booleanValue()) {
                return false;
            }

            return true;
        }

        return false;
    }

    @Override
    public String toString(boolean html) {
        String newLine = System.getProperty("line.separator");

        if (html) {
            newLine = "<br>";
        }

        StringBuilder output = new StringBuilder();
        Advocate advocate = getAlgorithm();
        output.append("# ------------------------------------------------------------------");
        output.append(newLine);
        output.append("# ").append(advocate.getName()).append(" Specific Parameters");
        output.append(newLine);
        output.append("# ------------------------------------------------------------------");
        output.append(newLine);
        output.append(newLine);

        output.append("MIN_PEP_LENGTH=");
        output.append(minPeptideLength);
        output.append(newLine);
        output.append("MAX_PEP_LENGTH=");
        output.append(maxPeptideLength);
        output.append(newLine);
        output.append("MIN_PRECURSOR_MASS=");
        output.append(minPrecursorMass);
        output.append(newLine);
        output.append("MAX_PRECURSOR_MASS=");
        output.append(maxPrecursorMass);
        output.append(newLine);
        output.append("MAX_VARIABLE_PTMS_PER_TYPE_PER_PEPTIDE=");
        output.append(maxVariablePtmsPerTypePerPeptide);
        output.append(newLine);
        output.append("MAX_VARIABLE_PTMS_PER_PEPTID=");
        output.append(maxVariablePtmsPerPeptide);
        output.append(newLine);
        output.append("DECOY_FORMAT=");
        output.append(decoyFormat);
        output.append(newLine);
        output.append("KEEP_TERMINAL_AMINO_ACIDS=");
        output.append(keepTerminalAminoAcids);
        output.append(newLine);
        output.append("DECOY_SEED=");
        output.append(decoySeed);
        output.append(newLine);
        output.append("OUTPUT_FOLDER_NAME=");
        output.append(outputFolderName);
        output.append(newLine);
        output.append("PRINT_PEPTIDES=");
        output.append(printPeptides);
        output.append(newLine);
        output.append("VERBOSITY=");
        output.append(verbosity);
        output.append(newLine);
        output.append("MONOISOTOPIC_PRECURSOR=");
        output.append(monoisotopicPrecursor);
        output.append(newLine);
        output.append("CLIP_NTERM_METHIONINE=");
        output.append(clipNtermMethionine);
        output.append(newLine);
        output.append("DIGESTION_TYPE=");
        output.append(digestionType);
        output.append(newLine);
        output.append("COMPUTE_SP_SCORE=");
        output.append(computeSpScore);
        output.append(newLine);
        output.append("NUMBER_SPECTRUM_MATCHES=");
        output.append(numberOfSpectrumMatches);
        output.append(newLine);
        output.append("COMPUTE_EXACT_P_VALUES=");
        output.append(computeExactPValues);
        output.append(newLine);
        output.append("MIN_SPECTRUM_MZ=");
        output.append(minSpectrumMz);
        output.append(newLine);
        output.append("MAX_SPECTRUM_MZ=");
        output.append(maxSpectrumMz);
        output.append(newLine);
        output.append("MIN_SPECTRUM_PEAKS=");
        output.append(minSpectrumPeaks);
        output.append(newLine);
        output.append("SPECTRUM_CHARGES=");
        output.append(spectrumCharges);
        output.append(newLine);
        output.append("REMOVE_PRECURSOR=");
        output.append(removePrecursor);
        output.append(newLine);
        output.append("REMOVE_PRECURSOR_TOLERANCE=");
        output.append(removePrecursorTolerance);
        output.append(newLine);
        output.append("PRINT_PROGRESS_INDICATOR_SIZE=");
        output.append(printProgressIndicatorSize);
        output.append(newLine);
        output.append("USE_FLANKING_PEAKS=");
        output.append(useFlankingPeaks);
        output.append(newLine);
        output.append("USE_NEUTRAL_LOSS_PEAKS=");
        output.append(useNeutralLossPeaks);
        output.append(newLine);
        output.append("MZ_BIN_WIDTH=");
        output.append(mzBinWidth);
        output.append(newLine);
        output.append("MZ_BIN_OFFSET=");
        output.append(mzBinOffset);
        output.append(newLine);
        output.append("CONCATENATE_TARGET_DECOY=");
        output.append(concatenateTargetDecoy);
        output.append(newLine);
        output.append("STORE_SPECTRA_FILE_NAME=");
        output.append(storeSpectraFileName);
        output.append(newLine);
        output.append("TEXT_OUTPUT=");
        output.append(textOutput);
        output.append(newLine);
        output.append("SQT_OUTPUT=");
        output.append(sqtOutput);
        output.append(newLine);
        output.append("PEPXML_OUTPUT=");
        output.append(pepXmlOutput);
        output.append(newLine);
        output.append("MZID_OUTPUT=");
        output.append(mzidOutput);
        output.append(newLine);
        output.append("PERCOLATOR_OUTPUT=");
        output.append(pinOutput);
        output.append(newLine);

        return output.toString();
    }

    /**
     * Returns the maximum number of variable modifications allowed on a single
     * peptide. Null if no upper limit is set.
     *
     * @return the maxVariablePtmsPerPeptide
     */
    public Integer getMaxVariablePtmsPerPeptide() {
        return maxVariablePtmsPerPeptide;
    }

    /**
     * Set the maximum number of variable modifications of allowed on a single
     * peptide. Set to null if no upper limit is to be used.
     *
     * @param maxVariablePtmsPerPeptide the maxVariablePtmsPerPeptide to set
     */
    public void setMaxVariablePtmsPerPeptide(Integer maxVariablePtmsPerPeptide) {
        this.maxVariablePtmsPerPeptide = maxVariablePtmsPerPeptide;
    }

    /**
     * Returns the maximum number of variable modifications of each modification
     * type allowed on a single peptide. Null if no upper limit is set.
     *
     * @return the maxVariablePtmsPerPeptide
     */
    public Integer getMaxVariablePtmsPerTypePerPeptide() {
        return maxVariablePtmsPerTypePerPeptide;
    }

    /**
     * Set the maximum number of variable modifications of each modification
     * type allowed on a single peptide. Set to null if no upper limit is to be
     * used.
     *
     * @param maxVariablePtmsPerTypePerPeptide the
     * maxVariablePtmsPerTypePerPeptide to set
     */
    public void setMaxVariablePtmsPerTypePerPeptide(Integer maxVariablePtmsPerTypePerPeptide) {
        this.maxVariablePtmsPerTypePerPeptide = maxVariablePtmsPerTypePerPeptide;
    }

    /**
     * Returns the maximal peptide length allowed.
     *
     * @return the maximal peptide length allowed
     */
    public Integer getMaxPeptideLength() {
        return maxPeptideLength;
    }

    /**
     * Sets the maximal peptide length allowed.
     *
     * @param maxPeptideLength the maximal peptide length allowed
     */
    public void setMaxPeptideLength(Integer maxPeptideLength) {
        this.maxPeptideLength = maxPeptideLength;
    }

    /**
     * Sets the minimal peptide length allowed.
     *
     * @return the minimal peptide length allowed
     */
    public Integer getMinPeptideLength() {
        return minPeptideLength;
    }

    /**
     * Sets the minimal peptide length allowed.
     *
     * @param minPeptideLength the minimal peptide length allowed
     */
    public void setMinPeptideLength(Integer minPeptideLength) {
        this.minPeptideLength = minPeptideLength;
    }

    /**
     * Returns the minimum precursor mass.
     *
     * @return the minimum precursor mass
     */
    public Double getMinPrecursorMass() {
        return minPrecursorMass;
    }

    /**
     * Sets the minimum precursor mass.
     *
     * @param minPrecursorMass the minPrecursorMass to set
     */
    public void setMinPrecursorMass(Double minPrecursorMass) {
        this.minPrecursorMass = minPrecursorMass;
    }

    /**
     * Returns the maxPrecursorMass precursor mass.
     *
     * @return the maximum precursor mass
     */
    public Double getMaxPrecursorMass() {
        return maxPrecursorMass;
    }

    /**
     * Sets the maximum precursor mass.
     *
     * @param maxPrecursorMass the maximum to set
     */
    public void setMaxPrecursorMass(Double maxPrecursorMass) {
        this.maxPrecursorMass = maxPrecursorMass;
    }

    /**
     * Returns the decoy format.
     *
     * @return the decoyFormat
     */
    public String getDecoyFormat() {
        return decoyFormat;
    }

    /**
     * Set the decoy format.
     *
     * @param decoyFormat the decoyFormat to set
     */
    public void setDecoyFormat(String decoyFormat) {
        this.decoyFormat = decoyFormat;
    }

    /**
     * Returns the option for keeping the terminal amino acids when generating
     * the decoys.
     *
     * @return the keepTerminalAminoAcids
     */
    public String getKeepTerminalAminoAcids() {
        return keepTerminalAminoAcids;
    }

    /**
     * Set the option for keeping the terminal amino acids when generating the
     * decoys
     *
     * @param keepTerminalAminoAcids the keepTerminalAminoAcids to set
     */
    public void setKeepTerminalAminoAcids(String keepTerminalAminoAcids) {
        this.keepTerminalAminoAcids = keepTerminalAminoAcids;
    }

    /**
     * Returns the decoy seed.
     *
     * @return the decoySeed
     */
    public Integer getDecoySeed() {
        return decoySeed;
    }

    /**
     * Set the decoy seed.
     *
     * @param decoySeed the decoySeed to set
     */
    public void setDecoySeed(Integer decoySeed) {
        this.decoySeed = decoySeed;
    }

    /**
     * Returns the name of the output folder.
     *
     * @return the outputFolderName
     */
    public String getOutputFolderName() {
        return outputFolderName;
    }

    /**
     * Set the name of the output folder.
     *
     * @param outputFolderName the outputFolderName to set
     */
    public void setOutputFolderName(String outputFolderName) {
        this.outputFolderName = outputFolderName;
    }

    /**
     * Returns true of a list of the peptides should be printed to the output
     * folder.
     *
     * @return the printPeptides
     */
    public Boolean getPrintPeptides() {
        return printPeptides;
    }

    /**
     * Set if a list of the peptides should be printed to the output folder.
     *
     * @param printPeptides the printPeptides to set
     */
    public void setPrintPeptides(Boolean printPeptides) {
        this.printPeptides = printPeptides;
    }

    /**
     * Returns the verbosity type of progress output.
     *
     * @return the verbosity
     */
    public Integer getVerbosity() {
        return verbosity;
    }

    /**
     * Set the verbosity type of progress output.
     *
     * @param verbosity the verbosity to set
     */
    public void setVerbosity(Integer verbosity) {
        this.verbosity = verbosity;
    }

    /**
     * Returns true if the precursor mass is monoisotopic, false if average.
     *
     * @return the monoisotopicPrecursor
     */
    public Boolean getMonoisotopicPrecursor() {
        return monoisotopicPrecursor;
    }

    /**
     * Set if the precursor mass is monoisotopic.
     *
     * @param monoisotopicPrecursor the monoisotopicPrecursor to set
     */
    public void setMonoisotopicPrecursor(Boolean monoisotopicPrecursor) {
        this.monoisotopicPrecursor = monoisotopicPrecursor;
    }

    /**
     * Returns if the starting methionine peptides will be included both with
     * and without the starting M.
     *
     * @return the clipNtermMethionine
     */
    public Boolean getClipNtermMethionine() {
        return clipNtermMethionine;
    }

    /**
     * Set if the starting methionine peptides will be included both with and
     * without the starting M.
     *
     * @param clipNtermMethionine the clipNtermMethionine to set
     */
    public void setClipNtermMethionine(Boolean clipNtermMethionine) {
        this.clipNtermMethionine = clipNtermMethionine;
    }

    /**
     * Returns the enzyme digestion type.
     *
     * @return the digestionType
     */
    public String getDigestionType() {
        return digestionType;
    }

    /**
     * Set the enzyme digestion type.
     *
     * @param digestionType the digestionType to set
     */
    public void setDigestionType(String digestionType) {
        this.digestionType = digestionType;
    }

    /**
     * Returns true of the SP score is to be computed.
     *
     * @return the computeSpScore
     */
    public Boolean getComputeSpScore() {
        return computeSpScore;
    }

    /**
     * Set if the SP score is to be computed.
     *
     * @param computeSpScore the computeSpScore to set
     */
    public void setComputeSpScore(Boolean computeSpScore) {
        this.computeSpScore = computeSpScore;
    }

    /**
     * Returns the number of PSMs to report per spectrum.
     *
     * @return the numberOfSpectrumMatches
     */
    public Integer getNumberOfSpectrumMatches() {
        return numberOfSpectrumMatches;
    }

    /**
     * Set the number of PSMs to report per spectrum.
     *
     * @param numberOfSpectrumMatches the numberOfSpectrumMatches to set
     */
    public void setNumberOfSpectrumMatches(Integer numberOfSpectrumMatches) {
        this.numberOfSpectrumMatches = numberOfSpectrumMatches;
    }

    /**
     * Returns true if the exact p-values are to be computed.
     *
     * @return the computeExactPValues
     */
    public Boolean getComputeExactPValues() {
        return computeExactPValues;
    }

    /**
     * Set if the exact p-values are to be computed.
     *
     * @param computeExactPValues the computeExactPValues to set
     */
    public void setComputeExactPValues(Boolean computeExactPValues) {
        this.computeExactPValues = computeExactPValues;
    }

    /**
     * Return the min spectrum m/z to search for.
     *
     * @return the minSpectrumMz
     */
    public Double getMinSpectrumMz() {
        return minSpectrumMz;
    }

    /**
     * Set the min spectrum m/z to search for.
     *
     * @param minSpectrumMz the minSpectrumMz to set
     */
    public void setMinSpectrumMz(Double minSpectrumMz) {
        this.minSpectrumMz = minSpectrumMz;
    }

    /**
     * Returns the max spectrum m/z to search for.
     *
     * @return the maxSpectrumMz
     */
    public Double getMaxSpectrumMz() {
        return maxSpectrumMz;
    }

    /**
     * Set the max spectrum m/z to search for.
     *
     * @param maxSpectrumMz the maxSpectrumMz to set
     */
    public void setMaxSpectrumMz(Double maxSpectrumMz) {
        this.maxSpectrumMz = maxSpectrumMz;
    }

    /**
     * Returns the min number of peaks in a spectrum.
     *
     * @return the minSpectrumPeaks
     */
    public Integer getMinSpectrumPeaks() {
        return minSpectrumPeaks;
    }

    /**
     * Set the min number of peaks in a spectrum.
     *
     * @param minSpectrumPeaks the minSpectrumPeaks to set
     */
    public void setMinSpectrumPeaks(Integer minSpectrumPeaks) {
        this.minSpectrumPeaks = minSpectrumPeaks;
    }

    /**
     * Returns the spectrum charges to search for.
     *
     * @return the spectrumCharges
     */
    public String getSpectrumCharges() {
        return spectrumCharges;
    }

    /**
     * Set the spectrum charges to search for.
     *
     * @param spectrumCharges the spectrumCharges to set
     */
    public void setSpectrumCharges(String spectrumCharges) {
        this.spectrumCharges = spectrumCharges;
    }

    /**
     * Returns true of the precursor peak is to be removed.
     *
     * @return the removePrecursor
     */
    public Boolean getRemovePrecursor() {
        return removePrecursor;
    }

    /**
     * Set if the precursor peak is to be removed.
     *
     * @param removePrecursor the removePrecursor to set
     */
    public void setRemovePrecursor(Boolean removePrecursor) {
        this.removePrecursor = removePrecursor;
    }

    /**
     * Returns the tolerance with which to remove the precursor peak.
     *
     * @return the removePrecursorTolerance
     */
    public Double getRemovePrecursorTolerance() {
        return removePrecursorTolerance;
    }

    /**
     * Set the tolerance with which to remove the precursor peak.
     *
     * @param removePrecursorTolerance the removePrecursorTolerance to set
     */
    public void setRemovePrecursorTolerance(Double removePrecursorTolerance) {
        this.removePrecursorTolerance = removePrecursorTolerance;
    }

    /**
     * Returns the progress indicator size.
     *
     * @return the printProgressIndicatorSize
     */
    public Integer getPrintProgressIndicatorSize() {
        return printProgressIndicatorSize;
    }

    /**
     * Set the progress indicator size.
     *
     * @param printProgressIndicatorSize the printProgressIndicatorSize to set
     */
    public void setPrintProgressIndicatorSize(Integer printProgressIndicatorSize) {
        this.printProgressIndicatorSize = printProgressIndicatorSize;
    }

    /**
     * Returns true if flanking peaks are to be used.
     *
     * @return the useFlankingPeaks
     */
    public Boolean getUseFlankingPeaks() {
        return useFlankingPeaks;
    }

    /**
     * Set if flanking peaks are to be used.
     *
     * @param useFlankingPeaks the useFlankingPeaks to set
     */
    public void setUseFlankingPeaks(Boolean useFlankingPeaks) {
        this.useFlankingPeaks = useFlankingPeaks;
    }

    /**
     * Returns true if the neutral loss peaks are to be used.
     *
     * @return the useNeutralLossPeaks
     */
    public Boolean getUseNeutralLossPeaks() {
        return useNeutralLossPeaks;
    }

    /**
     * Set if the neutral loss peaks are to be used.
     *
     * @param useNeutralLossPeaks the useNeutralLossPeaks to set
     */
    public void setUseNeutralLossPeaks(Boolean useNeutralLossPeaks) {
        this.useNeutralLossPeaks = useNeutralLossPeaks;
    }

    /**
     * Returns the m/z bin width.
     *
     * @return the mzBinWidth
     */
    public Double getMzBinWidth() {
        return mzBinWidth;
    }

    /**
     * Set the m/z bin width.
     *
     * @param mzBinWidth the mzBinWidth to set
     */
    public void setMzBinWidth(Double mzBinWidth) {
        this.mzBinWidth = mzBinWidth;
    }

    /**
     * Returns the m/z bin offset.
     *
     * @return the mzBinOffset
     */
    public Double getMzBinOffset() {
        return mzBinOffset;
    }

    /**
     * Set the m/z bin offset.
     *
     * @param mzBinOffset the mzBinOffset to set
     */
    public void setMzBinOffset(Double mzBinOffset) {
        this.mzBinOffset = mzBinOffset;
    }

    /**
     * Returns true if the target and decoy results are to be concatenated.
     *
     * @return the concatenatTargetDecoy
     */
    public Boolean getConcatenatTargetDecoy() {
        return concatenateTargetDecoy;
    }

    /**
     * Set if the target and decoy results are to be concatenated.
     *
     * @param concatenate the concatenateTargetDecoy to set
     */
    public void setConcatenatTargetDecoy(Boolean concatenate) {
        this.concatenateTargetDecoy = concatenate;
    }

    /**
     * Returns the name of the file where the binary spectra are to be stored.
     *
     * @return the storeSpectraFileName
     */
    public String getStoreSpectraFileName() {
        return storeSpectraFileName;
    }

    /**
     * Set the name of the file where the binary spectra are to be stored.
     *
     * @param storeSpectraFileName the storeSpectraFileName to set
     */
    public void setStoreSpectraFileName(String storeSpectraFileName) {
        this.storeSpectraFileName = storeSpectraFileName;
    }

    /**
     * Returns true if the text output is to be created.
     *
     * @return the textOutput
     */
    public Boolean getTextOutput() {
        return textOutput;
    }

    /**
     * Set if the text output is to be created.
     *
     * @param textOutput the textOutput to set
     */
    public void setTextOutput(Boolean textOutput) {
        this.textOutput = textOutput;
    }

    /**
     * Returns true if the sqt output is to be created.
     *
     * @return the sqtOutput
     */
    public Boolean getSqtOutput() {
        return sqtOutput;
    }

    /**
     * Set if the sqt output is to be created.
     *
     * @param sqtOutput the sqtOutput to set
     */
    public void setSqtOutput(Boolean sqtOutput) {
        this.sqtOutput = sqtOutput;
    }

    /**
     * Returns true if the pepxml output is to be created.
     *
     * @return the pepXmlOutput
     */
    public Boolean getPepXmlOutput() {
        return pepXmlOutput;
    }

    /**
     * Set if the pepxml output is to be created.
     *
     * @param pepXmlOutput the pepXmlOutput to set
     */
    public void setPepXmlOutput(Boolean pepXmlOutput) {
        this.pepXmlOutput = pepXmlOutput;
    }

    /**
     * Returns true if the mzid output is to be created.
     *
     * @return the mzidOutput
     */
    public Boolean getMzidOutput() {
        return mzidOutput;
    }

    /**
     * Set if the mzid output is to be created.
     *
     * @param mzidOutput the mzidOutput to set
     */
    public void setMzidOutput(Boolean mzidOutput) {
        this.mzidOutput = mzidOutput;
    }

    /**
     * Returns true if the Percolator output is to be created.
     *
     * @return the pinOutput
     */
    public Boolean getPinOutput() {
        return pinOutput;
    }

    /**
     * Set if the Percolator output is to be created.
     *
     * @param pinOutput the pinOutput to set
     */
    public void setPinOutput(Boolean pinOutput) {
        this.pinOutput = pinOutput;
    }

    /**
     * Returns the name of the FASTA index folder.
     *
     * @return the fastIndexFolderName
     */
    public String getFastIndexFolderName() {
        return fastIndexFolderName;
    }

    /**
     * Set the name of the FASTA index folder.
     *
     * @param fastIndexFolderName the fastIndexFolderName to set
     */
    public void setFastIndexFolderName(String fastIndexFolderName) {
        this.fastIndexFolderName = fastIndexFolderName;
    }

    /**
     * Returns true if the output and index folders are to be removed when the
     * search has completed.
     *
     * @return the removeTempFolders
     */
    public Boolean getRemoveTempFolders() {
        if (removeTempFolders == null) {
            removeTempFolders = true;
        }
        return removeTempFolders;
    }

    /**
     * Set if the output and index folders are to be removed when the search has
     * completed.
     *
     * @param removeTempFolders the removeTempFolders to set
     */
    public void setRemoveTempFolders(Boolean removeTempFolders) {
        this.removeTempFolders = removeTempFolders;
    }
}
