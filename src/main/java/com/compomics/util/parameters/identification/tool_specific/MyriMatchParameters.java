package com.compomics.util.parameters.identification.tool_specific;

import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;
import com.compomics.util.parameters.identification.search.DigestionParameters;

/**
 * The MyriMatch specific parameters.
 *
 * @author Harald Barsnes
 */
public class MyriMatchParameters extends ExperimentObject implements IdentificationAlgorithmParameter {

    /**
     * Version number for deserialization.
     */
    static final long serialVersionUID = 8755937399680481097L;
    /**
     * The minimum peptide length.
     */
    private Integer minPeptideLength = 8; // note that for MyriMath default is 6
    /**
     * The maximal peptide length.
     */
    private Integer maxPeptideLength = 30; // note that for MyriMath default is 40
    /**
     * The maximum number of spectrum matches.
     */
    private Integer numberOfSpectrumMatches = 10;
    /**
     * The TicCutoffPercentage.
     */
    private Double ticCutoffPercentage = 0.98;
    /**
     * The maximum number of variable modifications.
     */
    private Integer maxDynamicMods = 2;
    /**
     * By default, when generating peptides from the protein database, a peptide
     * must start and end at a valid cleavage site. Setting this parameter to 0
     * or 1 will reduce that requirement, so that neither terminus or only one
     * terminus of the peptide must match one of the cleavage rules specified in
     * the CleavageRules parameter. This parameter is useful to turn a tryptic
     * digest into a semi-tryptic digest.
     */
    private Integer minTerminiCleavages = 2;
    /**
     * The minimum precursor mass considered.
     */
    private Double minPrecursorMass = 600.0;
    /**
     * The maximum precursor mass considered.
     */
    private Double maxPrecursorMass = 5000.0;
    /**
     * If true, the UseSmartPlusThreeModel is used.
     */
    private Boolean useSmartPlusThreeModel = false;
    /**
     * If true, a Sequest-like cross correlation (xcorr) score will be
     * calculated for the top ranking hits in each spectrum’s result set.
     */
    private Boolean computeXCorr = false;
    /**
     * The number of intensity classes.
     */
    private Integer numIntensityClasses = 3;
    /**
     * The multiplier controlling the size of each intensity class relative to
     * the class above it.
     */
    private Integer classSizeMultiplier = 2;
    /**
     * The number of batches per node to strive for when using the MPI-based
     * parallelization features.
     */
    private Integer numberOfBatches = 50;
    /**
     * The lower isotope correction range.
     *
     * @deprecated now general search setting
     */
    private final Integer lowerIsotopeCorrection = -1; // note: the new default is 0
    /**
     * The upper isotope correction range.
     *
     * @deprecated now general search setting
     */
    private final Integer upperIsotopeCorrection = 2; // note: new default is 1
    /**
     * The fragmentation rules. CID (b, y), ETD (c, z*) or manual (user-defined
     * (a comma-separated list of [abcxyz] or z* (z+1), e.g. manual:b,y,z)
     */
    private String fragmentationRule = "CID";
    /**
     * The max number of peaks to use.
     */
    private Integer maxPeakCount = 300;
    /**
     * The output format.
     */
    private String outputFormat = "mzIdentML";

    /**
     * @TODO: parameters not currently supported:
     *
     * AvgPrecursorMzTolerance: "1.5mz", EstimateSearchTimeOnly: "0",
     * KeepUnadjustedPrecursorMz: "0", MaxFragmentChargeState: "0",
     * MaxPeptideVariants: "1000000", MinMatchedFragments: "5", MinResultScore:
     * "9.9999999999999995e-008", NumMzFidelityClasses: "3",
     * PreferIntenseComplements: "1", ProteinSamplingTime: "15",
     * ResultsPerBatch: "200000"
     */
    /**
     * Constructor.
     */
    public MyriMatchParameters() {
    }

    @Override
    public Advocate getAlgorithm() {
        return Advocate.myriMatch;
    }

    @Override
    public boolean equals(IdentificationAlgorithmParameter identificationAlgorithmParameter) {

        if (identificationAlgorithmParameter instanceof MyriMatchParameters) {
            MyriMatchParameters myriMatchParameters = (MyriMatchParameters) identificationAlgorithmParameter;
            if (!minPeptideLength.equals(myriMatchParameters.getMinPeptideLength())) {
                return false;
            }
            if (!maxPeptideLength.equals(myriMatchParameters.getMaxPeptideLength())) {
                return false;
            }
            if (!numberOfSpectrumMatches.equals(myriMatchParameters.getNumberOfSpectrumMatches())) {
                return false;
            }
            double diff = Math.abs(ticCutoffPercentage - myriMatchParameters.getTicCutoffPercentage());
            if (diff > 0.0000000000001) {
                return false;
            }
            if (!maxDynamicMods.equals(myriMatchParameters.getMaxDynamicMods())) {
                return false;
            }
            if (!minTerminiCleavages.equals(myriMatchParameters.getMinTerminiCleavages())) {
                return false;
            }
            diff = Math.abs(minPrecursorMass - myriMatchParameters.getMinPrecursorMass());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(maxPrecursorMass - myriMatchParameters.getMaxPrecursorMass());
            if (diff > 0.0000000000001) {
                return false;
            }
            if (useSmartPlusThreeModel != myriMatchParameters.getUseSmartPlusThreeModel()) {
                return false;
            }
            if (computeXCorr != myriMatchParameters.getComputeXCorr()) {
                return false;
            }
            if (!numIntensityClasses.equals(myriMatchParameters.getNumIntensityClasses())) {
                return false;
            }
            if (!classSizeMultiplier.equals(myriMatchParameters.getClassSizeMultiplier())) {
                return false;
            }
            if (!numberOfBatches.equals(myriMatchParameters.getNumberOfBatches())) {
                return false;
            }
            if (!fragmentationRule.equalsIgnoreCase(myriMatchParameters.getFragmentationRule())) {
                return false;
            }
            if (!maxPeakCount.equals(myriMatchParameters.getMaxPeakCount())) {
                return false;
            }
            if (!getOutputFormat().equalsIgnoreCase(myriMatchParameters.getOutputFormat())) {
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
        output.append("NUMBER_SPECTRUM_MATCHES=");
        output.append(numberOfSpectrumMatches);
        output.append(newLine);
        output.append("TIC_CUTOFF_PERCENTAGE=");
        output.append(ticCutoffPercentage);
        output.append(newLine);
        output.append("MAX_DYNMIC_MODS=");
        output.append(maxDynamicMods);
        output.append(newLine);
        output.append("MIN_TERMINI_CLEAVAGES=");
        output.append(minTerminiCleavages);
        output.append(newLine);
        output.append("MIN_PRECURSOR_MASS=");
        output.append(minPrecursorMass);
        output.append(newLine);
        output.append("MAX_PRECURSOR_MASS=");
        output.append(maxPrecursorMass);
        output.append(newLine);
        output.append("USE_SMART_PLUS_THREE_MODEL=");
        output.append(useSmartPlusThreeModel);
        output.append(newLine);
        output.append("COMPUTE_XCORR=");
        output.append(computeXCorr);
        output.append(newLine);
        output.append("NUM_INTENSITY_CLASSES=");
        output.append(numIntensityClasses);
        output.append(newLine);
        output.append("CLASS_SIZE_MULTIPLIER=");
        output.append(classSizeMultiplier);
        output.append(newLine);
        output.append("NUM_BATCHES=");
        output.append(numberOfBatches);
        output.append(newLine);
        output.append("FRAGMENTATION_RULE=");
        output.append(fragmentationRule);
        output.append(newLine);
        output.append("MAX_PEAK_COUNT=");
        output.append(maxPeakCount);
        output.append(newLine);
        output.append("OUTPUT_FORMAT=");
        output.append(outputFormat);
        output.append(newLine);

        return output.toString();
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
     * Returns the maximum number of spectrum matches.
     *
     * @return the numberOfSpectrumMarches
     */
    public Integer getNumberOfSpectrumMatches() {
        if (numberOfSpectrumMatches == null) {
            numberOfSpectrumMatches = 10;
        }
        return numberOfSpectrumMatches;
    }

    /**
     * Set the maximum number of spectrum matches.
     *
     * @param numberOfSpectrumMarches the numberOfSpectrumMarches to set
     */
    public void setNumberOfSpectrumMatches(Integer numberOfSpectrumMarches) {
        this.numberOfSpectrumMatches = numberOfSpectrumMarches;
    }

    /**
     * Returns the TicCutoffPercentage.
     *
     * @return the ticCutoffPercentage
     */
    public Double getTicCutoffPercentage() {
        return ticCutoffPercentage;
    }

    /**
     * Set the TicCutoffPercentage.
     *
     * @param ticCutoffPercentage the ticCutoffPercentage to set
     */
    public void setTicCutoffPercentage(Double ticCutoffPercentage) {
        this.ticCutoffPercentage = ticCutoffPercentage;
    }

    /**
     * Returns the maximum number of variable modifications.
     *
     * @return the maxDynamicMods
     */
    public Integer getMaxDynamicMods() {
        return maxDynamicMods;
    }

    /**
     * Set the maximum number of variable modifications.
     *
     * @param maxDynamicMods the maxDynamicMods to set
     */
    public void setMaxDynamicMods(Integer maxDynamicMods) {
        this.maxDynamicMods = maxDynamicMods;
    }

    /**
     * Returns the minimum number of termini cleavages.
     *
     * @return the maxDynamicMods
     */
    public Integer getMinTerminiCleavages() {
        return minTerminiCleavages;
    }

    /**
     * Set the minimum number of termini cleavages.
     *
     * @param minTerminiCleavages the minTerminiCleavages to set
     */
    public void setMinTerminiCleavages(Integer minTerminiCleavages) {
        this.minTerminiCleavages = minTerminiCleavages;
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
     * Returns true if the UseSmartPlusThreeModel is to be used.
     *
     * @return true if the UseSmartPlusThreeModel is to be used
     */
    public boolean getUseSmartPlusThreeModel() {
        return useSmartPlusThreeModel;
    }

    /**
     * Sets if the UseSmartPlusThreeModel is to be used.
     *
     * @param useSmartPlusThreeModel if the UseSmartPlusThreeModel is to be used
     */
    public void setUseSmartPlusThreeModel(boolean useSmartPlusThreeModel) {
        this.useSmartPlusThreeModel = useSmartPlusThreeModel;
    }

    /**
     * Returns true if a Sequest-like cross correlation score will be calculated
     * for the top ranking hits in each spectrum’s result set.
     *
     * @return true if the Sequest-like cross correlation score is to be
     * calculated
     */
    public boolean getComputeXCorr() {
        return computeXCorr;
    }

    /**
     * Sets if a Sequest-like cross correlation score will be calculated for the
     * top ranking hits in each spectrum’s result set.
     *
     * @param computeXCorr if the Sequest-like cross correlation score is to be
     * calculated
     */
    public void setComputeXCorr(boolean computeXCorr) {
        this.computeXCorr = computeXCorr;
    }

    /**
     * Returns the number of intensity classes.
     *
     * @return the number of intensity classes
     */
    public Integer getNumIntensityClasses() {
        return numIntensityClasses;
    }

    /**
     * Set the number of intensity classes.
     *
     * @param numIntensityClasses he number of intensity classes
     */
    public void setNumIntensityClasses(Integer numIntensityClasses) {
        this.numIntensityClasses = numIntensityClasses;
    }

    /**
     * Returns the intensity class size multiplier.
     *
     * @return the intensity class size multiplier
     */
    public Integer getClassSizeMultiplier() {
        return classSizeMultiplier;
    }

    /**
     * Set the intensity class size multiplier.
     *
     * @param classSizeMultiplier the intensity class size multiplier
     */
    public void setClassSizeMultiplier(Integer classSizeMultiplier) {
        this.classSizeMultiplier = classSizeMultiplier;
    }

    /**
     * Set the number of batches per node to strive for when using the MPI-based
     * parallelization features.
     *
     * @return the number of batches per node
     */
    public Integer getNumberOfBatches() {
        return numberOfBatches;
    }

    /**
     * Set the number of batches per node to strive for when using the MPI-based
     * parallelization features.
     *
     * @param numberOfBatches the number of batches per node
     */
    public void setNumberOfBatches(Integer numberOfBatches) {
        this.numberOfBatches = numberOfBatches;
    }

    /**
     * Returns the fragmentation rule.
     *
     * @return the fragmentation rule
     */
    public String getFragmentationRule() {
        return fragmentationRule;
    }

    /**
     * Set the fragmentation rule.
     *
     * @param fragmentationRule the fragmentation rule
     */
    public void setFragmentationRule(String fragmentationRule) {
        this.fragmentationRule = fragmentationRule;
    }

    /**
     * Returns the max peak count.
     *
     * @return the max peak count
     */
    public Integer getMaxPeakCount() {
        if (maxPeakCount == null) {
            maxPeakCount = 100;
        }
        return maxPeakCount;
    }

    /**
     * Set the max peak count.
     *
     * @param maxPeakCount the max peak count
     */
    public void setMaxPeakCount(Integer maxPeakCount) {
        this.maxPeakCount = maxPeakCount;
    }

    /**
     * Returns the output format.
     *
     * @return the outputFormat
     */
    public String getOutputFormat() {
        if (outputFormat == null) {
            outputFormat = "mzIdentML";
        }
        return outputFormat;
    }

    /**
     * Set the output format.
     *
     * @param outputFormat the outputFormat to set
     */
    public void setOutputFormat(String outputFormat) {
        this.outputFormat = outputFormat;
    }

    /**
     * Tries to map the utilities enzyme in the digestion preferences to the
     * enzymes supported by MyriMatch.
     *
     * @param digestionPreferences the digestion preferences
     *
     * @return the MyriMatch enzyme as a string, or null of no mapping is found
     */
    public static String enzymeMapping(DigestionParameters digestionPreferences) {

        // Try to map to one of the default Myrimatch enzymes
        if (digestionPreferences.getCleavageParameter() == DigestionParameters.CleavageParameter.unSpecific) {
            return "unspecific cleavage";
        }
        if (digestionPreferences.getCleavageParameter() == DigestionParameters.CleavageParameter.wholeProtein) {
            return "no cleavage";
        }
        if (digestionPreferences.getEnzymes().size() == 1) {
            Enzyme enzyme = digestionPreferences.getEnzymes().get(0);
            String enzymeName = enzyme.getName();
            if (enzymeName.equalsIgnoreCase("Trypsin")) {
                return "Trypsin";
            } else if (enzymeName.equalsIgnoreCase("Trypsin (no P rule)")) {
                return "Trypsin/P";
            } else if (enzymeName.equalsIgnoreCase("Chymotrypsin")) {
                return "Chymotrypsin";
            } else if (enzymeName.equalsIgnoreCase("Glu-C")) {
                return "glutamyl endopeptidase";
            } else if (enzymeName.equalsIgnoreCase("Arg-C")) {
                return "Arg-C";
            } else if (enzymeName.equalsIgnoreCase("Asp-N")) {
                return "Asp-N";
            } else if (enzymeName.equalsIgnoreCase("CNBr")) {
                return "CNBr";
            } else if (enzymeName.equalsIgnoreCase("Formic Acid")) {
                return "Formic_acid";
            } else if (enzymeName.equalsIgnoreCase("Lys-C")) {
                return "Lys-C/P";
            } else if (enzymeName.equalsIgnoreCase("Lys-C (no P rule)")) {
                return "Lys-C";
            } else if (enzymeName.equalsIgnoreCase("Pepsin A")) {
                return "Pepsin A";
            } else if (enzymeName.equalsIgnoreCase("Asp-N")) {
                return "Asp-N";
            }
        }

        // Make a custom cleavage
        return digestionPreferences.getMyriMatchFormat();
    }
}
