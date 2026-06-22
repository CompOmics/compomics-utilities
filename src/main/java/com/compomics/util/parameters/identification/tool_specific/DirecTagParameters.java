package com.compomics.util.parameters.identification.tool_specific;

import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;
import java.util.ArrayList;

/**
 * The DirecTag specific parameters.
 *
 * @author Thilo Muth
 * @author Harald Barsnes
 */
public class DirecTagParameters extends ExperimentObject implements IdentificationAlgorithmParameter {

    // @TODO: options not implemented: OutputSuffix, StartSpectraScanNum, EndSpectraScanNum 
    //        (advanced: ClassSizeMultiplier, NumBatches, ThreadCountMultiplier, UseMultipleProcessors)
    //
    /**
     * Serial number for backward compatibility.
     */
    static final long serialVersionUID = -3107089648450731330L;
    /**
     * In order to maximize the effectiveness of the MVH scoring algorithm, an
     * important step in preprocessing the experimental spectra is filtering out
     * noise peaks. Noise peaks are filtered out by sorting the original peaks
     * in descending order of intensity, and then picking peaks from that list
     * until the cumulative ion current of the picked peaks divided by the total
     * ion current (TIC) is greater than or equal to this parameter. Lower
     * percentages mean that less of the spectra total intensity will be
     * allowed to pass through preprocessing. See the section on Advanced Usage
     * for tips on how to use this parameter optimally.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double ticCutoffPercentage = 100;
    /**
     * Another way of increasing the effectiveness of the MVH scoring algorithm
     * when used for tagging is to set an upper bound on the number of peaks in
     * a spectrum before generating tags. This step tends to get rid of most
     * noise peaks and makes tagging much more feasible because so many fewer
     * false positives are generated.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int maxPeakCount = 100;
    /**
     * Before scoring any candidates, experimental spectra have their peaks
     * stratified into the number of intensity classes specified by this
     * parameter. Spectra that are very dense in peaks will likely benefit from
     * more intensity classes in order to best take advantage of the variation
     * in peak intensities. Spectra that are very sparse will not see much
     * benefit from using many intensity classes.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int numIntensityClasses = 3;
    /**
     * If true, the preprocessing step will correct the precursor mass by
     * adjusting it through a specified range in steps of a specified length,
     * finally choosing the optimal adjustment. The optimal adjustment is the
     * one that maximizes the sum of products of all complementary peaks in the
     * spectrum. (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private boolean adjustPrecursorMass = false;
    /**
     * When adjusting the precursor mass, this parameter sets the lower mass
     * limit of adjustment allowable from the original precursor mass, measured
     * in Daltons. (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double minPrecursorAdjustment = -0.5;
    /**
     * When adjusting the precursor mass, this parameter sets the upper mass
     * limit of adjustment allowable from the original precursor mass, measured
     * in Daltons. (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double maxPrecursorAdjustment = 1.5;
    /**
     * When adjusting the precursor mass, this parameter sets the size of the
     * steps between adjustments, measured in Daltons.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double precursorAdjustmentStep = 0.1;
    /**
     * Controls the number of charge states that DirecTag will handle during all
     * stages of the program.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int numChargeStates = 3;
    /**
     * The output of a DirecTag job will be a TAGS file for each input file. The
     * string specified by this parameter will be appended to each TAGS0.5
     * filename. It is useful for differentiating jobs within a single
     * directory. (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private String outputSuffix = "";
    /**
     * If true, DirecTag will use the charge state from the input data if it is
     * available. If false, or if charge state is not available from a
     * particular spectrum, DirecTag will use its internal algorithm to
     * determine charge state.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private boolean useChargeStateFromMS = true;
    /**
     * If DirecTag determines a spectrum to be multiply charged and this
     * parameter is true, the spectrum will be copied and treated as if it was
     * all possible charge states from +2 to +&lt;NumChargeStates&gt;. If this
     * parameter is false, the spectrum will simply be treated as a +2.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private boolean duplicateSpectra = true;
    /**
     * Deisotoping a spectrum (consolidating isotopic peak intensities into the
     * monoisotopic peaks intensity) during preprocessing will significantly
     * improve precursor adjustment, and it may be desirable to keep the
     * deisotoped spectrum around for candidate scoring as well. Set to 0, no
     * deisotoping will be used. Set to 1, deisotoping will be used for
     * precursor adjustment only. Set to 2, deisotoping will be used for both
     * precursor adjustment and for candidate scoring.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int deisotopingMode = 0;
    /**
     * When deisotoping a spectrum, an isotopic peak is one that is the mass of
     * a neutron higher than another peak, tolerating variation based on the
     * value of this parameter. Deisotoping actually traverses the spectrum at
     * multiple charge states, starting from the highest (NumChargeStates) and
     * ending at the lowest.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double isotopeMzTolerance = 0.25;
    /**
     * When adjusting the precursor mass, this parameter controls how much
     * tolerance there is on each side of the calculated m/z when looking for a
     * peaks complement.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double complementMzTolerance = 0.1;
    /**
     * A sequence tag is generated from the gaps between a number of peaks equal
     * to this parameter plus one. Longer tag lengths are more specific, but
     * harder to find because many consecutive ion fragments are rare.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int tagLength = 4;
    /**
     * This parameter sets the maximum number of modified residues that may be
     * in any candidate sequence.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int maxDynamicMods = 2;
    /**
     * This parameter sets the maximum number of sequence tags to report for
     * each spectrum. (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private int maxTagCount = 10;
    /**
     * This parameter controls how intensity scores are combined to form a total
     * score. DirecTag scores tags on the basis of their peak intensities. Tags
     * that contain intense peaks are more likely to be correct than those that
     * contain average peaks.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double intensityScoreWeight = 1.0;
    /**
     * This parameter controls how mzFidelity scores are combined to form a
     * total score. m/z fidelity for a tag can be characterized through SSE.
     * DirecTag evaluates the consistency of fragment ion m/z values for each
     * tag. (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double mzFidelityScoreWeight = 1.0;
    /**
     * This parameter controls how complement scores are combined to form a
     * total score. Peaks that match to complementary ions within the spectrum
     * are more trustworthy than other peaks. DirecTag assesses the number and
     * concordance of complementary ions for each tag.
     * (http://fenchurch.mc.vanderbilt.edu/bumbershoot/directag)
     */
    private double complementScoreWeight = 1.0;
    /**
     * The variable modifications in the order used.
     */
    private ArrayList<String> variableModifications = new ArrayList<>(0);

    /**
     * Constructor.
     */
    public DirecTagParameters() {
    }

    /**
     * Returns the number of charge states.
     *
     * @return numChargeStates the number of charge states
     */
    public int getNumChargeStates() {
        return numChargeStates;
    }

    /**
     * Sets the number of charge states.
     *
     * @param numChargeStates the number of charge states
     */
    public void setNumChargeStates(int numChargeStates) {
        this.numChargeStates = numChargeStates;
    }

    /**
     * Returns the output tags suffix.
     *
     * @return the output tags suffix.
     */
    public String getOutputSuffix() {
        return outputSuffix;
    }

    /**
     * Sets the output tags suffix.
     *
     * @param outputSuffix the output tags suffix.
     */
    public void setOutputSuffix(String outputSuffix) {
        this.outputSuffix = outputSuffix;
    }

    /**
     * Indicates whether the charge state from the spectrum should be used.
     *
     * @return boolean indicating whether the charge states from the spectrum
     * should be used
     */
    public boolean isUseChargeStateFromMS() {
        return useChargeStateFromMS;
    }

    /**
     * Sets whether the charge state from the spectrum should be used.
     *
     * @param useChargeStateFromMS boolean indicating whether the charge states
     * from the spectrum should be used
     */
    public void setUseChargeStateFromMS(boolean useChargeStateFromMS) {
        this.useChargeStateFromMS = useChargeStateFromMS;
    }

    /**
     * Indicates whether to use duplicate spectra option.
     *
     * @return boolean indicating whether to use duplicate spectra option
     */
    public boolean isDuplicateSpectra() {
        return duplicateSpectra;
    }

    /**
     * Sets whether to use duplicate spectra option.
     *
     * @param duplicateSpectra boolean indicating whether to use duplicate
     * spectra option
     */
    public void setDuplicateSpectra(boolean duplicateSpectra) {
        this.duplicateSpectra = duplicateSpectra;
    }

    /**
     * Returns the deisotoping flag value.
     *
     * @return deisotoping flag value
     */
    public int getDeisotopingMode() {
        return deisotopingMode;
    }

    /**
     * Sets the deisotoping flag value.
     *
     * @param deisotopingMode deisotoping flag value
     */
    public void setDeisotopingMode(int deisotopingMode) {
        this.deisotopingMode = deisotopingMode;
    }

    /**
     * Returns the deisotoping m/z tolerance.
     *
     * @return deisotoping m/z tolerance
     */
    public double getIsotopeMzTolerance() {
        return isotopeMzTolerance;
    }

    /**
     * Sets the deisotoping m/z tolerance.
     *
     * @param isotopeMzTolerance deisotoping m/z tolerance
     */
    public void setIsotopeMzTolerance(double isotopeMzTolerance) {
        this.isotopeMzTolerance = isotopeMzTolerance;
    }

    /**
     * Returns the complement m/z tolerance.
     *
     * @return the complement m/z tolerance
     */
    public double getComplementMzTolerance() {
        return complementMzTolerance;
    }

    /**
     * Sets the complement m/z tolerance.
     *
     * @param complementMzTolerance the complement m/z tolerance
     */
    public void setComplementMzTolerance(double complementMzTolerance) {
        this.complementMzTolerance = complementMzTolerance;
    }

    /**
     * Returns the tag length.
     *
     * @return the tag length
     */
    public int getTagLength() {
        return tagLength;
    }

    /**
     * Sets the tag length.
     *
     * @param tagLength the tag length
     */
    public void setTagLength(int tagLength) {
        this.tagLength = tagLength;
    }

    /**
     * Returns the maximum number of dynamic mods.
     *
     * @return the maximum number of dynamic mods
     */
    public int getMaxDynamicMods() {
        return maxDynamicMods;
    }

    /**
     * Sets the maximum number of dynamic mods.
     *
     * @param maxDynamicMods the maximum number of dynamic mods
     */
    public void setMaxDynamicMods(int maxDynamicMods) {
        this.maxDynamicMods = maxDynamicMods;
    }

    /**
     * Returns the maximum number of tags per spectrum to be generated.
     *
     * @return the maximum number of tags per spectrum to be generated
     */
    public int getMaxTagCount() {
        return maxTagCount;
    }

    /**
     * Sets the maximum number of tags per spectrum to be generated.
     *
     * @param maxTagCount the maximum number of tags per spectrum to be
     * generated
     */
    public void setMaxTagCount(int maxTagCount) {
        this.maxTagCount = maxTagCount;
    }

    /**
     * Returns the intensity score weight.
     *
     * @return the intensity score weight
     */
    public double getIntensityScoreWeight() {
        return intensityScoreWeight;
    }

    /**
     * Set the intensity score weight.
     *
     * @param intensityScoreWeight the intensity score weight
     */
    public void setIntensityScoreWeight(double intensityScoreWeight) {
        this.intensityScoreWeight = intensityScoreWeight;
    }

    /**
     * Returns the mzFidelity score weight.
     *
     * @return the mzFidelity score weight
     */
    public double getMzFidelityScoreWeight() {
        return mzFidelityScoreWeight;
    }

    /**
     * Set the mzFidelity score weight.
     *
     * @param mzFidelityScoreWeight the mzFidelity score weight
     */
    public void setMzFidelityScoreWeight(double mzFidelityScoreWeight) {
        this.mzFidelityScoreWeight = mzFidelityScoreWeight;
    }

    /**
     * Returns the complement score weight.
     *
     * @return the complement score weight
     */
    public double getComplementScoreWeight() {
        return complementScoreWeight;
    }

    /**
     * Set the complement score weight.
     *
     * @param complementScoreWeight the complement score weight
     */
    public void setComplementScoreWeight(double complementScoreWeight) {
        this.complementScoreWeight = complementScoreWeight;
    }

    /**
     * Returns the TIC cutoff in percent. Default is 85%.
     *
     * @return the ticCutoffPercentage
     */
    public double getTicCutoffPercentage() {
        return ticCutoffPercentage;
    }

    /**
     * Set the TIC cutoff in percent. Default is 85%.
     *
     * @param ticCutoffPercentage the ticCutoffPercentage to set
     */
    public void setTicCutoffPercentage(double ticCutoffPercentage) {
        this.ticCutoffPercentage = ticCutoffPercentage;
    }

    /**
     * Returns the maximum peak count.
     *
     * @return the maxPeakCount
     */
    public int getMaxPeakCount() {
        return maxPeakCount;
    }

    /**
     * Set the maximum peak count.
     *
     * @param maxPeakCount the maxPeakCount to set
     */
    public void setMaxPeakCount(int maxPeakCount) {
        this.maxPeakCount = maxPeakCount;
    }

    /**
     * Returns the number of intensity classes.
     *
     * @return the numIntensityClasses
     */
    public int getNumIntensityClasses() {
        return numIntensityClasses;
    }

    /**
     * Sets the number of intensity classes.
     *
     * @param numIntensityClasses the numIntensityClasses to set
     */
    public void setNumIntensityClasses(int numIntensityClasses) {
        this.numIntensityClasses = numIntensityClasses;
    }

    /**
     * Returns if the precursor mass is to be adjusted.
     *
     * @return the adjustPrecursorMass
     */
    public boolean isAdjustPrecursorMass() {
        if (adjustPrecursorMass) { // Default to false, DirecTag crashes otherwise...
            adjustPrecursorMass = false;
        }
        return adjustPrecursorMass;
    }

    /**
     * Set if the precursor mass is to be adjusted.
     *
     * @param adjustPrecursorMass the adjustPrecursorMass to set
     */
    public void setAdjustPrecursorMass(boolean adjustPrecursorMass) {
        this.adjustPrecursorMass = adjustPrecursorMass;
    }

    /**
     * Returns the minimum precursor adjustment.
     *
     * @return the minPrecursorAdjustment
     */
    public double getMinPrecursorAdjustment() {
        return minPrecursorAdjustment;
    }

    /**
     * Set the minimum precursor adjustment.
     *
     * @param minPrecursorAdjustment the minPrecursorAdjustment to set
     */
    public void setMinPrecursorAdjustment(double minPrecursorAdjustment) {
        this.minPrecursorAdjustment = minPrecursorAdjustment;
    }

    /**
     * Returns the maximum precursor adjustment.
     *
     * @return the maxPrecursorAdjustment
     */
    public double getMaxPrecursorAdjustment() {
        return maxPrecursorAdjustment;
    }

    /**
     * Set the maximum precursor adjustment.
     *
     * @param maxPrecursorAdjustment the maxPrecursorAdjustment to set
     */
    public void setMaxPrecursorAdjustment(double maxPrecursorAdjustment) {
        this.maxPrecursorAdjustment = maxPrecursorAdjustment;
    }

    /**
     * Returns the precursor adjustment step.
     *
     * @return the precursorAdjustmentStep
     */
    public double getPrecursorAdjustmentStep() {
        return precursorAdjustmentStep;
    }

    /**
     * Set the precursor adjustment step.
     *
     * @param PrecursorAdjustmentStep the precursorAdjustmentStep to set
     */
    public void setPrecursorAdjustmentStep(double PrecursorAdjustmentStep) {
        this.precursorAdjustmentStep = PrecursorAdjustmentStep;
    }

    /**
     * Sets the variable modifications searched. The order is the one used by DirecTag
     * and the name is the utilities one.
     *
     * @param variableModifications list of the names of the searched variable modifications
     */
    public void setModifications(ArrayList<String> variableModifications) {
        this.variableModifications = variableModifications;
    }

    public void setPtms(ArrayList<String> variableModifications) {
        setModifications(variableModifications);
    }

    /**
     * Returns the name of the modification indexed by the given index.
     *
     * @param index the index of the modification of interest
     *
     * @return the name of the modification of interest
     */
    public String getUtilitiesModificationName(int index) {
        if (variableModifications == null || variableModifications.isEmpty()) {
            throw new IllegalArgumentException("Variable modification index map not set for this DirecTag sequencing.");
        }
        if (index < 0 || index >= variableModifications.size()) {
            throw new IllegalArgumentException("Variable modification index " + index + " not found in mapping.");
        }
        return variableModifications.get(index);
    }

    @Override
    public Advocate getAlgorithm() {
        return Advocate.direcTag;
    }

    @Override
    public boolean equals(IdentificationAlgorithmParameter identificationAlgorithmParameter) {
        if (identificationAlgorithmParameter instanceof DirecTagParameters) {
            DirecTagParameters direcTagParameters = (DirecTagParameters) identificationAlgorithmParameter;
            if (ticCutoffPercentage != direcTagParameters.getTicCutoffPercentage()) {
                return false;
            }
            if (maxPeakCount != direcTagParameters.getMaxPeakCount()) {
                return false;
            }
            if (numIntensityClasses != direcTagParameters.getNumIntensityClasses()) {
                return false;
            }
            if (adjustPrecursorMass != direcTagParameters.isAdjustPrecursorMass()) {
                return false;
            }
            double diff = Math.abs(minPrecursorAdjustment - direcTagParameters.getMinPrecursorAdjustment());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(maxPrecursorAdjustment - direcTagParameters.getMaxPrecursorAdjustment());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(precursorAdjustmentStep - direcTagParameters.getPrecursorAdjustmentStep());
            if (diff > 0.0000000000001) {
                return false;
            }
            if (numChargeStates != direcTagParameters.getNumChargeStates()) {
                return false;
            }
            if (outputSuffix.equalsIgnoreCase(direcTagParameters.getOutputSuffix())) {
                if (!outputSuffix.isEmpty() && direcTagParameters.getOutputSuffix().isEmpty()) {
                    return false;
                }
            }
            if (useChargeStateFromMS != direcTagParameters.isUseChargeStateFromMS()) {
                return false;
            }
            if (duplicateSpectra != direcTagParameters.isDuplicateSpectra()) {
                return false;
            }
            if (deisotopingMode != direcTagParameters.getDeisotopingMode()) {
                return false;
            }
            diff = Math.abs(isotopeMzTolerance - direcTagParameters.getIsotopeMzTolerance());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(complementMzTolerance - direcTagParameters.getComplementMzTolerance());
            if (diff > 0.0000000000001) {
                return false;
            }
            if (tagLength != direcTagParameters.getTagLength()) {
                return false;
            }
            if (maxDynamicMods != direcTagParameters.getMaxDynamicMods()) {
                return false;
            }
            if (maxTagCount != direcTagParameters.getMaxTagCount()) {
                return false;
            }
            diff = Math.abs(intensityScoreWeight - direcTagParameters.getIntensityScoreWeight());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(mzFidelityScoreWeight - direcTagParameters.getMzFidelityScoreWeight());
            if (diff > 0.0000000000001) {
                return false;
            }
            diff = Math.abs(complementScoreWeight - direcTagParameters.getComplementScoreWeight());
            if (diff > 0.0000000000001) {
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

        output.append("TIC_CUTOFF_PERCENTAGE=");
        output.append(ticCutoffPercentage);
        output.append(newLine);
        output.append("MAX_PEAK_COUNT=");
        output.append(maxPeakCount);
        output.append(newLine);
        output.append("NUM_INTENSITY_CLASSES=");
        output.append(numIntensityClasses);
        output.append(newLine);
        output.append("ADJUST_PRECURSOR_MASS=");
        output.append(adjustPrecursorMass);
        output.append(newLine);
        output.append("MIN_PRECUSOR_ADJUSTMENT=");
        output.append(minPrecursorAdjustment);
        output.append(newLine);
        output.append("MAX_PRECUSOR_ADJUSTMENT=");
        output.append(maxPrecursorAdjustment);
        output.append(newLine);
        output.append("PRECUSOR_ADJUSTMENT_STEP=");
        output.append(precursorAdjustmentStep);
        output.append(newLine);
        output.append("NUM_CHARGE_STATES=");
        output.append(numChargeStates);
        output.append(newLine);
        output.append("OUTPUT_SUFFIX=");
        output.append(outputSuffix);
        output.append(newLine);
        output.append("USE_CHARGE_STATE_FROM_MS=");
        output.append(useChargeStateFromMS);
        output.append(newLine);
        output.append("DUPLICATE_SPECTRA=");
        output.append(duplicateSpectra);
        output.append(newLine);
        output.append("DEISOTOPING_MODE=");
        output.append(deisotopingMode);
        output.append(newLine);
        output.append("ISOTOPE_MZ_TOLERANCE=");
        output.append(isotopeMzTolerance);
        output.append(newLine);
        output.append("COMPLEMENT_MZ_TOLERANCE=");
        output.append(complementMzTolerance);
        output.append(newLine);
        output.append("TAG_LENGTH=");
        output.append(tagLength);
        output.append(newLine);
        output.append("MAX_DYNAMIC_MODS=");
        output.append(maxDynamicMods);
        output.append(newLine);
        output.append("MAX_TAG_COUNT=");
        output.append(maxTagCount);
        output.append(newLine);
        output.append("INTENSITY_SCORE_WEIGHT=");
        output.append(intensityScoreWeight);
        output.append(newLine);
        output.append("MZ_FIDELITY_SCORE_WEIGHT=");
        output.append(mzFidelityScoreWeight);
        output.append(newLine);
        output.append("COMPLEMENT_SCORE_WEIGHT=");
        output.append(complementScoreWeight);
        output.append(newLine);

        return output.toString();
    }
}
