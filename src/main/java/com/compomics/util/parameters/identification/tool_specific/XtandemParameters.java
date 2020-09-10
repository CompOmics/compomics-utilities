package com.compomics.util.parameters.identification.tool_specific;

import com.compomics.util.db.object.DbObject;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;

/**
 * The X!Tandem specific parameters.
 *
 * @author Marc Vaudel
 */
public class XtandemParameters extends DbObject implements IdentificationAlgorithmParameter {

    /**
     * Version number for deserialization.
     */
    static final long serialVersionUID = -5898951075262732261L;
    /**
     * Maximal e-value cut-off.
     */
    private double maxEValue = 0.01;
    /**
     * The output results filter: all, valid or stochastic.
     */
    private String outputResults = "all";
    /**
     * The dynamic range for spectrum filtering. When the highest peak is given
     * the dynamic range value peaks smaller than one are ignored. e.g. for 100
     * peaks with an intensity &lt;1% of the highest peak are ignored.
     */
    private double dynamicRange = 100.0;
    /**
     * The number of most intense peaks to consider.
     */
    private int nPeaks = 50;
    /**
     * The minimum precursor mass.
     */
    private double minPrecursorMass = 500.0;
    /**
     * The minimum fragment mass.
     */
    private double minFragmentMz = 200.0;
    /**
     * The minimum number of peaks per spectrum.
     */
    private int minPeaksPerSpectrum = 5;
    /**
     * Indicates whether the protein quick acetylation option should be
     * triggered.
     */
    private boolean proteinQuickAcetyl = true;
    /**
     * Indicates whether the quick pyrolidone option should be triggered.
     */
    private boolean quickPyrolidone = true;
    /**
     * Triggers the refinement process.
     */
    private boolean refine = true;
    /**
     * Sets whether semi enzymatic peptides should be search for during the
     * refinement process.
     */
    private boolean refineSemi = false;
    /**
     * Sets whether point mutations should be search for during the refinement
     * process.
     */
    private boolean refinePointMutations = false;
    /**
     * Sets whether the spectrum synthesis option should be used during the
     * refinement process.
     */
    private boolean refineSpectrumSynthesis = true;
    /**
     * Sets whether unexpected cleavages should be search for during the
     * refinement process.
     */
    private boolean refineUnanticipatedCleavages = true;
    /**
     * Indicates whether snAPs should be used during the refinement process.
     */
    private boolean refineSnaps = true;
    /**
     * The maximum expectation value for a hit to be considered during the
     * refinement process.
     */
    private double maximumExpectationValueRefinement = 0.01;
    /**
     * Sets the modifications to be used during the refinement process.
     */
    private boolean potentialModificationsForFullRefinment = false;
    /**
     * The skyline path.
     */
    private String skylinePath = "";
    /**
     * If true protein details will be exported to the to the result file.
     */
    private boolean outputProteins = true;
    /**
     * If true protein sequences will be added to the protein details to the
     * result file.
     */
    private boolean outputSequences = false;
    /**
     * If true spectra will be exported to the result file.
     */
    private boolean outputSpectra = true;
    /**
     * If true histograms will be exported to the result file.
     */
    private boolean outputHistograms = false;
    /**
     * Indicates whether the phospho stp bias option should be triggered.
     */
    private boolean stpBias = false;
    /**
     * Triggers the noise suppression function. Note: ignored in X!Tandem
     * VENGEANCE (2015.12.15) and newer.
     */
    private boolean useNoiseSuppression = false;
    /**
     * The value of the command "protein, ptm complexity" (C, a floating point
     * number 0.0–12.0) sets the maximum number of variable modification
     * alternatives that will be tested for a particular peptide. The number of
     * alternatives is 2.0C. If this number is not specified, the default value
     * C = 6.0 will be used.
     */
    private double proteinPtmComplexity = 6.0;
    /**
     * Sets whether the parent ion mass tolerance is expanded by opening up
     * multiple tolerance windows centered on the first and second 13C isotope
     * peaks for a peptide.
     */
    private boolean parentMonoisotopicMassIsotopeError = true;

    /**
     * Constructor.
     */
    public XtandemParameters() {

    }

    /**
     * Returns the dynamic range for spectrum filtering.
     *
     * @return the dynamic range for spectrum filtering
     */
    public double getDynamicRange() {
        return dynamicRange;
    }

    /**
     * Sets the dynamic range for spectrum filtering.
     *
     * @param dynamicRange the dynamic range for spectrum filtering
     */
    public void setDynamicRange(double dynamicRange) {
        this.dynamicRange = dynamicRange;
    }

    /**
     * Returns the number of most intense peaks to consider.
     *
     * @return the number of most intense peaks to consider
     */
    public int getnPeaks() {
        return nPeaks;
    }

    /**
     * Sets the number of most intense peaks to consider.
     *
     * @param nPeaks the number of most intense peaks to consider
     */
    public void setnPeaks(int nPeaks) {
        this.nPeaks = nPeaks;
    }

    /**
     * Returns the minimal precursor mass.
     *
     * @return the minimal precursor mass
     */
    public double getMinPrecursorMass() {
        return minPrecursorMass;
    }

    /**
     * Sets the minimal precursor mass.
     *
     * @param minPrecursorMass the minimal precursor mass
     */
    public void setMinPrecursorMass(double minPrecursorMass) {
        this.minPrecursorMass = minPrecursorMass;
    }

    /**
     * Returns the minimal fragment m/z.
     *
     * @return the minimal fragment m/z
     */
    public double getMinFragmentMz() {
        return minFragmentMz;
    }

    /**
     * Sets the minimal fragment m/z.
     *
     * @param minFragmentMz the minimal fragment m/z
     */
    public void setMinFragmentMz(double minFragmentMz) {
        this.minFragmentMz = minFragmentMz;
    }

    /**
     * Returns the minimal number of peaks per spectrum.
     *
     * @return the minimal number of peaks per spectrum
     */
    public int getMinPeaksPerSpectrum() {
        return minPeaksPerSpectrum;
    }

    /**
     * Sets the minimal number of peaks per spectrum.
     *
     * @param minPeaksPerSpectrum the minimal number of peaks per spectrum
     */
    public void setMinPeaksPerSpectrum(int minPeaksPerSpectrum) {
        this.minPeaksPerSpectrum = minPeaksPerSpectrum;
    }

    /**
     * Indicates whether the protein quick acetylation option should be
     * triggered.
     *
     * @return true if the protein quick acetylation option should be triggered
     */
    public boolean isProteinQuickAcetyl() {
        return proteinQuickAcetyl;
    }

    /**
     * Sets whether the protein quick acetylation option should be triggered.
     *
     * @param proteinQuickAcetyl true if the protein quick acetylation option
     * should be triggered
     */
    public void setProteinQuickAcetyl(boolean proteinQuickAcetyl) {
        this.proteinQuickAcetyl = proteinQuickAcetyl;
    }

    /**
     * Returns whether the quick pyrolidone option should be triggered.
     *
     * @return true if the quick pyrolidone option should be triggered
     */
    public boolean isQuickPyrolidone() {
        return quickPyrolidone;
    }

    /**
     * Sets whether the quick pyrolidone option should be triggered.
     *
     * @param quickPyrolidone the quick pyrolidone option should be triggered
     */
    public void setQuickPyrolidone(boolean quickPyrolidone) {
        this.quickPyrolidone = quickPyrolidone;
    }

    /**
     * Returns whether the second pass search should be triggered.
     *
     * @return true if the second pass search should be triggered
     */
    public boolean isRefine() {
        return refine;
    }

    /**
     * Sets whether the second pass search should be triggered.
     *
     * @param refine true if the second pass search should be triggered
     */
    public void setRefine(boolean refine) {
        this.refine = refine;
    }

    /**
     * Returns whether the stP bias should be triggered.
     *
     * @return true if the stP bias should be triggered
     */
    public boolean isStpBias() {
        return stpBias;
    }

    /**
     * Sets whether the stP bias should be triggered
     *
     * @param stpBias true if the stP bias should be triggered
     */
    public void setStpBias(boolean stpBias) {
        this.stpBias = stpBias;
    }

    /**
     * Returns the maximal e-value searched for.
     *
     * @return the maximal e-value searched for
     */
    public double getMaxEValue() {
        return maxEValue;
    }

    /**
     * Sets the maximal e-value searched for.
     *
     * @param maxEValue the maximal e-value searched for
     */
    public void setMaxEValue(double maxEValue) {
        this.maxEValue = maxEValue;
    }

    /**
     * Indicates whether the semi enzymatic option of the second pass search
     * should be triggered.
     *
     * @return true if the semi enzymatic option of the second pass search
     * should be triggered
     */
    public boolean isRefineSemi() {
        return refineSemi;
    }

    /**
     * Sets whether the semi enzymatic option of the second pass search should
     * be triggered.
     *
     * @param refineSemi true if the semi enzymatic option of the second pass
     * search should be triggered
     */
    public void setRefineSemi(boolean refineSemi) {
        this.refineSemi = refineSemi;
    }

    /**
     * Indicates whether point mutations should be looked for during the
     * refinement process.
     *
     * @return true if point mutations should be looked for during the
     * refinement process
     */
    public boolean isRefinePointMutations() {
        return refinePointMutations;
    }

    /**
     * Sets whether point mutations should be looked for during the refinement
     * process.
     *
     * @param refinePointMutations true if point mutations should be looked for
     * during the refinement process
     */
    public void setRefinePointMutations(boolean refinePointMutations) {
        this.refinePointMutations = refinePointMutations;
    }

    /**
     * Indicates whether the spectrum synthesis option should be used during the
     * refinement process.
     *
     * @return true if the spectrum synthesis option should be used during the
     * refinement process
     */
    public boolean isRefineSpectrumSynthesis() {
        return refineSpectrumSynthesis;
    }

    /**
     * Sets whether the spectrum synthesis option should be used during the
     * refinement process.
     *
     * @param refineSpectrumSynthesis true if the spectrum synthesis option
     * should be used during the refinement process
     */
    public void setRefineSpectrumSynthesis(boolean refineSpectrumSynthesis) {
        this.refineSpectrumSynthesis = refineSpectrumSynthesis;
    }

    /**
     * Returns whether the unanticipated cleavages option should be used during
     * the refinement process.
     *
     * @return true if the unanticipated cleavages option should be used during
     * the refinement process
     */
    public boolean isRefineUnanticipatedCleavages() {
        return refineUnanticipatedCleavages;
    }

    /**
     * Sets whether the unanticipated cleavages option should be used during the
     * refinement process.
     *
     * @param refineUnanticipatedCleavages true if the unanticipated cleavages
     * option should be used during the refinement process
     */
    public void setRefineUnanticipatedCleavages(boolean refineUnanticipatedCleavages) {
        this.refineUnanticipatedCleavages = refineUnanticipatedCleavages;
    }

    /**
     * Returns the maximum expectation value to use for refinement.
     *
     * @return the maximum expectation value to use for refinement
     */
    public double getMaximumExpectationValueRefinement() {
        return maximumExpectationValueRefinement;
    }

    /**
     * Sets the maximum expectation value to use for refinement.
     *
     * @param maximumExpectationValue the maximum expectation value to use for
     * refinement
     */
    public void setMaximumExpectationValueRefinement(double maximumExpectationValue) {
        this.maximumExpectationValueRefinement = maximumExpectationValue;
    }

    /**
     * Indicates whether the refinement modifications should be used for the
     * full refinement.
     *
     * @return true if the refinement modifications should be used for the full
     * refinement
     */
    public boolean isPotentialModificationsForFullRefinment() {
        return potentialModificationsForFullRefinment;
    }

    /**
     * Sets whether the refinement modifications should be used for the full
     * refinement
     *
     * @param potentialModificationsForFullRefinment true if the refinement
     * modifications should be used for the full refinement
     */
    public void setPotentialModificationsForFullRefinment(boolean potentialModificationsForFullRefinment) {
        this.potentialModificationsForFullRefinment = potentialModificationsForFullRefinment;
    }

    /**
     * Returns the skyline path.
     *
     * @return the skyline path
     */
    public String getSkylinePath() {
        return skylinePath;
    }

    /**
     * Sets the skyline path.
     *
     * @param skylinePath the skyline path
     */
    public void setSkylinePath(String skylinePath) {
        this.skylinePath = skylinePath;
    }

    /**
     * Indicates whether the protein bloc should be included in the export.
     *
     * @return true if the protein bloc should be included in the export
     */
    public boolean isOutputProteins() {
        return outputProteins;
    }

    /**
     * Sets whether the protein bloc should be included in the export.
     *
     * @param outputProteins the protein bloc should be included in the export
     */
    public void setOutputProteins(boolean outputProteins) {
        this.outputProteins = outputProteins;
    }

    /**
     * Returns whether the protein sequences should be included in the protein
     * block of the export.
     *
     * @return true if the protein sequences should be included in the protein
     * block of the export
     */
    public boolean isOutputSequences() {
        return outputSequences;
    }

    /**
     * Sets whether the protein sequences should be included in the protein
     * block of the export.
     *
     * @param outputSequences true if the protein sequences should be included
     * in the protein block of the export
     */
    public void setOutputSequences(boolean outputSequences) {
        this.outputSequences = outputSequences;
    }

    /**
     * Indicate whether the spectra should be exported in the result file.
     *
     * @return true if the spectra should be exported in the result file
     */
    public boolean isOutputSpectra() {
        return outputSpectra;
    }

    /**
     * Sets whether the spectra should be exported in the result file.
     *
     * @param outputSpectra true if the spectra should be exported in the result
     * file
     */
    public void setOutputSpectra(boolean outputSpectra) {
        this.outputSpectra = outputSpectra;
    }

    /**
     * Indicates whether histograms should be written in the result file.
     *
     * @return true if histograms should be written in the result file
     */
    public boolean isOutputHistograms() {
        return outputHistograms;
    }

    /**
     * Sets whether histograms should be written in the result file
     *
     * @param outputHistograms true if histograms should be written in the
     * result file
     */
    public void setOutputHistograms(boolean outputHistograms) {
        this.outputHistograms = outputHistograms;
    }

    /**
     * Indicates whether noise suppression should be used when importing
     * spectra. Note: ignored in X!Tandem VENGEANCE (2015.12.15) and newer
     *
     * @return true if noise suppression should be used when importing spectra
     */
    public boolean isUseNoiseSuppression() {
        return useNoiseSuppression;
    }

    /**
     * Sets whether noise suppression should be used when importing spectra.
     * Note: ignored in X!Tandem VENGEANCE (2015.12.15) and newer
     *
     * @param useNoiseSuppression true if noise suppression should be used when
     * importing spectra
     */
    public void setUseNoiseSuppression(boolean useNoiseSuppression) {
        this.useNoiseSuppression = useNoiseSuppression;
    }

    /**
     * Sets whether snAPs should be used during the refinement process.
     *
     * @return true if snAPs should be used during the refinement process
     */
    public boolean isRefineSnaps() {
        return refineSnaps;
    }

    /**
     * Sets whether snAPs should be used during the refinement process.
     *
     * @param refineSnaps true if snAPs should be used during the refinement
     * process
     */
    public void setRefineSnaps(boolean refineSnaps) {
        this.refineSnaps = refineSnaps;
    }

    /**
     * Returns the output results filter.
     *
     * @return the outputResults
     */
    public String getOutputResults() {
        return outputResults;
    }

    /**
     * Set the output results filter.
     *
     * @param outputResults the outputResults to set
     */
    public void setOutputResults(String outputResults) {
        this.outputResults = outputResults;
    }

    /**
     * Returns the proteinPtmComplexity. 
     * 
     * @return the proteinPtmComplexity
     */
    public double getProteinPtmComplexity() {
        return proteinPtmComplexity;
    }

    /**
     * Set the proteinPtmComplexity.
     * 
     * @param proteinPtmComplexity the proteinPtmComplexity to set
     */
    public void setProteinPtmComplexity(double proteinPtmComplexity) {
        this.proteinPtmComplexity = proteinPtmComplexity;
    }

    @Override
    public Advocate getAlgorithm() {
        return Advocate.xtandem;
    }

    @Override
    public boolean equals(IdentificationAlgorithmParameter identificationAlgorithmParameter) {

        if (identificationAlgorithmParameter instanceof XtandemParameters) {
            XtandemParameters xtandemParameters = (XtandemParameters) identificationAlgorithmParameter;
            if (maxEValue != xtandemParameters.getMaxEValue()) {
                return false;
            }
            if (dynamicRange != xtandemParameters.getDynamicRange()) {
                return false;
            }
            if (getnPeaks() != xtandemParameters.getnPeaks()) {
                return false;
            }
            if (minPrecursorMass != xtandemParameters.getMinPrecursorMass()) {
                return false;
            }
            if (minFragmentMz != xtandemParameters.getMinFragmentMz()) {
                return false;
            }
            if (getMinPeaksPerSpectrum() != xtandemParameters.getMinPeaksPerSpectrum()) {
                return false;
            }
            if (isProteinQuickAcetyl() != xtandemParameters.isProteinQuickAcetyl()) {
                return false;
            }
            if (isQuickPyrolidone() != xtandemParameters.isQuickPyrolidone()) {
                return false;
            }
            if (isRefine() != xtandemParameters.isRefine()) {
                return false;
            }
            if (isRefineSemi() != xtandemParameters.isRefineSemi()) {
                return false;
            }
            if (isRefinePointMutations() != xtandemParameters.isRefinePointMutations()) {
                return false;
            }
            if (isRefineSpectrumSynthesis() != xtandemParameters.isRefineSpectrumSynthesis()) {
                return false;
            }
            if (isRefineUnanticipatedCleavages() != xtandemParameters.isRefineUnanticipatedCleavages()) {
                return false;
            }
            if (isRefineSnaps() != xtandemParameters.isRefineSnaps()) {
                return false;
            }
            if (maximumExpectationValueRefinement != xtandemParameters.getMaximumExpectationValueRefinement()) {
                return false;
            }
            if (isPotentialModificationsForFullRefinment() != xtandemParameters.isPotentialModificationsForFullRefinment()) {
                return false;
            }
            if (!getSkylinePath().equals(xtandemParameters.getSkylinePath())) {
                return false;
            }
            if (isOutputProteins() != xtandemParameters.isOutputProteins()) {
                return false;
            }
            if (isOutputSpectra() != xtandemParameters.isOutputSpectra()) {
                return false;
            }
            if (isOutputSequences() != xtandemParameters.isOutputSequences()) {
                return false;
            }
            if (isOutputHistograms() != xtandemParameters.isOutputHistograms()) {
                return false;
            }
            if (isStpBias() != xtandemParameters.isStpBias()) {
                return false;
            }
            if (isUseNoiseSuppression() != xtandemParameters.isUseNoiseSuppression()) {
                return false;
            }
            if (!getOutputResults().equalsIgnoreCase(xtandemParameters.getOutputResults())) {
                return false;
            }
            if (getProteinPtmComplexity() != xtandemParameters.getProteinPtmComplexity()) {
                return false;
            }
            if (getParentMonoisotopicMassIsotopeError() != xtandemParameters.getParentMonoisotopicMassIsotopeError()) {
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

        output.append("DYNAMIC_RANGE=");
        output.append(dynamicRange);
        output.append(newLine);

        output.append("NUMBER_OF_PEAKS=");
        output.append(nPeaks);
        output.append(newLine);

        output.append("MIN_FRAG_MASS=");
        output.append(minFragmentMz);
        output.append(newLine);

        output.append("MIN_NUMBER_OF_PEAKS=");
        output.append(minPeaksPerSpectrum);
        output.append(newLine);

        output.append("NOISE_SUPPRESSION=");
        if (useNoiseSuppression) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);
        
        output.append("PARENT_MONOISOTOPIC_MASS_ISOTOPE_ERROR=");
        if (parentMonoisotopicMassIsotopeError) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("MIN_PREC_MASS=");
        output.append(minPrecursorMass);
        output.append(newLine);

        output.append("PROTEIN_QUICK_ACETYL=");
        if (proteinQuickAcetyl) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("QUICK_PYROLIDONE=");
        if (quickPyrolidone) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);
        
        output.append("PROTEIN_PTM_COMPLEXITY=");
        output.append(getProteinPtmComplexity());
        output.append(newLine);

        output.append("STP_BIAS=");
        if (stpBias) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE=");
        if (refine) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE_SEMI=");
        if (refineSemi) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE_POINT_MUTATIONS=");
        if (refinePointMutations) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE_SPECTRUM_SYNTHESIS=");
        if (refineSpectrumSynthesis) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE_UNANTICIPATED_CLEABAGES=");
        if (refineUnanticipatedCleavages) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE_SNAPS=");
        if (refineSnaps) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("REFINE_MAX_EVALUE=");
        output.append(maximumExpectationValueRefinement);
        output.append(newLine);

        output.append("POTENTIAL_MODIFICATIONS_FOR_FULL_REFINEMENT=");
        if (potentialModificationsForFullRefinment) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("EVALUE_CUTOFF=");
        output.append(maxEValue);
        output.append(newLine);

        output.append("SKYLINE_PATH=");
        output.append(skylinePath);
        output.append(newLine);

        output.append("OUTPUT_RESULTS=");
        output.append(getOutputResults());
        output.append(newLine);

        output.append("OUTPUT_PROTEINS=");
        if (outputProteins) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("OUTPUT_SEQUENCES=");
        if (outputSequences) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("OUTPUT_SPECTRA=");
        if (outputSpectra) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        output.append("OUTPUT_HISTOGRAMS=");
        if (outputHistograms) {
            output.append("YES");
        } else {
            output.append("NO");
        }
        output.append(newLine);

        return output.toString();
    }

    /**
     * Returns true if the parent ion mass tolerance is expanded by opening up
     * multiple tolerance windows centered on the first and second 13C isotope
     * peaks for a peptide.
     * 
     * @return the parentMonoisotopicMassIsotopeError
     */
    public boolean getParentMonoisotopicMassIsotopeError() {
        return parentMonoisotopicMassIsotopeError;
    }

    /**
     * Sets whether the parent ion mass tolerance is expanded by opening up
     * multiple tolerance windows centered on the first and second 13C isotope
     * peaks for a peptide.
     * 
     * @param parentMonoisotopicMassIsotopeError the parentMonoisotopicMassIsotopeError to set
     */
    public void setParentMonoisotopicMassIsotopeError(boolean parentMonoisotopicMassIsotopeError) {
        this.parentMonoisotopicMassIsotopeError = parentMonoisotopicMassIsotopeError;
    }
}
