package com.compomics.util.parameters.identification.advanced;

import com.compomics.util.experiment.personalization.ExperimentObject;

/**
 * The sequence matching options.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class SequenceMatchingParameters extends ExperimentObject {

    /**
     * The different types of amino acid matching.
     */
    public static enum MatchingType {

        /**
         * Matches character strings only.
         */
        string(0, "Character Sequence"),
        /**
         * Matches amino acids.
         */
        aminoAcid(1, "Amino Acids"),
        /**
         * Matches amino acids of indistinguishable masses (I/L).
         */
        indistiguishableAminoAcids(2, "Indistinguishable Amino Acids");

        /**
         * The index of the type as integer.
         */
        public final int index;
        /**
         * The description.
         */
        public final String description;

        /**
         * Constructor.
         *
         * @param index the index
         * @param description the description
         */
        private MatchingType(int index, String description) {
            this.index = index;
            this.description = description;
        }

        @Override
        public String toString() {
            return description;
        }

        /**
         * Returns the different matching types as command line options
         * description.
         *
         * @return the different matching types as command line options
         * description
         */
        public static String getCommandLineOptions() {

            StringBuilder optionsStringBuilder = new StringBuilder();

            for (MatchingType matchingType : values()) {

                if (optionsStringBuilder.length() != 0) {
                    optionsStringBuilder.append(", ");
                }

                optionsStringBuilder.append(matchingType.index).append(": ").append(matchingType.description);
            }

            return optionsStringBuilder.toString();
        }

        /**
         * Returns the matching type corresponding to the given index.
         *
         * @param index the index of the matching type
         *
         * @return the matching type
         */
        public static MatchingType getMatchingType(int index) {

            for (MatchingType matchingType : values()) {

                if (matchingType.index == index) {
                    return matchingType;
                }

            }

            throw new IllegalArgumentException("No matching type found for index " + index + ".");
        }
    }

    /**
     * The amino acid matching type.
     */
    private MatchingType sequenceMatchingType;
    /**
     * Limit the share of X's a match can contain, range [0.0-1.0].
     */
    private double limitX = 0.25;
    /**
     * Boolean deciding whether tags should only be mapped to enzymatic
     * peptides.
     */
    private Boolean enzymaticTagsOnly = false;
    /**
     * The maximum number of PTMs per peptide when mapping tags.
     */
    private int maxPtmsPerTagPeptide = 3;
    /**
     * The minimum amino acid score [0-100]. Used when converting Novor peptides
     * into tags.
     */
    private Integer minAminoAcidScore = 30;
    /**
     * The minimum tag length. Used when converting Novor peptides into tags.
     */
    private Integer minTagLength = 3;
    /**
     * Default string matching.
     */
    public static final SequenceMatchingParameters DEFAULT_STRING_MATCHING = getStringMatching();

    public static final SequenceMatchingParameters defaultStringMatching = DEFAULT_STRING_MATCHING;

    /**
     * Constructor for empty preferences.
     */
    public SequenceMatchingParameters() {

    }

    /**
     * Returns preferences for simple string matching.
     *
     * @return preferences for simple string matching
     */
    public static SequenceMatchingParameters getStringMatching() {

        SequenceMatchingParameters sequenceMatchingPreferences = new SequenceMatchingParameters();
        sequenceMatchingPreferences.setSequenceMatchingType(MatchingType.string);
        return sequenceMatchingPreferences;

    }

    /**
     * Returns default preferences from amino acid matching.
     *
     * @return default preferences from amino acid matching
     */
    public static SequenceMatchingParameters getDefaultSequenceMatching() {

        SequenceMatchingParameters sequenceMatchingPreferences = new SequenceMatchingParameters();
        sequenceMatchingPreferences.setSequenceMatchingType(MatchingType.indistiguishableAminoAcids);
        sequenceMatchingPreferences.setLimitX(0.25);
        sequenceMatchingPreferences.setMaxPtmsPerTagPeptide(3);
        return sequenceMatchingPreferences;

    }

    /**
     * Returns the sequence matching type.
     *
     * @return the sequence matching type
     */
    public MatchingType getSequenceMatchingType() {

        return sequenceMatchingType;

    }

    /**
     * Sets the sequence matching type.
     *
     * @param sequenceMatchingType the sequence matching type
     */
    public void setSequenceMatchingType(MatchingType sequenceMatchingType) {

        this.sequenceMatchingType = sequenceMatchingType;

    }

    /**
     * Returns the maximal share of X's a match can contain, range [0.0-1.0].
     * Null if not set.
     *
     * @return the maximal share of X's a match can contain
     */
    public double getLimitX() {

        return limitX;

    }

    /**
     * Sets the maximal share of X's a match can contain, range [0.0-1.0].
     *
     * @param limitX the maximal share of X's a match can contain
     */
    public void setLimitX(double limitX) {

        this.limitX = limitX;

    }

    /**
     * Indicates whether another sequence matching preferences is the same as
     * this one.
     *
     * @param sequenceMatchingPreferences the other sequence matching
     * preferences
     *
     * @return whether another sequence matching preferences is the same as this
     * one
     */
    public boolean isSameAs(SequenceMatchingParameters sequenceMatchingPreferences) {

        if (sequenceMatchingType != sequenceMatchingPreferences.getSequenceMatchingType()) {
            return false;
        }

        if (limitX != sequenceMatchingPreferences.getLimitX()) {
            return false;
        }

        if (isEnzymaticTagsOnly() != sequenceMatchingPreferences.isEnzymaticTagsOnly()) {
            return false;
        }

        if (maxPtmsPerTagPeptide != sequenceMatchingPreferences.getMaxPtmsPerTagPeptide()) {
            return false;
        }

        if (getMinAminoAcidScore() != sequenceMatchingPreferences.getMinAminoAcidScore()) {
            return false;
        }

        if (getMinTagLength() != sequenceMatchingPreferences.getMinTagLength()) {
            return false;
        }

        return true;
    }

    /**
     * Returns a short description of the parameters.
     *
     * @return a short description of the parameters
     */
    public String getShortDescription() {

        String newLine = System.getProperty("line.separator");

        StringBuilder output = new StringBuilder();

        output.append("Method: ").append(sequenceMatchingType).append(".").append(newLine);
        output.append("Max share of x's: ").append(limitX).append(".").append(newLine);
        output.append("Enzymatic tags matching: ").append(isEnzymaticTagsOnly()).append(".").append(newLine);
        output.append("Max PTMs per tag: ").append(getMaxPtmsPerTagPeptide()).append(".");
        output.append("Min amino acid score: ").append(getMinAminoAcidScore()).append(".");
        output.append("Min tag length: ").append(getMinTagLength()).append(".");

        return output.toString();
    }

    /**
     * Returns true if tags should only be mapped to enzymatic peptides.
     *
     * @return true if tags should only be mapped to enzymatic peptides
     */
    public boolean isEnzymaticTagsOnly() {

        if (enzymaticTagsOnly == null) {
            enzymaticTagsOnly = false;
        }

        return enzymaticTagsOnly;
    }

    /**
     * Sets whether tags should only be mapped to enzymatic peptides.
     *
     * @param enzymaticTagsOnly the enzymaticTagsOnly to set
     */
    public void setEnzymaticTagsOnly(boolean enzymaticTagsOnly) {

        this.enzymaticTagsOnly = enzymaticTagsOnly;

    }

    /**
     * Returns the maximum number of PTMs to consider when mapping tags to
     * protein sequences.
     *
     * @return the maximum number of PTMs to consider when mapping tags to
     * protein sequences
     */
    public int getMaxPtmsPerTagPeptide() {

        return maxPtmsPerTagPeptide;

    }

    /**
     * Sets the maximum number of PTMs to consider when mapping tags to protein
     * sequences.
     *
     * @param numberOfPtmsPerTagPeptide the maxPtmsPerTagPeptide to set
     */
    public void setMaxPtmsPerTagPeptide(int numberOfPtmsPerTagPeptide) {

        this.maxPtmsPerTagPeptide = numberOfPtmsPerTagPeptide;

    }

    /**
     * Returns the minimum amino acid score.
     *
     * @return the minAminoAcidScore
     */
    public int getMinAminoAcidScore() {

        if (minAminoAcidScore == null) {
            minAminoAcidScore = 30;
        }
        
        return minAminoAcidScore;

    }

    /**
     * Set the minimum amino acid score.
     *
     * @param minAminoAcidScore the minAminoAcidScore to set
     */
    public void setMinAminoAcidScore(int minAminoAcidScore) {

        this.minAminoAcidScore = minAminoAcidScore;

    }

    /**
     * Returns the minimum tag length.
     *
     * @return the minTagLength
     */
    public int getMinTagLength() {

        if (minAminoAcidScore == null) {
            minAminoAcidScore = 3;
        }
        
        return minTagLength;

    }

    /**
     * Set the minimum tag length.
     *
     * @param minTagLength the minTagLength to set
     */
    public void setMinTagLength(int minTagLength) {

        this.minTagLength = minTagLength;

    }
}
