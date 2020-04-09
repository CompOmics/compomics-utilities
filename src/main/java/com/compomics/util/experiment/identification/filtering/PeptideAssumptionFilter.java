package com.compomics.util.experiment.identification.filtering;

import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.Util;
import com.compomics.util.db.object.DbObject;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.utils.ProteinUtils;
import com.compomics.util.experiment.io.biology.protein.SequenceProvider;
import com.compomics.util.experiment.mass_spectrometry.SpectrumProvider;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import java.util.Arrays;
import java.util.TreeMap;

/**
 * This class filters peptide assumptions based on various properties.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class PeptideAssumptionFilter extends DbObject {

    /**
     * The minimal peptide length allowed.
     */
    private int minPepLength;
    /**
     * The maximal peptide length allowed.
     */
    private int maxPepLength;
    /**
     * The maximal m/z deviation allowed.
     */
    private double maxMassDeviation;
    /**
     * Boolean indicating the unit of the allowed m/z deviation (true: ppm,
     * false: Da).
     */
    private boolean isPpm;
    /**
     * Boolean indicating whether peptides presenting unknown modifications
     * should be ignored.
     */
    private boolean unknownModification;
    /**
     * The minimum number of missed cleavages allowed. Null means no lower
     * limit.
     */
    private Integer minMissedCleavages;
    /**
     * The maximum number of missed cleavages allowed. Null means no upper
     * limit.
     */
    private Integer maxMissedCleavages;
    /**
     * The minimum number of isotopes allowed. Null means no lower limit.
     */
    private Integer minIsotopes;
    /**
     * The maximum number of isotopes allowed. Null means no upper limit.
     */
    private Integer maxIsotopes;

    /**
     * Constructor with default settings.
     */
    public PeptideAssumptionFilter() {
        minPepLength = 8;
        maxPepLength = 30;
        maxMassDeviation = -1;
        isPpm = true;
        unknownModification = true;
        minMissedCleavages = null;
        maxMissedCleavages = null;
        minIsotopes = null;
        maxIsotopes = null;
    }

    /**
     * Constructor for an Identification filter.
     *
     * @param minPepLength the minimal peptide length allowed (0 or less for
     * disabled)
     * @param maxPepLength the maximal peptide length allowed (0 or less for
     * disabled)
     * @param maxMzDeviation the maximal m/z deviation allowed (0 or less for
     * disabled)
     * @param isPpm boolean indicating the unit of the allowed m/z deviation
     * (true: ppm, false: Da)
     * @param unknownModification shall peptides presenting unknown
     * modifications be removed
     * @param minMissedCleavages the minimum number of missed cleavages allowed
     * (null for disabled)
     * @param maxMissedCleavages the maximum number of missed cleavages allowed
     * (null for disabled)
     * @param minIsotopes the minimum number of isotopes allowed (null for
     * disabled)
     * @param maxIsotopes the maximum number of isotopes allowed (null for
     * disabled)
     */
    public PeptideAssumptionFilter(
            int minPepLength, 
            int maxPepLength, 
            double maxMzDeviation, 
            boolean isPpm, 
            boolean unknownModification, 
            Integer minMissedCleavages, 
            Integer maxMissedCleavages, 
            Integer minIsotopes, 
            Integer maxIsotopes
    ) {

        this.minPepLength = minPepLength;
        this.maxPepLength = maxPepLength;
        this.maxMassDeviation = maxMzDeviation;
        this.isPpm = isPpm;
        this.unknownModification = unknownModification;
        this.minMissedCleavages = minMissedCleavages;
        this.maxMissedCleavages = maxMissedCleavages;
        this.minIsotopes = minIsotopes;
        this.maxIsotopes = maxIsotopes;

    }

    /**
     * Updates the filter based on the search parameters.
     *
     * @param searchParameters the search parameters where to take the
     * information from
     */
    public void setFilterFromSearchParameters(
            SearchParameters searchParameters
    ) {
        writeDBMode();

        this.isPpm = searchParameters.isPrecursorAccuracyTypePpm();
        this.maxMassDeviation = searchParameters.getPrecursorAccuracy();
        this.minIsotopes = searchParameters.getMinIsotopicCorrection();
        this.maxIsotopes = searchParameters.getMaxIsotopicCorrection();
        this.unknownModification = true;

    }

    /**
     * Validates the peptide based on the peptide length, the share of X's in
     * the sequence and the allowed number of missed cleavages.
     *
     * @param peptide the peptide to validate
     * @param sequenceMatchingPreferences the sequence matching preferences
     * containing the maximal share of X's allowed
     * @param digestionPreferences the digestion preferences
     *
     * @return a boolean indicating whether the peptide passed the test
     */
    public boolean validatePeptide(
            Peptide peptide, 
            SequenceMatchingParameters sequenceMatchingPreferences, 
            DigestionParameters digestionPreferences
    ) {
        
        readDBMode();

        String peptideSequence = peptide.getSequence();
        int sequenceLength = peptideSequence.length();

        if ((maxPepLength > 0 && sequenceLength > maxPepLength)
                || (minPepLength > 0 && sequenceLength < minPepLength)) {

            return false;

        }

        double xShare = ((double) Util.getOccurrence(peptideSequence, 'X')) / sequenceLength;

        if (xShare > sequenceMatchingPreferences.getLimitX()) {

            return false;

        }

        if (minMissedCleavages != null || maxMissedCleavages != null) {

            int peptideMinMissedCleavages = peptide.getNMissedCleavages(digestionPreferences);

            if (minMissedCleavages != null && peptideMinMissedCleavages < minMissedCleavages) {

                return false;

            }

            if (maxMissedCleavages != null && peptideMinMissedCleavages > maxMissedCleavages) {

                return false;

            }
        }

        return true;
    }

    /**
     * Validates a peptide depending on its protein inference status.
     *
     * @param peptide the peptide
     * @param sequenceMatchingPreferences the sequence matching preferences
     * @param sequenceProvider a sequence provider
     *
     * @return a boolean indicating whether the peptide passed the test
     */
    public boolean validateProteins(
            Peptide peptide, 
            SequenceMatchingParameters sequenceMatchingPreferences, 
            SequenceProvider sequenceProvider
    ) {

        readDBMode();
        TreeMap<String, int[]> proteinMapping = peptide.getProteinMapping();

        if (proteinMapping != null && proteinMapping.size() > 1) {

            boolean target = false;
            boolean decoy = false;

            for (String accession : proteinMapping.navigableKeySet()) {

                if (ProteinUtils.isDecoy(accession, sequenceProvider)) {

                    decoy = true;

                } else {

                    target = true;

                }
            }

            if (target && decoy) {

                return false;

            }
        }

        return true;
    }

    /**
     * Verifies that the definition of every modification name is available.
     *
     * @param peptide the peptide of interest
     * @param sequenceMatchingPreferences the sequence matching preferences for
     * peptide to protein mapping
     * @param modificationSequenceMatchingPreferences the sequence matching
     * preferences for modification to peptide mapping
     * @param modificationProfile the modification profile of the identification
     *
     * @return a boolean indicating whether the peptide passed the test
     */
    public boolean validateModifications(
            Peptide peptide, 
            SequenceMatchingParameters sequenceMatchingPreferences,
            SequenceMatchingParameters modificationSequenceMatchingPreferences, 
            ModificationParameters modificationProfile
    ) {

        readDBMode();
        ModificationFactory modificationFactory = ModificationFactory.getInstance();

        // check if a modification could not be parsed
        if (unknownModification) {

            ModificationMatch[] modificationMatches = peptide.getVariableModifications();

            if (Arrays.stream(modificationMatches)
                    .map(ModificationMatch::getModification)
                    .anyMatch(modName -> !modificationFactory.containsModification(modName))) {

                return false;

            }
        }

        return true;
    }

    /**
     * Validates the mass deviation of a peptide assumption.
     *
     * @param assumption the considered peptide assumption
     * @param spectrumFile the file of the spectrum used to get the precursor
     * @param spectrumTitle the file of the spectrum used to get the precursor
     * @param spectrumProvider the spectrum provider
     * @param searchParameters the search parameters
     *
     * @return a boolean indicating whether the given assumption passes the
     * filter
     */
    public boolean validatePrecursor(
            PeptideAssumption assumption, 
            String spectrumFile, 
            String spectrumTitle, 
            SpectrumProvider spectrumProvider, 
            SearchParameters searchParameters
    ) {

        readDBMode();
        double precursorMz = spectrumProvider.getPrecursorMz(
                spectrumFile, 
                spectrumTitle
        );
        int isotopeNumber = assumption.getIsotopeNumber(
                precursorMz, 
                searchParameters.getMinIsotopicCorrection(), 
                searchParameters.getMaxIsotopicCorrection()
        );
        
        if (minIsotopes != null && isotopeNumber < minIsotopes) {
            return false;
        }
        
        if (maxIsotopes != null && isotopeNumber > maxIsotopes) {
            return false;
        }
        
        double mzDeviation = assumption.getDeltaMass(
                precursorMz, 
                isPpm, 
                searchParameters.getMinIsotopicCorrection(), 
                searchParameters.getMaxIsotopicCorrection()
        );
        
        return (maxMassDeviation <= 0 || Math.abs(mzDeviation) <= maxMassDeviation);
    }

    /**
     * Returns a boolean indicating whether unknown modifications shall be
     * removed.
     *
     * @return a boolean indicating whether unknown modifications shall be
     * removed
     */
    public boolean removeUnknownModifications() {
        readDBMode();

        return unknownModification;

    }

    /**
     * Set whether unknown modifications shall be removed.
     *
     * @param unknownModification whether unknown modifications shall be removed
     */
    public void setRemoveUnknownModifications(
            boolean unknownModification
    ) {
        writeDBMode();
        this.unknownModification = unknownModification;

    }

    /**
     * Indicates whether the mass tolerance is in ppm (true) or Dalton (false).
     *
     * @return a boolean indicating whether the mass tolerance is in ppm (true)
     * or Dalton (false)
     */
    public boolean isIsPpm() {
        readDBMode();
        return isPpm;
    }

    /**
     * Sets whether the mass tolerance is in ppm (true) or Dalton (false).
     *
     * @param isPpm a boolean indicating whether the mass tolerance is in ppm
     * (true) or Dalton (false)
     */
    public void setIsPpm(
            boolean isPpm
    ) {
        writeDBMode();
        this.isPpm = isPpm;
    }

    /**
     * Returns the maximal m/z deviation allowed.
     *
     * @return the maximal mass deviation allowed
     */
    public double getMaxMzDeviation() {
        readDBMode();
        return maxMassDeviation;
    }

    /**
     * Sets the maximal m/z deviation allowed.
     *
     * @param maxMzDeviation the maximal mass deviation allowed
     */
    public void setMaxMzDeviation(
            double maxMzDeviation
    ) {
        writeDBMode();
        this.maxMassDeviation = maxMzDeviation;
    }

    /**
     * Returns the maximal peptide length allowed.
     *
     * @return the maximal peptide length allowed
     */
    public int getMaxPepLength() {
        readDBMode();
        return maxPepLength;
    }

    /**
     * Sets the maximal peptide length allowed.
     *
     * @param maxPepLength the maximal peptide length allowed
     */
    public void setMaxPepLength(
            int maxPepLength
    ) {
        writeDBMode();
        this.maxPepLength = maxPepLength;
    }

    /**
     * Returns the maximal peptide length allowed.
     *
     * @return the maximal peptide length allowed
     */
    public int getMinPepLength() {
        readDBMode();
        return minPepLength;
    }

    /**
     * Sets the maximal peptide length allowed.
     *
     * @param minPepLength the maximal peptide length allowed
     */
    public void setMinPepLength(
            int minPepLength
    ) {
        writeDBMode();
        this.minPepLength = minPepLength;
    }

    /**
     * Returns the minimal number of isotopes allowed (inclusive).
     *
     * @return the minimal number of isotopes allowed
     */
    public Integer getMinIsotopes() {
        readDBMode();
        return minIsotopes;
    }

    /**
     * Sets the minimal number of isotopes allowed (inclusive).
     *
     * @param minIsotopes the minimal number of isotopes allowed
     */
    public void setMinIsotopes(
            Integer minIsotopes
    ) {
        writeDBMode();
        this.minIsotopes = minIsotopes;
    }

    /**
     * Returns the maximal number of isotopes allowed (inclusive).
     *
     * @return the maximal number of isotopes allowed
     */
    public Integer getMaxIsotopes() {
        readDBMode();
        return maxIsotopes;
    }

    /**
     * Sets the maximal number of isotopes allowed (inclusive).
     *
     * @param maxIsotopes the maximal number of isotopes allowed
     */
    public void setMaxIsotopes(
            Integer maxIsotopes
    ) {
        writeDBMode();
        this.maxIsotopes = maxIsotopes;
    }

    /**
     * Indicates whether this filter is the same as another one.
     *
     * @param anotherFilter another filter
     * 
     * @return a boolean indicating that the filters have the same parameters
     */
    public boolean isSameAs(
            PeptideAssumptionFilter anotherFilter
    ) {
        readDBMode();

        if (minMissedCleavages != null && anotherFilter.getMinMissedCleavages() != null) {
            if (!minMissedCleavages.equals(anotherFilter.getMinMissedCleavages())) {
                return false;
            }
        }
        if (minMissedCleavages != null && anotherFilter.getMinMissedCleavages() == null) {
            return false;
        }
        if (minMissedCleavages == null && anotherFilter.getMinMissedCleavages() != null) {
            return false;
        }
        if (maxMissedCleavages != null && anotherFilter.getMaxMissedCleavages() != null) {
            if (maxMissedCleavages.equals(anotherFilter.getMaxMissedCleavages())) {
                return false;
            }
        }
        if (maxMissedCleavages != null && anotherFilter.getMaxMissedCleavages() == null) {
            return false;
        }
        if (maxMissedCleavages == null && anotherFilter.getMaxMissedCleavages() != null) {
            return false;
        }

        if (minIsotopes != null && anotherFilter.getMinIsotopes() != null) {
            if (!minIsotopes.equals(anotherFilter.getMinIsotopes())) {
                return false;
            }
        }
        if (minIsotopes != null && anotherFilter.getMinIsotopes() == null) {
            return false;
        }
        if (minIsotopes == null && anotherFilter.getMinIsotopes() != null) {
            return false;
        }
        if (maxIsotopes != null && anotherFilter.getMaxIsotopes() != null) {
            if (!maxIsotopes.equals(anotherFilter.getMaxIsotopes())) {
                return false;
            }
        }
        if (maxIsotopes != null && anotherFilter.getMaxIsotopes() == null) {
            return false;
        }
        if (maxIsotopes == null && anotherFilter.getMaxIsotopes() != null) {
            return false;
        }

        return isPpm == anotherFilter.isPpm
                && unknownModification == anotherFilter.removeUnknownModifications()
                && minPepLength == anotherFilter.getMinPepLength()
                && maxPepLength == anotherFilter.getMaxPepLength()
                && maxMassDeviation == anotherFilter.getMaxMzDeviation();
    }

    /**
     * Returns a short description of the parameters.
     *
     * @return a short description of the parameters
     */
    public String getShortDescription() {
        readDBMode();

        String newLine = System.getProperty("line.separator");

        StringBuilder output = new StringBuilder();

        output.append("Peptide Length: ").append(minPepLength).append("-").append(maxPepLength).append(".").append(newLine);
        if (maxMassDeviation >= 0) {
            output.append("Precursor m/z Deviation: ").append(maxMassDeviation);
            if (isPpm) {
                output.append(" ppm.").append(newLine);
            } else {
                output.append(" Da.").append(newLine);
            }
        }
        output.append("Remove Unknown Modifications: ").append(unknownModification).append(".").append(newLine);

        if (minMissedCleavages != null || maxMissedCleavages != null) {

            output.append("Missed Cleavages: ");

            if (minMissedCleavages != null) {
                output.append(minMissedCleavages);
            } else {
                output.append("0");
            }

            output.append("-");

            if (maxMissedCleavages != null) {
                output.append(maxMissedCleavages);
            } else {
                output.append("n");
            }

            output.append(".").append(newLine);
        }

        if (minIsotopes != null || maxIsotopes != null) {

            output.append("Isotopes: ");

            if (minIsotopes != null) {
                output.append(minIsotopes);
            } else {
                output.append("n");
            }

            output.append("-");

            if (maxIsotopes != null) {
                output.append(maxIsotopes);
            } else {
                output.append("n");
            }

            output.append(".").append(newLine);
        }

        return output.toString();
    }

    /**
     * Returns the minimum number of missed cleavages. Null means no limit.
     *
     * @return the minMissedCleavages
     */
    public Integer getMinMissedCleavages() {
        readDBMode();
        return minMissedCleavages;
    }

    /**
     * Set the minimum number of missed cleavages. Null means no limit.
     *
     * @param minMissedCleavages the minMissedCleavages to set
     */
    public void setMinMissedCleavages(
            Integer minMissedCleavages
    ) {
        writeDBMode();
        this.minMissedCleavages = minMissedCleavages;
    }

    /**
     * Returns the maximum number of missed cleavages. Null means no limit.
     *
     * @return the maxMissedCleavages
     */
    public Integer getMaxMissedCleavages() {
        readDBMode();
        return maxMissedCleavages;
    }

    /**
     * Set the maximum number of missed cleavages. Null means no limit.
     *
     * @param maxMissedCleavages the maxMissedCleavages to set
     */
    public void setMaxMissedCleavages(
            Integer maxMissedCleavages
    ) {
        writeDBMode();
        this.maxMissedCleavages = maxMissedCleavages;
    }
}
